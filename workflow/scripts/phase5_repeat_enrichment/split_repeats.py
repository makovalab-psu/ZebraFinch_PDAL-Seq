#!/usr/bin/env python3
"""Split one repeat annotation into one BED file per repeat class.

Same job as the manuscript's Catagorize_CenSat.R, rewritten because ours has
to do three extra things that R script did not:

  * filter to the contigs that exist in our genome and clamp the coordinates
    to the contig length. `bedtools shuffle` aborts the whole job on a single
    out-of-range record, and these annotations were made against Linnea's
    matZ extraction rather than against our make_PDAL-Seq_fasta.sh extraction;

  * sort by .fai RANK, not ASCII. The .fai order is chr1_mat, chr1A_mat,
    chr2_mat, ... which is neither lexicographic nor numeric;

  * check the class names it finds against the list hardcoded in
    workflow/Snakefile, and fail loudly on a mismatch. Snakemake's outputs are
    that hardcoded list, so a class that appeared or was renamed upstream
    would otherwise show up as a missing-output error with no explanation.

Python rather than awk because this writes up to 54 output files from one
input, and the pipeline convention is that awk never does that -- BSD awk
silently drops output files when it has many open at once.

Usage:
  split_repeats.py --input in.bed --chrom-sizes sizes --outdir dir \\
                   --report report.txt --source TE \\
                   (--column 5 | --fixed-name CEN) [--all-name all_TE] \\
                   --expected c1,c2,...
"""

import argparse
import os
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True)
    p.add_argument("--chrom-sizes", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--report", required=True)
    p.add_argument("--source", required=True,
                   help="label written into the report, e.g. TE or Satellite")
    p.add_argument("--column", type=int, default=0,
                   help="1-based column holding the class name")
    p.add_argument("--fixed-name", default="",
                   help="put every record in one class of this name")
    p.add_argument("--all-name", default="",
                   help="also write every record to this extra class")
    p.add_argument("--expected", required=True,
                   help="comma separated class names the Snakefile declares")
    return p.parse_args()


def read_chrom_sizes(path):
    """contig -> (rank, length), in file order."""
    sizes = {}
    with open(path) as handle:
        for rank, line in enumerate(handle):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            sizes[fields[0]] = (rank, int(fields[1]))
    if not sizes:
        sys.exit("ERROR: {} held no contigs".format(path))
    return sizes


def main():
    args = parse_args()

    if (args.column > 0) == bool(args.fixed_name):
        sys.exit("ERROR: pass exactly one of --column and --fixed-name")

    sizes = read_chrom_sizes(args.chrom_sizes)
    expected = [c for c in args.expected.split(",") if c]
    if args.all_name and args.all_name not in expected:
        expected.append(args.all_name)
    expected_set = set(expected)

    # class -> list of (rank, chrom, start, end). Held in memory on purpose:
    # the largest of these annotations is ~290,000 records, which is a few
    # tens of MB, and sorting in memory avoids 54 concurrent sort processes.
    records = {}

    n_in = 0
    n_off_contig = 0
    n_malformed = 0
    n_empty_after_clamp = 0
    unexpected = {}

    with open(args.input) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            n_in += 1
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 3:
                n_malformed += 1
                continue

            chrom = fields[0]
            if chrom not in sizes:
                n_off_contig += 1
                continue

            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                n_malformed += 1
                continue

            if args.fixed_name:
                name = args.fixed_name
            else:
                if len(fields) < args.column:
                    n_malformed += 1
                    continue
                name = fields[args.column - 1].strip()
                if not name:
                    n_malformed += 1
                    continue

            rank, length = sizes[chrom]
            if start < 0:
                start = 0
            if end > length:
                end = length
            if end <= start:
                n_empty_after_clamp += 1
                continue

            if name not in expected_set:
                unexpected[name] = unexpected.get(name, 0) + 1
                continue

            records.setdefault(name, []).append((rank, chrom, start, end))
            if args.all_name:
                records.setdefault(args.all_name, []).append(
                    (rank, chrom, start, end))

    if unexpected:
        sys.stderr.write(
            "ERROR: {} holds {} class name(s) that workflow/Snakefile does "
            "not declare.\n".format(args.input, len(unexpected)))
        for name in sorted(unexpected, key=lambda k: -unexpected[k]):
            sys.stderr.write("       {:>10} records  {}\n".format(
                unexpected[name], name))
        sys.stderr.write(
            "       Add them to the list in workflow/Snakefile (or remove the "
            "ones that\n       have gone away) and re-run.\n")
        sys.exit(1)

    missing = sorted(expected_set - set(records))
    if missing:
        sys.stderr.write(
            "ERROR: workflow/Snakefile declares {} class(es) that survived "
            "nothing in {}:\n".format(len(missing), args.input))
        for name in missing:
            sys.stderr.write("       {}\n".format(name))
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(os.path.dirname(args.report) or ".", exist_ok=True)

    # Records are NOT merged. `bedtools intersect -a windows -b class.bed`
    # downstream emits one fragment per overlapping annotation, so a window
    # covered by many elements of a class contributes many sampling units --
    # which is what the methods' "each annotation type was subsampled at
    # random" means. Merging first would silently turn that into uniform
    # sampling over windows. See README section 07.
    with open(args.report, "w") as report:
        for name in sorted(records):
            rows = sorted(records[name], key=lambda r: (r[0], r[2], r[3]))
            path = os.path.join(args.outdir, name + ".bed")
            with open(path, "w") as out:
                for _, chrom, start, end in rows:
                    out.write("{}\t{}\t{}\n".format(chrom, start, end))
            contigs = len({r[1] for r in rows})
            bp = sum(r[3] - r[2] for r in rows)
            report.write("{}\t{}\t{}\t{}\t{}\n".format(
                args.source, name, len(rows), contigs, bp))

    sys.stderr.write(
        "{}: read {} records, wrote {} classes to {}\n".format(
            args.input, n_in, len(records), args.outdir))
    if n_off_contig:
        sys.stderr.write("  {} records were on contigs not in {}\n".format(
            n_off_contig, args.chrom_sizes))
    if n_empty_after_clamp:
        sys.stderr.write("  {} records were empty after clamping to the "
                         "contig length\n".format(n_empty_after_clamp))
    if n_malformed:
        sys.stderr.write("  {} records were malformed and skipped\n".format(
            n_malformed))


if __name__ == "__main__":
    main()
