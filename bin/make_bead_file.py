#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT

import argparse
import pathlib

import bioframe as bf
import pandas as pd


def make_cli() -> argparse.ArgumentParser:
    def existing_file(arg):
        if (file := pathlib.Path(arg)).exists():
            return file

        raise FileNotFoundError(arg)

    cli = argparse.ArgumentParser()

    cli.add_argument(
        "bedpe",
        type=existing_file,
        help="Path to a BEDPE with the list of significant interactions.",
    )
    cli.add_argument(
        "domains",
        type=existing_file,
        help="Path to a BED3+ file with the list of domains used to generate significant interactions.",
    )
    cli.add_argument(
        "chrom.sizes",
        type=existing_file,
        help="Path to the .chrom.sizes file for the reference genome used to call significant interactions.",
    )

    return cli


def import_beads(path_to_beads: pathlib.Path, chrom_sizes: pd.DataFrame) -> pd.DataFrame:
    beads = pd.read_table(path_to_beads, usecols=list(range(3)), names=["chrom", "start", "end"])

    beads = pd.concat([beads, bf.complement(beads, chrom_sizes)])
    return bf.sort_bedframe(beads)[["chrom", "start", "end"]]


def generate_gtrack(beads: pd.DataFrame, sig_interactions):
    records = {}

    for chrom1, start1, end1, chrom2, start2, end2 in sig_interactions.itertuples(index=False):
        r1 = (chrom1, start1, end1)
        r2 = (chrom2, start2, end2)
        records.setdefault(r1, []).append(r2)
        records.setdefault(r2, []).append(r1)

    for chrom, start, end in beads.itertuples(index=False):
        r = (chrom, start, end)
        if r not in records:
            records[r] = ()

    print("##gtrack version: 1.0")
    print("##track type: linked segments\n")
    print("###seqid\tstart\tend\tid\tradius\tperiphery\tedges")

    for k, v in records.items():
        k_bed = "\t".join((str(x) for x in k))
        k_ucsc = f"{k[0]}:{k[1]}-{k[2]}"

        if len(v) == 0:
            v = ["."]
        else:
            v = [f"{chrom}:{start}-{end}" for chrom, start, end in v]

        print("\t".join((k_bed, k_ucsc, "0.2", ".", ";".join(v))))


def main():
    args = vars(make_cli().parse_args())
    chrom_sizes = pd.read_table(args["chrom.sizes"], names=["chrom", "end"])
    chrom_sizes["start"] = 0

    sig_interactions = pd.read_table(
        args["bedpe"],
        usecols=list(range(6)),
        names=["chrom1", "start1", "end1", "chrom2", "start2", "end2"],
    )
    beads = import_beads(args["domains"], chrom_sizes)

    generate_gtrack(beads, sig_interactions)


if __name__ == "__main__":
    main()
