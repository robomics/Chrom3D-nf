#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import sys

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

    cli.add_argument(
        "--lads",
        type=existing_file,
        help="Path to a BED3+ file with the list of LADs",
    )

    return cli


def import_beads(path_to_beads: pathlib.Path, chrom_sizes: pd.DataFrame) -> pd.DataFrame:
    beads = pd.read_table(path_to_beads, skiprows=1, usecols=list(range(3)), names=["chrom", "start", "end"])
    return pd.concat([beads, bf.complement(beads, chrom_sizes)])[["chrom", "start", "end"]]


def intersect_with_lads(beads: pd.DataFrame, lads: pd.DataFrame) -> pd.DataFrame:
    df = bf.count_overlaps(beads, lads)
    df.loc[df["count"] != 0, "periphery"] = "1"
    return df[beads.columns.tolist()]


def generate_gtrack(beads: pd.DataFrame, sig_interactions: pd.DataFrame, lads: pathlib.Path):
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

    data = []
    for k, v in records.items():
        k_ucsc = f"{k[0]}:{k[1]}-{k[2]}"

        if len(v) == 0:
            v = ["."]
        else:
            v = [f"{chrom}:{start}-{end}" for chrom, start, end in v]

        data.append([*k, k_ucsc, "0.2", ".", ";".join(v)])

    beads = pd.DataFrame(data, columns=["chrom", "start", "end", "tid", "radius", "periphery", "edges"])

    if lads is not None:
        lads = pd.read_table(lads, usecols=list(range(3)), names=["chrom", "start", "end"])
        return intersect_with_lads(beads, lads)

    return beads


def main():
    args = vars(make_cli().parse_args())
    chrom_sizes = pd.read_table(args["chrom.sizes"], names=["chrom", "end"])
    chrom_sizes["start"] = 0

    sig_interactions = pd.read_table(
        args["bedpe"],
        usecols=list(range(6)),
    )
    beads = import_beads(args["domains"], chrom_sizes)

    beads = generate_gtrack(beads, sig_interactions, args["lads"])
    print(beads, file=sys.stderr)
    beads = bf.sort_bedframe(beads, chrom_sizes)

    print("##gtrack version: 1.0")
    print("##track type: linked segments")
    print("###seqid\tstart\tend\tid\tradius\tperiphery\tedges")

    beads.to_csv(sys.stdout, sep="\t", header=False, index=False)


if __name__ == "__main__":
    main()
