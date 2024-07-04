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
        "beads",
        type=existing_file,
        help="Path to a .gtrack file with the list of beads.",
    )
    cli.add_argument(
        "domains",
        type=existing_file,
        help="Path to a BED3+ file with the list of LADs.",
    )

    return cli


def main():
    args = vars(make_cli().parse_args())

    beads = pd.read_table(
        args["beads"], comment="#", names=["chrom", "start", "end", "tid", "radius", "periphery", "edges"]
    )
    lads = pd.read_table(args["domains"], usecols=list(range(3)), names=["chrom", "start", "end"])
    beads = bf.overlap(beads, lads)

    print("##gtrack version: 1.0")
    print("##track type: linked segments\n")
    print("###seqid\tstart\tend\tid\tradius\tperiphery\tedges")

    beads.to_csv(sys.stdout, sep="\t", header=False, index=False)


if __name__ == "__main__":
    main()
