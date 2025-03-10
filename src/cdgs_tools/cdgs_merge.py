from .cdgs_class import read_cdgs_file, CDGS

import argparse


def main():
    """
    this file is the script for merging cdgs files.
    """


    parser = argparse.ArgumentParser(
        description="merge cdgs files"
    )

    parser.add_argument("cdgs", nargs='+', type=str, help="enter two or more cdgs files to merge")

    args = parser.parse_args()

    if len(args.cdgs) < 2:
        parser.error("at least two cdgs files are required to merge")

    cdgs_filename = "merged.cdgs"
    vtk_filename = "cdgs_merged.vtk"

    cdgs_merge_object = CDGS()

    for cdgs_file in args.cdgs:

        cdgs_object = read_cdgs_file(cdgs_file)
        cdgs_merge_object.extend_cdgs(cdgs_object)

    cdgs_merge_object.write_cdgs(cdgs_filename)

    cdgs_merge_object.write_vtk(vtk_filename)

    return
