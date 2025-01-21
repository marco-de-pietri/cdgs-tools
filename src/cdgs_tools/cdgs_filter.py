from .cdgs_class import read_cdgs_file
import argparse


def main():
    """
    this file is the script for filtering the CDGS file with a stl.
    """

    # get the a cdgs file and a stl as inputs with argparse

    parser = argparse.ArgumentParser(description="Filter a CDGS file with a stl file")

    parser.add_argument("cdgs", type=str, help="The CDGS file to filter")
    parser.add_argument("stl", type=str, help="The STL file to filter with")

    args = parser.parse_args()

    cdgs_filename_split = args.cdgs.split(".")

    if len(cdgs_filename_split) == 1:
        cdgs_filename = cdgs_filename_split[0]
        cdgs_extension = "cdgs"

    elif len(cdgs_filename_split) == 2:
        cdgs_filename = cdgs_filename_split[0]
        cdgs_extension = cdgs_filename_split[1]
    else:
        raise ValueError("The CDGS file should not contain more than one dot")

    cdgs_object = read_cdgs_file(args.cdgs)

    cdgs_object.write_cdgs(f"{cdgs_filename}_test.{cdgs_extension}")

    #cdgs_object_filtered = cdgs_object.filter_with_stl(args.stl)
