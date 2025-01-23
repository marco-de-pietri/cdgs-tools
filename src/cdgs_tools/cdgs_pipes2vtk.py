from .cdgs_class import read_cdgs_file
import argparse


def main():
    """
    this file is the script for simply converting a cdgs with only cylinders to a vtk file.
    """

    # get the cdgs file as input

    parser = argparse.ArgumentParser(description="write a CDGS file to a vtk file")

    parser.add_argument("cdgs", type=str, help="The CDGS file to filter")
    parser.add_argument(
        "-o", "--output", default=None, type=str, help="The output vtk file"
    )

    args = parser.parse_args()

    cdgs_filename_split = args.cdgs.split(".")
    cdgs_filename = cdgs_filename_split[0]

    if args.output is None:
        output_name = f"{cdgs_filename}.vtk"
    else:
        output_name = args.output

    cdgs_object = read_cdgs_file(args.cdgs)

    cdgs_object.write_vtk(output_name)
