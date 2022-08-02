#!/usr/bin/env python

_version__ = "0.1"
__author__ = "Michi Jackson"
__contact__ = "noe.pozzan@stud.unibas.ch"
__doc__ = "rename files according to a mapping file"

# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

import sys
import os
import csv
from argparse import ArgumentParser, RawTextHelpFormatter


def main():
    """ Main function """

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter
    )

    parser.add_argument(
        "--map_file",
        dest="map_file",
        help="mapping",
        default=False,
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--dir",
        dest="dir",
        help="directory where files live",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--mode",
        dest="mode",
        help="mode in which the files are to be renamed ('rename' or 'reverse')",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        dest="verbose",
        default=False,
        required=False,
        help="Verbose"
    )

    # -------------------------------------------------------------------------
    # get the arguments
    # -------------------------------------------------------------------------

    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    path = options.dir
    with open(options.map_file, 'r') as file:
        my_reader = csv.reader(file, delimiter=',')
        for index, row in enumerate(my_reader):
            if int(index) > 0 and row:
                filename = row[0]
                new_filename = str(row[1]) + "_" + str(row[2]) + ".mzML"
                if options.mode == "rename":
                    os.rename(os.path.join(path, filename),
                              os.path.join(path, new_filename))
                if options.mode == "reverse":
                    os.rename(os.path.join(path, new_filename),
                              os.path.join(path, filename))

    # df = pd.read_csv(options.map_file)
    # path = options.dir
    # for index, row in df.iterrows():
    #     filename = row['file']
    #     new_filename = str(row['experiment']) + "_" + str(row['rep']) + ".mzML"
    #     if options.mode == "rename":
    #         os.rename(os.path.join(path, filename),
    #                   os.path.join(path, new_filename))
    #     if options.mode == "reverse":
    #         os.rename(os.path.join(path, new_filename),
    #                   os.path.join(path, filename))


    # file_map = []
    # path = options.dir
    # for count, filename in enumerate(os.listdir(path)):
    #     filename = filename.strip()
    #     if filename.endswith("sorted_indexed.bam.para.py"):
    #         new_tmp = "_".join(filename.split("_", 2)[:2])
    #         new_filename = new_tmp + "_1"
    #         str_to_app = "mapped {} to {} \n".format(filename, new_filename)
    #         os.rename(os.path.join(path, filename),
    #                   os.path.join(path, new_filename))
    #         file_map.append(str_to_app)
    #     if filename.endswith("folder_sorted_indexed_bam"):
    #         new_tmp = "_".join(filename.split("_", 2)[:2])
    #         new_filename = new_tmp + "_2"
    #         str_to_app = "mapped {} to {} \n".format(filename, new_filename)
    #         file_map.append(str_to_app)
    #         os.rename(os.path.join(path, filename),
    #                   os.path.join(path, new_filename))

    # if options.out:
    #     out_name = options.out_file
    # else:
    #     out_name = "file_map.txt"

    # with open(out_name, "w") as out_file:
    #     out_file.writelines(file_map)


# -----------------------------------------------------------------------------
# Call the Main function and catch Keyboard interrups
# -----------------------------------------------------------------------------


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!" + os.linesep)
        sys.exit(0)
