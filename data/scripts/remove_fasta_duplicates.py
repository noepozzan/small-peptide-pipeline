#!/usr/bin/env python

_version__ = "0.1"
__author__ = "Noe Pozzan"
__contact__ = "noe.pozzan@stud.unibas.ch"
__doc__ = """
            Parse a fasta file and remove duplicate entries.
          """

# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

import sys
import os
from argparse import ArgumentParser, RawTextHelpFormatter
from Bio import SeqIO

def main():
    """ Main function """

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter
    )

    parser.add_argument(
        "--input",
        dest="input",
        help="input fasta file with duplicate entries",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--output",
        dest="output",
        help="output fasta file withOUT duplicate entries",
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

    if options.verbose:
        sys.stdout.write(
            "Parsing params file: {} {}".format(
                options.input,
                os.linesep
            )
        )

    # if entry is already in seen list, don't put it to the records
    seen = []
    records = []

    for record in SeqIO.parse(options.input, "fasta"):
        if str(record.seq) not in seen:
            seen.append(str(record.seq))
            records.append(record)


    #writing to the output fasta file
    SeqIO.write(records, options.output, "fasta")

# -----------------------------------------------------------------------------
# Call the Main function and catch Keyboard interrups
# -----------------------------------------------------------------------------


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!" + os.linesep)
        sys.exit(0)
