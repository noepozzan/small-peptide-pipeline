#!/usr/bin/env python

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

from argparse import ArgumentParser, RawTextHelpFormatter
import sys
import os
from Bio import SeqIO

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------

def main():
    """ Main function """

    __doc__ = "Find overepresented sequences from a SAM file."
    __version__ = "0.1"

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter
    )

    parser.add_argument(
        "--fastq",
        dest="fastq",
        help="Input fastq file",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--oligos",
        dest="oligos",
        help="Input oligos file",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--out",
        dest="out",
        help="Output file",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        dest="verbose",
        default=False,
        required=False,
        help="Be verbose"
    )

    parser.add_argument(
        '--version',
        action='version',
        version=__version__
    )

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # get the arguments
    # -------------------------------------------------------------------------
    try:
        options = parser.parse_args()
    except Exception:
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    if options.verbose:
        sys.stdout.write("arsing {} {}".format(options.oligos, os.linesep))

    oligos = {}
    with open(options.oligos) as fp:
        for line in fp:
            oligo = line.strip()
            if oligo not in oligos:
                oligos[oligo] = 0

    if options.verbose:
        sys.stdout.write("Parsing {} {}".format(options.fastq, os.linesep))

    for line in SeqIO.parse(options.fastq, "fastq"):
        sequence = str(line.seq)
        for oligo, number_of_oligos in oligos.items():
            if oligo in sequence:
                oligos[oligo] += 1
    if options.verbose:
        sys.stdout.write("Writing {} {}".format(options.out, os.linesep))

    fh = open(options.out, "w")
    for w in sorted(oligos, key=oligos.get, reverse=True):
        fh.write(str(w) + "\t" + str(oligos[w]) + "\n")
    fh.close()

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Call the Main function and catch Keyboard interrups
# -----------------------------------------------------------------------------


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!" + os.linesep)
        sys.exit(0)
