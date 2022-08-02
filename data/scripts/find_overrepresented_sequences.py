#!/usr/bin/env python

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

from argparse import ArgumentParser, RawTextHelpFormatter
import sys
import os
import pysam

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
        "--sam",
        dest="sam",
        help="Input bam file (sorted and indexed)",
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
        sys.stdout.write("and parsing {} {}".format(options.sam, os.linesep))

    # Init dictionaries
    # observed_reads: Store observed_reads to count multimappers only once
    # sequence_counts: Number of times a sequence was observed

    # key:read id, value: read_id
    observed_reads = {}
    # key: sequence, value: count
    sequence_counts = {}

    # Open alignment file
    sam = pysam.AlignmentFile(options.sam, "rb")

    for read in sam.fetch():

        sequence = read.seq
        qname = read.qname

        if qname in observed_reads:
            continue
        else:
            observed_reads[qname] = qname

        if sequence not in sequence_counts:
            sequence_counts[sequence] = 1
        else:
            sequence_counts[sequence] += 1

    sam.close()

    fh = open(options.out, "w")
    for w in sorted(sequence_counts, key=sequence_counts.get, reverse=True):
        fh.write(str(w) + "\t" + str(sequence_counts[w]) + "\n")
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
