#!/usr/bin/env python

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

from argparse import ArgumentParser, RawTextHelpFormatter
import sys
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------

def main():
    """ Main function """

    __doc__ = "Plot length of reads"
    __version__ = "0.1"

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter
    )

    parser.add_argument(
        "--sam",
        dest="sam",
        help="Input sam file",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--outdir",
        dest="outdir",
        help="Output directory",
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
        sys.stdout.write("Creating output directory: {} {}".format(options.outdir, os.linesep))

    if not os.path.exists(options.outdir):
        os.mkdir(options.outdir)

    if options.verbose:
        sys.stdout.write("Parsing {} {}".format(options.sam, os.linesep))

    lengths = []
    with open(options.sam) as sam:
        for line in sam:
            if line[0] != "@":
                lengths.append(len(line.split("\t")[9]))

    fig = plt.figure()
    plt.hist(lengths, bins=30)
    plt.xlabel('Length')
    plt.ylabel('Number of cases')
    # plt.axis([0, 20000, None, None])
    plt.savefig(os.path.join(options.outdir, "read_length_histogram.pdf"))



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
