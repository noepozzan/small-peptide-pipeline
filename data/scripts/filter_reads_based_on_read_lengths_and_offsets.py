#!/usr/bin/env python

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

from argparse import ArgumentParser, RawTextHelpFormatter
import itertools
import sys
import os
import pysam
import json

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------

def main():
    """ Main function """

    __doc__ = "Determine p-site offset"
    __version__ = "0.1"

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter
    )

    parser.add_argument(
        "--bam",
        dest="bam",
        help="Input bam file (sorted and indexed)",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--p_site_offsets",
        dest="p_site_offsets",
        help="P-site offsets (dictionary of read length and offset)",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--bam_out",
        dest="bam_out",
        help="Output bam file",
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
        sys.stdout.write("Parsing {} {}".format(options.p_site_offsets, os.linesep))

    with open(options.p_site_offsets, 'r') as fp:
        p_site_offsets = json.load(fp)

    # cast dict to contain keys and values as ints
    p_site_offsets = {int(k):int(v) for k,v in p_site_offsets.items()}

    if options.verbose:
        sys.stdout.write("and parsing {} {}".format(options.bam, os.linesep))

    # Open alignment file
    bam = pysam.AlignmentFile(options.bam, "rb")

    bam_out = pysam.AlignmentFile(options.bam_out, "wb", header=bam.header)

    for read in bam.fetch():
        # keep only + strand
        if read.is_reverse:
            continue
        read_len = len(str(read.seq))

        if read_len in p_site_offsets:

            a_site_start = p_site_offsets[read_len] + 3

            read.reference_start += a_site_start
            read.cigar = [(0, 3)]
            read.seq = read.seq[a_site_start: 3]

            bam_out.write(read)

    # Close alignment files
    bam_out.close()
    bam.close()


    # alignment_offset = {}
    #
    # for key, value in readsDict.items():
    #     best_offset = 0
    #     best_offset_score = 0
    #     for offset in possible_offsets:
    #         raw_counts = readsDict[key][offset][0:3]
    #         if sum(raw_counts) > 0:
    #             if ((sum(raw_counts) > 200) \
    #             and (raw_counts[0] > raw_counts[1]) \
    #             and (raw_counts[0] > raw_counts[2])):
    #                 if ((raw_counts[0] > best_offset_score) \
    #                 and ((not options.filterReadsForPeriodicity) \
    #                 or (float(raw_counts[0])/sum(raw_counts) > 0.4))):
    #                     best_offset = offset
    #                     best_offset_score = raw_counts[0]
    #     if best_offset != 0:
    #         alignment_offset[key] = best_offset
    #         # self.filtered_alignments += (readsDict[rl]["reads"])
    #         sys.stderr.write("PEAK CALL:" + \
    #                          str(key) + \
    #                          ": peak found at position " + \
    #                          str(best_offset) \
    #                          + os.linesep)
    #     else:
    #         sys.stderr.write("PEAK CALL:" + \
    #                          str(key) + \
    #                          ": peak not found" + \
    #                          os.linesep)
    #
    # sys.stdout.write(str(alignment_offset) +  os.linesep)
    #
    # with open(os.path.join(options.outdir, "alignment_offset.json"), 'w') as file:
    #     file.write(json.dumps(alignment_offset))

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
