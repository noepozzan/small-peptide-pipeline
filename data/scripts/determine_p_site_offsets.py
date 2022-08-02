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
        "--cds_coordinates",
        dest="cds_coordinates",
        help="CDS coordinates tsv table",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--filterReadsForPeriodicity",
        dest="filterReadsForPeriodicity",
        help="Filter reads for periodicity",
        required=False,
        default=False,
        action='store_true'
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
        sys.stdout.write("Parsing {} {}".format(options.cds_coordinates, os.linesep))
        sys.stdout.write("and parsing {} {}".format(options.bam, os.linesep))

    # Open alignment file
    bam = pysam.AlignmentFile(options.bam, "rb")

    ### Calculating P-site offset
    possible_offsets = [10, 11, 12, 13, 14, 15, 16]
    readsDict = {}

    # loop over CDS coordinates
    with open(options.cds_coordinates) as cds_coordinates:
        for line in cds_coordinates:
            sp_line = line.strip().split("\t")

            transcript_id = sp_line[0]
            cds_start = int(sp_line[2]) - 1
            cds_end = int(sp_line[3])

            for read in bam.fetch(transcript_id,
                                  cds_start,
                                  cds_end):

                # keep only + strand
                if read.is_reverse:
                    continue

                read_len = len(str(read.seq))

                # Create a dictionary to score the different offsets for the distinct read lenghts
                if read_len not in readsDict.keys():
                    readsDict[read_len] = {}
                    for offset in possible_offsets:
                        readsDict[read_len][offset] = [0, 0, 0, 0, 0] #1st,2nd,3rd,UTR,CDS
                        readsDict[read_len]["reads"] = []

                read_start = read.reference_start

                for offset in possible_offsets:
                    psite_pos = (read_start + offset) #P-site
                    asite_pos = psite_pos + 3         #A-site
                    rel_pos_psite = psite_pos - cds_start
                    rel_pos_asite = asite_pos - cds_start - 3

                    if options.filterReadsForPeriodicity:
                        #There's NO peak at start codon
                        if rel_pos_asite < 0:
                            readsDict[read_len][offset][3] += 1
                        elif rel_pos_asite > (cds_end - cds_start):
                            readsDict[read_len][offset][3] += 1
                        else:
                            if rel_pos_asite > 60: #ignore first 20 codons
                                frame_rel_pos = rel_pos_asite % 3
                                readsDict[read_len][offset][frame_rel_pos] += 1
                    else:
                        #There's a peak at start codon
                        if rel_pos_psite < 0:
                            readsDict[read_len][offset][3] += 1
                        elif rel_pos_psite > 2:
                            readsDict[read_len][offset][4] += 1
                        elif rel_pos_psite in [0, 1, 2]: #anchor at start codon
                            readsDict[read_len][offset][rel_pos_psite] += 1

    # Close alignment file
    bam.close()

    alignment_offset = {}

    for key, value in readsDict.items():
        best_offset = 0
        best_offset_score = 0
        for offset in possible_offsets:
            raw_counts = readsDict[key][offset][0:3]
            if sum(raw_counts) > 0:
                if ((sum(raw_counts) > 200) \
                and (raw_counts[0] > raw_counts[1]) \
                and (raw_counts[0] > raw_counts[2])):
                    if ((raw_counts[0] > best_offset_score) \
                    and ((not options.filterReadsForPeriodicity) \
                    or (float(raw_counts[0])/sum(raw_counts) > 0.4))):
                        best_offset = offset
                        best_offset_score = raw_counts[0]
        if best_offset != 0:
            alignment_offset[key] = best_offset
            # self.filtered_alignments += (readsDict[rl]["reads"])
            sys.stderr.write("PEAK CALL:" + \
                             str(key) + \
                             ": peak found at position " + \
                             str(best_offset) \
                             + os.linesep)
        else:
            sys.stderr.write("PEAK CALL:" + \
                             str(key) + \
                             ": peak not found" + \
                             os.linesep)

    sys.stdout.write(str(alignment_offset) +  os.linesep)

    with open(os.path.join(options.outdir, "alignment_offset.json"), 'w') as file:
        json.dump(alignment_offset, file)
        # file.write(json.dumps(alignment_offset))

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
