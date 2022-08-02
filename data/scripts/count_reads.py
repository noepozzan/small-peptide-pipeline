#!/usr/bin/env python

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

from argparse import ArgumentParser, RawTextHelpFormatter
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

    __doc__ = "Count reads mapped to transcripts using BAM file and offsets from json file."
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
        "--json",
        dest="json",
        help="Input json file",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--tsv",
        dest="tsv",
        help="Input tsv file",
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
        sys.stdout.write("and parsing {} {}".format(options.bam, os.linesep))

    with open(options.json, 'r') as json_file:
        data = json_file.read()

    cds_start_ind = {}
    cds_end_ind = {}

    with open(options.tsv) as cds_coordinates:
        for line in cds_coordinates:
            sp_line = line.strip().split("\t")
            transcript_id = sp_line[0]
            cds_start_ind[transcript_id] = int(sp_line[2])
            cds_end_ind[transcript_id] = int(sp_line[3])

    count_discarded = 0
    all_transcripts = {}
    offset_obj = json.loads(data)
    bam = pysam.AlignmentFile(options.bam, "rb")  
    for read in bam.fetch():
        if read.is_reverse:
            continue
        transcript = read.reference_name
        if transcript not in all_transcripts:
           all_transcripts[transcript] = [0, 0, 0] 
        read_length = len(read.seq)
        if str(read_length) in offset_obj:
            read_start = read.reference_start
            psite_pos = (read_start + offset_obj[str(read_length)])
            if psite_pos >= cds_start_ind[transcript] and psite_pos <= cds_end_ind[transcript]:
                all_transcripts[transcript][1] = all_transcripts[transcript][1] + 1
            elif psite_pos < cds_start_ind[transcript]:
                all_transcripts[transcript][0] = all_transcripts[transcript][0] + 1
            else:
                all_transcripts[transcript][2] = all_transcripts[transcript][2] + 1
        else:
            count_discarded = count_discarded + 1
        
    w = open(os.path.join(options.outdir, "counts.tsv"), 'w')
    w.write("Transcript" + "\t" + "5' UTR" + "\t" + "CDS" + "\t" + "3' UTR" + "\t" + "CDS length" + "\n")
    for trcpt_key in all_transcripts:
        gene_length = cds_end_ind[trcpt_key]-cds_start_ind[trcpt_key]+1
        w.write(trcpt_key + "\t" + str(all_transcripts[trcpt_key][0]) + "\t" + str(all_transcripts[trcpt_key][1]) + "\t" + str(all_transcripts[trcpt_key][2]) + "\t" + str(gene_length) + "\n")	
    bam.close()
    w.close() 

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!" + os.linesep)
        sys.exit(0)
