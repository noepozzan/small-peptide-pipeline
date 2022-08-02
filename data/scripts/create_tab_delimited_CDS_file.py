#!/usr/bin/env python

__version__ = "0.1"
__author__ = "Foivos Gypas"
__contact__ = "foivos.gypas@unibas.ch"
__doc__ = "Parse a gtf file and a transcripts sequence file (as generated " + \
          "from gffread and create a tsv file containing the following: " + \
          "transcript id, gene id, start codon, stop codon)"

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

import sys
import os
import HTSeq
from Bio import SeqIO
from argparse import ArgumentParser, RawTextHelpFormatter

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Custom functions
# -----------------------------------------------------------------------------
def split_gff_fasta_header(header):

    header = header.strip().split(" ")

    transcript_id = header[0].strip("\'")
    CDS = header[2].strip("CDS=").split("-")

    CDS_start = CDS[0]
    CDS_stop = CDS[1]

    return transcript_id, CDS_start, CDS_stop

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------


def main():
    """ Main function """

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter
    )

    parser.add_argument(
        "--gtf",
        dest="gtf",
        help="Input GTF file in ENSEMBL format",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--fasta",
        dest="fasta",
        help="Selected transcripts fasta file",
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
        help="Verbose"
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
        sys.stdout.write(
            "Parsing gtf file: {} {}".format(
                options.gtf,
                os.linesep
            )
        )


    # dictionary of keys: transcript_id values: gene_id
    transcripts = {}

    # parse gtf file
    gtf_file = HTSeq.GFF_Reader(options.gtf)

    for gtf_line in gtf_file:

        if gtf_line.type == 'transcript':
            transcript_id = gtf_line.attr['transcript_id']
            gene_id = gtf_line.attr['gene_id']

            if transcript_id not in transcripts.keys():
                transcripts[transcript_id] = gene_id

    # parse fasta file
    w = open(options.out, 'w')
    for entry in SeqIO.parse(options.fasta, 'fasta'):
        transcript_id, CDS_start, CDS_stop = split_gff_fasta_header(entry.description)
        if transcript_id in transcripts.keys():
            w.write("\t".join([transcript_id,
                               transcripts[transcript_id],
                               CDS_start,
                               CDS_stop + os.linesep]))
    w.close()



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
