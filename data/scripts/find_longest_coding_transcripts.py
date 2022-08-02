#!/usr/bin/env python

__version__ = "0.1"
__author__ = "Foivos Gypas"
__contact__ = "foivos.gypas@unibas.ch"
__doc__ = "Parse a gtf file and find the longest expressed protein " + \
          "coding transcripts for each gene."

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

import sys
import os
import HTSeq
from argparse import ArgumentParser, RawTextHelpFormatter

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Custom classes
# -----------------------------------------------------------------------------


class Gene(object):

    """Gene object"""

    def __init__(self, gene_id):
        """Init object"""
        self.gene_id = gene_id
        self.transcripts = []

    def add_transcript(self, transcript):
        """Add transcript in the list of transcripts"""
        if transcript not in self.transcripts:
            self.transcripts.append(transcript)

    def get_all_transcripts(self):
        """Get all transcripts of that gene"""
        return self.transcripts

    def get_longest_coding_transcript(self):
        """"""
        selected_transcript = None
        for transcript in self.transcripts:
            if selected_transcript is None:
                selected_transcript = transcript
                continue
            if (transcript.get_CDS_length() > selected_transcript.get_CDS_length()):
                selected_transcript = transcript
        return selected_transcript


class Transcript(object):

    """Transcript object"""

    def __init__(self, transcript_id, gene_id):
        """Init object"""
        self.transcript_id = transcript_id
        self.gene_id = gene_id
        self.CDS_length = 0

    def get_CDS_length(self):
        """Get CDS length of the transcript"""
        return self.CDS_length

    def update_CDS_length(self, length):
        """Update CDS length"""
        self.CDS_length += length

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
        "--out",
        dest="out",
        help="GTF output file",
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
    except(Exception):
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

    # dictionary of genes
    genes = {}

    # dictionary of transcripts
    transcripts = {}

    # parse gtf file
    gtf_file = HTSeq.GFF_Reader(options.gtf)

    for gtf_line in gtf_file:

        if gtf_line.type == 'CDS':

            transcript_id = gtf_line.attr['transcript_id']
            gene_id = gtf_line.attr['gene_id']
            length = abs(gtf_line.iv.end - gtf_line.iv.start)

            if gene_id not in genes.keys():
                genes[gene_id] = Gene(gene_id=gene_id)

            if transcript_id not in transcripts.keys():
                transcripts[transcript_id] = Transcript(transcript_id, gene_id)
                genes[gene_id].transcripts.append(transcripts[transcript_id])

            transcripts[transcript_id].update_CDS_length(length)

    # keep the longest coding transcripts per gene
    selected_transcsripts_ids = {}

    for gene_id, gene_object in genes.items():

        if options.verbose:
            print(gene_object.gene_id)
            for tr in gene_object.transcripts:
                print(tr.transcript_id, tr.get_CDS_length())

        longest_transcript = gene_object.get_longest_coding_transcript()

        selected_transcsripts_ids[longest_transcript.transcript_id] = \
            longest_transcript.transcript_id

    # parse again the gtf file and keep only the longest transcripts that
    # were selected previously
    if options.verbose:
        sys.stdout.write(
            "Re-parse gtf file: {} and write out: {} {}".format(
                options.gtf,
                options.out,
                os.linesep
            )
        )

    w = open(options.out, 'w')
    for gtf_line in gtf_file:
        if (
            gtf_line.type == 'transcript' or
            gtf_line.type == 'exon' or
            gtf_line.type == 'CDS' or
            gtf_line.type == 'start_codon' or
            gtf_line.type == 'stop_codon'
        ):

            if (gtf_line.attr['transcript_id'] in selected_transcsripts_ids):
                w.write(gtf_line.get_gff_line())
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
