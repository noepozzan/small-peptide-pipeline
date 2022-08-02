#!/usr/bin/env python

__version__ = "0.1"
__author__ = "Meric Ataman"
__contact__ = "meric.ataman@unibas.ch"
__doc__ = "Parse a riboTISH file and map predicted sORFs to peptides."

# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

import sys
import os
import pdb
from Bio import SeqIO
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from argparse import ArgumentParser, RawTextHelpFormatter


def main():
    """ Main function """

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter
    )

    parser.add_argument(
        "--ribo_pred",
        dest="ribo_pred",
        help="riboTISH predict output files in txt format, in a list",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--fasta",
        dest="fasta",
        help="fasta file from gffread containing annotation",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--out",
        dest="out",
        help="fasta output file name",
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
            "Parsing riboTISH file: {} {}".format(
                options.ribo_pred,
                os.linesep
            )
        )

    ##########
    all_transcript_sequences = {}
    for entry in SeqIO.parse(options.fasta, 'fasta'):
        seq = str(entry.seq)
        trp = str(entry.id)
        if trp in all_transcript_sequences:
            disp("problem!")
            pdb.set_trace()
        all_transcript_sequences[trp] = seq

    # pdb.set_trace()
    aa_seq_small_peptides = []
    all_small_peptides = []

    # Go through each sample, identify the
    # -'UTR/dORF-3'UTR TisType and from them
    # generate peptides which are candidate small peptides.

    # for sam in options.ribo_pred:
    all_uorfs = pd.read_csv(options.ribo_pred, sep='\t')
    utr_uorf = all_uorfs.loc[(all_uorfs['TisType'] == "5'UTR")
                             | (all_uorfs['TisType'] == "3'UTR")]
    for ind, row in utr_uorf.iterrows():
        transcript_all = row['Tid']
        transcript_all = transcript_all.split('.')
        transcript = transcript_all[0]
        # just to check if we have the transcript sequence
        if transcript in all_transcript_sequences:
            transcript_seq = all_transcript_sequences[transcript]
        else:
            print(transcript)
            continue
        nuc_seq_uorf = Seq(transcript_seq[int(row['Start']):row['Stop'] + 1])
        if str(nuc_seq_uorf[0:3]).upper() != row['StartCodon']:
            # just to check if they match
            print(transcript)
            continue
        peptide = str(nuc_seq_uorf.translate())
        peptide = peptide.replace('*', '')
        if peptide in all_small_peptides:
            # just to check there are no duplicate entries
            continue
        all_small_peptides.append(peptide)
        aa_seq_uorf = SeqRecord(Seq(peptide), id=row['Symbol']
                                + '_' + transcript + '_Nuc_Pos:'
                                + str(row['Start']) + '-' + str(row['Stop'])
                                + '_' + row['TisType'],
                                description=row['Symbol']
                                + '_' + transcript + '_small_peptide_'
                                + str(row['Start']) + '-' + str(row['Stop'])
                                + '_' + row['TisType'])
        aa_seq_small_peptides.append(aa_seq_uorf)

    with open(options.out, "w") as output_handle:
        SeqIO.write(aa_seq_small_peptides, output_handle, "fasta")


# -----------------------------------------------------------------------------
# Call the Main function and catch Keyboard interrups
# -----------------------------------------------------------------------------


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!" + os.linesep)
        sys.exit(0)
