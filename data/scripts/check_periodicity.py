#!/usr/bin/env python

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

from argparse import ArgumentParser, RawTextHelpFormatter
import sys
import os
import pysam
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import json

def main():
    """ Main function """

    __doc__ = "Checks periodicity around start and stop codon."
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
        "--codnum",
        dest="codnum",
        help="Input codon coverage",
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

    codon_number_for_plot = int(options.codnum)
    bam = pysam.AlignmentFile(options.bam, "rb")

    with open(options.json, 'r') as json_file:
        data = json_file.read()
    offsets = json.loads(data)
    cds_start_ind = {}
    cds_end_ind = {}

    with open(options.tsv) as cds_coordinates:
        for line in cds_coordinates:
            sp_line = line.strip().split("\t")
            transcript_id = sp_line[0]
            cds_start_ind[transcript_id] = int(sp_line[2]) - 1
            cds_end_ind[transcript_id] = int(sp_line[3])

    p_site_signal_codon = {}
    p_site_signal_codon_stop = {}
    p_site_signal_transcipts = {}
    p_site_signal_transcipts_stop = {}
    for read in bam.fetch():
        #if read.reference_name not in candidate:
            #continue
        read_start = read.reference_start
        cds_start = cds_start_ind[read.reference_name]
        cds_end = cds_end_ind[read.reference_name]
        if read.is_reverse:
            continue
        if str(len(str(read.seq))) in offsets:
            codon_position = read_start - cds_start + offsets[str(len(read.seq))]
            codon_position_stop = read_start - cds_end + offsets[str(len(read.seq))]
            if codon_position in p_site_signal_codon:
                p_site_signal_codon[codon_position] = p_site_signal_codon[codon_position] + 1
                old_list = p_site_signal_transcipts[codon_position]
                old_list.append(read.reference_name)
                p_site_signal_transcipts[codon_position]=old_list
            else:
                p_site_signal_codon[codon_position] = 1
                p_site_signal_transcipts[codon_position] = (read.reference_name).split()

            if codon_position_stop in p_site_signal_codon_stop:
                p_site_signal_codon_stop[codon_position_stop] = p_site_signal_codon_stop[codon_position_stop] + 1
                old_list = p_site_signal_transcipts_stop[codon_position_stop]
                old_list.append(read.reference_name)
                p_site_signal_transcipts_stop[codon_position_stop]=old_list
            else:
                p_site_signal_codon_stop[codon_position_stop] = 1
                p_site_signal_transcipts_stop[codon_position_stop] = (read.reference_name).split()

    w = open(os.path.join(options.outdir, 'Periodicity_Analysis_Start_Ribo_Seq.txt'), "wt")
    for key, value in p_site_signal_codon.items():
       w.write(str(key) + "\t" + str(value) + "\n")
    w.close()
    bam.close()
    
    w1 = open(os.path.join(options.outdir, 'Periodicity_Analysis_Stop_Ribo_Seq.txt'), "wt")
    for key, value in p_site_signal_codon_stop.items():
       w1.write(str(key) + "\t" + str(value) + "\n")
    w1.close()

    codon_occupancy=[]
    codon_positions=list(range(-codon_number_for_plot, codon_number_for_plot+1))
    for i in codon_positions:
        if i in p_site_signal_codon:
            codon_occupancy.append(p_site_signal_codon[i])
        else:
            codon_occupancy.append(0)
    plt.figure(1)
    plt.plot(codon_positions, codon_occupancy)
    plt.title('Periodicity Around Start Codon')
    plt.xlabel('Codon Position')
    plt.ylabel('P-site Signal')
    y_limits=plt.gca().get_ylim()
    plt.savefig(os.path.join(options.outdir, 'periodicity_start.pdf'))
    codon_occupancy_stop=[]
    codon_positions_stop=list(range(-codon_number_for_plot, codon_number_for_plot+1))
    for i in codon_positions_stop:
        if i in p_site_signal_codon_stop:
            codon_occupancy_stop.append(p_site_signal_codon_stop[i])
        else:
            codon_occupancy_stop.append(0)
    plt.figure(2)
    plt.plot(codon_positions_stop, codon_occupancy_stop)
    plt.title('Periodicity Around Stop Codon')
    plt.xlabel('Codon Position')
    plt.ylabel('P-site Signal')
    plt.ylim(y_limits[0], y_limits[1]) 
    plt.savefig(os.path.join(options.outdir, 'periodicity_stop.pdf'))

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!" + os.linesep)
        sys.exit(0)
