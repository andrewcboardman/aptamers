from pysster.utils import annotate_structures
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile', type=str, action='store', dest='infile',	help='Input filename')
parser.add_argument('-o', '--outfile', type=str, action='store', dest='infile',	help='Output filename')
args = parser.parse_args()
annotate_structures(args.infile,args.outfile)