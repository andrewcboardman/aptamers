from Bio import SeqIO
import argparse

# Take name of read file as input at command line

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--infile', type=str, action='store', dest='infile',default='test.fastq'
	help='Name of raw reads file (FASTQ format)')
args = parser.parse_args()

# Select for reads starting with the primer sequence and trim them to remove the primer

def trim_primer(record, primer):
    if record.seq.startswith(primer):
        return record[len(primer):]
    else:
        return record

trimmed_reads = (trim_primer(record, "GATGACGGTGT") for record in SeqIO.parse(infile, "fastq"))

# Remove reads if any base has a quality lower than 10

good_reads = (rec for rec in trimmed_reads if min(rec.letter_annotations["phred_quality"]) >= 10)

# Write to a new file

SeqIO.write(good_reads, 'filter.fastq', "fastq")
