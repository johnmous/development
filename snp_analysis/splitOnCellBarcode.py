import pysam
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description='Split multiplexed files according to barcode')
parser.add_argument('--barcodes', type=str, help='Text file with all 96 cell barcodes')
parser.add_argument('--multiplexed_bam', type=str, help='Multiplexed bam file')
parser.add_argument('--output_dir', type=str, help='PAth where the demultiplexed files are stored')
args = parser.parse_args()

# Load the 96 barcodes into a list
f = open(args.barcodes, "r")
barcodes = []
for line in f:
    barcodes.append(line[:-1])

# Open one output file per barcode, put it in dict with barcode as key
template_file = pysam.AlignmentFile(args.multiplexed_bam, 'rb', threads=2)
brcd_to_file = {}
for barcode in barcodes:
    f = pysam.AlignmentFile(args.output_dir+'/'+barcode+'.bam', 'wb', template=template_file)
    brcd_to_file[barcode] = f

# Go through the reads, if the read has the CB tag, write it to the corresponding file
for read in template_file.fetch():
    if read.has_tag('CB'):
        cb = read.get_tag('CB')
        brcd_to_file[cb].write(read)

# Close all files opened
template_file.close()
for f in brcd_to_file.values():
    f.close()

Path(args.output_dir+'/demultiplexing_done.txt').touch()
