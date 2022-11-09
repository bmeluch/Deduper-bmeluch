#!/usr/bin/env python
import argparse
import bioinfo
import re

# assumption: the input SAM file is sorted

############################## FUNCTIONS ##############################

def pos_adjust(stranded: bool, lmap: int, cig: str) -> int:
    """Parses CIGAR string and returns a corrected leftmost start position
    based on strandedness and soft clipping."""

    # if read is on the - strand
    if stranded:
        # find numbers before each of the CIGAR letters and sum them
        M: int = sum([int(i) for i in re.findall(r'(\d+)M', cig)])
        D: int = sum([int(i) for i in re.findall(r'(\d+)D', cig)])
        N: int = sum([int(i) for i in re.findall(r'(\d+)N', cig)])
        S: int = 0
        # if there is soft clipping at the end, add that
        soft = re.search(r'(\d+)S\Z', cig)
        if soft:
            S = int(soft.group(1))
        return lmap + M + D + N + S
    # if read is on the + strand
    else:
        # if soft clipping at the beginning of the read, subtract it from pos
        soft = re.search(r'^(\d+)S', cig)
        if soft:
            return lmap - int(soft.group(1))
        else:
            return lmap

############################## GET INPUT FILES ##############################

parser = argparse.ArgumentParser(description="Get file paths for a SORTED SAM input file, output SAM file, and file containing the list of UMIs.")
parser.add_argument("-f", "--file", help="Absolute path to SORTED SAM input file", required=True)
parser.add_argument("-o", "--outfile", help="Absolute path to output SAM File", required=True)
parser.add_argument("-u", "--UMI", help="File containing the list of UMIs", required=True)
args = parser.parse_args()
inpath: str = args.file
outpath: str = args.outfile
umipath: str = args.UMI

# test files
# inpath = "/projects/bgmp/bmeluch/bioinfo/Bi624/Deduper-bmeluch/test.sam"
# outpath = "/projects/bgmp/bmeluch/bioinfo/Bi624/Deduper-bmeluch/Part3/test_out.sam"
# inpath = "/projects/bgmp/bmeluch/bioinfo/Bi624/Deduper-bmeluch/Part1/unit-test_input.sam"
# outpath = "/projects/bgmp/bmeluch/bioinfo/Bi624/Deduper-bmeluch/Part1/test_output.sam"
# umipath = "/projects/bgmp/bmeluch/bioinfo/Bi624/Deduper-bmeluch/Part1/unit-test_UMIs.txt"

############################## IMPORT UMIS ##############################

# Read file with list of known UMIs, save UMI sequences into a set
umis = set()
with open(umipath, 'r') as umifile:
    for line in umifile:
        i = line.strip()
        if bioinfo.validate_base_seq(i):
            umis.add(i)

############################## COUNTERS AND HOLDERS ##############################

# counts of things
headercount: int = 0
uniqcount: int = 0
wrongumicount: int = 0
dupsremoved: int = 0

# counts of unique reads per chrom
chroms: dict = {}

# set holding past read info
archive = set()

############################## MAIN LOGIC ##############################

with open(inpath, 'r') as infile, open(outpath, 'w') as outfile:
    for line in infile:

        # Write out header lines
        if line[0] == "@":
            outfile.write(line)
            headercount += 1

        # Once at real alignment lines, split by whitespace into component fields
        else:
            # Pull out important data from the resulting list
            alignment = line.strip().split()
            umi: str = alignment[0][-8:]
            stranded: bool = ((int(alignment[1]) and 16) == 16)
            chrom: str = alignment[2]
            lmap: int = int(alignment[3])
            cig: str = alignment[5]

            # If UMI is not recognized, do not output read, move to next line
            if umi not in umis:
                wrongumicount += 1
                continue

            # Adjust start position based on strandedness, CIGAR string
            lmap = pos_adjust(stranded, lmap, cig)

            # current read data as tuple
            current: tuple = (umi, stranded, chrom, lmap)

            # Duplicate check? UMI, strand, chromosome, leftmost mapping position
            if current in archive:
                dupsremoved += 1
            else:
                if chrom in chroms:
                    chroms[chrom] += 1
                else:
                    # when you hit a new chromosome, clear archive to save memory
                    chroms[chrom] = 1
                    archive = set()
                # write out read, count a unique read written out, save
                outfile.write(line)
                uniqcount += 1
                archive.add(current)

############################## COUNTS OUTPUT ##############################

print("Header lines:\t\t\t", headercount)
print("Total unique reads written out:\t", uniqcount)
print("Wrong UMIs:\t\t\t", wrongumicount)
print("Duplicates removed:\t\t", dupsremoved)
print("Unique reads per chromosome:")
print(chroms)