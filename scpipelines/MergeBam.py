#!/usr/bin/env python

import pysam 
from optparse import OptionParser

parser = OptionParser()
parser.add_option("--unmapped", dest = "unmapped_file",
                  help = "Path to unmapped trimmed, filtered, tagged BAM file. <samplename>_polyA_filtered.bam")
parser.add_option("--aligned", dest = "aligned_file", 
                  help = "Path to BAM file aligned by STAR, missing CB/UMI tags")
parser.add_option("-o", "--output", dest = "output", help = "Merged output file")

(options, args) = parser.parse_args()

unmapped = options.unmapped_file
mapped = options.aligned_file
output = options.output

mapped_bamfile = pysam.AlignmentFile(mapped, "rb", check_sq=False)
unmapped_bamfile = pysam.AlignmentFile(unmapped, "rb", check_sq=False)
merged_reads = pysam.AlignmentFile(output, "wb", template=mapped_bamfile)
    
for read_map, read_unmap in zip(mapped_bamfile.fetch(until_eof=True), unmapped_bamfile.fetch(until_eof=True)):

    barcode_tags = read_unmap.get_tags()

    existing_tags = read_map.get_tags()

    read_map.set_tags(existing_tags + barcode_tags)

    merged_reads.write(read_map)


merged_reads.close()
mapped_bamfile.close()
unmapped_bamfile.close()

