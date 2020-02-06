#!/usr/bin/env python

from optparse import OptionParser
import os

parser = OptionParser()
parser.add_option("--input", dest = "input_bam",
                  help = "Input bam file. Format velocyto.dir/<sample_name>.bam")
parser.add_option("--barcode_suffix", dest = "barcode_suffix", default = "_barcodes.tsv",
                  help = "Suffix of barcode whitelist file, [default= %default]")
parser.add_option("--config", dest = "config", default = "config_desc.xml" ,
                  help = "Config xml file, [default= %default]")

(options, args) = parser.parse_args()

input_bam = options.input_bam
barcode_suffix = options.barcode_suffix
generic_xml_file = options.config

sample_name = os.path.basename(input_bam).replace(".merged.aligned.coord.bam", "")
barcode_whitelist = "./" + sample_name + barcode_suffix
sample_config = "dropest.dir/" + sample_name + "_config_desc.xml"


with open(generic_xml_file, 'r+') as fin:
    with open(sample_config, "w") as fout:
        for line in fin:
            if "<barcodes_file>" in line:
                start = line.find(">") + 1
                end = line.find("<", start)
                changed = line[:start] + barcode_whitelist + line[end:]
                fout.write(changed)
            else:
                fout.write(line)



