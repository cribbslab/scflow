import os 
import sys
import pysam
import re
import cgatcore.experiment as E

def main(argv=sys.argv):

    ##################
    # Option parsing 
    ##################

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-i", "--in", dest="infile", type ="string" ,
                      help="BAM file without CB and UB tags. UMI and CB located in read name")

    parser.add_option("-o", "--out", dest="outfile", type ="string" ,
                      help="BAM file with CB and UB tags")


    (options, args) = E.start(parser) 

    bam_in = options.infile
    bam_out = options.outfile

    # Bam files open using pysam
    bamfile = pysam.AlignmentFile(bam_in, "rb")
    reads_barcodes_tags = pysam.AlignmentFile(bam_out, "wb", template=bamfile)
    
    # Write to new bam file
    for read in bamfile.fetch():
        reg = re.search('.*_(.*)_(.*)',str(read.query_name))
        cb = reg[1]
        ub = reg[2]
        read.tags = read.tags + [("CB",cb)] + [("UB", ub)]

        reads_barcodes_tags.write(read)

    bamfile.close()
    reads_barcodes_tags.close()

if __name__ == "__main__":
    sys.exit(main())
