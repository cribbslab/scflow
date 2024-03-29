import Bio.SeqIO.FastaIO as FastaIO
import cgatcore.iotools as iotools
import argparse
import sys
import logging
import re

# ########################################################################### #
# ###################### Set up the logging ################################# #
# ########################################################################### #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("tr2gene.py")


def main(argv=None):
    if argv is None:
        argv = sys.argv

parser = argparse.ArgumentParser()

parser.add_argument("--fasta", default=None, type=str,
                        help="An ensembl cdna fasta file with the properly formatted header")

parser.add_argument("--output", default=None, type=str,
                        help="output tr2gene file")


args = parser.parse_args()

L.info("args:")
print(args)

outf = iotools.open_file(args.output, "w")

with iotools.open_file(args.fasta, "r") as handle:
    for record in FastaIO.FastaIterator(handle):
        
        description = record.description
        trans = description.split(" ")[0]
        m = re.search('gene:(\S+)', description)
        gene = m.group(1)
        try:
            x = re.search('gene_symbol:(\S+)', description)
            symbol = x.group(1)
        except Exception:
            pass

        outf.write("%s\t%s\t%s\n" % (trans, gene, symbol))

outf.close()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
