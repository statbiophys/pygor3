#!/usr/bin/python3
import pygor3 as p3
from optparse import OptionParser

# parser = OptionParser(usage=usage)
# parser.add_option("-s", "--specie",
#                   dest="specie",
#                   help="IGoR defined specie")
# parser.add_option("-c", "--chain",
#                   dest="chain",
#                   help="IGoR defined chain")
# # parser.add_option("-f", "--filename",
# #                   metavar="FILE", help="write output to FILE")
# # parser.add_option("-m", "--mode",
# #                   default="intermediate",
# #                   help="interaction mode: novice, intermediate, "
# #                        "or expert [default: %default]")
#
# (options, args) = parser.parse_args()


import argparse
usage = "usage: %prog [options] arg1 arg2"

parser = argparse.ArgumentParser(usage=usage, description="Report the IGoR's models.")
parser.add_argument('integers', metavar='N', type=int, nargs='+',
                    help='an integer for the accumulator')
parser.add_argument('--sum', dest='accumulate', action='store_const',
                    const=sum, default=max,
                    help='sum the integers (default: find the max)')

args = parser.parse_args()


def main():
    """Plot a general report for a model"""
    igor_specie = "mouse"
    igor_chain = "tcr_beta"
    mdl0 = p3.IgorModel.load_default(igor_specie, igor_chain)
    mdl0.parms.plot_Graph()


if __name__ == '__main__':
    main()