import numpy as np
from mesh_main import main
import sys
import argparse as ap

parser = ap.ArgumentParser(description='Parse diff-2d.py command line \
                           arguments.')
parser.add_argument('inf',metavar='infile',help='Input file.')
parser.add_argument('outd', nargs='?', metavar='outdir',default='',help='Output directory\
                    (optional).')
parser.add_argument('-p','--print',dest='flag',action='store_true',help='Re-print input,\
                    skip over actual calculation of flux.')
args = parser.parse_args()
argd = vars(args)
inf  = argd['inf']
outd = argd['outd']
flag = argd['flag']

#try:
main(inf,outd,flag)
#except Exception as e:
    #print e
    #sys.exit(1)



