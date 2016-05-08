import numpy as np
from mesh_main import main
import sys
import argparse as ap
import matplotlib.pyplot as pp

parser = ap.ArgumentParser(description='Parse diff-2d.py command line \
                           arguments.')
parser.add_argument('inf',metavar='infile',help='Input file.')
parser.add_argument('outd', nargs='?', metavar='outdir',default='',help='Output directory\
                    (optional).')
parser.add_argument('-p','--print',dest='flag',action='store_true',help='Re-print input,\
                    skip over actual calculation of flux.')
parser.add_argument('-g','--graph',dest='gflag',action='store_true',help='Graph\
                    flux as heat map.')
parser.add_argument('-s','--sgraph',dest='sflag',action='store_true',help='Graph\
                    source strength as heat map.')
args = parser.parse_args()
argd = vars(args)
inf  = argd['inf']
outd = argd['outd']
flag = argd['flag']
sflag = argd['sflag']
gflag = argd['gflag']

try:
    main(inf,outd,flag,sflag,gflag)
except Exception as e:
    print e
    #raise
    #sys.exit(1)

if sflag or gflag:
    pp.show()

