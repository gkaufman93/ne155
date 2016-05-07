import numpy as np
from math import *
from mesh_comp import comp_flux
import datetime as dt

def main(inf, outd, flag):
    # Begin log_str string (output) and print first lines."
    version_number = 0
    log_str   = "2D Diffusion Solver V{:1.1f}\n".format(version_number)
    log_str  += "Authors: Gabriel Kaufman and Kyle Leverett\n"
    log_str  += "Input file: " + inf + '\n'
    print "2D Diffusion Solver V{:1.1f}\n".format(version_number)
    print "Authors: Gabriel Kaufman and Kyle Leverett\n"
    print "Input file: " + inf + '\n'

    mat_dict  = {}
    cell_dict = {}
    # Parse input file
    with open(inf,'r') as f:
        # Get title from first line of input file
        sr = (f.readline()).strip(' \t\n\r')
        log_str += "Input title: " + sr + '\n'
        print "Input title: " + sr + '\n'
        st = str(dt.datetime.now())
        log_str += "\nStart time: " + st + '\n'
        print "Start time: " + st + '\n'
        line_num = 2
        sr = f.readline().upper()
        while len(sr) > 0 and sr[0] == 'C':
            sr = f.readline().upper()
            line_num += 1
        if sr == '':
            raise InputFileError('Effectively empty file.')
        elif not (sr[0].isdigit() or sr[0] == 'M' or sr[0] == 'S'):
            raise InputFileError('Incorrect line in file, line\
                                 {}'.format(line_num))
        l = sr.split()
        if len(l) != 4:
            raise ParseError("System definition has incorrect number of \
                             arguments, line {}.".format(line_num))
        # Get (in order) width, Nx, height, Ny of system
        # Then calculate delta and epsilon 
        W  = float(l[0])
        Nx = int(float(l[1]))
        d  = W/Nx
        H  = float(l[2])
        Ny = int(float(l[3]))
        e  = H/Ny
        if not all(float(i) > 0 for i in l):
            raise LessThanZeroError("One or more arguments have values that are\
                                    less than zero, line {}".format(line_num))
        log_str += "\nMesh definition:\n"
        log_str += "\t\tWidth = {: >16.4f}".format(W)
        log_str += "\t\tNum. mesh points = {: >10}".format(Nx)
        log_str += "\n\t\tHeight = {: >15.4f}".format(H)
        log_str += "\t\tNum. mesh points = {: >10}".format(Ny)
        log_str += "\n\n"
        print "Mesh successfully parsed."
        # Initialize map of materials
        mat_map = [['0' for x in range(Nx+1)] for y in range(Ny+1)]
        mat_dict['0'] = (0,0)
        sr = f.readline().upper()
        line_num += 1
        if len(sr) == 0:
            raise InputFileError('End of file reached with only mesh\
                                 definitions.')

        # Get cell values in order: cell number, cell material, min x, min y,
        # max x, max y
        #
        # If cells overlap and have different materials, their properties
        # will be combined
        log_str += "Cells:\n"
        log_str += "Cell" + 4*' '+"Material" + 4*' '+"X Range" + 23*' '+\
                "Y Range" + 23*' ' + "\n\n"
        while len(sr) > 0 and sr[0] != 'M':
            while len(sr) > 0 and sr[0] == 'C':
                sr = f.readline().upper()
                line_num += 1
            if sr == '':
                raise InputFileError('End of file reached before material and\
                                     source definitions.')
            elif not (sr[0].isdigit() or sr[0] == 'M' or sr[0] == 'S'):
                raise InputFileError('Incorrect line in file, line\
                                     {}'.format(line_num))
            if (sr[0] == 'M' or sr[0] == 'S') and len(cell_dict) == 0:
                raise InputFileError('No cells defined in input file.')
            elif sr[0] == 'M' or sr[0] == 'S':
                break
            l = sr.split()
            if len(l) != 6:
                raise ParseError("System definition has incorrect number of \
                                 arguments, line {}.".format(line_num))
            cell_num = int(float(l[0]))
            cell_mat = l[1]
            cell_xmin = float(l[2])
            cell_ymin = float(l[3])
            cell_xmax = float(l[4])
            cell_ymax = float(l[5])
            print "Parsing cell {}".format(cell_num)+"..."
            if not all(float(i) >= 0 for i in l):
                raise LessThanZeroError("One or more arguments have values that are\
                                    less than zero, line {}".format(line_num))
            elif cell_xmax <= cell_xmin or cell_ymax <= cell_ymin:
                raise BoundError('Bounds are nonsensical, line\
                                 {}'.format(line_num))
            cell_dict[cell_num] = (cell_mat, cell_xmin, cell_ymin, cell_xmax,\
                                   cell_ymax)
            log_str += "{: <8}".format(cell_num)
            log_str += "{: <12}".format(cell_mat)
            log_str += "{: <15.4f}".format(cell_xmin)
            log_str += "{: <15.4f}".format(cell_xmax)
            log_str += "{: <15.4f}".format(cell_ymin)
            log_str += "{: <15.4f}".format(cell_ymax)
            log_str += '\n'
            cx_min = int(cell_xmin/d)
            cx_it  = int((cell_xmax - cell_xmin)/d)+1
            cy_min = int(cell_ymin/e)
            cy_it  = int((cell_ymax - cell_ymin)/e)+1
            # Fill material map for easy retrieval during computation
            for i in range(cx_it):
                for j in range(cy_it):
                    if mat_map[i][j] == '0':
                        mat_map[i][j] = cell_mat
                        if cell_mat not in mat_dict:
                            mat_dict[cell_mat] = (0,0)
                    elif cell_mat not in mat_map[i][j].split():
                        mat_map[i][j] += ' ' + cell_mat
                        if mat_map[i][j] not in mat_dict:
                            mat_dict[mat_map[i][j]] = (0,0)
            sr = f.readline().upper()
            line_num += 1
            if len(sr) == 0:
                raise InputFileError('End of file reached before material and\
                                     source definitions.')

        # Get materials and material properties (D, sig_a)
        print "Cells successfully parsed."
        mat_dict_temp = {}
        log_str += "\nMaterials:\n\n"
        log_str += "Material" + 4*' ' + "D" + 14*' ' + "Sig_a" + 10*' '\
                + "\n\n"
        while sr[0] != 'S':
            while len(sr) > 0 and sr[0] == 'C':
                sr = f.readline().upper()
                line_num += 1
            if sr == '':
                raise InputFileError('End of file reached before source\
                                     definitions.')
            elif not (sr[0] == 'M' or sr[0] == 'S'):
                raise InputFileError('Incorrect line in file, line\
                                     {}'.format(line_num))
            if sr[0] == 'S' and len(mat_dict_temp) == 0:
                raise InputFileError('Material definitions not found or\
                                     incorrectly organized in input file.')
            elif sr[0] == 'S':
                break
            l = sr.split()
            if len(l) != 3:
                raise ParseError("System definition has incorrect number of \
                                 arguments, line {}.".format(line_num))
            mat_num = l[0][1:]
            mat_D   = float(l[1])
            mat_Sig = float(l[2])
            l[0] = l[0][1:]
            print "Parsing material {}".format(mat_num) + "..."
            if not all(float(i) > 0 for i in l):
                raise LessThanZeroError("One or more arguments have values that are\
                                    less than zero, line {}".format(line_num))
            mat_dict_temp[mat_num] = (mat_D, mat_Sig)
            log_str += "{:<12}".format(mat_num)
            log_str += "{: <15.4f}".format(mat_D)
            log_str += "{: <15.4f}".format(mat_Sig)
            log_str += '\n'
            sr = f.readline().upper()
            line_num += 1

        print "Successfully parsed materials."
        print "Adjusting material dictionary..."
        # Set material properties in mat_dict. Combine if necessary. Raise
        # MaterialError if it is not described in precious block.
        mats = mat_dict.keys()
        for k in mats:
            l  = k.split() 
            mat_D = 0
            mat_sig = 0
            if l == ['0']:
                continue
            for m in l:
                try:
                    t = mat_dict_temp[m]
                except KeyError:
                    raise MaterialError("Material " + m + " not included in material \
                    definitions. Check cell definitions.")
                mat_D += 1/t[0]
                mat_sig += t[1]
            mat_D = 1/mat_D
            mat_dict[k] = (mat_D, mat_sig)

        # Iterate over source definitions. Sources are combined on elementwise
        # basis. However, only one key word can be defined per line. Number of
        # required arguments depends on the keyword; for example, a constant
        # keyword only requires one argument; 'COS_X' requires three keywords.
        # 
        # See README for more details, and usable functions.
        # Also possible to define domain for source with additional arguments.
        # 
        # NOTE: All functions, regardless of parameters, will always be
        # positive. An absolute value function is applied in order to keep 
        # source values sensical (i.e., source value cannout be less than zero).
        print "Material dictionary complete. Parsing source(s)..."
        S = np.zeros(((Nx+1)*(Ny+1),1))
        Smat = np.zeros((Nx+1, Ny+1))
        log_str += "Source types: \n\n"
        while len(sr) > 0:
            while len(sr) > 0 and sr[0] == 'C':
                sr = f.readline().upper()
                line_num += 1
            if sr == '' and log_str == "Source types: \n\n":
                raise InputFileError('No sources defined.')
            elif sr == '':
                break
            elif not sr[0] == 'S':
                raise InputFileError('Incorrect line in file, line\
                                     {}'.format(line_num))
            l = sr.split()
            l.pop(0)
            s_key = l.pop(0)
            log_str += s_key + "\n" 
            if s_key == "CONST":
                args = float(l.pop(0))
            elif s_key == "LIN_X":
                args = []
                for i in range(0,2):
                    args.append(l.pop(0))
            elif s_key == "LIN_Y":
                args = []
                for i in range(0,2):
                    args.append(l.pop(0))
            elif s_key == "QUAD_X":
                args = []
                for i in range(0,3):
                    args.append(l.pop(0))
            elif s_key == "QUAD_Y":
                args = []
                for i in range(0,3):
                    args.append(l.pop(0))
            elif s_key == "COS_X":
                args = []
                for i in range(0,4):
                    args.append(l.pop(0))
            elif s_key == "COS_Y":
                args = []
                for i in range(0,4):
                    args.append(l.pop(0))
            if len(l) > 0 and (not all(float(i) >= 0 for i in l)):
                raise LessThanZeroError("One or more arguments have values that are\
                                    less than zero, line {}".format(line_num))
            if len(l) == 0:
                sx_min = 0
                sy_min = 0
                sx_max = Nx
                sy_max = Ny
            elif len(l) == 1:
                sx_min = int(float(l[0])/d)
                sy_min = 0
                sx_max = Nx
                sy_max = Ny
            elif len(l) == 2:
                sx_min = int(float(l[0])/d)
                sy_min = int(float(l[1])/e)
                sx_max = Nx
                sy_max = Ny
            elif len(l) == 3:
                sx_min = int(float(l[0])/d)
                sy_min = int(float(l[1])/e)
                sx_max = int(float(l[2])/d)
                sy_max = Ny
            elif len(l) == 4:
                sx_min = int(float(l[0])/d)
                sy_min = int(float(l[1])/e)
                sx_max = int(float(l[2])/d)
                sy_max = int(float(l[3])/e)
            if (sx_min >= sx_max) or (sy_min >= sy_max):
                raise BoundError('Bounds are nonsensical, line\
                                 {}'.format(line_num))
            for i in xrange(sx_min, sx_max+1):
                for j in xrange(sy_min, sy_max+1):
                    if Smat[i,j] == 0:
                        Smat[i,j] = func(s_key, i*d, j*e, args)
                    else:
                        Smat[i,j] *= func(s_key, i*d, j*e, args)
            sr = f.readline().upper()
            line_num += 1

        print "Sources parsed. Creating source vector for flux calculation..."
        for ind in xrange(Nx+1,(Nx+1)*(Ny+1)):
            i = ind%(Nx+1)
            j = ind/(Nx+1)
            if i == 0:
                continue
            elif j != 0:
                S1 = Smat[i,j]
            else:
                S1 = 0
            if i != Nx and j != 0:
                S2 = Smat[i+1,j]
            else:
                S2 = 0
            if i != Nx and j != Ny:
                S3 = Smat[i+1,j+1]
            else:
                S3 = 0
            if j != Ny:
                S4 = Smat[i,j+1]
            else:
                S4 = 0
            S[i+j*(Nx+1)] = 0.25*d*e*(S1+S2+S3+S4)

    log_str += "Source output:\n\n"
    log_str += 'X' + 11*' ' + '\t\t' + 'Y' + 11*' ' + '\t\t' + 'Source' + 5*' ' + '\n\n'
    for ind in range(len(S)):
        px = (ind%(Nx+1))*d
        py = (ind/(Nx+1))*e
        log_str += "{:.6E}".format(px) + '\t\t'
        log_str += "{:.6E}".format(py) + '\t\t'
        log_str += "{:.6E}".format(float(S[ind])) + '\t\t'
        log_str += '\n'
    print "Source vector complete."

    log_str += '\n\n' + "End of file reached. Parse complete"
    log_str += '\n\n'
    print "Input file successfully parsed."
    if flag:
        log_str += "Input print mode flag, computation not carried out.\n\n"
    else:
        print "Computing flux..."
        cf = comp_flux(d, Nx, e, Ny, mat_map, mat_dict, S)
        F,numIt = cf.calcFlux()
        log_str += "Computation done!"
        log_str += "\nNumber of Gauss-Seidel iterations required: {}".format(numIt)
        log_str += "\n\nFlux output:\n\n"
        log_str += 'X' + 11*' ' + '\t\t' + 'Y' + 11*' ' + '\t\t' + 'Flux' + 8*' ' +\
                '\t\t' + '\n\n'
        for ind in range(len(F)):
            px = (ind%(Nx+1))*d
            py = (ind/(Nx+1))*e
            log_str += "{:.6E}".format(px) + '\t\t'
            log_str += "{:.6E}".format(py) + '\t\t'
            log_str += "{:.6E}".format(float(F[ind])) + '\t\t'
            log_str += '\n'
        log_str += '\n\n'
    if outd != '':
        if '.' in inf and inf.rfind('.') < inf.rfind('/'):
            inf = inf[:inf.rfind('.')]
        while '/' in inf:
            inf = inf[inf.find('/')+1:]
        outf = outd + inf + '.out'
    else:
        if '.' in inf:
            inf = inf[:inf.rfind('.')]
        outf = inf + '.out'
    log_str += "Process finished.\n" 
    et = str(dt.datetime.now())
    log_str += "End time: " + et + '\n'
    with open(outf,'w') as f:
        f.write(log_str)
    print "Data written to output file: " + outf
    print "Process finished, end time: " + et
    return;

def func(argkey, x, y, args):
    if argkey == "CONST":
        return abs(args)
    elif argkey == "LIN_X":
        return abs(args[0] + x*args[1])
    elif argkey == "LIN_Y":
        return abs(args[0] + y*args[1])
    elif argkey == "QUAD_X":
        return abs(args[0] + args[1]*x + args[2]*pow(x,2))
    elif argkey == "QUAD_Y":
        return abs(args[0] + args[1]*y + args[2]*pow(y,2))
    elif argkey == "COS_X":
        return abs(args[0]*cos(args[1]*x + args[2]) + args[3])
    elif argkey == "COS_Y":
        return abs(args[0]*cos(args[1]*y + args[2]) + args[3])
    else:
        raise SourceKeyError("Source keyword \""+argkey+"\" is not a possible argument.")

class ParseError(Exception):
    def __init__(self, value):
        self.value = value

class InputFileError(Exception):
    def __init__(self, value):
        self.value = value

class MaterialError(Exception):
    def __init__(self, value):
        self.value = value

class SourceKeyError(Exception):
    def __init__(self, value):
        self.value = value
                    
class LessThanZeroError(Exception):
    def __init__(self, value):
        self.value = value

class BoundError(Exception):
    def __init__(self, value):
        self.value = value
