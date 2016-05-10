2D Diffusion Equation Solver, Version 1.0

Authors: Gabriel Kaufman and Kyle Leverett

Synopsis: this code takes a two-dimensional system, defined by cells and materials as well as a source term, and calculates the neutron flux throughout the system with vacuum boundaries along the leftmost and bottom edges of the system and reflecting boundaries along the rightmost and top edges of the system, using assumptions made to derive the two-dimensional neutron diffusion equation as well as Gauss-Seidel iterative solutions to calculate the final flux. Note that the system must always be rectangular in nature, and failure to define the entire mesh may result in error.

Code usage: python diff-2d.py [-h] [-p] [-g] [-s] infile [outdir] 
File arguments:
  infile        Input file.
  outdir        Output directory for storing output file (optional).

Optional arguments:
  -h, --help    Show the help message and exit
  -p, --print   Re-print input, skip over actual calculation of flux. The output file will only contain parsed input file
                information unless "--sgraph" is flagged.
  -g, --graph   Graph flux as heat map. This plot will not appear until the code has completed.
  -s, --sgraph  Graph source strength as heat map. This plot will not appear until the code has completed. Additionally,
                the source vector used in the flux calculation will be written to file.

During runtime, if the code terminates due to an error or keyboard interrupt ("Ctrl+C") the output file will NOT be written. However, if a keyboard interrupt occurs during the Gauss-Seidel iteration stage, it will only terminate that process prematurely and continue with the rest of the code.

Expected input:
The input file resembles that for MCNP. However, given the more limited geometry, the cell and surface definitions are adjustedi to match a simple mesh. The file starts with a title; the following lines beginning with the character 'C' are considered comments, and will thus be ignored (capitalization, for this and following characters, is not important). The first non-commented line determines the total width, number of cells in the x-direction, height, and number of cells in the y-direction (in order). I.e.,

   Mesh definition format: W Nx H Ny

Non-commented lines following this input will be considered "macrocells." Since it is an inherently rectangular system, all macrocells are likewise rectangular in nature. As such, the lines start with the cell number, followed by the material number, and then the minimum and maximum x and y coordinates. 

   Cell definition format: Cell# Mat# x-min y-min x-max y-max

The user must keep the boundaries "sensical"; i.e., the max values should not surpass the total width of the system. Overlaps of different materials are considered combinations of their respective properties, while overlaps of the same material will simply be overwritten (i.e., it will combine the two regions). In order to leave the system in a solvable state, the total area of the cells should encompass the entire area of the mesh. Else, the program may prematurely end due to a singular coefficient matrix. 
The next set of inputs will be the materials. Lines giving material properties start with the character 'M', followed by the material number (as per MCNP material definitions), the value of the diffusion coefficient D, and the macroscopic cross section Sigma_a. 

   Material definition format: M[Mat#] D Sig_a

Material overlaps from cell definitions use the inverse sum of D^(-1) and the sum of Sigma_a for each material.

Finally, the input gives the user the option to define a source, designating the definition with a line starting with 'S'. This source can be defined as a constant or a function for a given x and y range. For example, using the keyword "COS_X" gives the source strength as a cosine function in the x-direction, and constant along the y-direction. The available functions are:

   CONST  A0          => A
   LIN_X  A0 A1       => A0 + A1*x
   LIN_Y  A0 A1       => A0 + A1*y
   QUAD_X A0 A1 A2    => A0 + A1*x + A2*x^2
   QUAD_Y A0 A1 A2    => A0 + A1*y + A2*y^2
   COS_X  A0 A1 A2 A3 => A0*cos(A1*pi*x + A2*pi) + A3
   COS_Y  A0 A1 A2 A3 => A0*cos(A1*pi*y + A2*pi) + A3

The function will always be taken as the absolute value of the desired function, such that no negative flux is given. After the function has been defined, the user has the option to determine the range. This takes the form of four numbers following the function parameters; if fewer than four numbers follow, the remaining dimensions are assumed to encompass the minimum/maximum dimensions for the grid. Multiple source definitions will combine as a product on an elementwise basis in the grid.

   Source definition fomrat: S FUNC A0 [A1 A2 A3] [x-min y-min x-max y-max]

Code status: fully operational.
