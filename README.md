Program Overview
This program is designed to identify and calculate the number of interfacial water molecules between two proteins.

Program Functionality
Protein Definition: The program assumes there are two proteins in the system, where the atomic index range for Protein A is from 1 to 2186, Protein B’s atomic index range is from 2187 to 4372, and the atomic index range for water molecules is from 4373 to 539389.

Interfacial Water Definition:

A 2 nm spherical region centered at the surface of each protein.

A cylindrical volume with a radius of Rh/2 nm, with its axis aligned along the line connecting the centers of mass of the two proteins.

Input Data:
Molecular dynamics trajectory file (md.xtc), which contains the coordinates of the system at different time steps.

Center-of-mass information file (mass.txt), which contains the positions of the proteins' centers of mass.

Output:
The program outputs the number of interfacial water molecules for each time step, saving the results in the file number_of_interfacial_water.dat.
