# DFT-Super-Cell-Input-File-Generator
A Python script that takes an existing input file for XCrysden or Quantum Espresso and scales the compound based on the parameters you input.
## Objective
Convert 3D Mxene to Supercell
Multiply cell by n times in x and y (ex. multiply cell 4 times in x and 3 times in y) 
Remove atoms between certain reduced z coordinates (ex. remove all atoms between z = 0.33333 & 0.66667) 

## Inputs 
* “How many times should the cell be multiplied in x?”
* “How many times should the cell be multiplied in y?”
* Original Input file

## Outputs
Supercell output file
