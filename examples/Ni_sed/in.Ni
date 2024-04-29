#LAMMPS

units		metal
boundary	p p p 
atom_style	atomic

lattice		fcc 3.52 origin 0.25 0.25 0.25 orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
region		box block 0 10 0 10 0 10 
create_box 	1 box
create_atoms 1 box

pair_style  eam
pair_coeff  * * Ni.eam

timestep 0.001
velocity all create 200 12345 mom yes rot no
compute 1 all saed 0.0251 Ni Kmax 1.70 Zone 0 0 0 c 0.056818 0.056818 0.056818 manual echo
fix 1 all saed/vtk 1 1 1 c_1 file Ni_10x10x10.000.saed.vtk
run 0