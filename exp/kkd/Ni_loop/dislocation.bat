atomsk Ni_crystal.lmp cfg
atomsk Ni_crystal.cfg -dislocation loop 0.5*box 0.5*box 0.5*box X 20 2.489 0 0 0.33 Ni_loop.cfg
atomsk Ni_loop.cfg lmp
@REM plane_normal radius_in_Angstrom Burgers_vector Poisson_ratio