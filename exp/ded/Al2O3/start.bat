..\..\..\bin\AAVDP_win.exe --ded .\Al2O3.vasp -q 0.8 -o .\Al2O3_spot.001.ded
..\..\..\bin\AAVDP_win.exe --ded .\Al2O3.vasp -q 0.8 -o .\Al2O3.001.ded -gauss
..\..\..\bin\AAVDP_win.exe --ded .\Al2O3.vasp -z 2 1 0 -fn 2 1 0 -q 0.8 -o .\Al2O3_spot.210.ded -rotate -rotate_x 0 0 1 -rotate_y -0.44721360 0.89442719 0
..\..\..\bin\AAVDP_win.exe --ded .\Al2O3.vasp -z 2 1 0 -fn 2 1 0 -q 0.8 -o .\Al2O3.210.ded -gauss -gauss_dx 0.00375 -gauss_sig 0.0075 -rotate -rotate_x 0 0 1 -rotate_y -0.44721360 0.89442719 0