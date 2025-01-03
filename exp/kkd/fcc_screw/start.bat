..\..\..\bin\AAVDP_win.exe --kkd .\screw.vasp -q 1.5 -o .\screw.ked3 --ked3
..\..\..\bin\AAVDP_win.exe --kkd .\screw.ked3 -z 1 0 0 -rx 0.075 -ry 0.075 -t 0.1 -px 1000 -py 1000 -o .\screw.100.kkd --rotate -x 0.00000000 0.47609129 -0.33664738 -y 0.00000000 0.33664738 0.47609129 --scale -i 0 350
..\..\..\bin\AAVDP_win.exe --kkd .\screw.ked3 -z 0 1 0 -rx 0.075 -ry 0.075 -t 0.1 -px 1000 -py 1000 -o .\screw.010.kkd --rotate -x 0 0 1 -y 1 0 0 --scale -i 0 130
..\..\..\bin\AAVDP_win.exe --kkd .\screw.ked3 -rx 0.075 -ry 0.075 -t 0.1 -px 1000 -py 1000 -o .\screw.001.kkd --scale -i 0 350