# aavdp

#### 介绍
Automatic Analysis of Virtual Diffraction Pattern （AAVDP）for artificial atomistic structures

#### 软件架构
1.  X-ray diffraction (XRD)
2.  Neutron diffraction (NED)
3.  Kinematic electron diffraction (KED)
4.  Kinematic Kikuchi diffraction (KKD)
5.  Dynamical electron diffraction (DED)
6.  Dynamical Kikuchi diffraction (DKD)
7.  Radial distribution function (RDF)
8.  Static structure factor (SSF)

#### 安装教程
1.  win:
	>>cp makefile.win makefile
	>>make
2.  mac:
	>>cp makefile.mac makefile
	>>make
3.  linux:
	>>cp makefile.linux makefile
	>>make

#### 使用说明
1.  xrd:
	>>./bin/AAVDP_win.exe --xrd ./exp/xrd/NaCl/NaCl.vasp
	>>./bin/AAVDP_win.exe --xrd ./exp/xrd/FeCo/FeCo.lmp -e Fe Co
	>>./bin/AAVDP_win.exe --xrd ./exp/xrd/FeCo/FeCo.lmp -e Fe Co -l 1.5400
	>>./bin/AAVDP_win.exe --xrd ./exp/xrd/FeCo/FeCo.lmp -e Fe Co -l 1.5400 -scherrer

	>>./bin/AAVDP_win.exe --xrd ./exp/xrd/NaCl/NaCl.vasp -dw 1.72 1.41 -l 1.54056 -2t 20 90 -o ./exp/xrd/NaCl/NaCl1.xrd
	>>./bin/AAVDP_win.exe --xrd ./exp/xrd/NaCl/NaCl.vasp -dw 1.72 1.41 -l 1.54439 -2t 20 90 -o ./exp/xrd/NaCl/NaCl2.xrd

	>>./bin/AAVDP_win.exe --xrd ./exp/xrd/ZnO/ZnO.vasp -2t 0 80 -o ./exp/xrd/ZnO/ZnO_line.xrd
	>>./bin/AAVDP_win.exe --xrd ./exp/xrd/ZnO/ZnO.vasp -2t 0 80 -o ./exp/xrd/ZnO/ZnO.xrd -scherrer -scherrer_m 0.5 -scherrer_d 294 -scherrer_d2t 0.02

	>>./bin/AAVDP_win.exe --xrd ./exp/xrd/FeCo/FeCo.lmp -e Fe Co -dw 0.5500 0.5500 -2t 0 150 -l 1.5400 -o ./exp/xrd/FeCo/FeCo.xrd
	>>./bin/AAVDP_win.exe --xrd ./exp/xrd/FeCo/FeCo_random.lmp -e Fe Co -dw 0.5500 0.5500 -2t 0 150 -l 1.5400 -c 10 10 10 -o ./exp/xrd/FeCo/FeCo_random.xrd
2.  ned:
	>>./bin/AAVDP_win.exe --ned ./exp/ned/LaCrGe3/LaCrGe3.vasp
	>>./bin/AAVDP_win.exe --ned ./exp/ned/FeCo/FeCo.lmp -e Fe Co
	>>./bin/AAVDP_win.exe --ned ./exp/ned/FeCo/FeCo.lmp -e Fe Co -l 1.5400
	>>./bin/AAVDP_win.exe --ned ./exp/ned/FeCo/FeCo.lmp -e Fe Co -l 1.5400 -scherrer

	>>./bin/AAVDP_win.exe --ned ./exp/ned/LaCrGe3/LaCrGe3.vasp -dw 0.7646 0.1028 0.2458 -2t 10 80 -l 2.43955 -o ./exp/ned/LaCrGe3/LaCrGe3_line.ned
	>>./bin/AAVDP_win.exe --ned ./exp/ned/LaCrGe3/LaCrGe3.vasp -dw 0.7646 0.1028 0.2458 -2t 10 80 -l 2.43955 -o ./exp/ned/LaCrGe3/LaCrGe3.ned -scherrer -scherrer_d 400 -scherrer_d2t 0.02
	>>./bin/AAVDP_win.exe --ned ./exp/ned/FeCo/FeCo.lmp -e Fe Co -dw 0.5500 0.5500 -2t 0 150 -l 1.5400 -o ./exp/ned/FeCo/FeCo.ned
	>>./bin/AAVDP_win.exe --ned ./exp/ned/FeCo/FeCo_random.lmp -e Fe Co -dw 0.5500 0.5500 -2t 0 150 -l 1.5400 -c 10 10 10 -o ./exp/ned/FeCo/FeCo_random.ned
3.  ked:
	>>./bin/AAVDP_win.exe --ked ./exp/ked/Al2O3/Al2O3.vasp
	>>./bin/AAVDP_win.exe --ked ./exp/ked/Cu_twin/Cu111_twin.lmp -e Cu
	>>./bin/AAVDP_win.exe --ked ./exp/ked/Cu_twin/Cu111_twin.lmp -e Cu -en 100
	>>./bin/AAVDP_win.exe --ked ./exp/ked/Cu_twin/Cu111_twin.lmp -e Cu -en 100 -gauss

	>>./bin/AAVDP_win.exe --ked ./exp/ked/Al2O3/Al2O3.vasp -q 0.8 -o ./exp/ked/Al2O3/Al2O3_spot.001.ked
	>>./bin/AAVDP_win.exe --ked ./exp/ked/Al2O3/Al2O3.vasp -q 0.8 -o ./exp/ked/Al2O3/Al2O3.001.ked -gauss
	>>./bin/AAVDP_win.exe --ked ./exp/ked/Al2O3/Al2O3.vasp -z 2 1 0 -q 0.8 -o ./exp/ked/Al2O3/Al2O3_spot.210.ked -rotate -rotate_x 0 0 1 -rotate_y -0.44721360 0.89442719 0
	>>./bin/AAVDP_win.exe --ked ./exp/ked/Al2O3/Al2O3.vasp -z 2 1 0 -q 0.8 -o ./exp/ked/Al2O3/Al2O3.210.ked -gauss -gauss_dx 0.00375 -gauss_sig 0.0075 -rotate -rotate_x 0 0 1 -rotate_y -0.44721360 0.89442719 0
	>>./bin/AAVDP_win.exe --ked ./exp/ked/Cu_twin/Cu111_twin.lmp -e Cu -z 0 1 0 -o ./exp/ked/Cu_twin/Cu111_twin.ked -gauss -gauss_sig 0.02 -rotate -rotate_x 0 0 1 -rotate_y 1 0 0
	>>./bin/AAVDP_win.exe --ked ./exp/ked/Fe_twin/Fe112_twin.lmp -e Fe -z 0 1 0 -o ./exp/ked/Fe_twin/Fe112_twin.ked -gauss -gauss_sig 0.02 -rotate -rotate_x 0 0 1 -rotate_y 1 0 0
4.  kkd:
	>>./bin/AAVDP_win.exe --kkd ./exp/kkd/bcc/bcc.vasp
	>>./bin/AAVDP_win.exe --kkd ./exp/kkd/bcc/bcc.vasp -en 100
	>>./bin/AAVDP_win.exe --kkd ./exp/kkd/bcc/bcc.vasp -en 100 -ked3
	>>./bin/AAVDP_win.exe --kkd ./AAVDP.ked3

	>>./bin/AAVDP_win.exe --kkd ./exp/kkd/bcc/bcc.vasp -q 1.5 -o ./exp/kkd/bcc/bcc.ked3 -ked3
	>>./bin/AAVDP_win.exe --kkd ./exp/kkd/bcc/bcc.ked3 -z 1 1 3 -rx 0.25 -ry 0.30 -t 0.3 -px 2000 -py 2400 -o ./exp/kkd/bcc/bcc.113.kkd -rotate -rotate_x 0 3 -1 -rotate_y 10 -1 -3
	>>./bin/AAVDP_win.exe --kkd ./exp/kkd/bcc/bcc.ked3 -t 0.3 -px 4000 -py 4000 -o ./exp/kkd/bcc/bcc.kkd

	>>./bin/AAVDP_win.exe --kkd ./exp/kkd/hcp/hcp.vasp -q 1.5 -t 0.3 -px 2000 -py 2000 -o ./exp/kkd/hcp/hcp.kkd
	
	>>./bin/AAVDP_win.exe --kkd ./exp/kkd/fcc/fcc.vasp -q 1.5 -o ./exp/kkd/fcc/fcc.ked3 -ked3
	>>./bin/AAVDP_win.exe --kkd ./exp/kkd/fcc/fcc.ked3 -z 1 0 0 -rx 0.075 -ry 0.075 -t 0.1 -px 1000 -py 1000 -o ./exp/kkd/fcc/fcc.100.kkd -rotate -rotate_x 0.00000000 0.47609129 -0.33664738 -rotate_y 0.00000000 0.33664738 0.47609129 -scale -scale_i 0 350
	>>./bin/AAVDP_win.exe --kkd ./exp/kkd/fcc/fcc.ked3 -z 0 1 0 -rx 0.075 -ry 0.075 -t 0.1 -px 1000 -py 1000 -o ./exp/kkd/fcc/fcc.010.kkd -rotate -rotate_x 0 0 1 -rotate_y 1 0 0 -scale -scale_i 0 130
	>>./bin/AAVDP_win.exe --kkd ./exp/kkd/fcc/fcc.ked3 -rx 0.075 -ry 0.075 -t 0.1 -px 1000 -py 1000 -o ./exp/kkd/fcc/fcc.001.kkd -scale -scale_i 0 350

	>>./bin/AAVDP_win.exe --kkd ./exp/kkd/fcc_screw/screw.vasp -q 1.5 -o ./exp/kkd/fcc_screw/screw.ked3 -ked3
	>>./bin/AAVDP_win.exe --kkd ./exp/kkd/fcc_screw/screw.ked3 -z 1 0 0 -rx 0.075 -ry 0.075 -t 0.1 -px 1000 -py 1000 -o ./exp/kkd/fcc_screw/screw.100.kkd -rotate -rotate_x 0.00000000 0.47609129 -0.33664738 -rotate_y 0.00000000 0.33664738 0.47609129 -scale -scale_i 0 350
	>>./bin/AAVDP_win.exe --kkd ./exp/kkd/fcc_screw/screw.ked3 -z 0 1 0 -rx 0.075 -ry 0.075 -t 0.1 -px 1000 -py 1000 -o ./exp/kkd/fcc_screw/screw.010.kkd -rotate -rotate_x 0 0 1 -rotate_y 1 0 0 -scale -scale_i 0 130
	>>./bin/AAVDP_win.exe --kkd ./exp/kkd/fcc_screw/screw.ked3 -rx 0.075 -ry 0.075 -t 0.1 -px 1000 -py 1000 -o ./exp/kkd/fcc_screw/screw.001.kkd -scale -scale_i 0 350

	>>./bin/AAVDP_win.exe --kkd ./exp/kkd/fcc_edge/edge.vasp -q 1.5 -o ./exp/kkd/fcc_edge/edge.ked3 -ked3
	>>./bin/AAVDP_win.exe --kkd ./exp/kkd/fcc_edge/edge.ked3 -z 1 0 0 -rx 0.075 -ry 0.075 -t 0.1 -px 1000 -py 1000 -o ./exp/kkd/fcc_edge/edge.100.kkd -rotate -rotate_x 0.00000000 0.47609129 -0.33664738 -rotate_y 0.00000000 0.33664738 0.47609129 -scale -scale_i 0 350
	>>./bin/AAVDP_win.exe --kkd ./exp/kkd/fcc_edge/edge.ked3 -z 0 1 0 -rx 0.075 -ry 0.075 -t 0.1 -px 1000 -py 1000 -o ./exp/kkd/fcc_edge/edge.010.kkd -rotate -rotate_x 0 0 1 -rotate_y 1 0 0 -scale -scale_i 0 130
	>>./bin/AAVDP_win.exe --kkd ./exp/kkd/fcc_edge/edge.ked3 -rx 0.075 -ry 0.075 -t 0.1 -px 1000 -py 1000 -o ./exp/kkd/fcc_edge/edge.001.kkd -scale -scale_i 0 350
5.  ded:
	>>./bin/AAVDP_win.exe --ded ./exp/ded/Al2O3/Al2O3.vasp
	>>./bin/AAVDP_win.exe --ded ./exp/ded/Al2O3/Al2O3.vasp -en 100
	>>./bin/AAVDP_win.exe --ded ./exp/ded/Al2O3/Al2O3.vasp -en 100 -gauss

	>>./bin/AAVDP_win.exe --ded ./exp/ded/Al2O3/Al2O3.vasp -q 0.8 -o ./exp/ded/Al2O3/Al2O3_spot.001.ded
	>>./bin/AAVDP_win.exe --ded ./exp/ded/Al2O3/Al2O3.vasp -q 0.8 -o ./exp/ded/Al2O3/Al2O3.001.ded -gauss
	>>./bin/AAVDP_win.exe --ded ./exp/ded/Al2O3/Al2O3.vasp -z 2 1 0 -fn 2 1 0 -q 0.8 -o ./exp/ded/Al2O3/Al2O3_spot.210.ded -rotate -rotate_x 0 0 1 -rotate_y -0.44721360 0.89442719 0
	>>./bin/AAVDP_win.exe --ded ./exp/ded/Al2O3/Al2O3.vasp -z 2 1 0 -fn 2 1 0 -q 0.8 -o ./exp/ded/Al2O3/Al2O3.210.ded -gauss -gauss_dx 0.00375 -gauss_sig 0.0075 -rotate -rotate_x 0 0 1 -rotate_y -0.44721360 0.89442719 0
6.  dkd:
	>>./bin/AAVDP_win.exe --dkd ./exp/dkd/bcc/bcc.vasp
	>>./bin/AAVDP_win.exe --dkd ./exp/dkd/bcc/bcc.vasp -en 100
	>>./bin/AAVDP_win.exe --dkd ./exp/dkd/bcc/bcc.vasp -monte

	>>./bin/AAVDP_win.exe --dkd ./exp/dkd/bcc/bcc.vasp -q 1.5 -z 1 1 3 -rx 0.25 -ry 0.30 -px 500 -py 600 -o ./exp/dkd/bcc/bcc.113.dkd -monte -monte_seed ./RandomSeeds.data -monte_o ./exp/dkd/bcc/bcc.mc -rotate -rotate_x 0 3 -1 -rotate_y 10 -1 -3
	>>./bin/AAVDP_win.exe --dkd ./exp/dkd/bcc/bcc.vasp -q 1.5 -px 600 -py 600 -o ./exp/dkd/bcc/bcc.dkd -monte -monte_seed ./RandomSeeds.data -monte_o ./exp/dkd/bcc/bcc.mc
	>>./bin/AAVDP_win.exe --dkd ./exp/dkd/bcc/bcc.vasp -q 1.5 -px 600 -py 600 -o ./exp/dkd/bcc/bcc.no_mc.dkd
	>>./bin/AAVDP_win.exe --dkd ./exp/dkd/hcp/hcp.vasp -q 1.5 -px 600 -py 600 -o ./exp/dkd/hcp/hcp.dkd -monte -monte_seed ./RandomSeeds.data -monte_o ./exp/dkd/hcp/hcp.mc
	>>./bin/AAVDP_win.exe --dkd ./exp/dkd/hcp/hcp.vasp -q 1.5 -px 600 -py 600 -o ./exp/dkd/hcp/hcp.no_mc.dkd
7.  rdf:
	>>./bin/AAVDP_win.exe --rdf ./exp/rdf/Cu_glass/Cu_glass.lammps

	>>./bin/AAVDP_win.exe --rdf ./exp/rdf/Cu_glass/Cu_glass.lammps -r 8 -n 160 -o ./exp/rdf/Cu_glass/Cu_glass.rdf
	>>./bin/AAVDP_win.exe --rdf ./exp/rdf/Cu_liquid/Cu_liquid.lammps -r 8 -n 160 -o ./exp/rdf/Cu_liquid/Cu_liquid.rdf
	>>./bin/AAVDP_win.exe --rdf ./exp/rdf/CuZr/Cu50Zr50.liquid.lammps -r 10 -n 100 -partial -o ./exp/rdf/CuZr/Cu50Zr50.liquid.rdf
8.  ssf:
	>>./bin/AAVDP_win.exe --ssf ./exp/ssf/Cu_glass/Cu_glass.lammps
	>>./bin/AAVDP_win.exe --ssf ./exp/ssf/Cu_glass/Cu_glass.lammps -rdf

	>>./bin/AAVDP_win.exe --ssf ./exp/ssf/Cu_glass/Cu_glass.lammps -q 8 -n 160 -o ./exp/ssf/Cu_glass/Cu_glass.ssf -rdf -rdf_r 7.5 -rdf_n 150
	>>./bin/AAVDP_win.exe --ssf ./exp/ssf/Cu_liquid/Cu_liquid.lammps -q 8 -n 160 -o ./exp/ssf/Cu_liquid/Cu_liquid.ssf -rdf -rdf_r 7.5 -rdf_n 150
	>>./bin/AAVDP_win.exe --ssf ./exp/ssf/CuZr/Cu50Zr50.liquid.lammps -q 10 -n 100 -partial -o ./exp/ssf/CuZr/Cu50Zr50.liquid.ssf
9.  other:
	>>./bin/AAVDP_win.exe -h/--help
	>>./bin/AAVDP_win.exe -v/--version

