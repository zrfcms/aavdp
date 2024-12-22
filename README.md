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
	>>./AAVDP.exe --xrd ./exp/xrd/ZnO/ZnO.vasp -2t 0 80 -o ./exp/xrd/ZnO/ZnO_line.xrd
	>>./AAVDP.exe --xrd ./exp/xrd/ZnO/ZnO.vasp -2t 0 80 -o ./exp/xrd/ZnO/ZnO.xrd --scherrer -m 0.5 -d 294 -d2t 0.02
	>>./AAVDP.exe --xrd ./exp/xrd/NaCl/NaCl.vasp -dw 1.72 1.41 -l 1.54056 -2t 20 90 -o ./exp/xrd/NaCl/NaCl1.xrd
	>>./AAVDP.exe --xrd ./exp/xrd/NaCl/NaCl.vasp -dw 1.72 1.41 -l 1.54439 -2t 20 90 -o ./exp/xrd/NaCl/NaCl2.xrd
	>>./AAVDP.exe --xrd ./exp/xrd/FeCo/FeCo.lmp -e Fe Co -dw 0.5500 0.5500 -2t 0 150 -l 1.5400 -o ./exp/xrd/FeCo/FeCo.xrd
	>>./AAVDP.exe --xrd ./exp/xrd/FeCo/FeCo_random.lmp -e Fe Co -dw 0.5500 0.5500 -2t 0 150 -l 1.5400 -c 10 10 10 -o ./exp/xrd/FeCo/FeCo_random.xrd
2.  ned:
	>>./AAVDP.exe --ned ./exp/ned/LaCrGe3/LaCrGe3.vasp -dw 0.7646 0.1028 0.2458 -2t 10 80 -l 2.43955 -o ./exp/ned/LaCrGe3/LaCrGe3_line.ned
	>>./AAVDP.exe --ned ./exp/ned/LaCrGe3/LaCrGe3.vasp -dw 0.7646 0.1028 0.2458 -2t 10 80 -l 2.43955 -o ./exp/ned/LaCrGe3/LaCrGe3.ned --scherrer -d 400 -d2t 0.02
	>>./AAVDP.exe --ned ./exp/ned/FeCo/FeCo.lmp -e Fe Co -dw 0.5500 0.5500 -2t 0 150 -l 1.5400 -o ./exp/ned/FeCo/FeCo.ned
	>>./AAVDP.exe --ned ./exp/ned/FeCo/FeCo_random.lmp -e Fe Co -dw 0.5500 0.5500 -2t 0 150 -l 1.5400 -c 10 10 10 -o ./exp/ned/FeCo/FeCo_random.ned
3.  ked:
	>>./AAVDP.exe --ked ./exp/ked/Al2O3/Al2O3.vasp -q 0.8 -o ./exp/ked/Al2O3/Al2O3_point.001.ked
	>>./AAVDP.exe --ked ./exp/ked/Al2O3/Al2O3.vasp -q 0.8 -o ./exp/ked/Al2O3/Al2O3.001.ked --gauss
	>>./AAVDP.exe --ked ./exp/ked/Al2O3/Al2O3.vasp -z 2 1 0 -q 0.8 -o ./exp/ked/Al2O3/Al2O3_point.210.ked --rotate -x 0 0 1 -y -0.44721360 0.89442719 0
	>>./AAVDP.exe --ked ./exp/ked/Al2O3/Al2O3.vasp -z 2 1 0 -q 0.8 -o ./exp/ked/Al2O3/Al2O3.210.ked --gauss -dx 0.00375 -sig 0.0075 --rotate -x 0 0 1 -y -0.44721360 0.89442719 0
	>>./AAVDP.exe --ked ./exp/ked/Cu_twin/Cu111_twin.lmp -e Cu -z 0 1 0 -o ./exp/ked/Cu_twin/Cu111_twin.ked --gauss -sig 0.02 --rotate -x 0 0 1 -y 1 0 0
	>>./AAVDP.exe --ked ./exp/ked/Fe_twin/Fe112_twin.lmp -e Fe -z 0 1 0 -o ./exp/ked/Fe_twin/Fe112_twin.ked --gauss -sig 0.02 --rotate -x 0 0 1 -y 1 0 0
4.  kkd:
	>>./AAVDP.exe --kkd ./exp/kkd/bcc/bcc.vasp -q 1.5 -o ./exp/kkd/bcc/bcc.ked3
	>>./AAVDP.exe --kkd ./exp/kkd/bcc/bcc.ked3 -x 0 3 -1 -y 10 -1 -3 -z 1 1 3 -rx 0.25 -ry 0.30 -t 0.3 -px 2000 -py 2400 -o ./exp/kkd/bcc/bcc.113.kkd
	>>./AAVDP.exe --kkd ./exp/kkd/bcc/bcc.ked3 -t 0.3 -px 4000 -py 4000 -o ./exp/kkd/bcc/bcc.kkd

	>>./AAVDP.exe --kkd ./exp/kkd/hcp/hcp.vasp -q 1.5 -o ./exp/kkd/hcp/hcp.ked3
	>>./AAVDP.exe --kkd ./exp/kkd/hcp/hcp.ked3 -t 0.3 -px 2000 -py 2000 -o ./exp/kkd/hcp/hcp.kkd
	
	>>./AAVDP.exe --kkd ./exp/kkd/fcc/fcc.vasp -q 1.5 -o ./exp/kkd/fcc/fcc.ked3
	>>./AAVDP.exe --kkd ./exp/kkd/fcc/fcc.ked3 -rx 0.075 -ry 0.075 -t 0.1 -px 1000 -py 1000 -o ./exp/kkd/fcc/fcc.001.kkd --scale -max 350
	>>./AAVDP.exe --kkd ./exp/kkd/fcc/fcc.ked3 -x 0.00000000 0.47609129 -0.33664738 -y 0.00000000 0.33664738 0.47609129 -z 1 0 0 -rx 0.075 -ry 0.075 -t 0.1 -px 1000 -py 1000 -o ./exp/kkd/fcc/fcc.100.kkd --scale -max 350
	>>./AAVDP.exe --kkd ./exp/kkd/fcc/fcc.ked3 -x 0 0 1 -y 1 0 0 -z 0 1 0 -rx 0.075 -ry 0.075 -t 0.1 -px 1000 -py 1000 -o ./exp/kkd/fcc/fcc.010.kkd --scale -max 130

	>>./AAVDP.exe --kkd ./exp/kkd/fcc_edge/edge.vasp -q 1.5 -o ./exp/kkd/fcc_edge/edge.ked3
	>>./AAVDP.exe --kkd ./exp/kkd/fcc_edge/edge.ked3 -rx 0.075 -ry 0.075 -t 0.1 -px 1000 -py 1000 -o ./exp/kkd/fcc_edge/edge.001.kkd --scale -max 350
	>>./AAVDP.exe --kkd ./exp/kkd/fcc_edge/edge.ked3 -x 0.00000000 0.47609129 -0.33664738 -y 0.00000000 0.33664738 0.47609129 -z 1 0 0 -rx 0.075 -ry 0.075 -t 0.1 -px 1000 -py 1000 -o ./exp/kkd/fcc_edge/edge.100.kkd --scale -max 350
	>>./AAVDP.exe --kkd ./exp/kkd/fcc_edge/edge.ked3 -x 0 0 1 -y 1 0 0 -z 0 1 0 -rx 0.075 -ry 0.075 -t 0.1 -px 1000 -py 1000 -o ./exp/kkd/fcc_edge/edge.010.kkd --scale -max 130

	>>./AAVDP.exe --kkd ./exp/kkd/fcc_screw/screw.vasp -q 1.5 -o ./exp/kkd/fcc_screw/screw.ked3
	>>./AAVDP.exe --kkd ./exp/kkd/fcc_screw/screw.ked3 -rx 0.075 -ry 0.075 -t 0.1 -px 1000 -py 1000 -o ./exp/kkd/fcc_screw/screw.001.kkd --scale -max 350
	>>./AAVDP.exe --kkd ./exp/kkd/fcc_screw/screw.ked3 -x 0.00000000 0.47609129 -0.33664738 -y 0.00000000 0.33664738 0.47609129 -z 1 0 0 -rx 0.075 -ry 0.075 -t 0.1 -px 1000 -py 1000 -o ./exp/kkd/fcc_screw/screw.100.kkd --scale -max 350
	>>./AAVDP.exe --kkd ./exp/kkd/fcc_screw/screw.ked3 -x 0 0 1 -y 1 0 0 -z 0 1 0 -rx 0.075 -ry 0.075 -t 0.1 -px 1000 -py 1000 -o ./exp/kkd/fcc_screw/screw.010.kkd --scale -max 130
5.  ded:
	>>./AAVDP.exe --ded ./exp/ded/Al2O3/Al2O3.vasp -q 0.8 -o ./exp/ded/Al2O3/Al2O3_point.001.ded
	>>./AAVDP.exe --ded ./exp/ded/Al2O3/Al2O3.vasp -q 0.8 -o ./exp/ded/Al2O3/Al2O3.001.ded --gauss
	>>./AAVDP.exe --ded ./exp/ded/Al2O3/Al2O3.vasp -z 2 1 0 -fn 2 1 0 -q 0.8 -o ./exp/ded/Al2O3/Al2O3_point.210.ded --rotate -x 0 0 1 -y -0.44721360 0.89442719 0
	>>./AAVDP.exe --ded ./exp/ded/Al2O3/Al2O3.vasp -z 2 1 0 -fn 2 1 0 -q 0.8 -o ./exp/ded/Al2O3/Al2O3.210.ded --gauss -dx 0.00375 -sig 0.0075 --rotate -x 0 0 1 -y -0.44721360 0.89442719 0
6.  dkd: 
	>>./AAVDP.exe --dkd ./exp/dkd/bcc/bcc.vasp -q 1.5 -x 0 3 -1 -y 10 -1 -3 -z 1 1 3 -rx 0.25 -ry 0.30 -px 500 -py 600 -o ./exp/dkd/bcc/bcc.113.dkd --monte -o ./exp/dkd/bcc/bcc.mc
	>>./AAVDP.exe --dkd ./exp/dkd/bcc/bcc.vasp -q 1.5 -px 600 -py 600 -o ./exp/dkd/bcc/bcc.dkd --monte -o ./exp/dkd/bcc/bcc.mc
	>>./AAVDP.exe --dkd ./exp/dkd/hcp/hcp.vasp -q 1.5 -px 600 -py 600 -o ./exp/dkd/hcp/hcp.dkd --monte -o ./exp/dkd/hcp/hcp.mc
	>>./AAVDP.exe --dkd ./exp/dkd/hcp/hcp.vasp -q 1.5 -px 600 -py 600 -o ./exp/dkd/hcp/hcp.no_mc.dkd
7.  rdf:
	>>./AAVDP.exe --rdf ./exp/rdf/Cu_glass/Cu_glass.lammps -r 8 -n 160 -o ./exp/rdf/Cu_glass/Cu_glass.rdf
	>>./AAVDP.exe --rdf ./exp/rdf/Cu_liquid/Cu_liquid.lammps -r 8 -n 160 -o ./exp/rdf/Cu_liquid/Cu_liquid.rdf
	>>./AAVDP.exe --rdf ./exp/rdf/CuZr/Cu50Zr50.liquid.lammps -r 10 -n 100 -partial -o ./exp/rdf/CuZr/Cu50Zr50.liquid.rdf
8.  ssf:
	>>./AAVDP.exe --ssf ./exp/ssf/Cu_glass/Cu_glass.lammps -q 8 -n 160 -o ./exp/ssf/Cu_glass/Cu_glass.ssf --rdf -r 7.5 -n 150
	>>./AAVDP.exe --ssf ./exp/ssf/Cu_liquid/Cu_liquid.lammps -q 8 -n 160 -o ./exp/ssf/Cu_liquid/Cu_liquid.ssf --rdf -r 7.5 -n 150
	>>./AAVDP.exe --ssf ./exp/ssf/CuZr/Cu50Zr50.liquid.lammps -q 10 -n 100 -partial -o ./exp/ssf/CuZr/Cu50Zr50.liquid.ssf
9.  other:
	>>./AAVDP.exe -h/--help
	>>./AAVDP.exe -h xrd
	>>./AAVDP.exe -v/--version

