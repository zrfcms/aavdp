# aavdp

#### 介绍
Automatic Analysis of Virtual Diffraction Pattern （AAVDP）for artificial atomistic structures

#### 软件架构
1.  X-ray diffraction (XRD)
2.  Neutron diffraction (NED)
3.  Kinematical electron diffraction (KED)
4.  Kinematical Kikuchi diffraction (KKD)
5.  Dynamatical electron diffraction (DED)
6.  Dynamatical Kikuchi diffraction (DKD)
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
	>>./AAVDP.exe --kkd ./exp/kkd/bcc/bcc.ked3 -x 0 3 -1 -y 10 -1 -3 -z 1 1 3 -rx 0.25 -ry 0.30 -px 1000 -py 1200 -o ./exp/kkd/bcc/bcc.113.kkd
5.  ded:
	>>./AAVDP.exe --ded ./exp/ded/Al2O3/Al2O3.vasp -q 0.8 -o ./exp/ded/Al2O3/Al2O3_point.001.ded
	>>./AAVDP.exe --ded ./exp/ded/Al2O3/Al2O3.vasp -q 0.8 -o ./exp/ded/Al2O3/Al2O3.001.ded --gauss
	>>./AAVDP.exe --ded ./exp/ded/Al2O3/Al2O3.vasp -z 2 1 0 -fn 2 1 0 -q 0.8 -o ./exp/ded/Al2O3/Al2O3_point.210.ded --rotate -x 0 0 1 -y -0.44721360 0.89442719 0
	>>./AAVDP.exe --ded ./exp/ded/Al2O3/Al2O3.vasp -z 2 1 0 -fn 2 1 0 -q 0.8 -o ./exp/ded/Al2O3/Al2O3.210.ded --gauss -dx 0.00375 -sig 0.0075 --rotate -x 0 0 1 -y -0.44721360 0.89442719 0
6.  dkd:
	>>./AAVDP.exe --dkd ./exp/dkd/bcc/bcc.vasp -o ./exp/dkd/bcc/bcc.dkd --monte
	
	>>./AAVDP.exe --dkd ./exp/bcc_dkd/bcc.vasp -x 0 3 -1 -y 10 -1 -3 -z 1 1 3 -rx 0.25 -ry 0.30 -px 1000 -py 1200 -o ./exp/bcc_dkd/bcc.113.dkd
	>>./AAVDP.exe --dkd ./exp/bcc_dkd/bcc.vasp -rx 0.1 -ry 0.1 -px 500 -py 500 -o ./exp/bcc_dkd/bcc.001.dkd
	>>./AAVDP.exe --kkd ./exp/hcp_kkd/hcp.vasp -q 1.5 -o ./exp/hcp_kkd/hcp.ked3
	>>./AAVDP.exe --kkd ./exp/hcp_kkd/hcp.ked3 -t 0.1 -px 3000 -py 3000 -ortho -o ./exp/hcp_kkd/hcp.001.kkd
7.  rdf:
	>>./AAVDP.exe --rdf ./exp/rdf/Cu_glass/Cu_glass.lammps -r 8 -n 160 -o ./exp/rdf/Cu_glass/Cu_glass.rdf
	>>./AAVDP.exe --rdf ./exp/rdf/Cu_liquid/Cu_liquid.lammps -r 8 -n 160 -o ./exp/rdf/Cu_liquid/Cu_liquid.rdf
	>>./AAVDP.exe --rdf ./exp/rdf/CuZr/Cu50Zr50.liquid.lammps -r 5.0 -n 100 -partial -o ./exp/rdf/CuZr/Cu50Zr50.liquid.rdf
8.  ssf:
	>>./AAVDP.exe --ssf ./exp/ssf/Cu_glass/Cu_glass.lammps -q 8 -n 160 -o ./exp/ssf/Cu_glass/Cu_glass.ssf --rdf -r 7.5 -n 150
	>>./AAVDP.exe --ssf ./exp/ssf/Cu_liquid/Cu_liquid.lammps -q 8 -n 160 -o ./exp/ssf/Cu_liquid/Cu_liquid.ssf --rdf -r 7.5 -n 150
	>>./AAVDP.exe --ssf ./exp/ssf/CuZr/Cu50Zr50.liquid.lammps -q 5.0 -n 100 -partial -o ./exp/ssf/CuZr/Cu50Zr50.liquid.ssf

