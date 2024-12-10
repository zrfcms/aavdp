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

1.  kkd:
	>>./AAVDP.exe --kkd ./exp/bcc_kkd/bcc.vasp -q 1.5 -o ./exp/bcc_kkd/bcc.ked3
	>>./AAVDP.exe --kkd ./exp/bcc_kkd/bcc.ked3 -x 0 3 -1 -y 10 -1 -3 -z 1 1 3 -rx 0.25 -ry 0.30 -px 1000 -py 1200 -o ./exp/bcc_kkd/bcc.113.kkd
	>>./AAVDP.exe --kkd ./exp/hcp_kkd/hcp.vasp -q 1.5 -o ./exp/hcp_kkd/hcp.ked3
	>>./AAVDP.exe --kkd ./exp/hcp_kkd/hcp.ked3 -t 0.1 -px 3000 -py 3000 -ortho -o ./exp/hcp_kkd/hcp.001.kkd
2.  dkd:
	>>./AAVDP.exe --dkd ./exp/bcc_dkd/bcc.vasp -x 0 3 -1 -y 10 -1 -3 -z 1 1 3 -rx 0.25 -ry 0.30 -px 1000 -py 1200 -o ./exp/bcc_dkd/bcc.113.dkd
	>>./AAVDP.exe --dkd ./exp/bcc_dkd/bcc.vasp -rx 0.1 -ry 0.1 -px 500 -py 500 -o ./exp/bcc_dkd/bcc.001.dkd
	>>./AAVDP.exe --kkd ./exp/hcp_kkd/hcp.vasp -q 1.5 -o ./exp/hcp_kkd/hcp.ked3
	>>./AAVDP.exe --kkd ./exp/hcp_kkd/hcp.ked3 -t 0.1 -px 3000 -py 3000 -ortho -o ./exp/hcp_kkd/hcp.001.kkd

