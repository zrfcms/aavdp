# aavdp

#### 介绍
Automatic Analysis of Virtual Diffraction Pattern （AAVDP）for artificial atomistic structures

#### 软件架构
1.  X-ray diffraction (XRD)
2.  Neutron diffraction (NED)
3.  Kinematical electron diffraction (KED)
4.  Kinematical Kikuchi diffraction (KKD)
5.  Dynamatical Kikuchi diffraction (DKD)
6.  Radial distribution function (RDF)
7.  Static structure factor (SSF)

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

1.  ./bin/AAVDP ./exp/dkd/Fe/Fe.lmp -e Fe -dw 0.0032


