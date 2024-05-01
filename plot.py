#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from pathlib import Path
import argparse
import numpy as np
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

rcParams['xtick.direction'] = 'in' 
rcParams['ytick.direction'] = 'in' 
rcParams['xtick.top'] = 'True'
rcParams['ytick.right'] = 'True'
rcParams['figure.frameon'] = 'False'
plt.rcParams['font.sans-serif'] = ['Times New Roman']
plt.rcParams['font.size'] = 18
plt.rcParams['axes.unicode_minus'] = False
class AAVDPPlot(object):

    def __init__(self, xlabel=r'2$\theta$', ylabel='Intensity (A.U.)', title='Simulated XRD pattern'):
        _, curve = plt.subplots(1, 1, figsize = (8, 5), dpi = 120, constrained_layout = True)
        curve.set_xlabel(xlabel)
        curve.set_ylabel(ylabel)
        curve.set_title(title)
        curve.tick_params(which = 'both', width = 1)
        curve.tick_params(which = 'major', length = 6)
        curve.tick_params(which = 'minor', length = 4)
        for dir in ['left', 'right', 'top', 'bottom']:
            curve.spines[dir].set_linewidth(1)
        self.curve = curve

    def cxrd(self, input, xlim = ()):
        def fmt(x):
            if x<0: return r'$\overline{' + str(abs(x)) + '}$'
            else: return str(x)
        data = np.loadtxt(input, dtype = np.float32, skiprows = 1)
        x, y = data[:, 3], data[:, 5]
        bar = self.curve.bar(x, y, width = 0.5, linewidth = 0.0)

        hkl = data[:, :3].astype('int')
        annotations = ['{0} {1} {2}'.format(fmt(hkl[i,0]), fmt(hkl[i,1]), fmt(hkl[i,2])) for i in range(hkl.shape[0])]
        self.curve.bar_label(bar, annotations,
                             padding=10, color='black', fontweight='bold', rotation='vertical', fontsize=8)
        
        self.curve.set_ylim(-5, 105)
        if xlim:
            self.curve.set_xlim(xlim)
        else:
            self.curve.set_xlim(np.floor(min(x)/10)*10,np.ceil(max(x)/10)*10)
        yminor, ymajor = 10, 20
        xminor, xmajor = 10, 20
        self.curve.yaxis.set_major_locator(MultipleLocator(ymajor))  
        self.curve.yaxis.set_minor_locator(MultipleLocator(yminor))
        self.curve.xaxis.set_major_locator(MultipleLocator(xmajor))
        self.curve.xaxis.set_minor_locator(MultipleLocator(xminor))
        plt.axhline(0, color='red')

    def rxrd(self, input, binlim = (10.0, 100.0), binnum = 250, xlim = ()):
        data = np.loadtxt(input, dtype = np.float32, skiprows = 3, usecols=[1,2])
        x, y = data[:, 0], data[:, 1]
        if binnum:
            avex, avey = np.zeros(binnum), np.zeros(binnum)
            binlen = (binlim[1]-binlim[0])/binnum
            for i in range(binnum):
                ilim = (binlim[0]+binlen*i, binlim[0]+binlen*(i+1))
                selectx = (x>=ilim[0]) & (x<=ilim[1])
                limy = y[selectx]
                print(limy)
                avex[i] = (ilim[0]+ilim[1])/2.0
                if limy.size>0: avey[i] = limy.sum()
                else: avey[i] = 0.0
            avey = avey[np.argsort(avex)]
            avey = avey/np.max(avey)*100.0
            avex.sort()
        else:
            avex, avey = data[:, 1], data[:, 2]
        self.curve.plot(avex, avey, linewidth = 1)
        # for i in range(binnum):
        #     print(avex[i], avey[i])
        
        self.curve.set_ylim(-5, 105)
        if xlim:
            self.curve.set_xlim(xlim)
        else:
            self.curve.set_xlim(np.floor(min(avex)/10)*10,np.ceil(max(avex)/10)*10)
        yminor, ymajor = 10, 20
        xminor, xmajor = 10, 20
        self.curve.yaxis.set_major_locator(MultipleLocator(ymajor))
        self.curve.yaxis.set_minor_locator(MultipleLocator(yminor))
        self.curve.xaxis.set_major_locator(MultipleLocator(xmajor))
        self.curve.xaxis.set_minor_locator(MultipleLocator(xminor))
        plt.axhline(0, color='red')

    def save(self, path):
        plt.savefig(path, bbox_inches = 'tight', dpi = 300)
        plt.close()

class Plot(object):
    def __init__(self, xlabel='r', ylabel='g(r)'):
        _, curve = plt.subplots(1, 1, figsize=(8, 5), dpi=300, constrained_layout=True)
        curve.set_xlabel(xlabel)
        curve.set_ylabel(ylabel)
        curve.tick_params(which='both', width=1)
        curve.tick_params(which='major', length=6)
        curve.tick_params(which='minor', length=4)
        for dir in ['left', 'right', 'top', 'bottom']:
            curve.spines[dir].set_linewidth(1)
        self.curve=curve
        self.n=0
        self.x, self.y=[], []

    def load(self, txt_path, skiprows=0, max_rows=400, usecols_x=(0, ), usecols_y=(1, )):
        self.x=np.loadtxt(txt_path, skiprows=skiprows, max_rows=max_rows, usecols=usecols_x, comments='#')
        self.y=np.loadtxt(txt_path, skiprows=skiprows, max_rows=max_rows, usecols=usecols_y, comments='#')
        nx, ny=len(usecols_x), len(usecols_y)
        self.n=max(nx, ny)
        if(nx==1 and ny>1):
            self.x=np.tile(self.x[:, 0], (1, ny))
        if(ny==1 and nx>1):
            self.y=np.tile(self.y[:, 0], (1, nx))

    def rdf(self, xlim=(), ylim=()):
        if(self.n==1):
            self.curve.plot(self.x, self.y, linewidth=1)
        else:
            for i in range(self.n):
                self.curve.plot(self.x[:, i], self.y[:, i], linewidth=1)
        if xlim:
            self.curve.set_xlim(xlim)
        else:
            self.curve.set_xlim(np.floor(min(self.x)/10)*10,np.ceil(max(self.x)/10)*10)
        if ylim:
            self.curve.set_ylim(ylim)
        else:
            self.curve.set_ylim(np.floor(min(self.y)/10)*10,np.ceil(max(self.y)/10)*10)
        xmajor, ymajor=1, 10
        xminor, yminor=xmajor/5, ymajor/5
        self.curve.yaxis.set_major_locator(MultipleLocator(ymajor))
        self.curve.yaxis.set_minor_locator(MultipleLocator(yminor))
        self.curve.xaxis.set_major_locator(MultipleLocator(xmajor))
        self.curve.xaxis.set_minor_locator(MultipleLocator(xminor))

    def save(self, path):
        plt.savefig(path, bbox_inches = 'tight', dpi = 300)
        plt.close()

if __name__ == "__main__":
    # parser = argparse.ArgumentParser(description='Check an interatomic potential')
    # parser.add_argument('pairs', nargs='+', help='A series of files required for an interatomic potential user input')
    # parser.add_argument('-regular', action='store_true', default=False, dest='regularFlag', help='Determine whether the input is regular expression.')
    # parser.add_argument('-c', required=True, dest='combo', nargs='+', help='Items list for checking, user-defined or default')
    # parser.add_argument('-e', required=True, dest='eles', nargs='+', help='A series of elements required for checking')
    # parser.add_argument('-o', default='.', dest='exportDir', help='Directory for saving checking results saved in json format, i.e., .chk.json, accompanying temporary folder')  
    # parser.add_argument('-t', default='', dest='pairType', help='Type of an interatomic potential user input, e.g. eam/alloy, eam/fs, meam, bop, tersoff, sw') 
    # parser.add_argument('-aniso', default=0, dest='aniso', type=float, help='Aniso pressure (in GPa)') 
    # parser.add_argument('-dump', default=0, dest='dumpStep', type=int, help='Steps for dumping simulation trajectories, i.e., 0/1, among which the default is 0') 
    # parser.add_argument('-lammps', default=_lammpsPath, dest='lammpsPath', help='Path of lammps.exe') 
    # parser.add_argument('-struct', default=_structDir, dest='structDir', help='Path of struct folder that includes the default structures')
    # parser.add_argument('-log', action='store_true', default=False, dest='log', help='Save the log file, i.e., log.lammps, specified for LAMMPS to write status information to') 
    # parser.add_argument('-clear', action='store_true', default=False, dest='clear', help='Clear the temporary folder that has existed')
    # parser.add_argument('-v', action='version', version='%(prog)s 1.0')
    # args = parser.parse_args()
    # p=Path("./examples/B4C_ssf/B4C_mp696746.ssf")
    # a=Plot()
    # a.load(p)
    # a.rdf()
    # a.save(p.with_suffix('.png'))
    p = Path(".\\examples\\CdS_xrd\\CdS.xrd")
    a = AAVDPPlot(title = p.stem)
    a.cxrd(p, xlim = (20, 180))
    a.save(p.with_suffix('.jpg'))