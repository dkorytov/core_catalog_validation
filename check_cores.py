#!/bin/bash/python

import numpy as np
import sys
sys.path.append('/home/prlarsen/usr/genericio/python')
import genericio as gio
import matplotlib.pyplot as plt

file1 = sys.argv[1]

def check_infall_masses(file1):
    try:
        infall_mass = gio.gio_read(file1,'infall_mass')
        criteria  = (np.min(infall_mass)>0)
        if criteria:
            print("minimum infall mass has positive value: "+str(np.min(infall_mass)/1.09086e9)+" particles")
        else:
            print("minimum infall mass not sensible: "+str(np.min(infall_mass)))
    except:
        print("can't find file")
    try:
        plt.figure()
        plt.hist(np.log10(infall_mass),bins=100)
        plt.show()
    except:
        print("can't make plot")
     

def check_infall_timesteps(file1):
    try:
        infall_step = gio.gio_read(file1,'infall_step')
        criteria  = (np.min(infall_step)>0)
        if criteria:
            print("minimum infall step has positive value: "+str(np.min(infall_step)))
        else:
            print("minimum infall step not sensible: "+str(np.min(infall_step)))
    except:
        print("can't find file")
    try:
        plt.figure()
        plt.hist(infall_step,bins=100)
        plt.show()
    except:
        print("can't make plot")
    try:
       print(np.unique(infall_step))
    
def check_positions(file1):
    try:
        x = gio.gio_read(file1,'x')
        y = gio.gio_read(file1,'y')
        z = gio.gio_read(file1,'z')
    except:
        print("Can't read file")
    try:
        plt.figure()
        plt.hist(x,bins=100,alpha=0.4)
        plt.hist(y,bins=100,alpha=0.4)
        plt.hist(z,bins=100,alpha=0.4)
        plt.show()
    except:
        print("can't make plot")
    try:
        plt.figure()
        plt.hist2d(x,y,bins=100)
        plt.show()
    except:
        print("can't make plot")

def check_central(file1):
    try:
        halo_tag = gio.gio_read(file1,'fof_halo_tag')

