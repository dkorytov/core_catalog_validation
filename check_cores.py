#!/bin/bash/python

import numpy as np
import sys
sys.path.append('/home/prlarsen/usr/genericio/python')
import genericio as gio
import matplotlib.pyplot as plt

file1 = sys.argv[1]

def check_infall_masses(file1):
    ''' Check how many particles make up the minimum infall mass, and plot the distribution '''
    try:
        infall_mass = gio.gio_read(file1,'infall_mass')
        criteria  = (np.min(infall_mass)>0)
        if criteria:
            print("minimum infall mass has positive value: "+str(np.min(infall_mass)/1.09086e9)+" particles")
        else:
            print("minimum infall mass not sensible: "+str(np.min(infall_mass)))
        plt.figure()
        plt.hist(np.log10(infall_mass),bins=100)
        plt.show()
    except:
        print("can't run test")
     

def check_infall_timesteps(file1):
    ''' Check the distribution of infall timesteps'''
    try:
        infall_step = gio.gio_read(file1,'infall_step')
        criteria  = (np.min(infall_step)>0)
        if criteria:
            print("minimum infall step has positive value: "+str(np.min(infall_step)))
        else:
            print("minimum infall step not sensible: "+str(np.min(infall_step)))
        plt.figure()
        plt.hist(infall_step,bins=100)
        plt.show()
        print(np.unique(infall_step))
    except:
        print("can't run test")

    
def check_positions(file1):
    ''' Check that the positions of the cores look sensible'''
    try:
        x = gio.gio_read(file1,'x')
        y = gio.gio_read(file1,'y')
        z = gio.gio_read(file1,'z')
        plt.figure()
        plt.hist(x,bins=100,alpha=0.4)
        plt.hist(y,bins=100,alpha=0.4)
        plt.hist(z,bins=100,alpha=0.4)
        plt.show()
        plt.figure()
        plt.hist2d(x,y,bins=100)
        plt.show()
    except:
        print("can't run test")

def check_central(file1):
    ''' Check that all halos have a central core. Note: Check is currently only for 499 - extend this '''
    try:
        halo_tag = gio.gio_read(file1,'fof_halo_tag')
        infall_step = gio.gio_read(file1,'infall_step')
        unique_tags, tag_counts = np.unique(halo_tag,return_counts=True)
        print("maximum number of cores in a halo is " +str(max(tag_counts)))
        print("number of halos is " + str(len(unique_tags)))
        print("number of centrals at 499 is " + str(np.sum(infall_step==499)))
        unique_tags_central, tag_counts_central = np.unique(halo_tag[infall_step==499],return_counts=True)
        print("All halos have a central = " + str((np.sort(unique_tags)==np.sort(unique_tags_central)).all()))
    except:
        print("test failed")
    
