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

def check_nan_in_positions(label, array):
    num_non_finite = np.sum(~np.isfinite(array))
    min_val = np.min(array)
    max_val = np.max(array)
    print("Position", label)
    print("\t min/max: {}/{}".format(min_val, max_val))
    print("\t non-finite vals: {}, {}%".format(num_non_finite, 100.0*num_non_finite/len(array)))

def check_positions(file1):
    ''' Check that the positions of the cores look sensible'''
    try:
        x = gio.gio_read(file1,'x')
        y = gio.gio_read(file1,'y')
        z = gio.gio_read(file1,'z')
        check_nan_in_position('x', x)
        check_nan_in_position('y', y)
        check_nan_in_position('z', z)
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

def check_central(file1,step):
    ''' Check that all halos have a central core. Note: Check is currently only for 499 - extend this '''
    try:
        halo_tag = gio.gio_read(file1,'fof_halo_tag')
        infall_step = gio.gio_read(file1,'infall_step')
        unique_tags, tag_counts = np.unique(halo_tag,return_counts=True)
        print("maximum number of cores in a halo is " +str(max(tag_counts)))
        print("number of halos is " + str(len(unique_tags)))
        print("number of centrals at "+str(step)+" is " + str(np.sum(infall_step==step)))
        unique_tags_central, tag_counts_central = np.unique(halo_tag[infall_step==step],return_counts=True)
        print("All halos have a central = " + str((np.sort(unique_tags)==np.sort(unique_tags_central)).all()))
    except:
        print("test failed")


def check_mt_files(file1,file_mt):
    ''' look at how many fof halos we have above threshold '''
    try:
        halo_tag = gio.gio_read(file_mt,'fof_halo_tag')
        halo_mass = gio.gio_read(file_mt,'fof_halo_mass')
        halo_tag_c = gio.gio_read(file1,'fof_halo_tag')
        infall_step = gio.gio_read(file1,'infall_step')
        halo_tag_c=halo_tag_c[infall_step==499]
        print(len(halo_tag),len(halo_tag_c))
    except:
        print("test_failed")
 
# positions of cores that exist at one timestep and not the next
=======
    

def check_unique_core_tags(file1):
    """
    """
    try:
        core_tag = gio.gio_read(file1, "core_tag")
        infall_tree_node_index = gio.gio_read(file1, "infall_tree_node_index")
        print("Unique core_tags: ", is_unique_array(core_tag))
        print("Unique infall_tree_node_index: ", is_unique_array(infall_tree_node_index))
    except:
        print("test failed")

def is_unique_array(array):
    """Returns true is all the values in the array are unique. False if
    there are repeating values.

    """
    u = np.unique(array)
    return len(u) == len(array)


if __name__ == "__main__":
    file1 = sys.argv[1]
    check_infall_masses(file1)
    check_infall_timestep(file1)
    check_positions(file1)
    check_central(file1)

