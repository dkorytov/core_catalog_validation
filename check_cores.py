#!/bin/bash/python

import numpy as np
import os
import sys
import matplotlib.pyplot as plt

cooley=False
if cooley:
    sys.path.append('/home/prlarsen/usr/genericio/python')
else:
    sys.path.append('/home/prlarsen/genericio2/genericio/python')
import genericio as gio

# global variable if you want visualisations
vis = True

def check_infall_masses(file1):
    ''' Check how many particles make up the minimum infall mass, and plot the distribution '''
    try:
        infall_mass = gio.gio_read(file1,'infall_mass')
        criteria  = (np.min(infall_mass)>0)
        if criteria:
            print("Minimum infall mass has positive value: "+str(round(np.min(infall_mass)/1.09086e9))+" particles")
            if vis:
                plt.figure()
                plt.hist(np.log10(infall_mass),bins=100)
                plt.xlabel('log10(infall mass)')
                plt.show()
            return 0
        else:
            print("minimum infall mass not sensible: "+str(np.min(infall_mass)))
            return 1
    except:
        print("can't run infall mass test")
        return 1 

     
def check_infall_timesteps(file1):
    ''' Check the distribution of infall timesteps'''
    try:
        infall_step = gio.gio_read(file1,'infall_step')
        criteria  = (np.min(infall_step)>0)
        print("List of infall steps:")
        print(np.unique(infall_step))
        if criteria:
            print("Minimum infall step has positive value: "+str(np.min(infall_step)))
            if vis:
                plt.figure()
                plt.hist(infall_step,bins=100)
                plt.xlabel('infall step')
                plt.show()
            return 0
        else:
            print("Minimum infall step not sensible: "+str(np.min(infall_step)))
            return 1
    except:
        print("can't run infall timestep test")
        return 1

def check_nan_in_positions(label, array):
    num_non_finite = np.sum(~np.isfinite(array))
    criteria = (num_non_finite==0)
    min_val = np.min(array)
    max_val = np.max(array)
    print("Position", label)
    print("\t min/max: {}/{}".format(min_val, max_val))
    print("\t non-finite vals: {}, {}%".format(num_non_finite, 100.0*num_non_finite/len(array)))
    if criteria:
        return 0 
    else:
        return 1

def check_positions(file1):
    ''' Check that the positions of the cores look sensible'''
    try:
        test = 0
        x = gio.gio_read(file1,'x').flatten()
        y = gio.gio_read(file1,'y').flatten()
        z = gio.gio_read(file1,'z').flatten()
        test += check_nan_in_positions('x', x)
        test += check_nan_in_positions('y', y)
        test += check_nan_in_positions('z', z)
        if vis:
            plt.figure()
            plt.hist(x,bins=100,alpha=0.4)
            plt.hist(y,bins=100,alpha=0.4)
            plt.hist(z,bins=100,alpha=0.4)
            plt.xlabel('position (Mpc)')
            plt.figure()
            plt.hist2d(x,y,bins=100)
            plt.xlabel('x (Mpc)')
            plt.ylabel('y (Mpc)')
            plt.show()
        return test
    except:
        print("can't run position test")

def check_central(file1,step):
    ''' Check that all halos have a central core. Note: Check is currently only for 499 - extend this '''
    try:
        halo_tag = gio.gio_read(file1,'fof_halo_tag')
        infall_step = gio.gio_read(file1,'infall_step')
        unique_tags, tag_counts = np.unique(halo_tag,return_counts=True)
        print("Maximum number of cores in a halo is " +str(max(tag_counts)))
        print("Number of halos is " + str(len(unique_tags)))
        unique_tags_central, tag_counts_central = np.unique(halo_tag[infall_step==step],return_counts=True)
        criteria = (np.sort(unique_tags)==np.sort(unique_tags_central)).all()
        if criteria:
            print("All halos have a central core")
            return 0
        else:
            print("Some halos are missing a central core")
            return 1
    except:
        print("test failed")


def check_mt_files(file1,file_mt,step):
    ''' look at how many fof halos we have above threshold '''
    try:
        plot_missing = False # useful to turn on if step = 499, so you can see if all missing objects are at mass resolution
        if (step!=499):
            print("Warning, merger tree files should only match properly at step 499 (at higher redshifts N_mergertree > N_centralcores due to fragmentation). As you are looking at a higher redshift, expect a ~5% level of missing cores with a range of masses (including 0 mass fragments).")
        halo_tag = gio.gio_read(file_mt,'fof_halo_tag')
        halo_mass = gio.gio_read(file_mt,'fof_halo_mass')
        halo_tag_c = gio.gio_read(file1,'fof_halo_tag')
        infall_step = gio.gio_read(file1,'infall_step')
        halo_tag_c = halo_tag_c[infall_step==step]
        print("Number of halos in merger tree file = "+str(len(halo_tag)))
        print("Number of central cores in core file = "+str(len(halo_tag_c)))
        diff_tot = (len(halo_tag)-len(halo_tag_c)+0.0)/len(halo_tag_c)*100.
        criteria = (diff_tot<10)
        if vis&plot_missing:
            print("Plotting missing cores - this may take a minute.")
            import random
            missing_cores = np.setdiff1d(halo_tag, halo_tag_c) 
            print("Number of missing objects = "+str(len(missing_cores)))
            random.shuffle(missing_cores)
            masses = []
            for i in range(10000):
                idx = np.where(halo_tag==missing_cores[i])[0][0]
                masses.append(halo_mass[idx]); 
            masses = np.array(masses)
            if 0 in masses:
                print("Fraction of zero-mass objects = "+str(np.sum(masses==0)/(len(masses)+0.0)))
            plt.figure()
            plt.hist(np.log10(masses[masses>0]),bins=100)
            plt.xlabel('non-zero masses (log10)')
            plt.show()
        if criteria:
            return 0 
        else: 
            return 1
    except:
        print("test_failed")
 
def check_core_evolution(file1,file2):
    ''' '''
    try:
        plot_missing = False
        core_tag1 = gio.gio_read(file1, "core_tag")
        core_tag2 = gio.gio_read(file2, "core_tag")
        # file1 should be later in time (higher step number) than file2
        missing_cores = np.setdiff1d(core_tag2,core_tag1,assume_unique=False)
        perc_missing = float(len(missing_cores))/len(core_tag2)*100
        criteria = (perc_missing<5)
        print("Percentage of cores present in earlier timestep missing in subsequent timestep = "+str(perc_missing))
        if vis&plot_missing:
            import random
            random.shuffle(missing_cores)
            print("Plotting positions of 10,000 missing cores - this may take a minute")
            x = gio.gio_read(file2,'x').flatten()
            y = gio.gio_read(file2,'y').flatten()
            z = gio.gio_read(file2,'z').flatten()
            xx = []
            yy = []
            zz = []
            for i in range(10000):
                idx = np.where(core_tag2==missing_cores[i])[0][0]
                xx.append(x[idx]); yy.append(y[idx]); zz.append(z[idx])
            xx = np.array(xx); yy = np.array(yy); zz = np.array(zz)
            plt.hist2d(xx,yy,bins=100)
            plt.title('location of missing cores')
            plt.figure()
            plt.hist2d(x,y,bins=100)
            plt.title('location of all cores - for comparison')
            plt.show()
        if criteria:
            return 0
        else:
            return 1
    except:
        print("test_failed")   

def check_unique_core_tags(file1):
    """
    """
    try:
        core_tag = gio.gio_read(file1, "core_tag")
        infall_tree_node_index = gio.gio_read(file1, "infall_tree_node_index")
        core_tag_unique = is_unique_array(core_tag)
        tree_node_unique = is_unique_array(infall_tree_node_index)
        if core_tag_unique:
            print("Core tag is unique")
        else:
            print("Core tag is not unique")
        if tree_node_unique:
            print("Infall tree node index is unique")
        else: 
            print("Infall tree node index is not unique")
        if core_tag_unique&tree_node_unique:
            return 0
        else: 
            return 1
    except:
        print("test failed")

def is_unique_array(array):
    """Returns true is all the values in the array are unique. False if
    there are repeating values.

    """
    u = np.unique(array)
    return len(u) == len(array)


if __name__ == "__main__":

    check_mt = True
    check_ev = True
    base_name = '/media/luna1/prlarsen/core_stuff/09_03_2019.AQ.' 
    base_name_mt = '/media/luna1/prlarsen/core_stuff/output/09_03_2019.AQ.'
    step = sys.argv[1]
    step2 = sys.argv[2] # note step2 should be earlier in time than step1
    file1 = base_name+step+'.coreproperties'
    file_mt = base_name_mt +step+'.treenodes'
    file_ev = base_name + step2 + '.coreproperties'

    test = 0
    test += check_infall_masses(file1)
    test += check_infall_timesteps(file1)
    test += check_positions(file1)
    test += check_central(file1,int(step))
    test += check_unique_core_tags(file1)
    if check_mt:
        test += check_mt_files(file1,file_mt,int(step))
    if check_ev:
        test += check_core_evolution(file1,file_ev)
    print("Total number of issues found = "+str(test))

