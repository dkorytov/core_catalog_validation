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
    ''' Check how many particles make up the minimum infall mass, and plot the distribution 
        criteria: minimum infall mass > 0 (this is effectively a sanity check on mass assignment)
        inputs: file1 = core catalog file 
        outputs: criteria (pass=0)
        plots: histogram of log10(mass)
    '''
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
    ''' Check the distribution of infall steps
        criteria: minimum infall step > 0 (sanity check on time step assignment)
        inputs: file1 = core catalog file 
        outputs: criteria (pass=0)
        plots: histogram of infall time step
        display: list of unique infall time steps
    '''
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

def _check_nan_in_positions(label, array):
    ''' Check for non-finite values in position arrays
        criteria: only finite values exist
        inputs: label - x,y,z label
                array - position array
        outputs: criteria (pass=0)
        display: minimum/maximum values and percentage of non-finite values
    '''
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
    ''' Check that the positions of the cores look sensible
        criteria: _check_nan_in_positions test passess for all position variables 
        inputs: file1 = core catalog file
        outputs: criteria (pass=0)
        plot: 1d histogram of x,y,z values, 2d histogram of (x,y) values
        display: output of _check_nan_in_positions test
    '''
    try:
        test = 0
        x = gio.gio_read(file1,'x').flatten()
        y = gio.gio_read(file1,'y').flatten()
        z = gio.gio_read(file1,'z').flatten()
        test += _check_nan_in_positions('x', x)
        test += _check_nan_in_positions('y', y)
        test += _check_nan_in_positions('z', z)
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
    ''' Check that all halos have a central core
        criteria: a core with infall step = current step exists for all fof halos
        inputs: file1 = core catalog file
                step = time step associated with core catalog file 
        outputs: criteria (pass=0)
        display: Maximum number of cores in a halo
                 Total number of fof halos in file
    '''
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
        
def _is_unique_array(array):
    ''' Check array for uniqueness
        outputs: True if array is unique, False otherwise
    '''
    u = np.unique(array)
    return len(u) == len(array)

        
def check_unique_core_tags(file1):
    ''' Check that core tags and infall tree nodes are unique
        criteria: no repeating core tags or infall tree nodes
        inputs: file1 = core catalog file
        outputs: criteria (pass=0)
        display: Uniqueness of core tag and infall tree node
    '''
    try:
        core_tag = gio.gio_read(file1, "core_tag")
        infall_tree_node_index = gio.gio_read(file1, "infall_tree_node_index")
        core_tag_unique = _is_unique_array(core_tag)
        tree_node_unique = _is_unique_array(infall_tree_node_index)
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




def check_mt_files(file1,file_mt,step):
    ''' Check that the number of halos in the merger tree file is consistent with the number of central cores
        criteria: number of halos in MT file and number of central cores agree to within 10%
        inputs: file1 = core catalog file
                file_mt = merger tree file 
                step = current time step
        outputs: criteria (pass=0)
        display: Number of halos in merger tree file 
                 Number of central cores in core file
        option: plot_missing
                if plot_missing, then print out the number of 'missing objects', check 
                how many of these have zero-masses, and plot the positions of the first 
                10,000 objects. (note at step=499 this shows boundary issues, at higher 
                redshifts this shows fragments)
    '''
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
    ''' Check that core tags present in an earlier timestep are retained in the later timestep
        criteria: less than 5% of the core tags present in the earlier timestep are missing in the later step
        inputs: file1 = core catalog file
                file2 = core catalog file at lower step number (higher redshift) 
        outputs: criteria (pass=0)
        display: Percentage of cores present in earlier timestep missing in subsequent timestep 
        option: plot_missing
                if plot_missing, then print out the number of 'missing objects', and plot the 
                positions of the first 10,000 objects alongside a core density map for comparison. 
    '''
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

