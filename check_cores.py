#!/bin/bash/python

import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import scipy
from scipy.stats import binned_statistic as bs

cooley=False
if cooley:
    sys.path.append('/home/prlarsen/usr/genericio/python')
else:
    sys.path.append('/home/prlarsen/genericio2/genericio/python')
import genericio as gio

# global variable if you want visualisations
vis = True


def norm_mass(file1):
    '''
    Checking that core abundance increases with host halo mass as a power law with index unity
    '''

    infall_mass = gio.gio_read(file1,'infall_mass').flatten()
    infall_step = gio.gio_read(file1,'infall_step').flatten()
    host_tag = gio.gio_read(file1,'fof_halo_tag').flatten()
    host_mass = infall_mass
    infall_mass=infall_mass[host_mass>=1.e12]
    host_tag = host_tag[host_mass>=1.e12]
    infall_step = infall_step[host_mass>=1.e12]
    host_tag_u = host_tag[infall_step==499]
    host_mass = host_mass[host_mass>=1.e12][infall_step==499]
    dict_values = dict(zip(host_tag_u, host_mass))
    def apply_dict(i):
        return dict_values[i]


    infall_mass = gio.gio_read(file1,'infall_mass').flatten()
    infall_step = gio.gio_read(file1,'infall_step').flatten()
    host_tag = gio.gio_read(file1,'fof_halo_tag').flatten()
    mask = (infall_mass>1.e12)&(infall_step!=499)
    tags,counts = np.unique(host_tag[mask],return_counts=True)
    mass_tot =  np.array(map(apply_dict,tags))

    mean_sats, edges, nums = bs(np.log10(mass_tot),counts-1.,statistic='mean',bins=100)
    bin_mean = (edges[1:]-edges[:-1])/2.+edges[:-1]

    idx_13 = np.where(np.abs(bin_mean-13.5)==np.min(np.abs(bin_mean-13.5)))
    idx_15 = np.where(np.abs(bin_mean-14.5)==np.min(np.abs(bin_mean-14.5)))
    diff_val = (np.log10(mean_sats[idx_15])-np.log10(mean_sats[idx_13]))

    criteria = (diff_val<1.2)&(diff_val>0.8)

    print((np.log10(mean_sats[idx_15])-np.log10(mean_sats[idx_13])))
    if vis:
        plt.figure()
        plt.plot(bin_mean,mean_sats)
        plt.ylabel('N_cores (>1.e12 mass)')
        plt.yscale('log')
        plt.xlabel('log10(host mass)')
        plt.xlim([12.5,15])
        plt.show()
    if criteria:
        print("power law index of core abundance as a function of host halo mass has an index of nearly unity")
        return 0
    else:
        print("core abundance as a function of host halo mass needs checking")
        return 1 

def dndmu(file1):
    print file1
    gio.gio_inspect(file1)
    infall_mass = gio.gio_read(file1,'infall_mass').flatten()
    infall_step = gio.gio_read(file1,'infall_step').flatten()
    host_tag = gio.gio_read(file1,'fof_halo_tag').flatten()
    host_mass = infall_mass
    infall_mass=infall_mass[host_mass>=1.e12]
    host_tag = host_tag[host_mass>=1.e12]
    infall_step = infall_step[host_mass>=1.e12]
    host_tag_u = host_tag[infall_step==499]
    host_mass = host_mass[host_mass>=1.e12][infall_step==499]
    dict_values = dict(zip(host_tag_u, host_mass))
    
    infall_mass = infall_mass[infall_step!=499]
    host_tag = host_tag[infall_step!=499]

    def apply_dict(i):
        return dict_values[i]

    print(host_tag[0])
    print(dict_values[host_tag[0]])
    print(apply_dict(host_tag[0]))
    mass_tot =  np.array(map(apply_dict,host_tag))
    ratio_tot = infall_mass/mass_tot
    print(len(mass_tot))
    print(len(ratio_tot)) 
    print(mass_tot[:100])
    '''
    print(len(host_tag_u))
    ratio_tot=np.array([]).astype(float)
    mass_tot=np.array([]).astype(float)
    for i in range(len(host_tag_u)):
        if i%1000==0:
            print(i)
        core_masses = infall_mass[(host_tag==host_tag_u[i])&(infall_step!=499)]
        ratio = core_masses/host_mass[i]
        mass = ratio*0.0 + host_mass[i]
        ratio_tot = np.concatenate((ratio_tot,ratio))
        mass_tot = np.concatenate((mass_tot,mass))
    print(len(host_mass))
    '''
    bin12 = (mass_tot>1.e12)&(mass_tot<5.e12)
    bin13 = (mass_tot>1.e13)&(mass_tot<5.e13)
    bin14 = (mass_tot>1.e14)&(mass_tot<5.e14)
    bin125 = (mass_tot>5.e12)&(mass_tot<1.e13)
    bin135 = (mass_tot>5.e13)&(mass_tot<1.e14)
    bin145 = (mass_tot>5.e14)&(mass_tot<1.e15)

    print(np.sum(bin12),np.sum(bin125),np.sum(bin13),np.sum(bin135),np.sum(bin14),np.sum(bin145))
    bins = np.logspace(-2,0,100)
    hist12, edges = np.histogram(ratio_tot[bin12], normed=True,bins = bins ) #np.logspace(-1,0,100))#np.linspace(0.0,1.0,100))
    hist13, edges = np.histogram(ratio_tot[bin13], normed=True,bins = bins ) #np.logspace(-1,0,100))#linspace(0.0,1.0,100))
    hist14, edges = np.histogram(ratio_tot[bin14], normed=False,bins = bins )#np.logspace(-1,0,100))#inspace(0.0,1.0,100))
    hist125, edges = np.histogram(ratio_tot[bin125], normed=True,bins = bins ) # np.logspace(-1,0,100))#np.linspace(0.0,1.0,100))
    hist135, edges = np.histogram(ratio_tot[bin135], normed=False,bins = bins )#np.logspace(-1,0,100))#linspace(0.0,1.0,100))
    hist145, edges = np.histogram(ratio_tot[bin145], normed=True,bins = bins)#np.logspace(-3,0,100))#inspace(0.0,1.0,100))

    #plt.plot(edges[1:], hist12,label='bin 12')
    #plt.plot(edges[1:], hist13,label = 'bin 13')
    plt.plot(edges[1:], hist14,label= 'bin 14')
    #plt.plot(edges[1:], hist125,label='bin 12.5')
    plt.plot(edges[1:], hist135,label = 'bin 13.5')
    #plt.plot(edges[1:], hist145,label= 'bin 14.5')
    plt.yscale('log')
    plt.xscale('log')
    plt.legend()
    plt.xlabel('mass ratio')
    plt.ylabel('N_cores (normalized)')
    plt.show()

    bins = np.logspace(-0.5,0,100)
    hist12, edges = np.histogram(ratio_tot[bin12], normed=True,bins = bins ) #np.logspace(-1,0,100))#np.linspace(0.0,1.0,100))
    hist13, edges = np.histogram(ratio_tot[bin13], normed=True,bins = bins ) #np.logspace(-1,0,100))#linspace(0.0,1.0,100))
    hist14, edges = np.histogram(ratio_tot[bin14], normed=True,bins = bins )#np.logspace(-1,0,100))#inspace(0.0,1.0,100))
    hist125, edges = np.histogram(ratio_tot[bin125], normed=True,bins = bins ) # np.logspace(-1,0,100))#np.linspace(0.0,1.0,100))
    hist135, edges = np.histogram(ratio_tot[bin135], normed=True,bins = bins )#np.logspace(-1,0,100))#linspace(0.0,1.0,100))
    hist145, edges = np.histogram(ratio_tot[bin145], normed=True,bins = bins)#np.logspace(-3,0,100))#inspace(0.0,1.0,100))
    plt.plot(edges[1:], hist12,label='bin 12')
    plt.plot(edges[1:], hist13,label = 'bin 13')
   # plt.plot(edges[1:], hist14,label= 'bin 14')
    plt.plot(edges[1:], hist125,label='bin 12.5')
    plt.plot(edges[1:], hist135,label = 'bin 13.5')
    #plt.plot(edges[1:], hist145,label= 'bin 14.5')
    plt.yscale('log')
    plt.xscale('log')
    plt.legend()
    plt.xlabel('mass ratio')
    plt.ylabel('N_cores (normalized)')
    plt.show()


   # except:
    #    print('test_failed')
    return


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
                plt.yscale('log')
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

def check_heirarchy(file1):
    '''
    Check that host core tag is within the same tree node index as the core tag. 
    Check that infall step of the host is higher than infall step of the core
    '''
    
    core_tag = gio.gio_read(file1, "core_tag")
    host_core_tag = gio.gio_read(file1, "host_core")
    vals,ar1,ar2 = np.intersect1d(core_tag, host_core_tag, assume_unique=False, return_indices=True)     
    tni = gio.gio_read(file1,'tree_node_index')
    infall_step = gio.gio_read(file1,'infall_step')
    if (tni[ar1]==tni[ar2]).all()&(infall_step[ar1]>infall_step[ar2]).all():
        return 
    else:
       print(np.sum(tni[ar1]!=tni[ar2]),np.sum(infall_step[ar1]<=infall_step[ar2]))
       return 1

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
        plot_missing = True # useful to turn on if step = 499, so you can see if all missing objects are at mass resolution
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
            for i in range(min(1000,len(missing_cores))):
                #print(i)
                idx = np.where(halo_tag==missing_cores[i])[0][0]
                #print(idx)
                masses.append(halo_mass[idx]); 
            masses = np.array(masses)
            if 0 in masses:
                print("Fraction of zero-mass objects = "+str(np.sum(masses==0)/(len(masses)+0.0)))
            plt.figure()
            plt.hist((masses[masses>0]/1.09e9),bins=100)
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


def _halo_central_flip_light(x,y,z,infall_step,infall_mass,core_tag,step,halo_tag):
    ''' For an individual halo, check proximity of massive cores to central core '''
    xc = x[infall_step==step]
    yc = y[infall_step==step]
    zc = z[infall_step==step]
    mc = infall_mass[infall_step==step]
    mask = (infall_step!=step)&(infall_mass>0.2*mc)

    dist = np.sqrt((x[mask]-xc)**2+(y[mask]-yc)**2+(z[mask]-zc)**2)
    mindist=0; flip = False

    if (len(dist)>0):
        if (np.min(dist)<0.03):
            flip = True; mindist = np.min(dist)
    return flip, mindist


def _halo_central_flip_detailed(x,y,z,infall_step,infall_mass,core_tag,step,x2,y2,z2,core_tag2,infall_mass2,halo_tag):
    ''' For an individual halo, check proximity of massive cores to central core '''
    xc = x[infall_step==step]
    yc = y[infall_step==step]
    zc = z[infall_step==step]
    mc = infall_mass[infall_step==step]

    mask = (infall_step!=step)&(infall_mass>0.2*mc)
    dist = np.sqrt((x[mask]-xc)**2+(y[mask]-yc)**2+(z[mask]-zc)**2)
    mindist=0; flip = False

    if (len(dist)>0):
        if (np.min(dist)<0.01):
            flip = True; mindist = np.min(dist)
            print('odd halo found, tag = ' + str(halo_tag))
            print('minimum distance = '+str(np.min(dist))+' Mpc')
            print(infall_step[mask][dist==np.min(dist)])
            print(infall_mass[mask][dist==np.min(dist)]/mc)
            x3=[];y3=[];z3=[];m3=[];
            for i in range(len(core_tag)):
                x3.append(x2[core_tag2==core_tag[i]])
                y3.append(y2[core_tag2==core_tag[i]])
                z3.append(z2[core_tag2==core_tag[i]])
                m3.append(infall_mass2[core_tag2==core_tag[i]])
            x3 = np.array(x3); y3 = np.array(y3); z3 = np.array(z3); m3 = np.array(m3)
            if vis:
                plt.figure()
                plt.scatter(x,y,c=np.log10(infall_mass),s=1.0)
                plt.colorbar()
                plt.title('later step')
                plt.figure()
                plt.scatter(x3,y3,c=np.log10(m3),s=1.0)
                plt.colorbar()
                plt.title('earlier step')
                plt.show()
    return flip, mindist   

def central_flip(file1,file2,step):
    ''' Check for flipping of centrals from one step to the next 
        This test is in development, and may need to read in the accum files to see if there is truly missing substructures
        
    '''
    try:
        x = gio.gio_read(file1,'x').flatten()
        y = gio.gio_read(file1,'y').flatten()
        z = gio.gio_read(file1,'z').flatten()
        halo_tag = gio.gio_read(file1,'fof_halo_tag').flatten()
        infall_step = gio.gio_read(file1,'infall_step').flatten()
        infall_mass = gio.gio_read(file1,'infall_mass').flatten()
        core_tag = gio.gio_read(file1, "core_tag").flatten()

        core_tag2 = gio.gio_read(file2, "core_tag").flatten()
        x2 = gio.gio_read(file2,'x').flatten()
        y2 = gio.gio_read(file2,'y').flatten()
        z2 = gio.gio_read(file2,'z').flatten()
        infall_mass2 = gio.gio_read(file2,'infall_mass').flatten()


        centrals = (infall_step==step)
        fof_mass_mask = (infall_mass[centrals]>500.*1.1e9) # at least 500 particles (let's not bother with low mass halos for now)
        nhalos = len(halo_tag[centrals][fof_mass_mask])
        print("Number of halos above mass cut = "+str(nhalos))
        count = 0 ; maxcount = min(nhalos,5000)
        distarr=[]
        for j in range(maxcount):
            fof_halo_tag = halo_tag[centrals][fof_mass_mask][j]
            mask_fof = (halo_tag==fof_halo_tag)
            flip,dist = _halo_central_flip_light(x[mask_fof],y[mask_fof],z[mask_fof],infall_step[mask_fof],infall_mass[mask_fof],core_tag[mask_fof],step,fof_halo_tag)
            # uncomment if you want for the worst cases to see comparison plots with the previous step
            #flip,dist = _halo_central_flip_detailed(x[mask_fof],y[mask_fof],z[mask_fof],infall_step[mask_fof],infall_mass[mask_fof],core_tag[mask_fof],step,x2,y2,z2,core_tag2,infall_mass2,fof_halo_tag)
            if flip:
                count+=1
                distarr.append(dist)
        distarr = np.array(distarr)
        np.savetxt('/media/luna1/prlarsen/dist1.txt',np.array(distarr))
        print("Number of pairs closer than 30kpc = "+str(np.sum(distarr<0.03))+" out of a possible "+str(maxcount))
    except:
        print("test_failed")
    return 
       
 
def relative_quantities():
    return #xc =        


if __name__ == "__main__":

    check_mt = False
    check_ev = False
    # old version (note merger tree files are on alcf systems for this):
    #base_name = '/media/luna1/dkorytov/data/AlphaQ/core_catalog5/07_13_17.AlphaQ.'
    # new version
    opt = sys.argv[3]
    if opt=='0.1':
        base_name = '/media/luna1/prlarsen/core_catalog_0.1/09_03_2019.AQ.' 
    elif opt=='old':
        base_name = '/media/luna1/dkorytov/data/AlphaQ/core_catalog5/07_13_17.AlphaQ.'
    elif opt=='0.2':
        base_name = '/media/luna1/prlarsen/core_catalog_0.2/09_03_2019.AQ.'
    elif opt=='new':
        base_name='/media/luna1/prlarsen/core_stuff/sep11/core_catalog/09_03_2019.AQ.'
    elif opt=='merge':
        base_name = '/media/luna1/prlarsen/core_stuff/oct3/core_catalog_merg/09_03_2019.AQ.'
    else:
        base_name = '/media/luna1/prlarsen/core_stuff/09_03_2019.AQ.' 
    base_name_mt = '/media/luna1/prlarsen/core_stuff/output/09_03_2019.AQ.'
    step = sys.argv[1]
    step2 = sys.argv[2] # note step2 should be earlier in time than step1
    file1 = base_name+step+'.coreproperties'
    file_mt = base_name_mt +step+'.treenodes'
    file_ev = base_name + step2 + '.coreproperties'
    check_heirarchy(file1)

    stop
    test = 0
    norm_mass(file1)
    dndmu(file1)
    test += check_infall_masses(file1)
    test += check_infall_timesteps(file1)
    test += check_positions(file1)
    test += check_central(file1,int(step))
    test += check_unique_core_tags(file1)
    if check_mt:
        test += check_mt_files(file1,file_mt,int(step))
    if check_ev:
        test += check_core_evolution(file1,file_ev)

    #central_flip(file1,file_ev,int(step))

    print("Total number of issues found = "+str(test))

