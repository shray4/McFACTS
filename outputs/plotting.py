from pylab import *
import sys, os, time, string, math, subprocess
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import AxesGrid

def get_data(path, filename):
    # read in the data and dump to arrays
    
    # open the file for reading
    file1 = open(path+filename, 'r')

    #read and strip(?) headerdata
    #read in the data line by line, split into float lists
    center_of_mass_list=[]
    total_mass_list=[]
    chi_eff_list=[]
    a_tot_list=[]
    spin_angle_list=[]
    m1_list=[]
    m2_list=[]
    a1_list=[]
    a2_list=[]
    theta1_list=[]
    theta2_list=[]
    gen1_list=[]
    gen2_list=[]
    tmerge_list=[]
    header1=file1.readline()
    for line in file1:
        line=line.strip()
        columns=line.split()
        center_of_mass_list.append(float(columns[0]))
        total_mass_list.append(float(columns[1]))
        chi_eff_list.append(float(columns[2]))
        a_tot_list.append(float(columns[3]))
        spin_angle_list.append(float(columns[4]))
        m1_list.append(float(columns[5]))
        m2_list.append(float(columns[6]))
        a1_list.append(float(columns[7]))
        a2_list.append(float(columns[8]))
        theta1_list.append(float(columns[9]))
        theta2_list.append(float(columns[10]))
        gen1_list.append(float(columns[11]))
        gen2_list.append(float(columns[12]))
        tmerge_list.append(float(columns[13]))

    # close file
    file1.close()

    # re-cast as arrays (from lists) for manipulation
    center_of_mass = np.array(center_of_mass_list)
    total_mass = np.array(total_mass_list)
    chi_eff = np.array(chi_eff_list)
    a_tot = np.array(a_tot_list)
    spin_angle = np.array(spin_angle_list)
    m1 = np.array(m1_list)
    m2 = np.array(m2_list)
    a1 = np.array(a1_list)
    a2 = np.array(a2_list)
    theta1 = np.array(theta1_list)
    theta2 = np.array(theta2_list)
    gen1 = np.array(gen1_list)
    gen2 = np.array(gen2_list)
    t_merge = np.array(tmerge_list)

    return center_of_mass, total_mass, chi_eff, a_tot, spin_angle, m1, m2, a1, a2, theta1, theta2, gen1, gen2, t_merge

def make_q_chi_eff():
    # make q-chi_eff plot from all mergers
    # read in masses, chi_eff from output_mergers.dat
    # compute q (is m1 always > m2 by definition or is it random?)

    path = ''
    inputfiles = ['output_mergers.dat']
    center_of_mass=[]
    total_mass=[]
    chi_eff=[]
    a_tot=[]
    spin_angle=[]
    m1=[]
    m2=[]
    a1=[]
    a2=[]
    theta1=[]
    theta2=[]
    gen1=[]
    gen2=[]
    t_merge=[]
    for i in range(len(inputfiles)):
        fn=inputfiles[i]
        cm, ms, chi, a_mag, a_ang, ms1, ms2, a_1, a_2, theta_1, theta_2, gen_1, gen_2, tmrg=get_data(path, fn)
        center_of_mass.append(cm)
        total_mass.append(ms)
        chi_eff.append(chi)
        a_tot.append(a_mag)
        spin_angle.append(a_ang)
        m1.append(ms1)
        m2.append(ms2)
        a1.append(a_1)
        a2.append(a_2)
        theta1.append(theta_1)
        theta2.append(theta_2)
        gen1.append(gen_1)
        gen2.append(gen_2)
        t_merge.append(tmrg)

    plot_done = False
    print(plot_done)

    # BEGIN DISPLAY INPUTS:
    # format=left, bottom, width, height
    rect1=0.12,0.12,0.85,0.85
        
    # make figure
    fig1=plt.figure(1)
    # add axes
    ax1=fig1.add_axes(rect1)
    # label them
    # ax1.yaxis.set_label_coords(1.0, 0.05)
    ax1.set_ylabel(r"$q$", fontsize=18)
    # ax1.xaxis.set_label_coords(-1.0, 1.0)
    ax1.set_xlabel(r"$\chi_{eff}$", fontsize=18)

    # set up range for x-axis
    plt.xlim(-1.0, 1.0)
    # range for y-axis
    plt.ylim(0.05, 1.0)
    #END DISPLAY INPUTS

    # set up y-axis variables
    primary = np.maximum(m1, m2)
    secondary = np.minimum(m1, m2)

    mass_ratio = secondary/primary

    ax1.scatter(chi_eff, mass_ratio, color='black')

    savefig('q_chi_eff.png')
    if (1==1):
        plot_done = True

    return plot_done

if __name__ == "__main__":
    print("In main")
    ran = make_q_chi_eff()
    print(ran)