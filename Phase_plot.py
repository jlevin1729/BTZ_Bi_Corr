#!/usr/bin/env python3
import numpy as np
import math as m
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as pat
import scipy.optimize as opt
import random as rand
import pickle

#retrieve dictionary of phases

file = open('dict_EQ','rb')
Dict_EQ = pickle.load(file)
file.close()

file = open('dict_EP','rb')
Dict_EP = pickle.load(file)
file.close()

#define function which plots the phase diagram for fixed A_center and rH

def phase_diagram(A_center, rH):

    fig, ax = plt.subplots()
    ax.set_xlim((0,40))
    ax.set_ylim((0,40))
    ax.set_aspect('equal')

    #ax.set_xlabel('A size [\u03C0/40 rad]')
    #ax.set_ylabel('B size [\u03C0/40 rad]')
    #ax.set_title(r'$E_P$' + ' phase diagram, r = 0.61, A center = 0')

    

    for k,v in Dict_EQ.items():

        if k[0] == A_center and k[3] == rH:

            if v == 1:

                phase = 'red'

            elif v == 2:

                phase = 'blue'

            elif v == 3:

                phase = 'green'

            elif v == 4:

                phase = 'yellow'

            elif v == 'EW disconn':

                phase = 'black'

            elif v == 'BH outside EW':

                phase = 'brown'

            ax.plot(k[1], k[2], '.', color = phase)

    for k,v in Dict_EP.items():

        if k[0] == A_center and k[3] == rH:

#            if v == 1:

#                ax.plot(k[1], k[2], '.', color = 'red')

            if v == 2:

                ax.plot(k[1], k[2], '.', color = 'orange')

#            elif v == 'EW disconn':

#                ax.plot(k[1], k[2], '.', color = 'black')

#            elif v == 'BH outside EW':

#                ax.plot(k[1], k[2], '.', color = 'brown')

    ax.tick_params(top = 'off', right = 'off', labeltop = 'off', labelbottom = 'off', labelright = 'off', labelleft = 'off' )

    plt.savefig('path' + str(rH) + '_A' + str(A_center) + '.png')
    plt.close()

    return

#loop over a parameter

#A_center = 0
#rH = 0.1

phase_diagram(3, 0.28)

phase_diagram(5, 0.28)

phase_diagram(6, 0.28)

phase_diagram(9, 0.28)

#for x in np.linspace(0, 40, 41):

    #phase_diagram(round(x,2), rH)

    
