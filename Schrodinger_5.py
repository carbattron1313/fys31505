#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Dec 4 2022

@author: Elena Muñoz Rivas, Alejandro Carballido Mantecón, Antonio Gómez Garrido
and David Martínez Hernández.
"""

import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa


#Creation of the settings:
h = 0.005
M = int(1.0/h);
dt = 2.5e-5


#---------------------------------------------------------------------
#Variables of the time
T1 = 0.008
N1 = int(T1/dt)

#Load of the bin
U_0 = pa.cx_cube()
U_0.load("cube.bin")
U_10 = pa.cx_cube()
U_10.load("cube1.bin")

#Creation of the probabilities for different velocities:
p_0_array = np.zeros(N1)
p_10_array = np.zeros(N1)
for n in range(N1):
    p_0 = 0.0
    p_10 = 0.0
    for i in range(M-2):
        for j in range(M-2):
            p_0 += np.real(np.conj(U_0[i, j, n])*U_0[i, j, n])
            p_10 += np.real(np.conj(U_10[i, j, n])*U_10[i, j, n])
    p_0_array[n] = p_0
    p_10_array[n] = p_10

#Create vector of the time and the desviation
t = np.linspace(0, T1, N1)
des_0 = abs(1 - p_0_array)
des_10 = abs(1 - p_10_array)


#Plot the deviation of the total probability from 1.0 as a function of time
#when v = 0, and when v = 1e+10 and a double-split
plt.figure()
plt.plot(t, des_0, ".", label = "$v_0 = 0$")
plt.plot(t, des_10, ".", label = "$v_0 = 10^{10}$")
plt.xlabel("Time")
plt.ylabel("|1 - $\sum_{i, j}\ p_{ij}$|")
plt.yscale("log")
plt.legend()
plt.savefig('desviation_prob.pdf')

plt.show()

#------------------------------------------------------------------
#Variables of the time
T2 = 0.002
N2 = int(T2/dt)

#Load of the bin
U_2 = pa.cx_cube()
U_2.load("cube2.bin")

#Create of the real and imaginary arrays
U_Re0 = np.zeros((M-2, M-2))
U_Re0001 = np.zeros((M-2, M-2))
U_Re0002 = np.zeros((M-2, M-2))

U_Im0 = np.zeros((M-2, M-2))
U_Im0001 = np.zeros((M-2, M-2))
U_Im0002 = np.zeros((M-2, M-2))

p_0 = np.zeros((M-2, M-2))
p_0001 = np.zeros((M-2, M-2))
p_0002 = np.zeros((M-2, M-2))

#Upload the values of the arrays
for i in range(M-2):
    for j in range(M-2):
        #Probabilities
        p_0[i, j] = np.real(np.conj(U_2[i, j, 0])*U_2[i, j, 0])
        p_0001[i, j] = (np.conj(U_2[i, j, int(N2/2)])*U_2[i, j, int(N2/2)]).real
        p_0002[i, j] = (np.conj(U_2[i, j, N2-1])*U_2[i, j, N2-1]).real

        #Real values
        U_Re0[i, j] = U_2[i, j, 0].real
        U_Re0001[i, j] = U_2[i, j, int(N2/2)].real
        U_Re0002[i, j] = U_2[i, j, N2-1].real

        #Imaginary values
        U_Im0[i, j] = U_2[i, j, 0].imag
        U_Im0001[i, j] = U_2[i, j, int(N2/2)].imag
        U_Im0002[i, j] = U_2[i, j, N2-1].imag


#Values of the different probabilities and cubes imaginaries and reals:
time_string = ["0", "0.001", "0.002"]
p_title = ["$p_0$", "$p_{0001}$", "$p_{0002}$"]
U_Re_title = ["Re($U_0$)", "Re($U_{0001}$)", "Re($U_{0002}$)"]
U_Im_title = ["Im($U_0$)", "Im($U_{0001}$)", "Im($U_{0002}$)"]

prob = [p_0, p_0001, p_0002]
U_Real = [U_Re0, U_Re0001, U_Re0002]
U_Imag = [U_Im0, U_Im0001, U_Im0002]

#Arrays of the axis for the colourmap plot:
x = np.linspace(0+h, 1-h, M-2);
y = np.linspace(0+h, 1-h, M-2);
X, Y = np.meshgrid(x,y)


#Colourmap plots that illustrate the time evolution of the 2D probability function
#The for is needed because of the different times that we have:
for i in range(len(prob)):
    plt.contourf(X, Y, np.sqrt(prob[i]), 20)
    plt.xlabel("x")
    plt.ylabel("y")
    cb = plt.colorbar()
    cb.set_label(label = "p(x, y; t = " + time_string[i] + "$)^{1/2}$")
    plt.savefig("prob" + time_string[i] + ".pdf")
    plt.show()


for i in range(len(prob)):
    plt.contourf(X, Y, U_Real[i], 20)
    plt.xlabel("x")
    plt.ylabel("y")
    cb = plt.colorbar()
    cb.set_label(label="Re(u(x, y, t = " + time_string[i] + ")")
    plt.savefig("real" + time_string[i] + ".pdf")
    plt.show()


for i in range(len(prob)):
    plt.contourf(X, Y, U_Imag[i], 20)
    plt.xlabel("x")
    plt.ylabel("y")
    cb = plt.colorbar()
    cb.set_label(label="Im(u(x, y, t = " + time_string[i] + ")")
    plt.savefig("imag" + time_string[i] + ".pdf")
    plt.show()

    #Load of the values of hte cubes for single bins, double bins and triple
    U_1s = pa.cx_cube()
    U_1s.load("cube2_1s.bin")
    U_3s = pa.cx_cube()
    U_3s.load("cube2_3s.bin")


    p_sing = np.zeros((M-2, M-2))
    p_double = np.zeros((M-2, M-2))
    p_trip = np.zeros((M-2, M-2))

    #Range to create the arrays of the probabilities for the three slit
    for i in range(M-2):
        for j in range(M-2):
            p_sing[i, j] = (np.conj(U_1s[i, j, N2-1])*U_1s[i, j, N2-1]).real
            p_double[i, j] = (np.conj(U_2[i, j, N2-1])*U_2[i, j, N2-1]).real
            p_trip[i, j] = (np.conj(U_3s[i, j, N2-1])*U_3s[i, j, N2-1]).real


    x_index = int(np.where(x==0.8)[0])

    #1D detection probability along a screen with x=0.8 at t=0.002 for the
    #different splits.
    plt.figure()
    plt.plot(y, p_sing[:, x_index], label = "Single-slit experiment")
    plt.plot(y, p_double[:, x_index], label = "Double-slit experiment")
    plt.plot(y, p_trip[:, x_index], label = "Triple-slit experiment")
    plt.xlabel("y")
    plt.ylabel("p(y | x = 0.8; t = 0.002)")
    plt.legend()
    plt.savefig('detection_prob.pdf')
    plt.show()