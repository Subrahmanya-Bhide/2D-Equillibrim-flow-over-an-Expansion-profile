# -*- coding: utf-8 -*-
"""
Created on Sun Apr  4 16:04:05 2021

@author: subbu
"""
import numpy as np
import matplotlib.pyplot as plt
from HTGD_TANEHILL_FUNCTIONS import T_prho
from HTGD_TANEHILL_FUNCTIONS import s_erho
from HTGD_TANEHILL_FUNCTIONS import reals_erho
from HTGD_TANEHILL_FUNCTIONS import a_ps
from HTGD_TANEHILL_FUNCTIONS import h_prho
from HTGD_TANEHILL_FUNCTIONS import e_ps
from HTGD_TANEHILL_FUNCTIONS import p_erho

#number of divisions of theta
n = 100

# Gas constant for Air
R = 287

#defining reference conditions 
p0 = 1.032*10**5 #reference pressure
T0 = 273.15      # reference temperature
rho0 = 1.2       #reference density
a0 = np.sqrt(1.4*R*T0) #reference speed of sound
e0 = 215000      #reference internal energy obtained from functional evaluations
s0 = reals_erho(e0,rho0,rho0,T0) #reference entropy
#angle of rotation of flow
THETA = 40
theta = THETA*np.pi/180 #in radians

#Incremental value of theta used in the calculations
dtheta = theta/n

#Array of thetas
Theta = np.arange(0,theta,dtheta)

#defining inlet conditions
#inlet pressure
p_in =  4*p0
p_stop = 10
n_p_step = 500
#inlet density
rho_in =  0.35*rho0
rho_stop = 0.00001
n_rho_step = 500
#inlet temperature
Tin = T_prho(p_in,rho_in,p0,rho0)
print(Tin)
# Array of internal energies for the search process
Es = np.linspace(10,215000*100,1000)

clist = []
for i in range(len(Es)):
    pcur = p_erho(Es[i],rho_in,T0,rho0)
    clist.append(abs(p_in - pcur))
c = clist.index(min(clist))
#inlet internal energy
e_in = Es[c]
#inlet Entropy
s_in = reals_erho(e_in,rho_in,rho0,T0)
#inlet enthalpy
hin = h_prho(p_in, rho_in, p0, rho0)
#speed of sound for the inlet
a_in = a_ps(p_in,s_in,p0,a0)
#inlet mach number
M_num = 1.8
#inlet velocity
Vin = a_in*M_num
#inlet total enthalpy
ht = hin + Vin**2/2

#function to find pressure given a enthalpy and entropy
def fin_p(h_true,s,p0,T0,s0):
    blist = []
    for i in range(len(Ps)):
        p_guess = Ps[i]
        e_cur = e_ps(p_guess,s,p0,T0,e0,s0)
        rho_cur = find_rho(e_cur,s,rhos)
        h_cur = h_prho(p_guess,rho_cur,p0,rho0)
        blist.append(abs(h_cur - h_true))
    b = blist.index(min(blist))
    return Ps[b]

#function to find density given an internal energy and entropy
def find_rho(e,s,rho_list):
    alist = []
    for i in range(len(rho_list)):
        alist.append(abs(s_erho(e,rho_list[i],rho0,T0,s)))
    a = alist.index(min(alist))
    return rho_list[a]

#List of density for search operation
rhos = np.linspace(rho_in,rho_stop,n_rho_step)
#List of pressure for search operation
Ps = np.linspace(p_in,p_stop,n_p_step)

#Velocity List
vel_array = [Vin]
#Enthalpy list
h = [hin]
#List for speed of sound
av = [a_in]
#List for pressure 
p = [p_in]
#Lsit for density
rho = [rho_in]
#List of temperature
T = [Tin]
#List for Mach number
M = [Vin/a_in]
#List for Mach angles
mu = [np.rad2deg(np.arcsin(1/M[0]))]

#Loop for the search operation and updation of values
for i in range(1,int(theta/dtheta)):
    sq = (((vel_array[i-1])/(av[i-1]))**2 - 1)**0.5
    vel_array.append(vel_array[i-1]  +   (dtheta*vel_array[i-1])/sq)
    h.append(ht - (vel_array[i]**2)/2)
    p.append(fin_p(h[i],s_in,p0,T0,s0))
    av.append(a_ps(p[i],s_in,p0,a0))
    rho.append(find_rho(e_ps(p[i],s_in,p0,T0,e0,s0),s_in,rhos))
    T.append(T_prho(p[i],rho[i],p0,rho0))
    M.append(vel_array[i]/av[i])
    mu.append(np.rad2deg(np.arcsin(1/M[i])))

#Plotting variation of properties with theta        
plt.plot(np.rad2deg(Theta),rho)
plt.xlabel('Theta (deg)')
plt.ylabel('Density') 

#plotting the Mach waves
X = np.array([0,1,2,3,4,5,6,7,8,9])
Y = np.zeros((n,10))

for i in range(n):
    for j in range(len(X)):
        Y[i,j] = np.tan(np.deg2rad(mu[i]))*X[j]
   
for i in range(n):
    plt.plot(X,Y[i])



#Function to give out values of properties of a given point
def prop(x0,y0):
    slope = np.rad2deg(np.arctan(y0/x0))
    for i in range(n-1):
        a = mu[i] - slope
        b = slope - mu[i+1]
        if a*b >= 0 :
            print('Pressure at the point is', p[i])
            print('Density at the point is', rho[i])
            print('Temperature at the point is', T[i])
            print('Velocity at the point is', vel_array[i])
            print('Speed of Sound at the point is', av[i])
            print('Mach Number at the point is', M[i])
            print('Enthalpy at the point is', h[i])
            break
        elif slope> mu[0]:
            print('Pressure at the point is', p[0])
            print('Density at the point is', rho[0])
            print('Temperature at the point is', T[0])
            print('Velocity at the point is', vel_array[0])
            print('Speed of Sound at the point is', av[0])
            print('Mach Number at the point is', M[0])
            print('Enthalpy at the point is', h[0])
            break
        elif slope < mu[-1]:
            print('Pressure at the point is', p[-1])
            print('Density at the point is', rho[-1])
            print('Temperature at the point is', T[-1])
            print('Velocity at the point is', vel_array[-1])
            print('Speed of Sound at the point is', av[-1])
            print('Mach Number at the point is', M[-1])
            print('Enthalpy at the point is', h[-1])
            break            
    return 1

#Example for the above function
prop(1,8)





 
    


















    
        
        
    

 
    
    




