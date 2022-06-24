##### preamble #####
from my_units import *
from secular_pt_param import *
import rebound
from scipy import stats
import os, sys

import numpy as np
import pandas as pd
import scipy as sp

from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from scipy.integrate import quad, nquad
from scipy.optimize import brentq

import matplotlib.pyplot as plt
import matplotlib as mp
from matplotlib import font_manager
from matplotlib import rcParams
from matplotlib import rc

plt.rcdefaults()
fontsize = 14
rcParams['font.family'] = 'sans-serif'
font_manager.findfont('serif', rebuild_if_missing=True)
rcParams.update({'font.size':fontsize})

from tqdm import tqdm, tqdm_notebook
from time import time as tictoc

##### orbital elements #####
def fn_e(h,k):
    return np.sqrt(h**2 + k**2)
def fn_omega_bar(h,k):
    return np.mod(np.arctan2(h,k),2*np.pi)
def fn_inc(p,q):
    """Inclination as a function of osculating elements p,q. Returns value between [0,np.pi/2], i.e. does not discriminate between prograde and retrograde orbits."""
    return np.arcsin(np.sqrt(p**2 + q**2)) 
def fn_Omega(p,q):
    return np.mod(np.arctan2(p,q),2*np.pi)
def fn_omega(h,k,p,q):
    return np.mod(fn_omega_bar(h,k)-fn_Omega(p,q),2*np.pi)

##### secular perturbations #####
def fn_alpha(a,a_i):
    return np.min([a/a_i,a_i/a])
def fn_alpha_bar(a,a_i):
    return np.min([a/a_i,1])

def fn_alpha_i(a,df_planets):
    return np.asarray([fn_alpha(a,a_i) for  a_i in df_planets['a']])
def fn_alpha_bar_i(a,df_planets):
    return np.asarray([fn_alpha_bar(a,a_i) for  a_i in df_planets['a']])

def fn_b(alpha,s,j):
    f = lambda psi,alpha,s,j : 1/np.pi * np.cos(j*psi) / (1-2*alpha*np.cos(psi)+alpha**2)**s
    return quad(f,0,2*np.pi,args=(alpha,s,j))[0]
def fn_b_i(a,s,j,df_planets):
    return np.asarray([fn_b(fn_alpha(a,a_i),s,j) for a_i in df_planets['a']])

def fn_secular_freq(a,df_planets,e_il,inc_il):
    alpha_i = fn_alpha_i(a,df_planets)
    alpha_bar_i = fn_alpha_bar_i(a,df_planets)
    b_i_1 = fn_b_i(a,3/2,1,df_planets)
    b_i_2 = fn_b_i(a,3/2,2,df_planets)
    n = a**(-3/2) # mean orbital angular velocity
    fac = np.asarray((n/4) * df_planets['m/M'] * fn_alpha_i(a,df_planets) * fn_alpha_bar_i(a,df_planets))
    A = np.sum(fac * fn_b_i(a,3/2,1,df_planets))
    A_i = -fac * fn_b_i(a,3/2,2,df_planets)
    B = -A
    B_i = fac * fn_b_i(a,3/2,1,df_planets)
    nu_l = A_i @ e_il
    mu_l = B_i @ inc_il
    return [A,B,A_i,B_i,nu_l,mu_l]

##### secular evolution #####
def fn_h_evol(a_0,e_0,inc_0,omega_0,Omega_0,t_0,t,df_planets,e_il,inc_il,g_l,f_l,beta_l,gamma_l):
    omega_bar_0 = Omega_0 + omega_0
    [A,B,A_i,B_i,nu_l,mu_l] = fn_secular_freq(a_0,df_planets,e_il,inc_il)
    return e_0 * np.sin(A*(t-t_0)+omega_bar_0) - np.sum( nu_l/(A-g_l) * ( np.sin(g_l*t+beta_l) - np.sin(A*(t-t_0)+g_l*t_0+beta_l) ) )
def fn_k_evol(a_0,e_0,inc_0,omega_0,Omega_0,t_0,t,df_planets,e_il,inc_il,g_l,f_l,beta_l,gamma_l):
    omega_bar_0 = Omega_0 + omega_0
    [A,B,A_i,B_i,nu_l,mu_l] = fn_secular_freq(a_0,df_planets,e_il,inc_il)
    return e_0 * np.cos(A*(t-t_0)+omega_bar_0) - np.sum( nu_l/(A-g_l) * ( np.cos(g_l*t+beta_l) - np.cos(A*(t-t_0)+g_l*t_0+beta_l) ) )
def fn_p_evol(a_0,e_0,inc_0,omega_0,Omega_0,t_0,t,df_planets,e_il,inc_il,g_l,f_l,beta_l,gamma_l):
    [A,B,A_i,B_i,nu_l,mu_l] = fn_secular_freq(a_0,df_planets,e_il,inc_il)
    return np.sin(inc_0) * np.sin(B*(t-t_0)+Omega_0) - np.sum( mu_l/(B-f_l) * ( np.sin(f_l*t+gamma_l) - np.sin(B*(t-t_0)+f_l*t_0+gamma_l) ) )
def fn_q_evol(a_0,e_0,inc_0,omega_0,Omega_0,t_0,t,df_planets,e_il,inc_il,g_l,f_l,beta_l,gamma_l):
    [A,B,A_i,B_i,nu_l,mu_l] = fn_secular_freq(a_0,df_planets,e_il,inc_il)
    return np.sin(inc_0) * np.cos(B*(t-t_0)+Omega_0) - np.sum( mu_l/(B-f_l) * ( np.cos(f_l*t+gamma_l) - np.cos(B*(t-t_0)+f_l*t_0+gamma_l) ) )

def fn_evol(e_0,inc_0,omega_0,Omega_0,a_0,t_0,t,df_planets,e_il,inc_il,g_l,f_l,beta_l,gamma_l):
    omega_bar_0 = Omega_0 + omega_0
    [A,B,A_i,B_i,nu_l,mu_l] = fn_secular_freq(a_0,df_planets,e_il,inc_il)
    if len(np.shape([t_0]))>1:
        h_evol = e_0 * np.sin(A*(t-t_0)+omega_bar_0) - np.sum(np.asarray([nu_l[i]/(A-g_l[i]) * ( np.sin(g_l[i]*t+beta_l[i]) - np.sin(A*(t-t_0)+g_l[i]*t_0+beta_l[i]) ) 
                                                                          for i in range(len(g_l))]),axis=0)
        k_evol = e_0 * np.cos(A*(t-t_0)+omega_bar_0) - np.sum(np.asarray([nu_l[i]/(A-g_l[i]) * ( np.cos(g_l[i]*t+beta_l[i]) - np.cos(A*(t-t_0)+g_l[i]*t_0+beta_l[i]) ) 
                                                                          for i in range(len(g_l))]),axis=0)
        p_evol = np.sin(inc_0) * np.sin(B*(t-t_0)+Omega_0) - np.sum(np.asarray([mu_l[i]/(B-f_l[i]) * ( np.sin(f_l[i]*t+gamma_l[i]) - np.sin(B*(t-t_0)+f_l[i]*t_0+gamma_l[i]) )
                                                                                for i in range(len(g_l))]),axis=0)
        q_evol = np.sin(inc_0) * np.cos(B*(t-t_0)+Omega_0) - np.sum(np.asarray([mu_l[i]/(B-f_l[i]) * ( np.cos(f_l[i]*t+gamma_l[i]) - np.cos(B*(t-t_0)+f_l[i]*t_0+gamma_l[i]) )
                                                                                for i in range(len(g_l))]),axis=0)
    else:
        h_evol = e_0 * np.sin(A*(t-t_0)+omega_bar_0) - np.sum( nu_l/(A-g_l) * ( np.sin(g_l*t+beta_l) - np.sin(A*(t-t_0)+g_l*t_0+beta_l) ) )
        k_evol = e_0 * np.cos(A*(t-t_0)+omega_bar_0) - np.sum( nu_l/(A-g_l) * ( np.cos(g_l*t+beta_l) - np.cos(A*(t-t_0)+g_l*t_0+beta_l) ) )
        p_evol = np.sin(inc_0) * np.sin(B*(t-t_0)+Omega_0) - np.sum( mu_l/(B-f_l) * ( np.sin(f_l*t+gamma_l) - np.sin(B*(t-t_0)+f_l*t_0+gamma_l) ) )
        q_evol = np.sin(inc_0) * np.cos(B*(t-t_0)+Omega_0) - np.sum( mu_l/(B-f_l) * ( np.cos(f_l*t+gamma_l) - np.cos(B*(t-t_0)+f_l*t_0+gamma_l) ) )
    e_evol = fn_e(h_evol,k_evol)
    inc_evol = fn_inc(p_evol,q_evol)
    omega_evol = fn_omega(h_evol,k_evol,p_evol,q_evol)
    Omega_evol = fn_Omega(p_evol,q_evol)
    return [e_evol,inc_evol,omega_evol,Omega_evol]

def fn_mat_rot_inv_E(t=0,df_planets=df_planets):
    inc = df_planets[df_planets['name']=='Earth']['inc'].to_numpy()[0]
    omega = df_planets[df_planets['name']=='Earth']['omega'].to_numpy()[0]
    Omega = df_planets[df_planets['name']=='Earth']['Omega'].to_numpy()[0]
    mat_omega = np.asarray([[np.cos(omega),np.sin(omega),0],[-np.sin(omega),np.cos(omega),0],[0,0,1]])
    mat_inc = np.asarray([[1,0,0],[0,np.cos(inc),np.sin(inc)],[0,-np.sin(inc),np.cos(inc)]])
    mat_Omega = np.asarray([[np.cos(Omega),np.sin(Omega),0],[-np.sin(Omega),np.cos(Omega),0],[0,0,1]])
    return mat_omega @ mat_inc @ mat_Omega

def fn_mat_rot(inc,omega,Omega):
    mat_omega = np.asarray([[np.cos(omega),-np.sin(omega),0],[np.sin(omega),np.cos(omega),0],[0,0,1]])
    mat_inc = np.asarray([[1,0,0],[0,np.cos(inc),-np.sin(inc)],[0,np.sin(inc),np.cos(inc)]])
    mat_Omega = np.asarray([[np.cos(Omega),-np.sin(Omega),0],[np.sin(Omega),np.cos(Omega),0],[0,0,1]])
    return mat_Omega @ mat_inc @ mat_omega

##### generate points along orbit ######
def fn_gen_pts_xyz(a,e,n_points):
    """Generate points along particle's orbit equidistant in time."""
    M = np.arange(0,2*np.pi,2*np.pi/n_points)
    theta = np.asarray([rebound.M_to_f(e,M[i]) for i in range(len(M))]) 
    r = a * (1-e**2) / (1+e*np.cos(theta))
    return np.asarray([r*np.cos(theta),r*np.sin(theta),0*np.sin(theta)])

def fn_gen_pts_xyz_E(e_0,inc_0,omega_0,Omega_0,a_0,t_0,t,df_planets,e_il,inc_il,g_l,f_l,beta_l,gamma_l,n_points=100):
    """Generate points along particle's orbit, evolve forward in time via secular perturbation theory, rotate to Earth's frame."""
    pts_xyz_E = np.zeros((len(t_0),n_points,3))
    pts_xyz = fn_gen_pts_xyz(a_0,e_0,n_points) # generate points along particle's orbit
    # possible forward evolutions, for different values of t_0
    [e_evol,inc_evol,omega_evol,Omega_evol] = fn_evol(e_0,inc_0,omega_0,Omega_0,a_0,t_0,t,df_planets,e_il,inc_il,g_l,f_l,beta_l,gamma_l)
    mat_rot_inv_E = fn_mat_rot_inv_E(t,df_planets) # rotation matrix to Earth's plane
    mat_rot = np.asarray([fn_mat_rot(inc_evol[i],omega_evol[i],Omega_evol[i]) for i in range(len(t_0))])
    for i in range(len(t_0)):
        pts_xyz_E[i] = np.transpose(mat_rot_inv_E @ mat_rot[i] @ pts_xyz)
    return pts_xyz_E.reshape(len(t_0)*n_points,3)