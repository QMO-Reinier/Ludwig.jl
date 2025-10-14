#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 11:49:24 2025

Code translated from Sr2RuO4.jl to python

@author: lion
"""

import numpy as np

# ===========================================
# Constants
# ===========================================

material = "Sr2RuO4"
c = 12.68e-10  # Interlayer distance (m)
a = 3.90e-10   # Lattice constant (m)

# γ band
t_gamma = 0.07735      # eV
tp_gamma = 0.0303212   # eV
mu_gamma = 0.114478    # eV

def ham_gamma(k):
    kx, ky = k
    return -2.0 * t_gamma * (np.cos(2*np.pi*kx) + np.cos(2*np.pi*ky)) \
           - 4 * tp_gamma * np.cos(2*np.pi*kx) * np.cos(2*np.pi*ky) - mu_gamma


alpha = 7.604
alpha_p = 7.604

# α and β bands
t_alpha = 0.099   # eV
t_beta = 0.1278   # eV
t3 = 0.08
t5 = 0.13
mu_alpha_beta = 1.08


# ===========================================
# Band structure functions
# ===========================================

def exz(k):
    kx, ky = k
    return -2.0 * np.cos(2*np.pi*kx) - 2.0 * t3 * np.cos(2*np.pi*ky) - mu_alpha_beta

def eyz(k):
    kx, ky = k
    return -2.0 * t3 * np.cos(2*np.pi*kx) - 2.0 * np.cos(2*np.pi*ky) - mu_alpha_beta

def V(k):
    kx, ky = k
    return 4.0 * t5 * np.sin(2*np.pi*kx) * np.sin(2*np.pi*ky)

def ham_alpha(k):
    x, y = exz(k), eyz(k)
    return 0.5 * ((x + y) - np.sqrt((x - y)**2 + 4 * V(k)**2)) * t_alpha

def ham_beta(k):
    x, y = exz(k), eyz(k)
    return 0.5 * ((x + y) + np.sqrt((x - y)**2 + 4 * V(k)**2)) * t_beta


# ===========================================
# Deformation potentials
# ===========================================

def dii_mu(k, i, mu, delta=0.0):
    kx, ky = k
    if mu in (1, 2):
        delta /= t_alpha if mu == 1 else t_beta
        x, y = exz(k), eyz(k)
        Delta = np.sqrt(0.25 * (x - y)**2 + V(k)**2)

        term1 = alpha * (1 + t3) * np.cos(2*np.pi*k[i])
        term2 = ((-1)**(mu+1)) * (((-1)**(i+1)) * (1 - t3) * 
                                  (alpha * (1 - t3) * np.cos(2*np.pi*k[i]) - delta / 2.0))
