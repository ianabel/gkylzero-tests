#!/usr/bin/env python3

import sys
sys.path.append("/home/ian/projects/postgkyl/")

import matplotlib
import matplotlib.pyplot as plt
import postgkyl as pg
import numpy as np
import tables

def read_data(pre, species, frame, gas_gamma):
    fluidData = pg.GData('%s-%s_%d.gkyl' % (pre, species, frame), mapc2p_name='%s-mapc2p.gkyl' % pre)
    dg = pg.GInterpModal(fluidData, 0, 'ms')
    gridx, gridy, gridz = dg.interpolateGrid()
    fluid5m = fluidData.getValues()

    rho = fluid5m[..., 0]
    rhoux = fluid5m[..., 1]
    rhouy = fluid5m[..., 2]
    rhouz = fluid5m[..., 3]
    e = fluid5m[..., 4] # 1/2 rho u^2 + p/(gas_gamma - 1)

    # gridx, gridy are vertex locations now
    plt.figure(1)
    plt.pcolormesh(gridx, gridy, rho[...,0])

    #center the grid values
    for d in range(len(gridx)):
        gridx[d] = 0.5*(gridx[d] + gridx[d+1])
        gridy[d] = 0.5*(gridy[d] + gridy[d+1])
        gridz[d] = 0.5*(gridz[d] + gridz[d+1])

    # gridx, gridy are cell center locations now
    plt.figure(2)
    plt.pcolormesh(gridx, gridy, rho[...,0] )        
    plt.show()
    
read_data('euler_taylorsedov_test', 'euler', 0, 5.0/3.0)
