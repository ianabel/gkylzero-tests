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


    r_min = 0.10;
    r_max = 1.38;
    NX = 64;
    dr = (r_max - r_min)/64;
    r_values = np.linspace(r_min + dr/2,r_max-dr/2,NX)
    plt.plot( r_values, rho[:,0,0] )


    plt.show()
    
read_data('euler_taylorsedov_test', 'euler', int(sys.argv[1]), 5.0/3.0)
