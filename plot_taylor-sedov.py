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
    #dg = pg.GInterpModal(fluidData, 0, 'ms')
    #gridx, gridy, gridz = dg.interpolateGrid()
    fluid5m = fluidData.getValues()
    grids = fluidData.getGrid()

    rho = fluid5m[..., 0]
    rhoux = fluid5m[..., 1]
    rhouy = fluid5m[..., 2]
    rhouz = fluid5m[..., 3]
    e = fluid5m[..., 4] # 1/2 rho u^2 + p/(gas_gamma - 1)


    r_min = 0.10;
    r_max = 0.78;
    NX = 64;
    dr = (r_max - r_min)/NX;
    r_values = np.linspace(r_min + dr/2,r_max - dr/2,NX)

    plt.plot( r_values, rho[:,0,0] , label = 'rho')
    plt.plot( r_values, rhoux[:,0,0]/rho[:,0,0], label = 'u' )

    p = (e - 0.5*rhoux*rhoux/rho)*(gas_gamma - 1.0)

    plt.plot( r_values, p[:,0,0], label = 'p' )

    plt.legend(loc='best')

    plt.show()

    print(rho[:,0,0])
    print("----")
    print(p[:,0,0])
    
read_data('euler_taylorsedov_test', 'euler', int(sys.argv[1]), 5.0/3.0)
