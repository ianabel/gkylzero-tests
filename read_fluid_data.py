import matplotlib
import matplotlib.pyplot as plt
import postgkyl as pg
import numpy as np
import tables

def read_data(pre, species, frame, gas_gamma):
    fluidData = pg.GData('%s-%s_%d.gkyl' % (pre, species, frame), mapc2p_name='%s-mapc2p.gkyl' % pre)
    dg = pg.GInterpModal(fluidData, 0, 'ms')
    gridx, gridy = dg.interpolateGrid()
    fluid5m = fluidData.getValues()

    rho = fluid5m[..., 0]
    rhoux = fluid5m[..., 1]
    rhouy = fluid5m[..., 2]
    rhouz = fluid5m[..., 3]
    e = fluid5m[..., 4] # 1/2 rho u^2 + p/(gas_gamma - 1)

    # gridx, gridy are vertex locations now
    plt.figure(1)
    plt.pcolormesh(gridx, gridy, rho)

    #center the grid values
    for d in range(len(gridx)):
        gridx[d] = 0.5*(gridx[d] + gridx[d])
        gridy[d] = 0.5*(gridy[d] + gridy[d])

    # gridx, gridy are cell center locations now
    plt.figure(2)
    plt.pcolormesh(gridx, gridy, rho)        
    plt.show()
    
read_data('euler_wedge_sodshock', 'euler', 1, 5.0/3.0)
