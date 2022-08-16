#!/usr/bin/env python3

import TaylorSedov as TS
import numpy as np
import sys

Time = float(sys.argv[1])
r    = float(sys.argv[2])
IE   = float(sys.argv[3])

tsProblem = TS.TaylorSedov( 1.0, IE, 5.0/3.0 )

print("At t = %f, the shock front is at r = %.8f" % (Time,tsProblem.ShockR(Time)))

print(" Rho = %.8f" % tsProblem.TaylorSedovRho(Time,r))
print(" u_r = %.8f" % tsProblem.TaylorSedovU(Time,r))
print(" P   = %.8f" % tsProblem.TaylorSedovP(Time,r))



