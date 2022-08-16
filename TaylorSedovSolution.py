#!/usr/bin/env python3

import TaylorSedov as TS
import numpy as np
import sys

Time = float(sys.argv[1])
rMin = float(sys.argv[2])
rMax = float(sys.argv[3])

NCELLS = int(sys.argv[4])
IE = float(sys.argv[5])

dr = (rMax - rMin)/(NCELLS)

answer = np.zeros( (NCELLS,4) )

tsProblem = TS.TaylorSedov( 1.0, IE, 5.0/3.0 )

print("At t = %f, the shock front is at r = %.8f" % (Time,tsProblem.ShockR(Time)))
print("# r\trho\tu\tp")
for i in range(NCELLS):
    r = rMin + (dr)*(i + 0.5)
    answer[i,0] = r
    answer[i,1:4] = tsProblem.Solution( Time, r )
    print("%.8f\t%.8f\t%.8f\t%.8f" % tuple(answer[i,:]))



