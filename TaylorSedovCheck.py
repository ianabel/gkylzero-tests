#!/usr/bin/env python3

import TaylorSedov as TS
import numpy

testProblem = TS.TaylorSedov(1.0,1.0,5.0/3.0)

print("Checking against stored data")

storedData = numpy.loadtxt('NormalisedData.dat')

NX = storedData.shape[0]
dXi = 1.0/(NX-1)

calculatedData = numpy.zeros(storedData.shape)

print("At xi = 1, we have, rho = %.8f, p = %.8f, u = %.8f" % (testProblem.rhoTilde(1.0),testProblem.pTilde(1.0),testProblem.uTilde(1.0)))

for i in range(NX):
    Xi = i * dXi
    calculatedData[i,0] = Xi
    calculatedData[i,1] = testProblem.pTilde( Xi ) / testProblem.pTilde( 1 )
    calculatedData[i,2] = testProblem.rhoTilde( Xi ) / testProblem.rhoTilde( 1 )
    calculatedData[i,3] = testProblem.uTilde( Xi ) / testProblem.uTilde( 1 )

print( "| p - p* |_Inf = %.8f" % numpy.linalg.norm(storedData[:,1] - calculatedData[:,1],numpy.inf))
print( "| u - u* |_Inf = %.8f" % numpy.linalg.norm(storedData[:,3] - calculatedData[:,3],numpy.inf))
print( "| rho - rho* |_Inf = %.8f" % numpy.linalg.norm(storedData[:,2] - calculatedData[:,2],numpy.inf))
