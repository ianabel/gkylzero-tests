#!/usr/bin/env python3
import ctypes

class TaylorSedov(ctypes.Structure):
	_fields_ = [("rhoZero",ctypes.c_double),("InjectedEnergy",ctypes.c_double),("gas_gamma",ctypes.c_double),("alpha",ctypes.c_double)]
	def __init__(self,rhoZ,IE,gg):
		self.rhoZero = rhoZ
		self.InjectedEnergy = IE
		self.gas_gamma = gg
		self.lib = ctypes.cdll.LoadLibrary("/home/ian/projects/gkylzero-tests/TaylorSedovSolution.so")
		self.lib.SetAlpha(ctypes.byref(self))
	
	def ShockR(self,t):
		return self.lib.TaylorSedovR(ctypes.c_double(t),ctypes.byref(self))



testProblem = TaylorSedov(1.0,1.0,5.0/3.0)

print(" Alpha =  %.8f" % testProblem.alpha)
print(" Shock location at t = %f is %.8f" % (0.1,testProblem.ShockR(0.1)))
