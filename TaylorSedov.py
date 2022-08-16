#!/usr/bin/env python3
import ctypes

class TaylorSedov(ctypes.Structure):
    _fields_ = [("rhoZero",ctypes.c_double),("InjectedEnergy",ctypes.c_double),("gas_gamma",ctypes.c_double),("alpha",ctypes.c_double),("dimension",ctypes.c_uint)]
    def __init__(self,rhoZ,IE,gg):
        self.rhoZero = rhoZ
        self.InjectedEnergy = IE
        self.gas_gamma = gg
        self.dimension = 3
        self.lib = ctypes.cdll.LoadLibrary("/home/ian/projects/gkylzero-tests/TaylorSedovSolution.so")
        self.lib.SetAlpha(ctypes.byref(self))
        self.lib.uTilde.restype = ctypes.c_double
        self.lib.pTilde.restype = ctypes.c_double
        self.lib.RhoTilde.restype = ctypes.c_double

        self.lib.TaylorSedovR.restype = ctypes.c_double
        self.lib.TaylorSedovU.restype = ctypes.c_double
        self.lib.TaylorSedovRho.restype = ctypes.c_double
        self.lib.TaylorSedovP.restype = ctypes.c_double

        self.lib.TaylorSedovR.argtypes = [ctypes.c_double,ctypes.c_void_p]
        self.lib.TaylorSedovU.argtypes = [ctypes.c_double,ctypes.c_double,ctypes.c_void_p]
        self.lib.TaylorSedovRho.argtypes = [ctypes.c_double,ctypes.c_double,ctypes.c_void_p]
        self.lib.TaylorSedovP.argtypes = [ctypes.c_double,ctypes.c_double,ctypes.c_void_p]

        self.lib.uTilde.argtypes = [ctypes.c_double,ctypes.c_void_p]
        self.lib.pTilde.argtypes = [ctypes.c_double,ctypes.c_void_p]
        self.lib.RhoTilde.argtypes = [ctypes.c_double,ctypes.c_void_p]

    def ShockR(self,t):
        return self.lib.TaylorSedovR(t,ctypes.byref(self))

    def uTilde(self,xi):
        return self.lib.uTilde(xi,ctypes.byref(self))

    def pTilde(self,xi):
        return self.lib.pTilde(xi,ctypes.byref(self))

    def rhoTilde(self,xi):
        return self.lib.RhoTilde(xi,ctypes.byref(self))

    def TaylorSedovRho( self, t, r ):
        return self.lib.TaylorSedovRho( t, r, ctypes.byref(self) )

    def TaylorSedovP( self, t, r ):
        return self.lib.TaylorSedovP( t, r, ctypes.byref(self) )

    def TaylorSedovU( self, t, r ):
        return self.lib.TaylorSedovU( t, r, ctypes.byref(self) )

    def Solution(self,t,r):
        return [self.TaylorSedovRho(t,r),self.TaylorSedovU(t,r),self.TaylorSedovP(t,r)]
