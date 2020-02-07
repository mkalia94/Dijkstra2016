from tps import *
def get_initvals(fm):
    fm.initvals = array([fm.NNai0, fm.NKi0, fm.NCli0, fm.m0, fm.h0, fm.n0, fm.Vi0,fm.Wi0])
    if fm.nogates:
        fm.initvals = fm.initvals + 1e-8*array([1,1,1,0,0,0,1,1])*fm.initvals
    else:    
        fm.initvals  = fm.initvals + 1e-5*fm.initvals
    
