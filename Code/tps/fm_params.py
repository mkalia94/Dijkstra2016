from tps import *


def parameters(p, dict_):
    p.__dict__.update(dict_)

    p.NNai0 = p.NaCi0*p.Wi0
    p.NKi0 = p.KCi0*p.Wi0
    p.NCli0 = p.ClCi0*p.Wi0


    
    # Impermeants and conserved quantities
    p.NAi = (- p.NCli0 + p.NKi0 + p.NNai0 - (p.C*p.Vi0)/p.F)
    p.ACe = (p.NaCi0+p.KCi0+p.ClCi0+p.NAi/p.Wi0 - p.NaCe0 - p.KCe0 - p.ClCe0)

    # Gates    
    p.alpham0 = 0.32*(p.Vi0+52)/(1-exp(-(p.Vi0+52)/4))
    p.betam0 = 0.28*(p.Vi0+25)/(exp((p.Vi0+25)/5)-1)
    p.alphah0 = 0.128*exp(-(p.Vi0+53)/18)
    p.betah0 = 4/(1+exp(-(p.Vi0+30)/5))
    p.alphan0 = 0.016*(p.Vi0+35)/(1-exp(-(p.Vi0+35)/5))
    p.betan0 = 0.25*exp(-(p.Vi0+50)/40)
    p.m0 = p.alpham0/(p.alpham0+p.betam0)
    p.h0 = p.alphah0/(p.alphah0+p.betah0)
    p.n0 = p.alphan0/(p.alphan0+p.betan0)
    
    # Neuronal leaks
    p.INaG0 = p.PNaG*(p.m0**3)*(p.h0)*(p.F**2)*(p.Vi0)/(
        p.R*p.T)*((p.NaCi0 -
                   p.NaCe0*exp(-(p.F*p.Vi0)/(p.R*p.T)))/(
                       1-exp(-(p.F*p.Vi0)/(p.R*p.T))))
    p.IKG0 = (p.PKG*(p.n0**2))*(p.F**2)*(p.Vi0)/(
        p.R*p.T)*((p.KCi0 -
                   p.KCe0*exp(-(p.F*p.Vi0)/(p.R*p.T)))/(
                       1-exp(-p.F*p.Vi0/(p.R*p.T))))
    
    p.IClG0 = p.PClG*1/(1+exp(-(p.Vi0+10)/10))*(p.F**2)*p.Vi0/(
        p.R*p.T)*((p.ClCi0 -
                   p.ClCe0*exp(p.F*p.Vi0/(p.R*p.T)))/(
                       1-exp(p.F*p.Vi0/(p.R*p.T))))
    p.INaL0 = (p.F**2)/(p.R*p.T)*p.Vi0*((
        p.NaCi0 -
        p.NaCe0*exp((-p.F*p.Vi0)/(p.R*p.T)))/(
            1-exp((-p.F*p.Vi0)/(p.R*p.T))))
    p.IKL0 = p.F**2/(p.R*p.T)*p.Vi0*((
        p.KCi0-p.KCe0*exp((-p.F*p.Vi0)/(p.R*p.T)))/(
            1-exp((-p.F*p.Vi0)/(p.R*p.T))))
    p.IClL0 = (p.F**2)/(p.R*p.T)*p.Vi0*((
        p.ClCi0-p.ClCe0*exp((p.F*p.Vi0)/(p.R*p.T)))/(
            1-exp((p.F*p.Vi0)/(p.R*p.T))))
    p.JKCl0 = p.UKCl*p.R*p.T/p.F*(log(p.KCi0) +
                                  log(p.ClCi0)-log(p.KCe0)-log(p.ClCe0))
    p.sigmapump = 1/7*(exp(p.NaCe0/67.3)-1)
    p.fpump = 1/(1+0.1245*exp(-0.1*p.F/p.R/p.T*p.Vi0) +
                 0.0365*p.sigmapump*exp(-p.F/p.R/p.T*p.Vi0))
    p.neurPump = p.pumpScaleNeuron*p.Qpump*p.fpump*(p.NaCi0**(1.5)/(
        p.NaCi0**(1.5)+p.nka_na**1.5))*(p.KCe0/(p.KCe0+p.nka_k))

    p.PNaL = (-p.INaG0 - 3*p.neurPump)/p.INaL0  # Estimated sodium leak
    #                                                    conductance in neuron
    p.PKL = (-p.IKG0 + 2*p.neurPump - p.F*p.JKCl0)/p.IKL0    # Estimated K leak conductance in neuron 
    p.PClL = (p.F*p.JKCl0 - p.IClG0)/p.IClL0                 # Estimated Cl leak conducatance in neuron
    

    #================================================================================
    #    CHECKS CHECKS CHECKS CHECKS CHECKS
    #=============================================================================
    
    ctr = 0
    p.err = 0
    
    if (p.NAi>0) & (p.ACe>0):
        disp('Quantity of impermeants....OK')
    else:
        disp('ERROR: Quantity of an impermeant is nonpositive')
        p.err = 10
    
    if (p.PNaL>0) & (p.PKL>0) & (p.PClL>0):
        disp('Leak cond. in neurons....OK')
    else:
        disp('ERROR: Sign error in neuron leak conductance')
        disp('PNaL (>0): {}'.format(p.PNaL))
        disp('PKL (>0): {}'.format(p.PKL))
        disp('PClL (>0): {}'.format(p.PClL))
        disp('----------------------------------')
        p.err = 10
