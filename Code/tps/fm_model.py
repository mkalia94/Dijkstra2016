from tps import *
def model(t, y, p, *args):
    if size(shape(y)) == 2:
        NNa = y[:, 0]
        NK = y[:, 1]
        NCl = y[:, 2]
        m = y[:, 3]
        h = y[:, 4]
        n = y[:, 5]
        Vtemp = y[:, 6]
        Wi = y[:, 7]
    else:
        NNa = y[0]
        NK = y[1]
        NCl = y[2]
        m = y[3]
        h = y[4]
        n = y[5]
        Vtemp = y[6]
        Wi = y[7]

    # Ionic amounts and concentrations
    # ECS
    NaCe = p.NaCe0
    KCe = p.KCe0
    ClCe = p.ClCe0
    # Neuron
    NaCi = NNa/Wi
    KCi = NK/Wi
    ClCi = NCl/Wi
    
    # Voltages
    if 'excite' in p.__dict__.keys():
        V = Vtemp
    else:
        V = p.F/p.C*(NNa+NK-NCl-p.NAi)
    Vi = V # Needed for plotting

    # ==========================================================================
    # --------------------NEURON-------------------------------------------------
    # ===========================================================================

    alpham = 0.32*(V+52)/(1-exp(-(V+52)/4))
    betam = 0.28*(V+25)/(exp((V+25)/5)-1)
    alphah = 0.128*exp(-(V+53)/18)
    betah = 4/(1+exp(-(V+30)/5))
    alphan = 0.016*(V+35)/(1-exp(-(V+35)/5))
    betan = 0.25*exp(-(V+50)/40)

    if p.nogates:
        gates_block = 0
        m = alpham/(alpham + betam)
        h = alphah/(alphah + betah)
        n = alphan/(alphan + betan)
    else:
        gates_block = 1   

    # Gated currents
    INaG = p.PNaG*(m**3)*(h)*(p.F**2)*(V)/(
       p.R*p.T)*((NaCi -
                  NaCe*exp(-(p.F*V)/(p.R*p.T)))/(
                     1-exp(-(p.F*V)/(p.R*p.T))))
    IKG = (p.PKG*(n**2))*(p.F**2)*(V)/(
       p.R*p.T)*((KCi -
                  KCe*exp(-(p.F*V)/(p.R*p.T)))/(
                     1-exp(-p.F*V/(p.R*p.T))))
    IClG = p.PClG*1/(1+exp(-(V+10)/10))*(
       p.F**2)*V/(p.R*p.T)*((ClCi -
                             ClCe*exp(p.F*V/(p.R*p.T)))/(
                                1-exp(p.F*V/(p.R*p.T))))

    # Leak currents
    INaL = p.PNaL*(p.F**2)/(
       p.R*p.T)*V*((NaCi -
                    NaCe*exp((-p.F*V)/(p.R*p.T)))/(
                       1-exp((-p.F*V)/(p.R*p.T))))
    IKL = p.PKL*p.F**2/(p.R*p.T)*V*((
       KCi -
       KCe*exp((-p.F*V)/(p.R*p.T)))/(
          1-exp((-p.F*V)/(p.R*p.T))))
    IClL = p.PClL*(p.F**2)/(
       p.R*p.T)*V*((ClCi -
                    ClCe*exp((p.F*V)/(p.R*p.T)))/(
                       1-exp((p.F*V)/(p.R*p.T))))

    # Blockade
    blockerExp = 1/(1+exp(p.beta1*(t-p.tstart))) + 1/(
       1+exp(-p.beta2*(t-p.tend)))
    blockerExpAlt = 1/(1+exp(p.beta1*(t-p.tstart))) + p.perc/(
       1+exp(-p.beta2*(t-p.tend)))
    blockerExp = p.perc + (1-p.perc)*blockerExp
    #blockerExp = blockerExpAlt

    # Na-K pump
    sigmapump = 1/7*(exp(NaCe/67.3)-1)
    fpump = 1/(1+0.1245*exp(-0.1*p.F/p.R/p.T*V) +
               0.0365*sigmapump*exp(-p.F/p.R/p.T*V))
    Ipump = blockerExp*p.pumpScaleNeuron*fpump*p.Qpump*(
                NaCi**(1.5)/(NaCi**(1.5)+p.nka_na**1.5))*(KCe/(KCe+p.nka_k))

    # KCl cotransport
    JKCl = p.UKCl*p.R*p.T/p.F*(log(KCi)+log(ClCi)-log(KCe)-log(ClCe))



    # =========================================================================
    # ----------------------------VOLUME DYNAMICS------------------------------
    # =========================================================================
    SCi = NaCi+KCi+ClCi+p.NAi/Wi
    SCe = NaCe+KCe+ClCe+p.ACe
    delpii = p.R*p.T*(SCi-SCe)
    fluxi = p.LH20i*(delpii)
    Voli = Wi/p.Wi0*100

    # ==========================================================================
    # ----------------------------INTERVENTIONS---------------------------------
    # ==========================================================================

    if 'block' in p.__dict__.keys():
        dict_ = p.block
        for key in dict_:
            val_ = dict_[key]
            blockOther = (1/(1+exp(100*(t-val_[0]))) +
                          1/(1+exp(-100*(t-val_[1]))))
            if key == 'INaG':
                INaG = INaG*blockOther
            elif key == 'IKG':
                IKG = IKG*blockOther
            elif key == 'IClG':
                IClG = IClG*blockOther
            elif key == 'JKCl':
                JKCl = JKCl*blockOther
            elif key == 'WaterN':
                fluxi = fluxi*blockOther
    if 'excite' in p.__dict__.keys():
        arg_excite = p.excite
        blocker_Excite = 1 - (1/(1+exp(100*(t-arg_excite[0]))) +
                              1/(1+exp(-100*(t-arg_excite[1]))))
        IExcite = blocker_Excite*arg_excite[2]/2/p.F*(1-signal.square(array(5*t),duty=arg_excite[3]))
        dur_ = arg_excite[3]
        duty_ = arg_excite[4]
        IExcite = blocker_Excite*arg_excite[2]/2/p.F*(1-signal.square(2*pi*array(t)*(1-duty_)/(dur_/60),duty=duty_))
        #IExcite = blocker_Excite*4.5/p.F
    else:
        IExcite = 0
    # ==========================================================================
    # ----------------------------FINAL MODEL-----------------------------------
    # ==========================================================================
    Ipumpi = Ipump
    
    ODEs = [  # Neuron
       -1/p.F*(INaG+INaL+3*Ipump),
       (-1/p.F*(IKG+IKL-2*Ipump)-JKCl),
       (1/p.F*(IClG+IClL)-JKCl),
       gates_block*(alpham*(1-m)-betam*m),
       gates_block*(alphah*(1-h)-betah*h),
       gates_block*(alphan*(1-n)-betan*n),
       0,
       # WATER
       fluxi]

    if 'excite' in p.__dict__.keys():
        ODEs[6] =  p.F/p.C*(ODEs[0]+ODEs[1]-ODEs[2] + IExcite)
    
    ODEs = array(ODEs)*60*1e3

    if args:
        if 'ax1' in args[0] or 'ax2' in args[0]:
            temp_ = args[0]
            return eval(temp_[3:])
        else:     
            return eval(args[0])
    else:
        return ODEs
