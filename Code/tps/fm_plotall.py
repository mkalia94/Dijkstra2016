from tps import *
def plotall(fm,t,y,fig,spec,case):
    if case == 1:
    #plt.rc('text', usetex=True)
        #plt.locator_params(axis='y', nbins=4)
        #plt.locator_params(axis='x', nbins=5)

        if 'excite' in fm.__dict__.keys():
            ax00 = fig.add_subplot(spec[0,0])
            ax00.spines['top'].set_visible(False)
            ax00.spines['bottom'].set_visible(False)
            ax00.spines['right'].set_visible(False)
            ax00.get_xaxis().set_ticks([])
            ax00.get_yaxis().set_ticks([0,fm.model(array(t),y,'max(p.F*IExcite)')])
            ax00.plot(t,fm.model(array(t),y,'p.F*IExcite'),color='black')
            ax00.xaxis.set_tick_params(labelsize=8)
            ax00.yaxis.set_tick_params(labelsize=8)
            tx1 = ax00.yaxis.get_offset_text()
            tx1.set_fontsize(2)
            ax00.set_xlim(fm.t0,fm.tfinal)
            ax00.set_ylim(-5+fm.model(array(t),y,'min(p.F*IExcite)'),5+fm.model(array(t),y,'max(p.F*IExcite)'))
            #plottwoaxes(fm,t,y,['blockerExp'],[],'',fig,ax00)
            ax01 = fig.add_subplot(spec[0,1])
            ax01.spines['top'].set_visible(False)
            ax01.spines['bottom'].set_visible(False)
            ax01.spines['right'].set_visible(False)
            ax01.get_xaxis().set_ticks([])
            ax01.get_yaxis().set_ticks([0,fm.model(array(t),y,'max(p.F*IExcite)')])
            ax01.plot(t,fm.model(array(t),y,'p.F*IExcite'),color='black')
            ax01.xaxis.set_tick_params(labelsize=8)
            ax01.yaxis.set_tick_params(labelsize=8)
            tx2 = ax00.yaxis.get_offset_text()
            tx2.set_fontsize(2)
            ax01.set_xlim(fm.t0,fm.tfinal)
            ax01.set_ylim(-5+fm.model(array(t),y,'min(p.F*IExcite)'),5+fm.model(array(t),y,'max(p.F*IExcite)'))
            #plottwoaxes(fm,t,y,['blockerExp'],[],'',fig,ax00)
        else:
            
            # blockerExp
            ax00 = fig.add_subplot(spec[0,0])
            ax00.spines['top'].set_visible(False)
            ax00.spines['bottom'].set_visible(False)
            ax00.spines['right'].set_visible(False)
            ax00.get_xaxis().set_ticks([])
            ax00.get_yaxis().set_ticks([0,100])
            ax00.plot(t,fm.model(array(t),y,'100*blockerExp'),color='black')
            ax00.plot([fm.t0,fm.tfinal],[fm.model(array(t),y,'100*min(blockerExp)'),fm.model(array(t),y,'100*min(blockerExp)')],'--',color='black',alpha=0.3)
            ax00.xaxis.set_tick_params(labelsize=8)
            ax00.yaxis.set_tick_params(labelsize=8)
            tx1 = ax00.yaxis.get_offset_text()
            tx1.set_fontsize(2)
            ax00.set_xlim(fm.t0,fm.tfinal)
            ax00.set_ylim(-5,105)
            #plottwoaxes(fm,t,y,['blockerExp'],[],'',fig,ax00)
            ax01 = fig.add_subplot(spec[0,1])
            ax01.spines['top'].set_visible(False)
            ax01.spines['bottom'].set_visible(False)
            ax01.spines['right'].set_visible(False)
            ax01.get_xaxis().set_ticks([])
            ax01.get_yaxis().set_ticks([0,100])
            ax01.plot(t,fm.model(array(t),y,'100*blockerExp'),color='black')
            ax01.plot([fm.t0,fm.tfinal],[fm.model(array(t),y,'100*min(blockerExp)'),fm.model(array(t),y,'100*min(blockerExp)')],'--',color='black',alpha=0.3)
            ax01.xaxis.set_tick_params(labelsize=8)
            ax01.yaxis.set_tick_params(labelsize=8)
            tx2 = ax00.yaxis.get_offset_text()
            tx2.set_fontsize(2)
            ax01.set_xlim(fm.t0,fm.tfinal)
            ax01.set_ylim(-5,105)
            #plottwoaxes(fm,t,y,['blockerExp'],[],'',fig,ax00)




        #plt.tick_params(axis='both', which='major', labelsize=2)
        #plt.tick_params(axis='both', which='minor', labelsize=2)

        ax10 = fig.add_subplot(spec[1,0])
        plottwoaxes(fm,t,y,['NaCi'],['NaCe+0*t'],r"[Na$^+$] (mM)",fig,ax10)

        ax10 = fig.add_subplot(spec[1,1])
        plottwoaxes(fm,t,y,['KCi'],['KCe+0*t'],r"[K$^+$] (mM)",fig,ax10)

        ax20 = fig.add_subplot(spec[2,0])
        plottwoaxes(fm,t,y,['ClCi'],['ClCe+0*t'],r"[Cl$^-$] (mM)",fig,ax20)


        ax31 = fig.add_subplot(spec[2,1])
        plottwoaxes(fm,t,y,['Voli'],[],r"Volume increase ($\%$)",fig,ax31)
        

        ax40 = fig.add_subplot(spec[3,0])
        plottwoaxes(fm,t,y,['Vi'],[],"Mem. Potential (mV)",fig,ax40)

        
        ax51 = fig.add_subplot(spec[3,1])
        plottwoaxes(fm,t,y,['Ipumpi'],[],"NKA current (pA)",fig,ax51)
        ax51.set_xlabel('time (min.)', fontdict={'fontsize': 8, 'fontweight': 'medium'})

    elif case == 2:

        if 'excite' in fm.__dict__.keys():
            ax00 = fig.add_subplot(spec[0,2])
            ax00.spines['top'].set_visible(False)
            ax00.spines['bottom'].set_visible(False)
            ax00.spines['right'].set_visible(False)
            ax00.get_xaxis().set_ticks([])
            ax00.get_yaxis().set_ticks([0,fm.model(array(t),y,'max(p.F*IExcite)')])
            ax00.plot(t,fm.model(array(t),y,'p.F*IExcite'),color='black')
            ax00.xaxis.set_tick_params(labelsize=8)
            ax00.yaxis.set_tick_params(labelsize=8)
            tx1 = ax00.yaxis.get_offset_text()
            tx1.set_fontsize(2)
            ax00.set_xlim(fm.t0,fm.tfinal)
            ax00.set_ylim(-5+fm.model(array(t),y,'min(p.F*IExcite)'),5+fm.model(array(t),y,'max(p.F*IExcite)'))
            #plottwoaxes(fm,t,y,['blockerExp'],[],'',fig,ax00)
            ax01 = fig.add_subplot(spec[0,3])
            ax01.spines['top'].set_visible(False)
            ax01.spines['bottom'].set_visible(False)
            ax01.spines['right'].set_visible(False)
            ax01.get_xaxis().set_ticks([])
            ax01.get_yaxis().set_ticks([0,fm.model(array(t),y,'max(p.F*IExcite)')])
            ax01.plot(t,fm.model(array(t),y,'p.F*IExcite'),color='black')
            ax01.xaxis.set_tick_params(labelsize=8)
            ax01.yaxis.set_tick_params(labelsize=8)
            tx2 = ax00.yaxis.get_offset_text()
            tx2.set_fontsize(2)
            ax01.set_xlim(fm.t0,fm.tfinal)
            ax01.set_ylim(-5+fm.model(array(t),y,'min(p.F*IExcite)'),5+fm.model(array(t),y,'max(p.F*IExcite)'))
            #plottwoaxes(fm,t,y,['blockerExp'],[],'',fig,ax00)
        else:
            
            # blockerExp
            ax00 = fig.add_subplot(spec[0,2])
            ax00.spines['top'].set_visible(False)
            ax00.spines['bottom'].set_visible(False)
            ax00.spines['right'].set_visible(False)
            ax00.get_xaxis().set_ticks([])
            ax00.get_yaxis().set_ticks([0,100])
            ax00.plot(t,fm.model(array(t),y,'100*blockerExp'),color='black')
            ax00.plot([fm.t0,fm.tfinal],[fm.model(array(t),y,'100*min(blockerExp)'),fm.model(array(t),y,'100*min(blockerExp)')],'--',color='black',alpha=0.3)
            ax00.xaxis.set_tick_params(labelsize=8)
            ax00.yaxis.set_tick_params(labelsize=8)
            tx1 = ax00.yaxis.get_offset_text()
            tx1.set_fontsize(2)
            ax00.set_xlim(fm.t0,fm.tfinal)
            ax00.set_ylim(-5,105)
            #plottwoaxes(fm,t,y,['blockerExp'],[],'',fig,ax00)
            ax01 = fig.add_subplot(spec[0,3])
            ax01.spines['top'].set_visible(False)
            ax01.spines['bottom'].set_visible(False)
            ax01.spines['right'].set_visible(False)
            ax01.get_xaxis().set_ticks([])
            ax01.get_yaxis().set_ticks([0,100])
            ax01.plot(t,fm.model(array(t),y,'100*blockerExp'),color='black')
            ax01.plot([fm.t0,fm.tfinal],[fm.model(array(t),y,'100*min(blockerExp)'),fm.model(array(t),y,'100*min(blockerExp)')],'--',color='black',alpha=0.3)
            ax01.xaxis.set_tick_params(labelsize=8)
            ax01.yaxis.set_tick_params(labelsize=8)
            tx2 = ax00.yaxis.get_offset_text()
            tx2.set_fontsize(2)
            ax01.set_xlim(fm.t0,fm.tfinal)
            ax01.set_ylim(-5,105)
            #plottwoaxes(fm,t,y,['blockerExp'],[],'',fig,ax00)





        #plt.tick_params(axis='both', which='major', labelsize=2)
        #plt.tick_params(axis='both', which='minor', labelsize=2)

        ax10 = fig.add_subplot(spec[1,2])
        plottwoaxes(fm,t,y,['NaCi'],['NaCe+0*t'],r"[Na$^+$] (mM)",fig,ax10)

        ax10 = fig.add_subplot(spec[1,3])
        plottwoaxes(fm,t,y,['KCi'],['KCe+0*t'],r"[K$^+$] (mM)",fig,ax10)

        ax20 = fig.add_subplot(spec[2,2])
        plottwoaxes(fm,t,y,['ClCi'],['ClCe+0*t'],r"[Cl$^-$] (mM)",fig,ax20)


        ax31 = fig.add_subplot(spec[2,3])
        plottwoaxes(fm,t,y,['Voli'],[],r"Volume increase ($\%$)",fig,ax31)
        

        ax40 = fig.add_subplot(spec[3,2])
        plottwoaxes(fm,t,y,['Vi'],[],"Mem. Potential (mV)",fig,ax40)


        ax51 = fig.add_subplot(spec[3,3])
        plottwoaxes(fm,t,y,['Ipumpi'],[],"NKA current (pA)",fig,ax51)
        ax51.set_xlabel('time (min.)', fontdict={'fontsize': 8, 'fontweight': 'medium'})

