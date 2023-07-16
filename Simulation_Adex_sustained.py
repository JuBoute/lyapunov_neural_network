#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: jboute
"""


import matplotlib.pyplot as plt
import numpy as np
from brian2 import *
import scipy.fftpack
from time import time
import pickle



#for NbS in range(0,1):
NbSim=5
N1=2000
N2=8000
prbC=0.05*2000/N1

print('simulation #'+str(NbSim))

def bin_array(array, BIN, time_array):
    N0 = int(BIN/(time_array[1]-time_array[0]))
    N1 = int((time_array[-1]-time_array[0])/BIN)
    return array[:N0*N1].reshape((N1,N0)).mean(axis=1)


def multi_Standard(input_freq,t_stim,N1,N2,prbC,NbSims):
    tic=time()
    print('New simulation, input freq = ', input_freq,'Hz')
    start_scope()

    ResultsTotal=[]
    

    for sim in NbSims:
        print('Simulation #{} out of {}'.format(sim+1,NbSim))
        DT=0.1
        defaultclock.dt = DT*ms
        
        TotTime=20000
        duration = TotTime*ms
        
        seed(sim)
        
        eqs='''
        dv/dt = (-GsynE*(v-Ee)-GsynI*(v-Ei)-gl*(v-El)+ gl*Dt*exp((v-Vt)/Dt)-w + Is)/Cm : volt (unless refractory)
        dw/dt = (a*(v-El)-w)/tau_w:ampere
        dGsynI/dt = -GsynI/Tsyn : siemens
        dGsynE/dt = -GsynE/Tsyn : siemens
        Is:ampere
        Cm:farad
        gl:siemens
        El:volt
        a:siemens
        tau_w:second
        Dt:volt
        Vt:volt
        Ee:volt
        Ei:volt
        Tsyn:second
        '''
        
        # Population 1 - FS
        b1 = 0.0*pA
        G1 = NeuronGroup(N1, eqs, threshold='v > -47.5*mV', reset='v = -65*mV', refractory='5*ms', method='heun')
        #init:
        G1.v = -55*mV
        G1.w = 0.0*pA
        G1.GsynI=0.0*nS
        G1.GsynE=0.0*nS
        #parameters
        G1.Cm = 200.*pF
        G1.gl = 8.*nS
        G1.El = -60.*mV
        G1.Vt = -50.*mV
        G1.Dt = 0.5*mV
        G1.tau_w = 1.0*ms
        G1.a = 0.0*nS
        G1.Is = 0.0
        
        G1.Ee=0.*mV
        G1.Ei=-80.*mV
        G1.Tsyn=5.*ms
        
        # Population 2 - RS
        b2 = 60.*pA
        G2 = NeuronGroup(N2, eqs, threshold='v > -40.0*mV', reset='v = -65*mV; w += b2', refractory='5*ms',  method='heun')
        G2.v = -42.*mV
        G2.w = 0.0*pA
        G2.GsynI=0.0*nS
        G2.GsynE=0.0*nS
        G2.Cm = 200.*pF
        G2.gl = 8.*nS
        G2.El = -60.*mV
        G2.Vt = -52.*mV
        G2.Dt = 2.*mV
        G2.tau_w = 100.*ms
        G2.a = 0.*nS
        G2.Is = 0.0*nA 
        
        G2.Ee=0.*mV
        G2.Ei=-80.*mV
        G2.Tsyn=5.*ms

        
        #Create neurons that will control metrics online
        
        W_tot = NeuronGroup(1,'W_tot : 1'  , method='heun') 

        # external drive--------------------------------------------------------------------------

        AmpStim=0. #92.
        plat = 1000
        def heaviside(x):
            return 0.5 * (1 + np.sign(x))
          
        t2 = np.arange(0, TotTime, DT)
        test_input = []
          
        for ji in t2:
            # print(i), print (input_rate(i))
            test_input.append(input_freq*(heaviside(ji)-heaviside(ji-t_stim)))
          
        stimulus=TimedArray(test_input*Hz, dt=DT*ms)
        #  print(max(test_input))
       
     
     
        P_ed=PoissonGroup(8000, rates='stimulus(t)')
            
        
        # connections-----------------------------------------------------------------------------
        
        Qi=5.0*nS
        Qe=1.5*nS
        
        prbC2=0.05
        
        S_12 = Synapses(G1, G2, on_pre='GsynI_post+=Qi') #'v_post -= 1.*mV')
        S_12.connect('i!=j', p=prbC2)
        
        S_11 = Synapses(G1, G1, on_pre='GsynI_post+=Qi')
        S_11.connect('i!=j',p=prbC2)
        
        S_21 = Synapses(G2, G1, on_pre='GsynE_post+=Qe')
        S_21.connect('i!=j',p=prbC)
        
        S_22 = Synapses(G2, G2, on_pre='GsynE_post+=Qe')
        S_22.connect('i!=j', p=prbC)
        
        
        
        S_ed_in = Synapses(P_ed, G1, on_pre='GsynE_post+=Qe')
        S_ed_in.connect(p=prbC2)
        
        S_ed_ex = Synapses(P_ed, G2, on_pre='GsynE_post+=Qe')
        S_ed_ex.connect(p=prbC2)#0.05)
        
        #Connect the control neurons to the populations. To do so, apply the interesting function as the post variable.

        S_W=Synapses(G2, W_tot, 'W_tot_post = w_pre : 1 (summed)')
        S_W.connect(p=1)

        
       
        ####################################################################################
        
                    ####        RECORDER GROUPS     ######
                    
        ####################################################################################
        
        
        dt_rec=1*ms
        
        M1G1 = SpikeMonitor(G1)
        FRG1 = PopulationRateMonitor(G1)
        
        M1G2 = SpikeMonitor(G2)
        FRG2 = PopulationRateMonitor(G2)

        
        #Control neurons
        
        MonW_tot=StateMonitor(W_tot, 'W_tot', record=0)
        
        print('--##Start simulation##--')
        run(duration)
        print('--##End simulation##--')
        tac=time()
        print('Simulation took {}s'.format(tac-tic))
        
       
        
        LfrG1=np.array(FRG1.smooth_rate(window='gaussian', width=500*ms)/Hz)
        
        LfrG2=np.array(FRG2.smooth_rate(window='gaussian', width=500*ms)/Hz)

        WG2_tot=MonW_tot.W_tot[0]
        
        
        #Save all
        
        FR_input=[LfrG1,LfrG2]
     
        ResultsTotal.append([FR_input,WG2_tot])
    
        fig,ax=plt.subplots()
        ax.plot(LfrG1)
        ax.plot(LfrG2)
                
    
    with open('Simulation_Adex_sustained_1',"wb") as f:
        pickle.dump(ResultsTotal,f)
    


input_freq=1.5
NbSim=[1]#1,3
multi_Standard(input_freq,100,N1,N2,prbC,NbSim)

