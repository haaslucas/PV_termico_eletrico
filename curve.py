from scipy.optimize import fsolve
import numpy as np

def curve (Voc,idl,Rp,Rs,Vt,Io,Ns,Ipv):
    '''
    generates the V-I characteristic curve for the inputs Rs, Rp and idl
    '''
    i=0
    V=[0]
    I=[]
    P=[]
    # V are the voltages between 0 and Voc.

    while V[i] <= Voc:
        V.append(V[i]+0.1)
        i += 1

    def f(cur):
        y = ( -cur+Ipv-Io*(np.exp((V[i]+Rs*cur)/(Vt*Ns*idl))-1)-(V[i]+Rs*cur)/Rp )
        return y


    ## calcula as correntes I pras tensões Vv entre 0<Vv<Vocn. 
    # "i" é o contador dos intervalos de tempo
    # "m" é o contador de quantas tensões há entre 0 e Vocn 
    i=0
    for volt in V:
        #incrementa Vv, até chegar a Vocn
        I.append(float((fsolve(f,1))))
        #método iterativo que calcula a correte I
        P.append(V[i]*I[i])  #[i] = (Ipv-Io*(np.exp((V[i]+I[i]*Rs)/Vt/Ns/idl)-1)-(V[i]+I[i]*Rs)/Rp)*V[i] #calcula a potencia
        i += 1
    Pmax = max(P) #Potência máxima da curva I-V
    return V,I,P,Pmax