from scipy.optimize import fsolve
import numpy as np

def R_shunt(Rs,Rp,Ns,Iscn,Ki,dT,Gi,Gn,Vocn,Kv,kb,q,Tn,Vt,Vmp,Imp,Pmaxe,Egap):
    '''
    calculate Rp for a given Rs and idl
    '''
    Ipvn = (Rs+Rp)/Rp * Iscn         # [A] Photogenerated current at STC
    idl = (Kv - Vocn/Tn) / ( Ns * Vt * ( Ki/Ipvn - 3/Tn - Egap/(kb*Tn**2) ) )
    Ipv = (Ipvn + Ki*dT) *Gi/Gn     # Photogenerated current
    Isc = (Iscn + Ki*dT) *Gi/Gn     # Short circuit current
    Voc = (Vocn + Kv*dT) *Gi/Gn*(1+(kb*Tn/q)*(np.log(Gi/Gn))/Vocn) # open circuit voltage
    Ion = (Ipv - Voc/Rp)/(np.exp(Voc/Vt/idl/Ns)-1)                 # Corrente de fuga
    Io = Ion                        # Corrente de fuga


    Rp = Vmp*(Vmp+Imp*Rs)/(Vmp*Ipv-Vmp*Io*np.exp((Vmp+Imp*Rs)/Vt/Ns/idl)+Vmp*Io-Pmaxe)

    return Rp,Voc,Ipv,Io,idl