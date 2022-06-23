# calcula a eficiência elétrica para uma Tc inicial.
from scipy.optimize import fsolve

Ipvn = (Rs+Rp)/Rp * Iscn  # Corrente fotogerada nas condições STC

#quando Gi=0 a eficiência não é calculada, portanto pula-se esse intervalo

while Gi[i] == 0:
    Gi[i] = float('nan')
    i = i + 1

Tk = np.zeros(78)
dT = np.zeros(78)
Ipv = np.zeros(78)
Isc = np.zeros(78)
Vocn = np.zeros(78)
P = np.zeros(78)
Pmax = np.zeros(78)
efic = np.zeros(78)
Ion = np.zeros(78)

Tk[i]=313 #Tk = Tc chute [K]
Vt  = kb*Tk[i]/q         # tensão térmica do diodo [V]
dT[i] = Tk[i]-Tn # Tk [K] é a temperatura da celula fotovoltaica temporária, Tn [K] é a temperatura STC 298.15K


Ipv[i] = (Ipvn + Ki*dT[i]) *Gi[i] / Gn # corrente fotogerada [A]
Isc[i] = (Iscn + Ki*dT[i]) *Gi[i]/Gn # corrente  de curto circuito [A]    
Vocn[i] = (Vocnn + Kv*dT[i]) *Gi[i]/Gn*(1+(kb*Tn/q)*(np.log(Gi[i]/Gn))/Vocnn) #tensão de circuito aberto [V]
Ion[i] = (Ipv[i] - Vocn[i]/Rp)/(np.exp(Vocn[i]/Vt/idl/Ns)-1) #corrente de fuga [A]
Io = Ion[i] #corrente de fuga [A]


j=1;
m=1;
Vv = np.zeros((78,800)) #Vv são as tensões da curva I-V característica. O cálculo da curva I-V começa com Vv=0.
I = np.zeros((78,800))
P = np.zeros((78,800))

#Vv[i,0]=3 
def f(cur):
    y = ( -cur+Ipv[i]-Io*(np.exp((Vv[i,m]+Rs*cur)/(Vt*Ns*idl))-1)-(Vv[i,m]+Rs*cur)/Rp )
    return y

## calcula as correntes I pras tensões Vv entre 0<Vv<Vocn. 
# "i" é o contador dos 78 intervalos de tempo
# "m" é o contador de quantas tensões há entre 0 e Vocn 
while Vv[i,m] < Vocn[i]:
    
    I[i,m] = fsolve(f,1)
    #método iterativo que calcula a correte I
    P[i,m] = (Ipv[i]-Io*(np.exp((Vv[i,m]+I[i,m]*Rs)/Vt/Ns/idl)-1)-(Vv[i,m]+I[i,m]*Rs)/Rp)*Vv[i,m] #calcula a potencia
    m=m+1
    Vv[i,m]= Vv[i,m-1] + 0.1 #incrementa Vv, até chegar a Vocn


    
Pmax[i] = max(P[i,:]); #Potência máxima da curva I-V
efic[i] = Pmax[i]/(area*Gi[i]); #eficiência elétrica
