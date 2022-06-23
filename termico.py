# Tchute inicial

#Tceu = 0.0552*Ta**1.5; # equação alternativa pra achar Tceu, não a usei.

Tceu = Ta -5
r=2 #contador
Tcc=np.zeros(30)
Tvv=np.zeros(30)
dif1=np.zeros(30)
dif2=np.zeros(30)
Tf=np.zeros(30)
K=np.zeros(30)
Pr=np.zeros(30)
v=np.zeros(30)
Re=np.zeros(30)
Nu=np.zeros(30)
hconv=np.zeros(30)
hrad=np.zeros(30)
Rrad=np.zeros(30)
Rconv=np.zeros(30)
Tcc[r] = Ta #1º chute: Tc temporária = temperatura do ar
Tvv[r] = Ta #1º chute: Tvidro temporária = temperatura do ar

# condutividade, Prandtl e viscosidade cinematica a 300K e 350K, a serem
# interpolados (Incropera tab. A4)
Ti = 300   #Temepratura = 300K
Tb = 350   #Temepratura = 350K
Ka = 26.3  #Condutividade ar a 300K
Kb = 30    #Condutividade ar a 350K
Pa = 0.707 #Prandtl a 300K
Pb = 0.7   #Prandtl a 350K
va = 15.89 #Viscosidade cinematica a 300K
vb = 20.92 #Viscosidade cinematica a 350K

f1=np.zeros(30)
f2=np.zeros(30)
erro=0.01

#valores iniciais pras variáveis não darem erro
dif1[1]=1
dif2[1]=1
dif1[2]=1
dif2[2]=1

#  método de Newton-Raphson para calcular Tv e Tc, metodologia do Chapra
while (abs(dif1[r-1])>=erro) or (abs(dif2[r-1])>=erro):
    Tf[r] = (Tvv[r]+Ta)/2                      #temperatura de filme
    K[r] = -10**-3*((Ka-Kb)*(Ti-Tf[r])/(Ti-Tb)-Ka) #interpolação pra r a condutividade ar
    Pr[r] = -1*((Pa-Pb)*(Ti-Tf[r])/(Ti-Tb)-Pa)    #interpolação pra achar o Prandtl
    v[r] = -10**-6*((va-vb)*(Ti-Tf[r])/(Ti-Tb)-va) #interpolação pra achar a viscosidade cinematica

    #convecção
    Re[r] = V*L/v[r] #numero de Reynolds
    Nu[r]=0.86*Re[r]**(0.5)*Pr[r]**(1/3) #numero de Nusselt
    hconv[r] = Nu[r]*K[r]/L #coeficiente convectivo
    
    #radiação
    hrad[r]= ev*sig*(Tvv[r]**2+Tceu**2)*(Tvv[r]+Tceu) #coeficiente radiativo
    
    #resistencias termicas do painel fotovoltaico 
    Rmodsup = Lv/kv +Leva1/keva
    Rmodinf = Lt/kt +Leva2/keva 
    Rrad[r] = 1/hrad[r]
    Rconv[r] = 1/hconv[r] 
    
        
    ## balanços em torno de Tv e Tc
    f1[r] = (Ta-Tvv[r])/Rconv[r] + (Tceu-Tvv[r])/Rrad[r] + afg*Gi[i] + (Tcc[r]-Tvv[r])/Rmodsup
    f2[r] = (Tvv[r]-Tcc[r])/Rmodsup + (Ta-Tcc[r])/(Rconv[r]+Rmodinf) + S[i] - Pmax[i]/area
    
    ## matriz jacobiana
    df1dTc = 1/Rmodsup                       #derivada de f1 em reação a Tc
    df1dTv = -1/Rconv[r] -1/Rrad[r] -1/Rmodsup     #derivada de f1 em reação a Tv 
    df2dTc = -1/Rmodsup -1/(Rconv[r]+Rmodinf)   #derivada de f2 em reação a Tc
    df2dTv = 1/Rmodsup                       #derivada de f2 em reação a Tv
    
    detJ = df1dTc*df2dTv - df1dTv*df2dTc     #determinante jacobiana
    
    Tcc[r+1] = Tcc[r] - (f1[r]*df2dTv-f2[r]*df1dTv)/detJ #Tcc é a temperatura Tc temporaria
    Tvv[r+1] = Tvv[r] - (f2[r]*df1dTc-f1[r]*df2dTc)/detJ #Tvv é a temperatura Tv temporaria
    
    ## estima a convergência
    dif1[r]=Tcc[r+1]-Tcc[r]
    dif2[r]=Tvv[r+1]-Tvv[r]
    r=r+1



Tc[i]=Tcc[r];
Tv[i]=Tvv[r];
aRe[i]=Re[r-1];
aNu[i]=Nu[r-1];
ahconv[i]=hconv[r-1];
