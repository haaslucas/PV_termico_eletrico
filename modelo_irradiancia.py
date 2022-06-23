
# modelo de irradiância de Erico



#Cálculo da declinação solar:
dec=23.45*np.sin(np.radians(360*((284+n)/365)))

#Equação do tempo:
B=(n-1)*(360/365);
E=229.2*(0.0000750+0.001868*np.cos(np.radians(B))-0.032077*np.sin(np.radians(B))-0.014615*np.cos(np.radians(2*B))-0.04089*np.sin(np.radians(2*B)));

#Ângulo horário máximo ("sunset hour angle"):
ahm = np.degrees( np.arccos( (-np.tan(np.radians(la))) * (np.tan(np.radians(dec))) ))

#Obtenção das médias de irradiância para cada intervalo de 10 min no 
#intervalo de 6h à 17h50min:
c= np.arange(300,1071,10)       #Marcador dos intervalos de tempo
G=np.zeros(len(c))              #Irradiância total na horizontal (kW/m^2)
Gd=np.zeros(len(c))             #Irradiância difusa na horizontal (kW/m^2)
Gi=np.zeros(len(c))             #Irradiância total na superfícia inclinada (kW/m^2)
S=np.zeros(len(c))              #Irradiância total absorvida nas células fotovoltaicas (kW/m^2)

#Coeficientes a e b:
a=0.409+0.5016*np.sin(np.radians(ahm-60)); 
b=0.6609-0.4767*np.sin(np.radians(ahm-60));



#antigo for
#Superfície horizontal:
hs=(c+5)/60+(abs(lp)-abs(ll))/15+E/60-12;   #Hora solar
ah=hs*15;                               #Ângulo horário 
rt=np.zeros(78)  
rd=np.zeros(78)
G=np.zeros(78)
Gd=np.zeros(78)
Rb=np.zeros(78)
caz=np.zeros(78)
Ob=np.zeros(78)
Kbb=np.zeros(78)
for i in range(0,77):
    if ((ah[i]-1.25)>=-ahm) and ((ah[i]+1.25)<=ahm):
        #rt = Fração horária da irradiação diária
        rt[i]=(np.pi/144)*(a+b*np.cos(np.radians(ah[i])))*((np.cos(np.radians(ah[i]))-np.cos(np.radians(ahm)))/(np.sin(np.radians(ahm))-((np.pi*ahm)/(180))*np.cos(np.radians(ahm))))
        #rd = Fração horária da irradiação diária
        rd[i]=(np.pi/144)*((np.cos(np.radians(ah[i])))-np.cos(np.radians(ahm)))/(np.sin(np.radians(ahm))-((np.pi*ahm)/(180))*np.cos(np.radians(ahm))) 
        #G=irradiação total [W/m²]
        G[i]=(H*rt[i])/600
        #G=irradiação difusa [W/m²]
        Gd[i]=(Hd*rd[i])/600
        #Fator de elevação para módulo localizado no hemisfério sul com face voltada para o norte  
        Rb[i]=(np.cos(np.radians(la+im))*np.cos(np.radians(dec))*np.cos(np.radians(ah[i]))+np.sin(np.radians(la+im))*np.sin(np.radians(dec)))/(np.cos(np.radians(la))*np.cos(np.radians(dec))*np.cos(np.radians(ah[i]))+np.sin(np.radians(la))*np.sin(np.radians(dec))) 
        #Irradiância total na superfície inclinada em condição de céu isotrópico
        Gi[i]=(G[i]-Gd[i])*Rb[i]+Gd[i]*((1+np.cos(np.radians(im)))/2)+G[i]*ref*((1-np.cos(np.radians(im)))/2) 
        #Cosseno do ângulo zenital
        caz[i]=np.cos(np.radians(la))*np.cos(np.radians(dec))*np.cos(np.radians(ah[i]))+np.sin(np.radians(la))*np.sin(np.radians(dec))
        Ob[i]=np.degrees(np.arccos((Rb[i]*caz[i])))
        Kbb[i]=1+bo*((1/np.cos(np.radians(Ob[i])))-1)

Od=59.68-0.1388*im+0.001497*(im**2)
Ogr=90-0.5788*im+0.002693*(im**2)
Kdd=1+bo*((1/np.cos(np.radians(Od)))-1)
Kgg=1+bo*((1/np.cos(np.radians(Ogr)))-1)
        
S=ta*(Rb*(G-Gd)*Kbb+Gd*Kdd*((1+np.cos(np.radians(im)))/2)+G*ref*Kgg*((1-np.cos(np.radians(im)))/2));

Gi=1000*Gi; #Conversão de kW/m2 para W/m2
S=1000*S;
