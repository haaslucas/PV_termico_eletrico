import numpy as np
import time
t1 = time.time()
Tc=np.zeros(78)
Tv=np.zeros(78)
aRe=np.zeros(78)
aNu=np.zeros(78)
ahconv=np.zeros(78)


# seleção do mês:
# #janeiro, fevereiro, marco, abril, maio, junho, julho, agosto, setembro, outubro, novembro, dezembro
execfile('abril.py')



execfile('data_UPSOLAR335.py') #carrega os parâmetros do painel fotovoltaico
execfile('parametros_termofisicos.py')  
execfile('modelo_irradiancia.py')   # calcula as irradiancias Gi a cada 10m 
i = 0
execfile('eletrico_Tchute.py')

 # calcula a 1a efic, a partir de uma temperatura inicial (Tk)
 
while i < 78:
    if Gi[i] == 0:
        
        i += 1
    else:
        execfile('termico.py')
        while abs(Tk[i]-Tc[i]) > 0.005: 
            Tk[i] = Tc[i]
            execfile('eletrico.py') #calcula nova efic(1) ( G(1), Tc(1) )
            execfile('termico.py') # calcula nova Tc(i)
        i += 1
        Tk[i] = Tk[i-1]
        efic[i] = efic[i-1]
        Pmax[i] = Pmax[i-1]
    

tempoExec = time.time() - t1

execfile('plotTE.py')


print("Tempo de execução: {} segundos".format(tempoExec))