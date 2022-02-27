import numpy as np
import matplotlib.pyplot as plt
def blockAverage(datos,tamaño_bloque_maximo=False,grafica=True):
    "Método para calcular incertidumes con datos correlacionados temporalmente, block-avarage. Sería ideal engadir o método de bootstraping para comparar"
    Nobs = len(datos) # número total de datos, dar un array ou lista
    tamaño_bloque_minimo=1 # mínimo tamaño de bloque, corresponde a non facer block_avarage
    if tamaño_bloque_maximo==False:
        tamaño_bloque_maximo=int(Nobs/4)#Criterio por se non se selecciona un dato óptimo de bloque
    Numero_de_bloques = tamaño_bloque_maximo-tamaño_bloque_minimo # total number of block sizes
    Media_bloque = np.zeros(Numero_de_bloques) 
    co = np.zeros(Numero_de_bloques) # definese no seguinte papper, pero é a varianza
    print('A información do método atópase en:  https://doi.org/10.1063/1.457480')
    bloqueCtr=0
    incertidume=np.zeros(Numero_de_bloques)
    barra_erros=np.zeros(Numero_de_bloques)
    for datos_bloque in range(tamaño_bloque_minimo,tamaño_bloque_maximo):
        Nbloque=int(Nobs/datos_bloque) # En que bloque estou según o numero de observacions
        almacenamento_temporal=np.zeros(Nbloque) # isto é para almacenar a media según en que bloque esté, é temporal
        # Loop to chop datastream into blocks
        # and take average
        for i in range(1,Nbloque+1):
            #Esta é a parte máis complicada do código pois (fago de 1 a N+1, para empregar logo as definicións matemáticas)
            #E como vou calculando as medias según o bloque no que me atope, e dicir o número de datos por bloque
            comeza = (i-1)*datos_bloque
            remata = comeza+datos_bloque
            almacenamento_temporal[i-1] = np.mean(datos[comeza:remata])
        Media_bloque[ bloqueCtr] = np.mean(almacenamento_temporal)
        co[bloqueCtr]  = np.var(almacenamento_temporal)
        incertidume[ bloqueCtr]=np.sqrt(np.var(almacenamento_temporal)/((Nbloque-1)))
        barra_erros[ bloqueCtr]=np.sqrt(np.var(almacenamento_temporal)/(Nbloque-1))*1/(np.sqrt(2*(Nbloque-1)))
        bloqueCtr=bloqueCtr+1
    tamaño_bloque= np.arange(tamaño_bloque_minimo,tamaño_bloque_maximo)
    if grafica:
        plt.errorbar(tamaño_bloque, incertidume,barra_erros, marker='o',ls='None',markersize=8, capsize=6)
        plt.xlabel('Tamaño do bloque (Número de datos por bloque, que deben ter todos o mesmo tamaño)')
        plt.ylabel('Desviación estándar')
        plt.show()
        plt.errorbar(tamaño_bloque, Media_bloque, incertidume,marker='o',ls='None',markersize=8, capsize=6,color='Orange')
        plt.ylabel('Media')
        plt.xlabel('Tamaño do bloque')
        plt.show()
        plt.tight_layout()
        plt.show()
        print(f"Valor medio = {Media_bloque[-1]:.3f} +/- {incertidume[-1]:.3f}")
    return tamaño_bloque, Media_bloque, incertidume
    
dato = np.loadtxt('output.txt',skiprows=2)
#print(dato)
E=dato[:,0]

blockAverage(E,tamaño_bloque_maximo=False,grafica=True)
