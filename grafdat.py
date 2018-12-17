""" declara la funcion para el ajuste """

#def func(x, A, B):
#   return A**2 * x / (B+x)**2

""" Librerias """ 

#from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
import pylab

"""Importar los datos """

dat=np.loadtxt('datos.txt',delimiter=',',usecols=(0,1),unpack=True)
x = dat[0] 
y = dat[1]


#erry=dat[2]
""" para el ajuste """

#Adat , Bdat  = 2.0 , 2.0
#(A,B), _ = curve_fit(func, x, y)
#print(A , B)  
#Y=func(x,A,B)
#pylab.plot(x,Y,'r',label='Ajuste')

""""Para hacer las barras de error """

#errx=0
#erry=0

#pylab.errorbar(x,y,erry,errx,'.',label='Datos ')
#pylab.title()
#pylab.ylabel(r'P(W)')
#pylab.xlabel(r'$R_L$($\Omega$)')
#pylab.legend()
#plt.semilogx(x,y,".")	
#pylab.grid()
plt.plot(x,y)	
plt.grid()
"""Guardar la imagen """

pylab.savefig('grafdat.pdf')

""" imprimir la figura """
plt.show()
#pylab.show()
