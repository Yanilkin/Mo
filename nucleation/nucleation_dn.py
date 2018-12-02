import matplotlib.pyplot as plt
import numpy as np
import math

import sys
sys.path.append('/home/yanilkin/spparks_code_generator/code')
import defect_properties as defpr

volume = 1.000000e-18 # in Litres
volume = volume/1000
a = 3.14e-10 # lattice constant in m
b = math.sqrt(3)/2*a # Burger's vector sqrt(3)/2*a
Va = (a*a*a)/2.0 # Atomic volume
na = volume/Va # Number of atoms
pi = 3.14 # Pui number

#Parameters
G = 5e-4
Ks = 1.2e13
Kd =5e13
Kiv = 4*3.14*(3.1e-10+2.9e-10)/Va
zi = 3.0
zv = 1.0
Ki = Ks+ zi*Kd
Kv = Ks+ zv*Kd

Ediss = 2.21 #v5c2
N = 4
v1C = 1e-4

def Difi(T):
	return 1e-8*np.exp(-0.4*11600/T)

def Difv(T):
	E_vac_mig = 1.3	#1.1	#eV
	d_vac = a*math.sqrt(3.)/2
	nu = 1.5e13#*55	#3.6e12
	return 8./6*d_vac**2*nu*np.exp(-E_vac_mig*11600/T)

def cv(T):
	Di = Difi(T)
	Dv = Difv(T)
	return np.sqrt(np.power(Ki*Di/(2*Kiv*(Di+Dv)),2)+G*Ki*Di/(Kiv*(Di+Dv)*Kv*Dv))-Ki*Di/(2*Kiv*(Di+Dv))

def ci(T):
	Di = Difi(T)
	Dv = Difv(T)
	return np.sqrt(np.power(Kv*Dv/(2*Kiv*(Di+Dv)),2)+G*Kv*Dv/(Kiv*(Di+Dv)*Ki*Di))-Kv*Dv/(2*Kiv*(Di+Dv))

def rate(T):
	vn = np.exp(Ediss*11600/T)*np.power(cv(T),N)*v1C
	J = 0.5*vn*(4*pi*(defpr.set_imp_vac_rad(5)+defpr.set_imp_vac_rad(1))*Difv(T)/Va)*cv(T)
	return J

def rate_vconst(T):
	v1 = 1e-5
	vn = np.exp(Ediss*11600/T)*np.power(v1,N)*v1C
	J = 0.5*vn*(4*pi*(defpr.set_imp_vac_rad(5)+defpr.set_imp_vac_rad(1))*Difv(T)/Va)*v1
	return J


print defpr.set_imp_vac_rad(5)+defpr.set_imp_vac_rad(1)

T = np.linspace(500,2000,100)
plt.plot(T,rate(T),label='cv=Eq')
plt.plot(T,rate_vconst(T), label='cv=1e-5')
plt.plot(T,cv(T),label='cv')
plt.xlabel('Temperature,K')
plt.ylabel('Nucleation rate, pores/at/s')
plt.yscale('log')
plt.legend(loc='best')
plt.savefig('Ediss_2.21.png')
plt.show()
