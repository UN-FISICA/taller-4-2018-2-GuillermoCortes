import scipy as sp
import numpy as np

class Derivada:
	def __init__(self,f,metodo="adelante",dx=0.001):
		self.f=f
		self.dx=dx
		if(metodo=="adelante"):
			self.metodo = "adelante"
		elif(metodo=="central"):
			self.metodo = "central"
		elif(metodo=="extrapolada"):
			self.metodo="extrapolada"
		elif(metodo=="segunda"):
			self.metodo="segunda"
	def calc(self,x):
		if(self.metodo=="adelante"):
			return (self.f(x+self.dx) - self.f(x))/self.dx
		elif(self.metodo=="central"):
			return (self.f(x+ (self.dx/2)) - self.f(x -(self.dx/2)))/self.dx
		elif(self.metodo=="extrapolada"):
			a=Derivada(self.f, "central", self.dx)
			b=Derivada(self.f, "central", (self.dx)/2)
			return (4*b.calc(x) - a.calc(x))/3
		elif(self.metodo=="segunda"):
			return (self.f(x + self.dx) + self.f(x - self.dx) -2*self.f(x))/((self.dx)*(self.dx))

class Zeros:
	def __init__(self,f,metodo,error=1e-4, max_iter=100):
		self.f=f
		self.error=error
		self.max_iter = max_iter
		if(metodo=="newton"):
			self.metodo="newton"
		elif(metodo=="bisectriz"):
			self.metodo="bisectriz"
		elif(metodo=="interpolacion"):
			self.metodo="interpolacion"
		elif(metodo=="newton-sp"):
			self.metodo="newton-sp"
		elif(metodo=="fsolve-sp"):
			self.metodo="fsolve-sp"
		elif(metodo=="brentq-sp"):
			self.metodo="brentq-sp"
	def zero(self,vi):
		if(self.metodo=="bisectriz"):
			a=vi[0]
			b=vi[1]
			c=(a+b)/2
			it=0
			while( self.f(c) > self.error or self.f(c) < -(self.error) or it<=self.max_iter):
				if(self.f(c) > 0):
					a=c
				elif(self.f(c) < 0):
					b=c
				c=(a+b)/2
				it=it+1
			return (c,it)
		if(self.metodo=="newton"):
			a=vi
			b=self.f(vi)
			c=Derivada(self.f, "adelante", 0.001)
			d= a - (b/c.calc(vi))
			it=1
			while(self.f(d) > (self.error) or self.f(d) < -(self.error) or it<=self.max_iter):
				a=d
				b=self.f(a)
				c=Derivada(self.f, "adelante", 0.001)
				d= a - (b/c.calc(a))
				it=it+1
			return (c,it)
		if(self.metodo=="interpolacion"):
			a=vi[0]
			b=vi[1]
			c = a - ((b-a)/(self.f(b) - self.f(a)))*self.f(a)
			it=1
			while(self.f(c)> (self.error) or self.f(c) < -(self.error) or it<=self.max_iter):
				if(self.f(c)>0):
					b=c
				elif(self.f(c)<0):
					a=c
				c = a - ((b-a)/(self.f(b) - self.f(a)))*self.f(a)
				it=it+1
			return (c,it)
		if(self.metodo=="newton-sp"):
			return sp.optimize.newton(f,vi)
		if(self.metodo=="brentq-sp"):
		   	return sp.optimize.brentq(f, vi[0],vi[1])
		if(self.metodo=="fsolve-sp"):
			return sp.optimize.newton(f,vi)

if __name__ == "__main__":
	a1=Derivada(np.cos,"adelante")
	da1=a1.calc((np.pi)/2)
	a2=Derivada(np.cos,"central")
	da2=a2.calc((np.pi)/2)
	a3=Derivada(np.cos,"extrapolada")
	da3=a3.calc((np.pi)/2)
	a4=Derivada(np.cos,"segunda")
	da4=a4.calc((np.pi)/2)
	print("Derivada de cos(x) en x = pi/2 usando el metodo de diferencias hacia adelante es de ",da1)
	print("Derivada de cos(x) en x = pi/2 usando el metodo de diferencia central es de ",da2)
	print("Derivada de cos(x) en x = pi/2 usando el metodo de diferencia extrapolada es de ",da3)
	print("Segunda derivada de cos(x) en x = pi/2 es de ", da4)
	c1=Zeros(np.cos,"newton")
	z1=c1.zero(1)
	c2=Zeros(np.cos,"bisectriz")
	z2=c2.zero((1,5))
	c3=Zeros(np.cos,"interpolacion")
	z3=c3.zero((1,5))
	print("El cero de cos(x) usando el metodo de Newton es ",z1)
	print("El cero de cos(x) en el intervalo (1,5), usando el metodo de la bisectriz, es ",z2)
	print("El cero de cos(x) en el intervalo (1,5), usando el metodo de la interpolacion es ",z3)
	c4=Zeros(np.cos,"newton-sp")
	z4=c4.zero(1)
	c5=Zeros(np.cos,"fsolve-sp")
	z5=c5.zero(1)
	c6=Zeros(np.cos,"brentq-sp")
	z6=c6.zero((1,5))
	print("El cero de cos(x), usando el metodo de Newton (Scipy) con punto inicial x=1, es ",z4)
	print("El cero de cos(x) en el intervalo (1,5), usando la funcion brentq de Scipy, es ",z6)
	print("El cero de cos(x), usando la funcion fsolve de Scipy con punto inicial x=1, es ",z5)		
