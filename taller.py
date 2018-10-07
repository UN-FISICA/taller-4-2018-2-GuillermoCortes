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

			
