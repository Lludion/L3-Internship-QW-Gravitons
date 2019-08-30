# -*- coding: utf-8 -*-
"""
Created on Fri Jun  12 17:17:12 2019

@author: Ulysse
"""



import numpy as np
import subprocess
import itertools
import matplotlib.pyplot as plt
import scipy.optimize as op
#import scipy.linalg

import time

import sys

def xfil(filename, globalz=None, localz=None):
	if globalz is None:
		globalz = sys._getframe(1).f_globals
	if localz is None:
		localz = sys._getframe(1).f_locals
	with open(filename, "r") as fh:
		exec(fh.read()+"\n", globalz, localz)

def polynom(p,x):
	
	s = 0
	n = len(p)
	for e in range(n):		
		s += p[n-1-e] * (x ** e)
	return s

def tofit(x, a=1,b=0,c=0):
	return a/(x*x)+b/x+c

def f2(r, c=500, R=1):
	return c - r/R - np.log(abs(r/R)-1)

def graph(a):
	xm, tm = a[0] + 1, len(a)
	T = np.arange(0,tm,1,int)
	#b = np.vectorize(lambda x : 1./x)(a)
	#p = np.polyfit(T,b,20)
	#a = a.astype(dtype=np.float32)
	#T = T.astype(dtype=np.float32)
	#popt,pcov = op.curve_fit(tofit,T,a,p0=[500,1])
	#print("!---------------------------------------------!")
	#print(popt)
	#print("!---------------------------------------------!")
	#c = np.vectorize(lambda x : 1/polynom(p,(x)))(T)
	plt.plot(T,a)
	#plt.plot(T,f2(T,popt[0],popt[1]))
	#plt.plot(T,c,color='orange')
	#plt.plot(T,b,color='red')
	#plt.ylabel('Position')
	#plt.xlabel('Temps')
	#plt.grid()
	#plt.show()
	#print(a)
"""
def geodesics(xm,tm):
	X = np.arange(0,xm,1,int)
	T = np.vectorize(schwarz(xm,xm/100))(X)
	z = T[-1]
	T = np.vectorize(lambda s : s - z)(T)
	plt.plot(T,X)"""


def gourgou(c,R=1):
	""" renvoie une fonction qui
	donne t en fonction de r dans une métrique de E. Gourgoulhon
	où la particule commence à la position c 
	et atteint un horizon des évènements en R
	"""
	return (lambda r :  c - r/R - np.log(abs(r/R)))

def schw(c,D):
	""" renvoie une fonction qui
	donne t en fonction de r dans une métrique radiale de schwarzschild
	où la particule commence à la position c 
	et atteint un horizon des évènements en R=1
	D = position de départ
	"""
	def f(x):
		try:
			return (D-x)/(c*(1-1/x))
		except ZeroDivisionError:
			return np.inf
	return f



def eV(xm,tm,c=1,e=0.1,D=None):
	if D == None:
		D = xm
	T = np.arange(0,tm,1,int)
	X = [D]
	Y = [D]
	r = D/10
	for t in range(1,tm,1):
		x0 = X[-1]
		dt = 1
		X.append( X[-1] - c * dt * (1 - r*r/((x0+r-1)*(x0+r-1))) )
	for t in range(1,tm,1):
		x0 = Y[-1]
		dt = 1
		Y.append( x0 - c * dt * (1 - r/(x0+r-1)))
	plt.plot(T,X,color="red")
	plt.plot(T,Y,color="black")
	S = np.arange(0,xm/100,0.01,int)
	#plt.plot(np.vectorize(gourgou(D,1))(S),S,color="grey")

def expectedValue(xm,tm,c=1,e=0.1,D=None):
	if D == None:
		D = xm
	print("Scaterred")#âååÂøÊ€ç±æþ 0±1=±1
	scatteredX = []
	scatteredT = []
	for x, t in np.ndindex((xm,tm)):
		if abs(x * x * x + (c * t - D) * x * x  -  c * t) < e:
			scatteredX.append(x)
			scatteredT.append(t)
	plt.plot(scatteredT,scatteredX)
	print("Scaterred")

def main(N=1):
	#a = [(0,1),(1,1)]
	#with open("./test.txt","w") as _: pass
	
	for _ in itertools.repeat(None, N):
		subprocess.call("./callOcaml.sh", shell=True)
		xfil("export.txt")
		xfil("export-hybrid.txt")
	plt.ylabel('Position')
	plt.xlabel('Temps')
	plt.grid()
	plt.show()


main(N=1)
print("Simulation terminée.")
