# -*- coding: utf-8 -*-
"""
Created on Fri Jun  12 17:17:12 2019

@author: Ulysse
-
Version légère avec des gravitons classiques représentés à part.
"""
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter#Used for 3d plotting

import numpy as np
import matplotlib.pyplot as plt
from cmath import cos
from cmath import sin
from cmath import pi
from cmath import polar
import scipy as sc
import scipy.linalg
import scipy.signal as signpy
from collections import defaultdict
from time import time
from math import sqrt as msqrt
from random import SystemRandom
def rand1():return SystemRandom().random()
I = np.identity(4, dtype = complex)
"""
#Quand il n'y a aucun graviton
T0 = np.array([
[1,0,0,0],
[0,0,1,0],
[0,1,0,0],
[0,0,0,-1],
])

Kt0 = np.kron(I,T0)


#Un seul graviton : la particule repart dans direction inverse.
T1 = np.array([
[1,0,0,0],
[0,1,0,0],
[0,0,1,0],
[0,0,0,-1],
])#demander à Giuseppe pour savoir si -1 sur dernière ligne?

Kt1 = np.kron(I,T1)

#DEUX gravitons : la particule repart dans direction inverse.
T2 = np.array([
[1,0,0,0],
[0,0,1,0],
[0,1,0,0],
[0,0,0,1],
])#demander à Giuseppe pour savoir si -1 sur dernière ligne?

Kt2 = np.kron(I,T2)"""


theta = 10e-2 * pi
c = np.cos(theta)
#print('%.30f' % c.real)
#print('%.30f' % np.cos(theta))
s = 1j * np.sin(theta)
#print(s)
#print(c)
#print(s*s+c*c)

U = np.array([
[1,0,0,0],
[0,c,s,0],
[0,-s,c,0],
[0,0,0,1],
])

#Ku = np.kron(I,U)

#print(U)
#print(np.dot(U,np.conj(U.T)))


def roundd(nar):
	return round(nar,6)

def abs(x):return np.abs(x)

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def graph(density,x ='Position', y ='Probabilité de présence',T=None):
	if T is None:
		plt.plot(density)
	else:
		plt.plot(T,density)	
	plt.ylabel(y)
	plt.xlabel(x)
	plt.grid()
	plt.show()

def normA(rep):return np.sqrt(np.sum(np.vectorize(lambda x : x*x)(rep)))
def normD(p):
	s = 0	
	for ((x,d),a) in p.items():
		s += abs(a)**2
	return msqrt(s)
def fusionLR(p):
	z = defaultdict(lambda: 0)
	for ((x,_),a) in p.items():
		z[x] += abs(a)**2
	for (x,a) in z.items():
		z[x] = msqrt(a)
	return z
def maximalPos(z):
	maxi = 0
	maxx = 0
	for (x,a) in z.items():
		if a > maxi:
			maxi = a
			maxx = x
	return maxx
def dict_to_array(p,size=1000,LR=True):
	if LR:
		flag = False
		rep = np.zeros((size,2))
		for ((x,d),a) in p.items():
			if a == 0 and flag:
				break
			elif a != 0 and not flag:
				flag = True
			if d == 0:
				rep[x,0] += abs(a) * abs(a)
			else:
				rep[x,1] += abs(a) * abs(a)
	else:
		rep = np.zeros(size)
		for ((x,d),a) in p.items():
			if abs(a) > 0: rep[x] += abs(a) * abs(a)
	rep = np.sqrt(rep)
	return rep
def inputt(st):
	try:
		return int(input(st))
	except:
		return 0

def basic_QW(N=None,it_number=None,th=None,r=None,halved=True,gauss=False,affreg=None,freq=None,disp=True,creat=True,ignore=False,affiche=True,d3=False):
	"""
	N = taille de la grille
	it_number = nombre d'itérations
	th = angle
	r = rayon du "trou noir" 
	halved : on n'affiche que la moitié intéressante de la QW
	gauss  : Donne t'on une gaussienne en entrée? 
	affreg : activer l'affichage régulier 
	freq   : fréquence des renormalisations et de l'affichage
	disp   : activer la disparition des gravitons
	creat  : activer la création des gravitons
	ignore : toujours ignorer les indications d'évènements
	affiche: affiche des informations sur l'évolution à la fin
	d3     : active l'affichage 3d
	"""
	
	if N is None:
		N = 1000#position de départ

	if it_number is None:
		it_number = 450#nombre d'itérations
	if r is None:
		r = N * 2
	if freq == None:
		freq = int(100)

	timing = False
	p = defaultdict(lambda: 0)
	if gauss:
		tau = int(N/10)+1
		ta = signpy.gaussian(tau, tau/10, True)
		print(ta)
		for i in range(tau):
			p[N + i,0] = ta[i]
		N += int(tau)#eviter debordements
	else:
		p[N,0] = 1
	print("Norme au temps"+str(0)+": " + str(normD(p)))
	cpt = -1
	Xm = []
	gr = defaultdict(lambda: False)#gravitons allant à droite (Right)
	gl = defaultdict(lambda: False)#Left-handed gravitons
	#pour qu'il y ait déjà des gravitons au début:
	gl[10] = True
	gr[0] = True
	gr[1]  = True
	gr[2]  = True
	gr[3]  = True
	gr[4]  = True
	gr[5]  = True
	gr[7]  = True
	gr[10] = True
	gr[12] = True
	gr[14] = True
	gr[16] = True
	gr[18] = True
	gr[100] = True
	#gl[3000] = True
	X = np.arange(it_number+N+1)
	threeDplot = []
	if halved:
		X = X[::2]
	for t in range(it_number):
		if timing:
			tcapt = time()
		if creat:
			gr[0] = True
			gr[-1] = True

		""" application des unitaires """
		L2 = defaultdict(lambda: 0)
		
		for ((x,d),a) in p.items():
			grx = gr[x]
			glx = gl[x]
			if abs(a) > 0.1:
				#print((glx,grx,x,d))
				if grx == True and not ignore:
					print("Il se passe un truc en x="+str(x))
					cpt = freq
			
			if d == 0:
				if grx == glx:# en log ! Donc, paradoxalement, mieux qu'avec une liste :
					q = (x - 1, 0)#la somme de tous les couts de recherche + mise à jour est en nlog n ! 
				else:
					q = (x + 1, 1)
			else:
				if grx == glx:
					q = (x + 1, 1)
				else:
					q = (x - 1, 0)
			L2[q] += a

		if timing:
			tim = [(time() - tcapt,time() - tcapt)]
		
		L = defaultdict(lambda: 0)
		for ((x,d),a) in L2.items():#U
			if gr[x-1] == 0 and gl[x-1] == 0 and gr[x+1] ==  0 and gl[x+1] == 0 :#d'où on vient -> est ce qu'un graviton va nous accompagner (et nous bloquer) ?
				L[x,d] += (a*c)
				L[x,1-d] += (a*s)
			else:
				L[x,d] += a## ? condition ci-dessus complexifiée pour que les matrices restent unitaires.
				
		"""mise à jour de p et renormalisations """
		
		if d3: 
			rep = dict_to_array(p,it_number+N+1,False)
			threeDplot.append(rep)
		else:
			Xm.append(maximalPos(fusionLR(p)))
		
		if timing:
			tim.append((time()-tcapt, time() - tim[0][1] - tcapt))

		cpt += 1
		if cpt >= freq:
			print("Norme au temps"+str(t)+": " + str(normD(p)))
			n = normD(L)
			p = defaultdict(lambda: 0)
			for ((x,d),a) in L.items():
				p[x,d] = a/n
		else:
			p = defaultdict(lambda: 0)
			for ((x,d),a) in L.items():
				p[x,d] = a
		
		if timing:
			tim.append((time()-tcapt, time() - tim[1][1] - tcapt))
		
		""" Mise à jour des gravitons """
		
		gl0 = defaultdict(lambda: False)
		gr0 = defaultdict(lambda: False)
		for k,v in gl.items():
			if k > 0:
				gl0[k - 1] = v#Mise à jour des gravitons
			#sinon, il disparaît
		for k,v in gr.items():
			if k < N + it_number + 1:
				gr0[k + 1] = v#Mise à jour des gravitons
		gr = gr0
		gl = gl0
		
		
		if disp:#disparition des gravitons  => !!! A CHANGER, pour changer les particules ! 
			for x,a in gr.items():
				try:
					if rand1() > r * (1 - (x+r)/(x+r+1)) and gr[x]:
						gr[x] = False
						#if p[x,0] > 0:
						a,phi = polar(p[x,1])
						b,psi = polar(p[x,0])
						if rand1() < a*a/(a*a + b*b):
							p[x, 0] = (np.sqrt(a*a + b*b)) * np.exp(1j * phi)#renormalisation -> mesure : il n'y a pas de gravitons ! 
						else:
							p[x, 0] = (np.sqrt(a*a + b*b)) * np.exp(1j * psi)
						p[x, 1] = 0
						#p[x,1] = p[x,0] - p[x,1]
						#p[x,0] = p[x,0] - p[x,1]
				except:
					pass
			for x,a in gl.items():
				if rand1() < 0.001:
					gl[x] = False

		if timing:
			tim.append((time()-tcapt, time() - tim[2][1] - tcapt))
			print(tim)

			
		""" affichages """
		if cpt >= freq:
			cpt -= freq
			if affreg:#affichage régulier
				rep = dict_to_array(p,it_number+N+1)
				if halved:
					rep = rep[::2]
				for x,a in gr.items():
					if a: plt.scatter(x,0.05,marker='4')
				for x,a in gl.items():
					if a: plt.scatter(x,0.05,marker='3',s=100)
				plt.plot(X, rep)
				plt.grid()
				plt.show()
				#graph(rep,T=X, x="Position", y="Densité linéique de Probabilisé de Présence")
			if not ignore:
				gDay = int(inputt("Combien d'étapes passer ?"))
				cpt -= gDay
				if gDay == 666:
					cpt += 666
					creat = False
					disp = False
					gr = defaultdict(lambda:0) 
				if gDay == 777:
					#cpt += 777
					ignore = True
				if gDay == -1:
					cpt -= 1
					timing = True 
			print(cpt)
			if timing:
				try:
					tim.append((time()-tcapt, time() - tim[3][1] - tcapt))
				except NameError:
					tim = []
		if timing:		
			print(tim)
			int(inputt("The time is displayed just over this line."))
	#for x in range(1000,1050):
	#	print(p[x,0])
	#	print(p[x,1])
	if affiche:
		print(affiche)
		print("Affichage", end = '')
		rep = dict_to_array(p,it_number+N+1)
		if halved:
			rep = rep[::2]
		plt.plot(X, rep)
		print("Somme : " + str(normA(rep)))
		print(".", end = '')
		plt.grid()
		for x,a in gr.items():
			if a:plt.scatter(x,0.1,marker='4')
		for x,a in gl.items():
			if a:plt.scatter(x,0.1,marker='3',s=100)
		print(".", end = '')
		print("Position de départ: "+str(N)+", Durée écoulée: "+str(it_number))
		plt.show()
		print(".")
		graph(Xm,x="Temps",y="Position",T=np.array(range(it_number)))
	if d3:
		Xm = np.empty(it_number)
		Xm[0] = N
	return Xm,threeDplot



def meaned_QW(mean,N=None,it_number=None,th=None,r=None,halved=True,gauss=False,affreg=False,freq=None,disp=True,creat=True,ignore=True,affiche=False,d3=True):
	
	
	X = [[] for _ in range(mean)]
	Z0 = [[] for _ in range(mean)]
	T = np.array(range(it_number))
	for m in range(mean):
		X[m],Z0[m] = basic_QW(N=N,it_number=it_number,th=th,r=r,halved=halved,gauss=gauss,affreg=affreg,freq=freq,disp=disp,creat=creat,ignore=ignore,affiche=affiche,d3=d3)
	Xa = np.empty(it_number)
	for t in range(it_number):
		Xa[t] = np.sum([X[i][t] for i in range(mean)])  / mean

	if d3:
		#print(Z0[0])
		Z = np.zeros((len(Z0[0]),len(Z0[0][0])))
		for k in range(len(Z0[0])):
			for t in range(len(Z0[0][0])):
				Z[k,t] = np.sum([Z0[i][k][t] for i in range(mean)])  / mean
		fig = plt.figure()
		#ax = fig.add_subplot(111, projection='3d')
		ax = fig.gca(projection='3d')
		X = np.arange(0,int(Xa[0])+len(Xa)+1,1)
		T = np.arange(0,len(Xa),1)
		print(np.shape(T))
		X, T = np.meshgrid(X, T)
		print(np.shape(X))
		print(np.shape(T))
		print(np.shape(Z))
		surf = ax.plot_surface(X, T, Z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
		ax.set_zlim(-0.01, 1.01)
		ax.zaxis.set_major_locator(LinearLocator(10))
		ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

		fig.colorbar(surf, shrink=0.5, aspect=5)
	else:
		#print(Xa)
		plt.plot(T,Xa)
		plt.grid()
	plt.show()	

	ffnorm = open("norms.txt","w")
	ffnorm.write("") 
	ffnorm.close() 
	ffnorm = open("norms.txt","a")
	ffnorm.write("Xh = "+str(Xa))
	ffnorm.close()

#print("La particule déterminée en position se disperse.")
meaned_QW(2,150,300,gauss=False,halved=False,affreg=False,disp=True,creat=True,d3=True) 

"""
p[0] = case_memoire[0] * state[j+1,0]#proba qu'il y ait |0000>
p[1] = case_memoire[0] * state[j+1,2]#proba qu'il y ait |0010>
p[2] = case_memoire[1] * state[j+1,0]#proba qu'il y ait |0100>
p[3] = case_memoire[1] * state[j+1,2]#proba qu'il y ait |0110>


p[4] = case_memoire[0] * state[j+1,1]#proba qu'il y ait |0001>
p[5] = case_memoire[0] * state[j+1,3]#proba qu'il y ait |0011>
p[6] = case_memoire[1] * state[j+1,1]
p[7] = case_memoire[1] * state[j+1,3]#|0111>


p[8] = case_memoire[2] * state[j+1,0]#proba qu'il y ait |1000>
p[9] = case_memoire[2] * state[j+1,2]#proba qu'il y ait |1010>
p[10] = case_memoire[3] * state[j+1,0]
p[11] = case_memoire[3] * state[j+1,2]			

p[12] = case_memoire[2] * state[j+1,1]#proba qu'il y ait |1001>
p[13] = case_memoire[2] * state[j+1,3]#proba qu'il y ait |1011>
p[14] = case_memoire[3] * state[j+1,1]
p[15] = case_memoire[3] * state[j+1,3]
"""
"""
p[0] = case_memoire[0] * case_memoire[2] *  state[j+1,0] *  state[j+1,2]#proba qu'il y ait |0000>
p[1] = case_memoire[0] * case_memoire[2] *  stat                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      e[j+1,1] *  state[j+1,2]#proba qu'il y ait |0010>
p[2] = case_memoire[0] * case_memoire[3] *  state[j+1,0] *  state[j+1,2]
p[3] = case_memoire[0] * case_memoire[3] *  state[j+1,1] *  state[j+1,2]


p[4] = case_memoire[0] * case_memoire[2] *  state[j+1,0] *  state[j+1,3]#proba qu'il y ait |0001>
p[5] = case_memoire[0] * case_memoire[2] *  state[j+1,1] *  state[j+1,3]#proba qu'il y ait |0011>
p[6] = case_memoire[0] * case_memoire[3] *  state[j+1,0] *  state[j+1,3]
p[7] = case_memoire[0] * case_memoire[3] *  state[j+1,1] *  state[j+1,3]


p[8] = case_memoire[1] * case_memoire[2] *  state[j+1,0] *  state[j+1,2]#proba qu'il y ait |1000>
p[9] = case_memoire[1] * case_memoire[2] *  state[j+1,1] *  state[j+1,2]#proba qu'il y ait |1010>
p[10] = case_memoire[1] * case_memoire[3] *  state[j+1,0] *  state[j+1,2]
p[11] = case_memoire[1] * case_memoire[3] *  state[j+1,1] *  state[j+1,2]			

p[12] = case_memoire[1] * case_memoire[2] *  state[j+1,0] *  state[j+1,3]#proba qu'il y ait |1001>
p[13] = case_memoire[1] * case_memoire[2] *  state[j+1,1] *  state[j+1,3]#proba qu'il y ait |1011>
p[14] = case_memoire[1] * case_memoire[3] *  state[j+1,0] *  state[j+1,3]
p[15] = case_memoire[1] * case_memoire[3] *  state[j+1,1] *  state[j+1,3]
"""
