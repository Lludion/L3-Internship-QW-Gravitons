# -*- coding: utf-8 -*-
"""
Created on Fri Jun  28 18:24:43 2019

@author: Ulysse
-
Version lourde avec gravitons et particules quantiques.
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
	for (_,a) in p.items():
		s += abs(a)**2
	return msqrt(s)
def fusionLR(p):
	z = defaultdict(float)
	for (c,a) in p.items():
		x1 = c[0]
		d1 = c[1]
		x2 = c[2]
		d2 = c[3]
		z[x1] += abs(a)**2
		z[x2] += abs(a)**2
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
def dict_to_array(p,size=1000,LR=False):
	if LR:
		flag = False
		rep = np.zeros((size,4))
		for (c,a) in p.items():
			x1 = c[0]
			d1 = c[1]
			x2 = c[2]
			d2 = c[3]
			if a == 0 and flag:
				break##############probablement a changer dans le futur
			elif a != 0 and not flag:
				flag = True
			print("man up")
	else:
		rep = np.zeros((size,2))
		for (c,a) in p.items():
			x1 = c[0]
			x2 = c[2]
			if abs(a) > 0:
				rep[x1,0] += abs(a) * abs(a)
				rep[x2,1] += abs(a) * abs(a)		
	rep = np.sqrt(rep)
	return rep
def inputt(st):
	try:
		return int(input(st))
	except:
		return 0

def ligrav(c):
	"""renvoie la liste des couples position direction des gravitons (non originellement couplés)"""
	n = len(c)
	t = n - 4
	gg = []
	for i in range(t//2):
		gg.append((c[2*i + 4], c[2*i + 5]))
	return gg

def lgrav(c):
	"""renvoie la liste des couples position direction des gravitons (originellement couplés)"""
	gg = []
	#print("!!"+str(c))
	for i in range(4,len(c)):
		gg.append(c[i])
	return gg

def reordonneG(gg):
	"""renvoie une liste de gravitons parfaitement reordonnee"""
	return gg#TODO - very important

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
	3d     : active l'affichage 3d

	Organisation du dictionnaire;
	p[x1,x1,x2,d2,(g,dg),(g2,dg2),(g3,dg3),....] =  . . . 
	"""

	
	theta = 5*10e-3 * pi
	c = np.cos(theta)
	s = 1j * np.sin(theta)
	psi = 10e-3 * pi
	k = np.cos(psi)
	z = 1j * np.sin(psi)

	if N is None:
		N = 50#position de départ
	taillemax = it_number+N+1+int(N/5)
	if it_number is None:
		it_number = 450#nombre d'itPythonérations
	if r is None:
		r = N * 2
	if freq == None:
		freq = int(1)
	timing = False


	p = defaultdict(float)
	p[N,1,N+int(N/5),0,(N+int(N/5)+10,0)] = 1
	print("Norme au temps"+str(0)+": " + str(normD(p)))
	cpt = -1
	Xm = []
	X = np.arange(taillemax)
	threeDplot = []
	if halved:
		X = X[::2]


	for t in range(it_number):
		#print("T T T T T T T T T "+str(t)+"T T T T T T T T T  ")
		if timing:
			tcapt = time()
		if creat:
			pass#ajout d'un truc dependant de la densité de probabilité de présence

		""" application des unitaires """
		L = defaultdict(float)
		
		for (ct,a) in p.items():
			"""on a un certain scenario c"""
			gg = lgrav(ct)
			
			x1 = ct[0]
			d1 = ct[1]
			x2 = ct[2]
			d2 = ct[3]
			#if abs(a) > 0.1:
			#	if not ignore and True:#rajouter une condition si l'envie me prend
			#		print("Il se passe un truc en x="+str(x))
			#		cpt = freq

			#ETAPE 1 CHANGER LA DIRECTION PERPENDICULAIREMENT AUX GRAVITONS 
			cpt1 = 0
			cpt2 = 0
			for (xg,dg) in gg:
				if xg == x1:
					cpt1 += 1
				if xg == x2:
					cpt2 += 1
			d1 = (d1 + cpt1) % 2#si le compteur est impair, on change de direction.
			d2 = (d2 + cpt2) % 2 
			
			""" ça n'a peut etre aucun sens d'un point de vue graviton quantique. Ceci est peut etre ABSURDE !!!"""
			#ETAPE 2 : ABSORPTION DE GRAVITONS
			""" dans cette partie il faut faire tr_s tr_s attention !!!!!!!!!!"""
			"""
			Qgraw = []
			if cpt1:
				gg1 = []
				for (x,d) in gg:
					if x != x1:
						gg1.append((x,d))
				Qgraw.append((gg1,1j* s))
				Qgraw.append((gg,c))
			"""
			#ETAPE 3: DEPLACEMENT DES PARTICULES
			
			x1 += -1 + 2 * d1
			x2 += -1 + 2 * d2
			
			
			""" attention aux doublons d'états ! -> imposer un ordre total sur les gravitons dans la clé du dictionnaire (gauche puis droite)"""
			""" et maintenir efficaement cet ordre"""
			#ETAPE 4: EMISSION DE GRAVITONS
			""" dans cette partie il faut faire tr_s tr_s attention !!!!!!!!!!"""
			Qgrav = []
			"""Qgrav sera une liste de tuples (liste de gravitons, amplitude de ce scénario), où la 
			liste de gravitons est une liste de tuples (position,direction)"""
			if cpt1 == 0:
				Qgrav.append((gg,c))
				Qgrav.append((gg+[(x1,d1)],1j*s))
			else:
				Qgrav.append((gg,1))
			
			gg = []#ici gg est une liste de liste de gravitons
			for (gh,amp) in Qgrav:
				if cpt2 == 0:
					gg.append((gh,amp*c))
					gg.append((gh+[(x2,d2)],1j*amp*s))#"""on ne risque pas de doublon ici par invariant + principe d'exclusion de Pauli"""
					
				else:
					gg.append((gh,amp))
			Qgrav = gg#ici Qgrav est la liste de liste de gravitons
			del gg
			for (gg,amp) in Qgrav:
				aamp = a * amp				
				#ETAPE 5: MOUVEMENT DES GRAVITONS 	
				#print(gg)
				#print(Qgrav)
				gh = [(x-1+2*d,d) for (x,d) in gg]
				gg = gh#pointer change O(1)
				"""                          - !!! -                        """
				#reordonnement des gravitons pour eviter les états doubles
				gg = reordonneG(gg)
				
				#ETAPE 6: DIFFUSION DES PARTICULES MASSIVES
				""" remarquons que si cpti est à 1 ou 2, on a été gelé par un graviton, on ne doit pas faire s'écouler le temps pour la particule i, donc empêcher la diffusion."""
				tg = tuple(gg)
				#print(tg)
				caa = (x1,d1,x2,d2) + tg
				
				C1 = []
				if cpt1 == 0:
					cba = (x1,1-d1,x2,d2) + tg
					
					#print("boolin"+str(cba))
					C1.append((caa,c*aamp))
					C1.append((cba,1j*s*aamp))
				else:
					C1.append((caa,aamp))
				
				for (cax,b) in C1:
					if cpt2 == 0:
						cab = (cax[0],cax[1],x2,d2) + tg
						
						cbb = (cax[0],cax[1],x2,1-d2) + tg
						
						L[cab] += b*k
						L[cbb] += 1j*b*z
					else:
						L[cax] += b
				#print("Straight up boolin"+str(L))
			
					
		if timing:
			tim = [(time() - tcapt,time() - tcapt)]
		
		
		#mise à jour de p et renormalisations 
		
		if d3: 
			rep = dict_to_array(p,taillemax,False)
			threeDplot.append(rep)
		else:
			Xm.append(maximalPos(fusionLR(p)))
		
		if timing:
			tim.append((time()-tcapt, time() - tim[0][1] - tcapt))

		cpt += 1
		if cpt >= freq:
			print("Norme au temps"+str(t)+": " + str(normD(p)))
			n = normD(L)
			p = defaultdict(float)
			for (ct,a) in L.items():
				p[ct] = a/n
		else:
			p = defaultdict(float)
			for (tc,a) in L.items():
				p[tc] = a
		
		if timing:
			tim.append((time()-tcapt, time() - tim[1][1] - tcapt))
		"""
		#Mise à jour des gravitons 
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
		"""
		if timing:
			tim.append((time()-tcapt, time() - tim[2][1] - tcapt))
			print(tim)

			
		""" affichages """
		if cpt >= freq:
			cpt -= freq
			if affreg:#affichage régulier
				rep = dict_to_array(p,taillemax)
				if halved:
					rep = rep[::2]
				#for x,a in gr.items():
				#	if a: plt.scatter(x,0.05,marker='4')
				#for x,a in gl.items():
				#	if a: plt.scatter(x,0.05,marker='3',s=100)
				plt.plot(X, rep)
				plt.grid()
				plt.show()
				#graph(rep,T=X, x="Position", y="Densité linéique de Probabilisé de Présence")
			"""if not ignore:
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
					timing = True """
			#print(cpt)
			if timing:
				try:
					tim.append((time()-tcapt, time() - tim[3][1] - tcapt))
				except NameError:
					tim = []
		if timing:		
			print(tim)
			int(inputt("Whassup my dude"))
	#for x in range(1000,1050):
	#	print(p[x,0])
	#	print(p[x,1])
	if affiche:
		print(affiche)
		print("Affichage", end = '')
		rep = dict_to_array(p,taillemax)
		if halved:
			rep = rep[::2]
		plt.plot(X, rep)
		#print("Somme : " + str(normA(rep)))
		print(".", end = '')
		plt.grid()
		#for x,a in gr.items():
		#	if a:plt.scatter(x,0.1,marker='4')
		#for x,a in gl.items():
		#	if a:plt.scatter(x,0.1,marker='3',s=100)
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
		p#rint(Xa)
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
basic_QW(30,6,300,gauss=False,halved=False,affreg=False,disp=True,creat=True,affiche=True) 

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
