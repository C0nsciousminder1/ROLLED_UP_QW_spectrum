# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 22:07:36 2019

@author: USUARIO
"""
#Librerias
import numpy as np
from scipy.integrate import simps
from scipy import interpolate 
from filon import cos_integral


def Arc_Lenght(ns,nc,NP,a,b):
    fi_0 = 0.0
    fi_f = 2*nc*np.pi
    beta = b/(2*np.pi)
    Long_Arc=np.zeros(ns)
    DS=np.zeros(NP)  # ds: longitud de Arco infinitesimal 
    X=np.zeros(NP)
    Phi = np.linspace(fi_0, fi_f, ns)
    for i in range(0,ns-1):
        fi_i = Phi[i] 
        fi_n = Phi[i+1]
        fi_grid = np.linspace(fi_i, fi_n, NP) 
        for j in range(1,NP+1):
            x = fi_grid[j-1]      # malla de phi
            ro   = a + beta*x       #\rho for spiral archimides
            X[j-1] = x              # arreglo de phi para entre dos nodos  
            DS[j-1] = np.sqrt(ro**2+beta**2)       #Arc Lenght to Archimedian      
        Long_Arc[i+1] = Long_Arc[i] + simps(DS,X)  #Total Arch Lenght
        #sf= Long_arc[-1]
    return Long_Arc,Phi

def fiss(s,Long_arc,phi):
    x_points = Long_arc.tolist() # el arreglo de numpy phi se pasa a ser una lista
    y_points = phi.tolist()  # el arreglo de numpy Lon_arc se pasa a ser una lista  
    tck = interpolate.splrep(x_points, y_points,s=0)
    return interpolate.splev(s, tck) #,der=0)

#  funcion Spline para angulo polar en funciòn de la longitud de arco (recorrido)
def ssfi(fi,Long_arc,phi):
    x_points = phi.tolist()  # el arreglo de numpy phi se pasa a ser una lista
    y_points = Long_arc.tolist()  # el arreglo de numpy Lon_arc se pasa a ser una lista  
    tck = interpolate.splrep(x_points, y_points,s=0)
    return interpolate.splev(fi, tck) #,der=0)

#  Coordenada radial ro en funciòn de la longitud de arco (recorrido)
def ross(a,b,fi):
    #fi = fiss(s,Long_arc,phi)
    rod = b/(2*np.pi)
    return a+rod*fi
#  Coordenada X en funciòn de la longitud de arco (recorrido)
def xss(a,b,fi):
    #fi= fiss(s,Long_arc,phi)
    ro= ross(a,b,fi)
    return ro*np.cos(fi)
#  Coordenada Y en funciónn de la longitud de arco (recorrido)
def yss(a,b,fi):
    #fi= fiss(s,Long_arc,phi)
    ro= ross(a,b,fi)
    return ro*np.sin(fi)

#Radio de curvatura en funciòn de la longitud de arco (recorrido)
def rads(a,b,fi):
#    fi=fiss(s,Long_arc,phi)
    rod=b/(2*np.pi)
    ro=a+rod*fi
    return  np.sqrt(ro**2+rod**2)
# Donor-electron separation
def reD(a,b,fi_e,fi_d):
    x = xss(a,b,fi_e) 
    y = yss(a,b,fi_e)
    XD = xss(a,b,fi_d)
    YD = yss(a,b,fi_d)
    return np.sqrt((x-XD)**2+(y-YD)**2 + 0.25)


#Parametros: Campo Eléctrico

#---Angulo de campo electrico entre 0 y 2pi?
#adimencional Field Camp_F
    
def CAMP(F_grid):
    nf= len(F_grid)
    Camp_F = np.zeros(nf)
    Camp_F[:]= 0.17*F_grid[:]
    return Camp_F


# Energía de potencial

def Vpot(s,sd,CAMP_F,teta,a,b,Long_arc,phi):
    Vpot = 0
    fi_e=fiss(s,Long_arc,phi)
    fi_d=fiss(sd,Long_arc,phi)
    Rs=rads(a,b,fi_e)
    ro=ross(a,b,fi_e)
    Red= reD(a,b,fi_e,fi_d)
    Vpot=CAMP_F*ro*np.cos(fi_e-teta)-0.25/(Rs**2)- 2.0/Red
    return Vpot

def int_vpot(sf,Nelectron,SD_grid,CAMP_F,teta,a,b,Long_arc,phi): 
    nf=len(CAMP_F)
    #----------- Donor position--------------   
    Ndonor= len(SD_grid)
    phi_D = np.zeros( Ndonor)
    SD_SF = np.zeros(Ndonor)
    #---------# Electron position-------------
    SE_grid = np.linspace(0, sf, Nelectron) 
    phi_e=np.zeros(Nelectron) # Electron position
    SD_SF= np.zeros(Ndonor)
    #-------------------   # Potential#  --------------------------------
    VS=np.zeros(shape=(nf,Ndonor,Nelectron))
    kindex=np.zeros(Nelectron)
    for f  in range(0, nf):
        for j  in range(1, Ndonor+1):
            sd = SD_grid[j-1]
            phi_D[j-1]= fiss(sd,Long_arc,phi)
            SD_SF[j-1] = sd/sf 
            for k in range( 0,Nelectron):
                se =SE_grid[k]
                phi_e[k]= fiss(se,Long_arc,phi)
                kindex[k]= (k*np.pi/sf)
                VS[f,j-1,k] = Vpot(se,sd,CAMP_F[f],teta,a,b,Long_arc,phi)       
    return VS,kindex,SE_grid,SD_SF,Ndonor,phi_e,nf


def coef(VS,kindex,sf,nf,Nelectron,Ndonor):
    ds= sf/(Nelectron)
    COEF=np.zeros((nf,Ndonor,Nelectron))
    for j in range(0,Ndonor):
        for i in range(0,nf):
            COEF[i,j,:]=(1/sf)*(cos_integral(VS[i,j,:],ds, kindex, x0=0.0, axis=0))
    return COEF

def schr(COEF,sf,nF,Nelectron,Ndonor):
    nd=int(Nelectron/2)
    AR=np.zeros(shape=(nd,nd)) # Empty Array "2Dimencional"
    eigenvalues = np.zeros(shape=(nF,Ndonor,nd))
    eigenvectors = np.zeros(shape=(nF,Ndonor,nd,nd))

    for j in range(0,nF):
        for s in range(0,Ndonor):
            for k in range(0,nd):
                for k1 in range(0,nd):
                    ka=abs(k-k1)
                    AR[k,k1]= 0.5*(COEF[j,s,ka] - COEF[j,s,k+k1])
                    if (k == k1):
                        AR[k,k1]=AR[k,k1]+(np.pi*k/sf)**2
            eigenvalues[j,s] ,eigenvectors[j,s] = np.linalg.eigh(AR)
    return eigenvalues,eigenvectors,nd

#  w : (…, M) ndarray  #The eigenvalues in ascending order, each repeated according to its multiplicity.

#  v : {(…, M, M) ndarray, (…, M, M) matrix}  #The column v[:, i] is the normalized eigenvector corresponding to the eigenvalue w[i]. Will return a matrix object if a is a matrix object


def DOS(ancho,Ee,Ndonor,nF,nd,points_plot):
    nEe= nd # #Number of eigenvalues
    Np = points_plot  #Number of points to graph 
    F = np.zeros((Np,Ndonor,nF)) 
    X = np.zeros(Np)
    Emin = -4.0 #np.min(Ee)
    Emax = 0.0 #np.max(Ee)
    Energy_grid = np.linspace(Emin, Emax,Np)
    s= ancho
    def fg(x,x0,s):
        return (s/((x-x0)**2+s**2))/np.pi
    
    for i in range(0,Np):
        x = Energy_grid[i]
        X[i]= x
        for j in range(0,nF):
            for k in range(0,nEe): 
                x0=Ee[j,:,k]
                F[i,:,j] = F[i,:,j] + fg(x,x0,s)
    return X,F

#--------- Wave Function: -----------
#VecR = Cn = Eigenvectors 
# 1 component : different electric fields  : nF: 20
# 2 component : different donor positions  : Ndonor 10 
# 3 component : nd  numbers of row egeinvectores  100
# 4 component : nd numbers of colums egeinvectores 100



