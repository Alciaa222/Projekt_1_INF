# -*- coding: utf-8 -*-
"""
Created on Sat Apr 29 20:52:34 2023

@author: andze&nikola&ala
"""

import numpy as np
from math import *
from argparse import ArgumentParser

class Transformacje:
    
     
    def __init__(self,model):
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "Elipsoida Krasowskiego":
            self.a = 6370245.0
            self.b = 6356863.01877
        else:
            raise NotImplementedError(f"{model} Ta operacja jest niemożliwa!")    
        self.sp = (self.a - self.b) / self.a
        self.e2 = (2 * self.sp - self.sp ** 2) 
        
           
   def hirvonen(self,X,Y,Z): #XYZ zamieniamy na flh
        p = np.sqrt(X**2+Y**2)
        fi = np.arctan(Z/(p*(1-self.e2)))
        while True: #pętla
            N=self.a/np.sqrt(1-self.e2*np.sin(fi)**2)
            h=p/np.cos(fi)-N
            fip=fi
            fi=np.arctan(Z/(p*(1-self.e2*N/(N+h))))
            if abs(fip-fi)<(0.000001/206265):
                break
        l=np.arctan2(Y,X) #lambda
        return(fi,l,h)     
        
    def flh2XYZ(self,fi,l,h):
        while True:
            N=self.a/np.sqrt(1-self.e2*np.sin(fi)**2)
            X=(N+h)*np.cos(fi)*np.cos(l)
            Xp=X
            Y=(N+h)*np.cos(fi)*np.sin(l)
            Z=(N*(1-self.e2)+h)*np.sin(fi)
            if abs(Xp-X)<(0.000001/206265):
                break
        return(X,Y,Z)
    
    
    def pl2000(self,fi,l,m=0.999923):
        l0=0 
        strefa = 0
        if l >np.deg2rad(13.5) and l < np.deg2rad(16.5):
            strefa = 5
            l0 = np.deg2rad(15)
        elif l >np.deg2rad(16.5) and l < np.deg2rad(19.5):
            strefa = 6
            l = np.deg2rad(18)
        elif l >np.deg2rad(19.5) and l < np.deg2rad(22.5):
            strefa =7
            l0 = np.deg2rad(21)
        elif l >np.deg2rad(22.5) and l < np.deg2rad(25.5):
            strefa = 8
            l0 = np.deg2rad(24)
        else:
            print("Punkt poza strefami odwzorowawczymi układu PL-2000")        
        
        # 1 parametry elipsoidy     
        b2 = self.a**2*(1-self.e2)
        ep2 = (self.a**2-b2)/b2
        # 2. Wielkosci pomocnicze     
        dell = l - l0
        t = np.tan(fi)
        ni2 = ep2*(np.cos(fi)**2)
        N=self.a/np.sqrt(1-self.e2*np.sin(fi)**2)
        
        # 3. Długosc luku poludnika 
        A0 = 1- (self.e2/4)-(3*self.e2**2/64)-(5*self.e2**3/256)
        A2 = (3/8)*(self.e2+(self.e2**2/4)+(15*self.e2**3/128))
        A4 = (15/256)*(self.e2**2+((3*self.e2**3)/4))
        A6 = (35*self.e2**3)/3072
        
        sigma = self.a *(A0*fi-A2*np.sin(2*fi)+A4*np.sin(4*fi)-A6*np.sin(6*fi))
        
        # wsolrzedne prostokatne lokalne na plaszczyznie gaussa-krugera
        
        xgk =  sigma    +    ( ((dell**2/2)*N*np.sin(fi)*np.cos(fi))    *    (1   +   ((dell**2/12)*(np.cos(fi)**2)*(5 - t**2 + 9*ni2 + 4*ni2**2))      +         ((dell**4/360)*(np.cos(fi)**4)*(61 - 58*t**2 + t**4 + 270*ni2 - 330*ni2*t**2))))
        
        ygk =  (dell* N * np.cos(fi))  *   ( 1 +  ((dell**2/6)   *   (np.cos(fi)**2)   *  (1 - t**2 + ni2))     +     (((dell**4/120)*(np.cos(fi)**4)) * (5 - (18*t**2) + t**4 + (14 * ni2) - (58*ni2*t**2))))
        
        x2000 = xgk * m 
        y2000 = ygk*m + (strefa *1000000) +500000
        return  x2000, y2000,xgk,ygk



    def XYZ2neu(self,s,alfa,z,fi,l):
            dneu=np.array([s*np.sin(z)*np.cos(alfa),
                           s*np.sin(z)*np.sin(alfa),
                           s*np.cos(z)])
            R=np.array([[-np.sin(fi)*np.cos(l),-np.sin(l),np.cos(fi)*np.cos(l)],
                        [-np.sin(fi)*np.sin(l),np.cos(l),np.cos(fi)*np.sin(l)],
                        [ np.cos(fi),            0.     ,np.sin(fi)]])
            return(R.T @ dneu)

        
 for xyz in dane:    
            if trans[args.trans]=="hirvonen":
                line = obiekt.hirvonen(xyz[0],xyz[1],xyz[2])
                result.append(line)
                ################## nikola dodanie wywolania modelu flh2XYZ
            if trans[args.trans]=="flh2XYZ":
                line = obiekt.flh2XYZ(xyz[0],xyz[1],xyz[2])
                result.append(line)
                ##################### ja dodanie wywolania modelu pl1992
            if trans[args.trans]=="pl1992":
                line = obiekt.pl1992(xyz[0],xyz[1])
                result.append(line)
                ############### Ala dodanie wywolania modelu pl2000
            if trans[args.trans]=="pl2000":
                line = obiekt.pl2000(xyz[0],xyz[1])
                result.append(line)
                ############### Nikola  dodanie wywolania modelu XYZ2neu
            if trans[args.trans]=="XYZ2neu":
                line = obiekt.XYZ2neu(xyz[0],xyz[1],xyz[2],xyz[3],xyz[4])
                result.append(line)






















    
