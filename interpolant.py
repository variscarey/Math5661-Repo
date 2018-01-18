# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 08:41:57 2016

@author: Undertrick
"""
import numpy as np
class Interpolant:

    
    def __init__(self,datax,datay,method='Newton'):
        self.datax=datax
        self.datay=datay
        self.method=method
        self.degree=np.size((datay))-1
        if self.method=='Lagrange':
            #Lagrange interpolant
            self.coeff=np.copy(self.datay)
            for i in range(self.degree+1):
                for j in range(self.degree+1):
                    if j!=i:
                        self.coeff[i]/=(self.datax[i]-self.datax[j])
        elif self.method=='Monomial':
            # build vandermonde matrix & solve for coeffs
            vander=np.ones((self.degree+1,self.degree+1),dtype=float)
            poly=self.datax
            for i in range(1,self.degree+1):
                vander[:,i]=vander[:,i-1]*poly
            self.coeff=np.linalg.solve(vander,self.datay)
        else:
            # Newton interpolant
            self.divdiff=np.zeros((self.degree+1,self.degree+1))
            self.divdiff[0,:]=self.datay
            for i in range(1,self.degree+1):
                for j in range(1,self.degree+2-i):
                    self.divdiff[i,j-1]=(self.divdiff[i-1,j]-self.divdiff[i-1,j-1])/(self.datax[j+i-1]-self.datax[j-1])
            self.coeff=self.divdiff[:,0]
            
    

    
    def eval(self,input):
        output=np.ones(np.shape(input)) #*self.coeff[0]
        if self.method=='Monomial':
            poly=np.ones(np.shape(input))
            output*=self.coeff[0]
            for i in range(1,self.degree+1):
                poly*=input
                output+=poly*self.coeff[i]
        elif self.method=='Lagrange':
            output=np.zeros(np.shape(input))
            for i in range(self.degree+1):
                temp=np.ones(np.shape(input))
                for j in range(self.degree+1):
                    if j!=i:
                        temp*=input-self.datax[j]
                output+=temp*self.coeff[i]
        else: #Newton 
            output*=self.coeff[0]
            poly=np.ones(np.shape(input))
            for i in range(1,self.degree+1):
                poly*=(input-self.datax[i-1])
                output+=poly*self.coeff[i]
        return output
            
