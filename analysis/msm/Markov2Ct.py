#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 13:57:54 2023

@author: albertsmith
"""

import numpy as np
from pyDR.MDtools import vft
from pyDR.MDtools.Ctcalc import Ctcalc
import matplotlib.pyplot as plt

class Markov2Ct():
    def __init__(self,tau:float,state=None,v=None,P=None,v_state=None):
        """
        Class for determining the NMR correlation functions for a Markov model.
        A few different possibilities exist for the input. One may give the
        Markov transition matrix directly (P), or one may provide a list of
        states (state) and the transition matrix will be solved for internally 
        
        *if Marc/Candide/Sereina, you have better algorithms for going from the
        list of states to the transition matrix, then perhaps better to put P
        in directly
        
        Next, we need the direction of the NMR interaction corresponding to each
        state. One may also provide these directly for each state (v_state), or
        one may provide a list of vectors (v, Nx3), and we will solve for v_state
        internally. Note that I will assume motion within each state of the
        Markov model yields an approximately symmetric tensor. For asymmetric
        tensors, I would start to say we should just go ahead and implement this
        for ROMANCE. Note that v/v_state should be provided in some reference
        frame (e.g. alignment of the CA-CB bond to z, CA-C in the xz-plane)
        
        
        Notes:
        Either state or P is required. P supercedes state if both provided
        
        If v is provided instead of v_state, then state is required, regardless
        of whether P is provided. v_state supercedes v.
        
        
        
        Changing the parameters after initialization will probably break this
        
        

        Parameters
        ----------
        tau : float
            Time step between states
        state : list-like, optional
            List of state indices (0 to N-1 for N states). The default is None.
        v : list-like, optional
            Nx3 list/array of vectors corresponding. The default is None.
        P : array, optional
            Transition matrix. The default is None.
        v_state : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        None.

        """
        
        assert state is not None or P is not None,"P or state is required"
        if v_state is None:
            assert v is not None,"v_state or v is required"
            assert state is not None,"If v_state is not provided, state is required"
         
        self.tau=tau
        
        #Eliminate empty states
        state=np.array(state)
        k=0
        while k < state.max()+1:
            if (state==k).sum()==0:
                state[state>k]-=1
                if v_state is not None:
                    v_state=np.concatenate((v_state[:k],v_state[k+1:]))
            else:
                k+=1
        
        self.state=None if state is None else np.array(state,dtype=int)
        self.v=None if v is None else (np.array(v).T/np.sqrt((v**2).sum(1))).T
        self._vstate=None if v_state is None else (np.array(v_state).T/np.sqrt((v_state**2).sum(1))).T
        self._P=None if P is None else np.array(P)
        
        self._Kex=None
        self._eig=None
        self._S2=None
        self._S2dir=None
        self._S2vstate=None
        self._P2=None
        self._i=None
        self._A=None
        self._A_unsorted=None
        self._delta=None
        self._eta=None
        self._pop=None
        self._Ctdir=None
        self._Ctvstate=None
        self._Ctmarkov=None
        
    @property
    def Nstates(self):
        """
        Returns the number of states in the Markov model

        Returns
        -------
        int

        """
        
        if self._P is None:
            return self.state.max()+1
        return self.P.shape[0]
    
    @property
    def P(self):
        """
        Returns the transition matrix. Performs calculation of the matrix if not
        provided.
        
        I assume P is set up for right multiplication:
            s1=s0@P

        Returns
        -------
        np.array

        """
        if self._P is None:
            P=np.zeros([self.Nstates,self.Nstates])
            for k in range(self.Nstates):
                i=(self.state==k)[:-1]
                nk=i.sum()
                for m in range(self.Nstates):
                    P[k,m]=(((self.state[1:][i]==m).sum()))/nk

            self._P=P
                
                        

        return self._P
    
    def _diagonalize(self):
        """
        Performs diagonalization of the transition matrix, as well as sorting
        of the eigenvalues and eigenvectors. Sorting is largest to smallest,
        but since the eigenvalues are negative, this is smallest magnitude to
        largest magnitude.

        Returns
        -------
        tuple
            eigenvalues and eigenvectors

        """
        if self._eig is None:
            out=np.linalg.eig(self.P.T)
            i=np.argsort(out[0])[::-1]
            self._eig=out[0][i],out[1][:,i]
        return self._eig
        
    
    @property
    def peq(self):
        self._diagonalize()
        return self.U[:,0].real/self.U[:,0].real.sum()
    
    @property
    def U(self):
        """
        Eigenvectors of the transition/exchange matrices

        Returns
        -------
        np.array
            Square matrix with eigenvectors on the columns.

        """
        self._diagonalize()
        return self._eig[1]
    
    @property
    def LambdaP(self):
        """
        Eigenvalues of the transition matrix

        Returns
        -------
        np.array
            1D array of the eigenvalues of the transition matrix.

        """
        self._diagonalize()
        return self._eig[0]
    
    @property
    def Lambda(self):
        """
        Eigenvalues of the exchange matrix

        Returns
        -------
        np.array
            1D array of the eigenvalues of the exchange matrix

        """
        return np.log(self.LambdaP)/self.tau
    
    @property
    def tc(self):
        """
        Correlation times occuring in the correlation function.
        Also referred to as implied timescales.
        
        Note that these are sorted from largest amplitude to smallest amplitude,
        and are not in the same order as the Lambda or LambdaP
        

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        if self._i is None:
            _=self.A
        return -1/self.Lambda[1:][self._i].real

    @property
    def tc_unsorted(self):
        """
        Correlation times occuring in the correlation function.
        Also referred to as implied timescales.
        
        These have the same order as the Lambda or LambdaP
        

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return -1/self.Lambda[1:].real


    @property
    def vstate(self):
        """
        Vectors pointing in the direction of the average dipole coupling for
        each state. Nx3

        Returns
        -------
        np.array

        """
        if self._vstate is None:
            euler=vft.getFrame(self.v.T)
            A0=vft.D2(*euler)*np.sqrt(3/2)
            Aavg=list()
            for k in range(self.Nstates):
                if self.pop[k]==0:
                    Aavg.append([0,0,np.sqrt(3/2),0,0])
                else:
                    Aavg.append(A0[:,self.state==k].mean(1))
            Aavg=np.array(Aavg)
            pars=vft.Spher2pars(Aavg.T)
            self._delta=pars[0]
            self._eta=pars[1]
            self._vstate=vft.R([0,0,1],*pars[2:]).T
            
        return self._vstate
    
    @property
    def delta(self):
        """
        Anisotropy of the averaged tensors due to motion within each Markov state

        Returns None if the vectors (v) have not been provided
        Returns
        -------
        np.array
            1D list of delta values

        """
        if self._delta is None:
            _=self.vstate
        return self._delta
    
    @property
    def eta(self):
        """
        Asymmetry of the averaged tensors due to motion within each Markov state

        Returns
        -------
        np.array
            1D list of eta values

        """
        if self._eta is None:
            _=self.vstate
        return self._eta
    
    @property
    def pop(self):
        """
        Returns a list of populations for each state

        Returns
        -------
        np.array

        """
        if self._pop is None:
            pop=[]
            for k in range(self.Nstates):
                pop.append((self.state==k).sum()/len(self.state))
            self._pop=np.array(pop)
        return self._pop
        
    
    @property
    def S2(self):
        """
        Order parameter for the exchange matrix

        Returns
        -------
        float

        """
        if self._S2 is None:
            pp=np.dot(np.atleast_2d(self.peq).T,np.atleast_2d(self.peq))
            self._S2=(pp*self.P2).sum()
        return self._S2
    
    @property
    def S2dir(self):
        """
        Order parameter calculated directly from the vectors

        Returns
        -------
        float

        """
        if self.v is None:
            return
        if self._S2dir is None:
            S2=0
            for a in self.v.T:
                for b in self.v.T:
                    S2+=((a*b).mean())**2
            S2*=3/2
            S2-=1/2
            self._S2dir=S2
        return self._S2dir
    
    @property
    def S2vstate(self):
        """
        Order parameter calculated from the list of states and the vectors for
        each state (vstate). Removes influence of motion within the states

        Returns
        -------
        float

        """
        if self.state is None:
            return
        if self._S2vstate is None:
            S2=0
            for a in self.vstate.T:
                for b in self.vstate.T:
                    S2+=((a*b*self.pop).sum())**2
            S2*=3/2
            S2-=1/2
            self._S2vstate=S2
        return self._S2vstate
            
    @property
    def P2(self):
        """
        Array of P2 values for all pairs of states in the Markov model

        Returns
        -------
        None.

        """
        if self._P2 is None:
            self._P2=-1/2+3/2*(np.sum([np.atleast_2d(q).T@np.atleast_2d(q) for q in self.vstate.T],0)**2)
        return self._P2
    
    @property
    def A(self):   
        """
        Returns the amplitudes of the individual correlation times.
        
        Ct=S^2+(1-S^2)sum(A_m*exp(-t/tau_m))
        
        returns peq,S2,tau,A
        """
        
        if self._A is None: 
            v=self.U
            
            vi=np.linalg.pinv(v)
            
            A=list()
            for vim,vm in zip(vi,v.T):
                A.append((np.dot(np.atleast_2d(vm).T,np.atleast_2d(vim*self.peq))*self.P2).sum())

            A=np.array(A).real
            
            A=A[1:]/(1-self.S2)
            
            self._i=np.argsort(A)[::-1]  #Sort the A from largest to smallest
            
            self._A_unsorted = A
            self._A=A[self._i]
        
        return self._A
    
    @property
    def t(self):
        if self.v is not None:
            N=len(self.v)
        elif self.state is not None:
            N=len(self.state)
        else:
            N=max(self.tc)*3/self.tau
            
        return np.arange(N).astype(float)*self.tau

    @property
    def Ctdir(self):
        if self._Ctdir is None:
            if self.v is None:
                return
            ctc=Ctcalc()
            for k in range(3):
                for j in range(k,3):
                    ctc.a=self.v[:,k]*self.v[:,j]
                    ctc.c=3/2 if j==k else 3
                    ctc.add()
            self._Ctdir=ctc.Return(-1/2)
        return self._Ctdir[0]
    
    @property
    def Ctvstate(self):
        if self._Ctvstate is None:
            if self.state is None:
                return
            
            v=np.zeros([self.state.shape[0],3])
            for k in range(3):
                v[:,k]=self.vstate[self.state,k]
                
            ctc=Ctcalc()
            for k in range(3):
                for j in range(k,3):
                    ctc.a=v[:,k]*v[:,j]
                    ctc.c=3/2 if j==k else 3
                    ctc.add()
            self._Ctvstate=ctc.Return(-1/2)
        return self._Ctvstate[0]
    
    @property
    def Ct(self):
        if self._Ctmarkov is None:
            self._Ctmarkov=np.ones(self.t.shape)*self.S2
            for A,tc in zip(self.A,self.tc):
                if np.isnan(A) or np.isnan(tc):
                    continue
                self._Ctmarkov+=(1-self.S2)*A*np.exp(-self.t/tc)
        return self._Ctmarkov
    
    def plotCt(self,ax=None,**kwargs):
        """
        Plots the various correlation functions

        Parameters
        ----------
        ax : TYPE, optional
            Provide an axis for plotting. The default is None.
        **kwargs : TYPE
            Parameters passed to matplotlib.

        Returns
        -------
        None.

        """
        if ax is None:ax=plt.subplots()[1]
        
        leg=[]
        t=self.t
        #t[0]=np.exp(np.log(t[1])-np.diff(np.log(t[1:3]))[0])
        
        hdl=list()
        if self.Ctdir is not None:
            ax.plot(t[1:],self.Ctdir[1:],**kwargs,label='Direct', color = 'black')
            #ax.plot(t[[0,-1]],np.ones(2)*self.S2dir,color='black',linestyle='-')
            leg.append('Direct')
        
        if self.Ctvstate is not None:
            ax.plot(t[1:],self.Ctvstate[1:],**kwargs,label='Constructed', color = 'royalblue', linestyle ='dotted')
            #ax.plot(t[[0,-1]],np.ones(2)*self.S2vstate,color=hdl[-1].get_color(),linestyle='--')
            leg.append('Constructed')
            
        ax.plot(t[1:],self.Ct[1:],**kwargs,label='Markov', color = 'darkorange')
        #ax.plot(t[[0,-1]],np.ones(2)*self.S2,color=hdl[-1].get_color(),linestyle='-.')
        leg.append('Markov')
        
        ax.set_xlabel('t / ns')
        ax.set_ylabel('C(t)')
        ax.legend(edgecolor='white', fontsize = 18)
        
        return ax
            
        
        
        