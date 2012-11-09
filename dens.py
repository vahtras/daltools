#!/usr/bin/env python
"""
Utility routines to generate density matrices
"""
#import numpy.oldnumeric as Numeric
import math
from daltools import one,sirifc
from util import full,blocked
def h1diag(nisht,nasht):
   """Generate density from diagonalized one-electron Hamiltonian"""
   h1=one.read('ONEHAMIL','AOONEINT').unpack().unblock()
   return fdiag(h1,nisht,nasht)
def fdiag(F,nisht,nasht):
   """Generate density from diagonalizing input Fock matrix"""
   C=cmo(F)
   return C2D(C,nisht,nasht)
def cmo(F,S=None,file='AOONEINT'):
   """Return MO coefficients from diagonalization of a provided Fock matrix F,
   an optional overlap matrix S. The returned MOs are normalized and ordered 
   with respect to eigenvalue"""
   nbast=F.shape[0]
   if S is not None:
      e,V=(F/S).eigvec()
   else:
      e,V=(F).eigvec()
      S=one.read('OVERLAP',file).unpack().unblock()
   #print "e,V",e,V
   #return V
   C=full.matrix((nbast,nbast))
   #ind=Numeric.argsort(e.T.data)[0]
   #ind=Numeric.argsort(e)
   #ind=e.argsort()
   #print "ind",ind
   N2=(V.T*S*V).diagonal()
   #print "N2",N2
   if 0:
      for i in range(nbast):
	 ii=ind[i]
	 Ni=1.0/math.sqrt(N2[ii])
	 for j in range(nbast):
	    C[j,i]=V[j,ii]*Ni
   else:
      for i in range(nbast):
	 Ni=1.0/math.sqrt(N2[i])
	 C[:,i]=V[:,i]*Ni
   #
   # Finish off with Gram Schmidt because degenerate orbitals are not
   # properly orthonormalized
   #
   C=C.GS(S)
   return C

def C2D(C,nisht,nasht):
   """Given orbitals and occupancy inactive/active return inactive/active densities"""
   Ci=C[:,:nisht]
   Ca=C[:,nisht:nisht+nasht]
   #print "C2D:Ci Ca",Ci,Ca
   Di=2*Ci*Ci.T
   Da=Ca*Ca.T
   #print "C2D",Di,Da
   return Di,Da
def C1D(C,n):
   """Given orbitals C and number of occupied n return density C*C.T"""
   Ci=C[:,:n]
   Di=Ci*Ci.T
   return Di

def ifc(filename='SIRIFC',ifc=None):
   """Return inactive/active densities from orbitals on interface file SIRIFC"""
   if ifc is None: ifc=sirifc.sirifc(filename)
   Di=blocked.matrix(ifc.nbas,ifc.nbas)
   for  isym in range(ifc.nsym):
      if ifc.nish[isym]:
         Ci=ifc.cmo.subblock[isym]
         Cocci=Ci[:,:ifc.nish[isym]]
         Di.subblock[isym]=2*Cocci*Cocci.T
   Di=Di.unblock()
   #print "ifc.nasht",ifc.nasht
   if ifc.nasht:
      cmoa=blocked.matrix(ifc.nbas,ifc.nash)
      for isym in range(ifc.nsym):
         if ifc.nash[isym]:
            nishi=ifc.nish[isym]
            nocci=ifc.nish[isym]+ifc.nash[isym]
            cmoa.subblock[isym]=ifc.cmo.subblock[isym][:,nishi:nocci]
      cmoa=cmoa.unblock()
      Dv=cmoa*ifc.dv.unpack()*cmoa.T
   else:
      Dv=full.matrix((ifc.nbast,ifc.nbast))
   return Di,Dv

def atomic(atomdirlist=[],atomfilelist=[]):
   print atomdirlist
   print atomfilelist
   assert atomdirlist or atomfilelist
   Dilist=[]
   Dalist=[]
   nbas=[]
   nbast=0
   if atomdirlist:
      for atomdir in atomdirlist:
         filename="%s/%s"%(atomdir,'SIRIFC')
         Di,Da=ifc(filename)
         nbas.append(Di.rdim)
         nbast+=Di.rdim
         Dilist.append(Di)
         Dalist.append(Da)
   else:
      for atomfile in atomfilelist:
         filename=atomfile
         Di,Da=ifc(filename)
         nbas.append(Di.rdim)
         nbast+=Di.rdim
         Dilist.append(Di)
         Dalist.append(Da)
   Di=blocked.matrix(nbas,nbas)
   Da=blocked.matrix(nbas,nbas)
   Di.subblock=Dilist
   Da.subblock=Dalist
   S=one.read().unblock().unpack()
   print S
   ibas=[]
   off=0
   SDS=full.matrix(nbast,nbast)
   for nb in nbas:
      ibas.append(off)
      off+=nb
   #print "ibas",ibas
   #print "nbas",nbas
   #print S[0:1,0:1]
   for A in range(len(nbas)):
      #print "A",range(ibas[A],ibas[A]+nbas[A])
      for B in range(len(nbas)):
         #print "B",range(ibas[B],ibas[B]+nbas[B])
         for C in range(len(nbas)):
            #print "C",range(ibas[C],ibas[C]+nbas[C])
            SAC=S[ibas[A]:ibas[A]+nbas[A],ibas[C]:ibas[C]+nbas[C]]
            #print "SAC",SAC
            SCB=S[ibas[C]:ibas[C]+nbas[C],ibas[B]:ibas[B]+nbas[B]]
            #print "SCB",SCB
            DC=Di.subblock[C]+Da.subblock[C]
            #print "DC",DC
            SDS[ibas[A]:ibas[A]+nbas[A],ibas[B]:ibas[B]+nbas[B]] += SAC*DC*SCB
   #D=Di.unblock()+Da.unblock()
   print "SDS,S",SDS,S
   print "SDSeig",(SDS/S).eig()
   C=cmo(-SDS)
   #print "D,S,C",D,S,C
   return C
   

if __name__ == "__main__":
   S=one.read('OVERLAP','AOONEINT').unpack().unblock()
   #
   # 
   #
   import sys
   try:
      nisht=int(sys.argv[1])
      nasht=int(sys.argv[2])
   except IndexError:
      nisht=1
      nasht=1
   print "Input orbital info: nisht=%d nasht=%d"%(nisht,nasht)
   if 1:
      Di,Da=h1diag(nisht,nasht)
      print "h1diag:Inactive/Active trace",S&Di,S&Da
   if 1:
      Di,Da=ifc()
      print "ifc:Inactive/Active trace",S&Di,S&Da
   if 0:
      Di,Da=C2D(atomic(["He","H"]),1,1)
      print "atomic:",Di,Da,S&Di,S&Da

   


