#!/usr/bin/env python
import numpy
from util import unformatted,blocked,full

class sirifc(unformatted.FortranBinary):
   """Read data from dalton interface file"""
   ifclabel = "SIR IPH "
   def __init__(self,name="SIRIFC"):
      unformatted.FortranBinary.__init__(self,name)
      self.getdata()
      self.close()
      self._cmo=None
      self._dv=None
      self._pv=None
      self._fock=None
      self._fc=None
      self._fv=None
   def getdata(self):
      self.find(sirifc.ifclabel)

      self.next()
      self.potnuc, self.emy, self.eactive, self.emcscf = self.readbuf(4,'d')
      self.istate, self.ispin, self.nactel, self.lsym  = self.readbuf(4,'i')

      self.next()
      self.nisht, self.nasht, self.nocct, self.norbt, self.nbast, \
      self.nconf, self.nwopt, self.nwoph, self.ncdets, self.ncmot, \
      self.nnashx, self.nnashy, self.nnorbt, self.n2orbt, self.nsym, \
      = self.readbuf(15,'i')

      self.muld2h = numpy.array(self.readbuf(64,'i')).reshape((8,8))

      self.nrhf = numpy.array(self.readbuf(8,'i'))
      self.nfro = numpy.array(self.readbuf(8,'i'))
      self.nish = numpy.array(self.readbuf(8,'i'))
      self.nash = numpy.array(self.readbuf(8,'i'))
      self.norb = numpy.array(self.readbuf(8,'i'))
      self.nbas = numpy.array(self.readbuf(8,'i'))

      self.nelmn1, self.nelmx1, self.nelmn3, self.nelmx3, self.mctype \
      = self.readbuf(5,'i')

      self.nas1 = numpy.array(self.readbuf(8,'i'))
      self.nas2 = numpy.array(self.readbuf(8,'i'))
      self.nas3 = numpy.array(self.readbuf(8,'i'))

      #
      # Skip readin
      #

      return

      self.next()
      ncmot4=max(self.ncmot,4)
      #print "ncmot4 len",ncmot4,len(self.data)
      dbl=self.readbuf(ncmot4,'d')
      #print "dbl:cmo...",dbl
      n=0
      self.cmo=blocked.matrix(self.nbas,self.norb)
      for isym in range(8):
         for mo in range(self.norb[isym]):
            for ao in range(self.nbas[isym]):
               self.cmo.subblock[isym][ao,mo]=dbl[n]
               n+=1
      assert(n == self.ncmot)

      self.next()
      mmashx=max(self.nnashx,4)
      #print "mmashx len",mmashx,len(self.data)
      dbl=self.readbuf(mmashx,'d')
      #print "dbl:dv...",dbl
      self.dv=full.triangular((self.nasht,self.nasht))
      for i in range(self.nnashx):self.dv[i,0]=dbl[i]

      self.next()
      self.next()
      m2orbt=max(self.n2orbt,4)
      #print "m2orbt,len",m2orbt,len(self.data)
      dbl=self.readbuf(m2orbt,'d')
      #print "dbl:fock...",dbl
      self.fock=blocked.matrix(self.norb,self.norb)
      n=0
      for isym in range(8):
         for i in range(self.norb[isym]):
            for j in range(self.norb[isym]):
               self.fock.subblock[isym][j,i] = dbl[n]
               n+=1
      assert (n == self.n2orbt)

      self.next()
      m2ashy=max(self.nnashx**2,4)
      dbl=self.readbuf(m2ashy,'d')
      self.pv=full.matrix((self.nnashx,self.nnashx))
      n=0
      for i in range(self.nnashx):
         for j in range(self.nnashx):
            self.pv[j,i]=dbl[n]
            n+=1
      assert(n == self.nnashx**2)

      self.next()
      mmorbt=max(self.nnorbt,4)
      dbl=self.readbuf(mmorbt,'d')
      self.fc=blocked.triangular(self.norb)
      n=0
      for isym in range(8):
         ij=0
         for i in range(self.norb[isym]):
            for j in range(i+1):
               self.fc.subblock[isym][i,j]=dbl[ij]
               ij+=1
         n+=ij
      assert (n == self.nnorbt)

      self.next()
      mmorbt=max(self.nnorbt,4)
      dbl=self.readbuf(mmorbt,'d')
      self.fv=blocked.triangular(self.norb)
      #print self.fc
      n=0
      #print "self.norb",self.norb
      for isym in range(8):
         ij=0
         for i in range(self.norb[isym]):
            for j in range(i+1):
               self.fv.subblock[isym][i,j]=dbl[ij]
               ij+=1
         n+=ij
      assert (n == self.nnorbt)
      return 
   def get_cmo(self):
      if self._cmo is None:
         #print  "in get_cmo"
         self.file=open(self.name,'rb')
         self.find(self.ifclabel)
         for i in range(3): self.next()
         self.file.close()
         ncmot4=max(self.ncmot,4)
         dbl=self.readbuf(ncmot4,'d')
         n=0
         self._cmo=blocked.matrix(self.nbas,self.norb)
         for isym in range(8):
            for mo in range(self.norb[isym]):
               for ao in range(self.nbas[isym]):
                  self._cmo.subblock[isym][ao,mo]=dbl[n]
                  n+=1
         assert(n == self.ncmot)
      return self._cmo

   def get_dv(self):
      if self._dv is None:
         self.file=open(self.name,'rb')
         self.find(self.ifclabel)
         for i in range(5): 
            self.next()
         mmashx=max(self.nnashx,4)
         dbl=self.readbuf(self.nnashx,'d')
         self._dv=full.triangular.init(dbl)
      return self._dv

   def get_pv(self):
      print "in get_pv"
      if self._pv is None:
         self.file=open(self.name,'rb')
         self.find(self.ifclabel)
         for i in range(7): self.next()
         self.file.close()
         m2ashy=max(self.nnashx**2,4)
         dbl=self.readbuf(m2ashy,'d')
         self._pv=full.matrix((self.nnashx,self.nnashx))
         n=0
         for i in range(self.nnashx):
            for j in range(self.nnashx):
               self._pv[j,i]=dbl[n]
               n+=1
         assert(n == self.nnashx**2)
      return self._pv

   def get_fock(self):
      if self._fock is None:
         self.file=open(self.name,'rb')
         self.find(self.ifclabel)
         for i in range(6): self.next()
         self.file.close()
         m2orbt=max(self.n2orbt,4)
         dbl=self.readbuf(m2orbt,'d')
         self._fock=blocked.matrix(self.norb,self.norb)
         n=0
         for isym in range(8):
            for i in range(self.norb[isym]):
               for j in range(self.norb[isym]):
                  self._fock.subblock[isym][j,i] = dbl[n]
                  n+=1
         assert (n == self.n2orbt)
      return self._fock
   def get_fc(self):
      if self._fc is None:
         self.file=open(self.name,'rb')
         self.find(self.ifclabel)
         for i in range(8): self.next()
         mmorbt=max(self.nnorbt,4)
         dbl=self.readbuf(mmorbt,'d')
         self._fc=blocked.triangular(self.norb)
         #print self.fc
         n=0
         #print "self.norb",self.norb
         for isym in range(8):
            ij=0
            for i in range(self.norb[isym]):
               for j in range(i+1):
                  self._fc.subblock[isym][i,j]=dbl[ij]
                  ij+=1
            n+=ij
         assert (n == self.nnorbt)
      return self._fc
   def get_fv(self):
      if self._fv is None:
         self.file=open(self.name,'rb')
         self.find(self.ifclabel)
         for i in range(9): self.next()
         mmorbt=max(self.nnorbt,4)
         dbl=self.readbuf(mmorbt,'d')
         self._fv=blocked.triangular(self.norb)
         #print self.fc
         n=0
         #print "self.norb",self.norb
         for isym in range(8):
            ij=0
            for i in range(self.norb[isym]):
               for j in range(i+1):
                  self._fv.subblock[isym][i,j]=dbl[ij]
                  ij+=1
            n+=ij
         assert (n == self.nnorbt)
      return self._fv
   dv = property(fget=get_dv)
   pv = property(fget=get_pv)
   cmo = property(fget=get_cmo)
   fock = property(fget=get_fock)
   fc = property(fget=get_fc)
   fv = property(fget=get_fv)
   def __repr__(self):
      retstr=""
      retstr+="Nuclear Potential Energy: %12.6f\n"%self.potnuc
      retstr+="Electronic energy       : %12.6f\n"%self.emy
      retstr+="Active energy           : %12.6f\n"%self.eactive
      retstr+="MCSCF energy            : %12.6f\n"%self.emcscf
      retstr+="State                   : %d\n"%self.istate
      retstr+="Spin                    : %d\n"%self.ispin
      retstr+="Active electrons        : %d\n"%self.nactel
      retstr+="Symmetry                : %d\n"%self.lsym
      retstr+="NISHT                   : %d\n"%self.nisht
      retstr+="NASHT                   : %d\n"%self.nasht
      retstr+="NOCCT                   : %d\n"%self.nocct
      retstr+="NORBT                   : %d\n"%self.norbt
      retstr+="NBAST                   : %d\n"%self.nbast
      retstr+="NCONF                   : %d\n"%self.nconf
      retstr+="NWOPT                   : %d\n"%self.nwopt
      retstr+="NWOPH                   : %d\n"%self.nwoph
      retstr+="NCDETS                  : %d\n"%self.ncdets
      retstr+="NCMOT                   : %d\n"%self.ncmot
      retstr+="NNASHX                  : %d\n"%self.nnashx
      retstr+="NNASHY                  : %d\n"%self.nnashy
      retstr+="NNORBT                  : %d\n"%self.nnorbt
      retstr+="N2ORBT                  : %d\n"%self.n2orbt
      retstr+="NSYM                    : %d\n"%self.nsym
      retstr+="MULD2H:\n"
      for i in range(8):
         retstr+="   "
         for j in range(8):
            retstr+=" %d"%self.muld2h[i,j]
         retstr+="\n"

      def strvec(lab,vec):
         locstr=lab + ":"
         for i in range(8): locstr+=" %d"%vec[i]
         return locstr + "\n"
      retstr+=strvec("NRHF",self.nrhf)
      retstr+=strvec("NFRO",self.nfro)
      retstr+=strvec("NISH",self.nish)
      retstr+=strvec("NASH",self.nash)
      retstr+=strvec("NORB",self.norb)
      retstr+=strvec("NBAS",self.nbas)
      retstr+="NELMN1                  : %d\n"%self.nelmn1
      retstr+="NELMX1                  : %d\n"%self.nelmx1
      retstr+="NELMN3                  : %d\n"%self.nelmn3
      retstr+="NELMX3                  : %d\n"%self.nelmx3
      retstr+="MCTYPE                  : %d\n"%self.mctype
      retstr+=strvec("NAS1",self.nas1)
      retstr+=strvec("NAS2",self.nas2)
      retstr+=strvec("NAS3",self.nas3)
      retstr+="CMO" + str(self.cmo)
      retstr+="DV\n" + str(self.dv) + "\n"
      retstr+="FOCK" + str(self.fock)
      retstr+="PV\n" + str(self.pv) + "\n"
      retstr+="FC" + str(self.fc)
      retstr+="FV" + str(self.fv)
      return retstr

if __name__ == "__main__":
    import sys
    try:
        filename = sys.argv[1]
        ifc = sirifc(name=filename)
        print ifc
    except(IndexError):
        print "Usage: %s [path]/SIRIFC"%sys.argv[0]
        sys.exit(1)
        
