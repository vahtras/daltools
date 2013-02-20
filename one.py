#!/usr/bin/env python
import struct,sys,numpy
from util import unformatted
def readhead(file="AOONEINT"):
   aooneint=unformatted.FortranBinary(file)
   aooneint.readrec()
   if len(aooneint.data) == 144:
      #
      # Newer versions: title in own record
      #
      title=aooneint.data
      #
      # Next record contains MAXREP...
      #
      aooneint.readrec()
   else:
      #
      # Older versions: title ans MAXREP in same record
      #
      title="".join(aooneint.readbuf(192,'c'))
   #
   # Buffer should now contains MAXREP data
   #
   nsym=aooneint.readbuf(1,'i')[0]
   #print "nsym",nsym
   naos=aooneint.readbuf(nsym,'i')
   #print "naos",naos
   potnuc=aooneint.readbuf(1,'d')[0]
   #print "naos",naos
   unlabeled={
         "ttitle":title,
         "nsym":nsym,
         "naos":naos,
         "potnuc":potnuc
         }
   aooneint.close()
   return unlabeled

def readisordk(file="AOONEINT"):
   aooneint=unformatted.FortranBinary(file)
   table1=aooneint.find("ISORDK")
   aooneint.readrec() #dummy
   aooneint.readrec()
   sizeofi=struct.calcsize('i')
   sizeofd=struct.calcsize('d')
   mxcent=(len(aooneint.data)-sizeofi)/(4*sizeofd)
   chrn=aooneint.readbuf(mxcent,'d')
   nucdep=aooneint.readbuf(1,'i')[0]
   cooo=aooneint.readbuf(3*mxcent,'d')
   isordk={
      "table":table1,
      "chrn":chrn,
      "nucdep":nucdep,
      "cooo":cooo,
      }
   aooneint.close()
   return isordk

def readscfinp(file="AOONEINT"):
   aooneint=unformatted.FortranBinary(file)
   table2=aooneint.find("SCFINP")
   #print table2
   scfinp={}
   scfinp["table"]=table2
   aooneint.readrec()
   if len(aooneint.data) == 144:
      title=aooneint.data
      aooneint.readrec()
   else:
      title=aooneint.readbuf(192,'c')
   scfinp["ttitle"]=title
   nsym=aooneint.readbuf(1,'i')[0]
   scfinp["nsym"]=nsym
   scfinp["naos"]=aooneint.readbuf(nsym,'i')
   scfinp["potnuc"]=aooneint.readbuf(1,'d')[0]
   kmax=aooneint.readbuf(1,'i')[0]
   scfinp["kmax"]=kmax
   scfinp["ncent"]=aooneint.readbuf(kmax,'i')
   nbasis=aooneint.readbuf(1,'i')[0]
   scfinp["nbasis"]=nbasis
   scfinp["jtran"]=aooneint.readbuf(nbasis,'i')
   scfinp["itran"]=aooneint.readbuf(8*nbasis,'i')
   scfinp["ctran"]=aooneint.readbuf(8*nbasis,'d')
   scfinp["nbasis"]=aooneint.readbuf(1,'i')[0]
   scfinp["inamn"]=aooneint.readbuf(nbasis,'i')
   scfinp["iptyp"]=aooneint.readbuf(nbasis,'i')
   scfinp["dpnuc"]=aooneint.readbuf(3,'d')
   nucdep=aooneint.readbuf(1,'i')[0]
   scfinp["nucdep"]=nucdep
   scfinp["cooo"]=aooneint.readbuf(3*nucdep,'d')
   scfinp["ifxyz"]=aooneint.readbuf(3,'i')
   scfinp["dummy"]=aooneint.readbuf(1,'d')[0]
   scfinp["qpol"]=aooneint.readbuf(6,'d')
   scfinp["qq"]=aooneint.readbuf(3,'d')
   scfinp["jfxyz"]=aooneint.readbuf(3,'i')
   #print "scfinp",scfinp; sys.exit(0)
   aooneint.close()
   return scfinp

def read(label="OVERLAP",filename="AOONEINT"):
   lbuf=600
   from util import full,blocked
   #
   # Initialize
   #
   unlabeled=readhead(filename)
   nsym=unlabeled["nsym"]
   nbas=unlabeled["naos"]
   nnbast=0
   for isym in range(nsym):
      nnbast+=nbas[isym]*(nbas[isym]+1)/2
   s=full.matrix((nnbast,))
   #
   # Open file, locate label
   #
   aooneint=unformatted.FortranBinary(filename)
   labinfo=aooneint.find(label)
   #
   # Loop over records
   #
   while True:
      try:
         aooneint.readrec()
         buf=aooneint.readbuf(lbuf,'d')
         ibuf=aooneint.readbuf(lbuf,'i')
         length=aooneint.readbuf(1,'i')[0]
         if length < 0: raise Exception("EOD")
         for i in range(length):
            ind=ibuf[i]-1
            s[ind]=buf[i]
      except Exception, inst:
         if inst[0] == "EOD":
            break
         else:
            raise
   S=blocked.triangular(nbas)
   off=0
   for isym in range(nsym):
      if 1:
         nbasi=nbas[isym]*(nbas[isym]+1)/2
         S.subblock[isym]=numpy.array(s[off:off+nbasi]).view(full.triangular)
      else:
         for i in range(nbas[isym]):
            for j in range(i+1):
               S.subblock[isym][i,j]=s[off]
               off+=1
   #aooneint.close()
   return S


if __name__ == "__main__":
   import sys,getopt
   from util.timing import timing
   gethead=False; getisordk=False; getscfinp=False; verbose=False;
   unpack=False;
   try:
      opt,arg=getopt.getopt(sys.argv[1:],'hisvu',['head','isordk','scfinp','verbose','unpack'])
      for o,v in opt:
         if o in ('-h','--head'):
            gethead=True
         if o in ('-i','--isordk'):
            getisordk=True
         if o in ('-s','--scfinp'):
            getscfinp=True
         if o in ('-v','--verbose'):
            verbose=True
         if o in ('-u','--unpack'):
            unpack=True
   except IndexError:
      print "Usage: %s integral" % sys.argv[0]
      sys.exit(1)

   if gethead: 
      t=timing('head')
      head=readhead()
      if head:
         print "Header on AOONEINT"
         for i in head.keys():
            print i,head[i]
      print t

   if getisordk:
      t=timing('getisrordk')
      isordk=readisordk()
      print "isordk table",isordk["table"]
      n=isordk["nucdep"]
      chrn=isordk["chrn"]
      cooo=isordk["cooo"]
      mxcent=len(chrn)
      print "nucdep=%i"%n,"mxcent=%i"%mxcent
      for i in range(n):
         print "Z[%i]=%g"%(i,chrn[i]),"r[%i]=(%g,%g,%g)"%(i,cooo[i],cooo[i+mxcent],cooo[i+mxcent*2])

      print t

   if getscfinp:
      t=timing('scfinp')
      scfinp=readscfinp()
      for i in scfinp.keys():
         print i,scfinp[i]
      print t

   for label in arg:
      tread=timing("read")
      s1=read(label)
      print tread
      if unpack:
         t=timing("unpack")
         s2=s1.unpack()
         print t
         t=timing("unblock")
         S=s2.unblock()
         print t
      else:
         S=s1
      if verbose:
         print label,S
