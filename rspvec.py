#!/usr/bin/env python
import struct
def read(property,propfile="RSPVEC",timing=False):
   import string,time,numpy
   from util import full,unformatted
   if timing:
      t0=time.clock()
   rspvec=unformatted.fortranbinary(propfile)
   rspvec.find(property)
   rspvec.readrec()
   kzyvar = rspvec.reclen / 8
   buffer_=rspvec.readbuf(kzyvar,'d')
   mat=numpy.array(buffer_).view(full.matrix)
   if timing:
      t1=time.clock()
      print "Time used in rspvec.read: %g"%(t1-t0)
   return mat


if __name__ == "__main__":
   import sys,string
   if len(sys.argv) != 4:
      print "Usage: %s ncdim property file" % sys.argv[0]
      sys.exit(1)
   kzyvar=string.atoi(sys.argv[1])
   prop=sys.argv[2]
   file=sys.argv[3]
   rvec=read(kzyvar,prop,file,timing=True)
   print rvec
