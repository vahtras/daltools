#!/usr/bin/env python
def read(n,property="OVERLAP",propfile="AOPROPER"):
   import numpy
   from daltools import one
   from util import unformatted,full
   ufile=unformatted.fortranbinary(propfile)
   ufile.find(property)
   line=ufile.data
   stars,data,symtype,label=(
	 line[:8],line[8:16],line[16:24],line[24:32]
	 )
   ufile.readrec()
   shape=(n,n)
   if symtype == "SQUARE  ":
      square=1
      matsize=n*n
      buffer=ufile.readbuf(matsize,'d')
      mat=full.matrix(shape)
   else:
      square=0
      matsize=n*(n+1)/2
      #mat=full.triangular(shape,anti=(symtype == "ANTISYMM"))
      buffer=ufile.readbuf(matsize,'d')
      mat=numpy.array(buffer).view(full.triangular)
      mat.anti = (symtype == "ANTISYMM")
   #print matsize,len(buffer)
   if 0:
      nn=0
      for i in range(n):
	 if square:
	    jtop=n
	 else:
	    jtop=i+1
	 for j in range(jtop):  
	    mat[i,j]=buffer[nn]
	    nn+=1
   return mat

if __name__ == "__main__":
   import sys,sirifc, timing,getopt
   from timing import timing
   t_all=timing("main")
   n=sirifc.sirifc().nbast
   #yo=sys.argv[1]
   unpack=False
   verbose=False
   try:
      opt,arg=getopt.getopt(sys.argv[1:],'uv',['unpack','verbose'])
      for o,v in opt:
	 if o in ('-v','--verbose'):
	    verbose=True
	 if o in ('-u','--unpack'):
	    unpack=True
      label=arg[0]
   except IndexError:
      print "Usage: %s property" % sys.argv[0]
      sys.exit(1)
   t_read=timing("prop.read")
   propint=read(n,property=label)
   print t_read
   if unpack:
      t_unpack=timing("unpack")
      propsq=propint.unpack()
      print t_unpack
   if verbose:
      if unpack:
	 print propsq
      else:
	 print propint
   print t_all
