DEBUG=0
def orbgrad(fc,fo,dc,do,S):
   import one,two
   g= (fc*dc+fo*do)*S - S*(dc*fc+do*fo)
   return g.T

def delta(S):
   if S == 1:
      return 1
   else:
      return 0

def e2n(N,hfx=1,Sg=1,Sv=1):
    """
    E[2]*N linear transformation:
     [k,F] = [k(pq)E(pq), F(rs)E(rs)]
           = k(pq)F(rs) (E(ps)d(rq) - E(rq)d(ps))
           = k(pq)F(ps)E(ps) - k(pq)F(rp)E(rq)
           = [k,F](pq)E(pq)
           = kF
    <[q,kF]> = <[E_{pq}, kF(rs)E(rs)]>
             = kF(rs) <E(ps)d(rq) - E(rq)d(ps)>
             = kF(qs)D(ps) - kF(rp)D(rq)
             = [D, kF.T](pq)
    kD = <[k, E(pq)]>
       = <[k(rs) E(rs), E(pq)]>
       = k(rs) (E(rq)d(ps) - E(ps)d(rq))
       = k(rp)D(rq) - k(qs)D(ps)
       = [k.T, D](p,q)
    Fk = F[kD]
    """
    import os
    import sirifc,one,two,dens

    from ConfigParser import SafeConfigParser
    parser = SafeConfigParser()
    parser.read('dalton.ini')
    tmp = parser.get('Dalton','tmp')

    SIRIFC = os.path.join(tmp, 'SIRIFC')
    AOONEINT = os.path.join(tmp, 'AOONEINT')
    AOTWOINT = os.path.join(tmp, 'AOTWOINT')
    LUINDF = os.path.join(tmp, 'LUINDF')

    ifc=sirifc.sirifc(SIRIFC)
    cmo=ifc.cmo.unblock()

    h=one.read('ONEHAMIL', filename=AOONEINT).unblock().unpack()
    S=one.read('OVERLAP',  filename=AOONEINT).unblock().unpack()


    dct,dot = dens.ifc(SIRIFC)
    da = dot + dct/2
    db = dct/2
    #print "da", da, "db", db

    kN = tomat(N, ifc, tmpdir = tmp).T #transpose  = (q, q+) to (q+/q)
    kn = cmo*kN*cmo.T

    dak = (kn.T*S*da - da*S*kn.T)
    dbk = (kn.T*S*db - db*S*kn.T)
    #print "dak",dak,"dbk",dbk


    fa, fb = two.fockab((da, db),  file=AOTWOINT)
    fa += h; fb += h
    #print "electronic energy", .5*(((h+fa)&da) + ((h+fb)&db))
    fka, fkb = two.fockab((dak, dbk), file=AOTWOINT)

    kfa=S*kn*fa - fa*kn*S
    kfb=S*kn*fb - fb*kn*S

    fat = fka + kfa
    fbt = fkb + kfb

    gao = S*(da*fat.T + db*fbt.T) - (fat.T*da + fbt.T*db)*S
    #print "gao",gao
    gm = cmo.T*gao*cmo
    #print "gm",gm

    # sign convention <[q,[k,F]]> = -E[2]*N
    gv = - tovec(gm, LUINDF)

    return gv

def s2n(N,hfx=1,Sg=1,Sv=1):
    """
    S[2]*N linear transformation:
    <[q,k]> = <[E_{ij}, k_{kl}, E_{kl}]> = 
            = k_{kl} [E_{il}d(kj) - E(kj)d(il)]
	    = k(jl)E(il) - k(ki)E(kj)
	    = Dk.T(ij) - k.TD(ij)
	    = [D, k.T](ij)
    """
    import os
    import sirifc,one,two,dens

    from ConfigParser import SafeConfigParser
    parser = SafeConfigParser()
    parser.read('dalton.ini')
    tmp = parser.get('Dalton','tmp')

    SIRIFC = os.path.join(tmp, 'SIRIFC')
    AOONEINT = os.path.join(tmp, 'AOONEINT')
    AOTWOINT = os.path.join(tmp, 'AOTWOINT')
    LUINDF = os.path.join(tmp, 'LUINDF')

    ifc=sirifc.sirifc(SIRIFC)
    cmo=ifc.cmo.unblock()

    S=one.read('OVERLAP',  filename=AOONEINT).unblock().unpack()


    dct,dot = dens.ifc(SIRIFC)
    da = dot + dct/2
    db = dct/2

    kN=tomat(N, ifc, tmpdir = tmp).T
    kn = cmo*kN*cmo.T

    dak = (kn.T*S*da - da*S*kn.T)
    dbk = (kn.T*S*db - db*S*kn.T)

    gv = -tovec(cmo.T*S*(dak+dbk)*S*cmo, LUINDF)

    return gv


def tovec(k,luindf="LUINDF"):
   from util import full
   wop=jwop(luindf=luindf)
   lwop=len(wop)
   kzyvar=2*lwop
   new=full.matrix((kzyvar,))
   for j in range(lwop):
      r=wop[j][0]-1
      s=wop[j][1]-1
      new[j]=k[r,s]
      new[j+lwop]=k[s,r]
   return new

def tomat(N,ifc,tmpdir='/tmp'):
   import os
   from util import full
   norbt=ifc.norbt
   new=full.matrix((norbt,norbt))
   wop=jwop(luindf=os.path.join(tmpdir,"LUINDF"))
   lwop=len(wop)
   #print "wop.shape",N.shape,lwop
   #print "wop",wop
   for j in range(lwop):
      #print "j wop",j,wop[j]
      r=wop[j][0]-1
      s=wop[j][1]-1
      #print "r s",r,s
      new[r,s]=N[j]
      new[s,r]=N[j+lwop]
   return new

def jwop(luindf="LUINDF"):
   from util import unformatted
   luindf=unformatted.fortranbinary(luindf)
   import one
   table=luindf.find("EXOPSYM1")
   #print "table",table
   luindf.readrec()
   ints=luindf.readbuf(8,'i')
   nwopt=ints[0]
   nwoph=ints[1]
   jwopsy=ints[2]
   nklwop=ints[3]
   #print "sirset",nwopt,nwoph,jwopsy,nklwop
   kzy=2*nwopt
   #print "kzy",kzy
   luindf.readrec()
   ints=luindf.readbuf(kzy,'i')
   wop=[]
   for i in range(0,kzy,2):
      wop.append((ints[i],ints[i+1]))
   return wop
      
   return ints

if __name__ == "__main__":
    import os
    import sirifc
    from util import full
    
    from ConfigParser import SafeConfigParser
    parser = SafeConfigParser()
    parser.read('dalton.ini')
    tmp = parser.get('Dalton','tmp')

    SIRIFC = os.path.join(tmp, 'SIRIFC')
    LUINDF = os.path.join(tmp, 'LUINDF')
    
    ifc=sirifc.sirifc(name=SIRIFC)
    wop=jwop(LUINDF)
    #print "wop",wop
    kzywop=2*len(wop)
    #print "kzywop",kzywop
    #
    # Unit trial vectors as in ABCHK
    #
    #N=full.matrix(kzywop,1).random()
    #X= E2N(N)
    #x= e2n(N)
    #print X,x
    #print "N",N
    #kappa=tomat(N,ifc)
    #print kappa, tovec(kappa,ifc)
    #print N-tovec(kappa,ifc)
    N=full.unit(kzywop)
    E2=full.matrix((kzywop,kzywop))
    S2=full.matrix((kzywop,kzywop))
    Est=full.matrix((kzywop,kzywop))
    Ets=full.matrix((kzywop,kzywop))
    for i in range(kzywop):
        n=N[:,i]
        _e2n = e2n(n)
	_s2n = s2n(n)
	#print n, 'e2', _e2n, 's2', _s2n
        E2[:,i]= _e2n
        S2[:,i]= _s2n
        #Est[:,i]=e2n(n,hfx=1.0,Sg=1,Sv=-1)
        #Ets[:,i]=e2n(n,hfx=1.0,Sg=-1,Sv=1)
    #print "E2",Est,Ets, Est-Ets.T
    #print "jwop",jwop()
    print E2, S2
    import shelve
    database = shelve.open("oli.ref")
    database["E2"] = E2
    database["S2"] = S2
    database.close()
      
   


