ll=('s','p','d','f','g','h','i')
def readin(inputfile):
   import sys
   import string
#
# Read inputfile into list lines
#
   molinp=open(inputfile,'r')
   lines=molinp.readlines()
#
# Check if input format is INTGRL, remove all lines before that
#
   start=lines.index('INTGRL\n')
   if start > 0: del lines[0:start]
#
# Two comment lines
#
   comments=lines[1:3]
#
# Symmetry info extracted form this line
#
   symline=lines[3]
   atom_types=string.atoi(symline[1:5])
   charge=string.atof(symline[5:8])
   symops=string.atoi(symline[8:10])
   symop=string.split(symline[10:19])
#
# The molecule is a list of atom types which are hashes with structure
#    atom{'charge'} charge of this atom type
#    atom{'label'}  list of atom label strings of this type
#    atom{'center'} list of atom label strings of this type
#    atom{'basis '} structure of form
#                   [angmom][block]{'prim':[], 'cont':[][]}
#
   molecule=[]
   for atype in range(atom_types):
      atom={}
      molecule.append(atom)
      atline=lines.pop(4)
   #print atline
      atfield=string.split(atline)
#
#   Information for this atom type
#
      atom['charge']=string.atof(atfield[0])
      unique =string.atoi(atfield[1])
      maxl=string.atoi(atfield[2])
      #print 'maxl',maxl
      atom['label'] =[]
      atom['center'] =[]
      atom['basis'] =[]
   #print atom
#
#   Information for all atoms of this atom type
#
      for at in range(unique):
         line=lines.pop(4)
         field=string.split(line)
         atom['label'].append([])
         atom['label'][at]=field[0]
         atom['center'].append([])
         for i in range(1,4): atom['center'][at].append(string.atof(field[i]))
#
#   Information for the basis set of this atom type
#
      for l in range(maxl):
         atom['basis'].append([])
      #print atom['basis']
         nbl=string.atoi(atfield[3+l])
         for block in range(nbl):
            atom['basis'][l].append([])
            atom['basis'][l][block]={}
            atom['basis'][l][block]['prim']=[]
            atom['basis'][l][block]['cont']=[]
            line=lines.pop(4)
            field=string.split(line)
            npr=string.atoi(field[1])
            nco=string.atoi(field[2])
            #print npr,nco
            for prim in range(npr):
               line=lines.pop(4)
               if nco > 3:
                  line+=lines.pop(4)
               #print 'prim',prim,line
               field=string.split(line)
               #print "field",field
               atom['basis'][l][block]['prim'].append(string.atof(field[0]))
               atom['basis'][l][block]['cont'].append([])
            #print atom['basis']
               for cont in range(nco):
                  atom['basis'][l][block]['cont'][prim].append(string.atof(field[cont+1]))
   return molecule

def printbasis(molecule):
   i=0
   for atom in molecule:
      i+=1
      print "Atom type %d charge %f" % (i,atom['charge'])
      for at in range(len(atom['label'])):
         print "center",at+1,atom['center'][at]
      for l in range(len(atom['basis'])):
         print "%s-functions" % ll[l]
         for block in range(len(atom['basis'][l])):
            prim=atom['basis'][l][block]['prim']
            cont=atom['basis'][l][block]['cont']
            for ip in range(len(prim)):
               print "   ",prim[ip],cont[ip][:]


#Test nuclear repulsion

def atoms_in(molecule):
   atomlist=[]
   ia=0
   for atype in molecule:
      for a in atype["center"]:
         atomlist.append({})
         Z=atype["charge"]
         atomlist[ia]["charge"]=Z
         atomlist[ia]["center"]=a
         atomlist[ia]["type"]=atype
         contracted=0
         for l in range(len(atype['basis'])):
            contracted_l=0
            for block in range(len(atype['basis'][l])):
               cont=atype['basis'][l][block]['cont']
               contracted+=len(cont[0])*(2*l+1)
               contracted_l+=len(cont[0])*(2*l+1)
            atomlist[ia][ll[l]]=contracted_l
         atomlist[ia]["contracted"]=contracted
         ia+=1
   return atomlist

def contracted_per_atom(molecule):
   atomlist=atoms_in(molecule)
   contracted=[]
   for atom in atomlist:
      contracted.append(atom["contracted"])
   return contracted

def contracted_per_atom_l(molecule):
   atomlist=atoms_in(molecule)
   contracted_l=[]
   for i in range(len(atomlist)):
      #print "i",i
      atom=atomlist[i]
      contracted_l.append([])
      for l in range(len(atom["type"]["basis"])):
         #print "l",l,ll[l]
         contracted_l[i].append(atom[ll[l]])
   return contracted_l

def occupied_per_atom(molecule):
   # assumes ANO type basis
   cpal=contracted_per_atom_l(molecule)
   atomlist=atoms_in(molecule)
   assert len(cpal) == len(atomlist)
   nocclist=[]
   for a in range(len(atomlist)):
      offset=0
      atom=atomlist[a]
      Z=int(atom["charge"]+0.5)
      nocclist.append([])
      for l in range(len(cpal[a])):
         #print "Z a l",Z,a,ll[l]
         if Z <= 1:
            if l == 0:
               nocclist[a].append(offset)
         elif Z<=3:
            if l == 0:
               nocclist[a].append(offset)
               nocclist[a].append(offset+1)
         elif Z<=8:
            if l == 0:
               nocclist[a].append(offset)
               nocclist[a].append(offset+1)
            elif l == 1:
               nocclist[a].append(offset)
               nocclist[a].append(offset+1)
               nocclist[a].append(offset+2)
         else:
            raise Exception( "occupied_per_atom:not implemented")
         offset+=cpal[a][l]
   return nocclist

def nuc(molecule):
   atomlist=atoms_in(molecule)
   Z=0.0
   for a in atomlist:
      ZA=a["charge"]
      RA=a["center"]
      ia=atomlist.index(a)
      for b in atomlist[:ia]:
         ZB=b["charge"]
         RB=b["center"]
         Z+=ZA*ZB/dist(RA,RB)
   return Z

def dist(A,B):
   import math
   d=0
   for i in range(3):
      d+=(A[i]-B[i])*(A[i]-B[i])
   return math.sqrt(d)

def overlap(molecule):
   S=property(molecule,["OVERLAP"])
   return S["OVERLAP"]

def diplen(molecule):
   label=["XDIPLEN","YDIPLEN","ZDIPLEN"]
   S=property(molecule,label)
   return S["XDIPLEN"],S["YDIPLEN"],S["ZDIPLEN"]

def property(molecule,oplabel):
   import mm
   khkt=(1,3,6,10)
   GA={}
   GB={}
   ABprim={}
   AB={}
   atomlist=[]
   if "NUCATTR" in oplabel: atomlist=atoms_in(molecule)
   #print "atoms_in->",atomlist
   for op in oplabel:
      AB[op]=[]
      ABprim[op]=[]
   ia=0
   for A in molecule:
      for a in range(len(A['label'])):
         for la in range(len(A['basis'])):
            ka=khkt[la] 
            for blocka in range(len(A['basis'][la])):
               Aprims=A['basis'][la][blocka]['prim']
               Aconts=A['basis'][la][blocka]['cont']
               for Aprim in Aprims:
                  for  op in oplabel:
                     for kka in range(ka):AB[op].append([])
                  ib=0
                  for B in molecule:
                     for b in range(len(B['label'])):
                        for lb in range(len(B['basis'])):
                           kb=khkt[lb]
                           for blockb in range(len(B['basis'][lb])):
                              Bprims=B['basis'][lb][blockb]['prim']
                              Bconts=B['basis'][lb][blockb]['cont']
                              #print "Aprims - Bprims",Aprims,"-",Bprims
                              for Bprim in Bprims:
                                 GA['center']=A['center'][a][:]
                                 GB['center']=B['center'][b][:]
                                 GA['exponent']=Aprim
                                 GB['exponent']=Bprim
                                 GA['angmom']=la
                                 GB['angmom']=lb
                                 ABprim=mm.matrix(GA,oplabel,GB,atomlist)
                                 for op in oplabel:
                                  for kka in range(ka):
                                    for kkb in range(kb):
                                       AB[op][ia+kka].append(ABprim[op][kka][kkb])
                                 ib+=kb
                  ia+=ka
   return AB

if __name__ == "__main__":
   #mo=readin("MOLECULE.INP")
   def header(title):
      print "\n%s\n%s\n"%(title,"_"*len(title))
   mo=readin("DALTON.BAS")
   if 1:
      header("printbasis")
      printbasis(mo)
   if 1:
      header("atoms_in")
      print atoms_in(mo)
   if 1:
      header("contracted_per_atom")
      print contracted_per_atom(mo)
      header("contracted_per_atom_l")
      print contracted_per_atom_l(mo)
   if 2:
      header("occupied_per_atom")
      print occupied_per_atom(mo)
   if 0:
      header("nuc")
      print nuc(mo)
   if 0:
      header("dist")
      print dist(mo[0]["center"][0],mo[1]["center"][0])
   if 0:
      header("OVERLAP")
      print overlap(mo)
