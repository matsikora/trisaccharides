import sys,os
import numpy as np
from subprocess import call
import subprocess
#TODO
# 0. add sub groups to pairs for prot, sugar, intermediate
# 1. fivetimes_protein
# 2. know what is prot what is not
class topology:
   
   def __init__(self,topfile):
      self.aanames= ['---','ala','asn','asp','arg','cys','gln','glu','gly','his','ile','leu','lys','met','pro','phe','ser','thr','trp','tyr','val','unk','ols','olt',
           'ALA','ASN','ASP','ARG','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PRO','PHE','SER','THR','TRP','TYR','VAL','UNK','OLS','OLT','OME']
      self.topfile=topfile
      self.entries=['[ atomtypes ]','[ moleculetype ]','[ atoms ]','[ bonds ]','[ pairs ]','[ angles ]','[ dihedrals ] ; propers','[ dihedrals ] ; impropers','[ system ]','[ molecules ]']
      self.onlyprotein=False # Flags to mark the pure system
      self.onlysugar=False
      prototypes=self.read(self.entries)
      self.comments,self.values,self.valuecomments={},{},{}
      # separate comments from values etc...
      for entry in self.entries:
         self.comments[entry],self.values[entry],self.valuecomments[entry]=self.polish_entry(prototypes[entry])
      # select which atoms belong to protein and sugar based on AA names.
      #~ print self.valuecomments
      
      self.assign_atom_origin()
      self.split_to_groups()
      #~ print self.values
      
   def read(self,entries):
      prototypes={}
      with open(self.topfile,'r') as f:
         buf=f.read()
         for entry in entries:
            #~ print entry
            #~ quit()
            prototypes[entry]=buf.split(entry)[1].split('[ ')[0]
         #~ p_atomtypes=buf.split('[ atomtypes ]')[1].split('[ ')[0]
         #~ p_moleculetype=buf.split('[ moleculetype ]')[1].split('[ ')[0]
         #~ p_atoms=buf.split('[ atoms ]')[1].split('[ ')[0]
         #~ p_bonds=buf.split('[ bonds ]')[1].split('[ ')[0]
         #~ p_pairs=buf.split('[ pairs ]')[1].split('[ ')[0]
         #~ p_angles=buf.split('[ angles ]')[1].split('[ ')[0]
         #~ p_dihedrals_prop=buf.split('[ dihedrals ] ; propers')[1].split('[ ')[0]
         #~ p_dihedrals_improp=buf.split('[ dihedrals ] ; impropers')[1].split('[ ')[0]
         #~ p_system=buf.split('[ system ]')[1].split('[ ')[0]
         #~ p_molecules=buf.split('[ molecules ]')[1].split('[ ')[0]
         #~ print self.polish_entry(p_atoms)
      return prototypes
   def polish_entry(self,entry):
      """
      split comments from values
      returns comments,values,valuecomments (when line ends with ";" and is 
      followed by comment)
      """
      comments=[]
      values=[]
      valuecomments=[]
      for line in entry.split('\n'):
         if len(line)>0:
            if ";" in line.split()[0]:
               comments.append(line)
            else:
               if (";") in line:
                  v=line.split(";")[0]
                  vc=line.split(";")[1]
               else:
                  v=line
                  vc=""
               values.append(v)
               valuecomments.append(vc)
               
      return comments,values,valuecomments
      
   def sugarize_atomtypes(self):
      """
      Add X before each atom definition.
      """
      new=[]
      for line in self.values['[ atomtypes ]']:
         ll=line.split()
         new.append('     '.join(["X"+ll[i] if i<2 else ll[i] for i in xrange(len(ll))]))
      self.values['[ atomtypes ]']=new
      
   def sugarize_atoms(self):
      """
      Add X before each atom name
      """
      new=[]
      for line in self.values['[ atoms ]']:
         ll=line.split()
         new.append('     '.join(["X"+ll[i] if i==1 else ll[i] for i in xrange(len(ll))]))
      self.values['[ atoms ]']=new
      
   def write(self,fname):
      """
      merge and write new topology
      """
      f=open(fname,'w')
      for entry in self.entries:
         f.write(entry+"\n")
         f.writelines([line+"\n" for line in self.comments[entry]])
         for line in xrange(len(self.values[entry])):
            txt=self.values[entry][line]+";"+self.valuecomments[entry][line]+"\n"
            f.write(txt)

      f.close()

   def sixtimes_sugar(self,zeroparams=True,pairtype='[ pairs ] ; sugar'):
      """
      place "[ pairs ]" entry for sugars 6x
      if zeroparams is true, put zeros instead of sigma and epsilon for the copies
      """
      
      new=[]
      new_c=[]
      pairs=self.values[pairtype]
      #~ print self.values.keys()
      
      pairs_c=self.valuecomments[pairtype]
      #~ new+=[pair for pair in pairs]
      #~ print pairtype,len(pairs),"!!"
      
      if len(pairs[0].split()) <5: # sanity check - if sigmas and epsilons were provided:
         raise BaseException("Provide sigma and epsilon values first!")
      for i in xrange(6):
         if zeroparams and i>0:
            new+=['     '.join(line.split()[:3]+['0.0','0.0']) for line in pairs] # zero epsilon and sigma.
         else:
            new+=[line for line in pairs]
         new_c+=[line for line in pairs_c]
         new_c[-1]+="\n;end of %d copy"%(i+1)
      #~ print new
      self.values[pairtype]=new
      self.valuecomments[pairtype]=new_c
   def fivetimes_protein(self,zeroparams=True,pairtype='[ pairs ] ; protein'):
      """
      place "[ pairs ]" entry for protein 5x
      if zeroparams is true, put zeros instead of sigma and epsilon for the copies
      """
      
      new=[]
      new_c=[]
      pairs=self.values[pairtype]
      pairs_c=self.valuecomments[pairtype]
      #~ new+=[pair for pair in pairs]
      if len(pairs[0].split()) <5: # sanity check - if sigmas and epsilons were provided:
         raise BaseException("Provide sigma and epsilon values first!")
      for i in xrange(5):
         if zeroparams and i>0:
            new+=['     '.join(line.split()[:3]+['0.0','0.0']) for line in pairs] # zero epsilon and sigma.
         else:
            new+=[line for line in pairs]
         new_c+=[line for line in pairs_c]
         new_c[-1]+="\n;end of %d copy"%(i+1)
      #~ print new
      self.values[pairtype]=new
      self.valuecomments[pairtype]=new_c
   
   def get_sigma_epsilon(self,pairtype='[ pairs ]',epsilonfactor=1.0):
      """
      calculates sigma and epsilon for 1-4 interaction based on the info stored in atomtypes
      the calculation assumes combination rule 2:
      sigma=1/2(sigma1+sigma2), epsilon=sqrt(epsilon1*epsilon2)
      """
      new=[]
      atoms=self.values['[ atoms ]']
      atomtypes=self.values['[ atomtypes ]']
      pairs=self.values[pairtype]
      
      atoms=np.array([line.split() for line in atoms])
      atomtypes=np.array([line.split() for line in atomtypes])

      
      for line in pairs:
         i,j=line.split()[:2]
         iname,jname=atoms[np.where(atoms[:,0]==i)][0][1],atoms[np.where(atoms[:,0]==j)][0][1]
         isigma,iepsilon=atomtypes[np.where(atomtypes[:,0]==iname)][0][5:7]
         jsigma,jepsilon=atomtypes[np.where(atomtypes[:,0]==jname)][0][5:7]
         newsigma=0.5*(float(isigma)+float(jsigma))
         newepsilon=epsilonfactor*np.sqrt(float(iepsilon)*float(jepsilon))
         new.append('     '.join([line,"%8.6f"%newsigma,"%8.6f"%newepsilon]))
         #~ print iname,jname,isigma,jsigma,newsigma,newepsilon
      #~ print new
      self.values[pairtype]=new
      
   def split_to_groups(self):
      """
      split pairs to groups, depending if they belong to sugar or protein.
      remove old instance of pairs
      """
      prot_like_pairs=[('CA','C1'),('HB2','C1'),('HB3','C1')]
      sugar_like_pairs=[('CB','H1'),('CB','C2'),('CB','O5'),('OG','H2'),('OG','O2'),('OG','C3'),('OG','C5')]

      pairs_sugar=[]
      pairs_protein=[]
      pairs_mixed_as_sugar=[]
      pairs_mixed_as_protein=[]
      
      valuecomments_sugar=[]
      valuecomments_protein=[]
      valuecomments_mixed_as_sugar=[]
      valuecomments_mixed_as_protein=[]
      
      pairs=self.values['[ pairs ]']
      valuecomments=self.valuecomments['[ pairs ]']
      atoms=self.values['[ atoms ]']
      #~ atomtypes=self.values['[ atomtypes ]']  
      atoms=np.array([line.split() for line in atoms])
      #~ atomtypes=np.array([line.split() for line in atomtypes])
      
      for npair in xrange(len(pairs)):
         pair=pairs[npair]
         valuecomment=valuecomments[npair]
         i=pair.split()[0]
         j=pair.split()[1]
         if i in self.sugaratoms and j in self.sugaratoms:
            pairs_sugar.append(pair)
            valuecomments_sugar.append(valuecomment)
         elif ( i in self.sugaratoms and j in self.proteinatoms) or ( j in self.sugaratoms and i in self.proteinatoms):
            # special case!
            iname,jname=atoms[np.where(atoms[:,0]==i)][0][1],atoms[np.where(atoms[:,0]==j)][0][1]
            pdbiname,pdbjname=atoms[np.where(atoms[:,0]==i)][0][4],atoms[np.where(atoms[:,0]==j)][0][4]
            #~ print i,j,iname,jname,pdbiname,pdbjname
            #~ ,atomtypes[np.where(atomtypes[:,0]==iname)][0],atomtypes[np.where(atomtypes[:,0]==jname)][0]
            pairs_mixed_as_sugar.append(pair) # For now all intermediate atoms are treated as sugar
            valuecomments_mixed_as_sugar.append(valuecomment)
         else:
            pairs_protein.append(pair)
            valuecomments_protein.append(valuecomment)
      
      if len(pairs_mixed_as_sugar)==0: # We have a clean system, either protein or sugar
         if len(pairs_sugar)==0:
            self.onlyprotein=True
         elif len(pairs_protein)==0:
            self.onlysugar=True
      pairposition=self.entries.index('[ pairs ]')
      self.entries.insert(pairposition,'[ pairs ] ; protein')
      self.entries.insert(pairposition,'[ pairs ] ; protein')
      self.entries.insert(pairposition,'[ pairs ] ; sugar')
      self.entries.insert(pairposition,'[ pairs ] ; mixed_as_sugar')
      
      self.comments['[ pairs ] ; protein']=[]
      self.comments['[ pairs ] ; sugar']=[]
      self.comments['[ pairs ] ; mixed_as_sugar']=[]
      
      self.values['[ pairs ] ; protein']=pairs_protein
      self.values['[ pairs ] ; sugar']=pairs_sugar
      self.values['[ pairs ] ; mixed_as_sugar']=pairs_mixed_as_sugar

      self.valuecomments['[ pairs ] ; protein']=valuecomments_protein
      self.valuecomments['[ pairs ] ; sugar']=valuecomments_sugar
      self.valuecomments['[ pairs ] ; mixed_as_sugar']=valuecomments_mixed_as_sugar
      
      self.comments['[ pairs ] ; protein']=self.comments['[ pairs ]']
      self.comments['[ pairs ] ; sugar']=self.comments['[ pairs ]']
      self.comments['[ pairs ] ; mixed_as_sugar']=self.comments['[ pairs ]']
      self.entries.remove('[ pairs ]')
      
      
      
   def assign_atom_origin(self):
      """
      read ATOMS and compare to amino acid list to split sugar from protein
      """
      self.proteinatoms=[]
      self.sugaratoms=[]
      
      atoms=self.values['[ atoms ]']
      #~ self.proteinresidues=[int(i.split()[0]) for i in atoms if i.split()[3] in self.aanames]
      for line in atoms:
         l=line.split()
         if l[3] in self.aanames:
            self.proteinatoms.append(l[0])
         else:
            self.sugaratoms.append(l[0])
   def sugarize(self):
      """
      perform actions to introduce proper 1-4 scaling for sugars.
      """
      if not self.onlyprotein:
         self.sugarize_atomtypes()
         self.sugarize_atoms()
         self.get_sigma_epsilon(pairtype='[ pairs ] ; sugar',epsilonfactor=1.0)
         self.sixtimes_sugar(zeroparams=True,pairtype='[ pairs ] ; sugar')
         if not self.onlysugar: # if the system is pure, no mixed pairs
            self.get_sigma_epsilon(pairtype='[ pairs ] ; mixed_as_sugar',epsilonfactor=1.0)
            self.sixtimes_sugar(zeroparams=True,pairtype='[ pairs ] ; mixed_as_sugar')
      return
   def proteinize(self):
      """
      perform actions to introduce proper 1-4 scaling for proteins.
      """
      if not self.onlysugar:
         self.get_sigma_epsilon(pairtype='[ pairs ] ; protein',epsilonfactor=0.5)
         self.fivetimes_protein(zeroparams=True,pairtype='[ pairs ] ; protein')
      return
   def remove_atom(self,atom_number):
      """
      remove an atom from an existing topology
      """
      for entry in self.entries:
         #~ if self.onlysugar:
            #~ if 'protein' in entry or 'mixed' in entry:
               #~ continue
         #~ print entry
         val=self.values[entry]
         com=self.valuecomments[entry]
         line_numbers=[0,0,0,0]
         if 'atoms' in entry:
            line_numbers=[0,5,0,0]
         elif 'bonds' in entry:
            line_numbers=[0,1,0,0]
         elif 'dihedrals' in entry:
            line_numbers=[0,1,2,3]
         elif 'angles' in entry:
            line_numbers=[0,1,2,0]
         elif 'pairs' in entry:
            line_numbers=[0,1,0,0]
         elif 'bonds' in entry:
            line_numbers=[0,1,0,0]
         else:
            continue
         val2=np.array([valline.split() for valline in val])
         #~ print entry
         #~ quit()
         newval=[]
         #~ print len(val2),entry
         if len(val2)>0: # Skip empty entries
            idx=np.where((val2[:,line_numbers[0]]==`atom_number`) | ((val2[:,line_numbers[1]]==`atom_number`)) | ((val2[:,line_numbers[2]]==`atom_number`)) | ((val2[:,line_numbers[3]]==`atom_number`)))[0]
            for line in val2:
               nr1=int(line[line_numbers[0]])
               nr2=int(line[line_numbers[1]])
               nr3=int(line[line_numbers[2]])
               nr4=int(line[line_numbers[3]])
               if nr1>atom_number:
                  line[line_numbers[0]]=str(nr1-1)
               if nr2>atom_number and line_numbers[1]!=line_numbers[0]:
                  line[line_numbers[1]]=str(nr2-1)
               if nr3>atom_number and line_numbers[2]!=line_numbers[0]:
                  line[line_numbers[2]]=str(nr3-1)
               if nr4>atom_number and line_numbers[3]!=line_numbers[0]:
                  line[line_numbers[3]]=str(nr4-1)
               newval.append(' '.join(line))
            self.values[entry]=newval
            for index in sorted(idx, reverse=True):
               del self.valuecomments[entry][index]
               del self.values[entry][index]
   def remove_atoms(self,atom_numbers):
      """
      remove an atom from an existing topology
      """
      atom_numbers=sorted(atom_numbers)
      for iaa in xrange(len(atom_numbers)):
         atom_number=atom_numbers[iaa]
         for entry in self.entries:

            val=self.values[entry]
            com=self.valuecomments[entry]
            line_numbers=[0,0,0,0]
            if 'atoms' in entry:
               line_numbers=[0,5,0,0]
            elif 'bonds' in entry:
               line_numbers=[0,1,0,0]
            elif 'dihedrals' in entry:
               line_numbers=[0,1,2,3]
            elif 'angles' in entry:
               line_numbers=[0,1,2,0]
            elif 'pairs' in entry:
               line_numbers=[0,1,0,0]
            elif 'bonds' in entry:
               line_numbers=[0,1,0,0]
            else:
               continue
            val2=np.array([valline.split() for valline in val])
            newval=[]
            #~ print len(val2),entry
            if len(val2)>0: # Skip empty entries
               idx=np.where((val2[:,line_numbers[0]]==`atom_number`) | ((val2[:,line_numbers[1]]==`atom_number`)) | ((val2[:,line_numbers[2]]==`atom_number`)) | ((val2[:,line_numbers[3]]==`atom_number`)))[0]
               for line in val2:
                  nr1=int(line[line_numbers[0]])
                  nr2=int(line[line_numbers[1]])
                  nr3=int(line[line_numbers[2]])
                  nr4=int(line[line_numbers[3]])
                  if nr1>atom_number:
                     line[line_numbers[0]]=str(nr1-1)
                  if nr2>atom_number and line_numbers[1]!=line_numbers[0]:
                     line[line_numbers[1]]=str(nr2-1)
                  if nr3>atom_number and line_numbers[2]!=line_numbers[0]:
                     line[line_numbers[2]]=str(nr3-1)
                  if nr4>atom_number and line_numbers[3]!=line_numbers[0]:
                     line[line_numbers[3]]=str(nr4-1)
                  newval.append(' '.join(line))
               self.values[entry]=newval
               for index in sorted(idx, reverse=True):
                  del self.valuecomments[entry][index]
                  del self.values[entry][index]
         atom_numbers=[i-1 for i in atom_numbers]
         #~ print atom_numbers
         #~ quit()
def remove_atoms_from_GRO(grofile,atomnumbers,outfile):
   """
   remove atom from a gro file
   """
   new=open(outfile,'w')
   atomnumbers=sorted(atomnumbers) # WAS WRONG, CORRECTED
   #~ newfile=[]
   #~ print 'aaa'
   cutoff=0
   with open(grofile,'r') as f:
      buf=f.readlines()
      for line in buf:
         l=line.split()
         #~ print 'a',l
         if len(l)==7:
            nr1=int(l[3])
            
            if nr1 in atomnumbers:
               #~ print atomnumbers.index(nr1)
               cutoff+=1
               continue
            else:
               
               l[3]=str(nr1-cutoff)
               #~ print l
            #~ new.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"%(int(l[0])))
            new.write("%5d  %-3s%5s%5d%8.3f%8.3f%8.3f\n"%(int(l[0]),l[1],l[2],int(l[3]),float(l[4]),float(l[5]),float(l[6])))
         elif len(l)==1:
            new.write("%s\n"%str(int(l[0])-len(atomnumbers)))
            #~ ll=' '.join(int(l)-1)
         else:
            
            ll=' '.join(l)
            new.write(ll+"\n")
   new.close()         
def remove_atom_from_GRO(grofile,atomnumber,outfile):
   """
   remove atom from a gro file
   """
   new=open(outfile,'w')
   #~ newfile=[]
   #~ print 'aaa'
   with open(grofile,'r') as f:
      buf=f.readlines()
      
      for line in buf:
         l=line.split()
         if len(l)==7:
            nr1=int(l[3])
            if nr1==atomnumber:
               continue
            else:
               l[3]=str(nr1-1)
               #~ print l
            #~ new.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"%(int(l[0])))
            new.write("%5d  %-3s%5s%5d%8.3f%8.3f%8.3f\n"%(int(l[0]),l[1],l[2],int(l[3]),float(l[4]),float(l[5]),float(l[6])))
         elif len(l)==1:
            new.write("%s\n"%str(int(l[0])-l))
            #~ ll=' '.join(int(l)-1)
         else:
            
            ll=' '.join(l)
            new.write(ll+"\n")
   new.close()
               
crdfile='LacNAc_4_15_2014_684_MIN.crd'
topfile='LacNAc_4_15_2014_684.top'
corfile='LacNAc_4_15_2014_684'
gmxgrofile='LacNAc_4_15_2014_684_GMX.gro'
gmxtopfile='LacNAc_4_15_2014_684_GMX.top'
import getopt
#~ opts, extraparams = getopt.getopt(sys.argv[1:],"t:s:",["topology","structure"])
deleteatoms=[]
opts,parms=getopt.getopt(sys.argv[1:],"t:s:o:d:",["topology=","structure=","output=","delete="]) 
for o,p in opts:
   if o in ['-t']:
      #~ print o,p
      gmxtopfile = p
   elif o in ['-s']:
      #~ print o,p
      gmxstructfile= p
   elif o in ['-d','delete']:
      deleteatoms=[int(i) for i in p.split(",")]
   elif o in ['-o','output']:
      core=p
#~ print gmxtopfile,gmxstructfile
#~ print deleteatoms
#~ quit()
#~ gmxtopfile=sys.argv[1]
   
t=topology(gmxtopfile)
#~ print deleteatoms
if len(deleteatoms)>0:
   #~ deleteatoms=[1,2,150]
   t.remove_atoms(deleteatoms)
   remove_atoms_from_GRO(gmxstructfile,deleteatoms,core+'.gro')
   #~ for at in deleteatoms:
      #~ t.remove_atom(at)
      #~ remove_atom_from_GRO(gmxstructfile,at,core+'.gro')
t.sugarize()
t.proteinize()
t.write(core+'.top')
