#!/usr/bin/env python3


# EDITABLE SECTIONS ARE MARKED WITH #@#


version = "1.6 based on martinize v2.6 for Python3"
authors = ["Djurre H. de Jong", "Jaakko J. Uusitalo", "Tsjerk A. Wassenaar", "Philipp S. Schmalhorst", "Mateusz Sikora"]

# This program has grown to be pretty complex.
# The routines have been organized in different files.
# For working versions, all files can be incorporated by using the catenate.py file. 
#
# Index of the program files:
#
#   1. Options and documentation                             @DOC.py
#   2. Description, options and command line parsing         @CMD.py
#   3. Helper functions and macros                           @FUNC.py
#   4. Finegrained to coarsegrained mapping                  @MAP.py
#   5. Secondary structure determination and interpretation  @SS.py
#   6. Elastic network                                       @ELN.py
#   7. Structure I/O                                         @IO.py
#   8. Topology generation                                   @TOP.py
#   9. Main                                                  @MAIN.py
#
#   Force field parameters are specified in the differentt forcefield modules,
#   e.g.: martini22_ff.py


################################
## 6 # FORCE FIELD PARAMETERS ##  -> @FF <-
################################

class martini21(object):
    ff = True
    def __init__(self):

        # parameters are defined here for the following (protein) forcefields:
        self.name = 'martini21'
        
        # Charged types:
        self.charges = {"Qd":1, "Qa":-1, "SQd":1, "SQa":-1, "RQd":1, "AQa":-1}                        #@#
        
        
        #----+---------------------+
        ## A | BACKBONE PARAMETERS |
        #----+---------------------+
        #
        # bbss  lists the one letter secondary structure code
        # bbdef lists the corresponding default backbone beads
        # bbtyp lists the corresponding residue specific backbone beads
        #
        # bbd   lists the structure specific backbone bond lengths
        # bbkb  lists the corresponding bond force constants
        #
        # bba   lists the structure specific angles
        # bbka  lists the corresponding angle force constants
        #
        # bbd   lists the structure specific dihedral angles
        # bbkd  lists the corresponding force constants
        #
        # -=NOTE=- 
        #  if the secondary structure types differ between bonded atoms
        #  the bond is assigned the lowest corresponding force constant 
        #
        # -=NOTE=-
        # if proline is anywhere in the helix, the BBB angle changes for 
        # all residues
        #
        
        ###############################################################################################
        ## BEADS ##                                                                         #                 
        #                          F     E     H     1     2     3     T     S     C        # SS one letter   
        self.bbdef    =    spl(" N0   Nda    N0    Nd    Na   Nda   Nda    P5    P5")  # Default beads   #@#
        self.bbtyp    = {                                                                   #                 #@#
                     "ALA": spl(" C5    N0    C5    N0    N0    N0    N0    P4    P4"),# ALA specific    #@#
                     "PRO": spl(" C5    N0    C5    N0    Na    N0    N0    Na    Na"),# PRO specific    #@#
                     "HYP": spl(" C5    N0    C5    N0    N0    N0    N0    Na    Na") # HYP specific    #@#
        }                                                                                   #                 #@#
        ## BONDS ##                                                                         #                 
        self.bbldef   =             (.365, .350, .350, .350, .350, .350, .350, .350, .350)  # BB bond lengths #@#
        self.bbkb     =             (1250, 1250, 1250, 1250, 1250, 1250,  500,  400,  400)  # BB bond kB      #@#
        self.bbltyp   = {}                                                                  #                 #@#
        self.bbkbtyp  = {}                                                                  #                 #@#
        ## ANGLES ##                                                                        #                 
        self.bbadef   =             ( 119.2,134,   96,   96,   96,   96,  100,  130,  127)  # BBB angles      #@#
        self.bbka     =             ( 150,   25,  700,  700,  700,  700,   25,   25,   25)  # BBB angle kB    #@#
        self.bbatyp   = {                                                                   #                 #@#
               "PRO":               ( 119.2,134,   98,   98,   98,   98,  100,  130,  127), # PRO specific    #@#
               "HYP":               ( 119.2,134,   98,   98,   98,   98,  100,  130,  127)  # PRO specific    #@#
        }                                                                                   #                 #@#
        self.bbkatyp  = {                                                                   #                 #@#
               "PRO":               ( 150,   25,  100,  100,  100,  100,   25,   25,   25), # PRO specific    #@#
               "HYP":               ( 150,   25,  100,  100,  100,  100,   25,   25,   25)  # PRO specific    #@#
        }                                                                                   #                 #@#
        ## DIHEDRALS ##                                                                     #                 
        self.bbddef   =             ( 90.7,   0, -120, -120, -120, -120)                    # BBBB dihedrals  #@#
        self.bbkd     =             ( 100,   10,  400,  400,  400,  400)                    # BBBB kB         #@#
        self.bbdmul   =             (   1,    1,    1,    1,    1,    1)                    # BBBB mltplcty   #@#
        self.bbdtyp   = {}                                                                  #                 #@#
        self.bbkdtyp  = {}                                                                  #                 #@#
                                                                                            #                 
        ###############################################################################################               
        
        # Some Forcefields use the Ca position to position the BB-bead (me like!)
        # martini 2.1 doesn't
        self.ca2bb = False 
        
        # BBS angle, equal for all ss types                                                         
        # Connects BB(i-1),BB(i),SC(i), except for first residue: BB(i+1),BB(i),SC(i)               
        #                 ANGLE   Ka                                                                
        self.bbsangle =      [   100,  25]                                                               #@#
        
        # Bonds for extended structures (more stable than using dihedrals)                          
        #               LENGTH FORCE                                                                
        self.ebonds   = {                                                                                #@#
               'short': [ .640, 2500],                                                              #@#
               'long' : [ .970, 2500]                                                               #@#
        }                                                                                           #@#
        
        
        #----+-----------------------+
        ## B | SIDE CHAIN PARAMETERS |
        #----+-----------------------+
        
        # To be compatible with Elnedyn, all parameters are explicitly defined, even if they are double.
        self.sidechains = {
            #RES#   BEADS                   BONDS                                                   ANGLES              DIHEDRALS
            #                               BB-SC          SC-SC                                        BB-SC-SC  SC-SC-SC
            "TRP": [spl("SC4 SP1 SC4 SC4"),[(0.300,5000)]+[(0.270,None) for i in range(5)],        [(210,50),(90,50),(90,50)], [(0,50),(0,200)]],
            "TYR": [spl("SC4 SC4 SP1"),    [(0.320,5000), (0.270,None), (0.270,None),(0.270,None)],[(150,50),(150,50)],        [(0,50)]],
            "PHE": [spl("SC4 SC4 SC4"),    [(0.310,7500), (0.270,None), (0.270,None),(0.270,None)],[(150,50),(150,50)],        [(0,50)]],
            "HIS": [spl("SC4 SP1 SP1"),    [(0.320,7500), (0.270,None), (0.270,None),(0.270,None)],[(150,50),(150,50)],        [(0,50)]],
            "HIH": [spl("SC4 SP1 SQd"),    [(0.320,7500), (0.270,None), (0.270,None),(0.270,None)],[(150,50),(150,50)],        [(0,50)]],
            "ARG": [spl("N0 Qd"),          [(0.330,5000), (0.340,5000)],                           [(180,25)]],
            "LYS": [spl("C3 Qd"),          [(0.330,5000), (0.280,5000)],                           [(180,25)]],
            "CYS": [spl("C5"),             [(0.310,7500)]],
            "ASP": [spl("Qa"),             [(0.320,7500)]],
            "GLU": [spl("Qa"),             [(0.400,5000)]],
            "ILE": [spl("AC1"),            [(0.310,None)]],
            "LEU": [spl("AC1"),            [(0.330,7500)]],
            "MET": [spl("C5"),             [(0.400,2500)]],
            "ASN": [spl("P5"),             [(0.320,5000)]],
            "PRO": [spl("AC2"),            [(0.300,7500)]],
            "HYP": [spl("P1"),             [(0.300,7500)]],
            "GLN": [spl("P4"),             [(0.400,5000)]],
            "SER": [spl("P1"),             [(0.250,7500)]],
            "THR": [spl("P1"),             [(0.260,None)]],
            "VAL": [spl("AC2"),            [(0.265,None)]],
            "ALA": [],
            "GLY": [],
            }
        
        # Not all (eg Elnedyn) forcefields use backbone-backbone-sidechain angles and BBBB-dihedrals.
        self.UseBBSAngles          = True 
        self.UseBBBBDihedrals      = True

        # Martini 2.2p has polar and charged residues with seperate charges.
        self.polar   = []
        self.charged = []

        # If masses or charged diverge from standard (45/72 and -/+1) they are defined here.
        self.mass_charge = {
        #RES   MASS               CHARGE
        }

        # Defines the connectivity between between beads
        self.connectivity = {
        #RES       BONDS                                   ANGLES             DIHEDRALS              V-SITE
        "TRP":     [[(0,1),(1,2),(1,3),(2,3),(2,4),(3,4)], [(0,1,2),(0,1,3)], [(0,2,3,1),(1,2,4,3)]],
        "TYR":     [[(0,1),(1,2),(1,3),(2,3)],             [(0,1,2),(0,1,3)], [(0,2,3,1)]],
        "PHE":     [[(0,1),(1,2),(1,3),(2,3)],             [(0,1,2),(0,1,3)], [(0,2,3,1)]],
        "HIS":     [[(0,1),(1,2),(1,3),(2,3)],             [(0,1,2),(0,1,3)], [(0,2,3,1)]],
        "HIH":     [[(0,1),(1,2),(1,3),(2,3)],             [(0,1,2),(0,1,3)], [(0,2,3,1)]],
        "GLN":     [[(0,1)]],
        "ASN":     [[(0,1)]],
        "SER":     [[(0,1)]],
        "THR":     [[(0,1)]],
        "ARG":     [[(0,1),(1,2)],                         [(0,1,2)]],
        "LYS":     [[(0,1),(1,2)],                         [(0,1,2)]],
        "ASP":     [[(0,1)]],
        "GLU":     [[(0,1)]],
        "CYS":     [[(0,1)]],
        "ILE":     [[(0,1)]],
        "LEU":     [[(0,1)]],
        "MET":     [[(0,1)]],
        "PRO":     [[(0,1)]],
        "HYP":     [[(0,1)]],
        "VAL":     [[(0,1)]],
        "ALA":     [],
        "GLY":     [],
        }
       
        #----+----------------+
        ## C | SPECIAL BONDS  |
        #----+----------------+
        
        self.special = {
            # Used for sulfur bridges
            # ATOM 1         ATOM 2          BOND LENGTH   FORCE CONSTANT
            (("SC1","CYS"), ("SC1","CYS")):     (0.24,         None),
            }
        
        # By default use an elastic network
        self.ElasticNetwork = False 

        # Elastic networks bond shouldn't lead to exclusions (type 6) 
        # But Elnedyn has been parametrized with type 1.
        self.EBondType = 6
        
        #----+----------------+
        ## D | INTERNAL STUFF |
        #----+----------------+
        
        
        ## BACKBONE BEAD TYPE ##                                                                    
        # Dictionary of default bead types (*D)                                                     
        self.bbBeadDictD  = hash(bbss,self.bbdef)                                                             
        # Dictionary of dictionaries of types for specific residues (*S)                            
        self.bbBeadDictS  = dict([(i,hash(bbss,self.bbtyp[i])) for i in list(self.bbtyp.keys())])                        
        
        ## BB BOND TYPE ##                                                                          
        # Dictionary of default abond types (*D)                                                    
        self.bbBondDictD = hash(bbss,list(zip(self.bbldef,self.bbkb)))                                                   
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbBondDictS = dict([(i,hash(bbss,list(zip(self.bbltyp[i],self.bbkbtyp[i])))) for i in list(self.bbltyp.keys())])       
        # This is tricky to read, but it gives the right bondlength/force constant
        
        ## BBB ANGLE TYPE ##                                                                        
        # Dictionary of default angle types (*D)                                                    
        self.bbAngleDictD = hash(bbss,list(zip(self.bbadef,self.bbka)))                                                  
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbAngleDictS = dict([(i,hash(bbss,list(zip(self.bbatyp[i],self.bbkatyp[i])))) for i in list(self.bbatyp.keys())])      
                    
        ## BBBB DIHEDRAL TYPE ##                                                                    
        # Dictionary of default dihedral types (*D)                                                 
        self.bbDihedDictD = hash(bbss,list(zip(self.bbddef,self.bbkd,self.bbdmul)))                                           
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbDihedDictS = dict([(i,hash(bbss,list(zip(self.bbdtyp[i],self.bbkdtyp[i])))) for i in list(self.bbdtyp.keys())])      
        
    # The following function returns the backbone bead for a given residue and                   
    # secondary structure type.                                                                 
    # 1. Look up the proper dictionary for the residue                                          
    # 2. Get the proper type from it for the secondary structure                                
    # If the residue is not in the dictionary of specials, use the default                      
    # If the secondary structure is not listed (in the residue specific                         
    # dictionary) revert to the default.                                                        
    def bbGetBead(self,r1,ss="C"):                                                                   
        return self.bbBeadDictS.get(r1,self.bbBeadDictD).get(ss,self.bbBeadDictD.get(ss))                      

    def bbGetBond(self,r,a,ss):
        # Retrieve parameters for each residue from table defined above
        b1 = self.bbBondDictS.get(r[0],self.bbBondDictD).get(ss[0],self.bbBondDictD.get(ss[0]))
        b2 = self.bbBondDictS.get(r[1],self.bbBondDictD).get(ss[1],self.bbBondDictD.get(ss[1]))
        # Determine which parameters to use for the bond
        # Python3 no longer returns a NoneType if None is compared to an integer.
        if b1[1] == None or b2[1] == None:
            return ( (b1[0]+b2[0])/2, None )
        else:
            return ( (b1[0]+b2[0])/2, min(b1[1],b2[1]) )
    
    def bbGetAngle(self,r,ca,ss):
        # PRO in helices is dominant
        if r[1] == "PRO" and ss[1] in "H123":
            return self.bbAngleDictS["PRO"].get(ss[1])
        else:
            # Retrieve parameters for each residue from table defined above
            a = [ self.bbAngleDictS.get(r[0],self.bbAngleDictD).get(ss[0],self.bbAngleDictD.get(ss[0])),
                  self.bbAngleDictS.get(r[1],self.bbAngleDictD).get(ss[1],self.bbAngleDictD.get(ss[1])),
                  self.bbAngleDictS.get(r[2],self.bbAngleDictD).get(ss[2],self.bbAngleDictD.get(ss[2])) ]
            # Sort according to force constant
            a.sort(key=lambda i: (i[1],i[0]))
            # This selects the set with the smallest force constant and the smallest angle
            return a[0]
        
    def messages(self):
        import logging 
        '''Prints any force-field specific logging messages.'''
        logging.info('Note: Cysteine bonds are 0.24 nm constraints, instead of the published 0.39nm/5000kJ/mol.')
    
    
################################
## 6 # FORCE FIELD PARAMETERS ##  -> @FF <-
################################

# New martini 2.2 parameters.
# Changed: 
#   Unstructured Pro backbone bead
#   Proline side chains
#   Phe sidechain
#   Trp sidechain
#   Helix BB-bonds to constraint      

class martini22(object):
    ff = True
    def __init__(self):

        # parameters are defined here for the following (protein) forcefields:
        self.name = 'martini22'
        
        # Charged types:
        self.charges = {"Qd":1, "Qa":-1, "SQd":1, "SQa":-1, "RQd":1, "AQa":-1}                                                           #@#
        
        
        #----+---------------------+
        ## A | BACKBONE PARAMETERS |
        #----+---------------------+
        #
        # bbss  lists the one letter secondary structure code
        # bbdef lists the corresponding default backbone beads
        # bbtyp lists the corresponding residue specific backbone beads
        #
        # bbd   lists the structure specific backbone bond lengths
        # bbkb  lists the corresponding bond force constants
        #
        # bba   lists the structure specific angles
        # bbka  lists the corresponding angle force constants
        #
        # bbd   lists the structure specific dihedral angles
        # bbkd  lists the corresponding force constants
        #
        # -=NOTE=- 
        #  if the secondary structure types differ between bonded atoms
        #  the bond is assigned the lowest corresponding force constant 
        #
        # -=NOTE=-
        # if proline is anywhere in the helix, the BBB angle changes for 
        # all residues
        #
        
        ###############################################################################################
        ## BEADS ##                                                                         #                 
        #                              F     E     H     1     2     3     T     S     C    # SS one letter   
        self.bbdef    =    spl(" N0   Nda    N0    Nd    Na   Nda   Nda    P5    P5")  # Default beads   #@#
        self.bbtyp    = {                                                                   #                 #@#
                    "ALA": spl(" C5    N0    C5    N0    N0    N0    N0    P4    P4"), # ALA specific    #@#
                    "PRO": spl(" C5    N0    C5    N0    Na    N0    N0    P4    P4"), # PRO specific    #@#
                    "HYP": spl(" C5    N0    C5    N0    N0    N0    N0    P4    P4")  # HYP specific    #@#
        }                                                                                   #                 #@#
        ## BONDS ##                                                                         #                 
        self.bbldef   =             (.365, .350, .310, .310, .310, .310, .350, .350, .350)  # BB bond lengths #@#
        self.bbkb     =             (1250, 1250, None, None, None, None, 1250, 1250, 1250)  # BB bond kB      #@#
        self.bbltyp   = {}                                                                  #                 #@#
        self.bbkbtyp  = {}                                                                  #                 #@#
        ## ANGLES ##                                                                        #                 
        self.bbadef   =             ( 119.2,134,   96,   96,   96,   96,  100,  130,  127)  # BBB angles      #@#
        self.bbka     =             ( 150,   25,  700,  700,  700,  700,   20,   20,   20)  # BBB angle kB    #@#
        self.bbatyp   = {                                                                   #                 #@#
               "PRO":               ( 119.2,134,   98,   98,   98,   98,  100,  130,  127), # PRO specific    #@#
               "HYP":               ( 119.2,134,   98,   98,   98,   98,  100,  130,  127)  # PRO specific    #@#
        }                                                                                   #                 #@#
        self.bbkatyp  = {                                                                   #                 #@#
               "PRO":               ( 150,   25,  100,  100,  100,  100,   25,   25,   25), # PRO specific    #@#
               "HYP":               ( 150,   25,  100,  100,  100,  100,   25,   25,   25)  # PRO specific    #@#
        }                                                                                   #                 #@#
        ## DIHEDRALS ##                                                                     #                 
        self.bbddef   =             ( 90.7,   0, -120, -120, -120, -120)                    # BBBB dihedrals  #@#
        self.bbkd     =             ( 100,   10,  400,  400,  400,  400)                    # BBBB kB         #@#
        self.bbdmul   =             (   1,    1,    1,    1,    1,    1)                    # BBBB mltplcty   #@#
        self.bbdtyp   = {}                                                                  #                 #@#
        self.bbkdtyp  = {}                                                                  #                 #@#
                                                                                            #                 
        ###############################################################################################               
        
        # Some Forcefields use the Ca position to position the BB-bead (me like!)
        # martini 2.1 doesn't
        self.ca2bb = False 
        
        # BBS angle, equal for all ss types                                                         
        # Connects BB(i-1),BB(i),SC(i), except for first residue: BB(i+1),BB(i),SC(i)               
        #                 ANGLE   Ka                                                                
        self.bbsangle =      [   100,  25]                                                               #@#
        
        # Bonds for extended structures (more stable than using dihedrals)                          
        #               LENGTH FORCE                                                                
        self.ebonds   = {                                                                                #@#
               'short': [ .640, 2500],                                                              #@#
               'long' : [ .970, 2500]                                                               #@#
        }                                                                                           #@#
        
        
        #----+-----------------------+
        ## B | SIDE CHAIN PARAMETERS |
        #----+-----------------------+
        
        # To be compatible with Elnedyn, all parameters are explicitly defined, even if they are double.
        self.sidechains = {
            #RES#   BEADS                   BONDS                                                   ANGLES              DIHEDRALS
            #                               BB-SC          SC-SC                                        BB-SC-SC  SC-SC-SC
            "TRP": [spl("SC4 SNd SC5 SC5"),[(0.300,5000)]+[(0.270,None) for i in range(5)],        [(210,50),(90,50),(90,50)], [(0,50),(0,200)]],
            "TYR": [spl("SC4 SC4 SP1"),    [(0.320,5000), (0.270,None), (0.270,None),(0.270,None)],[(150,50),(150,50)],        [(0,50)]],
            "PHE": [spl("SC5 SC5 SC5"),    [(0.310,7500), (0.270,None), (0.270,None),(0.270,None)],[(150,50),(150,50)],        [(0,50)]],
            "HIS": [spl("SC4 SP1 SP1"),    [(0.320,7500), (0.270,None), (0.270,None),(0.270,None)],[(150,50),(150,50)],        [(0,50)]],
            "HIH": [spl("SC4 SP1 SQd"),    [(0.320,7500), (0.270,None), (0.270,None),(0.270,None)],[(150,50),(150,50)],        [(0,50)]],
            "ARG": [spl("N0 Qd"),          [(0.330,5000), (0.340,5000)],                           [(180,25)]],
            "LYS": [spl("C3 Qd"),          [(0.330,5000), (0.280,5000)],                           [(180,25)]],
            "CYS": [spl("C5"),             [(0.310,7500)]],
            "ASP": [spl("Qa"),             [(0.320,7500)]],
            "GLU": [spl("Qa"),             [(0.400,5000)]],
            "ILE": [spl("AC1"),            [(0.310,None)]],
            "LEU": [spl("AC1"),            [(0.330,7500)]],
            "MET": [spl("C5"),             [(0.400,2500)]],
            "ASN": [spl("P5"),             [(0.320,5000)]],
            "PRO": [spl("C3"),             [(0.300,7500)]],
            "HYP": [spl("P1"),             [(0.300,7500)]],
            "GLN": [spl("P4"),             [(0.400,5000)]],
            "SER": [spl("P1"),             [(0.250,7500)]],
            "THR": [spl("P1"),             [(0.260,None)]],
            "VAL": [spl("AC2"),            [(0.265,None)]],
            "ALA": [],
            "GLY": [],
            }
        
        # Not all (eg Elnedyn) forcefields use backbone-backbone-sidechain angles and BBBB-dihedrals.
        self.UseBBSAngles          = True 
        self.UseBBBBDihedrals      = True

        # Martini 2.2p has polar and charged residues with seperate charges.
        self.polar   = []
        self.charged = []

        # If masses or charged diverge from standard (45/72 and -/+1) they are defined here.
        self.mass_charge = {
        #RES   MASS               CHARGE
        }

        # Defines the connectivity between between beads
        self.connectivity = {
        #RES       BONDS                                   ANGLES             DIHEDRALS              V-SITE
        "TRP":     [[(0,1),(1,2),(1,3),(2,3),(2,4),(3,4)], [(0,1,2),(0,1,3)], [(0,2,3,1),(1,2,4,3)]],  
        "TYR":     [[(0,1),(1,2),(1,3),(2,3)],             [(0,1,2),(0,1,3)], [(0,2,3,1)]], 
        "PHE":     [[(0,1),(1,2),(1,3),(2,3)],             [(0,1,2),(0,1,3)], [(0,2,3,1)]],
        "HIS":     [[(0,1),(1,2),(1,3),(2,3)],             [(0,1,2),(0,1,3)], [(0,2,3,1)]],
        "HIH":     [[(0,1),(1,2),(1,3),(2,3)],             [(0,1,2),(0,1,3)], [(0,2,3,1)]],
        "GLN":     [[(0,1)]],
        "ASN":     [[(0,1)]],
        "SER":     [[(0,1)]],
        "THR":     [[(0,1)]],
        "ARG":     [[(0,1),(1,2)],                         [(0,1,2)]],
        "LYS":     [[(0,1),(1,2)],                         [(0,1,2)]],
        "ASP":     [[(0,1)]],
        "GLU":     [[(0,1)]],
        "CYS":     [[(0,1)]],
        "ILE":     [[(0,1)]],
        "LEU":     [[(0,1)]],
        "MET":     [[(0,1)]],
        "PRO":     [[(0,1)]],
        "HYP":     [[(0,1)]],
        "VAL":     [[(0,1)]],
        "ALA":     [],
        "GLY":     [],
        }
        
        #----+----------------+
        ## C | SPECIAL BONDS  |
        #----+----------------+
        
        self.special = {
            # Used for sulfur bridges
            # ATOM 1         ATOM 2          BOND LENGTH   FORCE CONSTANT
            (("SC1","CYS"), ("SC1","CYS")):     (0.24,         None),
            }
        
        # By default use an elastic network
        self.ElasticNetwork = False 

        # Elastic networks bond shouldn't lead to exclusions (type 6) 
        # But Elnedyn has been parametrized with type 1.
        self.EBondType = 6
        
        #----+----------------+
        ## D | INTERNAL STUFF |
        #----+----------------+
        
        
        ## BACKBONE BEAD TYPE ##                                                                    
        # Dictionary of default bead types (*D)                                                     
        self.bbBeadDictD  = hash(bbss,self.bbdef)                                                             
        # Dictionary of dictionaries of types for specific residues (*S)                            
        self.bbBeadDictS  = dict([(i,hash(bbss,self.bbtyp[i])) for i in list(self.bbtyp.keys())])                        
        
        ## BB BOND TYPE ##                                                                          
        # Dictionary of default abond types (*D)                                                    
        self.bbBondDictD = hash(bbss,list(zip(self.bbldef,self.bbkb)))                                                   
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbBondDictS = dict([(i,hash(bbss,list(zip(self.bbltyp[i],self.bbkbtyp[i])))) for i in list(self.bbltyp.keys())])       
        # This is tricky to read, but it gives the right bondlength/force constant
        
        ## BBB ANGLE TYPE ##                                                                        
        # Dictionary of default angle types (*D)                                                    
        self.bbAngleDictD = hash(bbss,list(zip(self.bbadef,self.bbka)))                                                  
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbAngleDictS = dict([(i,hash(bbss,list(zip(self.bbatyp[i],self.bbkatyp[i])))) for i in list(self.bbatyp.keys())])      
                    
        ## BBBB DIHEDRAL TYPE ##                                                                    
        # Dictionary of default dihedral types (*D)                                                 
        self.bbDihedDictD = hash(bbss,list(zip(self.bbddef,self.bbkd,self.bbdmul)))                                           
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbDihedDictS = dict([(i,hash(bbss,list(zip(self.bbdtyp[i],self.bbkdtyp[i])))) for i in list(self.bbdtyp.keys())])      
        
    # The following function returns the backbone bead for a given residue and                   
    # secondary structure type.                                                                 
    # 1. Look up the proper dictionary for the residue                                          
    # 2. Get the proper type from it for the secondary structure                                
    # If the residue is not in the dictionary of specials, use the default                      
    # If the secondary structure is not listed (in the residue specific                         
    # dictionary) revert to the default.                                                        
    def bbGetBead(self,r1,ss="C"):                                                                   
        return self.bbBeadDictS.get(r1,self.bbBeadDictD).get(ss,self.bbBeadDictD.get(ss))                      
    
    
    def bbGetBond(self,r,a,ss):
        # Retrieve parameters for each residue from table defined above
        b1 = self.bbBondDictS.get(r[0],self.bbBondDictD).get(ss[0],self.bbBondDictD.get(ss[0]))
        b2 = self.bbBondDictS.get(r[1],self.bbBondDictD).get(ss[1],self.bbBondDictD.get(ss[1]))
        # Determine which parameters to use for the bond
        # Python3 no longer returns a NoneType if None is compared to an integer.
        if b1[1] == None or b2[1] == None:
            return ( (b1[0]+b2[0])/2, None )
        else:
            return ( (b1[0]+b2[0])/2, min(b1[1],b2[1]) )
    
    def bbGetAngle(self,r,ca,ss):
        # PRO in helices is dominant
        if r[1] == "PRO" and ss[1] in "H123":
            return self.bbAngleDictS["PRO"].get(ss[1])
        else:
            # Retrieve parameters for each residue from table defined above
            a = [ self.bbAngleDictS.get(r[0],self.bbAngleDictD).get(ss[0],self.bbAngleDictD.get(ss[0])),
                  self.bbAngleDictS.get(r[1],self.bbAngleDictD).get(ss[1],self.bbAngleDictD.get(ss[1])),
                  self.bbAngleDictS.get(r[2],self.bbAngleDictD).get(ss[2],self.bbAngleDictD.get(ss[2])) ]
            # Sort according to force constant
            a.sort(key=lambda i: (i[1],i[0]))
            # This selects the set with the smallest force constant and the smallest angle
            return a[0]
        
    def messages(self):
        '''Prints any force-field specific logging messages.'''
        import logging
        logging.info('Note: Cysteine bonds are 0.24 nm constraints, instead of the published 0.39nm/5000kJ/mol.')
################################
## 6 # FORCE FIELD PARAMETERS ##  -> @FF <-
################################

class martini21p(object):
    ff = True
    def __init__(self):

        # parameters are defined here for the following (protein) forcefields:
        self.name = 'martini21p'
        
        # Charged types:
        self.charges = {"Qd":1, "Qa":-1, "SQd":1, "SQa":-1, "RQd":1, "AQa":-1}                                                           #@#
        
        
        #----+---------------------+
        ## A | BACKBONE PARAMETERS |
        #----+---------------------+
        #
        # bbss  lists the one letter secondary structure code
        # bbdef lists the corresponding default backbone beads
        # bbtyp lists the corresponding residue specific backbone beads
        #
        # bbd   lists the structure specific backbone bond lengths
        # bbkb  lists the corresponding bond force constants
        #
        # bba   lists the structure specific angles
        # bbka  lists the corresponding angle force constants
        #
        # bbd   lists the structure specific dihedral angles
        # bbkd  lists the corresponding force constants
        #
        # -=NOTE=- 
        #  if the secondary structure types differ between bonded atoms
        #  the bond is assigned the lowest corresponding force constant 
        #
        # -=NOTE=-
        # if proline is anywhere in the helix, the BBB angle changes for 
        # all residues
        #
        
        ###############################################################################################
        ## BEADS ##                                                                         #                 
        #                          F     E     H     1     2     3     T     S     C        # SS one letter   
        self.bbdef    =    spl(" N0   Nda    N0    Nd    Na   Nda   Nda    P5    P5")  # Default beads   #@#
        self.bbtyp    = {                                                                   #                 #@#
                     "ALA": spl(" C5    N0    C5    N0    N0    N0    N0    P4    P4"),# ALA specific    #@#
                     "PRO": spl(" C5    N0    C5    N0    Na    N0    N0    Na    Na"),# PRO specific    #@#
                     "HYP": spl(" C5    N0    C5    N0    N0    N0    N0    Na    Na") # HYP specific    #@#
        }                                                                                   #                 #@#
        ## BONDS ##                                                                         #                 
        self.bbldef   =             (.365, .350, .350, .350, .350, .350, .350, .350, .350)  # BB bond lengths #@#
        self.bbkb     =             (1250, 1250, 1250, 1250, 1250, 1250,  500,  400,  400)  # BB bond kB      #@#
        self.bbltyp   = {}                                                                  #                 #@#
        self.bbkbtyp  = {}                                                                  #                 #@#
        ## ANGLES ##                                                                        #                 
        self.bbadef   =             ( 119.2,134,   96,   96,   96,   96,  100,  130,  127)  # BBB angles      #@#
        self.bbka     =             ( 150,   25,  700,  700,  700,  700,   25,   25,   25)  # BBB angle kB    #@#
        self.bbatyp   = {                                                                   #                 #@#
               "PRO":               ( 119.2,134,   98,   98,   98,   98,  100,  130,  127), # PRO specific    #@#
               "HYP":               ( 119.2,134,   98,   98,   98,   98,  100,  130,  127)  # PRO specific    #@#
        }                                                                                   #                 #@#
        self.bbkatyp  = {                                                                   #                 #@#
               "PRO":               ( 150,   25,  100,  100,  100,  100,   25,   25,   25), # PRO specific    #@#
               "HYP":               ( 150,   25,  100,  100,  100,  100,   25,   25,   25)  # PRO specific    #@#
        }                                                                                   #                 #@#
        ## DIHEDRALS ##                                                                     #                 
        self.bbddef   =             ( 90.7,   0, -120, -120, -120, -120)                    # BBBB dihedrals  #@#
        self.bbkd     =             ( 100,   10,  400,  400,  400,  400)                    # BBBB kB         #@#
        self.bbdmul   =             (   1,    1,    1,    1,    1,    1)                    # BBBB mltplcty   #@#
        self.bbdtyp   = {}                                                                  #                 #@#
        self.bbkdtyp  = {}                                                                  #                 #@#
                                                                                            #                 
        ###############################################################################################               
        
        # Some Forcefields use the Ca position to position the BB-bead (me like!)
        # martini 2.1 doesn't
        self.ca2bb = False 
        
        # BBS angle, equal for all ss types                                                         
        # Connects BB(i-1),BB(i),SC(i), except for first residue: BB(i+1),BB(i),SC(i)               
        #                 ANGLE   Ka                                                                
        self.bbsangle =      [   100,  25]                                                               #@#
        
        # Bonds for extended structures (more stable than using dihedrals)                          
        #               LENGTH FORCE                                                                
        self.ebonds   = {                                                                                #@#
               'short': [ .640, 2500],                                                              #@#
               'long' : [ .970, 2500]                                                               #@#
        }                                                                                           #@#
        
        
        #----+-----------------------+
        ## B | SIDE CHAIN PARAMETERS |
        #----+-----------------------+
        
        # To be compatible with Elnedyn, all parameters are explicitly defined, even if they are double.
        self.sidechains = {
            #RES#   BEADS                   BONDS                                                   ANGLES              DIHEDRALS
            #                               BB-SC          SC-SC                                        BB-SC-SC  SC-SC-SC
            "TRP": [spl("SC4 SP1 SC4 SC4"),[(0.300,5000)]+[(0.270,None) for i in range(5)],        [(210,50),(90,50),(90,50)], [(0,50),(0,200)]],
            "TYR": [spl("SC4 SC4 SP1"),    [(0.320,5000), (0.270,None), (0.270,None),(0.270,None)],[(150,50),(150,50)],        [(0,50)]],
            "PHE": [spl("SC4 SC4 SC4"),    [(0.310,7500), (0.270,None), (0.270,None),(0.270,None)],[(150,50),(150,50)],        [(0,50)]],
            "HIS": [spl("SC4 SP1 SP1"),    [(0.320,7500), (0.270,None), (0.270,None),(0.270,None)],[(150,50),(150,50)],        [(0,50)]],
            "HIH": [spl("SC4 SP1 SQd"),    [(0.320,7500), (0.270,None), (0.270,None),(0.270,None)],[(150,50),(150,50)],        [(0,50)]],
            "ARG": [spl("N0 Qd"),          [(0.330,5000), (0.340,5000)],                           [(180,25)]],
            "LYS": [spl("C3 Qd"),          [(0.330,5000), (0.280,5000)],                           [(180,25)]],
            "CYS": [spl("C5"),             [(0.310,7500)]],
            "ASP": [spl("Qa"),             [(0.320,7500)]],
            "GLU": [spl("Qa"),             [(0.400,5000)]],
            "ILE": [spl("C1"),            [(0.310,None)]],
            "LEU": [spl("C1"),            [(0.330,7500)]],
            "MET": [spl("C5"),             [(0.400,2500)]],
            "ASN": [spl("P5"),             [(0.320,5000)]],
            "PRO": [spl("C2"),            [(0.300,7500)]],
            "HYP": [spl("P1"),             [(0.300,7500)]],
            "GLN": [spl("P4"),             [(0.400,5000)]],
            "SER": [spl("P1"),             [(0.250,7500)]],
            "THR": [spl("P1"),             [(0.260,None)]],
            "VAL": [spl("C2"),            [(0.265,None)]],
            "ALA": [],
            "GLY": [],
            }
        
        # Not all (eg Elnedyn) forcefields use backbone-backbone-sidechain angles and BBBB-dihedrals.
        self.UseBBSAngles          = True 
        self.UseBBBBDihedrals      = True

        # Martini 2.2p has polar and charged residues with seperate charges.
        self.polar   = []
        self.charged = []

        # If masses or charged diverge from standard (45/72 and -/+1) they are defined here.
        self.mass_charge = {
        #RES   MASS               CHARGE
        }

        # Defines the connectivity between between beads
        self.connectivity = {
        #RES       BONDS                                   ANGLES             DIHEDRALS              V-SITE
        "TRP":     [[(0,1),(1,2),(1,3),(2,3),(2,4),(3,4)], [(0,1,2),(0,1,3)], [(0,2,3,1),(1,2,4,3)]],
        "TYR":     [[(0,1),(1,2),(1,3),(2,3)],             [(0,1,2),(0,1,3)], [(0,2,3,1)]],
        "PHE":     [[(0,1),(1,2),(1,3),(2,3)],             [(0,1,2),(0,1,3)], [(0,2,3,1)]],
        "HIS":     [[(0,1),(1,2),(1,3),(2,3)],             [(0,1,2),(0,1,3)], [(0,2,3,1)]],
        "HIH":     [[(0,1),(1,2),(1,3),(2,3)],             [(0,1,2),(0,1,3)], [(0,2,3,1)]],
        "GLN":     [[(0,1)]],
        "ASN":     [[(0,1)]],
        "SER":     [[(0,1)]],
        "THR":     [[(0,1)]],
        "ARG":     [[(0,1),(1,2)],                         [(0,1,2)]],
        "LYS":     [[(0,1),(1,2)],                         [(0,1,2)]],
        "ASP":     [[(0,1)]],
        "GLU":     [[(0,1)]],
        "CYS":     [[(0,1)]],
        "ILE":     [[(0,1)]],
        "LEU":     [[(0,1)]],
        "MET":     [[(0,1)]],
        "PRO":     [[(0,1)]],
        "HYP":     [[(0,1)]],
        "VAL":     [[(0,1)]],
        "ALA":     [],
        "GLY":     [],
        }

        #----+----------------+
        ## C | SPECIAL BONDS  |
        #----+----------------+
        
        self.special = {
            # Used for sulfur bridges
            # ATOM 1         ATOM 2          BOND LENGTH   FORCE CONSTANT
            (("SC1","CYS"), ("SC1","CYS")):     (0.24,         None),
            }
       
        # By default use an elastic network
        self.ElasticNetwork = False
 
        # Elastic networks bond shouldn't lead to exclusions (type 6) 
        # But Elnedyn has been parametrized with type 1.
        self.EBondType = 6
        
        #----+----------------+
        ## D | INTERNAL STUFF |
        #----+----------------+
        
        
        ## BACKBONE BEAD TYPE ##                                                                    
        # Dictionary of default bead types (*D)                                                     
        self.bbBeadDictD  = hash(bbss,self.bbdef)                                                             
        # Dictionary of dictionaries of types for specific residues (*S)                            
        self.bbBeadDictS  = dict([(i,hash(bbss,self.bbtyp[i])) for i in list(self.bbtyp.keys())])                        
        
        ## BB BOND TYPE ##                                                                          
        # Dictionary of default abond types (*D)                                                    
        self.bbBondDictD = hash(bbss,list(zip(self.bbldef,self.bbkb)))                                                   
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbBondDictS = dict([(i,hash(bbss,list(zip(self.bbltyp[i],self.bbkbtyp[i])))) for i in list(self.bbltyp.keys())])       
        # This is tricky to read, but it gives the right bondlength/force constant
        
        ## BBB ANGLE TYPE ##                                                                        
        # Dictionary of default angle types (*D)                                                    
        self.bbAngleDictD = hash(bbss,list(zip(self.bbadef,self.bbka)))                                                  
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbAngleDictS = dict([(i,hash(bbss,list(zip(self.bbatyp[i],self.bbkatyp[i])))) for i in list(self.bbatyp.keys())])      
                    
        ## BBBB DIHEDRAL TYPE ##                                                                    
        # Dictionary of default dihedral types (*D)                                                 
        self.bbDihedDictD = hash(bbss,list(zip(self.bbddef,self.bbkd,self.bbdmul)))                                           
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbDihedDictS = dict([(i,hash(bbss,list(zip(self.bbdtyp[i],self.bbkdtyp[i])))) for i in list(self.bbdtyp.keys())])      
        
    # The following function returns the backbone bead for a given residue and                   
    # secondary structure type.                                                                 
    # 1. Look up the proper dictionary for the residue                                          
    # 2. Get the proper type from it for the secondary structure                                
    # If the residue is not in the dictionary of specials, use the default                      
    # If the secondary structure is not listed (in the residue specific                         
    # dictionary) revert to the default.                                                        
    def bbGetBead(self,r1,ss="C"):                                                                   
        return self.bbBeadDictS.get(r1,self.bbBeadDictD).get(ss,self.bbBeadDictD.get(ss))                      
    
    
    def bbGetBond(self,r,a,ss):
        # Retrieve parameters for each residue from table defined above
        b1 = self.bbBondDictS.get(r[0],self.bbBondDictD).get(ss[0],self.bbBondDictD.get(ss[0]))
        b2 = self.bbBondDictS.get(r[1],self.bbBondDictD).get(ss[1],self.bbBondDictD.get(ss[1]))
        # Determine which parameters to use for the bond
        # Python3 no longer returns a NoneType if None is compared to an integer.
        if b1[1] == None or b2[1] == None:
            return ( (b1[0]+b2[0])/2, None )
        else:
            return ( (b1[0]+b2[0])/2, min(b1[1],b2[1]) )
    
    def bbGetAngle(self,r,ca,ss):
        # PRO in helices is dominant
        if r[1] == "PRO" and ss[1] in "H123":
            return self.bbAngleDictS["PRO"].get(ss[1])
        else:
            # Retrieve parameters for each residue from table defined above
            a = [ self.bbAngleDictS.get(r[0],self.bbAngleDictD).get(ss[0],self.bbAngleDictD.get(ss[0])),
                  self.bbAngleDictS.get(r[1],self.bbAngleDictD).get(ss[1],self.bbAngleDictD.get(ss[1])),
                  self.bbAngleDictS.get(r[2],self.bbAngleDictD).get(ss[2],self.bbAngleDictD.get(ss[2])) ]
            # Sort according to force constant
            a.sort(key=lambda i: (i[1],i[0]))
            # This selects the set with the smallest force constant and the smallest angle
            return a[0]
        
    def messages(self):
        '''Prints any force-field specific logging messages.'''
        logging.info('Note: Cysteine bonds are 0.24 nm constraints, instead of the published 0.39nm/5000kJ/mol.')

################################
## 6 # FORCE FIELD PARAMETERS ##  -> @FF <-
################################

# New martini 2.2p parameters.
# Changed: 
#   Unstructured Pro backbone bead
#   Proline side chains
#   Phe sidechain
#   Trp sidechain
#   Polar beads
#   Helix BB-bonds to constraint      
#   Helix BB-bond length

class martini22p(object):
    ff = True
    def __init__(self):

        # parameters are defined here for the following (protein) forcefields:
        self.name = 'martini22p'
        
        # Charged types:
        self.charges = {"Qd":1, "Qa":-1, "SQd":1, "SQa":-1, "RQd":1, "AQa":-1}    #@#
        
        #----+---------------------+
        ## A | BACKBONE PARAMETERS |
        #----+---------------------+
        #
        # bbss  lists the one letter secondary structure code
        # bbdef lists the corresponding default backbone beads
        # bbtyp lists the corresponding residue specific backbone beads
        #
        # bbd   lists the structure specific backbone bond lengths
        # bbkb  lists the corresponding bond force constants
        #
        # bba   lists the structure specific angles
        # bbka  lists the corresponding angle force constants
        #
        # bbd   lists the structure specific dihedral angles
        # bbkd  lists the corresponding force constants
        #
        # -=NOTE=- 
        #  if the secondary structure types differ between bonded atoms
        #  the bond is assigned the lowest corresponding force constant 
        #
        # -=NOTE=-
        # if proline is anywhere in the helix, the BBB angle changes for 
        # all residues
        #
        
        ###############################################################################################
        ## BEADS ##                                                                         #                 
        #                              F     E     H     1     2     3     T     S     C    # SS one letter   
        self.bbdef    =    spl(" N0   Nda    N0    Nd    Na   Nda   Nda    P5    P5")  # Default beads   #@#
        self.bbtyp    = {                                                                   #                 #@#
                    "ALA": spl(" C5    N0    C5    N0    N0    N0    N0    P4    P4"), # ALA specific    #@#
                    "PRO": spl(" C5    N0    C5    N0    Na    N0    N0    P4    P4"), # PRO specific    #@#
                    "HYP": spl(" C5    N0    C5    N0    N0    N0    N0    P4    P4")  # HYP specific    #@#
        }                                                                                   #                 #@#
        ## BONDS ##                                                                         #                 
        self.bbldef   =             (.365, .350, .310, .310, .310, .310, .350, .350, .350)  # BB bond lengths #@#
        self.bbkb     =             (1250, 1250, None, None, None, None, 1250, 1250, 1250)  # BB bond kB      #@#
        self.bbltyp   = {}                                                                  #                 #@#
        self.bbkbtyp  = {}                                                                  #                 #@#
        ## ANGLES ##                                                                        #                 
        self.bbadef   =             ( 119.2,134,   96,   96,   96,   96,  100,  130,  127)  # BBB angles      #@#
        self.bbka     =             ( 150,   25,  700,  700,  700,  700,   25,   25,   25)  # BBB angle kB    #@#
        self.bbatyp   = {                                                                   #                 #@#
               "PRO":               ( 119.2,134,   98,   98,   98,   98,  100,  130,  127), # PRO specific    #@#
               "HYP":               ( 119.2,134,   98,   98,   98,   98,  100,  130,  127)  # PRO specific    #@#
        }                                                                                   #                 #@#
        self.bbkatyp  = {                                                                   #                 #@#
               "PRO":               ( 150,   25,  100,  100,  100,  100,   25,   25,   25), # PRO specific    #@#
               "HYP":               ( 150,   25,  100,  100,  100,  100,   25,   25,   25)  # PRO specific    #@#
        }                                                                                   #                 #@#
        ## DIHEDRALS ##                                                                     #                 
        self.bbddef   =             ( 90.7,   0, -120, -120, -120, -120)                    # BBBB dihedrals  #@#
        self.bbkd     =             ( 100,   10,  400,  400,  400,  400)                    # BBBB kB         #@#
        self.bbdmul   =             (   1,    1,    1,    1,    1,    1)                    # BBBB mltplcty   #@#
        self.bbdtyp   = {}                                                                  #                 #@#
        self.bbkdtyp  = {}                                                                  #                 #@#
                                                                                            #                 
        ###############################################################################################               
        
        # Some Forcefields use the Ca position to position the BB-bead (me like!)
        # martini 2.1 doesn't
        self.ca2bb = False 
        
        # BBS angle, equal for all ss types                                                         
        # Connects BB(i-1),BB(i),SC(i), except for first residue: BB(i+1),BB(i),SC(i)               
        #                 ANGLE   Ka                                                                
        self.bbsangle =      [   100,  25]                                                               #@#
        
        # Bonds for extended structures (more stable than using dihedrals)                          
        #               LENGTH FORCE                                                                
        self.ebonds   = {                                                                                #@#
               'short': [ .640, 2500],                                                              #@#
               'long' : [ .970, 2500]                                                               #@#
        }                                                                                           #@#
        
        
        #----+-----------------------+
        ## B | SIDE CHAIN PARAMETERS |
        #----+-----------------------+
        
        # To be compatible with Elnedyn, all parameters are explicitly defined, even if they are double.
        self.sidechains = {
          #RES#   BEADS                       BONDS                                                                   ANGLES                      DIHEDRALS        V-SITES
          #                                   BB-SC          SC-SC                                                    BB-SC-SC  SC-SC-SC
          "TRP": [spl("SC4 SNd SC5 SC5"),[(0.300,5000)]+[(0.270,None) for i in range(5)],                    [(210,50),(90,50),(90,50)], [(0,50),(0,200)]],
          "TYR": [spl("SC4 SC4 SP1"),    [(0.320,5000), (0.270,None), (0.270,None),(0.270,None)],            [(150,50),(150,50)],        [(0,50)]],
          "PHE": [spl("SC5 SC5 SC5"),    [(0.310,7500), (0.270,None), (0.270,None),(0.270,None)],            [(150,50),(150,50)],        [(0,50)]],
          "HIS": [spl("SC4 SP1 SP1"),    [(0.320,7500), (0.270,None), (0.270,None),(0.270,None)],            [(150,50),(150,50)],        [(0,50)]],
          "HIH": [spl("SC4 SP1 SQd D"),  [(0.320,7500), (0.270,None), (0.270,None),(0.270,None),(0.11,None)],[(150,50),(150,50)],        [(0,50)]],
          "GLN": [spl("Nda D D"),        [(0.400,5000), (0.280,None)],                                       [],                         [],              [(0.5,)]],
          "ASN": [spl("Nda D D"),        [(0.320,5000), (0.280,None)],                                       [],                         [],              [(0.5,)]],
          "SER": [spl("N0 D D"),         [(0.250,7500), (0.280,None)],                                       [],                         [],              [(0.5,)]],
          "THR": [spl("N0 D D"),         [(0.260,9000), (0.280,None)],                                       [],                         [],              [(0.5,)]],
          "ARG": [spl("N0 Qd D"),        [(0.330,5000), (0.340,5000), (0.110,None)],                         [(180,25)]],
          "LYS": [spl("C3 Qd D"),        [(0.330,5000), (0.280,5000), (0.110,None)],                         [(180,25)]],
          "ASP": [spl("Qa D"),           [(0.320,7500), (0.110,None)]],
          "GLU": [spl("Qa D"),           [(0.400,5000), (0.110,None)]],
          "CYS": [spl("C5"),             [(0.310,7500)]],
          "ILE": [spl("C1"),             [(0.310,None)]],
          "LEU": [spl("C1"),             [(0.330,7500)]],
          "MET": [spl("C5"),             [(0.400,2500)]],
          "PRO": [spl("C3"),             [(0.300,7500)]],
          "HYP": [spl("P1"),             [(0.300,7500)]],
          "VAL": [spl("C2"),             [(0.265,None)]],
          "ALA": [],
          "GLY": [],
          }
        
        # Not all (eg Elnedyn) forcefields use backbone-backbone-sidechain angles and BBBB-dihedrals.
        self.UseBBSAngles          = True 
        self.UseBBBBDihedrals      = True

        # Martini 2.2p has polar and charged residues with separate charges.
        self.polar   = ["GLN","ASN","SER","THR"]
        self.charged = ["ARG","LYS","ASP","GLU","HIH"]

        # If masses or charged diverge from standard (45/72 and -/+1) they are defined here.
        self.mass_charge = {
        #RES   MASS               CHARGE
        "GLN":[[0,36,36],         [0,0.42,-0.42]], 
        "ASN":[[0,36,36],         [0,0.46,-0.46]], 
        "SER":[[0,36,36],         [0,0.40,-0.40]],
        "THR":[[0,36,36],         [0,0.36,-0.36]],
        "ARG":[[72,36,36],        [0,0,1]],
        "LYS":[[72,36,36],        [0,0,1]],
        "HIH":[[45,45,36,36],     [0,0,0,1]],
        "ASP":[[36,36],           [0,-1]],
        "GLU":[[36,36],           [0,-1]],
        }

        self.connectivity = {
        #RES       BONDS                                   ANGLES             DIHEDRALS              V-SITE
        "TRP":     [[(0,1),(1,2),(1,3),(2,3),(2,4),(3,4)], [(0,1,2),(0,1,3)], [(0,2,3,1),(1,2,4,3)]],  
        "TYR":     [[(0,1),(1,2),(1,3),(2,3)],             [(0,1,2),(0,1,3)], [(0,2,3,1)]], 
        "PHE":     [[(0,1),(1,2),(1,3),(2,3)],             [(0,1,2),(0,1,3)], [(0,2,3,1)]],
        "HIS":     [[(0,1),(1,2),(1,3),(2,3)],             [(0,1,2),(0,1,3)], [(0,2,3,1)]],
        "HIH":     [[(0,1),(1,2),(1,3),(2,3),(3,4)],       [(0,1,2),(0,1,3)], [(0,2,3,1)]],
        "GLN":     [[(0,1),(2,3)],                         [],                [],                    [(1,2,3)]],
        "ASN":     [[(0,1),(2,3)],                         [],                [],                    [(1,2,3)]],
        "SER":     [[(0,1),(2,3)],                         [],                [],                    [(1,2,3)]],
        "THR":     [[(0,1),(2,3)],                         [],                [],                    [(1,2,3)]],
        "ARG":     [[(0,1),(1,2),(2,3)],                   [(0,1,2)]],
        "LYS":     [[(0,1),(1,2),(2,3)],                   [(0,1,2)]],
        "ASP":     [[(0,1),(1,2)]],
        "GLU":     [[(0,1),(1,2)]],
        "CYS":     [[(0,1)]],
        "ILE":     [[(0,1)]],
        "LEU":     [[(0,1)]],
        "MET":     [[(0,1)]],
        "PRO":     [[(0,1)]],
        "HYP":     [[(0,1)]],
        "VAL":     [[(0,1)]],
        "ALA":     [],
        "GLY":     [],
        }
 
        #----+----------------+
        ## C | SPECIAL BONDS  |
        #----+----------------+
        
        self.special = {
            # Used for sulfur bridges
            # ATOM 1         ATOM 2          BOND LENGTH   FORCE CONSTANT
            (("SC1","CYS"), ("SC1","CYS")):     (0.24,         None),
            }
        
        # By default use an elastic network
        self.ElasticNetwork = False 

        # Elastic networks bond shouldn't lead to exclusions (type 6) 
        # But Elnedyn has been parametrized with type 1.
        self.EBondType = 6
        
        #----+----------------+
        ## D | INTERNAL STUFF |
        #----+----------------+
        
        
        ## BACKBONE BEAD TYPE ##                                                                    
        # Dictionary of default bead types (*D)                                                     
        self.bbBeadDictD  = hash(bbss,self.bbdef)                                                             
        # Dictionary of dictionaries of types for specific residues (*S)                            
        self.bbBeadDictS  = dict([(i,hash(bbss,self.bbtyp[i])) for i in list(self.bbtyp.keys())])                        
        
        ## BB BOND TYPE ##                                                                          
        # Dictionary of default abond types (*D)                                                    
        self.bbBondDictD = hash(bbss,list(zip(self.bbldef,self.bbkb)))                                                   
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbBondDictS = dict([(i,hash(bbss,list(zip(self.bbltyp[i],self.bbkbtyp[i])))) for i in list(self.bbltyp.keys())])       
        # This is tricky to read, but it gives the right bondlength/force constant
        
        ## BBB ANGLE TYPE ##                                                                        
        # Dictionary of default angle types (*D)                                                    
        self.bbAngleDictD = hash(bbss,list(zip(self.bbadef,self.bbka)))                                                  
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbAngleDictS = dict([(i,hash(bbss,list(zip(self.bbatyp[i],self.bbkatyp[i])))) for i in list(self.bbatyp.keys())])      
                    
        ## BBBB DIHEDRAL TYPE ##                                                                    
        # Dictionary of default dihedral types (*D)                                                 
        self.bbDihedDictD = hash(bbss,list(zip(self.bbddef,self.bbkd,self.bbdmul)))                                           
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbDihedDictS = dict([(i,hash(bbss,list(zip(self.bbdtyp[i],self.bbkdtyp[i])))) for i in list(self.bbdtyp.keys())])      
        
    # The following function returns the backbone bead for a given residue and                   
    # secondary structure type.                                                                 
    # 1. Look up the proper dictionary for the residue                                          
    # 2. Get the proper type from it for the secondary structure                                
    # If the residue is not in the dictionary of specials, use the default                      
    # If the secondary structure is not listed (in the residue specific                         
    # dictionary) revert to the default.                                                        
    def bbGetBead(self,r1,ss="C"):                                                                   
        return self.bbBeadDictS.get(r1,self.bbBeadDictD).get(ss,self.bbBeadDictD.get(ss))                      
    
    
    def bbGetBond(self,r,a,ss):
        # Retrieve parameters for each residue from table defined above
        b1 = self.bbBondDictS.get(r[0],self.bbBondDictD).get(ss[0],self.bbBondDictD.get(ss[0]))
        b2 = self.bbBondDictS.get(r[1],self.bbBondDictD).get(ss[1],self.bbBondDictD.get(ss[1]))
        # Determine which parameters to use for the bond
        # Python3 no longer returns a NoneType if None is compared to an integer.
        if b1[1] == None or b2[1] == None:
            return ( (b1[0]+b2[0])/2, None )
        else:
            return ( (b1[0]+b2[0])/2, min(b1[1],b2[1]) )

    def bbGetAngle(self,r,ca,ss):
        # PRO in helices is dominant
        if r[1] == "PRO" and ss[1] in "H123":
            return self.bbAngleDictS["PRO"].get(ss[1])
        else:
            # Retrieve parameters for each residue from table defined above
            a = [ self.bbAngleDictS.get(r[0],self.bbAngleDictD).get(ss[0],self.bbAngleDictD.get(ss[0])),
                  self.bbAngleDictS.get(r[1],self.bbAngleDictD).get(ss[1],self.bbAngleDictD.get(ss[1])),
                  self.bbAngleDictS.get(r[2],self.bbAngleDictD).get(ss[2],self.bbAngleDictD.get(ss[2])) ]
            # Sort according to force constant
            a.sort(key=lambda i: (i[1],i[0]))
            # This selects the set with the smallest force constant and the smallest angle
            return a[0]

    def messages(self):
        '''Prints any force-field specific logging messages.'''
        import logging
        logging.warning('Bead names of charges in sidechains differ between .top/.itp and .pdb.')
        logging.warning('Using names in topology, as Gromacs does, gives the correct result.')
        logging.info('Note: Cysteine bonds are 0.24 nm constraints, instead of the published 0.39nm/5000kJ/mol.')
################################
## 6 # FORCE FIELD PARAMETERS ##  -> @FF <-
################################

class elnedyn(object):
    ff = True
    def __init__(self):
        '''The forcefield has been implemented with some changes compared to the published parameters:
        - Backbone-Backbone bonds are constraints in stead of strong bonds.
        - Trp has an extra constrain added to the sidechain
        - The Backbone sidechain bonds with high force constants are replaced by constraints except for Trp.
        '''

        # parameters are defined here for the following (protein) forcefields:
        self.name = 'elnedyn'
        
        # Charged types:
        self.charges = {"Qd":1, "Qa":-1, "RQd":1, "AQa":-1}                                                           #@#
        
        
        #----+---------------------+
        ## A | BACKBONE PARAMETERS |
        #----+---------------------+
        #
        # bbss  lists the one letter secondary structure code
        # bbdef lists the corresponding default backbone beads
        # bbtyp lists the corresponding residue specific backbone beads
        #
        # bbd   lists the structure specific backbone bond lengths
        # bbkb  lists the corresponding bond force constants
        #
        # bba   lists the structure specific angles
        # bbka  lists the corresponding angle force constants
        #
        # bbd   lists the structure specific dihedral angles
        # bbkd  lists the corresponding force constants
        #
        # -=NOTE=- 
        #  if the secondary structure types differ between bonded atoms
        #  the bond is assigned the lowest corresponding force constant 
        #
        # -=NOTE=-
        # if proline is anywhere in the helix, the BBB angle changes for 
        # all residues
        #
        
        ###############################################################################################
        ## BEADS ##                                                                          #                 
        #                               F     E     H     1     2     3     T     S     C    # SS one letter   
        self.bbdef    =     spl(" N0   Nda    N0    Nd    Na   Nda   Nda    P5    P5")  # Default beads   #@#
        self.bbtyp    = {                                                                    #                 #@#
                     "ALA": spl(" C5    N0    C5    N0    N0    N0    N0    P4    P4"), # ALA specific    #@#
                     "PRO": spl(" C5    N0    C5    N0    Na    N0    N0    Na    Na"), # PRO specific    #@#
                     "HYP": spl(" C5    N0    C5    N0    N0    N0    N0    Na    Na")  # HYP specific    #@#
        }                                                                                    #                 #@#
        ## BONDS ##                                                                          #                 
        self.bbldef   =             (.365, .350, .350, .350, .350, .350, .350, .350, .350)   # BB bond lengths #@#
        self.bbkb     =             (1250, 1250, 1250, 1250, 1250, 1250,  500,  400,  400)   # BB bond kB      #@#
        self.bbltyp   = {}                                                                   #                 #@#
        self.bbkbtyp  = {}                                                                   #                 #@#
        ## ANGLES ##                                                                         #                 
        self.bbadef   =             (119.2, 134,   96,   96,   96,   96,  100,  130,  127)   # BBB angles      #@#
        self.bbka     =             ( 150,   25,  700,  700,  700,  700,   25,   25,   25)   # BBB angle kB    #@#
        self.bbatyp   = {                                                                    #                 #@#
                    "PRO":          ( 119.2,134,   98,   98,   98,   98,  100,  130,  127),  # PRO specific    #@#
                    "HYP":          ( 119.2,134,   98,   98,   98,   98,  100,  130,  127)   # PRO specific    #@#
        }                                                                                    #                 #@#
        self.bbkatyp  = {                                                                    #                 #@#
                    "PRO":          ( 150,   25,  100,  100,  100,  100,   25,   25,   25),  # PRO specific    #@#
                    "HYP":          ( 150,   25,  100,  100,  100,  100,   25,   25,   25)   # PRO specific    #@#
        }                                                                                    #                 #@#
        ## DIHEDRALS ##                                                                      #                 
        self.bbddef   =             (90.7,    0, -120, -120, -120, -120)                     # BBBB dihedrals  #@#
        self.bbkd     =             ( 100,   10,  400,  400,  400,  400)                     # BBBB kB         #@#
        self.bbdmul   =             (   1,    1,    1,    1,    1,    1)                     # BBBB mltplcty   #@#
        self.bbdtyp   = {}                                                                   #                 #@#
        self.bbkdtyp  = {}                                                                   #                 #@#
                                                                                             #                 
        ###############################################################################################               
        
        # Some Forcefields use the Ca position to position the BB-bead (me like!)
        self.ca2bb = True 
        
        # BBS angle, equal for all ss types                                                         
        # Connects BB(i-1),BB(i),SC(i), except for first residue: BB(i+1),BB(i),SC(i)               
        #                      ANGLE   Ka                                                                
        self.bbsangle =      [   100,  25]                                                          #@#
        
        # Bonds for extended structures (more stable than using dihedrals)                          
        #               LENGTH FORCE                                                                
        self.ebonds   = {                                                                           #@#
               'short': [ .640, 2500],                                                              #@#
               'long' : [ .970, 2500]                                                               #@#
        }                                                                                           #@#
        
        
        #----+-----------------------+
        ## B | SIDE CHAIN PARAMETERS |
        #----+-----------------------+
        
        # Sidechain parameters for Elnedyn. (read from cg-2.1.dat). 
        # For HIS the order of bonds is changed and a bond with fc=0 is added.
        # In the elnedyn2, TRP has an extra, cross-ring constraint
        self.sidechains = {
        #RES#   BEADS                      BONDS                                                                    ANGLES                          DIHEDRALS
        'TRP': [spl("SC4 SP1 SC4 SC4"), [(0.255,73000), (0.220,None), (0.250,None), (0.280,None), (0.255,None), (0.35454,None)], [(142,30), (143,20), (104,50)], [(180,200)]],
        'TYR': [spl("SC4 SC4 SP1"),     [(0.335, 6000), (0.335,6000), (0.240,None), (0.310,None), (0.310,None)], [(70,100), (130, 50)]],
        'PHE': [spl("SC4 SC4 SC4"),     [(0.340, 7500), (0.340,7500), (0.240,None), (0.240,None), (0.240,None)], [(70,100), (125,100)]],
        'HIS': [spl("SC4 SP1 SP1"),     [(0.195, None), (0.193,None), (0.295,None), (0.216,None)],               [(135,100),(115, 50)]],
        'ARG': [spl("N0 Qd"),           [(0.250,12500), (0.350,6200)],                                           [(150,15)]],
        'LYS': [spl("C3 Qd"),           [(0.250,12500), (0.300,9700)],                                           [(150,20)]],
        'CYS': [spl("C5"),              [(0.240, None)]],
        'ASP': [spl("Qa"),              [(0.255, None)]],
        'GLU': [spl("Qa"),              [(0.310, 2500)]],
        'ILE': [spl("C1"),              [(0.225,13250)]],
        'LEU': [spl("C1"),              [(0.265, None)]],
        'MET': [spl("C5"),              [(0.310, 2800)]],
        'ASN': [spl("P5"),              [(0.250, None)]],
        'PRO': [spl("C2"),              [(0.190, None)]],
        'GLN': [spl("P4"),              [(0.300, 2400)]],
        'SER': [spl("P1"),              [(0.195, None)]],
        'THR': [spl("P1"),              [(0.195, None)]],
        'VAL': [spl("C2"),              [(0.200, None)]],
        'GLY': [],
        'ALA': [],
        }
        
        # Not all (eg Elnedyn) forcefields use backbone-backbone-sidechain angles and BBBB-dihedrals.
        self.UseBBSAngles        = False 
        self.UseBBBBDihedrals    = False

        # Martini 2.2p has polar and charged residues with seperate charges.
        self.polar   = []
        self.charged = []

        # If masses or charged diverge from standard (45/72 and -/+1) they are defined here.
        self.mass_charge = {
        #RES   MASS               CHARGE
        }

        # Defines the connectivity between between beads
        # Connectivity records for Elnedyn (read from cg-2.1.dat). 
        # For HIS the order of bonds is changed and a bond with fc=0 is added.
        self.connectivity = {
        #RES       BONDS                                     ANGLES                            DIHEDRALS       V-SITE
        "TRP":     [[(0, 1), (1, 2), (2, 4), (4, 3), (3, 1), (1, 4)],[(0, 1, 2), (0, 1, 4), (0, 1, 3)],[(1, 2, 3, 4)]],
        "TYR":     [[(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)],[(0, 1, 2), (0, 1, 3)]],
        "PHE":     [[(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)],[(0, 1, 2), (0, 1, 3)]],
        "HIS":     [[(0, 1), (1, 2), (1, 3), (2, 3)],        [(0, 1, 2), (0, 1, 3)]],
        "GLN":     [[(0,1)]],
        "ASN":     [[(0,1)]],
        "SER":     [[(0,1)]],
        "THR":     [[(0,1)]],
        "ARG":     [[(0,1),(1,2)],                         [(0,1,2)]],
        "LYS":     [[(0,1),(1,2)],                         [(0,1,2)]],
        "ASP":     [[(0,1)]],
        "GLU":     [[(0,1)]],
        "CYS":     [[(0,1)]],
        "ILE":     [[(0,1)]],
        "LEU":     [[(0,1)]],
        "MET":     [[(0,1)]],
        "PRO":     [[(0,1)]],
        "HYP":     [[(0,1)]],
        "VAL":     [[(0,1)]],
        "ALA":     [],
        "GLY":     [],
        }
       
        #----+----------------+
        ## C | SPECIAL BONDS  |
        #----+----------------+
        
        self.special = {
            # Used for sulfur bridges
            # ATOM 1         ATOM 2          BOND LENGTH   FORCE CONSTANT
            (("SC1","CYS"), ("SC1","CYS")):     (0.24,         None),
            }
       
        # By default use an elastic network
        self.ElasticNetwork = True 

        # Elastic networks bond shouldn't lead to exclusions (type 6) 
        # But Elnedyn has been parametrized with type 1.
        self.EBondType = 1
        
        #----+----------------+
        ## D | INTERNAL STUFF |
        #----+----------------+
        
        
        ## BACKBONE BEAD TYPE ##                                                                    
        # Dictionary of default bead types (*D)                                                     
        self.bbBeadDictD  = hash(bbss,self.bbdef)                                                             
        # Dictionary of dictionaries of types for specific residues (*S)                            
        self.bbBeadDictS  = dict([(i,hash(bbss,self.bbtyp[i])) for i in list(self.bbtyp.keys())])                        
         
        ## BB BOND TYPE ##                                                                          
        # Dictionary of default abond types (*D)                                                    
        self.bbBondDictD = hash(bbss,list(zip(self.bbldef,self.bbkb)))                                                   
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbBondDictS = dict([(i,hash(bbss,list(zip(self.bbltyp[i],self.bbkbtyp[i])))) for i in list(self.bbltyp.keys())])       
        # This is tricky to read, but it gives the right bondlength/force constant

        ## BBB ANGLE TYPE ##                                                                        
        # Dictionary of default angle types (*D)                                                    
        self.bbAngleDictD = hash(bbss,list(zip(self.bbadef,self.bbka)))                                                  
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbAngleDictS = dict([(i,hash(bbss,list(zip(self.bbatyp[i],self.bbkatyp[i])))) for i in list(self.bbatyp.keys())])      
       
        ## BBBB DIHEDRAL TYPE ##                                                                    
        # Dictionary of default dihedral types (*D)                                                 
        self.bbDihedDictD = hash(bbss,list(zip(self.bbddef,self.bbkd,self.bbdmul)))                                           
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbDihedDictS = dict([(i,hash(bbss,list(zip(self.bbdtyp[i],self.bbkdtyp[i])))) for i in list(self.bbdtyp.keys())])      

    # The following function returns the backbone bead for a given residue and                   
    # secondary structure type.                                                                 
    # 1. Look up the proper dictionary for the residue                                          
    # 2. Get the proper type from it for the secondary structure                                
    # If the residue is not in the dictionary of specials, use the default                      
    # If the secondary structure is not listed (in the residue specific                         
    # dictionary) revert to the default.                                                        
    def bbGetBead(self,r1,ss="C"):                                                                   
        return self.bbBeadDictS.get(r1,self.bbBeadDictD).get(ss,self.bbBeadDictD.get(ss))                      
    
    # For Elnedyn we need something else to get the bond length (much simpler due to Ca position BB's)
    def bbGetBond(self,r,ca,ss):
        import math
        # The 150000 forceconstant gave an error message, turning to constraints would be better.
        return ( math.sqrt(distance2(ca[0],ca[1]))/10., 150000 )
    
    def bbGetAngle(self,r,ca,ss):
        import math
        # Elnedyn takes angles from structure, with fc=40
        return (math.acos(cos_angle([i-j for i,j in zip(ca[0],ca[1])],[i-j for i,j in zip(ca[2],ca[1])]))/d2r, 40)


    def messages(self):
        '''Prints any force-field specific logging messages.'''
        import logging
        logging.info('The Elnedyn forcefield has been implemented with some changes compared to the published parameters:')
        logging.info('- Backbone-Backbone bonds use high force constant bonds instead of constraints.')
        logging.info('- Trp has an extra constrain added to the sidechain.')
        logging.info('- The Backbone sidechain bonds with high force constants are replaced by constraints except for Trp and His.')
        logging.info('- Cysteine bonds are 0.24 nm constraints, instead of the published 0.39nm/5000kJ/mol.')
        logging.warning('Elnedyn topologies might not give numerical stable simulations with a 20fs timestep.')
        logging.warning('This can be solved by setting all S-type bead masses to 72amu.') 
        pass
################################
## 6 # FORCE FIELD PARAMETERS ##  -> @FF <-
################################

class elnedyn22(object):
    '''The forcefield has been implemented with some changes compared to the published parameters:
    - Trp has an extra constrain added to the sidechain
    - The Backbone-Sidechain bonds with high force constants are replaced by constraints except for Trp and His.
    '''
    ff = True
    def __init__(self):

        # parameters are defined here for the following (protein) forcefields:
        self.name = 'elnedyn22'
        
        # Charged types:
        self.charges = {"Qd":1, "Qa":-1, "SQd":1, "SQa":-1, "RQd":1, "AQa":-1}                                                           #@#
        
        
        #----+---------------------+
        ## A | BACKBONE PARAMETERS |
        #----+---------------------+
        #
        # bbss  lists the one letter secondary structure code
        # bbdef lists the corresponding default backbone beads
        # bbtyp lists the corresponding residue specific backbone beads
        #
        # bbd   lists the structure specific backbone bond lengths
        # bbkb  lists the corresponding bond force constants
        #
        # bba   lists the structure specific angles
        # bbka  lists the corresponding angle force constants
        #
        # bbd   lists the structure specific dihedral angles
        # bbkd  lists the corresponding force constants
        #
        # -=NOTE=- 
        #  if the secondary structure types differ between bonded atoms
        #  the bond is assigned the lowest corresponding force constant 
        #
        # -=NOTE=-
        # if proline is anywhere in the helix, the BBB angle changes for 
        # all residues
        #
        
        ###############################################################################################
        ## BEADS ##                                                                          #                 
        #                               F     E     H     1     2     3     T     S     C    # SS one letter   
        self.bbdef    =     spl(" N0   Nda    N0    Nd    Na   Nda   Nda    P5    P5")  # Default beads   #@#
        self.bbtyp    = {                                                                    #                 #@#
                     "ALA": spl(" C5    N0    C5    N0    N0    N0    N0    P4    P4"), # ALA specific    #@#
                     "PRO": spl(" C5    N0    C5    N0    Na    N0    N0    P4    P4"), # PRO specific    #@#
                     "HYP": spl(" C5    N0    C5    N0    N0    N0    N0    P4    P4")  # HYP specific    #@#
        }                                                                                    #                 #@#
        ## BONDS ##                                                                          #                 
        self.bbldef   =             (.365, .350, .350, .350, .350, .350, .350, .350, .350)   # BB bond lengths #@#
        self.bbkb     =             (1250, 1250, 1250, 1250, 1250, 1250,  500,  400,  400)   # BB bond kB      #@#
        self.bbltyp   = {}                                                                   #                 #@#
        self.bbkbtyp  = {}                                                                   #                 #@#
        ## ANGLES ##                                                                         #                 
        self.bbadef   =             (119.2, 134,   96,   96,   96,   96,  100,  130,  127)   # BBB angles      #@#
        self.bbka     =             ( 150,   25,  700,  700,  700,  700,   25,   25,   25)   # BBB angle kB    #@#
        self.bbatyp   = {                                                                    #                 #@#
                    "PRO":          ( 119.2,134,   98,   98,   98,   98,  100,  130,  127),  # PRO specific    #@#
                    "HYP":          ( 119.2,134,   98,   98,   98,   98,  100,  130,  127)   # PRO specific    #@#
        }                                                                                    #                 #@#
        self.bbkatyp  = {                                                                    #                 #@#
                    "PRO":          ( 150,   25,  100,  100,  100,  100,   25,   25,   25),  # PRO specific    #@#
                    "HYP":          ( 150,   25,  100,  100,  100,  100,   25,   25,   25)   # PRO specific    #@#
        }                                                                                    #                 #@#
        ## DIHEDRALS ##                                                                      #                 
        self.bbddef   =             (90.7,    0, -120, -120, -120, -120)                     # BBBB dihedrals  #@#
        self.bbkd     =             ( 100,   10,  400,  400,  400,  400)                     # BBBB kB         #@#
        self.bbdmul   =             (   1,    1,    1,    1,    1,    1)                     # BBBB mltplcty   #@#
        self.bbdtyp   = {}                                                                   #                 #@#
        self.bbkdtyp  = {}                                                                   #                 #@#
                                                                                             #                 
        ###############################################################################################               
        
        # Some Forcefields use the Ca position to position the BB-bead (me like!)
        self.ca2bb = True 
        
        # BBS angle, equal for all ss types                                                         
        # Connects BB(i-1),BB(i),SC(i), except for first residue: BB(i+1),BB(i),SC(i)               
        #                      ANGLE   Ka                                                                
        self.bbsangle =      [   100,  25]                                                          #@#
        
        # Bonds for extended structures (more stable than using dihedrals)                          
        #               LENGTH FORCE                                                                
        self.ebonds   = {                                                                           #@#
               'short': [ .640, 2500],                                                              #@#
               'long' : [ .970, 2500]                                                               #@#
        }                                                                                           #@#
        
        
        #----+-----------------------+
        ## B | SIDE CHAIN PARAMETERS |
        #----+-----------------------+
        
        # Sidechain parameters for Elnedyn. (read from cg-2.1.dat). 
        # For HIS the order of bonds is changed and a bond with fc=0 is added.
        # In the elnedyn2, TRP has an extra, cross-ring constraint
        self.sidechains = {
        #RES#   BEADS                      BONDS                                                                    ANGLES                          DIHEDRALS
        'TRP': [spl("SC4 SNd SC5 SC5"), [(0.255,73000), (0.220,None), (0.250,None), (0.280,None), (0.255,None), (0.35454,None)], [(142,30), (143,20), (104,50)], [(180,200)]],
        'TYR': [spl("SC4 SC4 SP1"),     [(0.335, 6000), (0.335,6000), (0.240,None), (0.310,None), (0.310,None)], [(70,100), (130, 50)]],
        'PHE': [spl("SC5 SC5 SC5"),     [(0.340, 7500), (0.340,7500), (0.240,None), (0.240,None), (0.240,None)], [(70,100), (125,100)]],
        'HIS': [spl("SC4 SP1 SP1"),     [(0.195, None), (0.193,None), (0.295,None), (0.216,None)],               [(135,100),(115, 50)]],
        'HIH': [spl("SC4 SP1 SP1"),     [(0.195, None), (0.193,None), (0.295,None), (0.216,None)],               [(135,100),(115, 50)]],
        'ARG': [spl("N0 Qd"),           [(0.250,12500), (0.350,6200)],                                           [(150,15)]],
        'LYS': [spl("C3 Qd"),           [(0.250,12500), (0.300,9700)],                                           [(150,20)]],
        'CYS': [spl("C5"),              [(0.240, None)]],
        'ASP': [spl("Qa"),              [(0.255, None)]],
        'GLU': [spl("Qa"),              [(0.310, 2500)]],
        'ILE': [spl("C1"),              [(0.225,13250)]],
        'LEU': [spl("C1"),              [(0.265, None)]],
        'MET': [spl("C5"),              [(0.310, 2800)]],
        'ASN': [spl("P5"),              [(0.250, None)]],
        'PRO': [spl("C3"),              [(0.190, None)]],
        'GLN': [spl("P4"),              [(0.300, 2400)]],
        'SER': [spl("P1"),              [(0.195, None)]],
        'THR': [spl("P1"),              [(0.195, None)]],
        'VAL': [spl("C2"),              [(0.200, None)]],
        'GLY': [],
        'ALA': [],
        }
        
        # Not all (eg Elnedyn) forcefields use backbone-backbone-sidechain angles and BBBB-dihedrals.
        self.UseBBSAngles        = False 
        self.UseBBBBDihedrals    = False

        # Martini 2.2p has polar and charged residues with seperate charges.
        self.polar   = []
        self.charged = []

        # If masses or charged diverge from standard (45/72 and -/+1) they are defined here.
        self.mass_charge = {
        #RES   MASS               CHARGE
        }

        # Defines the connectivity between between beads
        # Connectivity records for Elnedyn (read from cg-2.1.dat). 
        # For HIS the order of bonds is changed and a bond with fc=0 is added.
        self.connectivity = {
        #RES       BONDS                                             ANGLES                            DIHEDRALS       V-SITE
        "TRP":     [[(0, 1), (1, 2), (2, 4), (4, 3), (3, 1), (1, 4)],[(0, 1, 2), (0, 1, 4), (0, 1, 3)],[(1, 2, 3, 4)]],
        "TYR":     [[(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)],        [(0, 1, 2), (0, 1, 3)]],
        "PHE":     [[(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)],        [(0, 1, 2), (0, 1, 3)]],
        "HIS":     [[(0, 1), (1, 2), (1, 3), (2, 3)],        [(0, 1, 2), (0, 1, 3)]],
        "HIH":     [[(0, 1), (1, 2), (1, 3), (2, 3)],        [(0, 1, 2), (0, 1, 3)]],
        "GLN":     [[(0,1)]],
        "ASN":     [[(0,1)]],
        "SER":     [[(0,1)]],
        "THR":     [[(0,1)]],
        "ARG":     [[(0,1),(1,2)],                         [(0,1,2)]],
        "LYS":     [[(0,1),(1,2)],                         [(0,1,2)]],
        "ASP":     [[(0,1)]],
        "GLU":     [[(0,1)]],
        "CYS":     [[(0,1)]],
        "ILE":     [[(0,1)]],
        "LEU":     [[(0,1)]],
        "MET":     [[(0,1)]],
        "PRO":     [[(0,1)]],
        "HYP":     [[(0,1)]],
        "VAL":     [[(0,1)]],
        "ALA":     [],
        "GLY":     [],
        }
       
        #----+----------------+
        ## C | SPECIAL BONDS  |
        #----+----------------+
        
        self.special = {
            # Used for sulfur bridges
            # ATOM 1         ATOM 2          BOND LENGTH   FORCE CONSTANT
            (("SC1","CYS"), ("SC1","CYS")):     (0.24,         None),
            }
       
        # By default use an elastic network
        self.ElasticNetwork = True 

        # Elastic networks bond shouldn't lead to exclusions (type 6) 
        # But Elnedyn has been parametrized with type 1.
        self.EBondType = 1
        
        #----+----------------+
        ## D | INTERNAL STUFF |
        #----+----------------+
        
        
        ## BACKBONE BEAD TYPE ##                                                                    
        # Dictionary of default bead types (*D)                                                     
        self.bbBeadDictD  = hash(bbss,self.bbdef)                                                             
        # Dictionary of dictionaries of types for specific residues (*S)                            
        self.bbBeadDictS  = dict([(i,hash(bbss,self.bbtyp[i])) for i in list(self.bbtyp.keys())])                        
         
        ## BB BOND TYPE ##                                                                          
        # Dictionary of default abond types (*D)                                                    
        self.bbBondDictD = hash(bbss,list(zip(self.bbldef,self.bbkb)))                                                   
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbBondDictS = dict([(i,hash(bbss,list(zip(self.bbltyp[i],self.bbkbtyp[i])))) for i in list(self.bbltyp.keys())])       
        # This is tricky to read, but it gives the right bondlength/force constant

        ## BBB ANGLE TYPE ##                                                                        
        # Dictionary of default angle types (*D)                                                    
        self.bbAngleDictD = hash(bbss,list(zip(self.bbadef,self.bbka)))                                                  
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbAngleDictS = dict([(i,hash(bbss,list(zip(self.bbatyp[i],self.bbkatyp[i])))) for i in list(self.bbatyp.keys())])      
       
        ## BBBB DIHEDRAL TYPE ##                                                                    
        # Dictionary of default dihedral types (*D)                                                 
        self.bbDihedDictD = hash(bbss,list(zip(self.bbddef,self.bbkd,self.bbdmul)))                                           
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbDihedDictS = dict([(i,hash(bbss,list(zip(self.bbdtyp[i],self.bbkdtyp[i])))) for i in list(self.bbdtyp.keys())])      

    # The following function returns the backbone bead for a given residue and                   
    # secondary structure type.                                                                 
    # 1. Look up the proper dictionary for the residue                                          
    # 2. Get the proper type from it for the secondary structure                                
    # If the residue is not in the dictionary of specials, use the default                      
    # If the secondary structure is not listed (in the residue specific                         
    # dictionary) revert to the default.                                                        
    def bbGetBead(self,r1,ss="C"):                                                                   
        return self.bbBeadDictS.get(r1,self.bbBeadDictD).get(ss,self.bbBeadDictD.get(ss))                      
    
    # For Elnedyn we need something else to get the bond length (much simpler due to Ca position BB's)
    def bbGetBond(self,r,ca,ss):
        import math
        # The 150000 forceconstant sometimes gives an error message, turning to constraints could be better.
        return ( math.sqrt(distance2(ca[0],ca[1]))/10., 150000 )
    
    def bbGetAngle(self,r,ca,ss):
        import math
        # Elnedyn takes angles from structure, with fc=40
        return (math.acos(cos_angle([i-j for i,j in zip(ca[0],ca[1])],[i-j for i,j in zip(ca[2],ca[1])]))/d2r, 40)

    def messages(self):
        '''Prints any force-field specific logging messages.'''
        import logging
        logging.info('The elnedyn forcefield has been implemented with some changes compared to the published parameters:')
        #logging.info('- Backbone-Backbone bonds are constraints in stead of high force constant bonds.')
        logging.info('- Backbone-Backbone bonds use high force constant bonds instead of constraints.')
        logging.info('- Trp has an extra constrain added to the sidechain.')
        logging.info('- The Backbone sidechain bonds with high force constants are replaced by constraints except for Trp and His.')
        logging.info('- Cysteine bonds are 0.24 nm constraints, instead of the published 0.39nm/5000kJ/mol.')
        logging.warning('Elnedyn topologies might not give numerical stable simulations with a 20fs timestep.')
        logging.warning('This can be solved by setting all S-type bead masses to 72amu.')
        pass
################################
## 6 # FORCE FIELD PARAMETERS ##  -> @FF <-
################################

class elnedyn22p(object):
    ff = True
    def __init__(self):

        # parameters are defined here for the following (protein) forcefields:
        self.name = 'elnedyn22p'
        
        # Charged types:
        self.charges = {"Qd":1, "Qa":-1, "SQd":1, "SQa":-1, "RQd":1, "AQa":-1}                                                           #@#
        
        
        #----+---------------------+
        ## A | BACKBONE PARAMETERS |
        #----+---------------------+
        #
        # bbss  lists the one letter secondary structure code
        # bbdef lists the corresponding default backbone beads
        # bbtyp lists the corresponding residue specific backbone beads
        #
        # bbd   lists the structure specific backbone bond lengths
        # bbkb  lists the corresponding bond force constants
        #
        # bba   lists the structure specific angles
        # bbka  lists the corresponding angle force constants
        #
        # bbd   lists the structure specific dihedral angles
        # bbkd  lists the corresponding force constants
        #
        # -=NOTE=- 
        #  if the secondary structure types differ between bonded atoms
        #  the bond is assigned the lowest corresponding force constant 
        #
        # -=NOTE=-
        # if proline is anywhere in the helix, the BBB angle changes for 
        # all residues
        #
        
        ###############################################################################################
        ## BEADS ##                                                                          #                 
        #                               F     E     H     1     2     3     T     S     C    # SS one letter   
        self.bbdef    =     spl(" N0   Nda    N0    Nd    Na   Nda   Nda    P5    P5")  # Default beads   #@#
        self.bbtyp    = {                                                                    #                 #@#
                     "ALA": spl(" C5    N0    C5    N0    N0    N0    N0    P4    P4"), # ALA specific    #@#
                     "PRO": spl(" C5    N0    C5    N0    Na    N0    N0    P4    P4"), # PRO specific    #@#
                     "HYP": spl(" C5    N0    C5    N0    N0    N0    N0    P4    P4")  # HYP specific    #@#
        }                                                                                    #                 #@#
        ## BONDS ##                                                                          #                 
        self.bbldef   =             (.365, .350, .350, .350, .350, .350, .350, .350, .350)   # BB bond lengths #@#
        self.bbkb     =             (1250, 1250, 1250, 1250, 1250, 1250,  500,  400,  400)   # BB bond kB      #@#
        self.bbltyp   = {}                                                                   #                 #@#
        self.bbkbtyp  = {}                                                                   #                 #@#
        ## ANGLES ##                                                                         #                 
        self.bbadef   =             (119.2, 134,   96,   96,   96,   96,  100,  130,  127)   # BBB angles      #@#
        self.bbka     =             ( 150,   25,  700,  700,  700,  700,   25,   25,   25)   # BBB angle kB    #@#
        self.bbatyp   = {                                                                    #                 #@#
                    "PRO":          ( 119.2,134,   98,   98,   98,   98,  100,  130,  127),  # PRO specific    #@#
                    "HYP":          ( 119.2,134,   98,   98,   98,   98,  100,  130,  127)   # PRO specific    #@#
        }                                                                                    #                 #@#
        self.bbkatyp  = {                                                                    #                 #@#
                    "PRO":          ( 150,   25,  100,  100,  100,  100,   25,   25,   25),  # PRO specific    #@#
                    "HYP":          ( 150,   25,  100,  100,  100,  100,   25,   25,   25)   # PRO specific    #@#
        }                                                                                    #                 #@#
        ## DIHEDRALS ##                                                                      #                 
        self.bbddef   =             (90.7,    0, -120, -120, -120, -120)                     # BBBB dihedrals  #@#
        self.bbkd     =             ( 100,   10,  400,  400,  400,  400)                     # BBBB kB         #@#
        self.bbdmul   =             (   1,    1,    1,    1,    1,    1)                     # BBBB mltplcty   #@#
        self.bbdtyp   = {}                                                                   #                 #@#
        self.bbkdtyp  = {}                                                                   #                 #@#
                                                                                             #                 
        ###############################################################################################               
        
        # Some Forcefields use the Ca position to position the BB-bead (me like!)
        self.ca2bb = True 
        
        # BBS angle, equal for all ss types                                                         
        # Connects BB(i-1),BB(i),SC(i), except for first residue: BB(i+1),BB(i),SC(i)               
        #                      ANGLE   Ka                                                                
        self.bbsangle =      [   100,  25]                                                          #@#
        
        # Bonds for extended structures (more stable than using dihedrals)                          
        #               LENGTH FORCE                                                                
        self.ebonds   = {                                                                           #@#
               'short': [ .640, 2500],                                                              #@#
               'long' : [ .970, 2500]                                                               #@#
        }                                                                                           #@#
        
        
        #----+-----------------------+
        ## B | SIDE CHAIN PARAMETERS |
        #----+-----------------------+
        
        # Sidechain parameters for Elnedyn. (read from cg-2.1.dat). 
        # For HIS the order of bonds is changed and a bond with fc=0 is added.
        # In the elnedyn2, TRP has an extra, cross-ring constraint
        self.sidechains = {
        #RES#   BEADS                      BONDS                                                                                      ANGLES                          DIHEDRALS  V-SITES
        'TRP': [spl("SC4 SNd SC5 SC5"), [(0.255,73000), (0.220,None), (0.250,None), (0.280,None), (0.255,None), (0.35454,None)], [(142,30), (143,20), (104,50)], [(180,200)]],
        'TYR': [spl("SC4 SC4 SP1"),     [(0.335, 6000), (0.335,6000), (0.240,None), (0.310,None), (0.310,None)],                 [(70,100), (130, 50)]],
        'PHE': [spl("SC5 SC5 SC5"),     [(0.340, 7500), (0.340,7500), (0.240,None), (0.240,None), (0.240,None)],                 [(70,100), (125,100)]],
        'HIS': [spl("SC4 SP1 SP1"),     [(0.195, None), (0.193,None), (0.295,None), (0.216,None)],                               [(135,100),(115, 50)]],
        'HIH': [spl("SC4 SP1 SQd"),     [(0.195,94000), (0.193,None), (0.295,None), (0.216,None), (0.11,None)],                  [(135,100),(115, 50)]],
        'GLN': [spl("Nda D D"),         [(0.300, 2400), (0.280,None)],                                                           [],                             [],         [(0.5,)]],
        'ASN': [spl("Nda D D"),         [(0.250,61000), (0.280,None)],                                                           [],                             [],         [(0.5,)]],
        'SER': [spl("N0 D D"),          [(0.195,94000), (0.280,None)],                                                           [],                             [],         [(0.5,)]],
        'THR': [spl("N0 D D"),          [(0.195,94000), (0.280,None)],                                                           [],                             [],         [(0.5,)]],
        'ARG': [spl("N0 Qd D"),         [(0.250,12500), (0.350,6200), (0.110,None)],                                             [(150,15)]],
        'LYS': [spl("C3 Qd D"),         [(0.250,12500), (0.300,9700), (0.110,None)],                                             [(150,20)]],
        'ASP': [spl("Qa D"),            [(0.255, None), (0.110,None)]],
        'GLU': [spl("Qa D"),            [(0.310, 2500), (0.110,None)]],
        'CYS': [spl("C5"),              [(0.240, None)]],
        'ILE': [spl("C1"),              [(0.225,13250)]],
        'LEU': [spl("C1"),              [(0.265, None)]],
        'MET': [spl("C5"),              [(0.310, 2800)]],
        'PRO': [spl("C3"),              [(0.190, None)]],
        'HYP': [spl("P1"),              [(0.190, None)]],
        'VAL': [spl("C2"),              [(0.200, None)]],
        'GLY': [],
        'ALA': [],
        }
        
        # Not all (eg Elnedyn) forcefields use backbone-backbone-sidechain angles and BBBB-dihedrals.
        self.UseBBSAngles        = False 
        self.UseBBBBDihedrals    = False

        # Martini 2.2p has polar and charged residues with seperate charges.
        self.polar   = ["GLN","ASN","SER","THR"]
        self.charged = ["ARG","LYS","ASP","GLU","HIH"]

        # If masses or charged diverge from standard (45/72 and -/+1) they are defined here.
        self.mass_charge = {
        #RES   MASS               CHARGE
        "GLN":[[0,36,36],         [0,0.42,-0.42]],
        "ASN":[[0,36,36],         [0,0.46,-0.46]],
        "SER":[[0,36,36],         [0,0.40,-0.40]],
        "THR":[[0,36,36],         [0,0.36,-0.36]],
        "HIH":[[72,72,36,36],     [0,0,0,1]],
        "ARG":[[72,36,36],        [0,0,1]],
        "LYS":[[72,36,36],        [0,0,1]],
        "ASP":[[36,36],           [0,-1]],
        "GLU":[[36,36],           [0,-1]],
        }

        # Defines the connectivity between between beads
        # The polar sidechains have charged dummy beads, connected with a constraint
        # The charged sidechains have a charged dummy bead.
        self.connectivity = {
        #RES       BONDS                                              ANGLES                            DIHEDRALS       V-SITE
        "TRP":     [[(0, 1), (1, 2), (2, 4), (4, 3), (3, 1), (1, 4)],[(0, 1, 2), (0, 1, 4), (0, 1, 3)],[(1, 2, 3, 4)]],
        "TYR":     [[(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)],        [(0, 1, 2), (0, 1, 3)]],
        "PHE":     [[(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)],        [(0, 1, 2), (0, 1, 3)]],
        "HIS":     [[(0, 1), (1, 2), (1, 3), (2, 3)],                [(0, 1, 2), (0, 1, 3)]],
        "HIH":     [[(0, 1), (1, 2), (1, 3), (2, 3), (3, 4)],        [(0, 1, 2), (0, 1, 3)],           [(0, 2, 3, 1)]],
        "GLN":     [[(0, 1), (2, 3)],                                [],                               [],              [(1,2,3)]],
        "ASN":     [[(0, 1), (2, 3)],                                [],                               [],              [(1,2,3)]],
        "SER":     [[(0, 1), (2, 3)],                                [],                               [],              [(1,2,3)]],
        "THR":     [[(0, 1), (2, 3)],                                [],                               [],              [(1,2,3)]],
        "ARG":     [[(0, 1), (1, 2), (2, 3)],                        [(0,1,2)]],
        "LYS":     [[(0, 1), (1, 2), (2, 3)],                        [(0,1,2)]],
        "ASP":     [[(0, 1), (1, 2)]],
        "GLU":     [[(0, 1), (1, 2)]],
        "CYS":     [[(0, 1)]],
        "ILE":     [[(0, 1)]],
        "LEU":     [[(0, 1)]],
        "MET":     [[(0, 1)]],
        "PRO":     [[(0, 1)]],
        "HYP":     [[(0, 1)]],
        "VAL":     [[(0, 1)]],
        "ALA":     [],
        "GLY":     [],
        }
       
        #----+----------------+
        ## C | SPECIAL BONDS  |
        #----+----------------+
        
        self.special = {
            # Used for sulfur bridges
            # ATOM 1         ATOM 2          BOND LENGTH   FORCE CONSTANT
            (("SC1","CYS"), ("SC1","CYS")):     (0.24,         None),
            }
       
        # By default use an elastic network
        self.ElasticNetwork = True 

        # Elastic networks bond shouldn't lead to exclusions (type 6) 
        # But Elnedyn has been parametrized with type 1.
        self.EBondType = 1
        
        #----+----------------+
        ## D | INTERNAL STUFF |
        #----+----------------+
        
        
        ## BACKBONE BEAD TYPE ##                                                                    
        # Dictionary of default bead types (*D)                                                     
        self.bbBeadDictD  = hash(bbss,self.bbdef)                                                             
        # Dictionary of dictionaries of types for specific residues (*S)                            
        self.bbBeadDictS  = dict([(i,hash(bbss,self.bbtyp[i])) for i in list(self.bbtyp.keys())])                        
         
        ## BB BOND TYPE ##                                                                          
        # Dictionary of default abond types (*D)                                                    
        self.bbBondDictD = hash(bbss,list(zip(self.bbldef,self.bbkb)))                                                   
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbBondDictS = dict([(i,hash(bbss,list(zip(self.bbltyp[i],self.bbkbtyp[i])))) for i in list(self.bbltyp.keys())])       
        # This is tricky to read, but it gives the right bondlength/force constant

        ## BBB ANGLE TYPE ##                                                                        
        # Dictionary of default angle types (*D)                                                    
        self.bbAngleDictD = hash(bbss,list(zip(self.bbadef,self.bbka)))                                                  
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbAngleDictS = dict([(i,hash(bbss,list(zip(self.bbatyp[i],self.bbkatyp[i])))) for i in list(self.bbatyp.keys())])      
       
        ## BBBB DIHEDRAL TYPE ##                                                                    
        # Dictionary of default dihedral types (*D)                                                 
        self.bbDihedDictD = hash(bbss,list(zip(self.bbddef,self.bbkd,self.bbdmul)))                                           
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbDihedDictS = dict([(i,hash(bbss,list(zip(self.bbdtyp[i],self.bbkdtyp[i])))) for i in list(self.bbdtyp.keys())])      

    # The following function returns the backbone bead for a given residue and                   
    # secondary structure type.                                                                 
    # 1. Look up the proper dictionary for the residue                                          
    # 2. Get the proper type from it for the secondary structure                                
    # If the residue is not in the dictionary of specials, use the default                      
    # If the secondary structure is not listed (in the residue specific                         
    # dictionary) revert to the default.                                                        
    def bbGetBead(self,r1,ss="C"):                                                                   
        return self.bbBeadDictS.get(r1,self.bbBeadDictD).get(ss,self.bbBeadDictD.get(ss))                      
    
    # For Elnedyn we need something else to get the bond length (much simpler due to Ca position BB's)
    def bbGetBond(self,r,ca,ss):
        import math
        # The 150000 forceconstant gave an error message, turning to constraints would be better.
        return ( math.sqrt(distance2(ca[0],ca[1]))/10., 150000 )
    
    def bbGetAngle(self,r,ca,ss):
        import math
        # Elnedyn takes angles from structure, with fc=40
        return (math.acos(cos_angle([i-j for i,j in zip(ca[0],ca[1])],[i-j for i,j in zip(ca[2],ca[1])]))/d2r, 40)

    def messages(self):
        '''Prints any force-field specific logging messages.'''
        import logging
        logging.info('The elnedyn forcefield has been implemented with some changes compared to the published parameters:')
        #logging.info('- Backbone-Backbone bonds are constraints in stead of high force constant bonds.')
        logging.info('- Backbone-Backbone bonds use high force constant bonds instead of constraints.')
        logging.info('- Trp has an extra constrain added to the sidechain.')
        logging.info('- The Backbone sidechain bonds with high force constants are replaced by constraints except for Trp and His and the polar sidechains.')
        logging.info('- Cysteine bonds are 0.24 nm constraints, instead of the published 0.39nm/5000kJ/mol.')
        logging.warning('Elnedyn topologies might not give numerical stable simulations with a 20fs timestep.')
        logging.warning('This can be solved by setting all S-type bead masses to 72amu.')
        pass

###################################
## 1 # OPTIONS AND DOCUMENTATION ##  -> @DOC <-
###################################

class martini3(object):
    
    # Derived from martini22
    
    ff = True
    def __init__(self):

        # parameters are defined here for the following (protein) forcefields:
        self.name = 'martini3031_beta'
        
        # Charged types:
        self.charges = { 
             "Qd":1,  "Qa":-1,  # normal beads
            "SQd":1, "SQa":-1,  # small beads
            "TQd":1, "TQa":-1,  # tiny beads
            "RQd":1, "AQa":-1,  # are these included in v3?
            }

        
        #----+---------------------+
        ## A | BACKBONE PARAMETERS |
        #----+---------------------+
        #
        # bbss  lists the one letter secondary structure code
        # bbdef lists the corresponding default backbone beads
        # bbtyp lists the corresponding residue specific backbone beads
        #
        # bbd   lists the structure specific backbone bond lengths
        # bbkb  lists the corresponding bond force constants
        #
        # bba   lists the structure specific angles
        # bbka  lists the corresponding angle force constants
        #
        # bbd   lists the structure specific dihedral angles
        # bbkd  lists the corresponding force constants
        #
        # -=NOTE=- 
        #  if the secondary structure types differ between bonded atoms
        #  the bond is assigned the lowest corresponding force constant 
        #
        # -=NOTE=-
        # if proline is anywhere in the helix, the BBB angle changes for 
        # all residues
        #
        
        ###############################################################################################
        ## BEADS ##                                                                         #                 
        #                              F     E     H     1     2     3     T     S     C    # SS one letter   
        self.bbdef    =    spl(" N0   Nda    N0    Nd    Na   Nda   Nda    P5    P5")  # Default beads   #@#
        self.bbtyp    = {                                                                   #                 #@#
                    "ALA": spl(" C5    N0    C5    N0    N0    N0    N0    P4    P4"), # ALA specific    #@#
                    "PRO": spl(" C5    N0    C5    N0    Na    N0    N0    P4    P4"), # PRO specific    #@#
                    "HYP": spl(" C5    N0    C5    N0    N0    N0    N0    P4    P4")  # HYP specific    #@#
        }                                                                                   #                 #@#
        ## BONDS ##                                                                         #                 
        self.bbldef   =             (.365, .350, .310, .310, .310, .310, .350, .350, .350)  # BB bond lengths #@#
        self.bbkb     =             (1250, 1250, None, None, None, None, 1250, 1250, 1250)  # BB bond kB      #@#
        self.bbltyp   = {}                                                                  #                 #@#
        self.bbkbtyp  = {}                                                                  #                 #@#
        ## ANGLES ##                                                                        #                 
        self.bbadef   =             ( 119.2,134,   96,   96,   96,   96,  100,  130,  127)  # BBB angles      #@#
        self.bbka     =             ( 150,   25,  700,  700,  700,  700,   20,   20,   20)  # BBB angle kB    #@#
        self.bbatyp   = {                                                                   #                 #@#
               "PRO":               ( 119.2,134,   98,   98,   98,   98,  100,  130,  127), # PRO specific    #@#
               "HYP":               ( 119.2,134,   98,   98,   98,   98,  100,  130,  127)  # PRO specific    #@#
        }                                                                                   #                 #@#
        self.bbkatyp  = {                                                                   #                 #@#
               "PRO":               ( 150,   25,  100,  100,  100,  100,   25,   25,   25), # PRO specific    #@#
               "HYP":               ( 150,   25,  100,  100,  100,  100,   25,   25,   25)  # PRO specific    #@#
        }                                                                                   #                 #@#
        ## DIHEDRALS ##                                                                     #                 
        self.bbddef   =             ( 90.7,   0, -120, -120, -120, -120)                    # BBBB dihedrals  #@#
        self.bbkd     =             ( 100,   10,  400,  400,  400,  400)                    # BBBB kB         #@#
        self.bbdmul   =             (   1,    1,    1,    1,    1,    1)                    # BBBB mltplcty   #@#
        self.bbdtyp   = {}                                                                  #                 #@#
        self.bbkdtyp  = {}                                                                  #                 #@#
                                                                                            #                 
        ###############################################################################################               
        
        # Some Forcefields use the Ca position to position the BB-bead (me like!)
        # martini 2.1 doesn't
        self.ca2bb = False 
        
        # BBS angle, equal for all ss types                                                         
        # Connects BB(i-1),BB(i),SC(i), except for first residue: BB(i+1),BB(i),SC(i)               
        #                 ANGLE   Ka                                                                
        self.bbsangle =      [   100,  25]                                                               #@#
        
        # Bonds for extended structures (more stable than using dihedrals)                          
        #               LENGTH FORCE                                                                
        self.ebonds   = {                                                                                #@#
               'short': [ .640, 2500],                                                              #@#
               'long' : [ .970, 2500]                                                               #@#
        }                                                                                           #@#
        
        
        #----+-----------------------+
        ## B | SIDE CHAIN PARAMETERS |
        #----+-----------------------+
        
        ## This section needs to be updated with data from Paulo
        
        # To be compatible with Elnedyn, all parameters are explicitly defined, even if they are double.
        self.sidechains = {
            #RES#   BEADS                   BONDS                                                   ANGLES              DIHEDRALS
            #                               BB-SC          SC-SC                                        BB-SC-SC  SC-SC-SC
            "TRP": [spl("SC4 SNd SC5 SC5"),[(0.300,5000)]+[(0.270,None) for i in range(5)],        [(210,50),(90,50),(90,50)], [(0,50),(0,200)]],
            "TYR": [spl("SC4 SC4 SP1"),    [(0.320,5000), (0.270,None), (0.270,None),(0.270,None)],[(150,50),(150,50)],        [(0,50)]],
            "PHE": [spl("SC5 SC5 SC5"),    [(0.310,7500), (0.270,None), (0.270,None),(0.270,None)],[(150,50),(150,50)],        [(0,50)]],
            "HIS": [spl("SC4 SP1 SP1"),    [(0.320,7500), (0.270,None), (0.270,None),(0.270,None)],[(150,50),(150,50)],        [(0,50)]],
            "HIH": [spl("SC4 SP1 SQd"),    [(0.320,7500), (0.270,None), (0.270,None),(0.270,None)],[(150,50),(150,50)],        [(0,50)]],
            "ARG": [spl("N0 Qd"),          [(0.330,5000), (0.340,5000)],                           [(180,25)]],
            "LYS": [spl("C3 Qd"),          [(0.330,5000), (0.280,5000)],                           [(180,25)]],
            "CYS": [spl("C5"),             [(0.310,7500)]],
            "ASP": [spl("Qa"),             [(0.320,7500)]],
            "GLU": [spl("Qa"),             [(0.400,5000)]],
            "ILE": [spl("AC1"),            [(0.310,None)]],
            "LEU": [spl("AC1"),            [(0.330,7500)]],
            "MET": [spl("C5"),             [(0.400,2500)]],
            "ASN": [spl("P5"),             [(0.320,5000)]],
            "PRO": [spl("C3"),             [(0.300,7500)]],
            "HYP": [spl("P1"),             [(0.300,7500)]],
            "GLN": [spl("P4"),             [(0.400,5000)]],
            "SER": [spl("P1"),             [(0.250,7500)]],
            "THR": [spl("P1"),             [(0.260,None)]],
            "VAL": [spl("AC2"),            [(0.265,None)]],
            "ALA": [],
            "GLY": [],
            }
        
        # Not all (eg Elnedyn) forcefields use backbone-backbone-sidechain angles and BBBB-dihedrals.
        self.UseBBSAngles          = True 
        self.UseBBBBDihedrals      = True

        # Martini 2.2p has polar and charged residues with seperate charges.
        self.polar   = []
        self.charged = []

        # If masses or charged diverge from standard (45/72 and -/+1) they are defined here.
        self.mass_charge = {
        #RES   MASS               CHARGE
        }

        # Defines the connectivity between between beads
        self.connectivity = {
        #RES       BONDS                                   ANGLES             DIHEDRALS              V-SITE
        "TRP":     [[(0,1),(1,2),(1,3),(2,3),(2,4),(3,4)], [(0,1,2),(0,1,3)], [(0,2,3,1),(1,2,4,3)]],  
        "TYR":     [[(0,1),(1,2),(1,3),(2,3)],             [(0,1,2),(0,1,3)], [(0,2,3,1)]], 
        "PHE":     [[(0,1),(1,2),(1,3),(2,3)],             [(0,1,2),(0,1,3)], [(0,2,3,1)]],
        "HIS":     [[(0,1),(1,2),(1,3),(2,3)],             [(0,1,2),(0,1,3)], [(0,2,3,1)]],
        "HIH":     [[(0,1),(1,2),(1,3),(2,3)],             [(0,1,2),(0,1,3)], [(0,2,3,1)]],
        "GLN":     [[(0,1)]],
        "ASN":     [[(0,1)]],
        "SER":     [[(0,1)]],
        "THR":     [[(0,1)]],
        "ARG":     [[(0,1),(1,2)],                         [(0,1,2)]],
        "LYS":     [[(0,1),(1,2)],                         [(0,1,2)]],
        "ASP":     [[(0,1)]],
        "GLU":     [[(0,1)]],
        "CYS":     [[(0,1)]],
        "ILE":     [[(0,1)]],
        "LEU":     [[(0,1)]],
        "MET":     [[(0,1)]],
        "PRO":     [[(0,1)]],
        "HYP":     [[(0,1)]],
        "VAL":     [[(0,1)]],
        "ALA":     [],
        "GLY":     [],
        }
        
        #----+----------------+
        ## C | SPECIAL BONDS  |
        #----+----------------+
        
        self.special = {
            # Used for sulfur bridges
            # ATOM 1         ATOM 2          BOND LENGTH   FORCE CONSTANT
            (("SC1","CYS"), ("SC1","CYS")):     (0.24,         None),
            }
        
        # By default use an elastic network
        self.ElasticNetwork = False 

        # Elastic networks bond shouldn't lead to exclusions (type 6) 
        # But Elnedyn has been parametrized with type 1.
        self.EBondType = 6
        
        #----+----------------+
        ## D | INTERNAL STUFF |
        #----+----------------+
        
        
        ## BACKBONE BEAD TYPE ##                                                                    
        # Dictionary of default bead types (*D)                                                     
        self.bbBeadDictD  = hash(bbss,self.bbdef)                                                             
        # Dictionary of dictionaries of types for specific residues (*S)                            
        self.bbBeadDictS  = dict([(i,hash(bbss,self.bbtyp[i])) for i in list(self.bbtyp.keys())])                        
        
        ## BB BOND TYPE ##                                                                          
        # Dictionary of default abond types (*D)                                                    
        self.bbBondDictD = hash(bbss,list(zip(self.bbldef,self.bbkb)))                                                   
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbBondDictS = dict([(i,hash(bbss,list(zip(self.bbltyp[i],self.bbkbtyp[i])))) for i in list(self.bbltyp.keys())])       
        # This is tricky to read, but it gives the right bondlength/force constant
        
        ## BBB ANGLE TYPE ##                                                                        
        # Dictionary of default angle types (*D)                                                    
        self.bbAngleDictD = hash(bbss,list(zip(self.bbadef,self.bbka)))                                                  
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbAngleDictS = dict([(i,hash(bbss,list(zip(self.bbatyp[i],self.bbkatyp[i])))) for i in list(self.bbatyp.keys())])      
                    
        ## BBBB DIHEDRAL TYPE ##                                                                    
        # Dictionary of default dihedral types (*D)                                                 
        self.bbDihedDictD = hash(bbss,list(zip(self.bbddef,self.bbkd,self.bbdmul)))                                           
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbDihedDictS = dict([(i,hash(bbss,list(zip(self.bbdtyp[i],self.bbkdtyp[i])))) for i in list(self.bbdtyp.keys())]) 


        #----+--------+
        ## E | SUGARS |
        #----+--------+
        
        # Non-reducing end sugars. 
        
        self.sugarinternals = {
            #RES    BEADS                          BONDS                                                     ANGLES                      DIHEDRALS
            # terminal beta-D-GlcpNAc
            #             C12 C34 C56 2Ac            C12-C34       C34-C56       C12-C56       C12-2Ac       C34-C12-2Ac C56-C12-2Ac    
            "0YB":  [spl("TNd SP4 SP1 TNa"),     [(0.276, None),(0.307, None),(0.334,18000),(0.300,15000)],[  ( 90,200) , (120,200)],],
                 
            # 4-beta-D-GlcpNAc
            #             C12 C34 C56 2Ac            C12-C34       C34-C56       C12-C56       C12-2Ac       C34-C12-2Ac C56-C12-2Ac    
            "4YB":  [spl("TNd SP1 SP1 TNa"),     [(0.276, None),(0.307, None),(0.334,18000),(0.300,15000)],[  ( 90,200) , (120,200)],],
                 
            # 3,6-beta-D-Manp     
            #             C12 C34 C56                C12-C34       C34-C56       C12-C56                         
            "VMB":  [spl("TP1 SP1 SN0"    ),     [(0.276, None),(0.318,22000),(0.350,22000)],],
                 
            # 3,4,6-beta-D-Manp     
            #             C12 C34 C56                C12-C34       C34-C56       C12-C56                         
            "QMB":  [spl("TP1 SN0 SN0"    ),     [(0.276, None),(0.318,22000),(0.350,22000)],],
                 
            # 2-alpha-D-Manp     
            #             C12 C34 C56                C12-C34       C34-C56       C12-C56                        
            "2MA":  [spl("TN0 SP4 SP1"    ),     [(0.277, None),(0.322,16000),(0.350,13250)],],
     
            # 3-alpha-D-Manp     
            #             C12 C34 C56                C12-C34       C34-C56       C12-C56                        
            "3MA":  [spl("TP1 SP1 SP1"    ),     [(0.277, None),(0.322,16000),(0.350,13250)],],
     
            # terminal alpha-D-Manp     
            #             C12 C34 C56                C12-C34       C34-C56       C12-C56                        
            "0MA":  [spl("TP1 SP4 SP1"    ),     [(0.278, None),(0.295, None),(0.308, None)],],
     
            # 3,4-beta-D-GlcpNAc     
            #             C12 C34 C56 2Ac            C12-C34       C34-C56       C12-C56       C12-2Ac       C34-C12-2Ac C34-C12-2Ac    
            "WYB":  [spl("TNd SN0 SP1 TNa"),     [(0.276, None),(0.307, None),(0.334,18000),(0.300,15000)],[  ( 90,200) , (120,200)],],
     
            # terminal alpha-L-Fucp                 
            #             C12 C34 C56                C12-C34       C34-C56       C12-C56                        
            "0fA":  [spl("TP1 SP4 TN0"    ),     [(0.269, None),(0.285, None),(0.348, None)],],
     
            # 3-beta-D-Galp     
            #             C12 C34 C56                C12-C34       C34-C56       C12-C56                         
            "3LB":  [spl("TP1 SP1 SP1"    ),     [(0.276, None),(0.318,22000),(0.350,22000)],],

            # 6-beta-D-Galp     
            #             C12 C34 C56                C12-C34       C34-C56       C12-C56                         
            "6LB":  [spl("TP1 SP4 SN0"    ),     [(0.276, None),(0.318,22000),(0.350,22000)],],

            # terminal alpha-D-Neup5Ac
            #             C12 345  C67 NAc C89       C12-345       345-C67       C12-C67       345-NAc       C67-C89         (0,1,3),   (2,1,3)    (0,2,4),  ( 1,2,4)       
            "0SA":  [spl("SQa SNda SP1 SP2  P4"),[(0.336, None),(0.315, None),(0.382,15000),(0.357,15000),(0.294,10000)],[  (122,700), ( 70,100), (109,200), (157,350)],],

            # terminal alpha-D-Glcp     
            #             C12 C34 C56                C12-C34       C34-C56       C12-C56                        
            "0GA":  [spl("TP1 SP4 SP1"    ),     [(0.269, None),(0.375, None),(0.331, None)],],

            # terminal beta-D-Glcp     
            #             C12 C34 C56                C12-C34       C34-C56       C12-C56                        
            "0GB":  [spl("TP1 SP4 SP1"    ),     [(0.276, None),(0.323, None),(0.384, None)],],

            }
            
        
        self.redends = {
            
            # terminal alpha-D-Glcp     
            #             C12 C34 C56                C12-C34       C34-C56       C12-C56                        
            "0GA":  [spl("SP4 SP4 SP1"    ),     [(0.322, None),(0.375, None),(0.331, None)],],

            # terminal beta-D-Glcp     
            #             C12 C34 C56                C12-C34       C34-C56       C12-C56                        
            "0GB":  [spl("SP4 SP4 SP1"    ),     [(0.331, None),(0.323, None),(0.384, None)],],

            # 4-beta-D-GlcpNAc
            #             C12 C34 C56 2Ac            C12-C34       C34-C56       C12-C56       C12-2Ac       C34-C12-2Ac C56-C12-2Ac    
            "4YB":  [spl("SP1 SP1 SP1 TNa"),     [(0.331, None),(0.307, None),(0.334,18000),(0.300,15000)],[  ( 90,200) , (120,200)],],

            }

        
        
        
        self.sugarlinkpar = {
            #                                       BONDS                  ANGLES                                                  DIHEDRALS
            #                                         (0,1)             [[0],[1,0]], [[1,0],[1]] ,[[2,0],[1]]
            "DManpa1-OG1SER":                   [[(0.310, None)],      [(142,300),    ( 97,350),     ( 68,400)],                    []],
            #                                         (0,1)             [[0],[1,0]], [[1,0],[1]] ,[[2,0],[1]]
            "DManpa1-OG1THR":                   [[(0.355, None)],      [(108,150),    ( 95,200),     ( 67,200)],                    []],
            #                                         (0,1)             [[0],[1,0]], [[1,0],[1]] ,[[2,0],[1]]
            "DGlcpNAcb1-ND2ASN":                [[(0.377, 1800)],      [(150,100),   ( 98,170),     ( 63,150)],                    []],
            #                                         (0,1)             [[0],[1,0]], [[0],[1,2]], [[1,0],[1]] ,[[2,0],[1]]
            "DGlcpNAcb1-4DGlcpNAcb":            [[(0.433, 5000)],      [(154,400),   ( 72,250),     ( 86,600),   ( 50,300)],       []],
            #                                         (0,1)             [[0],[1,0]], [[0],[1,2]], [[1,0],[1]] ,[[2,0],[1]]
            "DManpb1-4DGlcpNAcb":               [[(0.338,15000)],      [(168,500),   ( 90, 50),     (131,700),   ( 78,200)],       []],
            #                                         (0,1)             [[0],[1,0]], [[0],[1,2]], [[1,0],[1]] ,[[2,0],[1]]
            "DManpa1-3DManpb":                  [[(0.372, 9000)],      [(105,200),   (169,600),     (110,300),   ( 86,200)],       []],
            #                                         (0,2)             [[0],[2,0]], [[0],[2,1]], [[1,0],[2]] ,[[2,0],[2]]
            "DManpa1-6DManpb":                  [[(0.354, 8000)],      [(90, 150),   (117,150),     (118, 400),  ( 86,100)],       []],
            #                                         (0,0)             [[0],[0,1]], [[0],[0,2]], [[1,0],[0]] ,[[2,0],[0]]
            "DGlcpNAcb1-2DManpa":               [[(0.370, 3250)],      [( 75, 80),   ( 90, 80),     ( 97,550),   ( 61,400)],       []],
            #                                         (0,1)             [[0],[1,0]], [[0],[1,2]], [[1,0],[1]] ,[[2,0],[1]]
            "DGalpb1-4DGlcpNAcb":               [[(0.365, None)],      [(167,700),   ( 78,250),     (112,500),   ( 71,350)],       []],
            #                                         (0,1)             [[0],[1,0]], [[0],[1,2]], [[1,0],[1]] ,[[2,0],[1]]
            "DNeup5Aca2-3DGalpb":               [[(0.335, 2000)],      [(108,180),   (137, 50),     (100,250),   ( 80,220)],       []],
            #                                         (0,1)             [[0],[1,0]], [[0],[1,2]], [[1,0],[1]] ,[[2,0],[1]]
            "LFucpa1-3DGlcpNAcb":               [[(0.365, None)],      [(167, 700),  ( 78, 250),    (112,500),   ( 71,350)],       []],
            }                                         
        # Defines the connectivity of sugars. Independent of anomeric state, D/L-configuration, but ring size. First letter
        # of the key encodes one-letter code (eg. "Y" = GlcNAc), second letter ring size ("p" = pyranose, "f"=furanose) 
        
        self.sugarinternalcon = {
        
            #RES       BONDS                                   ANGLES                             DIHEDRALS              V-SITES              EXCLUSIONS
            "Yp":     [[(0,1),(1,2),(0,2),(0,3)],              [(1,0,3),(2,0,3)],                 [],                    [],                  [(1,3),(2,3)]],  
            "Gp":     [[(0,1),(1,2),(0,2)],                    [],                                [],                    [],                  []],  
            "Mp":     [[(0,1),(1,2),(0,2)],                    [],                                [],                    [],                  []],  
            "fp":     [[(0,1),(1,2),(0,2)],                    [],                                [],                    [],                  []],  
            "Lp":     [[(0,1),(1,2),(0,2)],                    [],                                [],                    [],                  []],  
            "Sp":     [[(0,1),(1,2),(0,2),(1,3),(2,4)],        [(0,1,3),(2,1,3),(0,2,4),(1,2,4)], [],                    [],                  [(0,3),(0,4),(1,4),(2,3)]],  
            }

        self.sugarlinkcon = {
            #   Sugar A   Res B           BONDS                 ANGLES                                              DIHEDRALS
            #                             A,B                    A    B        A    B       A    B      A     B   etc...
            "DGlcpNAcb1-ND2ASN":       [[[0,1]],              [[[0],[1,0]],              [[1,0],[1]], [[2,0],[1]]], [],                    [],                  [[0,0],      [1,1],[2,1]]],
            "DManpa1-OG1SER":          [[[0,1]],              [[[0],[1,0]],              [[1,0],[1]], [[2,0],[1]]], [],                    [],                  [[0,0],      [1,1],[2,1]]],
            "DManpa1-OG1THR":          [[[0,1]],              [[[0],[1,0]],              [[1,0],[1]], [[2,0],[1]]], [],                    [],                  [[0,0],      [1,1],[2,1]]],
            "DGlcpNAcb1-4DGlcpNAcb":   [[[0,1]],              [[[0],[1,0]], [[0],[1,2]], [[1,0],[1]] ,[[2,0],[1]]], [],                    [],                  [[0,0],[0,2],[1,1],[2,1]]],
            "DManpb1-4DGlcpNAcb":      [[[0,1]],              [[[0],[1,0]], [[0],[1,2]], [[1,0],[1]] ,[[2,0],[1]]], [],                    [],                  [[0,0],[0,2],[1,1],[2,1]]],
            "DManpa1-3DManpb":         [[[0,1]],              [[[0],[1,0]], [[0],[1,2]], [[1,0],[1]] ,[[2,0],[1]]], [],                    [],                  [[0,0],[0,2],[1,1],[2,1]]],
            "DManpa1-6DManpb":         [[[0,2]],              [[[0],[2,0]], [[0],[2,1]], [[1,0],[2]] ,[[2,0],[2]]], [],                    [],                  [[0,0],[0,1],[1,2],[2,2]]],
            "DGlcpNAcb1-2DManpa":      [[[0,0]],              [[[0],[0,1]], [[0],[0,2]], [[1,0],[0]] ,[[2,0],[0]]], [],                    [],                  [[0,1],[0,2],[1,0],[2,0]]],
            "DGalpb1-4DGlcpNAcb":      [[[0,1]],              [[[0],[1,0]], [[0],[1,2]], [[1,0],[1]] ,[[2,0],[1]]], [],                    [],                  [[0,0],[0,2],[1,1],[2,1]]],
            "DNeup5Aca2-3DGalpb":      [[[0,1]],              [[[0],[1,0]], [[0],[1,2]], [[1,0],[1]] ,[[2,0],[1]]], [],                    [],                  [[0,0],[0,2],[1,1],[2,1]]],
            "LFucpa1-3DGlcpNAcb":      [[[0,1]],              [[[0],[1,0]], [[0],[1,2]], [[1,0],[1]] ,[[2,0],[1]]], [],                    [],                  [[0,0],[0,2],[1,1],[2,1]]],

            }


        
           
    # The following function returns the backbone bead for a given residue and                   
    # secondary structure type.                                                                 
    # 1. Look up the proper dictionary for the residue                                          
    # 2. Get the proper type from it for the secondary structure                                
    # If the residue is not in the dictionary of specials, use the default                      
    # If the secondary structure is not listed (in the residue specific                         
    # dictionary) revert to the default.                                                        
    def bbGetBead(self,r1,ss="C"):                                                                   
        return self.bbBeadDictS.get(r1,self.bbBeadDictD).get(ss,self.bbBeadDictD.get(ss))                      
    
    
    def bbGetBond(self,r,a,ss):
        # Retrieve parameters for each residue from table defined above
        b1 = self.bbBondDictS.get(r[0],self.bbBondDictD).get(ss[0],self.bbBondDictD.get(ss[0]))
        b2 = self.bbBondDictS.get(r[1],self.bbBondDictD).get(ss[1],self.bbBondDictD.get(ss[1]))
        # Determine which parameters to use for the bond
        # Python3 no longer returns a NoneType if None is compared to an integer.
        if b1[1] == None or b2[1] == None:
            return ( (b1[0]+b2[0])/2, None )
        else:
            return ( (b1[0]+b2[0])/2, min(b1[1],b2[1]) )
    
    def bbGetAngle(self,r,ca,ss):
        # PRO in helices is dominant
        if r[1] == "PRO" and ss[1] in "H123":
            return self.bbAngleDictS["PRO"].get(ss[1])
        else:
            # Retrieve parameters for each residue from table defined above
            a = [ self.bbAngleDictS.get(r[0],self.bbAngleDictD).get(ss[0],self.bbAngleDictD.get(ss[0])),
                  self.bbAngleDictS.get(r[1],self.bbAngleDictD).get(ss[1],self.bbAngleDictD.get(ss[1])),
                  self.bbAngleDictS.get(r[2],self.bbAngleDictD).get(ss[2],self.bbAngleDictD.get(ss[2])) ]
            # Sort according to force constant
            a.sort(key=lambda i: (i[1],i[0]))
            # This selects the set with the smallest force constant and the smallest angle
            return a[0]
        
    def messages(self):
        '''Prints any force-field specific logging messages.'''
        import logging
        bl, fc = self.special[(("SC1","CYS"), ("SC1","CYS"))] 
        logging.info('Note: Disulfide bonds have length {} nm and force constant {}, instead of the published (0.39 nm, 5000 kJ/mol) bond.'.format(bl, fc))


import types, os
    
# This is a simple and versatile option class that allows easy
# definition and parsing of options.
class Option(object):
    def __init__(self, func=str, num=1, default=None, description=""):
        self.func        = func
        self.num         = num
        self.value       = default
        self.description = description

    def __bool__(self): 
        if self.func == bool:
            return self.value != False
        return bool(self.value)

    def __str__(self):
        return self.value and str(self.value) or ""

    def setvalue(self, v):
        if len(v) == 1:
            self.value = self.func(v[0])
        else:
            self.value = [self.func(i) for i in v]
    
# Parameters can be defined for multiple forcefields
# We look for them within the script...
forcefields = [str(ff).split('.')[-1] for ff in list(globals().values()) if (type(ff) == type and hasattr(ff,"ff"))]
# ... in the local directory, ....
forcefields += [ff[:-3] for ff in os.listdir(".") if ff[-6:] == "_ff.py"]
# ... and in the GMXDATA dir.
if "GMXDATA" in os.environ:
    forcefields += [ff[:-3] for ff in os.listdir(os.environ["GMXDATA"]+"/top/") if ff[-6:] == "_ff.py"]

# Lists for gathering arguments to options that can be specified
# multiple times on the command line.
lists = {
    'cystines': [],
    'merges'  : [],
    'links'   : [],
    'multi'   : [],
    }

# List of Help text and options.
# This way we can simply print this list if the user wants help.
options = [
#   NOTE: Options marked with (+) can be given multiple times on the command line
#   option              type number default description
    """
MARTINIZE.py is a script to create Coarse Grain Martini input files of
proteins, ready for use in the molecular dynamics simulations package
Gromacs. For more information on the Martini forcefield, see:
www.cgmartini.nl
and read our papers:
Monticelli et al., J. Chem. Theory Comput., 2008, 4(5), 819-834
de Jong et al., J. Chem. Theory Comput., 2013, DOI:10.1021/ct300646g

Primary input/output
--------------------
The input file (-f) should be a coordinate file in PDB or GROMOS
format. The format is inferred from the structure of the file. The
input can also be provided through stdin, allowing piping of
structures. The input structure can have multiple frames/models. If an output
structure file (-x) is given, each frame will be coarse grained,
resulting in a multimodel output structure. Having multiple frames may
also affect the topology. If secondary structure is determined
internally, the structure will be averaged over the frames. Likewise,
interatomic distances, as used for backbone bond lengths in Elnedyn
and in elastic networks, are also averaged over the frames available.

If an output file (-o) is indicated for the topology, that file will
be used for the master topology, using #include statements to link the
moleculetype definitions, which are written to separate files. If no
output filename is given, the topology and the moleculetype
definitions are written to stdout.

Secondary structure
-------------------
The secondary structure plays a central role in the assignment of atom
types and bonded interactions in MARTINI. Martinize allows
specification of the secondary structure as a string (-ss), or as a
file containing a specification in GROMACS' ssdump format
(-ss). Alternatively, DSSP can be used for an on-the-fly assignment of
the secondary structure. For this, the option -dssp has to be used
giving the location of the executable as the argument.
The option -collagen will set the whole structure to collagen. If this
is not what you want (eg only part of the structure is collagen, you
can give a secondary structure file/string (-ss) and specifiy collagen
as "F". Parameters for collagen are taken from: Gautieri et al.,
J. Chem. Theory Comput., 2010, 6, 1210-1218.
With multimodel input files, the secondary structure as determined with
DSSP will be averaged over the frames. In this case, a cutoff
can be specified (-ssc) indicating the fraction of frames to match a
certain secondary structure type for designation.

Topology
--------
Several options are available to tune the resulting topology. By
default, termini are charged, and chain breaks are kept neutral. This
behaviour can be changed using -nt and -cb, respectively.

Disulphide bridges can be specified using -cys. This option can be
given multiple times on the command line. The argument is a pair of
cysteine residues, using the format
chain/resn/resi,chain/resn/resi.
It is also possible to let martinize detect cysteine pairs based on a
cut-off distance of 0.22nm, by giving the keyword 'auto' as argument to -cys.
Alternatively, a different cut-off distance can be specified, which
will also trigger a search of pairs satisfying the distance
criterion (eg: -cys 0.32).

In addition to cystine bridges, links between other atoms can be
specified using -link. This requires specification of the atoms, using
the format
chain/resi/resn/atom,chain/resi/resn/atom,bondlength,forceconstant.
If only two atoms are given, a constraint will be added with length
equal to the (average) distance in the coordinate file. If a bond
length is added, but no force constant, then the bondlength will be
used to set a constraint.

Linking atoms requires that the atoms are part of the same
moleculetype. Therefore any link between chains will cause the chains
to be merged. Merges can also be specified explicitly, using the
option -merge with a comma-separated list of chain identifiers to be
joined into one moleculetype. The option -merge can be used several
times. Note that specifying a chain in several merge groups will cause
all chains involved to be merged into a single moleculetype.

The moleculetype definitions are written to topology include (.itp)
files, using a name consisting of the molecule class (e.g. Protein)
and the chain identifier. With -name a name can be specified instead.
By default, martinize only writes a moleculetype for each unique
molecule, inferred from the sequence and the secondary structure
definition. It is possible to force writing a moleculetype definition
for every single molecule, using -sep.

The option -p can be used to write position restraints, using the
force constant specified with -pf, which is set to 1000 kJ/mol
by default.

For stability, elastic bonds are used to retain the structure of
extended strands. The option -ed causes dihedrals to be used
instead.

Different forcefields can be specified with -ff. All the parameters and
options belonging to that forcefield  will be set (eg. bonded interactions,
BB-bead positions, Elastic Network, etc.). By default martini 2.1 is
used.

Elastic network
---------------
Martinize can write an elastic network for atom pairs within a cutoff
distance. The force constant (-ef) and the upper distance bound (-eu)
can be speficied. If a force field with an intrinsic Elastic
network is specified (eg. Elnedyn) with -ff, -elastic in implied and
the default values for the force constant and upper cutoff are used.
However, these can be overwritten.

Multiscaling
------------
Martinize can process a structure to yield a multiscale system,
consisting of a coordinate file with atomistic parts and
corresponding, overlaid coarsegrained parts. For chains that are
multiscaled, rather than writing a full moleculetype definition,
additional [atoms] and [virtual_sitesn] sections are written, to
be appended to the atomistic moleculetype definitions.
The option -multi can be specified multiple times, and takes a chain
identifier as argument. Alternatively, the keyword 'all' can be given
as argument, causing all chains to be multiscaled.
========================================================================\n
""",
    ("-f",        Option(str,                      1,     None, "Input file (PDB|GRO)")),
    ("-o",        Option(str,                      1,     None, "Output topology (TOP)")),
    ("-x",        Option(str,                      1,     None, "Output coarse grained structure (PDB)")),
    ("-n",        Option(str,                      1,     None, "Output index file with CG (and multiscale) beads.")),
    ("-nmap",     Option(str,                      1,     None, "Output index file containing per bead mapping.")),
    ("-v",        Option(bool,                     0,    False, "Verbose. Be load and noisy.")),
    ("-h",        Option(bool,                     0,    False, "Display this help.")),
    ("-ss",       Option(str,                      1,     None, "Secondary structure (File or string)")),
    ("-ssc",      Option(float,                    1,      0.5, "Cutoff fraction for ss in case of ambiguity (default: 0.5).")),
    ("-dssp",     Option(str,                      1,     None, "DSSP executable for determining structure")),
#    ("-pymol",    Option(str,                      1,     None, "PyMOL executable for determining structure")),
    ("-collagen", Option(bool,                     0,    False, "Use collagen parameters")),
    ("-his",      Option(bool,                     0,    False, "Interactively set the charge of each His-residue.")),
    ("-nt",       Option(bool,                     0,    False, "Set neutral termini (charged is default)")),
    ("-cb",       Option(bool,                     0,    False, "Set charges at chain breaks (neutral is default)")),
    ("-cys",      Option(lists['cystines'].append, 1,     None, "Disulphide bond (+)")),
    ("-link",     Option(lists['links'].append,    1,     None, "Link (+)")),
    ("-merge",    Option(lists['merges'].append,   1,     None, "Merge chains: e.g. -merge A,B,C (+)")),
#    ("-mixed",    Option(bool,                     0,    False, "Allow chains of mixed type (default: False)")),
    ("-name",     Option(str,                      1,     None, "Moleculetype name")),
    ("-p",        Option(str,                      1,   'None', "Output position restraints (None/All/Backbone) (default: None)")),
    ("-pf",       Option(float,                    1,     1000, "Position restraints force constant (default: 1000 kJ/mol/nm^2)")),
    ("-ed",       Option(bool,                     0,    False, "Use dihedrals for extended regions rather than elastic bonds)")),
    ("-sep",      Option(bool,                     0,    False, "Write separate topologies for identical chains.")),
    ("-ff",       Option(str,                      1, 'martini3', "Which forcefield to use: "+' ,'.join(n for n in forcefields))),
# Fij = Fc exp( -a (rij - lo)**p )
    ("-elastic",  Option(bool,                     0,    False, "Write elastic bonds")),
    ("-ef",       Option(float,                    1,      500, "Elastic bond force constant Fc")),
    ("-el",       Option(float,                    1,        0, "Elastic bond lower cutoff: F = Fc if rij < lo")),
    ("-eu",       Option(float,                    1,     0.90, "Elastic bond upper cutoff: F = 0  if rij > up")),
    ("-ea",       Option(float,                    1,        0, "Elastic bond decay factor a")),
    ("-ep",       Option(float,                    1,        1, "Elastic bond decay power p")),
    ("-em",       Option(float,                    1,        0, "Remove elastic bonds with force constant lower than this")),
    ("-eb",       Option(str,                      1,     'BB', "Comma separated list of bead names for elastic bonds")),
#    ("-hetatm",   Option(bool,                     0,    False, "Include HETATM records from PDB file (Use with care!)")),
    ("-multi",    Option(lists['multi'].append,    1,     None, "Chain to be set up for multiscaling (+)")),
    ]

# Martini Quotes
martiniq = [
    ("Robert Benchley",
     "Why don't you get out of that wet coat and into a dry martini?"),
    ("James Thurber",
     "One martini is all right, two is two many, three is not enough"),
    ("Philip Larkin",
     "The chromatic scale is what you use to give the effect of drinking a quinine martini and having an enema simultaneously."),
    ("William Emerson, Jr.",
     "And when that first martini hits the liver like a silver bullet, there is a sigh of contentment that can be heard in Dubuque."),
    ("Alec Waugh",
     "I am prepared to believe that a dry martini slightly impairs the palate, but think what it does for the soul."),
    ("Gerald R. Ford",
     "The three-martini lunch is the epitome of American efficiency. Where else can you get an earful, a bellyful and a snootful at the same time?"),
    ("P. G. Wodehouse",
     "He was white and shaken, like a dry martini."),
    ]

desc = ""


def help():
    """Print help text and list of options and end the program."""
    import sys
    for item in options:
        if type(item) == str:
            print(item)
    for item in options:
        if type(item) != str:
            print("%10s  %s" % (item[0], item[1].description))
    print()
    sys.exit()
##############################
## 2 # COMMAND LINE PARSING ##  -> @CMD <-
##############################
import sys, logging, os, inspect


def str2atom(a):
    """ Helper function to parse atom strings given on the command line:
    resid, resname/resid, chain/resname/resid, resname/resid/atom,
    chain/resname/resid/atom, chain//resid, chain/resname/atom.
    Returns (atom,resname,resid,chain)"""
    a = a.split("/")
    if len(a) == 1:  # Only a residue number (123):
        return (None, None, int(a[0]), None)
    if len(a) == 2:  # Residue name and number (CYS/123):
        return (None, a[0], int(a[1]), None)
    if len(a) == 3:
        if a[2].isdigit():  # Chain, residue name, residue number (A/CYS/123)
            return (None, a[1], int(a[2]), a[0])
        else:  # Residue name, residue number, atom name
            return (a[2], a[0], int(a[1]), None)
    return (a[3], a[1], int(a[2]), a[0])


def option_parser(args, options, lists, version=0):

    # Check whether there is a request for help
    if '-h' in args or '--help' in args:
        help()

    # Convert the option list to a dictionary, discarding all comments
    options = dict([i for i in options if not type(i) == str])

    # This information we would like to print to some files,
    # so let's put it in our information class
    options['Version']             = version
    options['Arguments']           = args[:]

    while args:
        ar = args.pop(0)
        options[ar].setvalue([args.pop(0) for i in range(options[ar].num)])

    ## LOGGING ##
    # Set the log level and communicate which options are set and what is happening
    # If 'Verbose' is set, change the logger level
    logLevel = options["-v"] and logging.DEBUG or logging.INFO
    logging.basicConfig(format='%(levelname)-7s    %(message)s', level=logLevel)

    logging.info('MARTINIZE, script version %s'%version)
    logging.info('If you use this script please cite:')
    logging.info('de Jong et al., J. Chem. Theory Comput., 2013, DOI:10.1021/ct300646g')

    # To make the program flexible, the forcefield parameters are defined
    # for multiple forcefield. We first check for a FF file in the current directory.
    # Next we check for the FF in globals (for catenated scripts). 
    # Next we check in at the location of the script and the subdiretory FF.
    try:
        options['ForceField'] = globals()[options['-ff'].value.lower()]()
    except KeyError:
        try:
            _tmp = __import__(options['-ff'].value.lower()+"_ff")
            options['ForceField'] = getattr(_tmp, options['-ff'].value.lower())()
        except ImportError:
            try:
                # We add the directory where the script resides and a possible "ForceFields" directory to the search path
                # realpath() will make your script run, even if you symlink it :)
                cmd_folder = os.path.realpath(os.path.dirname(inspect.getfile(inspect.currentframe())))
                if cmd_folder not in sys.path:
                    sys.path.insert(0, cmd_folder)
                # use this if you want to include modules from a subfolder
                cmd_subfolder = os.path.realpath(os.path.dirname(inspect.getfile(inspect.currentframe()))) + "/ForceFields"
                if cmd_subfolder not in sys.path:
                     sys.path.insert(0, cmd_subfolder)
                _tmp = __import__(options['-ff'].value.lower()+"_ff")
                options['ForceField'] = getattr(_tmp, options['-ff'].value.lower())()
            except:
                logging.error("Forcefield '%s' can not be loaded." % (options['-ff']))
                sys.exit()
        
    #    if os.path.exists(options['-ff'].value.lower()+'_ff.py'):
    #        _tmp = __import__(options['-ff'].value.lower()+"_ff")
    #        options['ForceField'] = getattr(_tmp, options['-ff'].value.lower())()
    #    elif os.path.exists('ForceFields/'+options['-ff'].value.lower()+'_ff.py'):
    #        _tmp = __import__("ForceFields."+options['-ff'].value.lower()+'_ff',fromlist="ForceFields")
    #        options['ForceField'] = getattr(_tmp, options['-ff'].value.lower())()
    #    elif options['-ff'].value.lower() in globals():
    #        options['ForceField'] = globals()[options['-ff'].value.lower()]()
    #    else:
    #        logging.error("Forcefield '%s' can not be found."%(options['-ff']))
    #        sys.exit()
    #except:
    #    logging.error("Forcefield '%s' can not be loaded."%(options['-ff']))
    #    sys.exit()

    # Process the raw options from the command line
    # Boolean options are set to more intuitive variables
    options['Collagen']            = options['-collagen']
    options['chHIS']               = options['-his']
    options['ChargesAtBreaks']     = options['-cb']
    options['NeutralTermini']      = options['-nt']
    options['ExtendedDihedrals']   = options['-ed']
    options['RetainHETATM']        = False  # options['-hetatm']
    options['SeparateTop']         = options['-sep']
    options['MixedChains']         = False  # options['-mixed']
    options['ElasticNetwork']      = options['-elastic']

    # Parsing of some other options into variables
    options['ElasticMaximumForce'] = options['-ef'].value
    options['ElasticMinimumForce'] = options['-em'].value
    options['ElasticLowerBound']   = options['-el'].value
    options['ElasticUpperBound']   = options['-eu'].value
    options['ElasticDecayFactor']  = options['-ea'].value
    options['ElasticDecayPower']   = options['-ep'].value
    options['ElasticBeads']        = options['-eb'].value.split(',')
    options['PosResForce']         = options['-pf'].value

    options['PosRes']              = [i.lower() for i in options['-p'].value.split(",")]
    if "none"     in options['PosRes']: options['PosRes'] = []
    if "backbone" in options['PosRes']: options['PosRes'].append("BB")

    if options['ForceField'].ElasticNetwork:
        # Some forcefields, like elnedyn, always use an elatic network.
        # This is set in the forcefield file, with the parameter ElasticNetwork.
        options['ElasticNetwork'] = True

    # Merges, links and cystines
    options['mergeList'] = "all" in lists['merges'] and ["all"] or [i.split(",") for i in lists['merges']]

    # Process links
    linkList   = []
    linkListCG = []
    for i in lists['links']:
        ln     = i.split(",")
        # str2atom returns (atomname,resname,resid,chain)
        a, b   = str2atom(ln[0]), str2atom(ln[1])        
        if len(ln) > 3:  # Bond with given length and force constant
            bl, fc = (ln[2] and float(ln[2]) or None, float(ln[3]))
        elif len(a) == 3:  # Constraint at given distance
            bl, fc = float(a[2]), None
        else:  # Constraint at distance in structure
            bl, fc = None, None
        # Store the link, but do not list the atom name in the
        # atomistic link list. Otherwise it will not get noticed
        # as a valid link when checking for merging chains
        linkList.append(((None, a[1], a[2], a[3]), (None, b[1], b[2], b[3])))
        linkListCG.append((a, b, bl, fc))

    # Cystines
    # This should be done for all special bonds listed in the _special_ dictionary
    CystineCheckBonds = False   # By default, do not detect cystine bridges
    CystineMaxDist2   = (10*0.22)**2  # Maximum distance (A) for detection of SS bonds
    for i in lists['cystines']:
        if i.lower() == "auto":
            CystineCheckBonds = True
        elif i.replace(".", "").isdigit():
            CystineCheckBonds = True
            CystineMaxDist2   = (10*float(i))**2
        else:
            # This item should be a pair of cysteines
            cysA, cysB = [str2atom(j) for j in i.split(",")]
            # Internally we handle the residue number shifted by ord(' ')<<20.
            # We have to add this to the cys-residue numbers given here as well.
            constant = 32 << 20
            linkList.append((("SG", "CYS", cysA[2]+constant, cysA[3]),
                            ("SG", "CYS", cysB[2]+constant, cysB[3])))
            linkListCG.append((("SC1", "CYS", cysA[2]+constant, cysA[3]),
                              ("SC1", "CYS", cysB[2]+constant, cysB[3]), -1, -1))

    # Now we have done everything to it, we can add Link/cystine related stuff to options
    # 'multi' is not stored anywhere else, so that we also add
    options['linkList']          = linkList
    options['linkListCG']        = linkListCG
    options['CystineCheckBonds'] = CystineCheckBonds
    options['CystineMaxDist2']   = CystineMaxDist2
    options['multi']             = lists['multi']

    logging.info("Chain termini will%s be charged."%(options['NeutralTermini'] and " not" or ""))

    logging.info("Residues at chain breaks will%s be charged."%((not options['ChargesAtBreaks']) and " not" or ""))

    if 'ForceField' in options:
        logging.info("The %s forcefield will be used."%(options['ForceField'].name))
    else:
        logging.error("Forcefield '%s' has not been implemented."%(options['-ff']))
        sys.exit()

    if options['ExtendedDihedrals']:
        logging.info('Dihedrals will be used for extended regions. (Elastic bonds may be more stable)')
    else:
        logging.info('Local elastic bonds will be used for extended regions.')

    if options['PosRes']:
        logging.info("Position restraints will be generated.")
        logging.warning("Position restraints are only enabled if -DPOSRES is set in the MDP file")

    if options['MixedChains']:
        logging.warning("So far no parameters for mixed chains are available. This might crash the program!")

    if options['RetainHETATM']:
        logging.warning("I don't know how to handle HETATMs. This will probably crash the program.")

    return options
#################################################
## 3 # HELPER FUNCTIONS, CLASSES AND SHORTCUTS ##  -> @FUNC <-
#################################################

import math

#----+------------------+
## A | STRING FUNCTIONS |
#----+------------------+


# Split a string
def spl(x):
    return x.split()


# Split each argument in a list
def nsplit(*x):
    return [i.split() for i in x]


# Make a dictionary from two lists
def hash(x, y):
    return dict(list(zip(x, y)))


# Function to reformat pattern strings
def pat(x, c="."):
    return x.replace(c, "\x00").split()


# Function to generate formatted strings according to the argument type
def formatString(i):
    if type(i) == str:
        return i
    if type(i) == int:
        return "%5d" % i
    if type(i) == float:
        return "%8.5f" % i
    else:
        return str(i)


#----+----------------+
## B | MATH FUNCTIONS |
#----+----------------+


def cos_angle(a, b):
    p = sum([i*j for i, j in zip(a, b)])
    q = math.sqrt(sum([i*i for i in a])*sum([j*j for j in b]))
    return min(max(-1, p/q), 1)


def norm2(a):
    return sum([i*i for i in a])


def norm(a):
    return math.sqrt(norm2(a))


def distance2(a, b):
    return (a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2
##########################
## 4 # FG -> CG MAPPING ##  -> @MAP <-
##########################


dnares3 = " DA DC DG DT"
dnares1 = " dA dC dG dT"
rnares3 = "  A  C  G  U"
rnares1 = " rA rC rG rU"

## Sugar naming
# Common IUPAC 3-letter code 
sugarres3 = " Ara Lyx Rib Xyl All Alt Gal Glc Gul Ido Man Tal Fru Psi Sor Tag Fuc Qui Rha GalA GlcA IdoA GalNAc GlcNAc ManNAc Neu5Ac  KDN  KDO Neu5Gc"
# GLYCAM 1-letter code; here only upper case letters which indicate D configuration.
sugarres1 = "   A   D   R   X   N   E   L   G   K   I   M   T   C   P   B   J   F   Q   H    O    Z    U      V      Y      W      S (KN) (KO)   (SG)"

# Amino acid nucleic acid codes:
# The naming (AA and '3') is not strictly correct when adding DNA/RNA, but we keep it like this for consistency.
# HIE added (GLYCAM convention for epsilon-protonated HIS)
AA3     = spl("TRP TYR PHE HIE HIS HIH ARG LYS CYS ASP GLU ILE LEU MET ASN PRO HYP GLN SER THR VAL ALA GLY"+dnares3+rnares3+sugarres3) #@#
AA1     = spl("  W   Y   F   H   H   H   R   K   C   D   E   I   L   M   N   P   O   Q   S   T   V   A   G"+dnares1+rnares1+sugarres1) #@#

# Dictionaries for conversion from one letter code to three letter code v.v.
AA123, AA321 = hash(AA1, AA3), hash(AA3, AA1)

# Residue classes:
protein = AA3[:23]   # amino acids only
water   = spl("HOH SOL TIP")
lipids  = spl("DPP DHP DLP DMP DSP POP DOP DAP DUP DPP DHP DLP DMP DSP PPC DSM DSD DSS")
nucleic = spl("DAD DCY DGU DTH ADE CYT GUA THY URA DA DC DG DT")

# from GLYCAM documentation
glycamCode = dict(
    [(i,"GlcA")   for i in spl("0ZA 0ZB 1ZA 1ZB 2ZA 2ZB 3ZA 3ZB 4ZA 4ZB ZZA ZZB YZA YZB WZA WZB TZA TZB 0zA 0zB 1zA 1zB 2zA 2zB 3zA 3zB 4zA 4zB ZzA ZzB YzA YzB WzA WzB TzA TzB 0ZBP")] +
    [(i,"GlcN")   for i in spl("0YN 0Yn 0YNP 0YnP 0YS 0Ys 3YS 3Ys 4YS 4Ys 6YS 6Ys QYS QYs UYS UYs VYS VYs WYS WYs 0yS 0ys 3yS 3ys 4yS 4ys")] +
    [(i,"GlcNAc") for i in spl("0YA 0YB 1YA 1YB 3YA 3YB 4YA 4YB 6YA 6YB WYA WYB VYA VYB UYA UYB QYA QYB 0yA 0yB 1yA 1yB 3yA 3yB 4yA 4yB 6yA 6yB WyA WyB VyA VyB UyA UyB QyA QyB")] +    
    [(i,"Xyl")    for i in spl("0XA 0XB 1XA 1XB 2XA 2XB 3XA 3XB 4XA 4XB ZXA ZXB YXA YXB WXA WXB TXA TXB 0XD 0XU 1XD 1XU 2XD 2XU 3XD 3XU 5XD 5XU ZXD ZXU 0xA 0xB 1xA 1xB 2xA 2xB 3xA 3xB 4xA 4xB ZxA ZxB YxA YxB WxA WxB TxA TxB 0xD 0xU 1xD 1xU 2xD 2xU 3xD 3xU 5xD 5xU ZxD ZxU")] + # Xyl
    [(i,"ManNAc") for i in spl("0WA 0WB 1WA 1WB 3WA 3WB 4WA 4WB 6WA 6WB WWA WWB VWA VWB UWA UWB QWA QWB 0wA 0wB 1wA 1wB 3wA 3wB 4wA 4wB 6wA 6wB WwA WwB VwA VwB UwA UwB QwA QwB")] + # ManNAc
    [(i,"GalNAc") for i in spl("0VA 0VB 1VA 1VB 3VA 3VB 4VA 4VB 6VA 6VB WVA WVB VVA VVB UVA UVB QVA QVB 0vA 0vB 1vA 1vB 3vA 3vB 4vA 4vB 6vA 6vB WvA WvB VvA VvB UvA UvB QvA QvB")] + # GalNAc
    [(i,"IdoA")   for i in spl("0UA 0UB 1UA 1UB 2UA 2UB 3UA 3UB 4UA 4UB ZUA ZUB YUA YUB WUA WUB TUA TUB 0uA 0uB 1uA 1uB 2uA 2uB 3uA 3uB 4uA 4uB ZuA ZuB YuA YuB WuA WuB TuA TuB YuAP")] + # IdoA
    [(i,"Neu5Ac") for i in spl("0SA 0SB 4SA 4SB 7SA 7SB 8SA 8SB 9SA 9SB ASA ASB BSA BSB CSA CSB DSA DSB ESA ESB FSA FSB GSA GSB HSA HSB ISA ISB JSA JSB KSA KSB 0sA 0sB 4sA 4sB 7sA 7sB 8sA 8sB 9sA 9sB AsA AsB BsA BsB CsA CsB DsA DsB EsA EsB FsA FsB GsA GsB HsA HsB IsA IsB JsA JsB KsA KsB")] + # Neu5Ac
    [(i,"GalA")   for i in spl("0OA 0OB 1OA 1OB 2OA 2OB 3OA 3OB 4OA 4OB ZOA ZOB YOA YOB WOA WOB TOA TOB 0oA 0oB 1oA 1oB 2oA 2oB 3oA 3oB 4oA 4oB ZoA ZoB YoA YoB WoA WoB ToA ToB")] + # GalA
    [(i,"Man")    for i in spl("0MA 0MB 1MA 1MB 2MA 2MB 3MA 3MB 4MA 4MB 6MA 6MB ZMA ZMB YMA YMB XMA XMB WMA WMB VMA VMB UMA UMB TMA TMB SMA SMB RMA RMB QMA QMB PMA PMB 0mA 0mB 1mA 1mB 2mA 2mB 3mA 3mB 4mA 4mB 6mA 6mB ZmA ZmB YmA YmB XmA XmB WmA WmB VmA VmB UmA UmB TmA TmB SmA SmB RmA RmB QmA QmB PmA PmB")] + # Man
    [(i,"Gal")    for i in spl("0LA 0LB 1LA 1LB 2LA 2LB 3LA 3LB 4LA 4LB 6LA 6LB ZLA ZLB YLA YLB XLA XLB WLA WLB VLA VLB ULA ULB TLA TLB SLA SLB RLA RLB QLA QLB PLA PLB 0lA 0lB 1lA 1lB 2lA 2lB 3lA 3lB 4lA 4lB 6lA 6lB ZlA ZlB YlA YlB XlA XlB WlA WlB VlA VlB UlA UlB TlA TlB SlA SlB RlA RlB QlA QlB PlA PlB")] + # Gal
    [(i,"Glc")    for i in spl("0GA 0GB 1GA 1GB 2GA 2GB 3GA 3GB 4GA 4GB 6GA 6GB ZGA ZGB YGA YGB XGA XGB WGA WGB VGA VGB UGA UGB TGA TGB SGA SGB RGA RGB QGA QGB PGA PGB 0gA 0gB 1gA 1gB 2gA 2gB 3gA 3gB 4gA 4gB 6gA 6gB ZgA ZgB YgA YgB XgA XgB WgA WgB VgA VgB UgA UgB TgA TgB SgA SgB RgA RgB QgA QgB PgA PgB")] + # Glc
    [(i,"Fuc")    for i in spl("0FA 0FB 1FA 1FB 2FA 2FB 3FA 3FB 4FA 4FB ZFA ZFB YFA YFB WFA WFB TFA TFB 0fA 0fB 1fA 1fB 2fA 2fB 3fA 3fB 4fA 4fB ZfA ZfB YfA YfB WfA WfB TfA TfB")] + # Fuc
    [(i,"Fru")    for i in spl("0CA 0CB 1CA 1CB 2CA 2CB 3CA 3CB 4CA 4CB 5CA 5CB WCA WCB 0CD 0CU 1CD 1CU 2CD 2CU 3CD 3CU 4CD 4CU 6CD 6CU WCD WCU VCD VCU UCD UCU QCD QCU 0cA 0cB 1cA 1cB 2cA 2cB 3cA 3cB 4cA 4cB 5cA 5cB WcA WcB 0cD 0cU 1cD 1cU 2cD 2cU 3cD 3cU 4cD 4cU 6cD 6cU WcD WcU VcD VcU UcD UcU QcD QcU")] + # Fru
    [(i,"Ara")    for i in spl("0AA 0AB 1AA 1AB 2AA 2AB 3AA 3AB 4AA 4AB ZAA ZAB YAA YAB WAA WAB TAA TAB 0AD 0AU 1AD 1AU 2AD 2AU 3AD 3AU 5AD 5AU ZAD ZAU 0aA 0aB 1aA 1aB 2aA 2aB 3aA 3aB 4aA 4aB ZaA ZaB YaA YaB WaA WaB TaA TaB 0aD 0aU 1aD 1aU 2aD 2aU 3aD 3aU 5aD 5aU ZaD ZaU")] + # Ara
    [(i,"Redend") for i in spl("ROH")] # Reducing end OH group
    )


    
residueTypes = dict(
    [(i, "Protein") for i in protein ] +
    [(i, "Water")   for i in water   ] +
    [(i, "Lipid")   for i in lipids  ] +
    [(i, "Nucleic") for i in nucleic ] +
    [(i, "Glycan")  for i in glycamCode ]
    )


class CoarseGrained(object):
    # Class for mapping an atomistic residue list to a coarse-grained one
    # Should get an __init__ function taking a residuelist, atomlist, Pymol selection or ChemPy model
    # The result should be stored in a list-type attribute
    # The class should have pdbstr and grostr methods

    # Standard mapping groups
    # Protein backbone
    bb        = "N CA C O H H1 H2 H3 O1 O2"                                                                    #@#
    # Lipid tails
    palmitoyl1    = nsplit("C1B C1C C1D C1E", "C1F C1G C1H C1I", "C1J C1K C1L C1M", "C1N C1O C1P")              #@#
    palmitoyl2    = nsplit("C2B C2C C2D C2E", "C2F C2G C2H C2I", "C2J C2K C2L C2M", "C2N C2O C2P")              #@#
    oleyl1        = nsplit("C1B C1C C1D C1E", "C1F C1G C1H", "C1I C1J", "C1K C1L C1M C1N", "C1O C1P C1Q C1R")   #@#
    oleyl2        = nsplit("C2B C2C C2D C2E", "C2F C2G C2H", "C2I C2J", "C2K C2L C2M C2N", "C2O C2P C2Q C2R")   #@#
    #lauroyl1      = []
    #stearoyl1     = []
    #arachidonoyl1 = []
    #linoleyl1     = []
    #hexanoyl1     = []
    # Lipid head groups
    #phoshpatidylcholine      =
    phosphatydilethanolamine = nsplit("N H1 H2 H3 CA", "CB P OA OB OC OD", "CC CD OG C2A OH", "CE OE C1A OF")      #@#
    phosphatidylglycerol     = nsplit("H1 O1 CA H2 O2 CB", "CC P OA OB OC OD", "CD CE OG C2A OH", "CF OE C1A OF")  #@#
    #phosphatidylserine       =

    dna_bb = "P OP1 OP2 O5' O3'", "C5' O4' C4'", "C3' O3' C2' C1'"

    # This is the mapping dictionary
    # For each residue it returns a list, each element of which
    # lists the atom names to be mapped to the corresponding bead.
    # The order should be the standard order of the coarse grained
    # beads for the residue. Only atom names matching with those
    # present in the list of atoms for the residue will be used
    # to determine the bead position. This adds flexibility to the
    # approach, as a single definition can be used for different
    # states of a residue (e.g., GLU/GLUH).
    # For convenience, the list can be specified as a set of strings,
    # converted into a list of lists by 'nsplit' defined above.
    mapping = {
        "ALA":  nsplit(bb + " CB"),
        "CYS":  nsplit(bb, "CB SG"),
        "ASP":  nsplit(bb, "CB CG OD1 OD2"),
        "GLU":  nsplit(bb, "CB CG CD OE1 OE2"),
        "PHE":  nsplit(bb, "CB CG CD1 HD1", "CD2 HD2 CE2 HE2", "CE1 HE1 CZ HZ"),
        "GLY":  nsplit(bb),
        "HIS":  nsplit(bb, "CB CG", "CD2 HD2 NE2 HE2", "ND1 HD1 CE1 HE1"),
        "HIH":  nsplit(bb, "CB CG", "CD2 HD2 NE2 HE2", "ND1 HD1 CE1 HE1"),     # Charged Histidine.
        "ILE":  nsplit(bb, "CB CG1 CG2 CD CD1"),
        "LYS":  nsplit(bb, "CB CG CD", "CE NZ HZ1 HZ2 HZ3"),
        "LEU":  nsplit(bb, "CB CG CD1 CD2"),
        "MET":  nsplit(bb, "CB CG SD CE"),
        "ASN":  nsplit(bb, "CB CG ND1 ND2 OD1 OD2 HD11 HD12 HD21 HD22"),
        "PRO":  nsplit(bb, "CB CG CD"),
        "HYP":  nsplit(bb, "CB CG CD OD"),
        "GLN":  nsplit(bb, "CB CG CD OE1 OE2 NE1 NE2 HE11 HE12 HE21 HE22"),
        "ARG":  nsplit(bb, "CB CG CD", "NE HE CZ NH1 NH2 HH11 HH12 HH21 HH22"),
        "SER":  nsplit(bb, "CB OG HG"),
        "THR":  nsplit(bb, "CB OG1 HG1 CG2"),
        "VAL":  nsplit(bb, "CB CG1 CG2"),
        "TRP":  nsplit(bb, "CB CG CD2", "CD1 HD1 NE1 HE1 CE2", "CE3 HE3 CZ3 HZ3", "CZ2 HZ2 CH2 HH2"),
        "TYR":  nsplit(bb, "CB CG CD1 HD1", "CD2 HD2 CE2 HE2", "CE1 HE1 CZ OH HH"),
        "POPE": phosphatydilethanolamine + palmitoyl1 + oleyl2,
        "DOPE": phosphatydilethanolamine + oleyl1     + oleyl2,
        "DPPE": phosphatydilethanolamine + palmitoyl1 + palmitoyl2,
        "POPG": phosphatidylglycerol     + palmitoyl1 + oleyl2,
        "DOPG": phosphatidylglycerol     + oleyl1     + oleyl2,
        "DPPG": phosphatidylglycerol     + palmitoyl1 + palmitoyl2,
        "DA": nsplit("P OP1 OP2 O5' O3' O1P O2P", "C5' O4' C4'", "C3' C2' C1'", "N9 C4", "C8 N7 C5", "C6 N6 N1", "C2 N3"),
        "DG": nsplit("P OP1 OP2 O5' O3' O1P O2P", "C5' O4' C4'", "C3' C2' C1'", "N9 C4", "C8 N7 C5", "C6 O6 N1", "C2 N2 N3"),
        "DC": nsplit("P OP1 OP2 O5' O3' O1P O2P", "C5' O4' C4'", "C3' C2' C1'", "N1 C6", "C5 C4 N4", "N3 C2 O2"),
        "DT": nsplit("P OP1 OP2 O5' O3' O1P O2P", "C5' O4' C4'", "C3' C2' C1'", "N1 C6", "C5 C4 O4 C7 C5M", "N3 C2 O2"),
        }
    
    # Sugars
    sugarmapping = {
        "GlcNAc": nsplit("C1 H1 C2 H2 N2 H2N", "C3 H3 O3 H3O C4 H4 O4 H4O", "C5 H5 O5 C6 H61 H62 O6 H6O", "C2N O2N CME H1M H2M H3M"),
        "GalNAc": nsplit("C1 H1 C2 H2 N2 H2N", "C3 H3 O3 H3O C4 H4 O4 H4O", "C5 H5 O5 C6 H61 H62 O6 H6O", "C2N O2N CME H1M H2M H3M"),
        "Man":    nsplit("C1 H1 C2 H2 O2 H2O", "C3 H3 O3 H3O C4 H4 O4 H4O", "C5 H5 O5 C6 H61 H62 O6 H6O"),
        "Gal":    nsplit("C1 H1 C2 H2 O2 H2O", "C3 H3 O3 H3O C4 H4 O4 H4O", "C5 H5 O5 C6 H61 H62 O6 H6O"),
        "Glc":    nsplit("C1 H1 C2 H2 O2 H2O", "C3 H3 O3 H3O C4 H4 O4 H4O", "C5 H5 O5 C6 H61 H62 O6 H6O"),
        "Fuc":    nsplit("C1 H1 C2 H2 O2 H2O", "C3 H3 O3 H3O C4 H4 O4 H4O", "C5 H5 O5 C6 H61 H62 H63"),
        "Neu5Ac": nsplit("C1 O1A O1B C2", "C3 H3E H3A C4 H4 O4 H4O C5 H5", "C6 H6 O6 C7 H7 O7 H7O", "N5 H5N C5N O5N CME H1M H2M H3M", "C8 H8 O8 H8O C9 H9S H9R O9 H9O"),
        }
    
    redends = {
        "ROH":    nsplit("O1 HO1")
        }
    
    
    # Generic names for side chain beads
    residue_bead_names = spl("BB SC1 SC2 SC3 SC4")
    # Generic names for DNA beads
    residue_bead_names_dna = spl("BB1 BB2 BB3 SC1 SC2 SC3 SC4")

    # This dictionary contains the bead names for all residues,
    # following the order in 'mapping'
    names  = {
        "POPE": "NH3 PO4 GL1 GL2 C1A C2A C3A C4A C1B C2B D3B C4B C5B".split(),
        "POPG": "GLC PO4 GL1 GL2 C1A C2A C3A C4A C1B C2B D3B C4B C5B".split()
        }
    # Add default bead names for all amino acids
    names.update([(i, ("BB", "SC1", "SC2", "SC3", "SC4")) for i in AA3])

    # Add the default bead names for all DNA nucleic acids
    names.update([(i, ("BB1", "BB2", "BB3", "SC1", "SC2", "SC3", "SC4")) for i in nucleic])
    
    # Bead names for sugars
    beadnames_sugar = {
        "GlcNAc": "C12 C34 C56 2Ac".split(),
        "GalNAc": "C12 C34 C56 2Ac".split(),
        "Man":    "C12 C34 C56".split(),
        "Gal":    "C12 C34 C56".split(),
        "Glc":    "C12 C34 C56".split(),
        "Fuc":    "C12 C34 C56".split(),
        "Neu5Ac": "C12 345 C67 NAc C89".split(),
        }

    # This dictionary allows determining four letter residue names
    # for ones specified with three letters, e.g., resulting from
    # truncation to adhere to the PDB format.
    # Each entry returns a prototypical test, given as a string,
    # and the residue name to be applied if eval(test) is True.
    # This is particularly handy to determine lipid types.
    # The test assumes there is a local or global array 'atoms'
    # containing the atom names of the residue in correct order.
    restest = {
        "POP": [('atoms[0] == "CA"', "POPG"),
                ('atoms[0] == "N"',  "POPE")]
        }

    # Crude mass for weighted average. No consideration of united atoms.
    # This will probably give only minor deviations, while also giving less headache
    mass = {'H': 1, 'C': 12, 'N': 14, 'O': 16, 'S': 32, 'P': 31, 'M': 0}


# Determine average position for a set of weights and coordinates
# This is a rather specific function that requires a list of items
# [(m,(x,y,z),id),..] and returns the weighted average of the
# coordinates [x, y, z] and the list of ids mapped to this bead
def aver(b):
    mwx, ids = list(zip(*[((m*x, m*y, m*z), i) for m, (x, y, z), i in b]))      # Weighted coordinates
    tm  = sum(next(zip(*b)))                                                 # Sum of weights
    return [sum(i)/tm for i in zip(*mwx)], ids                            # Centre of mass


# Return the CG beads for an atomistic residue, using the mapping specified above
# The residue 'r' is simply a list of atoms, and each atom is a list:
# [ name, resname, resid, chain, x, y, z ]
def map(r, ca2bb = False):
    p = CoarseGrained.mapping[r[0][1]]     # Mapping for this residue
    if ca2bb: p[0] = ["CA"]                # Elnedyn maps BB to CA, ca2bb is False or True
    # Get the name, mass and coordinates for all atoms in the residue
    a = [(i[0], CoarseGrained.mass.get(i[0][0], 0), i[4:]) for i in r]
    # Store weight, coordinate and index for atoms that match a bead
    q = [[(m, coord, a.index((atom, m, coord))) for atom, m, coord in a if atom in i] for i in p]

    # Bead positions
    return list(zip(*[aver(i) for i in q]))

def sugarmap(sugar, red_end=False):
    
    # List of list of atom names to be mapped onto CG beads for this sugar 
    atlist = CoarseGrained.sugarmapping[sugar.name]

    # Get (atmname, mass, [x, y, z]) for all atoms in the residue. Mass lookup works because
    # the first letter of the atom name defines its chemical identity (C, H, O, N, S, P, ...)
    try:
        a = [(i[0], CoarseGrained.mass[i[0][0]], i[4:]) for i in sugar.residue]
    except KeyError:
        logging.warning('Atom detected with unknown mass! Check CoarseGrained.mass dictionary.')
    
    # For reducing end sugars red_end is the atom list (Residue object) of the ROH residue
    # (atom name, res name, res id, chain id, x, y, z)
    if red_end:
        
        # First bead needs to include the ROH atoms
        for atom in CoarseGrained.redends[red_end[0][1]][0]:            
            if atom not in atlist[0]:
                atlist[0].append(atom)
        # Add (atmname, mass, [x, y, z]) of all ROH atoms to the atom list of the sugar
        a.extend([(i[0], CoarseGrained.mass[i[0][0]], i[4:]) for i in red_end])
    
    # List of list of (mass, coordinates, index) of atoms (inner list) belonging to the beads (outer list) of the residue
    q = [[(m, coord, a.index((atom, m, coord))) for atom, m, coord in a if atom in bead] for bead in atlist]
    
    # aver(bead) returns ([comx, comy, comz], [atom indices])
    # Function returns a list of two tuples: lists of com coords and atom index lists for every bead in the residue
    #           bead 1                 bead2                        bead 1                     bead2
    # [([comx1, comy1, comz1], [comx2, comy2, comz2], ...), ([atom indices bead1], [atom indices bead2], ...)]
    return list(zip(*[aver(bead) for bead in q]))


# Mapping for index file
def mapIndex(r, ca2bb = False):
    p = CoarseGrained.mapping[r[0][1]]        # Mapping for this residue
    if ca2bb: p[0] = ["CA"]                   # Elnedyn maps BB to CA, ca2bb is False or True
    # Get the name, mass and coordinates for all atoms in the residue
    a = [(i[0], CoarseGrained.mass.get(i[0][0], 0), i[4:]) for i in r]
    # Store weight, coordinate and index for atoms that match a bead
    return [[(m, coord, a.index((atom, m, coord))) for atom, m, coord in a if atom in i] for i in p]
#############################
## 5 # SECONDARY STRUCTURE ##  -> @SS <-
#############################
import logging, os, sys
import subprocess as subp

#----+--------------------------------------+
## A | SECONDARY STRUCTURE TYPE DEFINITIONS |
#----+--------------------------------------+

# This table lists all coarse grained secondary structure types
# The following are matched lists. Make sure they stay matched.
# The lists do not need to be of the same length. The longer list
# will be truncated when combined with a shorter list, e.g. with
# dihedral definitions, which are not present for coil and termini
#
ss_names = {
 "F": "Collagenous Fiber",                                                                  #@#
 "E": "Extended structure (beta sheet)",                                                    #@#
 "H": "Helix structure",                                                                    #@#
 "1": "Helix start (H-bond donor)",                                                         #@#
 "2": "Helix end (H-bond acceptor)",                                                        #@#
 "3": "Ambivalent helix type (short helices)",                                              #@#
 "T": "Turn",                                                                               #@#
 "S": "Bend",                                                                               #@#
 "C": "Coil",                                                                               #@#
}

bbss = list(ss_names.keys())
bbss = spl("  F     E     H     1     2     3     T     S     C")  # SS one letter


# The following dictionary contains secondary structure types as assigned by
# different programs. The corresponding Martini secondary structure types are
# listed in cgss
#
# NOTE:
#  Each list of letters in the dictionary ss should exactly match the list
#  in cgss.
#
ssdefs = {
    "dssp":  list(".HGIBETSC~"),             # DSSP one letter secondary structure code     #@#
    "pymol": list(".H...S...L"),             # Pymol one letter secondary structure code    #@#
    "gmx":   list(".H...ETS.C"),             # Gromacs secondary structure dump code        #@#
    "self":  list("FHHHEETSCC")              # Internal CG secondary structure codes        #@#
}
cgss     =   list("FHHHEETSCC")              # Corresponding CG secondary structure types   #@#


#----+-------------------------------------------+
## B | SECONDARY STRUCTURE PATTERN SUBSTITUTIONS |
#----+-------------------------------------------+


# For all structure types specific dihedrals may be used if four or
# more consecutive residues are assigned that type.

# Helix start and end regions are special and require assignment of
# specific types. The following pattern substitutions are applied
# (in the given order). A dot matches any other type.
# Patterns can be added to the dictionaries. This only makes sense
# if for each key in patterns there is a matching key in pattypes.
patterns = {
    "H": pat(".H. .HH. .HHH. .HHHH. .HHHHH. .HHHHHH. .HHHHHHH. .HHHH HHHH.")                #@#
}
pattypes = {
    "H": pat(".3. .33. .333. .3333. .13332. .113322. .1113222. .1111 2222.")                #@#
}


#----+----------+
## C | INTERNAL |
#----+----------+


# Pymol Colors
#          F   E   H   1   2   3   T   S   C
ssnum  = (13,  4,  2,  2,  2,  2,  6, 22,  0)                                             #@#

# Dictionary returning a number for a given type of secondary structure
# This can be used for setting the b-factor field for coloring
ss2num = hash(bbss, ssnum)


# List of programs for which secondary structure definitions can be processed
programs = list(ssdefs.keys())


# Dictionaries mapping ss types to the CG ss types
ssd = dict([(i, hash(ssdefs[i], cgss)) for i in programs])


# From the secondary structure dictionaries we create translation tables
# with which all secondary structure types can be processed. Anything
# not listed above will be mapped to C (coil).
# Note, a translation table is a list of 256 characters to map standard
# ascii characters to.
def tt(program):
    return "".join([ssd[program].get(chr(i), "C") for i in range(256)])


# The translation table depends on the program used to obtain the
# secondary structure definitions
sstt = dict([(i, tt(i)) for i in programs])


# The following translation tables are used to identify stretches of
# a certain type of secondary structure. These translation tables have
# every character, except for the indicated secondary structure, set to
# \x00. This allows summing the sequences after processing to obtain
# a single sequence summarizing all the features.
null = "\x00"
sstd = dict([(i, ord(i)*null+i+(255-ord(i))*null) for i in cgss])


# Pattern substitutions
def typesub(seq, patterns, types):
    seq = null+seq+null
    for i, j in zip(patterns, types):
        seq = seq.replace(i, j)
    return seq[1:-1]


# The following function translates a string encoding the secondary structure
# to a string of corresponding Martini types, taking the origin of the
# secondary structure into account, and replacing termini if requested.
def ssClassification(ss, program="dssp"):
    # Translate dssp/pymol/gmx ss to Martini ss
    ss  = ss.translate(sstt[program])
    # Separate the different secondary structure types
    sep = dict([(i, ss.translate(sstd[i])) for i in list(sstd.keys())])
    # Do type substitutions based on patterns
    # If the ss type is not in the patterns lists, do not substitute
    # (use empty lists for substitutions)
    typ = [typesub(sep[i], patterns.get(i, []), pattypes.get(i, [])) for i in list(sstd.keys())]
    # Translate all types to numerical values
    typ = [[ord(j) for j in list(i)] for i in typ]
    # Sum characters back to get a full typed sequence
    typ = "".join([chr(sum(i)) for i in zip(*typ)])
    # Return both the actual as well as the fully typed sequence
    return ss, typ


# The following functions are for determination of secondary structure,
# given a list of atoms. The atom format is generic and can be written out
# as PDB or GRO. The coordinates are in Angstrom.
# NOTE: There is the *OLD* DSSP and the *NEW* DSSP, which require
# different calls. The old version uses '--' to indicate reading from stdin
# whereas the new version uses '-i /dev/stdin'
def call_dssp(chain, atomlist, executable='dsspcmbi'):
    '''Get the secondary structure, by calling to dssp'''
    ssdfile = 'chain_%s.ssd' % chain.id

    try:
        if os.system(executable+" -V 2>/dev/null"):
            logging.debug("New version of DSSP; Executing '%s -i /dev/stdin -o %s'" % (executable, ssdfile))
            p = subp.Popen([executable, "-i", "/dev/stdin", "-o", ssdfile], stderr=subp.PIPE, stdout=subp.PIPE, stdin=subp.PIPE)
        else:
            logging.debug("Old version of DSSP; Executing '%s -- %s'" % (executable, ssdfile))
            p = subp.Popen([executable, "--", ssdfile], stderr=subp.PIPE, stdout=subp.PIPE, stdin=subp.PIPE)
    except OSError:
        logging.error("A problem occured calling %s." % executable)
        sys.exit(1)

    for atom in atomlist:
        if atom[0][:2] == 'O1': atom = ('O',)+atom[1:]
        if atom[0][0] != 'H' and atom[0][:2] != 'O2': p.stdin.write(pdbOut(atom).encode('utf-8'))
    p.stdin.write('TER\n'.encode('utf-8'))
    data = p.communicate()
    p.wait()
    main, ss = 0, ''
    for line in open(ssdfile).readlines():
        if main and not line[13] == "!": ss += line[16]
        if line[:15] == '  #  RESIDUE AA': main = 1
    return ss

ssDetermination = {
    "dssp": call_dssp
    }
#########################
## 7 # ELASTIC NETWORK ##  -> @ELN <-
#########################
import math

## ELASTIC NETWORK ##

# Only the decay function is defined here, the network
# itself is set up through the Topology class


# The function to determine the decay scaling factor for the elastic network
# force constant, based on the distance and the parameters provided.
# This function is very versatile and can be fitted to most commonly used
# profiles, including a straight line (rate=0)
def decayFunction(distance, shift, rate, power):
    return math.exp(-rate*math.pow(distance-shift, power))


def rubberBands(atomList, lowerBound, upperBound, decayFactor, decayPower, forceConstant, minimumForce):
    out = []
    u2  = upperBound**2
    while len(atomList) > 3:
        bi, xi = atomList.pop(0)
        for bj, xj in atomList[2:]:
            # Mind the nm/A conversion -- This has to be standardized! Global use of nm?
            d2 = distance2(xi, xj)/100

            if d2 < u2:
                dij  = math.sqrt(d2)
                fscl = decayFunction(dij, lowerBound, decayFactor, decayPower)
                if fscl*forceConstant > minimumForce:
                    out.append({"atoms": (bi, bj), "parameters": (dij, "RUBBER_FC*%f" % fscl)})
    return out
#######################
## 8 # STRUCTURE I/O ##  -> @IO <-
#######################
import logging, math, random, sys

#----+---------+
## A | PDB I/O |
#----+---------+

d2r = 3.14159265358979323846264338327950288/180

# Reformatting of lines in structure file
pdbBoxLine  = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n"


def pdbBoxString(box):
    # Box vectors
    u, v, w  = box[0:3], box[3:6], box[6:9]

    # Box vector lengths
    nu, nv, nw = [math.sqrt(norm2(i)) for i in (u, v, w)]

    # Box vector angles
    alpha = nv*nw == 0 and 90 or math.acos(cos_angle(v, w))/d2r
    beta  = nu*nw == 0 and 90 or math.acos(cos_angle(u, w))/d2r
    gamma = nu*nv == 0 and 90 or math.acos(cos_angle(u, v))/d2r

    return pdbBoxLine % (10*norm(u), 10*norm(v), 10*norm(w), alpha, beta, gamma)


def pdbAtom(a):
    ##01234567890123456789012345678901234567890123456789012345678901234567890123456789
    ##ATOM   2155 HH11 ARG C 203     116.140  48.800   6.280  1.00  0.00
    if a.startswith("TER"):
        return 0
    # NOTE: The 27th field of an ATOM line in the PDB definition can contain an
    #       insertion code. We shift that 20 bits and add it to the residue number
    #       to ensure that residue numbers will be unique.
    # NOTE: Some glycan res names in GLYCAM have four characters, therefore res name
    #       is extended to include the 21st column.
    ## ===> atom name,        res name,         res id,                         chain,
    atom = [a[12:16].strip(), a[17:21].strip(), int(a[22:26])+(ord(a[26])<<20), a[21],
    ##            x,               y,               z
            float(a[30:38]), float(a[38:46]), float(a[46:54])]
    # If the chain identifier is empty, the chain is to None
    if atom[3].strip() == '':
        atom[3] = None
    return tuple(atom) 

def pdbLink(a):
    ##01234567890123456789012345678901234567890123456789012345678901234567890123456789
    ##LINK         ND2 ASN   132                 C1  4YB   211
    
    # NOTE: The 27/57th field of a LINK line in the PDB definition can contain an
    #       insertion code. We shift that 20 bits and add it to the residue number
    #       to ensure that residue numbers will be unique.
    # NOTE: Some glycan res names in GLYCAM have four characters, therefore res name
    #       is extended to include the 21st column.
    
    # Probably a bit of sloppiness in the GLYCAM pdb file: The space for the insertion 
    # code in the 57th field is not blank if not used, but contains the Line Feed char.
    if a[56] == u"\u000A": # = Line Feed
        field57 = ' '
    else:
        field57 = a[56] 
    
    ## ===> atom name 1,      res name 1,       res id 1,                       chain1    
    line = [a[12:16].strip(), a[17:21].strip(), int(a[22:26])+(ord(a[26])<<20), a[21],
    ##      atom name 2,      res name 2,       res id 2,                         chain2
            a[42:46].strip(), a[47:51].strip(), int(a[52:56])+(ord(field57)<<20), a[51]]

    # If the chain identifier is empty, the chainID is set to None
    for i in (3, 7):
        if line[i] == ' ':
            line[i] = None
    
    return tuple(line)
    
     


def pdbOut(atom, i=1, **kwargs):
    # insc contains the insertion code, shifted by 20-bitwise.
    # This means there are multiple residues with the same "resi",
    # which we circumvent by subtracting "insc" from "resi".
    # At other places this subtraction as to be inverted.
    insc = atom[2] >> 20
    resi = atom[2]-(insc << 20)
    if atom[3] == None:
        chain = ' '
    else:
        chain = atom[3]
    pdbline = "ATOM  %5i  %-3s %3s%2s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f           %1s \n"
    if "ssid" in kwargs and type(kwargs["ssid"]) == type(int()):
        occupancy = kwargs["ssid"]
    else:
        occupancy = 40
    return pdbline % ((i, atom[0][:3], atom[1], chain, resi, chr(insc)) + atom[4:] + (1, occupancy, atom[0][0]))


def isPdbAtom(a):
    return a.startswith("ATOM") or (options["-hetatm"] and a.startswith("HETATM")) or a.startswith("TER")


def pdbBoxRead(a):
    fa, fb, fc, aa, ab, ac = [float(i) for i in a.split()[1:7]]
    ca, cb, cg, sg         = math.cos(d2r*aa), math.cos(d2r*ab), math.cos(d2r*ac), math.sin(d2r*ac)
    wx, wy                 = 0.1*fc*cb, 0.1*fc*(ca-cb*cg)/sg
    wz                     = math.sqrt(0.01*fc*fc - wx*wx - wy*wy)
    return [0.1*fa, 0, 0, 0.1*fb*cg, 0.1*fb*sg, 0, wx, wy, wz]


# Function for splitting a PDB file in chains, based
# on chain identifiers and TER statements
def pdbChains(pdbAtomList):
    ''' Iterator whose argument is a list of either
    - 7-tuples (atom name, res name, res id, chain, x, y, z) like 
        ('N', 'GLN', 33554433, None, -64.929, 35.092, 45.292), or 
    - 0 (where has been a TER statement).
    Yields lists of atoms (7-tuples) with the same chain id.
    '''
    chain = []
    for atom in pdbAtomList:
        if not atom:  # Was a "TER" statement
            if chain:
                yield chain
            else:
                logging.info("Skipping empty chain definition.")
            chain = []
            continue
        if not chain or chain[-1][3] == atom[3]: # Fill chain with atoms of same chain id
            chain.append(atom)
        else:
            yield chain # old chain is finished.
            chain = [atom] # Start new chain
    if chain:
        yield chain


# Simple PDB iterator
def pdbFrameIterator(streamIterator):
    title, atoms, box = [], [], []
    for i in streamIterator:
        if i.startswith("ENDMDL"):
            yield ("".join(title), atoms, box)
            title, atoms, box = [], [], []
        elif i.startswith("TITLE"):
            title.append(i)
        elif i.startswith("CRYST1"):
            box = pdbBoxRead(i)
        elif i.startswith("ATOM") or i.startswith("HETATM"):
            atoms.append(pdbAtom(i))
    if atoms:
        yield ("".join(title), atoms, box)


def glycamFrameIterator(streamIterator):
    ''' Basically a copy of the pdbFrameIterator that yields
    in addition the LINK section.
    '''
    title, atoms, box, links = [], [], [], []
    for i in streamIterator:
        if i.startswith("ENDMDL"):
            yield ("".join(title), atoms, box, links)
            title, atoms, box, links = [], [], [], []
        elif i.startswith("TITLE"):
            title.append(i)
        elif i.startswith("CRYST1"):
            box = pdbBoxRead(i)
        elif i.startswith("LINK"): # Extract link section for GLYCAM pdbs
            links.append(pdbLink(i))
        elif i.startswith("ATOM") or i.startswith("HETATM"):
            atoms.append(pdbAtom(i))
    if atoms:
        yield ("".join(title), atoms, box, links)    


#----+---------+
## B | GRO I/O |
#----+---------+

groline = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"


def groBoxRead(a):
    b = [float(i) for i in a.split()] + 6*[0]                    # Padding for rectangular boxes
    return b[0], b[3], b[4], b[5], b[1], b[6], b[7], b[8], b[2]  # Return full definition xx,xy,xz,yx,yy,yz,zx,zy,zz


def groAtom(a):
    # In PDB files, there might by an insertion code. To handle this, we internally add
    # constant to all resids. To be consistent, we have to do the same for gro files.
    # 32 equal ord(' '), eg an empty insertion code
    constant = 32 << 20
    #012345678901234567890123456789012345678901234567890
    #    1PRN      N    1   4.168  11.132   5.291
    #  ===> atom name,        res name,          res id,             chain,
    return (a[10:15].strip(), a[5:10].strip(),   int(a[:5])+constant, None,
    #                x,                 y,                 z
            10*float(a[20:28]), 10*float(a[28:36]), 10*float(a[36:44]))


# Simple GRO iterator
def groFrameIterator(streamIterator):
    while True:
        try:
            title = next(streamIterator)
        except StopIteration:
            break
        natoms = next(streamIterator).strip()
        if not natoms:
            break
        natoms = int(natoms)
        atoms  = [groAtom(next(streamIterator)) for i in range(natoms)]
        box    = groBoxRead(next(streamIterator))
        yield (title, atoms, box)


#----+-------------+
## C | GENERAL I/O |
#----+-------------+

# It is not entirely clear where this fits in best.
# Called from main.
def getChargeType(resname, resid, choices):
    '''Get user input for the charge of residues, based on list with choices.'''
    print('Which %s type do you want for residue %s:' % (resname, resid+1))
    for i, choice in choices.items():
        print('%s. %s' % (i, choice))
    choice = None
    while choice not in list(choices.keys()):
        choice = eval(input('Type a number:'))
    return choices[choice]


# *NOTE*: This should probably be a CheckableStream class that
# reads in lines until either of a set of specified conditions
# is met, then setting the type and from thereon functioning as
# a normal stream.
def streamTag(stream):
    # Tag the stream with the type of structure file
    # If necessary, open the stream, taking care of
    # opening using gzip for gzipped files

    # First check whether we have have an open stream or a file
    # If it's a file, check whether it's zipped and open it
    if type(stream) == str:
        if stream.endswith("gz"):
            logging.info('Read input structure from zipped file.')
            s = gzip.open(stream)
        else:
            logging.info('Read input structure from file.')
            s = open(stream)
    else:
        logging.info('Read input structure from command-line')
        s = stream

    # Read a few lines, but save them
    x = [s.readline(), s.readline()]
    if x[-1].strip().isdigit():
        # If the second line consists of a single integer (ie. the total number of atoms) it must be a GRO file
        logging.info("Input structure is a GRO file. Chains will be labeled consecutively.")
        yield "GRO"
    elif x[0].split()[0] == "LINK":
        # If the file starts with a LINK section assume it is a PDB file from GLYCAM
        logging.info("Input structure is a PDB file from GLYCAM.")
        yield "GLYCAM"
    else:
        # Must be a PDB file then
        # Could wind further to see if we encounter an "ATOM" record
        logging.info("Input structure is a PDB file.")
        yield "PDB"

    # Hand over the lines that were stored
    for i in x:
        yield i

    # Now give the rest of the lines from the stream
    for i in s:
        yield i


#----+-----------------+
## D | STRUCTURE STUFF |
#----+-----------------+

# This list allows to retrieve atoms based on the name or the index
# If standard, dictionary type indexing is used, only exact matches are
# returned. Alternatively, partial matching can be achieved by setting
# a second 'True' argument.
class Residue(list):
    def __getitem__(self, tag):
        if type(tag) == int:
            # Call the parent class __getitem__
            return list.__getitem__(self, tag)
        if type(tag) == str:
            for i in self:
                if i[0] == tag:
                    return i
            else:
                return
        if tag[1]:
            return [i for i in self if tag[0] in i[0]]  # Return partial matches
        else:
            return [i for i in self if i[0] == tag[0]]  # Return exact matches only


def residues(atomList):
    ''' Argument is a list of 7-tuples (atom name, res name, res id, chain id, x, y, z).
    Yields Residue instances of lists of atoms (7-tuples) belonging to the same residue in a chain.    
    '''
    residue = [atomList[0]] # residue is a bit misleading, this is the first atom 7-tuple
    for atom in atomList[1:]:
        if (atom[1] == residue[-1][1] and  # Residue name check
            atom[2] == residue[-1][2] and  # Residue id check
            atom[3] == residue[-1][3]):    # Chain id check
            residue.append(atom)
        else:
            yield Residue(residue)
            residue = [atom]
    yield Residue(residue)


def residueDistance2(r1, r2):
    return min([distance2(i, j) for i in r1 for j in r2])


def breaks(residuelist, selection=("N", "CA", "C"), cutoff=2.5):
    # Extract backbone atoms coordinates
    bb = [[atom[4:] for atom in residue if atom[0] in selection] for residue in residuelist]
    # Needed to remove waters residues from mixed residues.
    bb = [res for res in bb if res != []]

    # We cannot rely on some standard order for the backbone atoms.
    # Therefore breaks are inferred from the minimal distance between
    # backbone atoms from adjacent residues.
    return [i+1 for i in range(len(bb)-1) if residueDistance2(bb[i], bb[i+1]) > cutoff]


def contacts(atoms, cutoff=5):
    rla = list(range(len(atoms)))
    crd = [atom[4:] for atom in atoms]
    return [(i, j) for i in rla[:-1] for j in rla[i+1:]
            if distance2(crd[i], crd[j]) < cutoff]


def add_dummy(beads, dist=0.11, n=2):
    # Generate a random vector in a sphere of -1 to +1, to add to the bead position
    v    = [random.random()*2.-1, random.random()*2.-1, random.random()*2.-1]
    # Calculated the length of the vector and divide by the final distance of the dummy bead
    norm_v = norm(v)/dist
    # Resize the vector
    vn   = [i/norm_v for i in v]
    # m sets the direction of the added vector, currently only works when adding one or two beads.
    m = 1
    for j in range(n):
        newName = 'SCD'
        newBead = (newName, tuple([i+(m*j) for i, j in zip(beads[-1][1], vn)]), beads[-1][2])
        beads.append(newBead)
        m *= -2
    return beads


def check_merge(chains, m_list=[], l_list=[], ss_cutoff=0):
    chainIndex = list(range(len(chains)))

    if 'all' in m_list:
        logging.info("All chains will be merged in a single moleculetype.")
        return chainIndex, [chainIndex]

    chainID = [chain.id for chain in chains]

    # Mark the combinations of chains that need to be merged
    merges = []
    if m_list:
        # Build a dictionary of chain IDs versus index
        # To give higher priority to top chains the lists are reversed
        # before building the dictionary
        chainIndex.reverse()
        chainID.reverse()
        dct = dict(list(zip(chainID, chainIndex)))
        chainIndex.reverse()
        # Convert chains in the merge_list to numeric, if necessary
        # NOTE The internal numbering is zero-based, while the
        # command line chain indexing is one-based. We have to add
        # one to the number in the dictionary to bring it on par with
        # the numbering from the command line, but then from the
        # result we need to subtract one again to make indexing
        # zero-based
        merges = [[(i.isdigit() and int(i) or dct[i]+1)-1 for i in j] for j in m_list]
        for i in merges:
            i.sort()

    # Rearrange merge list to a list of pairs
    pairs = [(i[j], i[k]) for i in merges for j in range(len(i)-1) for k in range(j+1, len(i))]

    # Check each combination of chains for connections based on
    # ss-bridges, links and distance restraints
    for i in chainIndex[:-1]:
        for j in chainIndex[i+1:]:
            if (i, j) in pairs:
                continue
            
            # Check here whether chain j is a glycan that is part of (glyco)protein i
            if isinstance(chains[j],Glycan) and chains[j].anchorres in chains[i]:
                logging.info('Glycan chain {} is linked to {} chain {}.'.format(j+1,chains[i]._type,i+1))
                pairs.append((i, j))
                continue
            #~ else:
                #~ logging.debug('Glycoprotein NOT detected! Pair (%d, %d)' % (i, j))
                #~ logging.debug('Chain '+str(j)+' is a glycan: '+str(isinstance(chains[j],Glycan)))
                #~ logging.debug('Anchor res: '+str(chains[j].anchorres)+' is in chain '+str(i)+': '+str(chains[j].anchorres in chains[i]))
                #~ logging.debug('Chain ID: '+str(chains[i].id))
                
                #print(chains[i].atoms())
                

            # Check whether any link links these two groups
            for a, b in l_list:
                if ((a in chains[i] and b in chains[j]) or (a in chains[j] and b in chains[i])):
                    logging.info("Merging chains %d and %d to allow link %s" % (i+1, j+1, str((a, b))))
                    pairs.append(i < j and (i, j) or (j, i))
                    break
            if (i, j) in pairs:
                continue
            # Check whether any cystine bond given links these two groups
            #for a,b in s_list:
            #    if ((a in chains[i] and b in chains[j]) or
            #        (a in chains[j] and b in chains[i])):
            #        logging.info("Merging chains %d and %d to allow cystine bridge"%(i+1,j+1))
            #        pairs.append( i<j and (i,j) or (j,i) )
            #        break
            #if (i,j) in pairs:
            #    continue
            # Check for cystine bridges based on distance
            if not ss_cutoff:
                continue
            # Get SG atoms from cysteines from either chain
            # Check this pair of chains
            for cysA in chains[i]["CYS"]:
                for cysB in chains[j]["CYS"]:
                    d2 = distance2(cysA["SG"][4:7], cysB["SG"][4:7])
                    if d2 <= ss_cutoff:
                        logging.info("Found SS contact linking chains %d and %d (%f nm)" % (i+1, j+1, math.sqrt(d2)/10))
                        pairs.append((i, j))
                    break
                if (i, j) in pairs:
                    break

    # Sort the combinations
    pairs.sort(reverse=True)

    merges = []
    while pairs:
        merges.append(set([pairs[-1][0]]))
        for i in range(len(pairs)-1, -1, -1):
            if pairs[i][0] in merges[-1]:
                merges[-1].add(pairs.pop(i)[1])
            elif pairs[i][1] in merges[-1]:
                merges[-1].add(pairs.pop(i)[0])
    merges = [list(i) for i in merges]
    for i in merges:
        i.sort()

    order = [j for i in merges for j in i]

    if merges:
        logging.warning("Merging chains.")
        logging.warning("This may change the order of atoms and will change the number of topology files.")
        logging.info("Merges: " + ", ".join([str([j+1 for j in i]) for i in merges]))

    if len(merges) == 1 and len(merges[0]) > 1 and set(merges[0]) == set(chainIndex):
        logging.info("All chains will be merged in a single moleculetype.")

    # Determine the order for writing; merged chains go first
    merges.extend([[j] for j in chainIndex if j not in order])
    order.extend([j for j in chainIndex if j not in order])

    return order, merges

def splitGlycans(glycan_links, reslist):
        '''Splits the LINK section of a GLYCAM PDB file into glycan chains.
        Returns a list of Glycan instances, each containing their Sugars.
        '''        
        glycans = []
        
        logging.debug('LINK section contains %s entries.' % (len(glycan_links)))
        # The LINK entries are sorted by ascending residue id of atom 1. Since glycan 
        # chains are appended to the protein ATOM section, the aa-sugar linkages appear 
        # first in the list.
        
        for line in glycan_links:
            
            newlink = Link(line, reslist)
            newsugar = newlink.sugar
            
            # Try to add Sugar to an existing Glycan...
            for i in glycans:
                #... if the linked-to atom of the new sugar is part of the glycan... 
                if newsugar.link_to in i:
                    i += newlink
                    i += newsugar
                    break            
            # ...otherwise start a new glycan chain.
            else:
                glycans.append(Glycan(newsugar))
                glycans[-1] += newlink
            
        return glycans



class Link:
    def __init__(self, linkline, reslist):
        # ('ND2',    'ASN',     132,   None,     'C1',    '4YB',     211,  None)
        atm1name, res1name, res1num, chain1, atm2name, res2name, res2num, chain2 = linkline        
        
        # Test whether sugar is known
        if res2name in glycamCode:
            self.sugar = Sugar(res2name, res2num, chain2)
        else:
            logging.error('Residue name {} was not found in the glycamCode dict.'.format(self.res2name))
            return
                    
        self.sugar.linkingatom = (atm2name, res2name, res2num, chain2 or "")
        self.sugar.red_end_carbon_atom = self.sugar.linkingatom[0]                      # most often 'C1'; 'C2' for sialic acids, fructose
        self.sugar.link_to     = (atm1name, res1name, res1num, chain1 or "")            # ASN, SER, THR typically, or another sugar, or ROH for a free reducing end
        
        # Dictionary of substituents
        self.sugar.subst = {self.sugar.red_end_carbon_atom: self.sugar.link_to[1:]}

        # Get atom names from Residue instance and the Residue instance itself
        self.sugar.residue = None
        for index, res in enumerate(reslist):
            if res[0][1:3] == self.sugar.linkingatom[1:3] and res[0][3] == chain2:
                self.sugar.atoms = [res[i] for i in res]
                self.sugar.residue = res
                self.reslistIndex = index
        if not self.sugar.residue:
            logging.error('Sugar {} {} does not appear in the residue list!'.format(self.sugar.resname, self.sugar.resnum - (32 << 20)))
            sys.exit(1)

        # Test for reducing end. Reducing end residue should be the one preceding the red end sugar.
        if self.sugar.link_to[1] in ('ROH'):
            self.sugar.redend = reslist[self.reslistIndex-1]

        # Residue A (self.sugar) is linked to residue B
        self.atomnameA = atm2name
        self.resnameA  = res2name
        self.resnumA   = res2num
        self.chainA    = chain2
        self.atidsA    = []

        self.atomnameB = atm1name
        self.resnameB  = res1name
        self.resnumB   = res1num
        self.chainB    = chain1
        self.atidsB    = []

    def __str__(self):
        # return human-readable, GLYCAM style like "DGlcpNAcb1-4DGlcpNAcb" or "DGlcpNAcb1-ND2ASN"
        # If res1 is a sugar use sugar name, atom number, otherwise resname, atomname  
        self.nameB = (self.resnameB in glycamCode) and Sugar(self.resnameB, self.resnumB, self.chainB).longname or self.resnameB
        self.linktoatomnum = (self.resnameB in glycamCode) and self.atomnameB[1:] or self.atomnameB
        return self.sugar.longname + self.sugar.red_end_carbon_atom[1:] + "-" + self.linktoatomnum + self.nameB

class Sugar(object):
    
    def __init__(self, resname, resnum, chain):  
        
        # Residue name, number and chain
        self.resname = resname
        self.resnum  = resnum
        # NOTE: Empty chain ids come as None from pdbLink. Changed here to 
        # '' to be consistent with amino acids
        self.chain   = chain or ''

        # Name (GlcNAc, Man, Gal, Neu5Ac,...)
        self.name = glycamCode[self.resname]
        
        # Anomeric state (alpha or beta) and ring size (pyranose or furanose)
        if self.resname[2] == 'A':
            self.anomer = 'alpha'
            self.ringsize = 'p'
        elif self.resname[2] == 'B':
            self.anomer = 'beta'
            self.ringsize = 'p'
        elif self.resname[2] == 'D':
            self.anomer = 'alpha'
            self.ringsize = 'f'
        elif self.resname[2] == 'U':
            self.anomer = 'beta'
            self.ringsize = 'f'
        elif self.resname[2].isupper():
            self.anomer = 'alpha'
            self.ringsize = '?' # check GLYCAM documentation
        elif self.resname[2].islower():
            self.anomer = 'alpha'
            self.ringsize = '?' # check GLYCAM documentation            
        else:    
            logging.warning('Could not determine the anomeric state of {}'.format(self.resname))

        # Configuration of the asymmetric carbon with highest number (D or L)
        if self.resname[1].islower():
            self.D_L = 'L'
        elif self.resname[1].isupper():
            self.D_L = 'D'
            
        # GLYCAM style name: DGlcpNAcb
        self.longname = ''.join([self.D_L, self.name[:3], self.ringsize, self.name[3:], self.anomer[0]])
        
        # By default sugar is not at reducing end
        self.redend = False                

    def __str__(self):
        return ('%s (%s %d) with substituents: %s' % (self.name, self.resname, self.resnum, self.subst))
        

## !! NOTE !! ##
## XXX The chain class needs to be simplified by extracting things to separate functions/classes
class Chain(object):
    # Attributes defining a chain
    # When copying a chain, or slicing, the attributes in this list have to
    # be handled accordingly.
    _attributes = ("residues", "sequence", "seq", "ss", "ssclass", "sstypes")

    def __init__(self, options, residuelist=[], name=None, multiscale=False):
        
        # residuelist is a list of Residue instances of atom lists (7-tuples (atom name, res name, res id, chain id, x, y, z))
        # belonging to the chain (same chain id or separated by TER)
        # Here Residue.__getitem__ behaves as list.__getitem__

        self.residues   = residuelist
        self._atoms     = [atom[:3] for residue in residuelist for atom in residue]
        self.sequence   = [residue[0][1] for residue in residuelist]
        # *NOTE*: Check for unknown residues and remove them if requested
        #         before proceeding.
        self.seq        = "".join([AA321.get(i, "X") for i in self.sequence])
        self.ss         = ""
        self.ssclass    = ""
        self.sstypes    = ""
        self.mapping    = []
        self.multiscale = multiscale
        self.options    = options

        # Unknown residues
        self.unknowns   = "X" in self.seq

        # Determine number of atoms
        self.natoms     = self._natoms()
        
        # Determine the type of chain
        self._type      = ""
        self.type()

        # BREAKS: List of indices of residues where a new fragment starts
        # Only when polymeric (protein, DNA, RNA, ...)
        # For now, let's remove it for the Nucleic acids...
        self.breaks     = self.type() in ("Protein", "Mixed") and breaks(self.residues) or []

        # LINKS:  List of pairs of pairs of indices of linked residues/atoms
        # This list is used for cysteine bridges and peptide bonds involving side chains
        # The list has items like ((#resi1, #atid1), (#resi2, #atid2))
        # When merging chains, the residue number needs ot be update, but the atom id
        # remains unchanged.
        # For the coarse grained system, it needs to be checked which beads the respective
        # atoms fall in, and bonded terms need to be added there.
        self.links      = []

        # Chain identifier; try to read from residue definition if no name is given
        self.id         = name or residuelist and residuelist[0][0][3] or ""

        # Container for coarse grained beads
        self._cg        = None

    def __len__(self):
        # Return the number of residues
        # DNA/RNA contain non-CAP d/r to indicate type. We remove those first.
        return len(''.join(i for i in self.seq if i.isupper()))

    def __add__(self, other):
        newchain = Chain(name=self.id+"+"+other.id)
        # Combine the chain items that can be simply added
        for attr in self._attributes:
            setattr(newchain, attr, getattr(self, attr) + getattr(other, attr))
        # Set chain items, shifting the residue numbers
        shift  = len(self)
        newchain.breaks     = self.breaks + [shift] + [i+shift for i in other.breaks]
        newchain.links      = self.links + [((i[0]+shift, i[1]), (j[0]+shift, j[1])) for i, j in other.links]
        newchain.natoms     = len(newchain.atoms())
        newchain.multiscale = self.multiscale or other.multiscale
        # Return the merged chain
        return newchain

    def __eq__(self, other):
        return (self.seq        == other.seq    and
                self.ss         == other.ss     and
                self.breaks     == other.breaks and
                self.links      == other.links  and
                self.multiscale == other.multiscale)

    # Extract a residue by number or the list of residues of a given type
    # This facilitates selecting residues for links, like chain["CYS"]
    def __getitem__(self, other):
        if type(other) == str:
            if other not in self.sequence:
                return []
            return [i for i in self.residues if i[0][1] == other]
        elif type(other) == tuple:
            # This functionality is set up for links
            # between coarse grained beads. So these are
            # checked first,
            for i in self.cg():
                if other == i[:4]:
                    return i
            else:
                for i in self.atoms():
                    if other[:3] == i[:3]:
                        return i
                else:
                    return []
        elif type(other) == slice:
            # This implements the old __getslice__ method 
            i, j = other.start, other.stop
            newchain = Chain(self.options, name=self.id)
            # Extract the slices from all lists
            for attr in self._attributes:
                setattr(newchain, attr, getattr(self, attr)[i:j])
            # Breaks that fall within the start and end of this chain need to be passed on.
            # Residue numbering is increased by 20 bits!!
            ch_sta, ch_end      = newchain.residues[0][0][2], newchain.residues[-1][0][2]
            newchain.breaks     = [crack for crack in self.breaks if ch_sta < (crack << 20) < ch_end]
            newchain.links      = [link for link in self.links if ch_sta < (link << 20) < ch_end]
            newchain.multiscale = self.multiscale
            newchain.natoms     = len(newchain.atoms())
            # newchain.type()
            # Return the chain slice
            return newchain
        return self.sequence[other]

    def _contains(self, atomlist, atom):
        atnm, resn, resi, chn = atom

        # If the chain does not match, bail out
        if chn != self.id:
            return False

        # Check if the whole tuple is in
        if atnm and resn and resi:
            return (atnm, resn, resi) in self.atoms()

        # Fetch atoms with matching residue id
        match = (not resi) and atomlist or [j for j in atomlist if j[2] == resi]
        if not match:
            return False

        # Select atoms with matching residue name
        match = (not resn) and match or [j for j in match if j[1] == resn]
        if not match:
            return False

        # Check whether the atom is given and listed
        if not atnm or [j for j in match if j[0] == atnm]:
            return True

        # It just is not in the list!
        return False

    def __contains__(self, other):
        '''other is of format (atnm, resn, resi, chn)'''
        return self._contains(self.atoms(), other) or self._contains(self.cg(), other)

    def __hash__(self):
        return id(self)

    def atoms(self):
        if not self._atoms:
            self._atoms = [atom[:3] for residue in self.residues for atom in residue]
        return self._atoms

    # Split a chain based on residue types; each subchain can have only one type
    def split(self):
        chains = []
        chainStart = 0
        for i in range(len(self.sequence)-1):
            if residueTypes.get(self.sequence[i], "Unknown") != residueTypes.get(self.sequence[i+1], "Unknown"):
                # Use the __getslice__ method to take a part of the chain.
                chains.append(self[chainStart:i+1])
                chainStart = i+1
        if chains:
            logging.debug('Splitting chain %s in %s chains' % (self.id, len(chains)+1))
        return chains + [self[chainStart:]]

    def getname(self, basename=None):
        name = []
        if basename:                      name.append(basename)
        if self.type() and not basename:  name.append(self.type())
        if self.id:
            if type(self.id) == int:
                name.append(chr(64+self.id))
            elif self.id.strip():
                name.append(str(self.id))
        return "_".join(name)

    def set_ss(self, ss, source="self"):
        if len(ss) == 1:
            self.ss = len(self)*ss
        else:
            self.ss = ss
        # Infer the Martini backbone secondary structure types
        self.ssclass, self.sstypes = ssClassification(self.ss, source)

    def dss(self, method=None, executable=None):
        # The method should take a list of atoms and return a
        # string of secondary structure classifications
        if self.type() == "Protein":
            if method:
                atomlist = [atom for residue in self.residues for atom in residue]
                self.set_ss(ssDetermination[method](self, atomlist, executable), source=method)
            else:
                self.set_ss(len(self)*"C")
        else:
            self.set_ss(len(self.sequence)*"-")
        return self.ss

    def type(self, other=None):
        if other:
            self._type = other
        elif not self._type and len(self):
            # Determine the type of chain
            self._type     = set([residueTypes.get(i, "Unknown") for i in set(self.sequence)])
            self._type     = len(self._type) > 1 and "Mixed" or list(self._type)[0]
        return self._type

    def _natoms(self):
        '''Return number of atoms in the chain.'''
        return len(self._atoms)


    # XXX The following (at least the greater part of it) should be made a separate function, put under "MAPPING"
    def cg(self, force=False, com=False):
        # Generate the coarse grained structure
        # Set the b-factor field to something that reflects the secondary structure

        # If the coarse grained structure is set already, just return,
        # unless regeneration is forced.
        if self._cg and not force:
            return self._cg
        self._cg = []
        atid     = 1
        bb       = [1]
        fail     = False
        previous = ''
        for residue, rss, resname in zip(self.residues, self.sstypes, self.sequence):
            # For DNA we need to get the O3' to the following residue when calculating COM
            # The force and com options ensure that this part does not affect itp generation or anything else
            if com:
                # Just an initialization, this should complain if it isn't updated in the loop
                store = 0
                for ind, i in enumerate(residue):
                    if i[0] == "O3'":
                        if previous != '':
                            residue[ind] = previous
                            previous = i
                        else:
                            store = ind
                            previous = i
                # We couldn't remove the O3' from the 5' end residue during the loop so we do it now
                if store > 0:
                    del residue[store]

            # Check if residues names has changed, for example because user has set residues interactively.
            residue = [(atom[0], resname)+atom[2:] for atom in residue]
            if residue[0][1] in ("SOL", "HOH", "TIP"):
                continue
            if not residue[0][1] in list(CoarseGrained.mapping.keys()):
                logging.warning("Skipped unknown residue %s\n" % residue[0][1])
                continue
            # Get the mapping for this residue
            # CG.map returns bead coordinates and mapped atoms
            # This will fail if there are (too many) atoms missing, which is
            # only problematic if a mapped structure is written; the topology
            # is inferred from the sequence. So this is the best place to raise
            # an error
            try:
                beads, ids = map(residue, ca2bb=self.options['ForceField'].ca2bb)
                beads      = list(zip(CoarseGrained.names[residue[0][1]], beads, ids))
                if residue[0][1] in self.options['ForceField'].polar:
                    beads = add_dummy(beads, dist=0.14, n=2)
                elif residue[0][1] in self.options['ForceField'].charged:
                    beads = add_dummy(beads, dist=0.11, n=1)
            except ValueError:
                logging.error("Too many atoms missing from residue %s %d(ch:%s):",
                              residue[0][1], residue[0][2]-(32 << 20), residue[0][3])
                logging.error(repr([i[0] for i in residue]))
                fail = True

            for name, (x, y, z), ids in beads:
                # Add the bead with coordinates and secondary structure id to the list
                self._cg.append((name, residue[0][1][:3], residue[0][2], residue[0][3], x, y, z, ss2num[rss]))
                # Add the ids to the list, after converting them to indices to the list of atoms
                self.mapping.append([atid+i for i in ids])

            # Increment the atom id; This pertains to the atoms that are included in the output.
            atid += len(residue)

            # Keep track of the numbers for CONECTing
            bb.append(bb[-1]+len(beads))

        if fail:
            logging.error("Unable to generate coarse grained structure due to missing atoms.")
            sys.exit(1)

        return self._cg

    def conect(self):
        # Return pairs of numbers that should be CONECTed
        # First extract the backbone IDs
        cg = self.cg()
        bb = [i+1 for i, j in zip(list(range(len(cg))), cg) if j[0] == "BB"]
        bb = list(zip(bb, bb[1:]+[len(bb)]))
        # Set the backbone CONECTs (check whether the distance is consistent with binding)
        conect = [(i, j) for i, j in bb[:-1] if distance2(cg[i-1][4:7], cg[j-1][4:7]) < 14]
        # Now add CONECTs for sidechains
        for i, j in bb:
            nsc = j-i-1
            
            
class Glycan(Chain):
    
    def __init__(self, sugar, options=[]):
        
        # residuelist is a list of Residue instances of atom lists (7-tuples (atom name, res name, res id, chain id, x, y, z))
        # belonging to the chain (same chain id or separated by TER)
        # Here Residue.__getitem__ behaves as list.__getitem__

        # List of Sugar instances
        self.sugars    = [sugar]
       
        # List of Residue instances
        self.residues   = [sugar.residue]
        
        # List of all atoms (atmname, resname, resid) in the glycan chain
        self._atoms     = self.atoms()

        # List of sugar res names ['4YB', '4YB', 'VMB', ...]
        self.sequence   = [sugar.resname]
        
        # String sequence in one-letter code
        self.seq        = self._seq()
        
        # TODO: Do we need a more readable sequence format, eg. that indicates
        # branching?
        
        self.mapping    = []
        
        # MULTISCALING currently not supported
        self.multiscale = False
        
        # This might be useful at some point
        self.options    = options

        # Determine number of atoms
        self.natoms     = self._natoms()
        
        # Set the type of chain, should always return 'Glycan'
        self._type      = ''
        self._type      = self.type()

        # BREAKS and LINKS are currently not supported for glycans
        self.breaks     = breaks(self.residues) or []
        self.links      = []

        # Chain identifier
        self.id         = self.residues and self.residues[0][0][3] or ''

        # Container for coarse grained beads
        self._cg        = None

        # Anchor of this glycan
        self.anchorres = sugar.link_to                # (atmname, resname, resnum, chainid)
        self.anchorresnum = self.anchorres[2]         # necessary to check for circular glycan    

        # Determine glycan name
        if self.anchorres[0:2] == ('ND2', 'ASN'):
            self.name = 'N-glycan ({})'.format(sugar.name)            
        elif self.anchorres[0][0] == 'O' and self.anchorres[1] in ('SER','THR'):
            self.name = 'O-glycan ({})'.format(sugar.name)
        elif self.anchorres[0:2] == ('O1', 'ROH'):
            self.name = 'Free glycan'
        else:
            self.name = 'Special glycan (missing link?)'

        # sugar-sugar linkages
        self.glyclinks = []
    
    def __contains__(self, other):
        '''other is of format (atnm, resn, resi, chn)'''
        return self._contains(self.atoms(), other) #or self._contains(self.cg(), other)


    def __iadd__(self, other):

        if isinstance(other, Sugar):
            sugar = other
            if sugar.linkingatom not in self:
                # add sugar to glycan
                self.sugars.append(sugar)
                self.residues.append(sugar.residue)
                self.sequence.append(sugar.resname)
                self.seq        = self._seq()
                self.natoms     = self._natoms()
                
                # update substituent dict of connecting sugar
                for i in self.sugars:
                    if (i.resname, i.resnum, i.chain) == sugar.link_to[1:]:
                        i.subst[sugar.link_to[0]] = (sugar.resname, sugar.resnum, sugar.chain)
            else:
                logging.debug('Attempted to add a sugar that is already part of glycan.')
        elif isinstance(other, Link):
            self.glyclinks.append(other)
        return self 


    def __str__(self):
        return ('{} anchored to {} with {} branches'.format(self.name, self.anchorres, self._nbranches()))
        
    def __len__(self):
        return len(self.sugars)
        
    def atoms(self):
        return [atom[:3] for residue in self.residues for atom in residue]
        
    def _natoms(self):
        return len(self.atoms())
        
    def _seq(self):
        return ''.join([AA321.get(i, "?") for i in self.sequence])
        
    def _nbranches(self):
        #logging.debug('Branch count over sugars: {}'.format([len(sugar.subst) for sugar in self.sugars]))
        # The number of terminal sugars = the number of branches
        return sum([1 for sugar in self.sugars if len(sugar.subst)==1])
    
    def cg(self, force=False, com=False):
        # Generate the coarse grained structure
        # Set the b-factor field to something that reflects the secondary structure

        # If the coarse grained structure is set already, just return,
        # unless regeneration is forced.
        if self._cg and not force:
            return self._cg
        self._cg = []
        atid     = 1
        fail     = False
        previous = ''
        
        #logging.debug('{} has {} sugars'.format(self.name, len(self)))
        for sugar in self.sugars:
            
            if not sugar.name in list(CoarseGrained.sugarmapping.keys()):
                logging.error("No mapping for sugar %s! Check/Update CoarseGrained.sugarmapping dict." % sugar.name)
                fail = True
                continue
            # Get the mapping for this residue
            # CG.map returns bead coordinates and mapped atoms
            # This will fail if there are (too many) atoms missing, which is
            # only problematic if a mapped structure is written; the topology
            # is inferred from the sequence. So this is the best place to raise
            # an error
            try:
                # beads is a list of center of mass coordinates for the CG beads in the sugar
                # ids is a list of list of atom indices for the CG beads in the sugar
                beads, ids = sugarmap(sugar, sugar.redend)
                
                # Convert beads to [(beadname_sugar, [comx, comy, comz], [atom indices]),  ...]
                beads      = list(zip(CoarseGrained.beadnames_sugar[sugar.name], beads, ids))
                
            except ValueError:
                logging.error("Too many atoms missing from sugar %s %s %d(ch:%s):",
                              sugar.name, sugar.resname, sugar.resnum-(32 << 20), sugar.residue[0][3])
                logging.error(repr([i[0] for i in sugar.residue]))
                fail = True

            for name, (x, y, z), ids in beads:
                # Add the bead with bead name, res name, resnum, chain id, coordinates and secondary structure id to the list
                # ss id is set to 0 (coil) for sugars
                self._cg.append((name, sugar.resname[:3], sugar.resnum, sugar.residue[0][3], x, y, z, 0))
                
                # Add the ids to the list, after converting them to indices to the list of atoms
                self.mapping.append([atid+i for i in ids])

            # Increment the atom id; This pertains to the atoms that are included in the output.
            atid += len(sugar.atoms)

            # Keep track of the numbers for CONECTing
            #bb.append(bb[-1]+len(beads))

        if fail:
            logging.error("Unable to generate coarse grained structure due to missing atoms and/or missing mapping.")
            sys.exit(1)

        return self._cg



##################
## 7 # TOPOLOGY ##  -> @TOP <-
##################
import logging, math


# This is a generic class for Topology Bonded Type definitions
class Bonded(object):
    # The init method is generic to the bonded types,
    # but may call the set method if atoms are given
    # as (ID, ResidueName, SecondaryStructure) tuples
    # The set method is specific to the different types.
    def __init__(self, other=None, options=None, **kwargs):
        self.atoms      = []
        self.type       = -1
        self.parameters = []
        self.comments   = []
        self.category   = None

        if options and type(options) == dict:
            self.options = options
        if other:
            # If other is given, then copy the attributes
            # if it is of the same class or set the
            # attributes according to the key names if
            # it is a dictionary
            if other.__class__ == self.__class__:
                for attr in dir(other):
                    if not attr[0] == "_":
                        setattr(self, attr, getattr(other, attr))
            elif type(other) == dict:
                for attr in list(other.keys()):
                    setattr(self, attr, other[attr])
            elif type(other) in (list, tuple):
                self.atoms = other

        # For every item in the kwargs keys, set the attribute
        # with the same name. This can be used to specify the
        # attributes directly or to override attributes
        # copied from the 'other' argument.
        for key in kwargs:
            setattr(self, key, kwargs[key])

        # If atoms are given as tuples of
        # (ID, ResidueName[, SecondaryStructure])
        # then determine the corresponding parameters
        # from the lists above
        if self.atoms and type(self.atoms[0]) == tuple:
            self.set(self.atoms, **kwargs)

    def __bool__(self):
        return bool(self.atoms)

    def __str__(self):
        if not self.atoms or not self.parameters:
            return ""
        s = ["%5d" % i for i in self.atoms]
        # For exclusions, no type is defined, which equals -1
        if self.type != -1: s.append(" %5d " % self.type)
        # Print integers and floats in proper format and neglect None terms
        s.extend([formatString(i) for i in self.parameters if i is not None])
        if self.comments:
            s.append(';')
            if type(self.comments) == str:
                s.append(self.comments)
            else:
                s.extend([str(i) for i in self.comments])
        return " ".join(s)

    def __iadd__(self, num):
        self.atoms = [i + int(num) for i in self.atoms]
        return self

    def __add__(self, num):
        out  = self.__class__(self)
        out += num
        return out

    def __eq__(self, other):
        if type(other) in (list, tuple):
            return self.atoms == other
        else:
            return self.atoms == other.atoms and self.type == other.type and self.parameters == other.parameters

    # This function needs to be overridden for descendents
    def set(self, atoms, **kwargs):
        pass


# The set method of this class will look up parameters for backbone beads
# Side chain bonds ought to be set directly, using the constructor
# providing atom numbers, bond type, and parameters
# Constraints are bonds with kb = None, which can be extracted
# using the category
class Bond(Bonded):
    def set(self, atoms, **kwargs):
        ids, r, ss, ca  = list(zip(*atoms))
        self.atoms      = ids
        self.type       = 1
        self.positionCa = ca
        self.comments   = "%s(%s)-%s(%s)" % (r[0], ss[0], r[1], ss[1])
        # The category can be used to keep bonds sorted
        self.category   = kwargs.get("category")

        self.parameters = self.options['ForceField'].bbGetBond(r, ca, ss)
        # Backbone bonds also can be constraints.
        # We could change the type further on, but this is more general.
        # Even better would be to add a new type: BB-Constraint
        if self.parameters[1] == None:
            self.category = 'Constraint'

    # Overriding __str__ method to suppress printing of bonds with Fc of 0
    def __str__(self):
        if len(self.parameters) > 1 and self.parameters[1] == 0:
            return ""
        return Bonded.__str__(self)


# Similar to the preceding class
class Angle(Bonded):
    def set(self, atoms, **kwargs):
        ids, r, ss, ca  = list(zip(*atoms))
        self.atoms      = ids
        self.type       = 2
        self.positionCa = ca
        self.comments   = "%s(%s)-%s(%s)-%s(%s)" % (r[0], ss[0], r[1], ss[1], r[2], ss[2])
        self.category   = kwargs.get("category")
        self.parameters = self.options['ForceField'].bbGetAngle(r, ca, ss)


# Similar to the preceding class
class Vsite(Bonded):
    def set(self, atoms, **kwargs):
        ids, r, ss, ca  = list(zip(*atoms))
        self.atoms      = ids
        self.type       = 1
        self.positionCa = ca
        self.comments   = "%s" % (r[0])
        self.category   = kwargs.get("category")
        self.parameters = kwargs.get("parameters")


# Similar to the preceding class
class Exclusion(Bonded):
    def set(self, atoms, **kwargs):
        ids, r, ss, ca  = list(zip(*atoms))
        self.atoms      = ids
        self.positionCa = ca
        self.comments   = "%s" % (r[0])
        self.category   = kwargs.get("category")
        self.parameters = kwargs.get("parameters")


# Similar to the preceding class
class Dihedral(Bonded):
    def set(self, atoms, **kwargs):
        ids, r, ss, ca  = list(zip(*atoms))
        self.atoms      = ids
        self.type       = 1
        self.positionCa = ca
        self.comments   = "%s(%s)-%s(%s)-%s(%s)-%s(%s)" % (r[0], ss[0], r[1], ss[1], r[2], ss[2], r[3], ss[3])
        self.category   = kwargs.get("category")

        if ''.join(i for i in ss) == 'FFFF':
            # Collagen
            self.parameters = self.options['ForceField'].bbDihedDictD['F']
        elif ''.join(i for i in ss) == 'EEEE' and self.options['ExtendedDihedrals']:
            # Use dihedrals
            self.parameters = self.options['ForceField'].bbDihedDictD['E']
        elif set(ss).issubset("H123"):
            # Helix
            self.parameters = self.options['ForceField'].bbDihedDictD['H']
        else:
            self.parameters = None


# This list allows to retrieve Bonded class items based on the category
# If standard, dictionary type indexing is used, only exact matches are
# returned. Alternatively, partial matching can be achieved by setting
# a second 'True' argument.
class CategorizedList(list):
    def __getitem__(self, tag):
        if type(tag) == int:
            # Call the parent class __getitem__
            return list.__getitem__(self, tag)

        if type(tag) == str:
            return [i for i in self if i.category == tag]

        if tag[1]:
            return [i for i in self if tag[0] in i.category]
        else:
            return [i for i in self if i.category == tag[0]]

class CategorizedAtomlist(list):
    
    def __init__(self):
        # Container for original residue numbers
        self.origResnum = []
    
    def __getitem__(self, tag):
        if type(tag) == int:
            # Call the parent class __getitem__
            return list.__getitem__(self, tag)

        if type(tag) == str:
            return [i for i in self if i.category == tag]

        #                              resname            original resnum
        if type(tag) == tuple and type(tag[0]) == str and type(tag[1]) == int:
            return [i for i, j in zip(self, self.origResnum) if (i[3], j) == tag]
  
        if tag[1]:
            return [i for i in self if tag[0] in i.category]
        else:
            return [i for i in self if i.category == tag[0]]


class Topology(object):
    def __init__(self, other=None, options=None, name=""):
        self.name        = name or ''
        self.nrexcl      = 1
        self.atoms       = CategorizedAtomlist()
        self.vsites      = CategorizedList()
        self.exclusions  = CategorizedList()
        self.bonds       = CategorizedList()
        self.angles      = CategorizedList()
        self.dihedrals   = CategorizedList()
        self.impropers   = CategorizedList()
        self.constraints = CategorizedList()
        self.posres      = CategorizedList()
        self.sequence    = []
        self.secstruc    = ""
        # Okay, this is sort of funny; we will add a
        #   #define mapping virtual_sitesn
        # to the topology file, followed by a header
        #   [ mapping ]
        self.mapping     = []
        # For multiscaling we have to keep track of the number of
        # real atoms that correspond to the beads in the topology
        self.natoms      = 0
        self.multiscale  = options['multi']

        self.options     = options or {}
        self.interchainlinks = []
        if not other:
            # Returning an empty instance
            return
        elif isinstance(other, Topology):
            for attrib in ["atoms", "vsites", "bonds", "angles", "dihedrals", "impropers", "constraints", "posres", "exclusions"]:
                setattr(self, attrib, getattr(other, attrib, []))
        elif isinstance(other, Chain):
            self.type      = other.type()
            if self.type   == "Protein":
                self.fromAminoAcidSequence(other)
            elif self.type == "Nucleic":
                # Currently there are no Martini Nucleic Acids
                self.fromNucleicAcidSequence(other)
            elif self.type == "Glycan":
                self.fromGlycanSequence(other)            
            elif sef.type  == "Mixed":
                logging.warning('Mixed type chains are not yet implemented!')
                # How can you have a mixed chain?
                # Well, you could get a covalently bound lipid or piece of DNA to a protein :S
                # But how to deal with that?
                # Probably one should separate the chains into blocks of specified type,
                # determine the locations of links, then construct the topologies for the
                # blocks and combine them according to the links.
                pass
            else:
                # This chain should not be polymeric, but a collection of molecules
                # For each unique residue type fetch the proper moleculetype
                self.fromMoleculeList(other)

    def __iadd__(self, other):
        if not isinstance(other, Topology):
            other = Topology(other)
        shift     = len(self.atoms)
        last      = self.atoms[-1]
        # The following used to work: zip>list expansions>zip back, but that only works if
        # all the tuples in the original list are of equal length. With masses and charges
        # that is not necessarily the case.
        if shift != len(self.atoms.origResnum):
            logging.error('Original residue number list should contain {} instead of {} entries.'.format(shift, len(self.atoms.origResnum)))
            logging.error('This will lead to wrong interchainlink bonds.')        
        for atom, origResnum in zip(other.atoms, other.atoms.origResnum):
            atom = list(atom)
            atom[0] += shift    # Update atom numbers
            atom[2] += last[2]  # Update residue numbers
            atom[5] += last[5]  # Update charge group numbers
            self.atoms.append(tuple(atom))
            self.atoms.origResnum.append(origResnum)
        for attrib in ["bonds", "vsites", "angles", "dihedrals", "impropers", "constraints", "posres", "exclusions"]:
            getattr(self, attrib).extend([source + shift for source in getattr(other, attrib)])
        # Add interchain links from other topology
        self.interchainlinks.extend(other.interchainlinks)
        return self

    def __add__(self, other):
        out = Topology(self)
        if not isinstance(other, Topology):
            other = Topology(other)
        out += other
        return out

    def __str__(self):
        if self.multiscale:
            out = ['; MARTINI (%s) Multiscale virtual sites topology section for "%s"' % (self.options['ForceField'].name, self.name)]
        else:
            string = '; MARTINI (%s) Coarse Grained topology file for "%s"' % (self.options['ForceField'].name, self.name)
            string += '\n; Created by martinize_glycam.py version %s \n; Using the following options:  ' % (self.options['Version'])
            string += ' '.join(self.options['Arguments'])
            out = [string]
        if self.type == 'Glycan':
            # prevent error for glycans where AA321 keys are not residue names but sugar names
            out += [
                '; Sequence:',
                '; ' + ''.join([AA321.get(AA) for AA in [glycamCode[res] for res in self.sequence]]),
                ]            
            seq = '-'.join(self.sequence)
        elif self.sequence: 
            out += [
                '; Sequence:',
                '; ' + ''.join([AA321.get(AA) for AA in self.sequence]),
                '; Secondary Structure:',
                '; ' + self.secstruc,
                ]

        # Do not print a molecule name when multiscaling
        # In that case, the topology created here needs to be appended
        # at the end of an atomistic moleculetype
        if not self.multiscale:
            out += ['\n[ moleculetype ]',
                    '; Name         Exclusions',
                    '%-15s %3d' % (self.name, self.nrexcl)]

        out.append('\n[ atoms ]')

        # For virtual sites and dummy beads we have to be able to specify the mass.
        # Thus we need two different format strings:
        fs8 = '%5d %5s %5d %5s %5s %5d %7.4f ; %s'
        fs9 = '%5d %5s %5d %5s %5s %5d %7.4f %7.4f ; %s'
        out.extend([len(i) == 9 and fs9 % i or fs8 % i for i in self.atoms])

        # Print out the vsites only if they exist. Right now it can only be type 1 virual sites.
        vsites = [str(i) for i in self.vsites]
        if vsites:
            out.append('\n[ virtual_sites2 ]')
            out.extend(vsites)

        # Print out the exclusions only if they exist.
        exclusions = [str(i) for i in self.exclusions]
        if exclusions:
            out.append('\n[ exclusions ]')
            out.extend(exclusions)

        if self.multiscale:
            out += ['\n;\n; Coarse grained to atomistic mapping\n;',
                    '#define mapping virtual_sitesn',
                    '[ mapping ]']
            for i, j in self.mapping:
                out.append(("%5d     2 " % i)+" ".join(["%5d" % k for k in j]))

            logging.info('Created virtual sites section for multiscaled topology')
            return "\n".join(out)

        # Bonds in order: backbone, backbone-sidechain, sidechain, short elastic, long elastic
        out.append("\n[ bonds ]")
        # Backbone-backbone
        bonds = [str(i) for i in self.bonds["BB"]]
        if bonds:
            out.append("; Backbone bonds")
            out.extend(bonds)
        # Rubber Bands
        bonds = [str(i) for i in self.bonds["Rubber", True]]
        if bonds:
            # Add a CPP style directive to allow control over the elastic network
            out.append("#ifndef NO_RUBBER_BANDS")
            out.append("#ifndef RUBBER_FC\n#define RUBBER_FC %f\n#endif" % self.options['ElasticMaximumForce'])
            out.extend(bonds)
            out.append("#endif")
        # Backbone-Sidechain/Sidechain-Sidechain
        bonds = [str(i) for i in self.bonds["SC"]]
        if bonds:
            out.append("; Sidechain bonds")
            out.extend(bonds)
        # Short elastic/Long elastic
        bonds = [str(i) for i in self.bonds["Elastic short"]]
        if bonds:
            out.append("; Short elastic bonds for extended regions")
            out.extend(bonds)
        bonds = [str(i) for i in self.bonds["Elastic long"]]
        if bonds:
            out.append("; Long elastic bonds for extended regions")
            out.extend(bonds)
        # Cystine bridges
        bonds = [str(i) for i in self.bonds["Cystine"]]
        if bonds:
            out.append("; Cystine bridges")
            out.extend(bonds)
        # Glycans    
        bonds = [str(i) for i in self.bonds["Glycan"]]
        if bonds:
            out.append("; Glycans")
            out.extend(bonds)        
        # Other links
        bonds = [str(i) for i in self.bonds["Link"]]
        if bonds:
            out.append("; Links")
            out.extend(bonds)

        # Constraints
        out.append("\n[ constraints ]")
        out.extend([str(i) for i in self.bonds["Constraint"]])

        # Angles
        out.append("\n[ angles ]")
        out.append("; Backbone angles")
        out.extend([str(i) for i in self.angles["BBB"]])
        out.append("; Backbone-sidechain angles")
        out.extend([str(i) for i in self.angles["BBS"]])
        out.append("; Sidechain angles")
        out.extend([str(i) for i in self.angles["SC"]])
        # Glycans    
        angles = [str(i) for i in self.angles["Glycan"]]
        if angles:
            out.append("; Glycans")
            out.extend(angles)        


        # Dihedrals
        out.append("\n[ dihedrals ]")
        out.append("; Backbone dihedrals")
        out.extend([str(i) for i in self.dihedrals["BBBB"] if i.parameters])
        out.append("; Sidechain improper dihedrals")
        out.extend([str(i) for i in self.dihedrals["SC"] if i.parameters])
        # Glycans    
        dihedrals = [str(i) for i in self.dihedrals["Glycan"]]
        if dihedrals:
            out.append("; Glycans")
            out.extend(dihedrals)        


        # Postition Restraints
        if self.posres:
            out.append("\n#ifdef POSRES")
            out.append("#ifndef POSRES_FC\n#define POSRES_FC %.2f\n#endif" % self.options['PosResForce'])
            out.append(" [ position_restraints ]")
            out.extend(['  %5d    1    POSRES_FC    POSRES_FC    POSRES_FC' % i for i in self.posres])
            out.append("#endif")

        logging.info('Created coarse-grained topology.')
        return "\n".join(out)

    def fromAminoAcidSequence(self, sequence, secstruc=None, links=None,
                              breaks=None, mapping=None, rubber=False,
                              multi=False):
        '''The sequence function can be used to generate the topology for
           a sequence :) either given as sequence or as chain'''

        # Shift for the atom numbers of the atomistic part in a chain
        # that is being multiscaled
        shift = 0
        # First check if we get a sequence or a Chain instance
        if isinstance(sequence, Chain):
            chain         = sequence
            # Collect original residue numbers
            resnums       = [res[0][2] for res in chain.residues]
            links         = chain.links
            breaks        = chain.breaks
            # If the mapping is not specified, the actual mapping is taken,
            # used to construct the coarse grained system from the atomistic one.
            # The function argument "mapping" could be used to use a default
            # mapping scheme in stead, like the mapping for the GROMOS96 force field.
            mapping       = mapping or chain.mapping
            multi         = self.options['multi'] or chain.multiscale
            self.secstruc = chain.sstypes or len(chain)*"C"
            self.sequence = chain.sequence
            # If anything hints towards multiscaling, do multiscaling
            self.multiscale = self.multiscale or chain.multiscale or multi
            if self.multiscale:
                shift        = self.natoms
                self.natoms += len(chain.atoms())
        elif not secstruc:
            # If no secondary structure is provided, set all to coil
            chain         = None
            self.secstruc = len(self.sequence)*"C"
        else:
            # If a secondary structure is provided, use that. chain is none.
            chain         = None
            self.secstruc = secstruc

        logging.debug('Secondary structure ({},...), length {}'.format(self.secstruc[:10], len(self.secstruc)))
        logging.debug('Sequence ({},...), length {}'.format(self.sequence[:10], len(self.sequence)))

        # Fetch the sidechains
        # Pad with empty lists for atoms, bonds, angles
        # and dihedrals, and take the first four lists out
        # This will avoid errors for residues for which
        # these are not defined.

        sc = [(self.options['ForceField'].sidechains[res]+5*[[]])[:5] for res in self.sequence]

        # ID of the first atom/residue
        # The atom number and residue number follow from the last
        # atom c.q. residue id in the list processed in the topology
        # thus far. In the case of multiscaling, the real atoms need
        # also be accounted for.
        startAtom = self.natoms + 1
        startResi = self.atoms and self.atoms[-1][2]+1 or 1

        # Backbone bead atom IDs
        bbid = [startAtom]
        for i in next(zip(*sc)):
            bbid.append(bbid[-1]+len(i)+1)

        # Calpha positions, to get Elnedyn BBB-angles and BB-bond lengths
        # positionCa = [residue[1][4:] for residue in chain.residues]
        # The old method (line above) assumed no hydrogens: Ca would always be
        # the second atom of the residue. Now we look at the name.
        positionCa = []
        for residue in chain.residues:
            for atom in residue:
                if atom[0] == "CA":
                    positionCa.append(atom[4:])

        # Residue numbers for this moleculetype topology
        resid = list(range(startResi, startResi+len(self.sequence)))

        # This contains the information for deriving backbone bead types,
        # bb bond types, bbb/bbs angle types, and bbbb dihedral types and
        # Elnedyn BB-bondlength BBB-angles
        seqss = list(zip(bbid, self.sequence, self.secstruc, positionCa))

        # Fetch the proper backbone beads
        bb = [self.options['ForceField'].bbGetBead(res, typ) for num, res, typ, Ca in seqss]

        # If termini need to be charged, change the bead types
        if not self.options['NeutralTermini']:
            bb[0]  = "Qd"
            bb[-1] = "Qa"

        # If breaks need to be charged, change the bead types
        if self.options['ChargesAtBreaks']:
            for i in breaks:
                bb[i]   = "Qd"
                bb[i-1] = "Qa"

        # For backbone parameters, iterate over fragments, inferred from breaks
        for i, j in zip([0]+breaks, breaks+[-1]):
            # Extract the fragment
            frg = j == -1 and seqss[i:] or seqss[i:j]

            # Iterate over backbone bonds
            self.bonds.extend([Bond(pair, category="BB", options=self.options,) for pair in zip(frg, frg[1:])])

            # Iterate over backbone angles
            # Don't skip the first and last residue in the fragment
            self.angles.extend([Angle(triple, options=self.options, category="BBB") for triple in zip(frg, frg[1:], frg[2:])])

            # Get backbone quadruples
            quadruples = list(zip(frg, frg[1:], frg[2:], frg[3:]))

            # No i-1,i,i+1,i+2 interactions defined for Elnedyn
            if self.options['ForceField'].UseBBBBDihedrals:
                # Process dihedrals
                for q in quadruples:
                    id, rn, ss, ca = list(zip(*q))
                    # Maybe do local elastic networks
                    if ss == ("E", "E", "E", "E") and not self.options['ExtendedDihedrals']:
                        # This one may already be listed as the 2-4 bond of a previous one
                        if not (id[0], id[2]) in self.bonds:
                            self.bonds.append(Bond(
                                options    = self.options,
                                atoms      = (id[0], id[2]),
                                parameters = self.options['ForceField'].ebonds['short'],
                                type       = 1,
                                comments   = "%s(%s)-%s(%s) 1-3" % (rn[0], id[0], rn[2], id[2]),
                                category   = "Elastic short"))
                        self.bonds.append(Bond(
                                options    = self.options,
                                atoms      = (id[1], id[3]),
                                parameters = self.options['ForceField'].ebonds['short'],
                                type       = 1,
                                comments   = "%s(%s)-%s(%s) 2-4" % (rn[1], id[1], rn[3], id[3]),
                                category   = "Elastic short"))
                        self.bonds.append(Bond(
                                options    = self.options,
                                atoms      = (id[0], id[3]),
                                parameters = self.options['ForceField'].ebonds['long'],
                                type       = 1,
                                comments   = "%s(%s)-%s(%s) 1-4" % (rn[0], id[0], rn[3], id[3]),
                                category   = "Elastic long"))
                    else:
                        # Since dihedrals can return None, we first collect them separately and then
                        # add the non-None ones to the list
                        dihed = Dihedral(q, options=self.options, category="BBBB")
                        if dihed:
                            self.dihedrals.append(dihed)

            # Elnedyn does not use backbone-backbone-sidechain-angles
            if self.options['ForceField'].UseBBSAngles:
                # Backbone-Backbone-Sidechain angles
                # If the first residue has a sidechain, we take SBB, otherwise we skip it
                # For other sidechains, we 'just' take BBS
                if len(frg) > 1 and frg[1][0]-frg[0][0] > 1:
                    self.angles.append(Angle(
                        options    = self.options,
                        atoms      = (frg[0][0] + 1, frg[0][0], frg[1][0]),
                        parameters = self.options['ForceField'].bbsangle,
                        type       = 2,
                        comments   = "%s(%s)-%s(%s) SBB" % (frg[0][1], frg[0][2], frg[1][1], frg[1][2]),
                        category   = "BBS"))

                # Start from first residue: connects sidechain of second residue
                for (ai, ni, si, ci), (aj, nj, sj, cj), s in zip(frg[0:], frg[1:], sc[1:]):
                    if s[0]:
                        self.angles.append(Angle(
                            options    = self.options,
                            atoms      = (ai, aj, aj+1),
                            parameters = self.options['ForceField'].bbsangle,
                            type       = 2,
                            comments   = "%s(%s)-%s(%s) SBB" % (ni, si, nj, sj),
                            category   = "BBS"))

        # Now do the atom list, and take the sidechains along
        #
        # AtomID AtomType ResidueID ResidueName AtomName ChargeGroup Charge ; Comments
        atid = startAtom
        for resi, origResnum, resname, bbb, sidechn, ss in zip(resid, resnums, self.sequence, bb, sc, self.secstruc):
            scatoms, bon_par, ang_par, dih_par, vsite_par = sidechn

            # Side chain bonded terms
            # Collect bond, angle and dihedral connectivity
            bon_con, ang_con, dih_con, vsite_con = (self.options['ForceField'].connectivity[resname]+4*[[]])[:4]

            # Side Chain Bonds/Constraints
            for atids, par in zip(bon_con, bon_par):
                if par[1] == None:
                    self.bonds.append(Bond(
                        options    = self.options,
                        atoms      = atids,
                        parameters = [par[0]],
                        type       = 1,
                        comments   = resname,
                        category   = "Constraint"))
                else:
                    self.bonds.append(Bond(
                        options    = self.options,
                        atoms      = atids,
                        parameters = par,
                        type       = 1,
                        comments   = resname,
                        category   = "SC"))
                # Shift the atom numbers
                self.bonds[-1] += atid

            # Side Chain Angles
            for atids, par in zip(ang_con, ang_par):
                self.angles.append(Angle(
                    options    = self.options,
                    atoms      = atids,
                    parameters = par,
                    type       = 2,
                    comments   = resname,
                    category   = "SC"))
                # Shift the atom numbers
                self.angles[-1] += atid

            # Side Chain Dihedrals
            for atids, par in zip(dih_con, dih_par):
                self.dihedrals.append(Dihedral(
                    options    = self.options,
                    atoms      = atids,
                    parameters = par,
                    type       = 2,
                    comments   = resname,
                    category   = "SC"))
                # Shift the atom numbers
                self.dihedrals[-1] += atid

            # Side Chain V-Sites
            for atids, par in zip(vsite_con, vsite_par):
                self.vsites.append(Vsite(
                    options    = self.options,
                    atoms      = atids,
                    parameters = par,
                    type       = 1,
                    comments   = resname,
                    category   = "SC"))
                # Shift the atom numbers
                self.vsites[-1] += atid

            # Side Chain exclusions
            # The polarizable forcefield give problems with the charges in the sidechain,
            # if the backbone is also charged.
            # To avoid that, we add explicit exclusions
            if bbb in list(self.options['ForceField'].charges.keys()) and resname in list(self.options['ForceField'].mass_charge.keys()):
                for i in [j for j, d in enumerate(scatoms) if d == 'D']:
                    self.exclusions.append(Exclusion(
                        options    = self.options,
                        atoms      = (atid, i+atid+1),
                        comments   = '%s(%s)' % (resname, resi),
                        parameters = (None, )))

            # All residue atoms
            counter = 0  # Counts over beads
            for atype, aname in zip([bbb] + list(scatoms), CoarseGrained.residue_bead_names):
                if self.multiscale:
                    atype, aname = "v" + atype, "v" + aname
                # If mass or charge diverse, we adopt it here.
                # We don't want to do this for BB beads because of charged termini.
                if resname in list(self.options['ForceField'].mass_charge.keys()) and counter != 0:
                    M, Q = self.options['ForceField'].mass_charge[resname]
                    aname = Q[counter-1] > 0 and 'SCP' or Q[counter-1] < 0 and 'SCN' or aname
                    self.atoms.append((atid, atype, resi, resname, aname, atid,
                                       Q[counter-1], M[counter-1], ss))
                else:
                    self.atoms.append((atid, atype, resi, resname, aname, atid,
                                       self.options['ForceField'].charges.get(atype, 0), ss))
                
                # To allow later mapping of resi to original resnum
                self.atoms.origResnum.append(origResnum)

                # Doing this here save going over all the atoms onesmore.
                # Generate position restraints for all atoms or Backbone beads only.
                if 'all' in self.options['PosRes']:
                    self.posres.append((atid))
                elif aname in self.options['PosRes']:
                    self.posres.append((atid))
                if mapping:
                    self.mapping.append((atid, [i + shift for i in mapping[counter]]))
                atid    += 1
                counter += 1

        # The rubber bands are best applied outside of the chain class, as that gives
        # more control when chains need to be merged. The possibility to do it on the
        # chain level is retained to allow building a complete chain topology in
        # a straightforward manner after importing this script as module.
        if rubber and chain:
            rubberList = rubberBands(
                [(i[0], j[4:7]) for i, j in zip(self.atoms, chain.cg()) if i[4] in ElasticBeads],
                ElasticLowerBound, ElasticUpperBound,
                ElasticDecayFactor, ElasticDecayPower,
                ElasticMaximumForce, ElasticMinimumForce)
            self.bonds.extend([Bond(i, options=self.options, type=6,
                                    category="Rubber band") for i in rubberList])

        # Note the equivalent of atomistic atoms that have been processed
        if chain and self.multiscale:
            self.natoms += len(chain.atoms())

    def fromNucleicAcidSequence(self, sequence, secstruc=None, links=None, breaks=None,
                                mapping=None, rubber=False, multi=False):

        # Shift for the atom numbers of the atomistic part in a chain
        # that is being multiscaled
        shift = 0
        # First check if we get a sequence or a Chain instance
        if isinstance(sequence, Chain):
            chain         = sequence
            links         = chain.links
            breaks        = chain.breaks
            # If the mapping is not specified, the actual mapping is taken,
            # used to construct the coarse grained system from the atomistic one.
            # The function argument "mapping" could be used to use a default
            # mapping scheme in stead, like the mapping for the GROMOS96 force field.
            mapping       = mapping or chain.mapping
            multi         = self.options['multi'] or chain.multiscale
            self.secstruc = chain.sstypes or len(chain)*"C"
            self.sequence = chain.sequence
            # If anything hints towards multiscaling, do multiscaling
            self.multiscale = self.multiscale or chain.multiscale or multi
            if self.multiscale:
                shift        = self.natoms
                self.natoms += len(chain.atoms())
        elif not secstruc:
            # If no secondary structure is provided, set all to coil
            chain         = None
            self.secstruc = len(self.sequence)*"C"
        else:
            # If a secondary structure is provided, use that. chain is none.
            chain         = None
            self.secstruc = secstruc

        logging.debug(self.secstruc)
        logging.debug(self.sequence)

        # Fetch the base information
        # Pad with empty lists for atoms, bonds, angles
        # and dihedrals, and take the first five lists out
        # This will avoid errors for residues for which
        # these are not defined.
        sc = [(self.options['ForceField'].bases[res]+6*[[]])[:6] for res in self.sequence]

        # ID of the first atom/residue
        # The atom number and residue number follow from the last
        # atom c.q. residue id in the list processed in the topology
        # thus far. In the case of multiscaling, the real atoms need
        # also be accounted for.
        startAtom = self.natoms + 1
        startResi = self.atoms and self.atoms[-1][2]+1 or 1

        # Backbone bead atom IDs
        bbid = [[startAtom, startAtom+1, startAtom+2]]
        for i in next(zip(*sc)):
            bbid1 = bbid[-1][0]+len(i)+3
            bbid.append([bbid1, bbid1+1, bbid1+2])

        # Residue numbers for this moleculetype topology
        resid = list(range(startResi, startResi+len(self.sequence)))

        # This contains the information for deriving backbone bead types,
        # bb bond types, bbb/bbs angle types, and bbbb dihedral types.
        seqss = list(zip(bbid, self.sequence, self.secstruc))

        # Fetch the proper backbone beads
        # Since there are three beads we need to split these to the list
        bb = [self.options['ForceField'].bbGetBead(res, typ) for num, res, typ in seqss]
        bb3 = [i for j in bb for i in j]

        # This is going to be usefull for the type of the last backbone bead.
        # If termini need to be charged, change the bead types
        # if not self.options['NeutralTermini']:
        #    bb[0]  ="Qd"
        #    bb[-1] = "Qa"

        # If breaks need to be charged, change the bead types
        # if self.options['ChargesAtBreaks']:
        #    for i in breaks:
        #        bb[i]   = "Qd"
        #        bb[i-1] = "Qa"

        # For backbone parameters, iterate over fragments, inferred from breaks
        for i, j in zip([0]+breaks, breaks+[-1]):
            # Extract the fragment
            frg = j == -1 and seqss[i:] or seqss[i:j]
            # Expand the 3 bb beads per residue into one long list
            # Resulting list contains three tuples per residue
            # We use the useless ca parameter to get the correct backbone bond from bbGetBond
            frg = [(j[0][i], j[1], j[2], i) for j in frg for i in range(len(j[0]))]

            # Iterate over backbone bonds
            self.bonds.extend([Bond(pair, category="BB", options=self.options,) for pair in zip(frg, frg[1:])])

            # Iterate over backbone angles
            # Don't skip the first and last residue in the fragment
            self.angles.extend([Angle(triple, options=self.options, category="BBB") for triple in zip(frg, frg[1:], frg[2:])])

            # Get backbone quadruples
            quadruples = list(zip(frg, frg[1:], frg[2:], frg[3:]))

            # No i-1,i,i+1,i+2 interactions defined for Elnedyn
            # Process dihedrals
            for q in quadruples:
                id, rn, ss, ca = list(zip(*q))
                # Since dihedrals can return None, we first collect them separately and then
                # add the non-None ones to the list
                dihed = Dihedral(q, options=self.options, category="BBBB")
                if dihed:
                    self.dihedrals.append(dihed)

        # Now do the atom list, and take the sidechains along
        #
        atid  = startAtom
        # We need to do some trickery to get all 3 bb beads in to these lists
        # This adds each element to a list three times, feel free to shorten up
        resid3 = [i for i in resid for j in range(3)]
        sequence3 = [i for i in self.sequence for j in range(3)]
        sc3 = [i for i in sc for j in range(3)]
        secstruc3 = [i for i in self.secstruc for j in range(3)]
        count = 0
        for resi, resname, bbb, sidechn, ss in zip(resid3, sequence3, bb3, sc3, secstruc3):
            # We only want one side chain per three backbone beads so this skips the others
            if (count % 3) == 0:
                # Note added impropers in contrast to aa
                scatoms, bon_par, ang_par, dih_par, imp_par, vsite_par = sidechn

                # Side chain bonded terms
                # Collect bond, angle and dihedral connectivity
                # Impropers needed to be added here for DNA
                bon_con, ang_con, dih_con, imp_con, vsite_con = (self.options['ForceField'].connectivity[resname]+5*[[]])[:5]

                # Side Chain Bonds/Constraints
                for atids, par in zip(bon_con, bon_par):
                    if par[1] == None:
                        self.bonds.append(Bond(
                            options    = self.options,
                            atoms      = atids,
                            parameters = [par[0]],
                            type       = 1,
                            comments   = resname,
                            category   = "Constraint"))
                    else:
                        self.bonds.append(Bond(
                            options    = self.options,
                            atoms      = atids,
                            parameters = par,
                            type       = 1,
                            comments   = resname,
                            category   = "SC"))
                    # Shift the atom numbers
                    self.bonds[-1] += atid

                # Side Chain Angles
                for atids, par in zip(ang_con, ang_par):
                    self.angles.append(Angle(
                        options    = self.options,
                        atoms      = atids,
                        parameters = par,
                        type       = 2,
                        comments   = resname,
                        category   = "SC"))
                    # Shift the atom numbers
                    self.angles[-1] += atid

                # Side Chain Dihedrals
                for atids, par in zip(dih_con, dih_par):
                    self.dihedrals.append(Dihedral(
                        options    = self.options,
                        atoms      = atids,
                        parameters = par,
                        type       = 1,
                        comments   = resname,
                        category   = "BSC"))
                    # Shift the atom numbers
                    self.dihedrals[-1] += atid

                # Side Chain Impropers
                for atids, par in zip(imp_con, imp_par):
                    self.dihedrals.append(Dihedral(
                        options    = self.options,
                        atoms      = atids,
                        parameters = par,
                        type       = 2,
                        comments   = resname,
                        category   = "SC"))
                    # Shift the atom numbers
                    self.dihedrals[-1] += atid

                # Side Chain V-Sites
                for atids, par in zip(vsite_con, vsite_par):
                    self.vsites.append(Vsite(
                        options    = self.options,
                        atoms      = atids,
                        parameters = par,
                        type       = 1,
                        comments   = resname,
                        category   = "SC"))
                    # Shift the atom numbers
                    self.vsites[-1] += atid

                # Currently DNA needs exclusions for the base
                # The loop runs over the first backbone bead so 3 needs to be added to the indices
                for i in range(len(scatoms)):
                    for j in range(i+1, len(scatoms)):
                        self.exclusions.append(Exclusion(
                            options    = self.options,
                            atoms      = (i+atid+3, j+atid+3),
                            comments   = '%s(%s)' % (resname, resi),
                            parameters = (None, )))

                # All residue atoms
                counter = 0  # Counts over beads
                # Need to tweak this to get all the backbone beads to the list with the side chain
                bbbset = [bb3[count], bb3[count+1], bb3[count+2]]
                for atype, aname in zip(bbbset+list(scatoms), CoarseGrained.residue_bead_names_dna):
                    if self.multiscale:
                        atype, aname = "v"+atype, "v"+aname
                    self.atoms.append((atid, atype, resi, resname, aname, atid,
                                       self.options['ForceField'].charges.get(atype, 0), ss))
                    # Doing this here saves going over all the atoms onesmore.
                    # Generate position restraints for all atoms or Backbone beads only.
                    if 'all' in self.options['PosRes']:
                        self.posres.append((atid))
                    elif aname in self.options['PosRes']:
                        self.posres.append((atid))
                    if mapping:
                        self.mapping.append((atid, [i+shift for i in mapping[counter]]))
                    atid    += 1
                    counter += 1
            count += 1

        # One more thing, we need to remove dihedrals (2) and an angle (1)  that reach beyond the 3' end
        # This is stupid to do now but the total number of atoms seems not to be available before
        # This iterate the list in reverse order so that removals don't affect later checks
        for i in range(len(self.dihedrals)-1, -1, -1):
            if (max(self.dihedrals[i].atoms) > self.atoms[-1][0]):
                del self.dihedrals[i]
        for i in range(len(self.angles)-1, -1, -1):
            if (max(self.angles[i].atoms) > self.atoms[-1][0]):
                del self.angles[i]

    def fromGlycanSequence(self, sequence, secstruc=None, links=None,
                              breaks=None, mapping=None, rubber=False,
                              multi=False):
        '''The sequence function can be used to generate the topology for
           a sequence :) either given as sequence or as chain'''

        # Shift for the atom numbers of the atomistic part in a chain
        # that is being multiscaled
        shift = 0
        # First check if we get a sequence or a Chain instance
        if isinstance(sequence, Glycan):
            chain         = sequence
            links         = chain.links
            breaks        = chain.breaks
            # If the mapping is not specified, the actual mapping is taken,
            # used to construct the coarse grained system from the atomistic one.
            # The function argument "mapping" could be used to use a default
            # mapping scheme instead, like the mapping for the GROMOS96 force field.
            mapping       = mapping or chain.mapping
            multi         = self.options['multi'] or chain.multiscale
            self.sequence = chain.sequence
            self.seq      = chain.seq
            self.sugars   = chain.sugars
            self.glyclinks= chain.glyclinks

            # If anything hints towards multiscaling, do multiscaling
            self.multiscale = self.multiscale or chain.multiscale or multi
            if self.multiscale:
                shift        = self.natoms
                self.natoms += len(chain.atoms())
        elif not secstruc:
            # If no secondary structure is provided, set all to coil
            chain         = None
            self.secstruc = len(self.sequence)*"C"
        else:
            # If a secondary structure is provided, use that. chain is none.
            chain         = None
            self.secstruc = secstruc

        #logging.debug('Sequence ({},...), length {}'.format(self.sequence[:10], len(self.sequence)))

        # Fetch the sidechains
        # Pad with empty lists for atoms, bonds, angles
        # and dihedrals, and take the first four lists out
        # This will avoid errors for residues for which
        # these are not defined.
        logging.debug('Generating glycan topology...')
        
        # List of list of lists of CG bead type, bond, angle, dihedral tuples            
        # "4YB":  [spl("TNd SP1 SP1 TNa"),     [(0.276, None),(0.307, None),(0.334,18000),(0.300,15000)],[  ( 90,200) , (120,200)],],
        sugars = []
        for res, sugar in zip(self.sequence, self.sugars):
            if not sugar.redend:
                sugars.append((self.options['ForceField'].sugarinternals[res]+5*[[]])[:5])
            # Reducing end sugars have slightly different patterns
            else:
                sugars.append((self.options['ForceField'].redends[res]+5*[[]])[:5])
                
        
        # ID of the first atom/residue
        # The atom number and residue number follow from the last
        # atom c.q. residue id in the list processed in the topology
        # thus far. In the case of multiscaling, the real atoms need
        # also be accounted for.
        startAtom = self.natoms + 1
        startResi = self.atoms and self.atoms[-1][2]+1 or 1

        # Residue numbers for this moleculetype topology
        resids = list(range(startResi, startResi+len(self.sequence)))

        # Now do the atom list
        #
        # AtomID AtomType ResidueID ResidueName AtomName ChargeGroup Charge ; Comments
        atid = startAtom
        for resi, sugar, sugar_par, glyclink in zip(resids, self.sugars, sugars, self.glyclinks):
            cgBeadTypes, bon_par, ang_par, dih_par, vsite_par = sugar_par
                        
            # Internal bonded terms 
            # Collect bond, angle, dihedral, vsite and exclusion connectivity
            bon_con, ang_con, dih_con, vsite_con, excl_con = (self.options['ForceField'].sugarinternalcon[sugar.resname[1]+sugar.ringsize]+5*[[]])[:5]

            # Bonds/Constraints
            for atids, par in zip(bon_con, bon_par):
                if par[1] == None:
                    self.bonds.append(Bond(
                        options    = self.options,
                        atoms      = atids,
                        parameters = [par[0]],          # Only the bond length
                        type       = 1,
                        comments   = '{} {} {}'.format(sugar.resname, sugar.resnum - (32 << 20), str(atids)),
                        category   = "Constraint"))
                else:
                    self.bonds.append(Bond(
                        options    = self.options,
                        atoms      = atids,
                        parameters = par,
                        type       = 1,
                        comments   = '{} {} {}'.format(sugar.resname, sugar.resnum - (32 << 20), str(atids)),
                        category   = "Glycan"))
                # Shift the atom numbers
                self.bonds[-1] += atid

            # Angles
            for atids, par in zip(ang_con, ang_par):
                self.angles.append(Angle(
                    options    = self.options,
                    atoms      = atids,
                    parameters = par,
                    type       = 2,
                    comments   = '{} {} {}'.format(sugar.resname, sugar.resnum - (32 << 20), str(atids)),
                    category   = "Glycan"))
                # Shift the atom numbers
                self.angles[-1] += atid

            # Dihedrals
            for atids, par in zip(dih_con, dih_par):
                self.dihedrals.append(Dihedral(
                    options    = self.options,
                    atoms      = atids,
                    parameters = par,
                    type       = 2,
                    comments   = '{} {} {}'.format(sugar.resname, sugar.resnum - (32 << 20), str(atids)),
                    category   = "Glycan"))
                # Shift the atom numbers
                self.dihedrals[-1] += atid

            # V-Sites
            for atids, par in zip(vsite_con, vsite_par):
                self.vsites.append(Vsite(
                    options    = self.options,
                    atoms      = atids,
                    parameters = par,
                    type       = 1,
                    comments   = '{} {} {}'.format(sugar.resname, sugar.resnum - (32 << 20), str(atids)),
                    category   = "Glycan"))
                # Shift the atom numbers
                self.vsites[-1] += atid

            # Exclusions
            for atids in excl_con:
                self.exclusions.append(Exclusion(
                    options    = self.options,
                    atoms      = atids,
                    comments   = '{} {} {}'.format(sugar.resname, sugar.resnum - (32 << 20), str(atids)),
                    parameters = (None, )))
            
            # Assignment of the atids to the glycan links works only if the entries in self.sugars and self.glyclink
            # match. So we test that and exit if they don't.
            
            if (glyclink.resnameA, glyclink.resnumA) == (sugar.resname, sugar.resnum):
                
                # Add CG beads in this residue to the atom list
                counter = 0  # Counts over beads
                for atype, aname in zip(list(cgBeadTypes), CoarseGrained.beadnames_sugar[sugar.name]):
                    if self.multiscale:
                        atype, aname = "v" + atype, "v" + aname
                    
                    self.atoms.append((atid, atype, resi, sugar.resname, aname, atid,
                                           self.options['ForceField'].charges.get(atype, 0), ''))
                    
                    # The residue names from the source PDB file (sugar.resname) will NOT be kept but be replaced 
                    # by consecutively numbered residue numbers (resi). Add the original residue numbers here
                    # to the Topology.atoms list to be able to map resi and sugar.resname for later use (interchain link
                    # atid assignment)
                    self.atoms.origResnum.append(sugar.resnum)
                    
                    # Add the atids of sugar A to the Link objects
                    glyclink.atidsA.append(atid)
                    # Doing this here saves going over all the atoms once more.
                    # Generate position restraints for all atoms or Backbone beads only.
                    if 'all' in self.options['PosRes']:
                        self.posres.append((atid))
                    elif aname in self.options['PosRes']:
                        self.posres.append((atid))
                    if mapping:
                        
                        self.mapping.append((atid, [i + shift for i in mapping[counter]]))
                    atid    += 1
                    counter += 1
            else:
                logging.error('The glycan link list and the sugar list are not matching. Exiting.')
                sys.exit()
    
        # Now that all the atom numbers (atid) of the chain are assigned, try to fill in the atids B for the links. This should work in all
        # cases except for the first sugar if it links outside the chain (an interchain link), eg. to an amino acid
        for gl in self.glyclinks:
                
                l = [link.atidsA for link in self.glyclinks if (gl.resnameB, gl.resnumB) == (link.resnameA, link.resnumA)]
                if l:
                    gl.atidsB = l[0] 
                    try:
                        
                        bon_con, ang_con, dih_con, vsite_con, excl_con = (self.options['ForceField'].sugarlinkcon[str(gl)]+5*[[]])[:5]
                        bon_par, ang_par, dih_par, vsite_par           = (self.options['ForceField'].sugarlinkpar[str(gl)]+4*[[]])[:4]
                    except KeyError:
                        logging.warning('Parameters missing for link {}. Update sugarlinkcon/-par dictionaries.'.format(gl))
                        continue
                        
                    for con, par in zip(bon_con, bon_par):
                        atidsA = [gl.atidsA[con[0]]]
                        atidsB = [gl.atidsB[con[1]]]
                        #logging.debug('{}, atidsA: {}, atidsB: {}'.format(gl, gl.atidsA, gl.atidsB))

                        if par[1] == None:
                            self.bonds.append(Bond(
                                options    = self.options,
                                atoms      = tuple(atidsA + atidsB),
                                parameters = [par[0]],          # Only the bond length
                                type       = 1,
                                comments   = '{}, {} {}'.format(str(gl), con[0], con[1]),
                                category   = "Constraint"))
                        else:
                            self.bonds.append(Bond(
                                options    = self.options,
                                atoms      = tuple(atidsA + atidsB),
                                parameters = par,
                                type       = 1,
                                comments   = '{}, {} {}'.format(str(gl), con[0], con[1]),
                                category   = "Glycan"))

                        # Angles
                        for con, par in zip(ang_con, ang_par):
                            atidsA = [gl.atidsA[i] for i in con[0]]
                            atidsB = [gl.atidsB[i] for i in con[1]]
                            self.angles.append(Angle(
                                options    = self.options,
                                atoms      = tuple(atidsA + atidsB),
                                parameters = par,
                                type       = 2,
                                comments   = '{}, {} {}'.format(str(gl), con[0], con[1]),
                                category   = "Glycan"))
                        
                        # Dihedrals
                        for con, par in zip(dih_con, dih_par):
                            atidsA = [gl.atidsA[i] for i in con[0]]
                            atidsB = [gl.atidsB[i] for i in con[1]]
                            self.dihedrals.append(Dihedral(
                                options    = self.options,
                                atoms      = tuple(atidsA + atidsB),
                                parameters = par,
                                type       = 2,
                                comments   = '{}, {} {}'.format(str(gl), con[0], con[1]),
                                category   = "Glycan"))
                        
                        # V-Sites
                        for con, par in zip(vsite_con, vsite_par):
                            atidsA = [gl.atidsA[i] for i in con[0]]
                            atidsB = [gl.atidsB[i] for i in con[1]]

                            self.vsites.append(Vsite(
                                options    = self.options,
                                atoms      = tuple(atidsA + atidsB),
                                parameters = par,
                                type       = 1,
                                comments   = '{}, {} {}'.format(str(gl), con[0], con[1]),
                                category   = "Glycan"))
                        
                        # Exclusions
                        for con in excl_con:
                            atidsA = [gl.atidsA[con[0]]]
                            atidsB = [gl.atidsB[con[1]]]
                            self.exclusions.append(Exclusion(
                                options    = self.options,
                                atoms      = tuple(atidsA + atidsB),
                                comments   = '{}, {} {}'.format(str(gl), con[0], con[1]),
                                parameters = (None, )))
                        
                elif gl.resnameB not in CoarseGrained.redends:
                    # add Link to interchainlinks for later retrieval of atids B
                    logging.debug('Found interchain link {}.'.format(str(gl)))
                    self.interchainlinks.append(gl)
                    
                
                

        # The rubber bands are best applied outside of the chain class, as that gives
        # more control when chains need to be merged. The possibility to do it on the
        # chain level is retained to allow building a complete chain topology in
        # a straightforward manner after importing this script as module.
        if rubber and chain:
            rubberList = rubberBands(
                [(i[0], j[4:7]) for i, j in zip(self.atoms, chain.cg()) if i[4] in ElasticBeads],
                ElasticLowerBound, ElasticUpperBound,
                ElasticDecayFactor, ElasticDecayPower,
                ElasticMaximumForce, ElasticMinimumForce)
            self.bonds.extend([Bond(i, options=self.options, type=6,
                                    category="Rubber band") for i in rubberList])

        # Note the equivalent of atomistic atoms that have been processed
        if chain and self.multiscale:
            self.natoms += len(chain.atoms())


    def fromMoleculeList(self, other):
        pass
        
    def getlinkcon(atidsA, atidsB, con):
        # con is a 2-tuple of tuples containing the indices for the atidsA and B respectively
        indicesA, indicesB = con
        return tuple([atidsA[i] for i in indicesA] + [atidsB[i] for i in indicesB])
#############
## 8 # MAIN #  -> @MAIN <-
#############
import sys, logging, random, math, os, re


def main(options):
    # Check whether to read from a gro/pdb file or from stdin
    # We use an iterator to wrap around the stream to allow
    # inferring the file type, without consuming lines already
    inStream = streamTag(options["-f"] and options["-f"].value or sys.stdin)

    # The streamTag iterator first yields the file type, which
    # is used to specify the function for reading frames
    fileType = next(inStream)
    if fileType == "GRO":
        frameIterator = groFrameIterator
    elif fileType == "GLYCAM":
        frameIterator = glycamFrameIterator
    else:
        frameIterator = pdbFrameIterator

    # ITERATE OVER FRAMES IN STRUCTURE FILE #

    # Now iterate over the frames in the stream
    # This should become a StructureFile class with a nice .next method
    model     = 1
    cgOutPDB  = None
    ssTotal   = []
    cysteines = []
    for args in frameIterator(inStream):

        if fileType == "PDB":
            title, atoms, box = args
            
            # The PDB file can have multiple chains (different chain id), TER statements are also interpreted as chain separators.
            # The pdbChains iterator yields the individual chains as lists of atoms (7-tuples).
            # The residues iterator yields Residue instances of lists of atoms (7-tuples) belonging to the same residue in a chain.             
            # A chain may have breaks in which case the breaking residues are flagged
            chains = [Chain(options, [i for i in residues(chain)]) for chain in pdbChains(atoms)]
                        
        elif fileType == "GLYCAM":
            title, atoms, box, glycam_links = args
            
            # do the chain thing as for conventional PDB file
            chains = [Chain(options, [i for i in residues(chain)]) for chain in pdbChains(atoms)]
            
            # The last chain in chains is the one containing the sugars
            sugarreslist = chains[-1].residues
            
        else:
            title, atoms, box = args
            # The GRO file does not define chains. Here breaks in the backbone are
            # interpreted as chain separators.  
            residuelist = [residue for residue in residues(atoms)]
            # The breaks are indices to residues
            broken = breaks(residuelist)
            # Reorder, such that each chain is specified with (i,j,k)
            # where i and j are the start and end of the chain, and
            # k is a chain identifier
            chains = list(zip([0]+broken, broken+[len(residuelist)], list(range(len(broken)+1))))
            chains = [Chain(options, residuelist[i:j], name=chr(65+k)) for i, j, k in chains]

        for chain in chains:
            chain.multiscale = "all" in options['multi'] or chain.id in options['multi']

        # Check the chain identifiers
        if model == 1 and len(chains) != len(set([i.id for i in chains])):
            # Ending down here means that non-consecutive blocks of atoms in the
            # PDB file have the same chain ID. The warning pertains to PDB files only,
            # since chains from GRO files get a unique chain identifier assigned.
            logging.warning("Several chains have identical chain identifiers in the PDB file.")

        # Check if chains are of mixed type. If so, split them.
        # Note that in some cases HETATM residues are part of a
        # chain. This will get problematic. But we cannot cover
        # all, probably.
        if not options['MixedChains']:
            demixedChains = []
            for chain in chains:
                demixedChains.extend(chain.split())
            chains = demixedChains
        
        
        # In case of a GLYCAM file, identify and split glycans,
        # based on LINK section. 
        # Extend chains with glycans.
        if fileType == "GLYCAM":
            
            glycans = splitGlycans(glycam_links, sugarreslist)
            chains.extend(glycans)
            
            # Verbose info on identified glycans
            for i in glycans:
                logging.debug('Found {} with residues:'.format(i))
                for n,j in enumerate(i.sugars):
                    logging.debug('    %2d: %6s (%s %d) linked with:' % (n+1, j.longname, j.resname, j.resnum-(32<<20)))                        
                    for k in j.subst:
                        logging.debug('          %s to %s %d' % (k, j.subst[k][0], j.subst[k][1]-(32<<20)))
                            

        # Check all chains
        keep = []
        for chain in chains:
            if chain.type() == "Water":
                logging.info("Removing %d water molecules (chain %s)." % (len(chain), chain.id))
            elif chain.type() in ("Protein", "Nucleic") or isinstance(chain, Glycan):
                keep.append(chain)
            # This is currently not active:
            elif options['RetainHETATM']:
                keep.append(chain)
            else:
                logging.info("Removing HETATM chain %s consisting of %d residues." % (chain.id, len(chain)))
        chains = keep
        
        logging.info("Keeping %d chains:" % len(chains))
        for n, chain in enumerate(chains):
            logging.info("  %2d:   %s (%s), %d atoms in %d residues." % (n+1, chain.id, chain._type, chain.natoms, len(chain)))
            
        
        # Here we interactively check the charge state of residues
        # Can be easily expanded to residues other than HIS
        for chain in chains:
            if not isinstance(chain, Glycan):
                for i, resname in enumerate(chain.sequence):
                    if resname == 'HIS' and options['chHIS']:
                        choices = {0: 'HIH', 1: 'HIS'}
                        choice = getChargeType(resname, i, choices)
                        chain.sequence[i] = choice
                    
                    # epsilon-protonated histidine in GLYCAM
                    elif resname == 'HIE':
                        chain.sequence[i] = 'HIH'   # charged histidine
                        logging.info('Replaced epsilon-charged histidine HIE {} from GLYCAM with HIH.'.format(chain.residues[i][0][2]-(32<<20)))
                        

        # Check which chains need merging
        if model == 1:
            # order is a list of chain indices
            # merge is a list of list of chain indices to be merged
            order, merge = check_merge(chains, options['mergeList'], options['linkList'], options['CystineCheckBonds'] and options['CystineMaxDist2'])
            
        # Get the total length of the sequence
        seqlength = sum([len(chain) for chain in chains])
        logging.info('Total size of the system: %s residues.' % seqlength)

        ## SECONDARY STRUCTURE
        ss = ''
        if options['Collagen']:
            for chain in chains:
                if not isinstance(chain, Glycan):
                   chain.set_ss("F")
                   ss += chain.ss
        elif options["-ss"]:
            # XXX We need error-catching here,
            # in case the file doesn't exist, or the string contains bogus.
            # If the string given for the sequence consists strictly of upper case letters
            # and does not appear to be a file, assume it is the secondary structure
            ss = options["-ss"].value.replace('~', 'L').replace(' ', 'L')
            if ss.isalnum() and ss.isupper() and not os.path.exists(options["-ss"].value):
                ss = options["-ss"].value
                logging.info('Secondary structure read from command line:\n'+ss)
            else:
                # There ought to be a file with the name specified
                ssfile = [i.strip() for i in open(options["-ss"].value)]

                # Try to read the file as a Gromacs Secondary Structure Dump
                # Those have an integer as first line
                if ssfile[0].isdigit():
                    logging.info('Will read secondary structure from file (assuming Gromacs ssdump).')
                    ss = "".join([i for i in ssfile[1:]])
                else:
                    # Get the secondary structure type from DSSP output
                    logging.info('Will read secondary structure from file (assuming DSSP output).')
                    pss = re.compile(r"^([ 0-9]{4}[0-9]){2}")
                    ss  = "".join([i[16] for i in open(options["-ss"].value) if re.match(pss, i)])

            # Now set the secondary structure for each of the chains
            sstmp = ss
            for chain in chains:
                if not isinstance(chain, Glycan):
                    ln = min(len(sstmp), len(chain))
                    chain.set_ss(sstmp[:ln])
                    sstmp = ss[:ln]
        else:
            if options["-dssp"]:
                method, executable = "dssp", options["-dssp"].value
            #elif options["-pymol"]:
            #    method, executable = "pymol", options["-pymol"].value
            else:
                logging.warning("No secondary structure or determination method speficied. Protein chains will be set to 'COIL'.")
                method, executable = None, None

            for chain in chains:
                if not isinstance(chain, Glycan):
                    ss += chain.dss(method, executable)

            # Used to be: if method in ("dssp","pymol"): but pymol is not supported
            if method in ["dssp"]:
                logging.debug('%s determined secondary structure:\n' % method.upper()+ss)

        # Collect the secondary structure classifications for different frames
        ssTotal.append(ss)
        
        # Check whether glycans belong to proteins
        
        

        # Write the coarse grained structure if requested
        if options["-x"].value:
            logging.info("Writing coarse grained structure file '{}'.".format(options["-x"].value))
            if cgOutPDB is None:
                cgOutPDB = open(options["-x"].value, "w")
            cgOutPDB.write("MODEL %8d\n" % model)
            cgOutPDB.write(title)
            cgOutPDB.write(pdbBoxString(box))
            atid = 1
            for i in order:
                ci = chains[i]
                if ci.multiscale:
                    for r in ci.residues:
                        for name, resn, resi, chain, x, y, z in r:
                            cgOutPDB.write(pdbOut((name, resn[:3], resi, chain, x, y, z),i=atid))
                            atid += 1
                coarseGrained = ci.cg(com=True)

                if coarseGrained:
                    for name, resn, resi, chain, x, y, z, ssid in coarseGrained:
                        if ci.multiscale:
                            name = "v"+name
                        cgOutPDB.write(pdbOut((name, resn[:3], resi, chain, x, y, z),i=atid,ssid=ssid))
                        atid += 1
                    cgOutPDB.write("TER\n")
                else:
                    logging.warning("No mapping for coarse graining chain %s (%s); chain is skipped." % (ci.id, ci.type()))
            cgOutPDB.write("ENDMDL\n")

        # Gather cysteine sulphur coordinates
        cyslist = [cys["SG"] for chain in chains for cys in chain["CYS"]]
        cysteines.append([cys for cys in cyslist if cys])

        model += 1

    # Write the index file if requested.
    # Mainly of interest for multiscaling.
    # Could be improved by adding separate groups for BB, SC, etc.
    if options["-n"].value:
        logging.info("Writing index file.")
        # Lists for All-atom, Virtual sites and Coarse Grain.
        NAA, NVZ, NCG = [], [], []
        atid = 1
        for i in order:
            ci = chains[i]
            coarseGrained = ci.cg(force=True)
            if ci.multiscale:
                NAA.extend([" %5d" % (a+atid) for a in range(ci.natoms)])
                atid += ci.natoms
            if coarseGrained:
                if ci.multiscale:
                    NVZ.extend([" %5d" % (a+atid) for a in range(len(coarseGrained))])
                else:
                    NCG.extend([" %5d" % (a+atid) for a in range(len(coarseGrained))])
                atid += len(coarseGrained)
        outNDX = open(options["-n"].value, "w")
        outNDX.write("\n[ AA ]\n"+"\n".join([" ".join(NAA[i:i+15]) for i in range(0, len(NAA), 15)]))
        outNDX.write("\n[ VZ ]\n"+"\n".join([" ".join(NVZ[i:i+15]) for i in range(0, len(NVZ), 15)]))
        outNDX.write("\n[ CG ]\n"+"\n".join([" ".join(NCG[i:i+15]) for i in range(0, len(NCG), 15)]))
        outNDX.close()

    # Write the index file for mapping AA trajectory if requested
    if options["-nmap"].value:
        logging.info("Writing trajectory index file.")
        atid = 1
        outNDX = open(options["-nmap"].value, "w")
        # Get all AA atoms as lists of atoms in residues
        # First we skip hetatoms and unknowns then iterate over beads
        # In DNA the O3' atom is mapped together with atoms from the next residue
        # This stores it until we get to the next residue
        o3_shift = ''
        for i_count, i in enumerate(residues(atoms)):
            if i[0][1] in ("SOL", "HOH", "TIP"):
                continue
            if not i[0][1] in list(CoarseGrained.mapping.keys()):
                continue
            nra = 0
            names = [j[0] for j in i]
            # This gives out a list of atoms in residue, each tuple has other
            # stuff in it that's needed elsewhere so we just take the last
            # element which is the atom index (in that residue)
            for j_count, j in enumerate(mapIndex(i)):
                outNDX.write('[ Bead %i of residue %i ]\n' % (j_count+1, i_count+1))
                line = ''
                for k in j:
                    if names[k[2]] == "O3'":
                        line += '%s ' % (str(o3_shift))
                        o3_shift = k[2]+atid
                    else:
                        line += '%i ' % (k[2]+atid)
                line += '\n'
                nra += len(j)
                outNDX.write(line)
            atid += nra

    # Evertything below here we only need, if we need to write a Topology
    if options['-o']:

        # Collect the secondary structure stuff and decide what to do with it
        # First rearrange by the residue
        ssTotal = list(zip(*ssTotal))
        ssAver  = []
        for i in ssTotal:
            si = list(set(i))
            if len(si) == 1:
                # Only one type -- consensus
                ssAver.append(si[0])
            else:
                # Transitions between secondary structure types
                i = list(i)
                si = [(1.0*i.count(j)/len(i), j) for j in si]
                si.sort()
                if si[-1][0] > options["-ssc"].value:
                    ssAver.append(si[-1][1])
                else:
                    ssAver.append(" ")

        ssAver = "".join(ssAver)
        logging.info('(Average) Secondary structure has been determined (see head of .itp-file).')

        # Divide the secondary structure according to the division in chains
        # This will set the secondary structure types to be used for the
        # topology.
        for chain in chains:
            chain.set_ss(ssAver[:len(chain)])
            ssAver = ssAver[len(chain):]

        # Now the chains are complete, each consisting of a residuelist,
        # and a secondary structure designation if the chain is of type 'Protein'.
        # There may be mixed chains, there may be HETATM things.
        # Water has been discarded. Maybe this has to be changed at some point.
        # The order in the coarse grained files matches the order in the set of chains.
        #
        # If there are no merges to be done, i.e. no global Elnedyn network, no
        # disulphide bridges, no links, no distance restraints and no explicit merges,
        # then we can write out the topology, which will match the coarse grained file.
        #
        # If there are merges to be done, the order of things may be changed, in which
        # case the coarse grained structure will not match with the topology...

        # CYSTINE BRIDGES #
        # Extract the cysteine coordinates (for all frames) and the cysteine identifiers
        if options['CystineCheckBonds']:
            logging.info("Checking for cystine bridges, based on sulphur (SG) atoms lying closer than %.4f nm" % math.sqrt(options['CystineMaxDist2']/100))

            cyscoord  = list(zip(*[[j[4:7] for j in i] for i in cysteines]))
            cysteines = [i[:4] for i in cysteines[0]]

            bl, kb    = options['ForceField'].special[(("SC1", "CYS"), ("SC1", "CYS"))]

            # Check the distances and add the cysteines to the link list if the
            # SG atoms have a distance smaller than the cutoff.
            rlc = list(range(len(cysteines)))
            for i in rlc[:-1]:
                for j in rlc[i+1:]:
                    # Checking the minimum distance over all frames
                    # But we could also take the maximum, or the mean
                    d2 = min([distance2(a, b) for a, b in zip(cyscoord[i], cyscoord[j])])
                    if d2 <= options['CystineMaxDist2']:
                        a, b = cysteines[i], cysteines[j]
                        options['linkListCG'].append((("SC1", "CYS", a[2], a[3]), ("SC1", "CYS", b[2], b[3]), bl, kb))
                        a, b = (a[0], a[1], a[2]-(32 << 20), a[3]), (b[0], b[1], b[2]-(32 << 20), b[3])
                        logging.info("Detected SS bridge between %s and %s (%f nm)" % (a, b, math.sqrt(d2)/10))

        # REAL ITP STUFF #
        # Check whether we have identical chains, in which case we
        # only write the ITP for one...
        # This means making a distinction between chains and
        # moleculetypes.
        # merge is a list of lists of chain indices to be merged
        molecules = [tuple([chains[i] for i in j]) for j in merge]

        # At this point we should have a list or dictionary of chains
        # Each chain should be given a unique name, based on the value
        # of options["-o"] combined with the chain identifier and possibly
        # a number if there are chains with identical identifiers.
        # For each chain we then write an ITP file using the name for
        # moleculetype and name + ".itp" for the topology include file.
        # In addition we write a master topology file, using the value of
        # options["-o"], with an added extension ".top" if not given.

        # XXX *NOTE*: This should probably be gathered in a 'Universe' class
        itp = 0
        moleculeTypes = {}
        for mi, mol in enumerate(molecules):
            # Check if the moleculetype is already listed
            # If not, generate the topology from the chain definition
            if mol not in moleculeTypes or options['SeparateTop']:
                # Name of the moleculetype
                # XXX: The naming should be changed; now it becomes Protein_X+Protein_Y+...
                name = "+".join([chain.getname(options['-name'].value) for chain in mol])
                moleculeTypes[mol] = name

                # Write the molecule type topology. Generate topology for first chain in 
                # a molecule individually, then add the other topologies
                top = Topology(mol[0], options=options, name=name)
                for m in mol[1:]:
                    top += Topology(m, options=options)
                    
                # Have to add the connections, like the connecting network
                # Gather coordinates
                mcg, coords = list(zip(*[(j[:4], j[4:7]) for m in mol for j in m.cg(force=True)]))
                mcg         = list(mcg)

                # Run through the link list and add connections (links = cys bridges or hand specified links)
                for atomA, atomB, bondlength, forceconst in options['linkListCG']:
                    if bondlength == -1 and forceconst == -1:
                        bondlength, forceconst = options['ForceField'].special[(atomA[:2], atomB[:2])]
                    # Check whether this link applies to this group
                    atomA = atomA in mcg and mcg.index(atomA)+1
                    atomB = atomB in mcg and mcg.index(atomB)+1
                    if atomA and atomB:
                        cat = (forceconst is None) and "Constraint" or "Link"
                        top.bonds.append(Bond(
                            (atomA, atomB),
                            options    = options,
                            type       = 1,
                            parameters = (bondlength, forceconst),
                            category   = cat,
                            comments   = "Cys-bonds/special link"))
                
                # Now add bonded terms of the interchain links
                #~ for i,j in zip(top.atoms, top.atoms.origResnum):
                    #~ logging.debug('{} orig resnum {}'.format(str(i), str(j - (32 << 20))))

                for gl in top.interchainlinks:
                    atomAlist = top.atoms[(gl.resnameA, gl.resnumA)]
                    atomBlist = top.atoms[(gl.resnameB, gl.resnumB)]
                    #~ logging.debug('Interchain link {} binds to {} {} which is one of these beads {}'.format(gl, gl.resnameB, gl.resnumB - (32 << 20), atomBlist))
                    
                    try:
                        bon_con, ang_con, dih_con, vsite_con, excl_con = (options['ForceField'].sugarlinkcon[str(gl)]+5*[[]])[:5]
                        bon_par, ang_par, dih_par, vsite_par           = (options['ForceField'].sugarlinkpar[str(gl)]+4*[[]])[:4]
                    except KeyError:
                        logging.warning('Parameters missing for link {}. Update sugarlinkcon/-par dictionaries.'.format(gl))
                        continue
                        
                    for con, par in zip(bon_con, bon_par):
                        atidsA = [atomAlist[con[0]][0]]
                        atidsB = [atomBlist[con[1]][0]]
                        #logging.debug('{}, atidsA: {}, atidsB: {}'.format(gl, gl.atidsA, gl.atidsB))

                        if par[1] == None:
                            top.bonds.append(Bond(
                                options    = options,
                                atoms      = tuple(atidsA + atidsB),
                                parameters = [par[0]],          # Only the bond length
                                type       = 1,
                                comments   = '{}, {} {}'.format(str(gl), con[0], con[1]),
                                category   = "Constraint"))
                        else:
                            top.bonds.append(Bond(
                                options    = options,
                                atoms      = tuple(atidsA + atidsB),
                                parameters = par,
                                type       = 1,
                                comments   = '{}, {} {}'.format(str(gl), con[0], con[1]),
                                category   = "Glycan"))

                        # Angles
                        for con, par in zip(ang_con, ang_par):
                            atidsA = [atomAlist[i][0] for i in con[0]]
                            atidsB = [atomBlist[i][0] for i in con[1]]
                            top.angles.append(Angle(
                                options    = options,
                                atoms      = tuple(atidsA + atidsB),
                                parameters = par,
                                type       = 2,
                                comments   = '{}, {} {}'.format(str(gl), con[0], con[1]),
                                category   = "Glycan"))
                        
                        # Dihedrals
                        for con, par in zip(dih_con, dih_par):
                            atidsA = [atomAlist[i][0] for i in con[0]]
                            atidsB = [atomBlist[i][0] for i in con[1]]
                            top.dihedrals.append(Dihedral(
                                options    = options,
                                atoms      = tuple(atidsA + atidsB),
                                parameters = par,
                                type       = 2,
                                comments   = '{}, {} {}'.format(str(gl), con[0], con[1]),
                                category   = "Glycan"))
                        
                        # V-Sites
                        for con, par in zip(vsite_con, vsite_par):
                            atidsA = [atomAlist[i][0] for i in con[0]]
                            atidsB = [atomBlist[i][0] for i in con[1]]

                            top.vsites.append(Vsite(
                                options    = options,
                                atoms      = tuple(atidsA + atidsB),
                                parameters = par,
                                type       = 1,
                                comments   = '{}, {} {}'.format(str(gl), con[0], con[1]),
                                category   = "Glycan"))
                        
                        # Exclusions
                        for con in excl_con:
                            atidsA = [atomAlist[con[0]][0]]
                            atidsB = [atomBlist[con[1]][0]]
                            top.exclusions.append(Exclusion(
                                options    = options,
                                atoms      = tuple(atidsA + atidsB),
                                comments   = '{}, {} {}'.format(str(gl), con[0], con[1]),
                                parameters = (None, )))

                # Elastic Network
                # The elastic network is added after the topology is constructed, since that
                # is where the correct atom list with numbering and the full set of
                # coordinates for the merged chains are available.
                if options['ElasticNetwork']:
                    rubberType = options['ForceField'].EBondType
                    rubberList = rubberBands(
                        [(i[0], j) for i, j in zip(top.atoms, coords) if i[4] in options['ElasticBeads']],
                        options['ElasticLowerBound'], options['ElasticUpperBound'],
                        options['ElasticDecayFactor'], options['ElasticDecayPower'],
                        options['ElasticMaximumForce'], options['ElasticMinimumForce'])
                    top.bonds.extend([Bond(i, options=options, type=rubberType, category="Rubber band") for i in rubberList])

                # Write out the MoleculeType topology
                destination = options["-o"] and open(moleculeTypes[mol]+".itp", 'w') or sys.stdout
                destination.write(str(top))

                # Build conect dictionary if CG structure output is desired
                if options["-x"].value:
                    logging.info("Appending CONECT section to coarse-grained structure file '{}'.".format(options["-x"].value))
    
                    conect = {}
                    for bond in top.bonds:
                        # Do not display elastic network
                        if bond.category != 'Rubber band':    
                            try:
                                conect[bond.atoms[0]].append(bond.atoms[1])
                            except KeyError:
                                conect[bond.atoms[0]] = [bond.atoms[1]]
        
                            try:
                                conect[bond.atoms[1]].append(bond.atoms[0])
                            except KeyError:
                                conect[bond.atoms[1]] = [bond.atoms[0]]
                            
                    # Append CONECT section to the CG structure file if it exists
                    for i in conect:
                        cgOutPDB.write('CONECT' + '{0:5d}'.format(i) + "".join(['{0:5d}'.format(k) for k in conect[i]])+'\n')
                    cgOutPDB.close()
                itp += 1

            # Check whether other chains are equal to this one
            # Skip this step if we are to write all chains to separate moleculetypes
            if not options['SeparateTop']:
                for j in range(mi+1, len(molecules)):
                    if not molecules[j] in moleculeTypes and mol == molecules[j]:
                        # Molecule j is equal to a molecule mi
                        # Set the name of the moleculetype to the one of that molecule
                        moleculeTypes[molecules[j]] = moleculeTypes[mol]

        logging.info('Written %d ITP file%s.' % (itp, itp > 1 and "s" or ""))

        # WRITING THE MASTER TOPOLOGY
        # Output stream
        top  = options["-o"] and open(options['-o'].value, 'w') or sys.stdout

        # ITP file listing
        itps = '\n'.join(['#include "%s.itp"' % molecule for molecule in set(moleculeTypes.values())])

        # Molecule listing
        logging.info("Output contains %d molecules:" % len(molecules))
        n = 1
        for molecule in molecules:
            chainInfo = (n, moleculeTypes[molecule], len(molecule) > 1 and "s" or " ", " ".join([i.id for i in molecule]))
            logging.info("  %2d->  %s (chain%s %s)" % chainInfo)
            n += 1
        molecules   = '\n'.join(['%s \t 1' % moleculeTypes[molecule] for molecule in molecules])

        # Set a define if we are to use rubber bands
        useRubber   = options['ElasticNetwork'] and "#define RUBBER_BANDS" or ""

        # XXX Specify a better, version specific base-itp name.
        # Do not set a define for position restraints here, as people are more used to do it in mdp file?
        top.write(
'''#include "martini.itp"

%s

%s

[ system ]
; name
Martini system from %s

[ molecules ]
; name        number
%s''' % (useRubber, itps, options["-f"] and options["-f"].value or "stdin", molecules))

        logging.info('Written topology files')

    # Maybe there are forcefield specific log messages?
    options['ForceField'].messages()

    # The following lines are always printed (if no errors occur).
    print("\n\tThere you are. One MARTINI. Shaken, not stirred.\n")
    Q = martiniq.pop(random.randint(0, len(martiniq)-1))
    print("\n", Q[1], "\n%80s" % ("--"+Q[0]), "\n")
    
if __name__ == '__main__':
    import sys, logging
    args = sys.argv[1:]
    # Get the possible commandline arguments arguments and help text. 
    options, lists = options, lists
    # Parse commandline options.
    options = option_parser(args, options, lists, version)

    main(options)
