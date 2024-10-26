#! /usr/bin/env python
########################################################################
#            PHASE 3 (FOR MAIN_CHAIN and SIDE_CHAIN ONLY)
#PROGRAM NAME: HB_tracker v.1.0
#1. calculate the DHA angle...
#2. calculate HAAA angle...
#3. calculate DAAj angle...
#4. calculate DA and HA distances. . . 
#5. <outfile to be gennerated>
##########################################################################
from __future__ import with_statement #in some python version "with" statement does not work, so this line helps to overcome, but must be written at VERY FIRST of CODE

import timeit
start = timeit.default_timer()
import os
import math
import itertools
inp = open("/run/media/abhisek/Earth/Users/abhisek/Desktop/New_folder/x/4g78FH.atom",'r').read().strip().split('\n')
out = open("/run/media/abhisek/Earth/Users/abhisek/Desktop/New_folder/x/4g78fh.DA_cal",'w')
out2 = open("/run/media/abhisek/Earth/Users/abhisek/Desktop/New_folder/x/4g78fh.DA",'w')
hb_out = open("/run/media/abhisek/Earth/Users/abhisek/Desktop/New_folder/x/4g78fh.hbinfo",'w')  #for writing all finally calculated distance & angle values...

donor = []
H = []
acceptor = []
adjacent = []
donor_list = []
H_atom_list = []


for line in map(str.split,inp):
    atomName = line[2]
    residueName = line[3]
    chainName = line[4]
    residueNumber = line[5]
    xAxis = line[6]
    yAxis = line[7]
    zAxis = line[8]
    if atomName == 'N' and residueName[-3:] != 'PRO':
        donor.append(line)
    elif atomName == 'H':
        H.append(line)
    elif atomName == 'O':
        acceptor.append(line)
    elif atomName == 'C':
        adjacent.append(line)


#######################################################################################################################################################
#MODULE NAME:    RAT_h[rESONATING aTOM tRACKER for HISTIDINE] v.15.01.2014....
#.........completed on 15th Jan 2014.................
#PROGRAM TO CATEGORIZE DONOR AND ACCEPTOR ATOM INFORMATION FROM RESONATING IMMIDAZOLE RING OF HISTIDINE.
#NB: Histidine is the only amino acid which can act as both donor and acceptor depending upon the resonating double bond character.
#    i.e. ND1 & NE2 nitrogens while have correspondong H-atom then they act as donor, but if no H-atom is attached then they act as acceptor atom.
########################################################################################################################################################


#initiating array for accumulation of GRAND DATA cominf after each module. . .
side_donor = []
side_H = []
side_acceptor = []
side_adjacent = []

result = {"N" : [], "H" : [], "ace" : []} #initiating listed dictionary.
h_his = []
d_h_his = []   #d_ lists contains duplicate list items, created to avoid calculation easier.
donor_his = []
d_donor_his = []
acceptor_his = []
d_acceptor_his = []
his_inp_file = []
adjacent_his = [] #adjacent atoms of acceptor for histidine.
with open("/run/media/abhisek/Earth/Users/abhisek/Desktop/New_folder/x/4g78FH.atom",'r') as f:
    categories = {} #initiating a dictionary.
    for line in f:
        splitline = line.strip().split()
        #considering only HIS residue and some specified atom of interest.
        if splitline[3] == 'HIS' or splitline[3] == 'AHIS':
            if splitline[2] == "ND1" or splitline[2] == "NE2" or splitline[2] == "HE2" or splitline[2] == "HD1"  :
                new_line = splitline[2]+splitline[5] #joining selected columns for future use as a key
                identifier = "".join(new_line)[1:]  #making identifier...
                if identifier not in categories:
                    categories[identifier] = 1     #starting the counter as 1...
                else:
                    categories[identifier] = 2     #if the identifier is already there, the value becomes 2. i.e. if pair is present it becomes 2.
    #to go for second passs to create the lists.    
    f.seek(0)
    
    for line in f:
        splitline = line.strip().split()
        if splitline[3] == 'HIS' or splitline[3] == 'AHIS':
            if splitline[2] == "ND1" or splitline[2] == "NE2" or splitline[2] == "HE2" or splitline[2] == "HD1"  :
                #print splitline[-1]
                new_line = splitline[2]+splitline[5]
                identifier = "".join(new_line)[1:]
                if categories[identifier] == 1:    #means the stem just occured once.
                    if splitline[-1] == 'N':
                        result["ace"].append(line.strip())
                else:
                    #splitline[-1] indicates H or N and add it to the lists as result[H] or result[N].
                    result[splitline[-1]].append(line.strip())

    

#for type in result:
#    print (type + ' list:')+'\n'+( '\n'.join(result[type]))+'\n'
    
##############################################################################################################
#actual RAT_h program ends here.. what follows represent the procedings
#for H-bond calculations while "HIS" as side chain molecule......
#*******************************************************************<<<<<<<<
##############################################################################################################
    
for type in result:
    if type == 'H':
        for line in (result[type]):
            h_his.append(line.split())
    elif type == 'N':
        for line in (result[type]):
            donor_his.append(line.split())
    elif type == 'ace':
        for line in (result[type]):
            acceptor_his.append(line.split())


from itertools import repeat
#doubling acceptor atom list entry now...
acceptor_new = [i for item in acceptor_his for i in repeat(item,2)]
for items in acceptor_new:
    d_acceptor_his.append(items)

#Gathering adjacent atoms now... !!

inp = open("/run/media/abhisek/Earth/Users/abhisek/Desktop/New_folder/x/4g78FH.atom",'r').read().strip().split('\n')
for line in map(str.split,inp):
    residueName = line[3][-3:] #just taking last 3 alphabets so to avoid 'AHIS' instead of 'HIS'...
    if residueName == 'HIS':
        his_inp_file.append(line)
    
###############################################################
#creating a search string...
###############################################################        
for line in acceptor_his:
    if line[2] == 'ND1':
        search_string_a = 'CG'+line[3]+line[4]+line[5] #combining all columns to make searching easier in the in_file.
        search_string_b = 'CE1'+line[3]+line[4]+line[5]
        #print search_string_a
        for case in his_inp_file:
            target_string = case[2]+case[3]+case[4]+case[5]
            #print target_string
            if search_string_a == target_string:
                adjacent_his.append(case)
            elif search_string_b == target_string:
                adjacent_his.append(case)
    elif line[2] == 'NE2':
        search_string_a = 'CD2'+line[3]+line[4]+line[5]
        search_string_b = 'CE1'+line[3]+line[4]+line[5]
        for case in his_inp_file:
            target_string = case[2]+case[3]+case[4]+case[5]
            #print target_string
            if search_string_a == target_string:
                adjacent_his.append(case)
            elif search_string_b == target_string:
                adjacent_his.append(case)


#for don, hy, ace, adj in zip(d_donor_his, d_h_his, d_acceptor_his, adjacent_his):
#    print don, hy, ace, adj

#INCLUDING IN MOTHER LIST . . . > > > > > >
side_donor.extend(donor_his)
side_H.extend(h_his)
side_acceptor.extend(d_acceptor_his)
side_adjacent.extend(adjacent_his)


################################################################
#including LYSINE residue...[DONOR]
# it has 1 NZ and 3 HZ atoms. So gonna triple NZ to get zippable array structure...
################################################################
donor_lysine = []
donor_lysine_tripled = []
H_lysine = []
lysine_inp_file = []
inp = open("/run/media/abhisek/Earth/Users/abhisek/Desktop/New_folder/x/4g78FH.atom",'r').read().strip().split('\n')
for line in map(str.split,inp):
    residueName = line[3][-3:] #just taking last 3 alphabets so to avoid 'AHIS' instead of 'HIS'...
    if residueName == 'LYS':
        lysine_inp_file.append(line)

for line in lysine_inp_file:
    if line[2] == 'NZ' :
        donor_lysine.append(line)
    elif line [2] == 'HZ1' or line[2] == 'HZ2' or line[2] == 'HZ3':
        H_lysine.append(line)
#tripling the donor_lysine entry to get hold via *zip...
donor_lysine_new = [i for item in donor_lysine for i in repeat(item,3)]
for items in donor_lysine_new:
    donor_lysine_tripled.append(items)    

#for don, h in zip(donor_lysine_tripled, H_lysine):
#    print don,h
#categorization for LYSINE is completed. . .>>>>>put it in a list for calculations. . . 

#INCLUDING IN MOTHER LIST . . . > > > > > >
side_donor.extend(donor_lysine_tripled)
side_H.extend(H_lysine)

################################################################
#including ARGININE residue . . .[DONOR]
#it has 1 NE & 1 HE atom and NH1 with HH11 & HH12 and NH2 with HH21 & HH22. . .
#So, let put NE & HE in a list, and double NH1,NH2 to match with H-atom counterparts...
###################################################################
donor_arginine = []
donor_arginine_doubled = []
H_arginine = []
arginine_inp_file = []
for line in map(str.split,inp):
    residueName = line[3][-3:] #just taking last 3 alphabets so to avoid 'AARG' instead of 'ARG'...
    if residueName == 'ARG':
        arginine_inp_file.append(line)

for line in arginine_inp_file:
    if line[2] == 'NH1' or line[2] == 'NH2':
        donor_arginine.append(line)
    elif line[2] == 'HH11' or line[2] == 'HH12' or line[2] == 'HH21' or line[2] == 'HH22':
        H_arginine.append(line)
#doubling donor_arginine so far added...
donor_arginine_new = [i for item in donor_arginine for i in repeat(item,2)]
for items in donor_arginine_new:
    donor_arginine_doubled.append(items)

#now including NE and HE atoms to corresponding lists. . .
for line in arginine_inp_file:
    if line[2] == 'NE':
        donor_arginine_doubled.append(line)
    elif line[2] == 'HE':
        H_arginine.append(line)

#for don, h in zip(donor_arginine_doubled, H_arginine):
#   print don,h
#categorization for ARGININE is completed. . .>>>>>put it in a list for calculations. . .

#INCLUDING IN MOTHER LIST . . . > > > > > >
side_donor.extend(donor_arginine_doubled)
side_H.extend(H_arginine)
#for don, h in zip(side_donor, side_H): print don, h
############################################################################
#including ASPARTATE residue . . .[ACCEPTOR]
#it has 1OD1, 1OD2 and 1CG . . .
#So, putting OD1 & OD2 in a acceptor and CG in adjacent follwed by doubling adjacent list. . .
############################################################################
acceptor_asp = []
adjacent_asp = []
adjacent_asp_doubled = []
asp_inp_file = [] 
for line in map(str.split, inp):
    residueName = line[3][-3:] #just taking last 3 alphabets so to avoid 'AASP' instead of 'ASP'...
    if residueName == 'ASP':
        asp_inp_file.append(line)

for line in asp_inp_file:
    if line[2] == 'OD1' or line[2] == 'OD2':
        acceptor_asp.append(line)
    elif line[2] == 'CG':
        adjacent_asp.append(line)
#doubling adjacent CG atoms from adjacent_asp list.  . .
adjacent_asp_new = [i for item in adjacent_asp for i in repeat(item,2)]
for items in adjacent_asp_new:
    adjacent_asp_doubled.append(items)

#for ace, adj in zip(acceptor_asp, adjacent_asp_doubled):
#    print ace, adj
#categorization of ASPARTATE completed . . . >>[NOTE]>>>put it in a list for calculations. . .

#INCLUDING IN MOTHER LIST . . . > > > > > >
side_acceptor.extend(acceptor_asp)
side_adjacent.extend(adjacent_asp_doubled)

###############################################################################
#including GLUTAMATE residue. . .[ACCEPTOR]
#it has 1OE1, 1OE2 and 1 CD. . .
#So, putting OE1 & OE2 in a acceptor and CD in adjacent followed by doubling adjacent list. . .
###############################################################################
acceptor_glu = []
adjacent_glu= []
adjacent_glu_doubled = []
glu_inp_file = [] 

for line in map(str.split, inp):
    residueName = line[3][-3:] #just taking last 3 alphabets so to avoid 'AGLU' instead of 'GLU...
    if residueName == 'GLU':
        glu_inp_file.append(line)

for line in glu_inp_file:
    if line[2] == 'OE1' or line[2] == 'OE2':
        acceptor_glu.append(line)
    elif line[2] == 'CD':
        adjacent_glu.append(line)
#doubling adjacent CD atoms from adjacent_glu list.  . .
adjacent_glu_new = [i for item in adjacent_glu for i in repeat(item,2)]
for items in adjacent_glu_new:
    adjacent_glu_doubled.append(items)

#for ace, adj in zip(acceptor_glu, adjacent_glu_doubled):
#    print ace, adj
#categorization of GLUTAMATE completed. . . >>[NOTE]>>>put it in a list for calculations

#INCLUDING IN MOTHER LIST . . . > > > > > >
side_acceptor.extend(acceptor_glu)
side_adjacent.extend(adjacent_glu_doubled)

##################################################################################
#including METHIONINE residue. . .[ACCEPTOR]
#it has 1SD(as an acceptor) and adjacent CG & CE. . .
#So, putting CE & CG in adjacent list and SD in acceptor followed by doubling the acceptor list. . .
###################################################################################
acceptor_met = []
acceptor_met_doubled = []
adjacent_met = []
met_inp_file = []

for line in map(str.split, inp):
    residueName = line[3][-3:] #just taking last 3 alphabets so to avoid 'AMET instead of 'MET'...
    if residueName == 'MET':
        met_inp_file.append(line)

for line in met_inp_file:
    if line[2] == 'SD':
        acceptor_met.append(line)
    elif line[2] == 'CG' or line[2] == 'CE':
        adjacent_met.append(line)
#doubling acceptor SD atoms from acceptor_met list. . .
acceptor_met_new = [i for item in acceptor_met for i in repeat(item,2)]
for items in acceptor_met_new:
    acceptor_met_doubled.append(items)

#for ace, adj in zip(acceptor_met_doubled, adjacent_met):
#    print ace, adj
#categorization of METHIIONINE completed. . .>>[NOTE]>>>put it in a list for calculations

#INCLUDING IN MOTHER LIST . . . > > > > > >
side_acceptor.extend(acceptor_met_doubled)
side_adjacent.extend(adjacent_met)
#for ace,adj in zip(side_acceptor, side_adjacent): print ace, adj

#####################################################################################
#including ASPARAGINE(ASN) residue. . .[ACCEPTOR & DONOR]
#it's OD1 act as an Acceptor with adjacent CG atom.. and it's ND2 act as donor with HD21 & HD22 H-atoms. . .
#So, storing OD1 & CG in required lists and doubling ND2 to get zipping enabled with HD21 & HD22. . .
#####################################################################################
acceptor_asn = []
adjacent_asn = []
donor_asn = []
donor_asn_doubled = []
H_asn = []
asn_inp_file = []

for line in map(str.split, inp):
    residueName = line[3][-3:] #just taking last 3 alphabets so to avoid 'AASN' instead of 'ASN'...
    if residueName == 'ASN':
        asn_inp_file.append(line)

for line in asn_inp_file:
    if line[2] == 'OD1':
        acceptor_asn.append(line)
    elif line[2] == 'CG':
        adjacent_asn.append(line)
        
    elif line[2] == 'ND2':
        donor_asn.append(line)
    elif line[2] == 'HD21' or line[2] == 'HD22':
        H_asn.append(line)
#doubling donor 'ND2' atoms from donor_asn list. . .
donor_asn_new = [i for item in donor_asn for i in repeat(item, 2)]
for items in donor_asn_new:
    donor_asn_doubled.append(items)

#for ace, adj in zip(acceptor_asn, adjacent_asn):
#    print ace, adj
#categorization of ASPARAGINE as Acceptor  completed . . . >>[NOTE]>>> put it in subsequent list for calculations
#for don,h in zip(donor_asn_doubled, H_asn):
#    print don,h
#categorization of ASPARAGINE as Donor completed . . . >>[NOTE]>>> put it in subsequent list for calculations

#INCLUDING IN MOTHER LIST . . . > > > > > >
side_donor.extend(donor_asn_doubled)
side_H.extend(H_asn)
side_acceptor.extend(acceptor_asn)
side_adjacent.extend(adjacent_asn)
#for don, h in zip(side_donor, side_H): print don, h
#for ace, adj in zip(side_acceptor, side_adjacent): print ace, adj

######################################################################################
#including GLUTAMINE(GLN) residue . . .[DONOR & ACCEPTOR]
#it's OD1 act as an Acceptor with adjacent CG atom.. and it's ND2 act as donor with HD21 & HD22 H-atoms. . .
#So, storing OD1 & CG in required lists and doubling ND2 to get zipping enabled with HD21 & HD22. . .
#####################################################################################
acceptor_gln = []
adjacent_gln = []
donor_gln = []
donor_gln_doubled = []
H_gln = []
gln_inp_file = []

for line in map(str.split, inp):
    residueName = line[3][-3:] #just taking last 3 alphabets so to avoid 'AASN' instead of 'ASN'...
    if residueName == 'GLN':
        gln_inp_file.append(line)

for line in gln_inp_file:
    if line[2] == 'OE1':
        acceptor_gln.append(line)
    elif line[2] == 'CD':
        adjacent_gln.append(line)
        
    elif line[2] == 'NE2':
        donor_gln.append(line)
    elif line[2] == 'HE21' or line[2] == 'HE22':
        H_gln.append(line)
#doubling donor 'NE2' atoms from donor_gln list. . .
donor_gln_new = [i for item in donor_gln for i in repeat(item, 2)]
for items in donor_gln_new:
    donor_gln_doubled.append(items)

#for ace, adj in zip(acceptor_gln, adjacent_gln):
 #   print ace, adj
#categorization of ASPARAGINE as Acceptor  completed . . . >>[NOTE]>>> put it in subsequent list for calculations
#for don,h in zip(donor_gln_doubled, H_gln):
#    print don,h
#categorization of ASPARAGINE as Donor completed . . . >>[NOTE]>>> put it in subsequent list for calculations

#INCLUDING IN MOTHER LIST . . . > > > > > >
side_donor.extend(donor_gln_doubled)
side_H.extend(H_gln)
side_acceptor.extend(acceptor_gln)
side_adjacent.extend(adjacent_gln)
#for don, h in zip(side_donor, side_H): print don, h
#for ace, adj in zip(side_acceptor, side_adjacent): print ace, adj

################################################################################################
#including SERINE(SER) residue . . . [DONOR & ACCEPTOR]
#it's OG acts as both acceptor and donor atom. . .
#So, alloting OG as donor & HG as H-atom and OG as acceptor & CB as adjacent atom. . .
#################################################################################################
acceptor_ser = []
adjacent_ser = []
donor_ser = []
H_ser = []
ser_inp_file = []

for line in map(str.split, inp):
    residueName = line[3][-3:]    ##just taking last 3 alphabets so to avoid 'ASER' instead of 'SER'...
    if residueName == 'SER':
        ser_inp_file.append(line)

for line in ser_inp_file:
    if line[2] == 'OG':
        acceptor_ser.append(line)
        donor_ser.append(line)
    elif line[2] == 'CB':
        adjacent_ser.append(line)
    elif line[2] == 'HG':
        H_ser.append(line)

#for ace, adj in zip(acceptor_ser, adjacent_ser):
 #  print ace, adj
#categorization of SERINE as Acceptor  completed . . . >>[NOTE]>>> put it in subsequent list for calculations
#for don,h in zip(donor_ser, H_ser):
#  print don,h
#categorization of SERINE as Donor completed . . . >>[NOTE]>>> put it in subsequent list for calculations

#INCLUDING IN MOTHER LIST . . . > > > > > >
side_donor.extend(donor_ser)
side_H.extend(H_ser)
side_acceptor.extend(acceptor_ser)
side_adjacent.extend(adjacent_ser)
#for don, h in zip(side_donor, side_H): print don, h
#for ace, adj in zip(side_acceptor, side_adjacent): print ace, adj

        
###################################################################################################
#including THREONINE(THR) residue. . . [DONOR & ACCEPTOR]
#it's OG1 atom acts as both donor and acceptor. . .
#So, alocating OG1 in donor & HG1 in H-atom list and OG1 in acceptor & CB in adjacent list. . .
####################################################################################################
acceptor_thr = []
adjacent_thr = []
donor_thr = []
H_thr = []
thr_inp_file = []

for line in map(str.split, inp):
    residueName = line[3][-3:]    ##just taking last 3 alphabets so to avoid 'ATHR' instead of 'THR'...
    if residueName == 'THR':
        thr_inp_file.append(line)

for line in thr_inp_file:
    if  line[2] == 'OG1':
        acceptor_thr.append(line)
        donor_thr.append(line)
    elif line[2] == 'CB':
        adjacent_thr.append(line)
    elif line[2] == 'HG1':
        H_thr.append(line)

#for ace, adj in zip(acceptor_thr, adjacent_thr):
#  print ace, adj
#categorization of THREONINE as Acceptor  completed . . . >>[NOTE]>>> put it in subsequent list for calculations
#for don,h in zip(donor_thr, H_thr):
 # print don,h
#categorization of THREONINE as Donor completed . . . >>[NOTE]>>> put it in subsequent list for calculations

#INCLUDING IN MOTHER LIST . . . > > > > > >
side_donor.extend(donor_thr)
side_H.extend(H_thr)
side_acceptor.extend(acceptor_thr)
side_adjacent.extend(adjacent_thr)
#for don, h in zip(side_donor, side_H): print don, h
#for ace, adj in zip(side_acceptor, side_adjacent): print ace, adj


####################################################################################################
#including TYROSINE(TYR) residue. . . [DONOR & ACCEPTOR]
#its OH atom acts as both donor and acceptor. . .
#So, allocating OH in donor, HH in H-atom and OH in acceptor, CZ in adjacent atom list. . .
####################################################################################################
acceptor_tyr = []
adjacent_tyr = []
donor_tyr = []
H_tyr = []
tyr_inp_file = []

for line in map(str.split, inp):
    residueName = line[3][-3:]    ##just taking last 3 alphabets so to avoid 'ATYR' instead of 'TYR'...
    if residueName == 'TYR':
        tyr_inp_file.append(line)

for line in tyr_inp_file:
    if line[2] == 'OH':
        acceptor_tyr.append(line)
        donor_tyr.append(line)
    elif line[2] == 'CZ':
        adjacent_tyr.append(line)
    elif line[2] == 'HH':
        H_tyr.append(line)

#for ace, adj in zip(acceptor_tyr, adjacent_tyr):
#  print ace, adj
#categorization of TYROSINE as Acceptor  completed . . . >>[NOTE]>>> put it in subsequent list for calculations
#for don,h in zip(donor_tyr, H_tyr):
#    print don,h
#categorization of TYROSINE as Donor completed . . . >>[NOTE]>>> put it in subsequent list for calculations

#INCLUDING IN MOTHER LIST . . . > > > > > >
side_donor.extend(donor_tyr)
side_H.extend(H_tyr)
side_acceptor.extend(acceptor_tyr)
side_adjacent.extend(adjacent_tyr)
#for don, h in zip(side_donor, side_H): print don, h
#for ace, adj in zip(side_acceptor, side_adjacent): print ace, adj


#####################################################################################################
#including TRYPTOPHAN(TRP) residue . . . [DONOR]
#it's NE1 atom as donor and HE1 is corresponding H-atom. . .
#So, allocating NE1 in donor & HE1 in H-atom list. . .
#####################################################################################################
donor_trp = []
H_trp = []
trp_inp_file = []

for line in map(str.split, inp):
    residueName = line[3][-3:]    ##just taking last 3 alphabets so to avoid 'ATRP' instead of 'TRP'...
    if residueName == 'TRP':
        trp_inp_file.append(line)

for line in trp_inp_file:
    print line
    if line[2] == 'NE1':
        donor_trp.append(line)
    elif line[2] == 'HE1':
        H_trp.append(line)

#for don,h in zip(donor_trp, H_trp):
#    print don,h
#categorization of TYPTOPHAN as Donor completed . . . >>[NOTE]>>> put it in subsequent list for calculations

#INCLUDING IN MOTHER LIST . . . > > > > > >
side_donor.extend(donor_trp)
side_H.extend(H_trp)
#for don, h,ace,adj in zip(side_donor, side_H,side_acceptor, side_adjacent): print don, h, ace, adj

###########################################################
#>>>> EXTENDING MAIN_CHAIN_COORDINATES WITH SIDE_CHAIN_COORDINATES... >>>>
donor.extend(side_donor)
H.extend(side_H)
acceptor.extend(side_acceptor)
adjacent.extend(side_adjacent)
        
##################################################################################################################
#PREPARING FOR CALCULATION NOW. . .
#......... :)              :)           :) ..................... >>>>>
###################################################################################################################
for don,h in  zip(donor,H):
    #print don
    don_xAxis = don[6]
    don_yAxis = don[7]
    don_zAxis = don[8]
    h_xAxis = h[6]
    h_yAxis = h[7]
    h_zAxis = h[8]
    out.write('\n')
    #dx = 'DONOR', ' '*20 , 'ACCEPTOR', ' '*15, 'D-A'
    #header = ''.join(map(''.join,dx))
        
    donorAtom = ('donor: ' + don[2]),don[3],don[4],don[5],don[6],don[7],don[8]
    final_donorAtom = ' '.join(map(''.join,donorAtom))
       
    hAtom = ('H-atom: ' + h[2]),h[3],h[4],h[5],h[6],h[7],h[8]
    final_hAtom = ' '.join(map(''.join,hAtom))
   #print header
    #out.write(header)
    #print final_donorAtom
        
    
    for ace,adj in zip(acceptor, adjacent):
        ace_xAxis = ace[6]
        ace_yAxis = ace[7]
        ace_zAxis = ace[8]
        xy = ('acceptor: '+ace[2]),ace[3],ace[4],ace[5],ace[6],ace[7],ace[8]
        DA = ' '.join(map(''.join,xy))
        adj = ('adjacent : '+adj[2]),adj[3],adj[4],adj[5],adj[6],adj[7],adj[8]
        adjFinal = ' '.join(map(''.join,adj))
        ##########################################################
        #Generating file in a format to ease angle calculation(*.DA_cal)...
        ###########################################################
        out.write(final_donorAtom) #getting donor atom
        out.write(' ')
        out.write(final_hAtom)     #getting H-atom
        out.write(' ')
        out.write(DA)              #getting acceptor atom
        out.write(' ')
        out.write(adjFinal)        #getting acceptor-adjacent atom
        out.write('\n')
        ##########################################################
        #now generating same file but in READABLE format...
        ##########################################################
        out2.write('\n')
        out2.write(final_donorAtom) #getting donor atom
        out2.write('\n')
        out2.write(final_hAtom)     #getting H-atom
        out2.write('\n')
        out2.write(DA)              #getting acceptor atom
        out2.write('\n')
        out2.write(adjFinal)        #getting acceptor-adjacent atom
        out2.write('\n')
          
print " D-A distance generated successfully..."
print ".DA_cal and .DA outfile created..."
out.close()
out2.close()


#distance list for storage . . .
DA = []
HA = []
#angle list for storage...
DHA_list = []   #contain all calculated DHA angle values.
HAAj_list = []  #contain all calculated HAAj angle values.
DAAj_list = []  #contain all calculated DAAj angle values.
inp2_list = []
inp2_list_new = []   #list made after removing trailing spaces. . . 
inp2 = open("/run/media/abhisek/Earth/Users/abhisek/Desktop/New_folder/x/4g78fh.DA_cal",'r').read().strip().split('\n')
for lines in map(str.split,inp2):
    inp2_list.append(lines)
inp2_list_new = filter(None, inp2_list) #removing extra new_line characters to avoid index error. . . 
for lines in inp2_list_new:
    donor_atom = lines[1]
    donor_res_name = lines[2]
    donor_chain = lines[3]
    donor_res_no = lines[4]
    donor_xAxis = float(lines[5])
    donor_yAxis = float(lines[6])
    donor_zAxis = float(lines[7])

    hy_atom = lines[9]
    hy_res_name = lines[10]
    hy_chain = lines[11]
    hy_res_no = lines[12]
    hy_xAxis = float(lines[13])
    hy_yAxis = float(lines[14])
    hy_zAxis = float(lines[15])

    ace_atom = lines[17]
    ace_res_name = lines[18]
    ace_chain = lines[19]
    ace_res_no = lines[20]
    ace_xAxis = float(lines[21])
    ace_yAxis = float(lines[22])
    ace_zAxis = float(lines[23])

    adj_atom = lines[26]
    adj_res_name = lines[27]
    adj_chain = lines[28]
    adj_res_no = lines[29]
    adj_xAxis = float(lines[30])
    adj_yAxis = float(lines[31])
    adj_zAxis = float(lines[32])
    

    #print donor_xAxis,hy_xAxis,ace_xAxis,adj_xAxis
    if donor_atom == ace_atom and donor_res_name == ace_res_name and donor_chain == ace_chain and donor_res_no == ace_res_no:
        pass
    else:
        ###########################################################################
        #>>> STARTING DISTANCE CALCULATIONS. . .
        #Distance_1: DA . . .
        ###########################################################################
        xD = (ace_xAxis) - (donor_xAxis)
        yD = (ace_yAxis) - (donor_yAxis)
        zD = (ace_zAxis) - (donor_zAxis)
        DA_distance = math.sqrt(xD*xD + yD*yD + zD*zD) #distance b/w donor-acceptor.
        DA.append(DA_distance)

        ##########################################################################
        #Distance_2: HA. . .
        ##########################################################################
        xD1 = (ace_xAxis) - (hy_xAxis)
        yD1 = (ace_yAxis) - (hy_yAxis)
        zD1 = (ace_zAxis) - (hy_zAxis)
        HA_distance = math.sqrt(xD1*xD1 +yD1*yD1 + zD1*zD1)
        HA.append(HA_distance)

        ##########################################################################
        #>>> STARTING ANGLE CALCULATIONS. . . 
        #Angle_1: DHA...
        #Calculating magnitude of the two vectors HD & HA i.e. DHA ANGLE . . .!!!
        ##########################################################################
        #Calculating vector HD . . .(NB: in vector HD & DH is not same. Draw picture to explain to urself. HERE calculating vectors keeping 'H' as central point.)
        x1 = (donor_xAxis - hy_xAxis)
        y1 = (donor_yAxis - hy_yAxis)
        z1 = (donor_zAxis - hy_zAxis)
        #Calculating vector HA . . .
        x2 = (ace_xAxis - hy_xAxis)
        y2 = (ace_yAxis - hy_yAxis)
        z2 = (ace_zAxis - hy_zAxis)
        #Calculating cosine and rev-cosine function. . . 
        numerator = ((x1*x2)+(y1*y2)+(z1*z2))
        #print numerator
        denominator = ((math.sqrt(math.pow(x1,2)+math.pow(y1,2)+math.pow(z1,2)))*(math.sqrt(math.pow(x2,2)+math.pow(y2,2)+math.pow(z2,2))))
        #print denominator
        cosine = numerator/denominator
        #print cosine
        theta_rad = math.acos(cosine) #inverse cos function gives angle measure in radian...
        DHA_angle = math.degrees(theta_rad) #converted to degree...
        #if DHA_angle >= 111 and DHA_angle <= 112:
            #print donor_atom,donor_res_name,hy_atom,hy_res_name,ace_atom,ace_res_name,adj_atom,adj_res_name,numerator,denominator, DHA_angle
        DHA_list.append(DHA_angle)
        #print 'DHA = '+ str(DHA_angle)
      

        ##########################################################################
        #Angle_2: HAAj...(Aj = adjacent atom to Acceptor atom)
        #Calculating magnitude of the two vectors AH & Aj` i.e. HAAj ANGLE . . .!!!
        ##########################################################################
        #Calculating vector AH. . .
        x3 = (hy_xAxis - ace_xAxis)
        y3 = (hy_yAxis - ace_yAxis)
        z3 = (hy_zAxis - ace_zAxis)
        #Calculating vector AA'. . .
        x4 = (adj_xAxis - ace_xAxis)
        y4 = (adj_yAxis - ace_yAxis)
        z4 = (adj_zAxis - ace_zAxis) 
        #Calculating cosine and rev-cosine function. . .
        numerator2 = ((x3*x4)+(y3*y4)+(z3*z4))
        denominator2 = ((math.sqrt(math.pow(x3,2)+math.pow(y3,2)+math.pow(z3,2)))*(math.sqrt(math.pow(x4,2)+math.pow(y4,2)+math.pow(z4,2))))
        cosine2 = numerator2/denominator2
        theta_rad2 = math.acos(cosine2)
        HAAj_angle = math.degrees(theta_rad2)
        HAAj_list.append(HAAj_angle)
        #print 'HAAj = ' + str(HAAj_angle)
        
        ##########################################################################
        #Angle_3: DAAj...(Aj = adjacent atom to Acceptor atom)
        #Calculating magnitude of the two vectors AD & Aj` i.e. DAAj ANGLE . . .!!!
        ##########################################################################
        #Calculating vector AD. . .
        x5 = (donor_xAxis - ace_xAxis)
        y5 = (donor_yAxis - ace_yAxis)
        z5 = (donor_zAxis - ace_zAxis)
    
        #Calculating vector AA'. . .
        x6 = (adj_xAxis - ace_xAxis)
        y6 = (adj_yAxis - ace_yAxis)
        z6 = (adj_zAxis - ace_zAxis)
        #Calculating cosine and rev-cosine function. . .
        numerator3 = ((x5*x6)+(y5*y6)+(z5*z6))
        #print numerator3
        denominator3 = ((math.sqrt(math.pow(x5,2)+math.pow(y5,2)+math.pow(z5,2)))*(math.sqrt(math.pow(x6,2)+math.pow(y6,2)+math.pow(z6,2))))
        #print denominator3
        cosine3 = numerator3/denominator3
        theta_rad3 = math.acos(cosine3)
        DAAj_angle = math.degrees(theta_rad3)
        DAAj_list.append(DAAj_angle)
        #print 'DAAj = ' + str(DAAj_angle)

############################################################       
#>> checkpoint inactivated for now...                          #
final_list = []
dft_input_list = []    #CONTAINS DATA FOR NWCHEM INPUT.. >>>>>>
count = 0
h1 = str("_"*150)  #Adding a header line for viewing pleasure...
hb_out.write(h1)
hb_out.write('\n')
header_line = str( "DONOR"+" "*20+"HYDROGEN-ATOM"+" "*12+"ACCEPTOR ATOM"+" "*12+"ADJACENT ATOM"+" "*10+"DA"+" "*6+"HA"+" "*7+"DHA"+" "*7+"HAAj"+" "*9+"DAAj")   #headers line. . .
hb_out.write(header_line)
hb_out.write('\n')
hb_out.write(h1)
hb_out.write('\n')
HB_output_list=[]  #storing the replica of hydrogen bond output, for later ease of modification...@5.1.15
for lines,da,ha,items,case,entry in zip(inp2_list_new,DA,HA,DHA_list,HAAj_list,DAAj_list):
    if da <= 3.9 and ha <= 2.5 and items > 90 and case >90 and entry > 90:
        #print lines,da,ha,items,case,entry
        count += 1        #print lines, ('DA: '+ str(da)), ('HA: '+ str(ha)), ('DHA: ' + str(items)), ('HAAj: ' + str(case)), ('DAAj: ' + str(entry))
        #output = lines[1], lines[2], lines[3],lines[4],'    ', lines[9],lines[10],lines[11],lines[12],'    ',lines[17],lines[18],lines[19],lines[20],'    ',lines[26],lines[27],lines[28],lines[29],'    ', ('DA: %.4f ' % da),'    ',('HA: %.4f' % ha),'    ', ('DHA: %.4f' % items),'    ', ('HAAj: %.4f' % case),'    ', ('DAAj: %.4f' % entry)       

        donor_lines =str(str(lines[1])+" "*(6-len(lines[1]))+str(lines[2])+" "*(6-len(lines[2]))+ str(lines[3])+" "*(4-len(lines[3]))+str(lines[4]))
        dft_donor_lines = str(str(lines[1][0])+' '+str(lines[5])+' '+str(lines[6])+' '+str(lines[7]))    #to be used as DFT input data for DONOR_ATOM...
        #donor_lines_list.append(donor_lines)
        hydrogen_lines = str(" "*(8-len(lines[4]))+str(lines[9])+" "*(6-len(lines[9]))+str(lines[10])+" "*(6-len(lines[10]))+str(lines[11])+" "*(4-len(lines[11]))+lines[12])
        dft_hydrogen_lines = str(str(lines[9][0])+' '+str(lines[13])+' '+str(lines[14])+' '+str(lines[15])) #to be used as DFT input data for HYDROGEN_ATOM...
        #hydrogen_lines_list.append(hydrogen_lines)
        acceptor_lines = str(" "*(8-len(lines[12]))+str(lines[17])+" "*(6-len(lines[17]))+str(lines[18])+" "*(6-len(lines[18]))+str(lines[19])+" "*(4-len(lines[19]))+str(lines[20]))
        dft_acceptor_lines = str(str(lines[17][0])+' '+str(lines[21])+' '+str(lines[22])+' '+str(lines[23]))  #to be used as DFT input data for ACCEPTOR_ATOM...
        #acceptor_lines_list.append(acceptor_lines)
        adjacent_lines = str(" "*(8-len(lines[20]))+str(lines[26])+" "*(6-len(lines[26]))+str(lines[27])+" "*(6-len(lines[27]))+str(lines[28])+" "*(4-len(lines[28]))+str(lines[29]))
        dft_adjacent_lines = str(str(lines[26][0])+' '+str(lines[30])+ ' '+str(lines[31])+' '+str(lines[32]))  #to be used as DFT input data for ADJACENT_ATOMS...
        
        total_dft_input = dft_donor_lines+ ';' +dft_hydrogen_lines+';'+dft_acceptor_lines+';'+dft_adjacent_lines
        dft_input_list.append(total_dft_input)         #storing in the list for mass usage later on...
        #adjacent_lines_list.append(adjacent_lines)
        
        da_distance = str(" "*(6-len(lines[29]))+'DA:%.2f' % (da))
        ha_distance = str(" "*(6-len(da_distance))+'HA:%.2f' % (ha))
        dha_angle =  str(" "*(6-len(ha_distance))+'DHA:%.2f' % (items))
        #print dha_angle
        haaj_angle = str(" "*(6-len(dha_angle)) +'HAAj:%.2f' % (case))
        daaj_angle = str(" "*(11-len(haaj_angle))+'DAAj:%.2f' % (entry))
        
        zz =  donor_lines, hydrogen_lines, acceptor_lines, adjacent_lines,da_distance,ha_distance,dha_angle,haaj_angle,daaj_angle
        HB_output_list.append(zz)        
        zzz = ' '.join(zz)
        hb_out.write(zzz)
        hb_out.write('\n')

        #formatted_output = (' '.join((output)))
        #z = str(formatted_output)
        #final_list.append(z)
#for items in HB_output_list:
#    print items
'''
################################################################################################################        
#################################################################################################################
#       ...............PREPARING FOR DFT INPUT CONFIGURATION...............................
#               <<<data are stored in list named "dft_input_list">>>>
#################################################################################################################
import sys
import subprocess
dft_energy = []
for items in dft_input_list:
    #split_items = items.split()
    #if split_items[0]=='O':
    #    print 'its an O items'
    #elif split_items[0]=='N':
    #    print 'its a N item'
    f = open("/run/media/abhisek/Earth/Users/abhisek/Desktop/New_folder/dft/e_hydrogen_bond_template.nw",'r')
    filedata=f.read()
    f.close()
    newdata = filedata.replace("  N##",("  "+items))
    z = open("/run/media/abhisek/Earth/Users/abhisek/Desktop/New_folder/dft/e_hydrogen_bond_exe.nw",'w')
    z.write(newdata)
    z.close()
    dft_out = open("/run/media/abhisek/Earth/Users/abhisek/Desktop/New_folder/dft/e_hydrogen_bond.data",'w')
    p =subprocess.Popen(['nwchem ' +"/run/media/abhisek/Earth/Users/abhisek/Desktop/New_folder/dft/e_hydrogen_bond_exe.nw"], stdout=subprocess.PIPE, shell=True)
    output, err = p.communicate()
    for line in output:
        dft_out.write(line)
    dft_out.close()
    dft_file = open("/run/media/abhisek/Earth/Users/abhisek/Desktop/New_folder/dft/e_hydrogen_bond.data",'r').read().strip().split('\n')
    total_lines = []      #storing all nwchem generated DFT & BSSE energies...
    useful_energies = []  #filtering from list "total_lines" and storing only 'isolated' and 'BSSE corrected energies'...
    for lines in map(str.split,dft_file):
        if lines != []:
            if lines[0] == 'Total' and lines[1] == 'DFT' and lines[2] == 'energy':
                total_lines.append(lines)
            else:
                if lines[0] == 'Corrected' and lines[1] == 'energy':
                    total_lines.append(lines)
    isolated_energy1 = float(total_lines[1][-1])
    isolated_energy2 = float(total_lines[3][-1])
    BSSE_energy = float(total_lines[5][-1])
    net_energy = (BSSE_energy - (isolated_energy1 + isolated_energy2))  # generates energy in Hartree...
    net_energy_kcal = net_energy * 627.509 #Converting energy into kcal/mol unit...
    #print ('DFT Energy = %f' % net_energy_kcal +' kcal/mol')
    dft_energy.append(net_energy_kcal)

#MODULE REARRANGEMENT NECESSARY...    
for items in dft_energy:
    print items
'''
print "VOILA!! you have successfully computed all types of possible distances and angles. . ."
print (str(count) +" probable hydrogen bond calculated")
################################################################
stop = timeit.default_timer()
print ("Total time taken for the program to complete: %.3f sec" % (stop - start))
print "OUTPUT FILE GENERATED SUCCESSFULLY . . ."
