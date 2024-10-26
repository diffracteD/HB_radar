#! /usr/bin/env python

'''
################################################################################################################
                     PHASE 6: Generating all H-bond info.
                    PROGRAM NAME: HB_radar_ClusterClaculator_m062x.v.2.0 (05.05.2018)
_______________________________________________________________________________________________________________
1. calculate the DHA angle...
2. calculate HAAA angle...
3. calculate DAAj angle...
4. calculate DA and HA distances. . . 
5. all distance & angle outfile generated as .hbinfo
6. {MODIFICATION}: It calculates only CLUSTER (given certain radius) energy using DFT
    PS: It does not calculates BSSE or interaction energy...
7. Calculates Cluster energy with dummy/ghost DH-AA' atoms
8. Subtract ghost cluster energy from Cluster Energy...
9. dft energy + .hbinfo is outfiled as .HBradar
##########################################################################
New Func v.3.1: Multiple Exception Handler...
ALL ATOM HANDLER... :)
'''
from __future__ import with_statement #in some python version "with" statement does not work, so this line helps to overcome, 
                                        # but must be written at VERY FIRST of CODE


def show_me_hbs(path):

    from time import sleep
    import sys
    import time
    import timeit
    start = timeit.default_timer()
    import math
    import itertools
    import os
    from os.path import basename
    import glob  

    #EMPLOYING THE WALKING MODULE...
    for root, dirs, files in os.walk(path):
        list_of_file = glob.iglob(os.path.join(root,'*.atom'))
        listOfOpenedFiles = []
        print "            ENTERING DIRECTORY: " + root
        for atom_file_name in list_of_file:
            #count += 1
            if os.stat(atom_file_name)[6] != 0:
                filename = basename(atom_file_name)
                listOfOpenedFiles.append(filename)

                #Now opening the .atom file as input...
                inp = open(atom_file_name,'r').read().strip().split('\n')


                print "> PROCESSING FILE: " + filename  #shows current file under attack...

                loop_clock_start = timeit.default_timer() #Strarting loop timer $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$                
                #-----------------------------------------------------------------------------------
                #Opening oufiles to be written in later as data comes...[special walking_module format]
                #Opening DA_cal file while on walk...
                da_cal = filename.replace('atom','DA_cal')
                da_cal_out = open(os.path.join(root,da_cal),'w')

                #Opening DA file while on walk...
                da = filename.replace('atom','DA')
                da_out = open(os.path.join(root, da),'w')


                #--------------------------------------------------------------------------                               


                #inp = open("/home/saumen/abhisek/dft_cal/xpt_zone_3_fine-opt/4g78FH.atom",'r').read().strip().split('\n')
                #out = open("/home/saumen/abhisek/dft_cal/xpt_zone_3_fine-opt/4g78FH.DA_cal",'w')
                #out2 = open("/home/saumen/abhisek/dft_cal/xpt_zone_3_fine-opt/4g78FH.DA",'w')
                #hb_out = open("/home/saumen/abhisek/dft_cal/xpt_zone_3_fine-opt/4g78FH.hbinfo",'w')  #for writing all finally calculated distance & angle values...

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

                    #making-out donor n hydrogen atoms...
                    if atomName =='N' and residueName[-3] != 'PRO':
                        a = residueNumber
                        x = line
                    elif atomName == 'H':
                        b = residueNumber
                        y = line
                        if a == b:  #making sure we are taking pair of same residue i.e. N3-H3 not N3-H4
                            donor.append(x)
                            H.append(y)

                    #making-out acceptor n adjacent atoms...
                    #if atomName == 'O':
                    #    acceptor.append(line)
                    #elif atomName == 'C':
                    #    adjacent.append(line)

                    #appending donor-adjacent atoms...
                for line in map(str.split,inp):
                    if line[2]=='O':
                        acceptor_search_str_a = 'C'+line[3]+line[4]+line[5]
                        for item in map(str.split,inp):
                            if item[2]=='C':
                                adjacent_search_str_b = item[2]+item[3]+item[4]+item[5]
                                if acceptor_search_str_a == adjacent_search_str_b:
                                    #print line, item
                                    acceptor.append(line)
                                    adjacent.append(item)




                '''
                ==================================================================================================================================================
                MODULE NAME:    RAT_h[rESONATING aTOM tRACKER for HISTIDINE] v.15.01.2014....
                .........completed on 15th Jan 2014.................
                PROGRAM TO CATEGORIZE DONOR AND ACCEPTOR ATOM INFORMATION FROM RESONATING IMMIDAZOLE RING OF HISTIDINE.
                NB: Histidine is the only amino acid which can act as both donor and acceptor depending upon the resonating double bond character.
                i.e. ND1 & NE2 nitrogens while have correspondong H-atom then they act as donor, but if no H-atom is attached then they act as acceptor atom.
                ==================================================================================================================================================
                '''

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
                with open(atom_file_name,'r') as f:
                    categories = {} #initiating a dictionary.
                    for line in f:
                        splitline = line.strip().split()
                        #considering only HIS residue and some specified atom of interest.
                        if splitline[3] == 'HIS' or splitline[3] == 'AHIS':
                            if splitline[2] == "ND1" or splitline[2] == "NE2" or splitline[2] == "HE2" or splitline[2] == "HD1":
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

                #-----------------------------------------------------------------------------------------
                #actual RAT_h program ends here.. what follows represent the procedings
                #for H-bond calculations while "HIS" as side chain molecule......
                #*******************************************************************<<<<<<<<
                #-----------------------------------------------------------------------------------------

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
                #acceptor_new = [i for item in acceptor_his for i in repeat(item,2)]


                #Gathering adjacent atoms now... !!

                inp = open(atom_file_name,'r').read().strip().split('\n')
                for line in map(str.split,inp):
                    residueName = line[3][-3:] #just taking last 3 alphabets so to avoid 'AHIS' instead of 'HIS'...
                    if residueName == 'HIS':
                        his_inp_file.append(line)

                #---------------------------------------------------------------------
                #creating a search string...
                #---------------------------------------------------------------------        
                for item in his_inp_file:
                    if item[2] == 'ND1':
                        search_string_a = 'CG'+item[3]+item[4]+item[5] #combining all columns to make searching easier in the in_file.
                        search_string_b = 'CE1'+item[3]+item[4]+item[5]
                        #print search_string_a
                        for case in his_inp_file:
                            target_string_x = case[2]+case[3]+case[4]+case[5]
                            #print target_string
                            if search_string_a == target_string_x:
                                adjacent_his.append(case)
                                d_acceptor_his.append(item)
                            elif search_string_b == target_string_x:
                                adjacent_his.append(case)
                                d_acceptor_his.append(item)
                    elif item[2] == 'NE2':
                        search_string_c = 'CD2'+item[3]+item[4]+item[5]
                        search_string_d = 'CE1'+item[3]+item[4]+item[5]
                        for case in his_inp_file:
                            target_string_y = case[2]+case[3]+case[4]+case[5]
                            #print target_string
                            if search_string_c == target_string_y:
                                adjacent_his.append(case)
                                d_acceptor_his.append(item)
                            elif search_string_d == target_string_y:
                                adjacent_his.append(case)
                                d_acceptor_his.append(item)


                #for don, hy, ace, adj in zip(d_donor_his, d_h_his, d_acceptor_his, adjacent_his):
                #    print don, hy, ace, adj

                #INCLUDING IN MOTHER LIST . . . > > > > > >
                side_donor.extend(donor_his)
                side_H.extend(h_his)
                side_acceptor.extend(d_acceptor_his)
                side_adjacent.extend(adjacent_his)


                #------------------------------------------------------------------------------------
                #including LYSINE residue...[DONOR]
                # it has 1 NZ and 3 HZ atoms. So gonna triple NZ to get zippable array structure...
                #------------------------------------------------------------------------------------
                donor_lysine = []
                donor_lysine_tripled = []
                H_lysine = []
                lysine_inp_file = []
                inp = open(atom_file_name,'r').read().strip().split('\n')
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

                #---------------------------------------------------------------------------------------
                #including ARGININE residue . . .[DONOR]
                #it has 1 NE & 1 HE atom and NH1 with HH11 & HH12 and NH2 with HH21 & HH22. . .
                #So, lets put NE & HE in a list, and double NH1,NH2 to match with H-atom counterparts...
                #---------------------------------------------------------------------------------------
                donor_arginine = []
                donor_arginine_doubled = []
                H_arginine = []
                arginine_inp_file = []
                for line in map(str.split,inp):
                    residueName = line[3][-3:] #just taking last 3 alphabets so to avoid 'AARG' instead of 'ARG'...
                    if residueName == 'ARG':
                        arginine_inp_file.append(line)

                for line in arginine_inp_file:
                    residueNumber = line[5] #has to redefine the residue number as reopening...
                    if line[2] == 'NH1' or line[2] == 'NH2':
                        p = residueNumber
                        #print line, p
                        arg_x = line
                        donor_arginine.append(arg_x)
                    elif line[2] == 'HH11' or line[2] == 'HH12' or line[2] == 'HH21' or line[2] == 'HH22' or line[2] == 'HH1' or line[2] =='HH2':
                        q = residueNumber
                        arg_y = line
                        H_arginine.append(arg_y)

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


                #INCLUDING IN MOTHER LIST . . . > > > > > >
                side_donor.extend(donor_arginine_doubled)
                side_H.extend(H_arginine)
                #for don, h in zip(side_donor, side_H): print don, h

                #---------------------------------------------------------------------------------------------
                #including ASPARTATE residue . . .[ACCEPTOR]
                #it has 1OD1, 1OD2 and 1CG . . .
                #So, putting OD1 & OD2 in a acceptor and CG in adjacent follwed by doubling adjacent list. . .
                #----------------------------------------------------------------------------------------------
                acceptor_asp = []
                adjacent_asp = []
                adjacent_asp_doubled = []
                asp_inp_file = [] 
                for line in map(str.split, inp):
                    residueName = line[3][-3:] #just taking last 3 alphabets so to avoid 'AASP' instead of 'ASP'...
                    if residueName == 'ASP':
                        asp_inp_file.append(line)

                for line in asp_inp_file:
                    if line[2] == 'OD1':
                        asp_search_string_a = 'CG'+line[3]+line[4]+line[5]
                        for item in asp_inp_file:
                            if item[2] == 'CG':
                                asp_target_string = item[2]+item[3]+item[4]+item[5]
                                if asp_search_string_a == asp_target_string:
                                    adjacent_asp.append(item)
                                    acceptor_asp.append(line)

                    elif line[2] == 'OD2':
                        asp_search_string_b = 'CG'+line[3]+line[4]+line[5]
                        for item in asp_inp_file:
                            if item[2] == 'CG':
                                asp_target_string = item[2]+item[3]+item[4]+item[5]
                                if asp_search_string_b == asp_target_string:
                                    adjacent_asp.append(item)   
                                    acceptor_asp.append(line)


                    #resNum = line[5]
                    #if line[2] == 'OD1' or line[2] == 'OD2':
                    #    asp_x = line[5]  #residue number
                    #    asp_ace_line = line
                    #    acceptor_asp.append(asp_ace_line)
                    #elif line[2] == 'CG':
                    #    asp_y = line[5]  #residue number
                    #    asp_adj_line = line
                    #    adjacent_asp.append(asp_adj_line)

                #doubling adjacent CG atoms from adjacent_asp list.  . .
                #adjacent_asp_new = [i for item in adjacent_asp for i in repeat(item,2)]
                #for items in adjacent_asp_new:
                    #adjacent_asp_doubled.append(items)

                #for ace, adj in zip(acceptor_asp, adjacent_asp_doubled):
                #    print ace, adj
                #categorization of ASPARTATE completed . . . >>[NOTE]>>>put it in a list for calculations. . .

                #INCLUDING IN MOTHER LIST . . . > > > > > >
                side_acceptor.extend(acceptor_asp)
                side_adjacent.extend(adjacent_asp)

                #------------------------------------------------------------------------------------------------
                #including GLUTAMATE residue. . .[ACCEPTOR]
                #it has 1OE1, 1OE2 and 1 CD. . .
                #So, putting OE1 & OE2 in a acceptor and CD in adjacent followed by doubling adjacent list. . .
                #-------------------------------------------------------------------------------------------------
                acceptor_glu = []
                adjacent_glu= []
                adjacent_glu_doubled = []
                glu_inp_file = [] 

                for line in map(str.split, inp):
                    residueName = line[3][-3:] #just taking last 3 alphabets so to avoid 'AGLU' instead of 'GLU...
                    if residueName == 'GLU':
                        glu_inp_file.append(line)

                for line in glu_inp_file:
                    #print line
                    if line[2] == 'OE1':
                        glu_search_string_a = 'CD'+line[3]+line[4]+line[5]
                        for item in glu_inp_file:
                            if item[2] == 'CD':
                                glu_target_string_b = item[2]+item[3]+item[4]+item[5]
                                if glu_search_string_a == glu_target_string_b:
                                    adjacent_glu.append(item)
                                    acceptor_glu.append(line)

                    elif line[2] == 'OE2':
                        glu_search_string_c = 'CD'+line[3]+line[4]+line[5]
                        for item in glu_inp_file:
                            if item[2] == 'CD':
                                glu_target_string_d = item[2]+item[3]+item[4]+item[5]
                                if glu_search_string_c == glu_target_string_d:
                                    adjacent_glu.append(item)
                                    acceptor_glu.append(line)






                #for ace, adj in zip(acceptor_glu, adjacent_glu_doubled):
                #    print ace, adj
                #categorization of GLUTAMATE completed. . . >>[NOTE]>>>put it in a list for calculations

                #INCLUDING IN MOTHER LIST . . . > > > > > >
                side_acceptor.extend(acceptor_glu)
                side_adjacent.extend(adjacent_glu)
                #print acceptor_glu,adjacent_glu

                #--------------------------------------------------------------------------------------------------------
                #including METHIONINE residue. . .[ACCEPTOR]
                #it has 1SD(as an acceptor) and adjacent CG & CE. . .
                #So, putting CE & CG in adjacent list and SD in acceptor followed by doubling the acceptor list. . .
                #---------------------------------------------------------------------------------------------------------
                acceptor_met = []
                #acceptor_met_doubled = []
                adjacent_met = []
                met_inp_file = []

                for line in map(str.split, inp):
                    residueName = line[3][-3:] #just taking last 3 alphabets so to avoid 'AMET instead of 'MET'...
                    if residueName == 'MET':
                        met_inp_file.append(line)

                for line in met_inp_file:
                    #a new type of loop is used to efficiently treat 2 ADJACENT atoms per 1 ACCEPTOR atom...
                    if line[2] == 'SD':
                        search_string_a = 'CG'+line[3]+line[4]+line[5]
                        search_string_b = 'CE'+line[3]+line[4]+line[5]
                        for case in met_inp_file:
                            target_string_x = case[2]+case[3]+case[4]+case[5]
                            if search_string_a == target_string_x:
                                acceptor_met.append(line)
                                adjacent_met.append(case)
                        for case in met_inp_file:
                            target_string_y == case[2]+case[3]+case[4]+case[5]
                            if search_string_b == target_string_y:
                                acceptor_met.append(line)
                                adjacent_met.append(case)


                #for ace, adj in zip(acceptor_met_doubled, adjacent_met):
                #    print ace, adj
                #categorization of METHIIONINE completed. . .>>[NOTE]>>>put it in a list for calculations

                #INCLUDING IN MOTHER LIST . . . > > > > > >
                side_acceptor.extend(acceptor_met)
                side_adjacent.extend(adjacent_met)
                #for ace,adj in zip(side_acceptor, side_adjacent): print ace, adj

                #------------------------------------------------------------------------------------------------------------
                #including ASPARAGINE(ASN) residue. . .[ACCEPTOR & DONOR]
                #it's OD1 act as an Acceptor with adjacent CG atom.. and it's ND2 act as donor with HD21 & HD22 H-atoms. . .
                #So, storing OD1 & CG in required lists and doubling ND2 to get zipping enabled with HD21 & HD22. . .
                #-------------------------------------------------------------------------------------------------------------
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
                        #print line
                        asn_search_string_a = 'CG'+line[3]+line[4]+line[5]
                        for item in asn_inp_file:
                            if item[2] == 'CG':
                                #print item
                                asn_target_string = item[2]+item[3]+item[4]+item[5]
                                if asn_search_string_a == asn_target_string:
                                    adjacent_asn.append(item)
                                    acceptor_asn.append(line)
                                    #print asn_ace, item

                    #if line[2] == 'OD1':
                    #    acceptor_asn.append(line)
                    #elif line[2] == 'CG':
                    #    adjacent_asn.append(line)

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

                #------------------------------------------------------------------------------------------------------------
                #including GLUTAMINE(GLN) residue . . .[DONOR & ACCEPTOR]
                #it's OD1 act as an Acceptor with adjacent CG atom.. and it's ND2 act as donor with HD21 & HD22 H-atoms. . .
                #So, storing OD1 & CG in required lists and doubling ND2 to get zipping enabled with HD21 & HD22. . .
                #------------------------------------------------------------------------------------------------------------
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
                        #gln_ace = line
                        gln_search_string_a = 'CD'+line[3]+line[4]+line[5]
                        for item in gln_inp_file:
                            if item[2] == 'CD':
                                gln_target_string = item[2]+item[3]+item[4]+item[5]
                                if gln_search_string_a == gln_target_string:
                                    adjacent_gln.append(item)  
                                    acceptor_gln.append(line)
                    #***********
                    #if line[2] == 'OE1':
                        #acceptor_gln.append(line)
                    #elif line[2] == 'CD':
                        #adjacent_gln.append(line)
                    #**********
                    elif line[2] == 'NE2':
                        donor_gln.append(line)
                    elif line[2] == 'HE21' or line[2] == 'HE22' or line[2] == 'HE2':
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

                #----------------------------------------------------------------------------------------------
                #including SERINE(SER) residue . . . [DONOR & ACCEPTOR]
                #it's OG acts as both acceptor and donor atom. . .
                #So, alloting OG as donor & HG as H-atom and OG as acceptor & CB as adjacent atom. . .
                #-----------------------------------------------------------------------------------------------
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
                        donor_ser.append(line)
                    elif line[2] == 'HG':
                        H_ser.append(line)

                #Sometimes adjacent atom goes missing, so creating such loops helps maintain the whole data string hole-free.        
                for line in ser_inp_file:
                    if line[2] == 'OG':
                        ser_ace = line
                        ser_search_str_a = 'CB'+line[3]+line[4]+line[5]
                        for item in ser_inp_file:
                            if item[2] == 'CB':
                                ser_search_str_b = item[2]+item[3]+item[4]+item[5]
                                if ser_search_str_a == ser_search_str_b:
                                    adjacent_ser.append(item)
                                    acceptor_ser.append(ser_ace)


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


                #---------------------------------------------------------------------------------------------------
                #including THREONINE(THR) residue. . . [DONOR & ACCEPTOR]
                #it's OG1 atom acts as both donor and acceptor. . .
                #So, alocating OG1 in donor & HG1 in H-atom list and OG1 in acceptor & CB in adjacent list. . .
                #----------------------------------------------------------------------------------------------------
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
                        donor_thr.append(line)
                    elif line[2] == 'HG1':
                        H_thr.append(line)

                #Sometimes adjacent atom goes missing, so creating such loops helps maintain the whole data string hole-free.        
                for line in thr_inp_file:
                    if line[2] == 'OG1':
                        thr_ace = line
                        thr_search_str_a = 'CB'+line[3]+line[4]+line[5]
                        for item in thr_inp_file:
                            if item[2] == 'CB':
                                thr_search_str_b = item[2]+item[3]+item[4]+item[5]
                                if thr_search_str_a == thr_search_str_b:
                                    adjacent_thr.append(item)
                                    acceptor_thr.append(thr_ace)


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


                #---------------------------------------------------------------------------------------------------
                #including TYROSINE(TYR) residue. . . [DONOR & ACCEPTOR]
                #its OH atom acts as both donor and acceptor. . .
                #So, allocating OH in donor, HH in H-atom and OH in acceptor, CZ in adjacent atom list. . .
                #---------------------------------------------------------------------------------------------------
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
                        donor_tyr.append(line)
                    elif line[2] == 'HH':
                        H_tyr.append(line)

                #Sometimes adjacent atom goes missing, so creating such loops helps maintain the whole data string hole-free.        
                for line in tyr_inp_file:
                    if line[2] == 'OH':
                        tyr_ace = line
                        tyr_search_str_a = 'CZ'+line[3]+line[4]+line[5]
                        for item in tyr_inp_file:
                            if item[2] == 'CZ':
                                tyr_search_str_b = item[2]+item[3]+item[4]+item[5]
                                if tyr_search_str_a == tyr_search_str_b:
                                    adjacent_tyr.append(item)
                                    acceptor_tyr.append(tyr_ace)

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


                #----------------------------------------------------------------------------------------------------
                #including TRYPTOPHAN(TRP) residue . . . [DONOR]
                #it's NE1 atom as donor and HE1 is corresponding H-atom. . .
                #So, allocating NE1 in donor & HE1 in H-atom list. . .
                #----------------------------------------------------------------------------------------------------
                donor_trp = []
                H_trp = []
                trp_inp_file = []

                for line in map(str.split, inp):
                    residueName = line[3][-3:]    ##just taking last 3 alphabets so to avoid 'ATRP' instead of 'TRP'...
                    if residueName == 'TRP':
                        trp_inp_file.append(line)

                for line in trp_inp_file:
                    if line[2] == 'NE1':
                        donor_trp.append(line)
                    elif line[2] == 'HE1':
                        H_trp.append(line)

                #for don,h in zip(donor_trp, H_trp):
                #    print don,h
                #categorization of TYPTOPHAN as Donor completed . . . >>[NOTE]>>> put it in subsequent list for calculations

                #INCLUDING IN MOTHER LIST . . . > > > > > >
                #print donor_trp, H_trp
                #print "----" *20
                side_donor.extend(donor_trp)
                side_H.extend(H_trp)
                #for don, h,ace,adj in zip(side_donor, side_H,side_acceptor, side_adjacent): print don, h, ace, adj

                #--------------------------------------------------------------------------
                #>>>> EXTENDING MAIN_CHAIN_COORDINATES WITH SIDE_CHAIN_COORDINATES... >>>>
                donor.extend(side_donor)
                H.extend(side_H)
                acceptor.extend(side_acceptor)
                adjacent.extend(side_adjacent)

                #================================================================
                #PREPARING FOR CALCULATION NOW. . .
                #......... :)              :)           :) ..................... >>>>>
                #================================================================
                for don,h in  zip(donor,H):
                    #print don, h
                    #print "---"*20

                    '''
                    #--------------------------------------------------------------------------
                    #Constructing a Progressbar to minor loop progress...                    
                    toolbar_width = len(donor)
                    sys.stdout.write("[%s]" % (" "*toolbar_width))
                    sys.stdout.flush()
                    sys.stdout.write("\b"*(toolbar_width+1)) #return to the startline after '['
                    #----------------------------------------------------------------------------
                    '''
                    #print don
                    don_xAxis = don[6]
                    don_yAxis = don[7]
                    don_zAxis = don[8]
                    h_xAxis = h[6]
                    h_yAxis = h[7]
                    h_zAxis = h[8]
                    da_cal_out.write('\n')
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
                        #print ace, adj
                        ace_xAxis = ace[6]
                        ace_yAxis = ace[7]
                        ace_zAxis = ace[8]
                        xy = ('acceptor: '+ace[2]),ace[3],ace[4],ace[5],ace[6],ace[7],ace[8]
                        DA = ' '.join(map(''.join,xy))
                        adj = ('adjacent : '+adj[2]),adj[3],adj[4],adj[5],adj[6],adj[7],adj[8]
                        adjFinal = ' '.join(map(''.join,adj))
                        #-----------------------------------------------------------------
                        #Generating file in a format to ease angle calculation(*.DA_cal)...
                        #-----------------------------------------------------------------
                        da_cal_out.write(final_donorAtom) #getting donor atom
                        da_cal_out.write(' ')
                        da_cal_out.write(final_hAtom)     #getting H-atom
                        da_cal_out.write(' ')
                        da_cal_out.write(DA)              #getting acceptor atom
                        da_cal_out.write(' ')
                        da_cal_out.write(adjFinal)        #getting acceptor-adjacent atom
                        da_cal_out.write('\n')
                        #---------------------------------------------------------------
                        #now generating same file but in READABLE format(.DA)...
                        #---------------------------------------------------------------
                        da_out.write('\n')
                        da_out.write(final_donorAtom) #getting donor atom
                        da_out.write('\n')
                        da_out.write(final_hAtom)     #getting H-atom
                        da_out.write('\n')
                        da_out.write(DA)              #getting acceptor atom
                        da_out.write('\n')
                        da_out.write(adjFinal)        #getting acceptor-adjacent atom
                        da_out.write('\n')

                loop_clock_stop = timeit.default_timer() #Terminating loop timer $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                loop_time = str(loop_clock_stop - loop_clock_start)



                da_cal_out.close()
                da_out.close()
                print "  - Primary Scan Complete, .DA_cal and .DA outfile created for " + filename 

        listOfOpenedFilesOut = open(path+"listOfOpenedFiles.db", 'w')
        for item in listOfOpenedFiles: 
            listOfOpenedFilesOut.write(item)
            listOfOpenedFilesOut.write('\n')
        listOfOpenedFilesOut.close()

    print "  -- .DA_cal & .DA file generated for all .atom files..."
    print "  -- Searching dir for .DA_cal files to generate .hbinfo file..."
    #------------------------------------------------------------------------------
    #Now as the calculation ready format .DA_cal has been created.
    #We gonna search th directory for .DA-cal and process each onw by one...
    #-------------------------------------------------------------------------------

    #Now have to open .DA_cal files one by one to make further processing...
    #inp2 = open("/home/saumen/abhisek/dft_cal/xpt_zone_3_fine-opt/4g78FH.DA_cal",'r').read().strip().split('\n')
    #OPENING A NEW SEARCH STRING...
    for root, dirs, files in os.walk(path):
        f_count = 0 #initializing count, want to count total no of files...
        list_of_DA_cal_files = glob.iglob(os.path.join(root,'*.DA_cal'))

        for DA_cal_file in list_of_DA_cal_files:
            #distance list for storage . . .
            DA = []
            HA = []
            #angle list for storage...
            DHA_list = []   #contain all calculated DHA angle values.
            HAAj_list = []  #contain all calculated HAAj angle values.
            DAAj_list = []  #contain all calculated DAAj angle values.
            inp2_list = []
            inp2_list_new = []   #list made after removing trailing spaces. . . 

            #Putting a file marker for logical-checkpoint...
            DA_filename = basename(DA_cal_file)


            inp2 = open(DA_cal_file,'r').read().strip().split('\n')
            print "                                                 "
            print "    >> Now Processing file: " + DA_filename

            #----------------------------------------------------
            #Opening .hbinfo file for each opened .DA_cal file..
            hbinfo = DA_filename.replace('DA_cal','hbinfo') 
            hbinfo_out = open(os.path.join(root,hbinfo),'w')


            #Opening .Hbradar file for each open .DA-cal file...
            hbradar = DA_filename.replace('DA_cal','HBradar')
            hbradar_out = open(os.path.join(root,hbradar),'w')   

            #-----------------------------------------------------
            #print "[debug-step] blank files created..."    

            #Following list would be used for calculation of surrounding residues...
            #It needed to be oepend earlier to avoid repeated file opening caused by for loops...
            inp4 = []
            inp3 = open(path + "listOfOpenedFiles.db",'r').read().strip().split('\n') #CHANGE PATH as per REQUIRED
            for fileName in map(str.split,inp3):
                strFilename = ''.join(fileName)
                inp4surround = open(path + strFilename,'r').read().strip().split('\n')   #CHANGE PATH as per REQUIRED
                for lines in map(str.split,inp4surround):
                    inp4.append(lines)	    
            #####################################################################################

            for lines in map(str.split,inp2):
                #print lines
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
                if donor_atom == ace_atom and donor_res_name == ace_res_name and donor_chain == ace_chain and donor_res_no == ace_res_no and donor_xAxis == ace_xAxis:
                    break
                else:
                    #print lines
                    #--------------------------------------------------------------------------
                    #>>> STARTING DISTANCE CALCULATIONS. . .
                    #Distance_1: DA . . .
                    #--------------------------------------------------------------------------
                    xD = (ace_xAxis) - (donor_xAxis)
                    yD = (ace_yAxis) - (donor_yAxis)
                    zD = (ace_zAxis) - (donor_zAxis)
                    DA_distance = math.sqrt(xD*xD + yD*yD + zD*zD) #distance b/w donor-acceptor.
                    DA.append(DA_distance)

                    #-------------------------------------------------------------------------
                    #Distance_2: HA. . .
                    #-------------------------------------------------------------------------
                    xD1 = (ace_xAxis) - (hy_xAxis)
                    yD1 = (ace_yAxis) - (hy_yAxis)
                    zD1 = (ace_zAxis) - (hy_zAxis)
                    HA_distance = math.sqrt(xD1*xD1 +yD1*yD1 + zD1*zD1)
                    HA.append(HA_distance)

                    #--------------------------------------------------------------------------
                    #>>> STARTING ANGLE CALCULATIONS. . . 
                    #Angle_1: DHA...
                    #Calculating magnitude of the two vectors HD & HA i.e. DHA ANGLE . . .!!!
                    #-------------------------------------------------------------------------
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
                    denominator = ((math.sqrt(math.pow(x1,2)+math.pow(y1,2)+math.pow(z1,2)))*(math.sqrt(math.pow(x2,2)+math.pow(y2,2)+math.pow(z2,2))))
                    cosine = numerator/denominator
                    theta_rad = math.acos(cosine) #inverse cos function gives angle measure in radian...
                    DHA_angle = math.degrees(theta_rad) #converted to degree...
                    DHA_list.append(DHA_angle)
                    #print 'DHA = '+ str(DHA_angle)


                    #--------------------------------------------------------------------------
                    #Angle_2: HAAj...(Aj = adjacent atom to Acceptor atom)
                    #Calculating magnitude of the two vectors AH & Aj` i.e. HAAj ANGLE . . .!!!
                    #--------------------------------------------------------------------------
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

                    #---------------------------------------------------------------------------
                    #Angle_3: DAAj...(Aj = adjacent atom to Acceptor atom)
                    #Calculating magnitude of the two vectors AD & Aj` i.e. DAAj ANGLE . . .!!!
                    #---------------------------------------------------------------------------
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

            #Now Preparing for DFT calculations...
            #.HBradar file already been opened above as 'hbradar_out'...
            #--------------------------------------------------------------------
            print "        Generating_outfile: "+ hbinfo
            final_list = []
            dft_input_list = []    #CONTAINS DATA FOR NWCHEM INPUT.. >>>>>>

            h1 = str("_"*200)   #adding a header line for viewing pleasure...
            hbinfo_out.write(h1)
            hbinfo_out.write('\n')
            header_line = str( "DONOR"+" "*20+"HYDROGEN-ATOM"+" "*12+"ACCEPTOR ATOM"+" "*12+"ADJACENT ATOM"+" "*10+"DA"+" "*6+"HA"+" "*7+"DHA"+" "*7+"HAAj"+" "*9+"DAAj"+" "*9+"Cluster(m06-2x)"+" "*14+"DummyClusterE"+" "*10+"DiffEnergy(Cluster-dummyCluster)")   #headers line. . .

            hbinfo_out.write(header_line)
            hbradar_out.write(header_line)
            hbinfo_out.write('\n')
            hbradar_out.write('\n')
            hbinfo_out.write(h1)
            hbradar_out.write(h1)
            hbinfo_out.write('\n')
            hbradar_out.write('\n')
            #print "[debug-step]header written..."

            radial_distance_list=[]
            midpoint_list = [] 
            HB_coord_list=[] #Contains all d coordiantes of HB partners.
            HB_output_list=[] #to store output items forming hydrogen bond, in a particular format, for later ease of usage...
            for lines,da,ha,items,case,entry in zip(inp2_list_new,DA,HA,DHA_list,HAAj_list,DAAj_list):
                #print lines,da,ha,items,case,entry
                if da <= 3.9 and ha <= 2.5 and items > 90 and case >90 and entry > 90:
                    #print lines
                    HB_coord_list.append(lines)
                    #print lines, ('DA: '+ str(da)), ('HA: '+ str(ha)), ('DHA: ' + str(items)), ('HAAj: ' + str(case)), ('DAAj: ' + str(entry))
                    #output = lines[1], lines[2], lines[3],lines[4],'    ', lines[9],lines[10],lines[11],lines[12],'    ',lines[17],lines[18],lines[19],lines[20],'    ',lines[26],lines[27],lines[28],lines[29],'    ', ('DA: %.4f ' % da),'    ',('HA: %.4f' % ha),'    ', ('DHA: %.4f' % items),'    ', ('HAAj: %.4f' % case),'    ', ('DAAj: %.4f' % entry)       

                    donor_lines =str(str(lines[1])+" "*(6-len(lines[1]))+str(lines[2])+" "*(6-len(lines[2]))+
                                     str(lines[3])+" "*(4-len(lines[3]))+str(lines[4]))
                    dft_donor_lines = str(str(lines[1][0])+' '+str(lines[5])+' '+str(lines[6])+' '+
                                          str(lines[7]))    #to be used as DFT input data for DONOR_ATOM...
                    #donor_lines_list.append(donor_lines)
                    hydrogen_lines = str(" "*(8-len(lines[4]))+str(lines[9])+" "*(6-len(lines[9]))+str(lines[10])+
                                         " "*(6-len(lines[10]))+str(lines[11])+" "*(4-len(lines[11]))+lines[12])
                    dft_hydrogen_lines = str(str(lines[9][0])+' '+str(lines[13])+' '+str(lines[14])+' '+
                                             str(lines[15])) #to be used as DFT input data for HYDROGEN_ATOM...
                    #hydrogen_lines_list.append(hydrogen_lines)
                    acceptor_lines = str(" "*(8-len(lines[12]))+str(lines[17])+" "*(6-len(lines[17]))+
                                         str(lines[18])+" "*(6-len(lines[18]))+str(lines[19])+" "*(4-len(lines[19]))+str(lines[20]))
                    dft_acceptor_lines = str(str(lines[17][0])+' '+str(lines[21])+' '+str(lines[22])+' '+
                                             str(lines[23]))  #to be used as DFT input data for ACCEPTOR_ATOM...
                    #acceptor_lines_list.append(acceptor_lines)
                    adjacent_lines = str(" "*(8-len(lines[20]))+str(lines[26])+" "*(6-len(lines[26]))+
                                         str(lines[27])+" "*(6-len(lines[27]))+str(lines[28])+" "*(4-len(lines[28]))+str(lines[29]))
                    dft_adjacent_lines = str(str(lines[26][0])+' '+str(lines[30])+ ' '+str(lines[31])+' '+
                                             str(lines[32]))  #to be used as DFT input data for ADJACENT_ATOMS...

                    total_dft_input = dft_donor_lines+ ';' +dft_hydrogen_lines+';'+dft_acceptor_lines+';'+dft_adjacent_lines
                    dft_input_list.append(total_dft_input)         #storing in the list for mass usage later on...
                    #adjacent_lines_list.append(adjacent_lines)

                    da_distance = str(" "*(6-len(lines[29]))+'DA:%.2f' % (da))
                    ha_distance = str(" "*(6-len(da_distance))+'HA:%.2f' % (ha))
                    dha_angle =  str(" "*(6-len(ha_distance))+'DHA:%.2f' % (items))
                    haaj_angle = str(" "*(6-len(dha_angle)) +'HAAj:%.2f' % (case))
                    daaj_angle = str(" "*(11-len(haaj_angle))+'DAAj:%.2f' % (entry))

                    zz =  donor_lines, hydrogen_lines, acceptor_lines, adjacent_lines,da_distance,ha_distance,dha_angle,haaj_angle,daaj_angle
                    #string_zz = str(zz)  #converting zz-tuple to string
                    zzz = ' '.join(zz)
                    HB_output_list.append(zzz)
                    hbinfo_out.write(zzz)
                    hbinfo_out.write('\n')
                    #print zzz

                    #formatted_output = (' '.join((output)))
                    #z = str(formatted_output)
                    #final_list.append(z)

            for line in HB_coord_list:
                midpoint = (float(line[13])+float(line[21]))/2,(float(line[14])+float(line[22]))/2,\
                    (float(line[15])+float(line[23]))/2
                midpoint_list.append(midpoint)

            final_list=[]    #CONTAINS ALL THE COORDINATES FOR DFT/MP2 CALCULATIONs...
            ghost_final_list=[] #Contains all the coordinate info but hydrogen bond forming atoms are ghost here
            atomCountInEachCluster=[]

            for midpoint,hb_coord in zip(midpoint_list,HB_coord_list):
                print hb_coord #Prints coord info about the concerned H bond

                #The Donor, Hydrogen, Acceptor, Adjacent coordinates are here......
                x = [[hb_coord[1][0],hb_coord[5],hb_coord[6],hb_coord[7]],';',[hb_coord[9][0],hb_coord[13],hb_coord[14],\
                                                                               hb_coord[15]],';',[hb_coord[17][0],hb_coord[21],hb_coord[22],hb_coord[23]],';',\
                     [hb_coord[26][0],hb_coord[30],hb_coord[31],hb_coord[32],';']]

                #Adding 'bq' to the atom names, which to be declared as dummy during DFT calculations...
                firstAtom = 'bq'+hb_coord[1][0]
                secondAtom = 'bq'+hb_coord[9][0]
                thirdAtom = 'bq'+hb_coord[17][0]
                fourthAtom = 'bq'+hb_coord[26][0]

                #Group of atoms which will be appended as dummy atom in the dft file
                #ghost_x = [[firstAtom,hb_coord[5],hb_coord[6],hb_coord[7]],';',[secondAtom,hb_coord[13],hb_coord[14],\
                #           hb_coord[15]],';',[thirdAtom,hb_coord[21],hb_coord[22],hb_coord[23]],';',\
                #           [fourthAtom,hb_coord[30],hb_coord[31],hb_coord[32],';']]

                ghost_x = [[hb_coord[1][0],hb_coord[5],hb_coord[6],hb_coord[7]],';',[hb_coord[9][0],hb_coord[13],hb_coord[14],\
                                                                                     hb_coord[15]],';', [hb_coord[26][0],hb_coord[30],hb_coord[31],hb_coord[32],';']]                     #Deleted 3rd atom to break the H bond...


                #print x
                #print ghost_x
                #Inp4 was opened earlier for better calling purpose...

                #print lines
                #print midpoint[0], lines[6]
                for lines in inp4:
                    x_dist = float(midpoint[0]) - float(lines[6])
                    y_dist = float(midpoint[1]) - float(lines[7])
                    z_dist = float(midpoint[2]) - float(lines[8])
                    radial_distance = math.sqrt((x_dist**2)+(y_dist**2)+(z_dist**2))
                    #print radial_distance
                    #DETERMINING Search Radius...
                    if radial_distance <= 3.0:
                        #print hb_coord[4], lines
                        f = [[lines[2][0],lines[6], lines[7],lines[8],';']]
                        #print f
                        if lines[6] != hb_coord[5] and lines[6] != hb_coord[13] and lines[6] != hb_coord[21] and lines[6] != hb_coord[30]:
                            x.extend(f)
                            ghost_x.extend(f) #Surrounding residues will be the same
                            print "atoms found  in radial search" + str(lines) + str(radial_distance)
                atomCountInEachCluster.append(len(x)-3)
                #print x
                print len(x)-3 #Prints number of atoms in each cluster... -3 means, three ';' are deducted from len counting...
                                #x.extend(lines)

                final_list.append(x)
                ghost_final_list.append(ghost_x)
            #print final_list
            #for line in final_list: print (line)
            #for item in atomCountInEachCluster: print item

            hbinfo_out.close()

            print "        --.hbinfo file generated successfully      " 
            #""" Following part is active in HPC version, laptops wont be able to run this below part, so inactivated[no bu]
            #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$       
            #---------------------------------------------------------------------------------------------------------------
            #       ...............PREPARING FOR DFT INPUT CONFIGURATION...............................
            #               <<<data are stored in list named "dft_input_list">>>>
            #================================================================================================================
            print "  -- Starting the DFT calculations now"
            print "                      please wait, this may take a while. . ."
            import sys
            import subprocess
            m06_2x = [] #list to store dft energy coming from m06-2x function...
            m05_2x = []
            m06_HF = []
            ghost_m06_2x = [] #List to store DFT energy coming from calculations with ghost/dummy atoms...


            #--------------------------------------------------------------------------------------------------------------
            #CREATING FUNCTION TO ANALYZE DFT_code i.e. <analyze_dft_code>
            #_________________________________________________________________________________________________________________
            def analyze_dft_code(dft_code_exe,xc,func_name): #xc2 is for dummy atom list creation only...
                dft_out = open(path + "e_hydroden_bond.data",'w')
                print 'starting subprocess....'
                p =subprocess.Popen(['mpirun --mca btl self,sm,openib -np 64 -hostfile int_node /app/NWChem/bin/nwchem ' + dft_code_exe], stdout=subprocess.PIPE, shell=True)
                #Use the following line if wish to run in a SINGLE NODE....
                #p =subprocess.Popen(['mpirun -np 16 -hostfile $PBS_NODEFILE nwchem ' + dft_code_exe], stdout=subprocess.PIPE, shell=True)
                output, err = p.communicate()
                print "communication successful with nwchem"
                for line in output:
                    dft_out.write(line)
                dft_out.close()
                dft_file = open(path + "e_hydroden_bond.data",'r').read().strip().split('\n')
                total_lines = []     #storing all nwchem generated DFT & BSSE energies...
                useful_energies = []  #filtering from list "total_lines" and storing only 'isolated' and 'BSSE corrected energies'...

                #Loop to point out energy terms in DFT_script...
                for lines in map(str.split,dft_file):
                    if lines != []:
                        if lines[0] == 'Total' and lines[1] == 'DFT' and lines[2] == 'energy':
                            total_lines.append(lines)
                        else:
                            if lines[0] == 'Corrected' and lines[1] == 'energy':
                                total_lines.append(lines)			


                try:
                    isolated_energy1 = float(total_lines[1][-1])
                    isolated_energy2 = float(total_lines[3][-1])
                    BSSE_energy = float(total_lines[5][-1])
                    print isolated_energy1, isolated_energy2, BSSE_energy
                    net_energy = (BSSE_energy - (isolated_energy1 + isolated_energy2))  # generates energy in Hartree...
                    net_energy_kcal = net_energy * 627.509 #Converting energy into kcal/mol unit...
                    print ('                     DFT Energy(%s) = %f'  %(func_name, net_energy_kcal) +' kcal/mol')
                    xc.append(net_energy_kcal)
                    #print "COMMUNICATION COMPLETE WITH NWCHEM >>>>>>>>>>>>>>>>>>>>>>"

                except IndexError:
                    print ('                     DFT Energy = ') + 'Failed to calculate for this coordinate...'
                    xc.append('00.00000')
                return
                ##############################################################################################################


            #--------------------------------------------------------------------------------------------------------------
            #CREATING FUNCTION TO ANALYZE GHOST/DUMMY - DFT_code i.e. <analyze_dft_code>
            #_________________________________________________________________________________________________________________		
            def analyze_GHOST_dft_code(ghost_dft_code_exe,ghost_xc,func_name): #xc2 is for dummy atom list creation only...
                dft_out = open(path + "e_hydroden_bond.data",'w')
                p =subprocess.Popen(['mpirun --mca btl self,sm,openib -np 64 -hostfile int_node /app/NWChem/bin/nwchem ' + ghost_dft_code_exe], stdout=subprocess.PIPE, shell=True)
                #Use the following line if wish to run in a SINGLE NODE....
                #p =subprocess.Popen(['mpirun -np 16 -hostfile $PBS_NODEFILE nwchem ' + ghost_dft_code_exe], stdout=subprocess.PIPE, shell=True)
                output, err = p.communicate()
                print "communication successful with nwchem"
                for line in output:
                    dft_out.write(line)
                dft_out.close()
                dft_file = open(path + "e_hydroden_bond.data",'r').read().strip().split('\n')
                total_lines = []     #storing all nwchem generated DFT & BSSE energies...
                useful_energies = []  #filtering from list "total_lines" and storing only 'isolated' and 'BSSE corrected energies'...

                #Loop to point out energy terms in DFT_script...
                for lines in map(str.split,dft_file):
                    if lines != []:
                        if lines[0] == 'Total' and lines[1] == 'DFT' and lines[2] == 'energy':
                            total_lines.append(lines)
                        else:
                            if lines[0] == 'Corrected' and lines[1] == 'energy':
                                total_lines.append(lines)			

                try:
                    isolated_energy1 = float(total_lines[1][-1])
                    isolated_energy2 = float(total_lines[3][-1])
                    BSSE_energy = float(total_lines[5][-1])
                    #print isolated_energy1, isolated_energy2, BSSE_energy
                    net_energy = (BSSE_energy - (isolated_energy1 + isolated_energy2))  # generates energy in Hartree...
                    net_energy_kcal = net_energy * 627.509 #Converting energy into kcal/mol unit...
                    print ('                     DFT Energy(%s) = %f'  %(func_name, net_energy_kcal) +' kcal/mol')
                    ghost_xc.append(net_energy_kcal)
                    #print "COMMUNICATION COMPLETE WITH NWCHEM >>>>>>>>>>>>>>>>>>>>>>"
                except IndexError:
                    print ('                     DFT Energy = ') + 'Failed to calculate for this coordinate...'
                    ghost_xc.append('00.00000')
                return		


            #_________________________________________________________________________________________________________________



            import itertools


            #print final_list[-1][-1][-1] = ; #for info
            for items in final_list:
                #print items
                del items[-1][-1]
                electron_num = []  #Electron number of non-dummy atom group
                bsse_mon_2 = []
                #print items
                for i in items:
                    if i[0] =='N':
                        electron_num.append(7)
                    elif i[0]=='O':
                        electron_num.append(8)
                    elif i[0]=='C':
                        electron_num.append(6)
                    elif i[0] == 'H':
                        electron_num.append(1)



                total_electron_number = sum(electron_num)
                for i in range (3, len(electron_num)+1): bsse_mon_2.append(str(i))


                dft_string = list(itertools.chain.from_iterable(items))
                final_dft_string = ' '.join(dft_string) #to be appended in dft calculations
                final_bsse_mon_2 = ' '.join(bsse_mon_2) #to be appended as bsse mon 2
                #print final_dft_string
                #print final_bsse_mon_2


                '''
                 - calculating electron number based on atoms involved...
                 - As because, in case of ODD number of electrons of participating atoms, the calculation must go in ODFT-format
                 - But for EVEN no. of electrons, preferred is simple DFT-protocol.
                 -- Now First checking if the no. of electron is ODD/EVEN, then approaching further likewise...

                '''


                if (total_electron_number % 2) == 0:   #i.e. number of electron is even, calculation will go on with 'dft protocol'...
                    #print items

                    '''
                    ** New to v.3.1
                    ________________
                    Writing an EXCEPTION handler (in case the convergence fails in 6-31g(medium grid), DFT func gonna try 
                    another shot of calculation in 6-311g(fine-grid) and it it fails too then DFT gonna switch to
                    6-311g(coarse grid) template)
                    > So, the exception handler order is:
                          6-31g(medium grid) > 6-311g (fine grid) > 6-311g (coarse grid)
                    '''

                    try:
                        f_m06_2x = open(path + "dft_m06-2x_631g_template.nw",'r')
                        filedata=f_m06_2x.read()
                        f_m06_2x.close()
                        newdata = filedata.replace("  N##",("  "+final_dft_string)).replace("  #mon#","  mon second "+final_bsse_mon_2)

                        dft_code = open(path + "e_hydrogen_bond_exe.nw",'w')
                        dft_code.write(newdata)
                        dft_code.close()
                        #Now calling dft_analyzer function created earlier...
                        analyze_dft_code(dft_code_exe=path + "e_hydrogen_bond_exe.nw", xc=m06_2x, func_name='m06-2x')

                    except IndexError:   #could not converge in xfine so going for xfine now... 
                        #Going for 6-311g(fine-grid) exception ...
                        try:
                            print "Seems problem in Convergence, switching to 6-311g basis-set(grid: fine)..."
                            f_m06_2x = open(path + "dft_m06-2x_6311g-fine_template.nw",'r')
                            filedata=f_m06_2x.read()
                            f_m06_2x.close()
                            newdata = filedata.replace("  N##",("  "+final_dft_string)).replace("  #mon#","  mon second "+final_bsse_mon_2)

                            dft_code = open(path + "e_hydrogen_bond_exe.nw",'w')
                            dft_code.write(newdata)
                            dft_code.close()
                            #Now calling dft_analyzer function created earlier...
                            analyze_dft_code(dft_code_exe= path + "e_hydrogen_bond_exe.nw", xc=m06_2x, func_name='m06-2x')
                        #Going for 6-311g(Coarse-grid) exception ...
                        except IndexError:
                            try:
                                print "Seems a problem with convergence, switching to 6-311g basis set(grid: coarse)..."
                                f_m06_2x = open(path + "dft_m06-2x_6311g-coarse_template.nw",'r')
                                filedata=f_m06_2x.read()
                                f_m06_2x.close()
                                newdata = filedata.replace("  N##",("  "+final_dft_string)).replace("  #mon#","  mon second "+final_bsse_mon_2)

                                dft_code = open(path + "e_hydrogen_bond_exe.nw",'w')
                                dft_code.write(newdata)
                                dft_code.close()
                                #Now calling dft_analyzer function created earlier...
                                analyze_dft_code(dft_code_exe=path + "e_hydrogen_bond_exe.nw", xc=m06_2x, func_name='m06-2x')
                            except IndexError:
                                print "Seems it failed again ! Switching to 6-311g (grid: xcoarse)..."
                                f_m06_2x = open(path + "dft_m06-2x_6311g-xcoarse_template.nw",'r')
                                filedata=f_m06_2x.read()
                                f_m06_2x.close()
                                newdata = filedata.replace("  N##",("  "+final_dft_string)).replace("  #mon#","  mon second "+final_bsse_mon_2)

                                dft_code = open(path + "e_hydrogen_bond_exe.nw",'w')
                                dft_code.write(newdata)
                                dft_code.close()
                                #Now calling dft_analyzer function created earlier...
                                analyze_dft_code(dft_code_exe=path + "e_hydrogen_bond_exe.nw", xc=m06_2x, func_name='m06-2x')



                elif (total_electron_number % 2) > 0:  #i.e. no. of electron is ODD, Switching to ODFT-protocols...
                    #print items
                    #Going into the exception looping...
                    try:
                        f_m06_2x = open(path + "odft_m06-2x_631g_template.nw",'r')
                        filedata=f_m06_2x.read()
                        f_m06_2x.close()
                        #> Now if atom 3 is N them we might get multiplicity error, so if-else loop required...
                        if items[0][0] != 'N':
                            newdata = filedata.replace("  O##",("  "+final_dft_string)).replace("  #mult@@","  mult 2").replace("  #mon#","  mon second "+final_bsse_mon_2)
                        elif items[0][0] == 'N' and (sum(electron_num[2:]) % 2) > 0:
                            newdata = filedata.replace("  O##",("  "+final_dft_string)).replace("  #mult$$","  mult 2").replace("  #mon#","  mon second "+final_bsse_mon_2)

                        dft_code_exe = open(path + "e_hydrogen_bond_exe.nw",'w')
                        dft_code_exe.write(newdata)
                        dft_code_exe.close()
                        #Now calling dft_analyzer function created earlier...
                        analyze_dft_code(dft_code_exe=path + "e_hydrogen_bond_exe.nw", xc=m06_2x, func_name='m06-2x')
                    except IndexError:
                        #Going for 6-311g(fine-grid) exception ...
                        try:
                            print "Seems a problem with convergence, switching to 6-311g (grid: fine)..."
                            f_m06_2x = open(path + "odft_m06-2x_6311g-fine_template.nw",'r')
                            filedata=f_m06_2x.read()
                            f_m06_2x.close()
                            #> Now if atom 3 is N them we might get multiplicity error, so if-else loop required...
                            if items[0][0] != 'N':
                                newdata = filedata.replace("  O##",("  "+final_dft_string)).replace("  #mult@@","  mult 2").replace("  #mon#","  mon second "+final_bsse_mon_2)
                            elif items[4][0] == 'N' and (sum(electron_num[2:]) % 2) > 0:
                                newdata = filedata.replace("  O##",("  "+final_dft_string)).replace("  #mult$$","  mult 2").replace("  #mon#","  mon second "+final_bsse_mon_2)

                            dft_code_exe = open(path + "e_hydrogen_bond_exe.nw",'w')
                            dft_code_exe.write(newdata)
                            dft_code_exe.close()
                            #Now calling dft_analyzer function created earlier...
                            analyze_dft_code(dft_code_exe=path + "e_hydrogen_bond_exe.nw", xc=m06_2x, func_name='m06-2x')
                        #Going for 6-311g(coarse-grid) exception ...
                        except IndexError:
                            try:
                                print "Seems the convergence persists, switching to 6-311g (grid: coarse)..."
                                f_m06_2x = open(path + "odft_m06-2x_6311g-coarse_template.nw",'r')
                                filedata=f_m06_2x.read()
                                f_m06_2x.close()
                                #> Now if atom 3 is N them we might get multiplicity error, so if-else loop required...
                                if items[0][0] != 'N':
                                    newdata = filedata.replace("  O##",("  "+final_dft_string)).replace("  #mult@@","  mult 2").replace("  #mon#","  mon second "+final_bsse_mon_2)
                                elif items[4][0] == 'N' and (sum(electron_num[2:]) % 2) > 0:
                                    newdata = filedata.replace("  O##",("  "+final_dft_string)).replace("  #mult$$","  mult 2").replace("  #mon#","  mon second "+final_bsse_mon_2)

                                dft_code_exe = open(path + "e_hydrogen_bond_exe.nw",'w')
                                dft_code_exe.write(newdata)
                                dft_code_exe.close()
                                #Now calling dft_analyzer function created earlier...
                                analyze_dft_code(dft_code_exe=path + "e_hydrogen_bond_exe.nw", xc=m06_2x, func_name='m06-2x')
                            except IndexError:
                                print "oh faied again! switching to 6-311g (grid: xcoarse)..."
                                f_m06_2x = open(path + "odft_m06-2x_6311g-xcoarse_template.nw",'r')
                                filedata=f_m06_2x.read()
                                f_m06_2x.close()
                                #> Now if atom 3 is N them we might get multiplicity error, so if-else loop required...
                                if items[0][0] != 'N':
                                    newdata = filedata.replace("  O##",("  "+final_dft_string)).replace("  #mult@@","  mult 2").replace("  #mon#","  mon second "+final_bsse_mon_2)
                                elif items[4][0] == 'N' and (sum(electron_num[2:]) % 2) > 0:
                                    newdata = filedata.replace("  O##",("  "+final_dft_string)).replace("  #mult$$","  mult 2").replace("  #mon#","  mon second "+final_bsse_mon_2)

                                dft_code_exe = open(path + "e_hydrogen_bond_exe.nw",'w')
                                dft_code_exe.write(newdata)
                                dft_code_exe.close()
                                #Now calling dft_analyzer function created earlier...
                                analyze_dft_code(dft_code_exe=path + "e_hydrogen_bond_exe.nw", xc=m06_2x, func_name='m06-2x')





            #CHANGING THE STORED ENERGY VALUES INTO string...    
            m06_2x_string_dft_energy=[]   #converting the energy value in string and storing here for extending the HB-output_list later...
            #m05_2x_string_dft_energy=[]
            #m06_HF_string_dft_energy=[]
            for m06_2x_energy in m06_2x:
                m06_2x_string_value=str(m06_2x_energy)
                m06_2x_string_dft_energy.append(m06_2x_string_value)




            ############################################################################################################################
            #Now, starting calculation for GHOST/DUMMY atoms....
            ############################################################################################################################
            print "Starting caluclation with GHOST/DUMMY ATOMS..."
            for case in ghost_final_list:
                #del case[-1][-1]   #Not required in this case....
                ghost_electron_num = [] #Electon number of dummy atom group
                bsse_ghost = [] #bsse related to ghost atom lists...
                for j in case:
                    if j[0] =='N':
                        ghost_electron_num.append(7)
                    elif j[0]=='O':
                        ghost_electron_num.append(8)
                    elif j[0]=='C':
                        ghost_electron_num.append(6)
                    elif j[0] =='H':
                        ghost_electron_num.append(1)


                total_ghost_electron_number = sum(ghost_electron_num)
                for i in range (3, len(ghost_electron_num)+1): bsse_ghost.append(str(i))

                #For ghost...
                ghost_dft_string = list(itertools.chain.from_iterable(case))
                ghost_final_dft_string = ' '.join(ghost_dft_string) #to be appended in GHOST dft calculations
                ghost_final_bsse_ghost = ' '.join(bsse_ghost) #to be appended as bsse mon 2  [bsse_mon_2 is bsse_ghost, in case of ghost atomic calculations...]


                '''
                 - calculating electron number based on atoms involved...
                 - As because, in case of ODD number of electrons of participating atoms, the calculation must go in ODFT-format
                 - But for EVEN no. of electrons, preferred is simple DFT-protocol.
                 -- Now First checking if the no. of electron is ODD/EVEN, then approaching further likewise...

                '''


                if (total_ghost_electron_number % 2) == 0:   #i.e. number of electron is even, calculation will go on with 'dft protocol'...
                    #print items

                    '''
                    ** New to v.3.1
                    ________________
                    Writing an EXCEPTION handler (in case the convergence fails in 6-31g(medium grid), DFT func gonna try 
                    another shot of calculation in 6-311g(fine-grid) and it it fails too then DFT gonna switch to
                    6-311g(coarse grid) template)
                    > So, the exception handler order is:
                          6-31g(medium grid) > 6-311g (fine grid) > 6-311g (coarse grid)
                    '''

                    try:
                        f_m06_2x = open(path + "dft_m06-2x_631g_template.nw",'r')
                        filedata=f_m06_2x.read()
                        f_m06_2x.close()
                        newdata = filedata.replace("  N##",("  "+ghost_final_dft_string)).replace("  #mon#","  mon second "+ghost_final_bsse_ghost)

                        dft_code = open(path + "e_hydrogen_bond_exe.nw",'w')
                        dft_code.write(newdata)
                        dft_code.close()
                        #Now calling dft_analyzer function created earlier...
                        analyze_GHOST_dft_code(ghost_dft_code_exe=path + "e_hydrogen_bond_exe.nw", ghost_xc=ghost_m06_2x,func_name='m06-2x')

                    except IndexError:   #could not converge in xfine so going for xfine now... 
                        #Going for 6-311g(fine-grid) exception ...
                        try:
                            print "Seems problem in Convergence, switching to 6-311g basis-set(grid: fine)..."
                            f_m06_2x = open(path + "dft_m06-2x_6311g-fine_template.nw",'r')
                            filedata=f_m06_2x.read()
                            f_m06_2x.close()
                            newdata = filedata.replace("  N##",("  "+ghost_final_dft_string)).replace("  #mon#","  mon second "+ghost_final_bsse_ghost)

                            dft_code = open(path + "e_hydrogen_bond_exe.nw",'w')
                            dft_code.write(newdata)
                            dft_code.close()
                            #Now calling dft_analyzer function created earlier...
                            analyze_GHOST_dft_code(ghost_dft_code_exe= path + "e_hydrogen_bond_exe.nw", ghost_xc=ghost_m06_2x,func_name='m06-2x')
                        #Going for 6-311g(Coarse-grid) exception ...
                        except IndexError:
                            try:
                                print "Seems a problem with convergence, switching to 6-311g basis set(grid: coarse)..."
                                f_m06_2x = open(path + "dft_m06-2x_6311g-coarse_template.nw",'r')
                                filedata=f_m06_2x.read()
                                f_m06_2x.close()
                                newdata = filedata.replace("  N##",("  "+ghost_final_dft_string)).replace("  #mon#","  mon second "+ghost_final_bsse_ghost)

                                dft_code = open(path + "e_hydrogen_bond_exe.nw",'w')
                                dft_code.write(newdata)
                                dft_code.close()
                                #Now calling dft_analyzer function created earlier...
                                analyze_GHOST_dft_code(ghost_dft_code_exe=path + "e_hydrogen_bond_exe.nw", ghost_xc=ghost_m06_2x,func_name='m06-2x')
                            except IndexError:
                                print "Seems it failed again ! Switching to 6-311g (grid: xcoarse)..."
                                f_m06_2x = open(path + "dft_m06-2x_6311g-xcoarse_template.nw",'r')
                                filedata=f_m06_2x.read()
                                f_m06_2x.close()
                                newdata = filedata.replace("  N##",("  "+ghost_final_dft_string)).replace("  #mon#","  mon second "+ghost_final_bsse_ghost)

                                dft_code = open(path + "e_hydrogen_bond_exe.nw",'w')
                                dft_code.write(newdata)
                                dft_code.close()
                                #Now calling dft_analyzer function created earlier...
                                analyze_GHOST_dft_code(ghost_dft_code_exe=path + "e_hydrogen_bond_exe.nw", ghost_xc=ghost_m06_2x,func_name='m06-2x')



                elif (total_ghost_electron_number % 2) > 0:  #i.e. no. of electron is ODD, Switching to ODFT-protocols...
                    #print case
                    #Going into the exception looping...
                    try:
                        f_m06_2x = open(path + "odft_m06-2x_631g_template.nw",'r')
                        filedata=f_m06_2x.read()
                        f_m06_2x.close()
                        #> Now if atom 3 is N them we might get multiplicity error, so if-else loop required...
                        if case[0][0] != 'N':
                            newdata = filedata.replace("  O##",("  "+ghost_final_dft_string)).replace("  #mult@@","  mult 2").replace("  #mon#","  mon second "+ghost_final_bsse_ghost)
                        elif case[0][0] == 'N' and (sum(ghost_electron_num[2:]) % 2) > 0:
                            newdata = filedata.replace("  O##",("  "+ghost_final_dft_string)).replace("  #mult$$","  mult 2").replace("  #mon#","  mon second "+ghost_final_bsse_ghost)

                        dft_code_exe = open(path + "e_hydrogen_bond_exe.nw",'w')
                        dft_code_exe.write(newdata)
                        dft_code_exe.close()
                        #Now calling dft_analyzer function created earlier...
                        analyze_GHOST_dft_code(ghost_dft_code_exe=path + "e_hydrogen_bond_exe.nw", ghost_xc=ghost_m06_2x,func_name='m06-2x')
                    except IndexError:
                        #Going for 6-311g(fine-grid) exception ...
                        try:
                            print "Seems a problem with convergence, switching to 6-311g (grid: fine)..."
                            f_m06_2x = open(path + "odft_m06-2x_6311g-fine_template.nw",'r')
                            filedata=f_m06_2x.read()
                            f_m06_2x.close()
                            #> Now if atom 3 is N them we might get multiplicity error, so if-else loop required...
                            if case[0][0] != 'N':
                                newdata = filedata.replace("  O##",("  "+ghost_final_dft_string)).replace("  #mult@@","  mult 2").replace("  #mon#","  mon second "+ghost_final_bsse_ghost)
                            elif case[4][0] == 'N' and (sum(ghost_electron_num[2:]) % 2) > 0:
                                newdata = filedata.replace("  O##",("  "+ghost_final_dft_string)).replace("  #mult$$","  mult 2").replace("  #mon#","  mon second "+ghost_final_bsse_ghost)

                            dft_code_exe = open(path + "e_hydrogen_bond_exe.nw",'w')
                            dft_code_exe.write(newdata)
                            dft_code_exe.close()
                            #Now calling dft_analyzer function created earlier...
                            analyze_GHOST_dft_code(ghost_dft_code_exe=path + "e_hydrogen_bond_exe.nw", ghost_xc=ghost_m06_2x,func_name='m06-2x')
                        #Going for 6-311g(coarse-grid) exception ...
                        except IndexError:
                            try:
                                print "Seems the convergence persists, switching to 6-311g (grid: coarse)..."
                                f_m06_2x = open(path + "odft_m06-2x_6311g-coarse_template.nw",'r')
                                filedata=f_m06_2x.read()
                                f_m06_2x.close()
                                #> Now if atom 3 is N them we might get multiplicity error, so if-else loop required...
                                if case[0][0] != 'N':
                                    newdata = filedata.replace("  O##",("  "+ghost_final_dft_string)).replace("  #mult@@","  mult 2").replace("  #mon#","  mon second "+ghost_final_bsse_ghost)
                                elif case[4][0] == 'N' and (sum(ghost_electron_num[2:]) % 2) > 0:
                                    newdata = filedata.replace("  O##",("  "+ghost_final_dft_string)).replace("  #mult$$","  mult 2").replace("  #mon#","  mon second "+ghost_final_bsse_ghost)

                                dft_code_exe = open(path + "e_hydrogen_bond_exe.nw",'w')
                                dft_code_exe.write(newdata)
                                dft_code_exe.close()
                                #Now calling dft_analyzer function created earlier...
                                analyze_GHOST_dft_code(ghost_dft_code_exe=path + "e_hydrogen_bond_exe.nw", ghost_xc=ghost_m06_2x,func_name='m06-2x')
                            except IndexError:
                                print "oh faied again! switching to 6-311g (grid: xcoarse)..."
                                f_m06_2x = open(path + "odft_m06-2x_6311g-xcoarse_template.nw",'r')
                                filedata=f_m06_2x.read()
                                f_m06_2x.close()
                                #> Now if atom 3 is N them we might get multiplicity error, so if-else loop required...
                                if case[0][0] != 'N':
                                    newdata = filedata.replace("  O##",("  "+ghost_final_dft_string)).replace("  #mult@@","  mult 2").replace("  #mon#","  mon second "+ghost_final_bsse_ghost)
                                elif case[4][0] == 'N' and (sum(ghost_electron_num[2:]) % 2) > 0:
                                    newdata = filedata.replace("  O##",("  "+ghost_final_dft_string)).replace("  #mult$$","  mult 2").replace("  #mon#","  mon second "+ghost_final_bsse_ghost)

                                dft_code_exe = open(path + "e_hydrogen_bond_exe.nw",'w')
                                dft_code_exe.write(newdata)
                                dft_code_exe.close()
                                #Now calling dft_analyzer function created earlier...
                                analyze_GHOST_dft_code(ghost_dft_code_exe=path + "e_hydrogen_bond_exe.nw", ghost_xc=ghost_m06_2x,func_name='m06-2x')



            #CHANGING THE STORED ENERGY VALUES INTO string...    
            ghost_m06_2x_string_dft_energy=[]   #converting the energy value in string and storing here for extending the HB-output_list later...
            for ghost_m06_2x_energy in ghost_m06_2x:
                ghost_m06_2x_string_value=str(ghost_m06_2x_energy)
                ghost_m06_2x_string_dft_energy.append(ghost_m06_2x_string_value)	    		


            subtractedEnergy = []  #subtracting cluster energy and Ghost-cluster energy to get the contribution of hydrogen bond forming atomic energies...
            for m, n in zip(m06_2x, ghost_m06_2x):
                diffEnergy = float(m) - float(n)  #Subtracting ...
                z = str(diffEnergy)
                subtractedEnergy.append(z)


            #--------------------------------------------------------------------------------------
            #Generating the .HBradar outfile with all info. and hydrogen bond_energy...
            #including xc= m06-2x, m05-2x, m06-HF
            #--------------------------------------------------------------------------------------
            for items, m06_2x_value, ghost_m06_2x_value, e, clusterCount in zip(HB_output_list, m06_2x_string_dft_energy,ghost_m06_2x_string_dft_energy,subtractedEnergy,atomCountInEachCluster):
                #print items, m06_2x_value,m06_2x_value, m06_HF_value
                new_line=items +"    "+ m06_2x_value + "    "+ ghost_m06_2x_value +"    "+ e +"    "+ clusterCount
                hbradar_out.write(new_line)
                hbradar_out.write('\n')
                #print new_line
            hbradar_out.close()
            #"""

            #print "COMPLETED SUCCESSFULLY. . ."
            #print (str(h_count) +" probable hydrogen bond calculated")
            ################################################################
            #stop = timeit.default_timer()
        #print ("Total time taken to analyze this file: %.3f sec" % (stop - start))
        #print "OUTPUT FILE GENERATED SUCCESSFULLY AS: " + hbradar

        print ">>> ALL DIRECTORIES SCANNED, TOTAL FILE PROCESSED: " + str(f_count)   
        print "              Thanks for using me !"
        print "               Have a nice day !"

    return #'coz func-loop has to run till the path has files to be processed...

#Opening the file...
import time
#person_name = raw_input("Hello sir !, Please enter your name to proceed: ")
print "*********************************************************************************************************"
print "                  Hello sir! You are executing HB_radar algorithm on "+ time.asctime()+" IST "
print "                          PROGRAM NAME: HB_radar_Cluster_Calculator_m062x.v.3.0.py  "
print "                                             designed by => Abhisek Mondal"
print "                                                 PI => Dr. Saumen Datta"
print "                                               CSIR-IICB, KOLKATA, INDIA"
print "**********************************************************************************************************"
#file_path = "/home/saumen/abhisek/project_HBf/DFT_section/bacteria/segment_1/"
file_path = raw_input("ENTER THE PATH FOR THE FILE (viz. /home/usr/): ")
show_me_hbs(path = file_path)