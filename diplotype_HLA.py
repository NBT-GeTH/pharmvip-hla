#!/usr/bin/env python

import sys, getopt
from typing import BinaryIO, Dict, Any
from pprint import pprint
import pandas as pd
import numpy as np
import io
import re
import itertools
from itertools import combinations

def create_diplotype_table(df,haplotype_mapping,sampleid):
#    print (df)
	########################df(hla_result.txt)############################
	#         hla_res_tool hla_res_original hla_res_2digit hla_res_4digit hla_res_6digit  hla_res_8digit hla_res_type
	#    0       ATHLATES    A*24:02:01:14           A*24        A*24:02     A*24:02:01   A*24:02:01:14            A
	#    1       ATHLATES       A*02:03:01           A*02        A*02:03     A*02:03:01            None            A
	#    2       ATHLATES   A*24:02:01:02L           A*24        A*24:02     A*24:02:01  A*24:02:01:02L            A
	#    3       ATHLATES    A*24:02:01:03           A*24        A*24:02     A*24:02:01   A*24:02:01:03            A
	#    4       ATHLATES    A*24:02:01:04           A*24        A*24:02     A*24:02:01   A*24:02:01:04            A
	
	#######################haplotype_mapping_HLA.txt#############################
#    gene    diplotype
#    HLA-A    HLA-A*31:01:02
#    HLA-B    HLA-B*15:02:01
#    HLA-B    HLA-B*57:01:01
#    HLA-B    HLA-B*58:01

	next(haplotype_mapping)
	data_haplotype ={}
	#### create dictionary of haplotype_mapping
	for line in haplotype_mapping:
		line = line.strip()
		genename = line.split("\t")[0]
		haplotype = line.split("\t")[1]

		if genename in data_haplotype :
			data_haplotype[genename] = data_haplotype[genename]+","+haplotype
		else:
			data_haplotype[genename] = haplotype

	#pprint(data_haplotype)
	listgenehla = data_haplotype.keys() #list gene that have guideline
	guidedip_tmp ={}
	printdip_tmp ={}
	diplotype_guidedip={}
	diplotype_printdip ={}
	tool_HLA_all = ["ATHLATES","HLA-HD","KOURAMI"]

	#create df of each tool and match haplotype_mapping
	for key, value in data_haplotype.items(): #HLA-A, HLA-B
	    listhap = value.split(",")
	    dfhaplotype = df[df['hla_res_type'] == key]
	    dfhaplotype_all = dfhaplotype.filter(['hla_res_original','hla_res_type','hla_res_tool'])
	#     print (dfhaplotype_all)
	    tool_HLA = dfhaplotype_all['hla_res_tool'].unique()
	#     print (tool_HLA)
	    guidedip_tmp[key]={}
	    printdip_tmp[key]={}
	    diplotype_guidedip[key]={}
	    diplotype_printdip[key]={}
	    
	#     print ("result0") 
	#     pprint (guidedip_tmp)
	#     pprint (printdip_tmp)
	    for t in tool_HLA: #ATHLATES/ HLA-HD/ KOURAMI
	#         print (t)
	        dftool = dfhaplotype_all.loc[dfhaplotype_all['hla_res_tool'] == t]
	        dftool_original = dftool.filter(['hla_res_original','hla_res_type'])
	#         print(dftool_original)

	        for index, row in dftool_original.iterrows(): #each row that HLA-B
	            for l in listhap:
	                haplotype = row['hla_res_original']
	                lreg = l.replace('*', '\*')
	                regex = "^"+lreg
	#                 print (regex,haplotype)
	                x = re.search(regex, haplotype)
	                if x: #match
	                    haplotype_result = l
	                else: #no match
	                    haplotype_result = "Other"

	                if t in guidedip_tmp[key] :
	                    guidedip_tmp[key][t] = guidedip_tmp[key][t]+","+haplotype_result
	                else:
	                    guidedip_tmp[key][t] = haplotype_result

	                if t in printdip_tmp[key] :
	                    printdip_tmp[key][t] = printdip_tmp[key][t]+","+haplotype
	                else:
	                    printdip_tmp[key][t] = haplotype
	         
	        guidedip_tmp_list = guidedip_tmp[key][t].split(",")
	        guidedip_tmp_list_unique = list(set(guidedip_tmp_list))
	        guidedip_tmp_list_unique.sort()
	        diplotype_guidedip[key][t] = guidedip_tmp_list_unique

	        printdip_tmp_list = printdip_tmp[key][t].split(",")
	        printdip_tmp_list_unique = list(set(printdip_tmp_list))
	        printdip_tmp_list_unique.sort()
	        diplotype_printdip[key][t] = printdip_tmp_list_unique
			
	f= open(f"{sampleid}_diplotype_hla.tsv","w+")
	#f.write(f"sampleid\tgene\ttool\tguide_diplotype\tprint_diplotype\n")
	f.write(f"sampleid\tgene\ttool_detail\tguide_diplotype_detail\tprint_diplotype_detail\ttool\tguide_diplotype\tprint_diplotype\n")
	guidelist = {} 
	printlist = {}
	for l in listgenehla: #HLA-A, HLA-B
	    guidelist[l] = {}
	    printlist[l] = {}
	    for t in tool_HLA_all: #ATHLATES HLA-HD KOURAMI
	        if t in diplotype_printdip[l]:
	            #print (diplotype_printdip[l][t])
	            tmp_printlist = diplotype_printdip[l][t]
	            tmp_printlist.sort()
	            if len(tmp_printlist) ==1:
	                alter_diplotype_printdip = f"{tmp_printlist[0]}/{tmp_printlist[0]}"
	                if t in printlist[l] :
	                    printlist[l][t] = printlist[l][t]+","+alter_diplotype_printdip
	                else:
	                    printlist[l][t] = alter_diplotype_printdip
	            if len(tmp_printlist) ==2:
	                alter_diplotype_printdip = f"{tmp_printlist[0]}/{tmp_printlist[1]}"
	                if t in printlist[l] :
	                    printlist[l][t] = printlist[l][t]+","+alter_diplotype_printdip
	                else:
	                    printlist[l][t] = alter_diplotype_printdip
	            if len(tmp_printlist) >2:
	                comb = combinations(tmp_printlist, 2)
	                for i in list(comb):
	                    alter_diplotype_printdip = f"{i[0]}/{i[1]}"
	                    if t in printlist[l] :
	                        printlist[l][t] = printlist[l][t]+","+alter_diplotype_printdip
	                    else:
	                        printlist[l][t] = alter_diplotype_printdip

	        else:
	            alter_diplotype_printdip = f"?/?"
	           #printlist_alltool.append(alter_diplotype_printdip)
	            if t in printlist[l] :
	                printlist[l][t] = printlist[l][t]+","+alter_diplotype_printdip
	            else:
	                printlist[l][t] = alter_diplotype_printdip


	        if t in diplotype_guidedip[l]:
	            #print (diplotype_guidedip[t][l])
	            tmp_guidelist = diplotype_guidedip[l][t]
	            tmp_guidelist.sort()
	            if len(tmp_guidelist) ==1:
	                alter_diplotype_guidedip = f"{tmp_guidelist[0]}/{tmp_guidelist[0]}"
	                if t in guidelist[l] :
	                    guidelist[l][t] = guidelist[l][t]+", "+alter_diplotype_guidedip
	                else:
	                    guidelist[l][t]= alter_diplotype_guidedip
	            if len(tmp_guidelist) ==2:
	                alter_diplotype_guidedip = f"{tmp_guidelist[0]}/{tmp_guidelist[1]}"
	                if t in guidelist[l] :
	                    guidelist[l][t] = guidelist[l][t]+", "+alter_diplotype_guidedip
	                else:
	                    guidelist[l][t]= alter_diplotype_guidedip
	            if len(tmp_guidelist) >2:
	                for i in range(len(tmp_guidelist)-1) :
	                    alter_diplotype_guidedip = f"{tmp_guidelist[i]}/{tmp_guidelist[-1]}"
	#                     print (alter_diplotype_guidedip)
	                    if t in guidelist[l] :
	                        guidelist[l][t] = guidelist[l][t]+", "+alter_diplotype_guidedip
	                    else:
	                        guidelist[l][t] = alter_diplotype_guidedip

	#         else: #no ?/? in guidelist
	#             alter_diplotype_guidedip = f"?/?"
	#             if t in guidelist[l] :
	#                 guidelist[l][t] = guidelist[l][t]+", "+alter_diplotype_guidedip
	#             else:
	#                 guidelist[l][t]= alter_diplotype_guidedip

	# pprint (guidelist)
	# pprint (printlist)
			
	######print to file#########
	for l in listgenehla:
		f.write(f"{sampleid}\t")
		f.write(f"{l}\t") #gene in haplotype_mapping
		listtool =[]
		listguide = []
		listprint = []
		listtoolsum =[]
		listguidesum = []
		listprintsum = []
		
		##guidelist
		for tool, guidedip in guidelist[l].items(): #each tool
			listtool.append(f"{tool}")
			listguide.append(f"{guidedip}")
			printdip = printlist[l][tool]
			listprint.append(f"{printdip}")
		guidesum = {}
		for i in range(len(listguide)):
			if listguide[i] in guidesum:
				guidesum[listguide[i]] = guidesum[listguide[i]]+","+listtool[i]
			else:
				guidesum[listguide[i]] = listtool[i]

		for g,t in guidesum.items():
				listtoolsum.append(f"{t}")
				listguidesum.append(f"{g}")
				listprintsum.append(f"{g}")
		#listToStr_tool = ','.join([str(elem) for elem in listtool])
		#listToStr_guide = ','.join([str(elem) for elem in listguide])
	   # listToStr_print = ','.join([str(elem) for elem in listprint])
		f.write(f"{listtool}\t")
		f.write(f"{listguide}\t")
		f.write(f"{listprint}\t")
		f.write(f"{listtoolsum}\t")
		f.write(f"{listguidesum}\t")
		f.write(f"{listprintsum}\n")
		
	f.close()

	return


#input = hla_result.txt
hladata = sys.argv[1]
mapping = sys.argv[2]
sampleid = sys.argv[3]
df = pd.read_csv(hladata,sep='\t',header=(0))
mappingfile = open(mapping, "r")

create_diplotype_table(df,mappingfile,sampleid)
