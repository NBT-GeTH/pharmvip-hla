#!/usr/bin/env python

import sys, getopt
from typing import BinaryIO, Dict, Any
from pprint import pprint
import pandas as pd
import numpy as np
import io
import re


def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print('HLA.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('HLA.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
     #   elif opt in ("-o", "--ofile"):
      #      outputfile = arg
    #print('Input file is "', inputfile)
    #print('Output file is "', outputfile)
    return inputfile, outputfile

def create_result_table(dict_HLA):
    pprint (dict_HLA)

    f = open('hla_results.txt', 'w')
    df = pd.DataFrame(columns=['hla_res_tool','hla_res_original','hla_res_2digit','hla_res_4digit','hla_res_6digit','hla_res_8digit','hla_res_type'])
    f.write("hla_res_tool\thla_res_original\thla_res_2digit\thla_res_4digit\thla_res_6digit\thla_res_8digit\thla_res_type\n")
    for types in dict_HLA:
        for tool in dict_HLA[types]:
            HLAtyping = dict_HLA[types][tool].split(',')
            for t in HLAtyping:
                r2 = t
                if t[-1] == 'G':
                    r2 = t[:-1] 

                digit2 = "None"
                digit4 = "None"
                digit6 = "None"
                digit8 = "None"

                if r2.count(':') == 0:
                    digit2 = r2

                if r2.count(':') == 1:
                    digit2 = r2.split(':')[0]
                    digit4 = ':'.join(r2.split(':')[0:2])
                if r2.count(':') == 2:
                    digit2 = r2.split(':')[0]
                    digit4 = ':'.join(r2.split(':')[0:2])
                    digit6 = ':'.join(r2.split(':')[0:3])
                if r2.count(':') == 3:
                    digit2 = r2.split(':')[0]
                    digit4 = ':'.join(r2.split(':')[0:2])
                    digit6 = ':'.join(r2.split(':')[0:3])
                    digit8 = ':'.join(r2.split(':')[0:4])
                    
                 #   print(tool, r, digit2, digit4, digit6, digit8, type, sep='\t')
                f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (tool, t, digit2, digit4, digit6, digit8, types))
                df = df.append({'hla_res_tool': tool,'hla_res_original': t,'hla_res_2digit': digit2,'hla_res_4digit': digit4,'hla_res_6digit': digit6,'hla_res_8digit': digit8,'hla_res_type': types}, ignore_index=True)
    f.close 

    return df


def create_summary_table(df):

    type_digit = {'hla_res_2digit', 'hla_res_4digit', 'hla_res_6digit', 'hla_res_8digit'}
    #print (df)
    #####################################################
    #         hla_res_tool hla_res_original hla_res_2digit hla_res_4digit hla_res_6digit  hla_res_8digit hla_res_type
    #    0       ATHLATES    A*24:02:01:14           A*24        A*24:02     A*24:02:01   A*24:02:01:14            A
    #    1       ATHLATES       A*02:03:01           A*02        A*02:03     A*02:03:01            None            A
    #    2       ATHLATES   A*24:02:01:02L           A*24        A*24:02     A*24:02:01  A*24:02:01:02L            A
    #    3       ATHLATES    A*24:02:01:03           A*24        A*24:02     A*24:02:01   A*24:02:01:03            A
    #    4       ATHLATES    A*24:02:01:04           A*24        A*24:02     A*24:02:01   A*24:02:01:04            A
    
    df = df.sort_values(by=['hla_res_type'])
    type_HLA = df.hla_res_type.unique() # all type HLA unique sort

    f = open('hla_summary.txt', 'w')
    f.write("hla_sum_type\thla_sum_calculate\thla_sum_digit\thla_sum_tool\thla_sum_count\n")

    
    for t in type_HLA:
        #print ("type\t"+t)
        for d in type_digit:
            #print (d)
            df_tmp = df.loc[(df['hla_res_type'] == t) & (df[d].str.contains('None') == False)]
            df_type_filter = df_tmp[['hla_res_type',d,'hla_res_tool']]
            #print (df_type_filter)
            HLA_calculate = df_type_filter[d].unique()
         #   print(HLA_calculate)
            for c in HLA_calculate:
                #print (t,d,c)
                df_tool = df_type_filter.loc[(df_type_filter[d] == c)]
                #print (df_tool)
                HLA_tool = df_tool['hla_res_tool'].unique()
                HLA_sum_count =  df_tool['hla_res_tool'].nunique()
                list_HLA_tool = ','.join(HLA_tool)
                m = re.findall(r'\d', d) 
                digit_data = ' '.join(map(str, m)) 
                f.write("%s\t%s\t%s\t%s\t%s\n" % (t, c, digit_data, list_HLA_tool, HLA_sum_count))

    f.close
    return



    # for type in results_type:
    #     typeHLA = type[0]
    #     for d in type_digit:
    #         sql_each_type = "SELECT DISTINCT %s FROM HLA_result WHERE type = '%s' and %s <> 'None' and project ='%s' " % (d, typeHLA, d, last_project)
    #         cursor.execute(sql_each_type)
    #         results_each_type = cursor.fetchall()
    #         for row in results_each_type:
    #             HLA_cal = row[0]
    #             results_tool = None
    #             sql_each_cal = "SELECT DISTINCT tool FROM HLA_result WHERE %s = '%s' and project ='%s'" % (d, HLA_cal, last_project)
    #             cursor.execute(sql_each_cal)
    #             results_each_cal = cursor.fetchall()
    #             for row2 in results_each_cal:
    #                 if results_tool is None:
    #                     results_tool = row2[0]
    #                 else:
    #                     results_tool = results_tool + "," + row2[0]

    #             count_tool = results_tool.count(',') + 1
    #            # print(typeHLA, HLA_cal, results_tool, count_tool)
    #             sql = "INSERT INTO HLA_summary(project,user,type,HLA_calculate,tool,count ) " \
    #                   "VALUES ( '%s','1','%s', '%s', '%s', '%d')" % (last_project, typeHLA, HLA_cal, results_tool, count_tool)
    #             try:
    #                 cursor.execute(sql)
    #                 db.commit()
    #             except:
    #                 db.rollback()


    # db.close()
    # return

if __name__ == "__main__":
    inputfile, outputfile = main(sys.argv[1:])
    # Open a file
    inputHLA = inputfile.split(',')
    print(inputHLA)
    dict_HLA = {}
    classI ={"A","B","C","E","F","G","H","J","K","L","N","P","S","T","U","V","W","Y"}
    classII = {"DRA","DRB","DQA1","DQA2","DQB1","DPA1","DPA2","DPB1","DPB2","DMA","DMB","DOA","DOB","DRB1","DRB2","DRB3","DRB4","DRB5","DRB6","DRB7","DRB8","DRB9"}
    classother = {"HFE","MICA","MICB","TAP1","TAP2"}
    
    for x in inputHLA:
       # print(x)
        if x[0:5] == 'HLAHD':
            print('get HLAHD')
            with open(x, "r") as f:
                for line in f:
                    line = line.strip()
                    if line != "Couldn't read result file.": #### found some error line from HLAHD process
                        HLAHD_data = line.split('\t')
                        HLAHDtmptype = HLAHD_data[0]
                        HLAHD_type = HLAHDtmptype
                        if HLAHDtmptype in classI or classII:
                            HLAHD_type = "HLA-"+HLAHDtmptype

                        if not HLAHD_type in dict_HLA:
                            dict_HLA[HLAHD_type] = {}

                        if '*' in HLAHD_data[1]:
                            HLAHD_A = HLAHD_data[1]
                            dict_HLA[HLAHD_type]['HLA-HD'] = HLAHD_A

                        if '*' in HLAHD_data[2]:
                            HLAHD_B = HLAHD_data[2]
                            if not HLAHD_B in dict_HLA[HLAHD_type]['HLA-HD']:
                                dict_HLA[HLAHD_type]['HLA-HD'] = dict_HLA[HLAHD_type]['HLA-HD'] + ',' + HLAHD_B
            f.close()
        if x[0:5] == 'ATHLA':
            print('get ATHLATES')
            dict_check = {}
            with open(x, "r") as f:
                #next(f)
                for line in f:
                  #  print (line)
                    line = line.strip()
                    ATHLA_data = line.split()
                    ATHLA_A = ATHLA_data[0]
                    ATHLA_B = ATHLA_data[1]
                    dict_check[ATHLA_A] = 1
                    dict_check[ATHLA_B] = 1
            #pprint(dict_check)
            for d in dict_check:
                ATHLAtmptype = d.split('*')[0]
                ATHLA_type = ATHLAtmptype
                if ATHLAtmptype in classI or classII:
                    ATHLA_type = "HLA-"+ATHLAtmptype
                    d = "HLA-"+d

                if not ATHLA_type in dict_HLA:
                    dict_HLA[ATHLA_type] = {}
                    dict_HLA[ATHLA_type]['ATHLATES'] = d
                if not d in dict_HLA[ATHLA_type]['ATHLATES']:
                    dict_HLA[ATHLA_type]['ATHLATES'] = dict_HLA[ATHLA_type]['ATHLATES'] + ',' + d


            #pprint(dict_HLA)
            f.close()

        if x[0:5] == 'KOURA':
            print('get KOURA')
            first = "false"
            with open(x, "r") as f:
                for line in f:
                    line = line.strip()
                    KOURA_data = line.split('\t')
                    KOURAtmptype = KOURA_data[0].split('*')[0]
                    KOURA_type = KOURAtmptype
                    stringadd =''
                    if KOURAtmptype in classI or classII:
                        KOURA_type = "HLA-"+KOURAtmptype
                        stringadd = 'HLA-'

                    if not KOURA_type in dict_HLA:
                        dict_HLA[KOURA_type] = {}
                    if first == "true":
                        first = "false"
                        tmpKOURA_B = KOURA_data[0]
                        if ";" in tmpKOURA_B:
                            tmpKOURA_B = tmpKOURA_B.replace(";", ",")
                            ttB = tmpKOURA_B.split(',')
                            listKOURA_B = [stringadd + x for x in ttB]
                            KOURA_B = ','.join([str(elem) for elem in listKOURA_B])
                        else:
                            KOURA_B = stringadd + tmpKOURA_B

                        if not KOURA_B in dict_HLA[KOURA_type]['KOURAMI']:
                            dict_HLA[KOURA_type]['KOURAMI'] = dict_HLA[KOURA_type]['KOURAMI'] + ',' + KOURA_B
                        

                    else:
                        first = "true"
                        tmpKOURA_A = KOURA_data[0]
                        if ";" in tmpKOURA_A:
                            ttKOURA_A = tmpKOURA_A.replace(";", ",")
                            ttA = ttKOURA_A.split(',')
                            listKOURA_A = [stringadd + x for x in ttA]
                            KOURA_A = ','.join([str(elem) for elem in listKOURA_A]) 
                        else:
                            KOURA_A = stringadd+tmpKOURA_A
                        print ("Add"+KOURA_A)    
                        dict_HLA[KOURA_type]['KOURAMI'] = KOURA_A

            f.close()

        # if x[0:5] == 'arcas':
        #     print('get arcas')
        #     with open(x, "r") as f:
        #         next(f)
        #         for line in f:
        #             line = line.strip()
        #             arcas_ori = line.split('\t')
        #             for arcas_data in arcas_ori:
        #                 if "*" in arcas_data :
        #                     arcas_type = arcas_data.split('*')[0]
        #                     #print (arcas_type)
        #                     if not arcas_type in dict_HLA:
        #                         dict_HLA[arcas_type] = {}
        #                     #print (arcas_A)
        #                     if not 'arcasHLA' in dict_HLA[arcas_type]:
        #                         dict_HLA[arcas_type]['arcasHLA'] = arcas_data
        #                     if not arcas_data in dict_HLA[arcas_type]['arcasHLA']:
        #                         dict_HLA[arcas_type]['arcasHLA'] = dict_HLA[arcas_type]['arcasHLA'] + ',' + arcas_data

        #     f.close()

    df = create_result_table(dict_HLA)
    create_summary_table(df)
