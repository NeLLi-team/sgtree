import glob
import os
import datetime
import subprocess
import pandas as pd
import random
import copy 
from Bio import SeqIO
import itertools
import argparse
import math
import sys
import shutil
import numpy as np


# python3 sampler.py Class Bacteria 100 /global/projectb/scratch/jwhittak/jwhittak --exclude Acidobacteria

parser = argparse.ArgumentParser()

parser.add_argument("lvl", help="",
                    type=str)
parser.add_argument("name", help="",
                    type=str)
parser.add_argument("num", help="",
                    type=int)
parser.add_argument("dir", help="directory to save your sampled genomes to", type=str)
parser.add_argument("--table", help="table to sample from", default="/global/projectb/scratch/jwhittak/jwhittak/taxontable107697_22-feb-2019.xls", type=str)
parser.add_argument('--exclude', nargs='+', help='<Required> Set flag to exclude certain phyla/clades')
args = parser.parse_args()


tx_lvl = args.lvl
name = args.name
lst_ecl = args.exclude

if tx_lvl == "Class": 
    tx_num = 2

if tx_lvl == "Phylum": 
    tx_num = 1
    
if tx_lvl == "Order": 
    tx_num = 3
    
if tx_lvl == "Family": 
    tx_num = 4

if tx_lvl == "Domain": 
    tx_num = 0

if tx_lvl == "Genus": 
    tx_num = 5
    
currentdir = args.dir
print(name, currentdir, lst_ecl)
    
df = pd.read_csv(args.table, sep="\t", 
                 usecols=[0,1,4,6,7,8,9,10,11,12,13,14], low_memory=False)



df = df.loc[df['High Quality'] == "Yes"]
df = df.loc[df['Is Public'] == "Yes"]
# df = df.applymap(lambda stri: stri.replace(" ", "_") if type(stri) == "str" else stri)
df_bigcopy = copy.deepcopy(df)

df["fullname"] = df["Domain"]+"|"+ df["Phylum"] +"|"+ df["Class"] +"|"+ df["Order"] +"|"+ df["Family"] +"|"+ df["Genus"] +"|"+ df["Species"]

df_bigcopy = copy.deepcopy(df)

print(df["Phylum"].drop_duplicates())

# indicies=[]
# for row in df_bigcopy.iterrows(): 
#     if "candidatus" in row[1]["fullname"]: 
#         #print("ELIMINATING", row[1]["fullname"], row[1]["taxon_oid"], row[0])
#         indicies.append(row[0])
#     elif "Candidatus" in row[1]["fullname"]: 
#         #print("ELIMINATING", row[1]["fullname"], row[1]["taxon_oid"], row[0])
#         indicies.append(row[0])
#     elif "unknown" in row[1]["fullname"]: 
#         #print("ELIMINATING", row[1]["fullname"], row[1]["taxon_oid"], row[0])
#         indicies.append(row[0])
#     elif "unclassified|unclassified" in row[1]["fullname"]: 
#         print("ELIMINATING", row[1]["fullname"], row[1]["taxon_oid"], row[0])
#         indicies.append(row[0])
#     elif "Actinobacteria|Actinobacteria" in row[1]["fullname"]:
#         df[row[0]["fullname"]] = row[1]["fullname"].replace("Actinobacteria|Actinobacteria", "Actinobacteria|Actinobacteriia")
#         print("REPLACED", df.loc[row[0]])

# orig_stdout = sys.stdout
# f = open(currentdir + "/final_scripts/02data/" + name + '_out.txt', 'w')
# sys.stdout = f

df_bigcopy = df_bigcopy[df_bigcopy["fullname"].str.contains(name, case=True)]
if lst_ecl is not None: 
    df_bigcopy = df_bigcopy[~df_bigcopy["fullname"].str.contains('|'.join(lst_ecl), case=True)]
#print(df_bigcopy, tx_lvl)

lst_tx1 = {}
for i, each in df_bigcopy[tx_lvl].iteritems(): 
    if each not in lst_tx1.keys(): 
        lst_tx1[each] = [df_bigcopy["taxon_oid"][i]]
    else: 
        lst_tx1[each].append(df_bigcopy["taxon_oid"][i])
        
            
# sys.stdout = orig_stdout
# f.close()
            
# print(df_bigcopy[df_bigcopy.Domain.isin(["archaea"])])



directory1 = currentdir + "/data/" + name + "_by_" + tx_lvl
dest = directory1
            
if not os.path.exists(dest):
    os.makedirs(dest)
else: 
    shutil.rmtree(dest)
    os.makedirs(dest)

    
print(lst_tx1, "LEN", len(lst_tx1), lst_tx1.keys())


if len(lst_tx1.keys()) > 1: 
    orig_stdout = sys.stdout
    f = open(currentdir + "/data/" + name + "_by_" + tx_lvl + '_out.txt', 'w')
    sys.stdout = f

    seen_ids = []

    def one_round_sample(lst_tx1):
        """samples one genome from all of the the DomPhyClassOFG within the group to be sampled (Bacteria for example)"""
        for key, value_list in lst_tx1.items(): 

            if len(os.listdir(dest)) >= args.num: 
                print("we hit the number of genomes we desire")
                break

            elif all([x in seen_ids for x in value_list]): ## check whether every value from this list has been taken. 
                print("\n","###","skipped " + key + " we already sampled everything from it")
                print("SAMPLED IDS", seen_ids, "AND VALUE LIST", value_list,"\n","###",)
                continue

            value = random.sample(value_list, 1)[0]
            while value in seen_ids: 
                print(value, " in seen")
                value = random.sample(value_list, 1)[0]
            seen_ids.append(value)

            directroypath = "/global/projectb/sandbox/IMG/isolateRepo/" + str(value)
            # print(key, value, value_list, directroypath)
            for genomefaa in glob.glob(directroypath + "/*.faa"):
                if("IMG" + str(value) + ".faa" not in os.listdir(dest)):
                    for record in SeqIO.parse(genomefaa, "fasta"):
                        genomename = "IMG" + genomefaa.split(".")[0].split("/")[-1]
                        seqheader = genomename + "|" + record.id.split()[0]
                        record.id = seqheader
                        record.description = ""
                        with open(dest + "/IMG" + str(value) + ".faa", "a") as output_handle:
                            SeqIO.write(record, output_handle, "fasta")
                print(str(value), "\t", key, "\n")

    index = 0
    while len(os.listdir(dest)) < args.num: 
        """while we have less genomes in our destination directory than desired"""
        index += 1
        print("SEEN", seen_ids, "\n")
        print(len(os.listdir(dest)), args.num, "\n")
        # print(lst_tx1, "\n")
        one_round_sample(lst_tx1)
        ## if weve done more than 50 rounds theres something wrong
        if index >= 50: 
            print(name, " doesnt exist or there is another error/bug")
            break
    sys.stdout = orig_stdout
    f.close()
else: 
    print("something is wrong, couldnt find any valid ", tx_lvl, " for ", name, " --check your command line arguments")
    


