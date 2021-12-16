import pandas as pd
import numpy as np
import sys
import os
import inspect
import subprocess
import glob
formats_name = [
    "PED_MAP",
   	 "VCF",
	"RAW",
	"BED_BIM_FAM",
	"GEN_SAMPLE",
	"BGEN_SAMPLE",
	"_23andme",
	"AncestryDNA",
 "HAP_SAMPLE_LEGEND",
]
tools_name = [
    "Samtools",
    "Plink",
    "Gtool",
]
functionnames = []

def PED_MAPtoVCF(input_file):
    os.system("./plink --file "+input_file+" --recode vcf --out "+input_file.split(".")[0])    
 
    
    
def PED_MAPtoRAW(input_file):
    os.system("./plink --file "+input_file+" --recodeA --out "+input_file.split(".")[0])    
    print(inspect.stack()[1][3])
    
def PED_MAPtoBED_BIM_FAM(input_file):
    print(inspect.stack()[1][3])
    os.system("./plink --file "+input_file+" --make-bed --out "+input_file.split(".")[0])    
    
def PED_MAPtoGEN_SAMPLE(input_file):
    #os.system("./gtool  -P --ped "+input_file+".ped --map "+input_file+".map --og "+input_file+".gen --os "+input_file+".sample")
    os.system("./plink2 --ped  "+input_file+" dosage=DS --export oxford --out "+input_file)
    
def PED_MAPtoBGEN_SAMPLE(input_file):
    print("Not Implemented yet")
    
def PED_MAPto_23andme(input_file):
    PED_MAPtoVCF(input_file)
    VCFto_23andme(input_file+".vcf")
    
def PED_MAPtoAncestryDNA(input_file):
    PED_MAPtoVCF(input_file)
    VCFtoAncestryDNA(input_file+".vcf")
    
def PED_MAPtoHAP_SAMPLE_LEGEND(input_file):
    PED_MAPtoVCF(input_file)
    VCFtoHAP_SAMPLE_LEGEND(input_file+".vcf")
     
    

    
def VCFtoPED_MAP(input_file):
    print(inspect.stack()[1][3])
    # Include full file name. For example test.vcf
    os.system("./plink --vcf "+input_file+" --recode --out "+input_file.split(".")[0])    
  
def VCFtoRAW(input_file):
    # Include full file name. For example test.vcf

    os.system("./plink --vcf "+input_file+" --recodeA --out "+input_file.split(".")[0])    
    print(inspect.stack()[1][3])

def VCFtoBED_BIM_FAM(input_file):
    # Include full file name. For example test.vcf
    os.system("./plink --vcf "+input_file+" --make-bed --out "+input_file.split(".")[0])    

def VCFtoGEN_SAMPLE(input_file):
    #VCFtoPED_MAP(input_file)
    #os.system("./gtool  -P --ped "+input_file.split(".")[0]+".ped --map "+input_file.split(".")[0]+".map --og "+input_file.split(".")[0]+".gen --os "+input_file.split(".")[0]+".sample")
    # Include full file name. For example test.vcf
    os.system("./plink --vcf "+input_file+" --export oxford --out "+input_file.split(".")[0])
    print(inspect.stack()[1][3])

def VCFtoBGEN_SAMPLE(input_file):
    #os.system("./plink --vcf "+input_file+" dosage=DS --export oxford --out "+input_file.split(".")[0])
    print(inspect.stack()[1][3])

def VCFto_23andme(input_file):
    # Include full file name without extension. For example test

    if not os.path.isdir("23andme"):
        os.mkdir("23andme")
        
    VCFtoBED_BIM_FAM(input_file)
    os.system("bcftools query -l "+input_file+" > ./23andme/temp_samples.txt")
    f = open("./23andme/temp_samples.txt", "r")
    for x in f:
        temp = open("./23andme/temp.txt", "w")
        print("X")
        temp.write(x.strip('\n').split("_")[0] +"  "+x.strip('\n').split("_")[1])
        temp.close()
        os.system("./plink --bfile "+input_file.split(".")[0]+" --keep ./23andme/temp.txt --recode 23 --snps-only --out ./23andme/"+x.strip('\n')+".23andme")
        #print("./plink --bfile "+input_file.split(".")[0]+" --keep ./23andme/temp.txt --make-bed --out ./23andme/"+x.strip('\n'))  
    #os.system("./plink --vcf "+input_file+" --snps-only --recode 23 --out ./23andme/")


def VCFtoAncestryDNA(input_file):
  VCFto_23andme(input_file)
  if not os.path.isdir("AncestryDNA"):
    os.mkdir("AncestryDNA")
    
  _23andmefiles  = os.listdir('./23andme')
  for files in _23andmefiles:
    if ".23andme.txt" in files and "temp" not in files:
      if os.stat("./23andme"+os.sep+files).st_size == 0:
        continue
      else:
        data = pd.read_csv("./23andme"+os.sep+files,sep="\t",skiprows=8)
        new = pd.DataFrame()
        new['Rsid'] = data['# rsid'].values
        new['Chromosome'] = data['chromosome'].values
        new['position'] = data['position'].values
        new['allele1'] = data['genotype'].str[0]
        new['allele2'] =data['genotype'].str[1]
        new['Chromosome'] = new['Chromosome'].replace(23, 'X')
        new['Chromosome'] = new['Chromosome'].replace(24, 'Y')
        new['Chromosome'] = new['Chromosome'].replace(25, 'XY')
        new['Chromosome'] = new['Chromosome'].replace(26, 'MT')
        files = files.replace("23andme","ancestry")         
        new.to_csv("./AncestryDNA"+os.sep+files, sep="\t")
        

  
def VCFtoHAP_SAMPLE_LEGEND(input_file):
    os.system("bcftools convert "+input_file+"  -h  "+input_file.split(".")[0])



    
    
def RAWtoPED_MAP(input_file):
    
    def converter(x):
      x = x.fillna('0 0')
      x.astype(int, errors='ignore')
      ref = x.name[-1]
      if ref=="G":
        x = x.replace(0, "G G")
        x = x.replace(1, "G C")
        x = x.replace(2, "C C")
      
      if ref=="C":
        x = x.replace(0, "C C")
        x = x.replace(1, "C G")
        x = x.replace(2, "G G")
        
      if ref=="T":
        x = x.replace(0, "T T")
        x = x.replace(1, "T A")
        x = x.replace(2, "A A")
        
      if ref=="A":
        x = x.replace(0, "A A")
        x = x.replace(1, "A T")
        x = x.replace(2, "T T")
      
      
      return x
      
    # Extract SNPs names, which is in this format SNP_REFAllele
    #os.system("cat "+input_file+" | head -n 1  >> snps.txt")
    print("cat "+input_file+" | head -n 1  >> snps.txt")
    
    data = pd.read_csv("snps.txt",index_col=None,header=None,sep="\s+").loc[:, 6:].T
    data["ref"] = data[0].str[-1]    
    data["alt"] = data["ref"].replace({"C": "T", "T": "C","A": "G", "G": "A"})    
    if not os.path.isdir("Chunks"):
      os.mkdir("Chunks")   
    
    #Same file
    maps = pd.DataFrame()
    maps[0] = [0]*len(data)
    maps[1] = data[0].values
    maps[2] = [0]*len(data)
    maps[3] = [0]*len(data)
    maps.to_csv("final.map",sep="\t",header=False,index=False)
    _smallraw  = os.listdir('./Chunks')
    count=0
    _smallraw = sorted(_smallraw)
   
    os.system("split "+input_file+" ./Chunks/ ")
    
    for files in _smallraw:
      if ".txt" not in files:
        if count==0:
          count=1
          data2 = pd.read_csv("Chunks"+os.sep+files,sep="\s+")
          data2[list(data[0].values)] = data2[list(data[0].values)].apply(converter)
          data2.to_csv("Chunks"+os.sep+files+".txt",sep="\t",index=False,header=False)

        else:
          data2 = pd.read_csv("Chunks"+os.sep+files,sep="\s+",names=list(data2.columns.values))
          data2[list(data[0].values)] = data2[list(data[0].values)].apply(converter)
          data2.to_csv("Chunks"+os.sep+files+".txt",sep="\t",index=False,header=False)
    
    
    final = pd.DataFrame()
    #Merge the existing files
    for files in _smallraw:
      if ".txt" in files:
        if count==0:
          count=1
          final = pd.read_csv("Chunks"+os.sep+files,sep="\t",index_col=None,low_memory=False,header=None)
        else:
          data2 = pd.read_csv("Chunks"+os.sep+files,sep="\t",header=None,index_col=None,low_memory=False)
          final = final.append(data2, ignore_index=True)
          del data2
    final.to_csv("final.ped",sep="\t",index=False,header=None)    
 

def RAWtoVCF(input_file):
    #RAWtoPED_MAP(input_file)
    os.system("./plink --file final --recode vcf --out output_file")    
    print(inspect.stack()[1][3])
def RAWtoBED_BIM_FAM(input_file):
    #RAWtoPED_MAP(input_file)
    os.system("./plink --file final  --make-bed --out output_file")  
    print(inspect.stack()[1][3])
def RAWtoGEN_SAMPLE(input_file):
    #RAWtoPED_MAP(input_file)
    os.system("./plink --file final  --export oxford --out output_file")  
    print(inspect.stack()[1][3])
 
def RAWto_23andme(input_file):
    #RAWtoPED_MAP(input_file)
    os.system("./plink --file final  --make-bed --out output_file")  
    BED_BIM_FAMto_23andme("output_file")
 
    
def RAWtoAncestryDNA(input_file):
    #RAWtoPED_MAP(input_file)
    os.system("./plink --file final  --make-bed --out output_file")  
    BED_BIM_FAMto_23andme("output_file")
    _23andmeAncestryDNA("23andme")
 

def RAWtoHAP_SAMPLE_LEGEND(input_file):
    #RAWtoPED_MAP(input_file)
    os.system("./plink --file final --recode vcf --out output_file")
    os.system("bcftools convert output_file.vcf  -h  output_file2")
    print(inspect.stack()[1][3])



def BED_BIM_FAMtoPED_MAP(input_file):
    os.system("./plink --bfile "+input_file+" --recode --out "+input_file.split(".")[0])    
    print(inspect.stack()[1][3])
def BED_BIM_FAMtoVCF(input_file):
    os.system("./plink --bfile "+input_file+" --recode vcf --out "+input_file.split(".")[0])    
    print(inspect.stack()[1][3])
def BED_BIM_FAMtoRAW(input_file):
    os.system("./plink --bfile "+input_file+" --recodeA --out "+input_file.split(".")[0])    
    print(inspect.stack()[1][3])
    
def BED_BIM_FAMtoGEN_SAMPLE(input_file):
    os.system("./plink --bfile "+input_file+"  --export oxford --out "+input_file)
    print(inspect.stack()[1][3])
def BED_BIM_FAMtoBGEN_SAMPLE(input_file):
    print("Not implemented yet!")
    
    print(inspect.stack()[1][3])
def BED_BIM_FAMto_23andme(input_file):
    if not os.path.isdir("23andme"):
        os.mkdir("23andme")
    data = pd.read_csv(input_file+".fam",sep="\s+",header=None)
    print(data)
    data = data [[0,1]]
    data.to_csv("./23andme/temp_samples.txt",header=False,index=False,sep=" ")
    f = open("./23andme/temp_samples.txt", "r")
    for x in f:
        temp = open("./23andme/temp.txt", "w")
        temp.write(x)
        temp.close()
        os.system("./plink --bfile "+input_file.split(".")[0]+" --keep ./23andme/temp.txt --recode 23 --snps-only --out ./23andme/"+x.split(" ")[0]+"_"+x.split(" ")[0])
        
def BED_BIM_FAMtoAncestryDNA(input_file):
    print("Not implemented yet!")
    print(inspect.stack()[1][3])
    
def BED_BIM_FAMtoHAP_SAMPLE_LEGEND(input_file):
    #BED_BIM_FAMtoVCF(input_file)
    print("X")
    os.system("bcftools convert "+input_file+".vcf  -h  "+input_file)
    print(inspect.stack()[1][3])


def GEN_SAMPLEtoPED_MAP(input_file):
    os.system("./gtool  -G --g "+input_file+".gen --s "+input_file+".sample --ped "+input_file+".ped --map "+input_file+".map")
def GEN_SAMPLEtoVCF(input_file):
    GEN_SAMPLEtoPED_MAP(input_file)
    PED_MAPtoVCF(input_file)
    
def GEN_SAMPLEtoRAW(input_file):
    GEN_SAMPLEtoPED_MAP(input_file)
    PED_MAPtoRAW(input_file)
    
def GEN_SAMPLEtoBED_BIM_FAM(input_file):
    GEN_SAMPLEtoPED_MAP(input_file)
    PED_MAPtoBED_BIM_FAM(input_file)
    
 
    
def GEN_SAMPLEto_23andme(input_file):
    GEN_SAMPLEtoPED_MAP(input_file)
    PED_MAPto_23andme(input_file)
    
def GEN_SAMPLEtoAncestryDNA(input_file):
    GEN_SAMPLEtoPED_MAP(input_file)
    PED_MAPtoAncestryDNA(input_file)
    
def GEN_SAMPLEtoHAP_SAMPLE_LEGEND(input_file):
    GEN_SAMPLEtoPED_MAP(input_file)
    PED_MAPtoHAP_SAMPLE_LEGEND(input_file)    
    
 


def BGEN_SAMPLEtoPED_MAP(input_file):

    print(inspect.stack()[1][3])
def BGEN_SAMPLEtoVCF(input_file):
    print(inspect.stack()[1][3])
def BGEN_SAMPLEtoRAW(input_file):
    print(inspect.stack()[1][3])
def BGEN_SAMPLEtoBED_BIM_FAM(input_file):
    print(inspect.stack()[1][3])
def BGEN_SAMPLEtoGEN_SAMPLE(input_file):
    print(inspect.stack()[1][3])
def BGEN_SAMPLEto_23andme(input_file):
    print(inspect.stack()[1][3])
def BGEN_SAMPLEtoAncestryDNA(input_file):
    print(inspect.stack()[1][3])
def BGEN_SAMPLEtoHAP_SAMPLE_LEGEND(input_file):
    print(inspect.stack()[1][3])


def _23andmetoPED_MAP(input_file):
    #_23andmetoBED_BIM_FAM(input_file)
    BED_BIM_FAMtoPED_MAP("23andmetoBED")

    print(inspect.stack()[1][3])
def _23andmetoVCF(input_file):
    #_23andmetoBED_BIM_FAM(input_file)
    BED_BIM_FAMtoVCF("23andmetoBED")


    print(inspect.stack()[1][3])
def _23andmetoRAW(input_file):
    #_23andmetoBED_BIM_FAM(input_file)
    BED_BIM_FAMtoRAW("23andmetoBED")

def _23andmetoBED_BIM_FAM(input_file):
    # Place 23andme files in a new directory. 
    # input_file is actually input directory.
    
    if not os.path.isdir(input_file):
        print("Directory "+input_file+" does not exists...! Kindly place 23andme files in that directory.")
        exit(0)
    files = []
    allfiles = os.listdir("./"+input_file+"/")
    personname = [] 
    sexinfo = "0"
    for loop in allfiles:
      if ".23andme" in loop:
      
        data = loop.split(".")[0].split("_")
        personname.append(data[0])
        print(data)
        if data[5] =="XX":
          sexinfo = "2"
        elif data[5] =="XY":
          sexinfo = "1"
        else:
          sexinfo = "0" 
        
        os.system("./plink --23file ./"+input_file+os.sep+loop+" --snps-only --make-bed  --out ./"+input_file+os.sep+data[0])
        if os.path.exists("./"+input_file+os.sep+data[0]+".fam"):
          data2 = pd.read_csv("./"+input_file+os.sep+data[0]+".fam",header=None, sep="\s+")
          data2[0] = data[0]
          data2[1] = data[0]
          data2[4] = sexinfo    
          data2.to_csv("./"+input_file+os.sep+data[0]+".fam",sep="\t",header=False,index=False)
    allfiles =  os.listdir("./"+input_file+"/")
    count=0
    files=[]
    for loop in allfiles:
        if ".txt" in loop and ".bed" not in loop and ".fam" not in loop and ".bim" not in loop:
          print(loop)
          x = loop.split("_")[0]
          x = x + ".fam"
          me = os.path.exists("./"+input_file+os.sep+x)
          if me==True:
            x = x.split(".")[0]
            x = "./"+input_file+os.sep+x +".bed " + "./"+input_file+os.sep+x + ".bim " +  "./"+input_file+os.sep+x + ".fam"	
            files.append(x)
          else:
            count=count+1
          print(count," People removed due to missing fam file")
    with open("./"+input_file+os.sep+"All.txt", "w") as filehandle:
      for listitem in files:
        filehandle.write('%s\n' % listitem)
    
    os.system("./plink --merge-list ./"+input_file+os.sep+"/All.txt --make-bed --out 23andmetoBED")
    
    if os.path.exists("23andmetoBED.bed"):
      exit(0)
    else:
      allfiles =  os.listdir("./"+input_file+"/")
      count=0
      for loop in allfiles:
        if ".bed" in loop:
          x = loop
          x = x.split(".")[0]
    
          command = "./plink --bfile ./"+input_file+os.sep+x+" --exclude 23andmetoBED-merge.missnp --make-bed --out ./"+input_file+os.sep + x
          os.system(command)
      allfiles =  os.listdir("./"+input_file+"/")
      files=[]
      
      for loop in allfiles:
          if ".txt" in loop and ".bed" not in loop and ".fam" not in loop and ".bim" not in loop:
            print(loop)
            x = loop.split("_")[0]
            x = x + ".fam"
            me = os.path.exists("./"+input_file+os.sep+x)
            if me==True:
              x = x.split(".")[0]
              x = "./"+input_file+os.sep+x +".bed " + "./"+input_file+os.sep+x + ".bim " +  "./"+input_file+os.sep+x + ".fam"	
              files.append(x)
            else:
              count=count+1
            print(count," People removed due to missing fam file")
      with open("./"+input_file+os.sep+"All.txt", "w") as filehandle:
        for listitem in files:
          filehandle.write('%s\n' % listitem)      
      os.system("./plink --merge-list ./"+input_file+os.sep+"/All.txt --make-bed --out 23andmetoBED")    
      
 
def _23andmetoGEN_SAMPLE(input_file):
    #_23andmetoBED_BIM_FAM(input_file)
    BED_BIM_FAMtoGEN_SAMPLE("23andmetoBED")
    print(inspect.stack()[1][3])
def _23andmetoBGEN_SAMPLE(input_file):
    print(inspect.stack()[1][3])

def _23andmetoAncestryDNA(input_file):  
  if not os.path.isdir("AncestryDNA"):
    os.mkdir("AncestryDNA")
  _23andmefiles  = os.listdir("./"+input_file)
  for files in _23andmefiles:
    if "23andme.txt" in files and "temp" not in files:
      if os.stat("./"+input_file+os.sep+files).st_size == 0:
        continue
      else:
        print(files)
        data = pd.read_csv("./"+input_file+os.sep+files,sep="\t", comment='#',header=None,low_memory=False)
        new = pd.DataFrame()
        new['Rsid'] = data[0].values
        new['Chromosome'] = data[1].values
        new['position'] = data[2].values
        new['allele1'] = data[3].str[0]
        new['allele2'] =data[3].str[1]
        new['Chromosome'] = new['Chromosome'].replace(23, 'X')
        new['Chromosome'] = new['Chromosome'].replace(24, 'Y')
        new['Chromosome'] = new['Chromosome'].replace(25, 'XY')
        new['Chromosome'] = new['Chromosome'].replace(26, 'MT')
        files = files.replace("23andme","ancestry")
        new.to_csv("./AncestryDNA"+os.sep+files, sep="\t",index=False)
    
def _23andmetoHAP_SAMPLE_LEGEND(input_file):
    #_23andmetoBED_BIM_FAM(input_file)
    BED_BIM_FAMtoVCF("23andmetoBED")
    VCFtoHAP_SAMPLE_LEGEND("23andmetoBED.vcf")


def AncestryDNAtoPED_MAP(input_file):
    #AncestryDNAto_23andme(input_file)
    _23andmetoPED_MAP("23andme")
    
 
def AncestryDNAtoVCF(input_file):
    #AncestryDNAto_23andme(input_file)
    _23andmetoVCF("23andme")
 
def AncestryDNAtoRAW(input_file):
    #AncestryDNAto_23andme(input_file)
    _23andmetoRAW("23andme")
 
def AncestryDNAtoBED_BIM_FAM(input_file):
    #AncestryDNAto_23andme(input_file)
    _23andmetoBED_BIM_FAM("23andme")
 
def AncestryDNAtoGEN_SAMPLE(input_file):
    #AncestryDNAto_23andme(input_file)
    _23andmetoGEN_SAMPLE("23andme")
 
def AncestryDNAtoBGEN_SAMPLE(input_file):
    #AncestryDNAto_23andme(input_file)
    
    print(inspect.stack()[1][3])
def AncestryDNAto_23andme(input_file):
  if not os.path.isdir("23andme"):
    os.mkdir("23andme")
    
  _ancestryfiles  = os.listdir('./'+input_file)
  for files in _ancestryfiles:
    if ".txt" in files and "temp" not in files:
      if os.stat("./"+input_file+os.sep+files).st_size == 0:
        continue
      else:
        data = pd.read_csv("./"+input_file+os.sep+files,sep="\t",comment='#',low_memory=False)
        new = pd.DataFrame()
        new['Rsid'] = data['rsid'].values
        new['Chromosome'] = data['chromosome'].values
        new['position'] = data['position'].values
        new['genotype'] = data['allele1']+data['allele1']
        new['Chromosome'] = new['Chromosome'].replace(23, 'X')
        new['Chromosome'] = new['Chromosome'].replace(24, 'Y')
        new['Chromosome'] = new['Chromosome'].replace(25, 'XY')
        new['Chromosome'] = new['Chromosome'].replace(26, 'MT')
        files = files.replace("ancestry","23andme")        
        new.to_csv("./23andme"+os.sep+files, sep="\t",index=False,header=False)    

def AncestryDNAtoHAP_SAMPLE_LEGEND(input_file):
    #AncestryDNAto_23andme(input_file)
    _23andmetoHAP_SAMPLE_LEGEND("23andme")
 

def HAP_SAMPLE_LEGENDtoPED_MAP(input_file):
    #HAP_SAMPLE_LEGENDtoVCF(input_file)
    VCFtoPED_MAP("output_file.vcf")    
    
def HAP_SAMPLE_LEGENDtoVCF(input_file):
    
    data = pd.read_csv(input_file+".legend",index_col=False, sep="\s+")
    #Create a new ID column for legend file
    data["ID"] = "1:"+data['position'].astype(str)+"_"+data["a0"]+"_"+data["a1"]
    data =data[['ID','position','a0','a1']]
    data.to_csv(input_file+".legend",index=False, sep="\t")
    
    #Zip the legend file
    os.system("gzip "+input_file+".legend")
    #Rename input_file.haps.gz to input_file+.hap.gz
    os.rename(input_file+".haps.gz",input_file+".hap.gz")  
    
    #Rename input_file.sample to input_file.samples
    os.rename(input_file+".sample",input_file+".samples")
    
    os.system("bcftools convert  --haplegendsample2vcf  "+input_file+"  -o  output_file.vcf")
   
def HAP_SAMPLE_LEGENDtoRAW(input_file):
    #HAP_SAMPLE_LEGENDtoVCF(input_file)    
    VCFtoRAW("output_file.vcf")       

def HAP_SAMPLE_LEGENDtoBED_BIM_FAM(input_file):
    #HAP_SAMPLE_LEGENDtoVCF(input_file)    
    VCFtoBED_BIM_FAM("output_file.vcf")   
 

def HAP_SAMPLE_LEGENDtoGEN_SAMPLE(input_file):
    #HAP_SAMPLE_LEGENDtoVCF(input_file)    
    VCFtoGEN_SAMPLE("output_file.vcf")   
    
   
def HAP_SAMPLE_LEGENDto_23andme(input_file):
    #HAP_SAMPLE_LEGENDtoVCF(input_file)    
    if not os.path.isdir("23andme"):
        os.mkdir("23andme")
        
    VCFtoBED_BIM_FAM("output_file.vcf")
    os.system("bcftools query -l output_file.vcf > ./23andme/temp_samples.txt")
    f = open("./23andme/temp_samples.txt", "r")
    for x in f:
        temp = open("./23andme/temp.txt", "w")
        print("X")
        temp.write(x.strip('\n')+" "+x.strip('\n'))
        temp.write('\n')
        
        temp.close()
        os.system("./plink --bfile output_file --keep ./23andme/temp.txt --recode 23 --snps-only --out ./23andme/"+x.strip('\n')+".23andme")

    
    
    
def HAP_SAMPLE_LEGENDtoAncestryDNA(input_file):
    #HAP_SAMPLE_LEGENDto_23andme(input_file)        
    if not os.path.isdir("AncestryDNA"):
      os.mkdir("AncestryDNA")
    
    _23andmefiles  = os.listdir('./23andme')
    for files in _23andmefiles:
      if ".23andme.txt" in files and "temp" not in files:
        if os.stat("./23andme"+os.sep+files).st_size == 0:
          continue
        else:
          data = pd.read_csv("./23andme"+os.sep+files,sep="\t",skiprows=8)
          new = pd.DataFrame()
          new['Rsid'] = data['# rsid'].values
          new['Chromosome'] = data['chromosome'].values
          new['position'] = data['position'].values
          new['allele1'] = data['genotype'].str[0]
          new['allele2'] =data['genotype'].str[1]
          new['Chromosome'] = new['Chromosome'].replace(23, 'X')
          new['Chromosome'] = new['Chromosome'].replace(24, 'Y')
          new['Chromosome'] = new['Chromosome'].replace(25, 'XY')
          new['Chromosome'] = new['Chromosome'].replace(26, 'MT')
          files = files.replace("23andme","ancestry")          
          new.to_csv("./AncestryDNA"+os.sep+files, sep="\t",index=False)

if sys.argv[1]=="PED_MAPtoVCF":
    PED_MAPtoVCF(sys.argv[2])
if sys.argv[1]=="PED_MAPtoRAW":
    PED_MAPtoRAW(sys.argv[2])
if sys.argv[1]=="PED_MAPtoBED_BIM_FAM":
    PED_MAPtoBED_BIM_FAM(sys.argv[2])
if sys.argv[1]=="PED_MAPtoGEN_SAMPLE":
    PED_MAPtoGEN_SAMPLE(sys.argv[2])
if sys.argv[1]=="PED_MAPtoBGEN_SAMPLE":
    PED_MAPtoBGEN_SAMPLE(sys.argv[2])
if sys.argv[1]=="PED_MAPto_23andme":
    PED_MAPto_23andme(sys.argv[2])
if sys.argv[1]=="PED_MAPtoAncestryDNA":
    PED_MAPtoAncestryDNA(sys.argv[2])
if sys.argv[1]=="PED_MAPtoHAP_SAMPLE_LEGEND":
    PED_MAPtoHAP_SAMPLE_LEGEND(sys.argv[2])    
    
if sys.argv[1]=="VCFtoPED_MAP":
    VCFtoPED_MAP(sys.argv[2])
if sys.argv[1]=="VCFtoRAW":
    VCFtoRAW(sys.argv[2])
if sys.argv[1]=="VCFtoBED_BIM_FAM":
    VCFtoBED_BIM_FAM(sys.argv[2])
if sys.argv[1]=="VCFtoGEN_SAMPLE":
    VCFtoGEN_SAMPLE(sys.argv[2])
if sys.argv[1]=="VCFtoBGEN_SAMPLE":
    VCFtoBGEN_SAMPLE(sys.argv[2])
if sys.argv[1]=="VCFto_23andme":
    VCFto_23andme(sys.argv[2])
if sys.argv[1]=="VCFtoAncestryDNA":
    VCFtoAncestryDNA(sys.argv[2])
if sys.argv[1]=="VCFtoHAP_SAMPLE_LEGEND":
    VCFtoHAP_SAMPLE_LEGEND(sys.argv[2])    
    
    
    
if sys.argv[1]=="RAWtoPED_MAP":
    RAWtoPED_MAP(sys.argv[2])
if sys.argv[1]=="RAWtoVCF":
    RAWtoVCF(sys.argv[2])
if sys.argv[1]=="RAWtoBED_BIM_FAM":
    RAWtoBED_BIM_FAM(sys.argv[2])
if sys.argv[1]=="RAWtoGEN_SAMPLE":
    RAWtoGEN_SAMPLE(sys.argv[2])
if sys.argv[1]=="RAWtoBGEN_SAMPLE":
    RAWtoBGEN_SAMPLE(sys.argv[2])
if sys.argv[1]=="RAWto_23andme":
    RAWto_23andme(sys.argv[2])
if sys.argv[1]=="RAWtoAncestryDNA":
    RAWtoAncestryDNA(sys.argv[2])
if sys.argv[1]=="RAWtoHAP_SAMPLE_LEGEND":
    RAWtoHAP_SAMPLE_LEGEND(sys.argv[2])    
    
    
if sys.argv[1]=="BED_BIM_FAMtoPED_MAP":
    BED_BIM_FAMtoPED_MAP(sys.argv[2])
if sys.argv[1]=="BED_BIM_FAMtoVCF":
    BED_BIM_FAMtoVCF(sys.argv[2])
if sys.argv[1]=="BED_BIM_FAMtoRAW":
    BED_BIM_FAMtoRAW(sys.argv[2])
if sys.argv[1]=="BED_BIM_FAMtoGEN_SAMPLE":
    BED_BIM_FAMtoGEN_SAMPLE(sys.argv[2])
if sys.argv[1]=="BED_BIM_FAMtoBGEN_SAMPLE":
    BED_BIM_FAMtoBGEN_SAMPLE(sys.argv[2])
if sys.argv[1]=="BED_BIM_FAMto_23andme":
    BED_BIM_FAMto_23andme(sys.argv[2])
if sys.argv[1]=="BED_BIM_FAMtoAncestryDNA":
    BED_BIM_FAMtoAncestryDNA(sys.argv[2])
if sys.argv[1]=="BED_BIM_FAMtoHAP_SAMPLE_LEGEND":
    BED_BIM_FAMtoHAP_SAMPLE_LEGEND(sys.argv[2])

if sys.argv[1]=="GEN_SAMPLEtoPED_MAP":
    GEN_SAMPLEtoPED_MAP(sys.argv[2])
if sys.argv[1]=="GEN_SAMPLEtoVCF":
    GEN_SAMPLEtoVCF(sys.argv[2])
if sys.argv[1]=="GEN_SAMPLEtoRAW":
    GEN_SAMPLEtoRAW(sys.argv[2])
if sys.argv[1]=="GEN_SAMPLEtoBED_BIM_FAM":
    GEN_SAMPLEtoBED_BIM_FAM(sys.argv[2])
if sys.argv[1]=="GEN_SAMPLEtoBGEN_SAMPLE":
    GEN_SAMPLEtoBGEN_SAMPLE(sys.argv[2])
if sys.argv[1]=="GEN_SAMPLEto_23andme":
    GEN_SAMPLEto_23andme(sys.argv[2])
if sys.argv[1]=="GEN_SAMPLEtoAncestryDNA":
    GEN_SAMPLEtoAncestryDNA(sys.argv[2])
    
if sys.argv[1]=="GEN_SAMPLEtoHAP_SAMPLE_LEGEND":
    GEN_SAMPLEtoHAP_SAMPLE_LEGEND(sys.argv[2])
    
    
if sys.argv[1]=="BGEN_SAMPLEtoPED_MAP":
    BGEN_SAMPLEtoPED_MAP(sys.argv[2])
if sys.argv[1]=="BGEN_SAMPLEtoVCF":
    BGEN_SAMPLEtoVCF(sys.argv[2])
if sys.argv[1]=="BGEN_SAMPLEtoRAW":
    BGEN_SAMPLEtoRAW(sys.argv[2])
if sys.argv[1]=="BGEN_SAMPLEtoBED_BIM_FAM":
    BGEN_SAMPLEtoBED_BIM_FAM(sys.argv[2])
if sys.argv[1]=="BGEN_SAMPLEtoGEN_SAMPLE":
    BGEN_SAMPLEtoGEN_SAMPLE(sys.argv[2])
if sys.argv[1]=="BGEN_SAMPLEto_23andme":
    BGEN_SAMPLEto_23andme(sys.argv[2])
if sys.argv[1]=="BGEN_SAMPLEtoAncestryDNA":
    BGEN_SAMPLEtoAncestryDNA(sys.argv[2])
if sys.argv[1]=="_23andmetoPED_MAP":
    _23andmetoPED_MAP(sys.argv[2])
if sys.argv[1]=="_23andmetoVCF":
    _23andmetoVCF(sys.argv[2])
if sys.argv[1]=="_23andmetoRAW":
    _23andmetoRAW(sys.argv[2])
if sys.argv[1]=="_23andmetoBED_BIM_FAM":
    _23andmetoBED_BIM_FAM(sys.argv[2])
if sys.argv[1]=="_23andmetoGEN_SAMPLE":
    _23andmetoGEN_SAMPLE(sys.argv[2])
if sys.argv[1]=="_23andmetoBGEN_SAMPLE":
    _23andmetoBGEN_SAMPLE(sys.argv[2])
if sys.argv[1]=="_23andmetoAncestryDNA":
    _23andmetoAncestryDNA(sys.argv[2])
    
if sys.argv[1]=="_23andmetoHAP_SAMPLE_LEGEND":
    _23andmetoHAP_SAMPLE_LEGEND(sys.argv[2])    
    
if sys.argv[1]=="AncestryDNAtoPED_MAP":
    AncestryDNAtoPED_MAP(sys.argv[2])
if sys.argv[1]=="AncestryDNAtoVCF":
    AncestryDNAtoVCF(sys.argv[2])
if sys.argv[1]=="AncestryDNAtoRAW":
    AncestryDNAtoRAW(sys.argv[2])
if sys.argv[1]=="AncestryDNAtoBED_BIM_FAM":
    AncestryDNAtoBED_BIM_FAM(sys.argv[2])
if sys.argv[1]=="AncestryDNAtoGEN_SAMPLE":
    AncestryDNAtoGEN_SAMPLE(sys.argv[2])
if sys.argv[1]=="AncestryDNAtoBGEN_SAMPLE":
    AncestryDNAtoBGEN_SAMPLE(sys.argv[2])
if sys.argv[1]=="AncestryDNAto_23andme":
    AncestryDNAto_23andme(sys.argv[2])

if sys.argv[1]=="AncestryDNAtoHAP_SAMPLE_LEGEND":
    AncestryDNAtoHAP_SAMPLE_LEGEND(sys.argv[2])

if sys.argv[1]== "HAP_SAMPLE_LEGENDtoPED_MAP":
    HAP_SAMPLE_LEGENDtoPED_MAP(sys.argv[2])
if sys.argv[1]== "HAP_SAMPLE_LEGENDtoVCF":
    HAP_SAMPLE_LEGENDtoVCF(sys.argv[2])
if sys.argv[1]== "HAP_SAMPLE_LEGENDtoRAW":
    HAP_SAMPLE_LEGENDtoRAW(sys.argv[2])
if sys.argv[1]== "HAP_SAMPLE_LEGENDtoBED_BIM_FAM":
    HAP_SAMPLE_LEGENDtoBED_BIM_FAM(sys.argv[2])
if sys.argv[1]== "HAP_SAMPLE_LEGENDtoGEN_SAMPLE":
    HAP_SAMPLE_LEGENDtoGEN_SAMPLE(sys.argv[2])
if sys.argv[1]== "HAP_SAMPLE_LEGENDtoBGEN_SAMPLE":
    HAP_SAMPLE_LEGENDtoBGEN_SAMPLE(sys.argv[2])
if sys.argv[1]== "HAP_SAMPLE_LEGENDto_23andme":
    HAP_SAMPLE_LEGENDto_23andme(sys.argv[2])
if sys.argv[1]== "HAP_SAMPLE_LEGENDtoAncestryDNA":
    HAP_SAMPLE_LEGENDtoAncestryDNA(sys.argv[2])
    
###############################################################################################################Parser Functions
### Make functions
for loop in formats_name:
    for loop2 in formats_name:
        if loop2!=loop:
            #print("def "+loop+"to"+loop2+"(input_file):")
            functionnames.append(loop+"to"+loop2)
            #print("    print(inspect.stack()[1][3])")

 
################################################################################################################Execution
if sys.argv[1] not in functionnames:
    print("The function you called does not exists")


