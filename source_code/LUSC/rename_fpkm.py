#!/usr/bin/env python
#rename TCGA file param1 is sample_sheet_file
import sys
import re
import os
wd = os.getcwd()
param1=sys.argv[1]
fi = open(param1,"r")
os.mkdir('Expression_Count')
next(fi)
for fl in fi:
    fine = fl.strip().split('\t')
    os.system("cp {wd}/{fileid}/{filename} {wd}/{fileid}/{sampleid}.FPKM.{type}.txt.gz;cp {wd}/{fileid}/{sampleid}.FPKM.{type}.txt.gz ./Expression_Count".format(wd=wd,fileid=fine[0],filename=fine[1],sampleid=fine[5],type=fine[7].replace(' ','.')))
fi.close()
    
    
