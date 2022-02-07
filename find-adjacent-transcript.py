#!/usr/bin/env python3
# coding:utf-8
import argparse
import os
import sys
parser = argparse.ArgumentParser(description="find adjacent fusion transcript")
parser.add_argument('-gtf', "--gtf", required=True,  help="gft file for adjacent fusion transcript")
parser.add_argument('-ref', "--reference", required=True,  help="reference gtf/gff file")
parser.add_argument('-o', "--out", required=True,  help="out dir,please note: outdir end with /")
args = parser.parse_args()
out_dir = os.path.dirname(args.out)
if out_dir and not os.path.exists(out_dir):
    os.makedirs(out_dir)
def FGversion(x,y,out_dir):#

    os.system("gffcompare -r %s -o %ssense %s"%(x,out_dir,y))
    os.system("gffcompare -r %s -o %sre %s"%(y,out_dir,x))
#FGversion(args.gtf,args.reference,args.out)
def splitGffTo(y,z,out_dir):#
    tranid=[]
    same=['=','c','e','j','k','o']
    same2=['=','c','e','j','k']
    for line in open('%sre.tracking'%out_dir):
        if line.split('\t')[3] in same:
            tranid.append( (line.split('\t')[4]).split('|')[1])
    tranid.reverse();OtoO=[]
    myfile1=open('%sAll.txt'%out_dir,'w')
    while len(tranid)>0:
        print(len(tranid))
        pop=tranid.pop()
        if 'PB.%s'%(pop.split('.')[1]) in OtoO:
            pass
        else:
            myfile=open('%s%s.gtf'%(out_dir,pop),'w')
            for line in open(z):
                if '%s"'%pop in line:
                    myfile.write(line)
            myfile.close()
            os.system("gffcompare -r %s%s.gtf -o %s%s %s"%(out_dir,pop,out_dir,pop,y))
            for line1 in open('%s%s.tracking'%(out_dir,pop)):#
                if line1.split('\t')[3] in same2:#
                    myfile1.write(line1)
            os.system("rm %sPB.*"%out_dir)
#splitGffTo('FGRAMPH1-re.tracking','gff-FGRAMPH1.gff3','FilterPT-Hy-nofusion.collapsed.gff')
#splitGffTo(args.reference,args.gtf,args.out)

def FindFgene(out_dir):#
    tran=[];gene=[]
    for line in open('%sAll.txt'%out_dir):
        tran.append((line.split('\t')[2]).split('|')[1])
        gene.append((line.split('\t')[2]).split('|')[0])
    tran2=list(set(tran))
    PT=[]
    gene=list(set(gene))
    print('All gene %s'%(len(gene)))
    for i in tran2:
        num=tran.count(i)
        if num>1:
            PT.append('PB.%s'%(i.split('.')[1]))
    PT=list(set(PT))
    print('PT gene %s'%(len(PT)))
    FGPT=[];FusionFG=[];remaind=[]
    while len(gene)>0:
        #print len(gene)
        pop=gene.pop();FG=[]
        if pop in PT:
            pass
        else:
            FG=[]
            for line in open('%sAll.txt'%out_dir):
                if '\t%s|'%pop in line:
                    FG.append(((line.split('\t')[4]).split('.')[0]).split(':')[1])
            FG=list(set(FG))
            if len(FG)==1:
                FGPT.append(pop)
            else:
                FusionFG.append(pop)
    print('1to1 %s'%(len(FGPT)))
    print('F gene %s'%(len(FusionFG)))
    print(FusionFG)
    myfile=open('%sFinal-adjacent-transcript-list.txt'%out_dir,'w')
    myfile.write('Gtf-gene\tGtf-transcript\tRef\tNumber of reference gene\n')
    while len(FusionFG)>0:
        #print len(FusionFG)
        n=1
        pop=FusionFG.pop();FGfusion=[]
        print(pop)
        for line in  open('%sAll.txt'%out_dir):
            if '\t%s|'%pop in line:
                FGfusion.append(((line.split('\t')[4]).split('.')[0]).split(':')[1])
        FGfusion=list(set(FGfusion))
        while len(FGfusion)>0:
            FGpop=FGfusion.pop()
            for line in open('%sAll.txt'%out_dir):
                if ':%s.path1'%FGpop in line:
                    #print (line.split('\t')[2]).split('|')[1],FGpop,n
                    myfile.write('%s\t%s\t%s\t%s\n'%(pop,(line.split('\t')[2]).split('|')[1],FGpop,n))
            n+=1 
    myfile.close()
#
FindFgene(args.out)