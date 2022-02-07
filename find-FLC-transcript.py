#!/usr/bin/env python3
# coding:utf-8
import argparse
import os
import sys
parser = argparse.ArgumentParser(description="predict full length chimeric sequences")
parser.add_argument('-fasta', "--fasta", required=True,  help="fasta file for predict full length chimeric sequences")
parser.add_argument('-genome', "--genome", required=True,  help="genome fasta file")
parser.add_argument('-percent', "--percent", required=True,  type=float,default=50.0,help="Pecentage of full length chimeric sequences coverage itself, def=50.0")
parser.add_argument('-o', "--out", required=True,  help="out dir,please note: outdir end with /")
args = parser.parse_args()
out_dir = os.path.dirname(args.out)
if out_dir and not os.path.exists(out_dir):
    os.makedirs(out_dir)
def findflc(out_dir,Percent,fastaFile):
    myfile=open('%sflc1-fT.pslx'%out_dir,'w')
    for line in open('%sflc.pslx'%out_dir):
        T=line.split('\t')[21].replace(',','')
        Tlist=T.upper()
        if (float(Tlist.count('T'))/float(len(Tlist)))*100>70:
            pass
        elif (float(Tlist.count('A'))/float(len(Tlist)))*100>70:
            pass
        elif (int(line.split('\t')[0]))<=float(Percent):##
            pass
        else:
            myfile.write(line)
    myfile.close()
    print('Done filter polyA')
    tranidnum=[]
    for line in open('%sflc1-fT.pslx'%out_dir):
        tranidnum.append(line.split('\t')[9])
    tranidnum2=list(set(tranidnum))##
    tranid3=[]#
    tranid3.reverse()
    for i in tranidnum2:
        if tranidnum.count(i)==1:
            pass
        else:
            tranid3.append(i)
            print(i)
    print('filter one transcript complete,file name flc1-fT.pslx')

    uniq=[]#
    most=[]#

    while len(tranid3)>0:
        pop=tranid3.pop()
        print(len(tranid3))
        ran=[]#
        for line in open('%sflc1-fT.pslx'%out_dir):
            if pop in line:
                ran.append(line.split('\t')[13])
        ran2=list(set(ran))
        if len(ran2)==1:
            uniq.append(pop)
    myfile=open('%sflc1-fT-uniqchr.pslx'%out_dir,'w')

    for line in open('%sflc1-fT.pslx'%out_dir):
        if line.split('\t')[9] in uniq:
            myfile.write(line)
    myfile.close()
    print('uniq chr complete,file name flc1-fT-uniqchr.pslx')

    #
    myfile=open('%sflc1-fT-uniqchr-flc.pslx'%out_dir,'w')
    tranidlist=[]
    for line in open('%sflc1-fT-uniqchr.pslx'%out_dir):
        tranidlist.append(line.split('\t')[9])
    tranidlist2=list(set(tranidlist))
    while len(tranidlist2)>0:
        print(len(tranidlist2))
        SS=[]#
        pop=tranidlist2.pop()
        tranid=pop##
        for line in open('%sflc1-fT-uniqchr.pslx'%out_dir):#
            if pop in line:
                A=[int(line.split('\t')[15]),int(line.split('\t')[16])]
                SS.append(A)##
        maxtran=[]#
        for i in SS:
            maxtran.append(i[1]-i[0])
        #print maxtran
        lenSS=len(SS)##
        total=[]
        while len(SS)>0:##
            pop=SS.pop()
            for i in range(pop[0],pop[1]+1):
                total.append(i)     ##
        total2=list(set(total));mubiao=[]
        for i in total2:
            if total.count(i)==lenSS:
                mubiao.append(tranid)
        if (float(len(mubiao)-1))/float(max(maxtran))*100>float(Percent):#
            mubiao2=list(set(mubiao))
            for line in open('%sflc1-fT-uniqchr.pslx'%out_dir):##
                if line.split('\t')[9] in mubiao2:
                    myfile.write(line)
        else:pass
    myfile.close()
    print('flc filter complete flc1-fT-uniqchr-flc.pslx')
    ##

    tranidself=[]
    #for line in open(r'C:/Users/Mu/Desktop/flc1-fT-uniqchr-flc.pslx'):
    for line in open('%sflc1-fT-uniqchr-flc.pslx'%out_dir):
        tranidself.append(line.split('\t')[9])
    tranidself2=list(set(tranidself))##
    # print(fastaFile)
    myfile=open('%sFinal-File-Filter.fasta'%(out_dir),'w')
    while len(tranidself2)>0:
        pop=tranidself2.pop()
        print(len(tranidself2))
        n=[tranidself.count(pop)];n=list(set(n))
        m=[x for x in range(1,n[0]+1)];m.reverse()
        #for line in open(r'C:/Users/Mu/Desktop/flc1-fT-uniqchr-flc.pslx'):
        for line in open('%sflc1-fT-uniqchr-flc.pslx'%out_dir):
            if pop in line:
                mpop=m.pop()
                myfile.write('>%s-%s \n%s'%(line.split('\t')[9],mpop,line.split('\t')[22].replace(',','').upper()))
    myfile.close()
#findflc('Hy1-hlq-flag-filter-fusion.pslx')
def PTpredictFromPB2(out_dir,fastaFile,genome,Percent):
    # #1 
    # print('mapping fasta to genome')
    # os.system('blat %s %s -q=rna -dots=100 -maxIntron=4000 -out=pslx %sflc.pslx'%(genome,fastaFile,out_dir))
    #2
    findflc(out_dir,Percent,fastaFile)
PTpredictFromPB2(args.out,args.fasta,args.genome,args.percent)