#!/usr/bin/env python3
# coding:utf-8
import argparse
import os
import sys
parser = argparse.ArgumentParser(description="predict polycistronic transcript")
parser.add_argument('-gtf', "--gtf", required=True,  help="gft file for predict polycistronic transcript")
parser.add_argument('-ref', "--reference", required=True,  help="reference gtf/gff file with CDS")
parser.add_argument('-genome', "--genome", required=True,  help="genome fasta file")
parser.add_argument('-percent', "--percent", required=True,  type=float,default=50.0,help="Pecentage of polycistronic transcript coverage reference CDS, def=50.0")
parser.add_argument('-o', "--out", required=True,  help="out dir,please note: outdir end with /")
args = parser.parse_args()
out_dir = os.path.dirname(args.out)

if out_dir and not os.path.exists(out_dir):
    os.makedirs(out_dir)

def gff3togtf(x):
    myfile=open(r'tofu-%s.gtf'%x,'w')#
    for line in open(r'%s'%x):
        if 'mRNA\t' in line:
            a=line.split('ID')[0]#
            tranid=(line.split('=')[1]).split('.mrna')[0]
            myfile.write('%sgene_id "%s"; transcript_id "%s";\n'%(a,tranid,tranid))
        if 'exon\t' in line:
            a=line.split('ID')[0]#
            tranid=(line.split('=')[1]).split('.mrna')[0]
            myfile.write('%sgene_id "%s"; transcript_id "%s";\n'%(a,tranid,tranid))
    myfile.close()
def findPT(out_dir):
    tranidlist=[]
    for line in open('%sORF.gff3'%out_dir):
        if 'mRNA\t' in line:
            tranidlist.append((line.split('.p')[0]).split('=')[1])
    tranid=list(set(tranidlist))
    myfile=open(r'%spolycistronic-transcript.txt'%out_dir,'w')
    myfile2=open(r'%spolycistronic-transcript+start-stop.txt'%out_dir,'w')
    while len(tranid)>0:
        print('remain tran %s'%len(tranid))
        SS=[]#
        pop=tranid.pop()
        tran=pop
        for line in open('%sORF.gff3'%out_dir):
            if 'mRNA\t' in line:
                if '%s.p'%pop in line:
                    A=[int(line.split('\t')[3]),int(line.split('\t')[4])]#
                    SS.append(A)
        SS2=SS#
        SSall=[]#
        for i in SS:
            for j in SS2:
                if i==j:
                    pass
                else:
                    SSall.append([i,j])#
        SSover=[]
        for i in SSall:
            for j in range(i[0][0],i[0][1]):
                if j in range(i[1][0],i[1][1]):
                    SSover.append(i)#
        SSover2=[]#
        for i in SSover:
            if i in SSover2:
                pass
            else:
                SSover2.append(i)
        final=[]#
        for i in SSover2:
            if i in SSall:
                SSall.remove(i)
        SSall2=[]
        for i in SSall:
            if i in SSall2:
                pass
            else:
                SSall2.append(i)
        if len(SSall)>0:
            myfile.write('%s\n'%(tran))#
            myfile2.write('%s\t%s\n'%(tran,SSall2))#
def PTgtf(out_dir):
    PBtran=[]
    for line1 in open('%spolycistronic-transcript.txt'%out_dir):
        PBtran.append(line1.strip())
    PBtran=list(set(PBtran))

    myfile1=open('%spolycistronic-transcript.gtf'%out_dir,'w')
    for line2 in open('%sORF.gff3'%out_dir):
        if 'mRNA\t' in line2:
            if (line2.split('.p')[0]).split('=')[1] in PBtran:
                myfile1.write(line2)
        if 'exon\t' in line2:
            if (line2.split('.p')[0]).split('=')[1] in PBtran:
                myfile1.write(line2)
    myfile1.close()
def CompareToFG(x,y,out_dir):
    os.system(r'gffcompare -r %spolycistronic-transcript.gtf -o %spolycistronic-transcript-rev %s'%(out_dir,out_dir,y))
    os.system(r'gffcompare -r %s -o %spolycistronic-transcript %spolycistronic-transcript.gtf'%(y,out_dir,out_dir))

def PTfromGffcompare(x,y,out_dir):
    same=['=','c','k','j','e','o']
    tranid=[]
    for line in open('%spolycistronic-transcript.txt'%out_dir):
        tranid.append(line.strip())
    tranid2=list(set(tranid));tranid3=list(set(tranid2))#
    myfile3=open('%sPT-1-n.txt'%out_dir,'w')
    while len(tranid2)>0:
        print(len(tranid2))
        FGPB=[]
        pop=tranid2.pop()
        for line in open('%spolycistronic-transcript-rev.tracking'%out_dir):
            if '%s.p'%pop in line:
                if line.split('\t')[3] in same:
                    a=[(line.split('\t')[2]).split('.p')[0],((line.split('\t')[4]).split('|')[0]).split(':')[1]]#
                    b='%s\t%s'%(a[0],a[1])#
                    FGPB.append(b)
        FGPB2=list(set(FGPB))
        if len(FGPB)<2:
            pass
        elif FGPB2[0].split('\t')[0]==FGPB2[1].split('\t')[0]:
            for i in FGPB2:
                myfile3.write('%s\n'%i)

    while len (tranid3)>0:
        print(len(tranid3))
        FGPB=[]
        pop=tranid3.pop()
        for line in open('%spolycistronic-transcript.tracking'%out_dir):
            if '%s.p'%pop in line:
                if line.split('\t')[3] in same:
                    a=[(((line.split('\t')[4]).split('|')[0]).split(':')[1]).split('.p')[0],(line.split('\t')[2]).split('|')[0]]
                    b='%s\t%s'%(a[0],a[1])
                    FGPB.append(b)
        FGPB2=list(set(FGPB))
        if len(FGPB2)<2:
            pass
        elif FGPB2[0].split('\t')[0]==FGPB2[1].split('\t')[0]:
            for i in FGPB2:
                myfile3.write('%s\n'%i)
    myfile3.close()
    PTtran=[];PTgene=[]
    for line in open('%sPT-1-n.txt'%out_dir):
        PTtran.append(line.split('\t')[0])
        PTgene.append(line.split('.')[1])
    PTtran=list(set(PTtran));PTgene=list(set(PTgene))
    # print 'FGmapping tran\t%s'%(len(PTtran))
    # print 'FGmapping gene\t%s'%(len(PTgene))
def PBinfasta(x,y,z,a,b):

    PTtran=[]
    for line in open(r'%s-PT-1-n.txt'%(x.split('.')[0])):
        if line.split('\t')[0] in PTtran:
            pass
        else:
            PTtran.append(line.split('\t')[0])
    for line in open(r'%s-PT-1-n.txt'%(y.split('.')[0])):
        if line.split('\t')[0] in PTtran:
            pass
        else:
            PTtran.append(line.split('\t')[0])
    PTreads=[]
    for line in open(r'%s'%z):
        if line.split('\t')[0] in PTtran:
            PB=((line.strip()).replace(',','\t')).split('\t')[1:]
            for i in PB:
                PTreads.append(i)

    myfile4=open(r'PT-%s'%a,'w')
    myfile5=open(r'filter-PT-%s'%a,'w')
    title=[]
    fasta=[]
    for line in open(r'%s'%a):
        if '>' in line:
            title.append(line)
        else:
            fasta.append(line)
    while len(title)>0:
        print(len(title))
        titlepop=title.pop()
        fastapop=fasta.pop()
        if (titlepop.split(' ')[0]).split('>')[1] in PTreads:
            myfile4.write(titlepop)
            myfile4.write(fastapop)
        else:
            myfile5.write(titlepop)
            myfile5.write(fastapop)
    myfile4.close()
    myfile5.close()

    myfile6=open(r'filter-PT-%s'%b,'w')
    myfile7=open(r'PT-%s'%b,'w')
    for line in open(r'%s'%b):
        if line.split('\t')[0] in PTreads:
            myfile7.write(line)
        else:
            myfile6.write(line)
    myfile6.close()
    myfile7.close()
def LongORFinPT(out_dir):
    Already=[]
    same=['=','c','e','j','k','o']
    FGRAMlongORF=open('%slongORF-PT-1-n.txt'%out_dir,'w')
    print('%slongORF-PT-1-n.txt'%out_dir)
    FGRAMlongORF.write('PBtran\tToFu\tlongORFPT\tchr\tstrand\tleft\tright\tlen\n')
    for line5 in open('%sPT-1-n.txt'%out_dir):#
        PBtran=line5.split('\t')[0]
        FGRAMgene=(line5.split('\t')[1]).strip()
        PTorf=[]#
        for line6 in open('%spolycistronic-transcript.tracking'%out_dir):
            if line6.split('\t')[3] in same:
                if FGRAMgene in line6:
                    if '%s.'%PBtran in line6:
                        orf=((line6.split('\t')[4]).split('|')[1]).split('.mrna')[0]
                        if orf in PTorf:
                            pass
                        else:
                            PTorf.append(orf)
        for line7 in open('%spolycistronic-transcript-rev.tracking'%out_dir):
            if line7.split('\t')[3] in same:
                if FGRAMgene in line7:
                    if '%s.'%PBtran in line7:
                        orf=(line7.split('\t')[2]).split('.path')[0]
                        if orf in PTorf:
                            pass
                        else:
                            PTorf.append(orf)
        strand=[];CHR=[];left=[];right=[];orfLen=[]
        for i in PTorf:
            for line8 in open('%spolycistronic-transcript.gtf'%out_dir):
                if 'mRNA\t' in line8:
                    if '%s.'%i in line8:
                        #print i,line8
                        CHR.append(line8.split('\t')[0])#
                        strand.append(line8.split('\t')[6])#
                        left.append(int(line8.split('\t')[3]))#
                        right.append(int(line8.split('\t')[4]))#
                        orfLen.append(int(line8.split('\t')[4])-int(line8.split('\t')[3]))
        if len(orfLen)==1:
            #pass#print PBtran,FGSGgene,PTorf,orfLen,orfLen[0],PTorf[0]
            #print PBtran,FGSGgene,PTorf[0],CHR[0],strand[0],left[0],right[0],orfLen[0]
            if '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(PBtran,FGRAMgene,PTorf[0],CHR[0],strand[0],left[0],right[0],orfLen[0]) not in Already:
                FGRAMlongORF.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(PBtran,FGRAMgene,PTorf[0],CHR[0],strand[0],left[0],right[0],orfLen[0]))
                Already.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(PBtran,FGRAMgene,PTorf[0],CHR[0],strand[0],left[0],right[0],orfLen[0]))
        else:
            #print PBtran,FGSGgene,PTorf[orfLen.index(max(orfLen))],CHR[0],strand[0],left[0],right[0],max(orfLen)
            if '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(PBtran,FGRAMgene,PTorf[orfLen.index(max(orfLen))],CHR[0],strand[0],left[orfLen.index(max(orfLen))],right[orfLen.index(max(orfLen))],max(orfLen)) not in Already:
                FGRAMlongORF.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(PBtran,FGRAMgene,PTorf[orfLen.index(max(orfLen))],CHR[0],strand[0],left[orfLen.index(max(orfLen))],right[orfLen.index(max(orfLen))],max(orfLen)))
                Already.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(PBtran,FGRAMgene,PTorf[orfLen.index(max(orfLen))],CHR[0],strand[0],left[orfLen.index(max(orfLen))],right[orfLen.index(max(orfLen))],max(orfLen)))
    FGRAMlongORF.close()
    print('longORF complete!')
#LongORFinPT('polycistronic-transcript.gtf','polycistronic-transcript.tracking','polycistronic-transcript-rev.tracking')
def longORFmapSameFG(out_dir):#
    myfile1=open('%sover2FG.txt'%out_dir,'w')
    myfile2=open('%sfilter2FGorf.txt'%out_dir,'w')
    PT=[]#
    for line1 in open('%slongORF-PT-1-n.txt'%out_dir):
        if 'PB.' in line1:
            if line1.split('\t')[0] in PT:
                pass
            else:
                PT.append(line1.split('\t')[0])
    ORFoverTwoFG=[]
    for i in PT:
        PTorf=[]#
        for line2 in open('%slongORF-PT-1-n.txt'%out_dir):
            if '%s\t'%i in line2:
                PTorf.append(line2.split('\t')[2])
        PTorf2=list(set(PTorf))
        if len(PTorf)>len(PTorf2):#
            ORFoverTwoFG.append(i)
    for line3 in open('%slongORF-PT-1-n.txt'%out_dir):
        if line3.split('\t')[0] in ORFoverTwoFG:
            print(line3)
            myfile1.write(line3)
        else:
            myfile2.write(line3)
    myfile1.close()
    myfile2.close()
def longORFcoverFGpercent(x,out_dir,Percent):#
    myfile1=open('%sOver-percent.txt'%out_dir,'w')#
    myfile2=open('%sLess-percent.txt'%out_dir,'w')#
    myfile1.write('PBtran\tRef\tlongORFPT\tchr\tstrand\tleft\tright\tlen\tCoverage\n')
    myfile2.write('PBtran\tRef\tlongORFPT\tchr\tstrand\tleft\tright\tlen\tCoverage\n')
    less50=[]#
    for linePTFGSG in open('%sfilter2FGorf.txt'%out_dir):
        #print('%sfilter2FGorf.txt'%out_dir)
        if 'PB.' in linePTFGSG:
            PTorf= linePTFGSG.split('\t')[2];FGSG=linePTFGSG.split('\t')[1];PTleft=int(linePTFGSG.split('\t')[5]);PTright=int(linePTFGSG.split('\t')[6]);FGSGsite=[]#存放orf及基因的位点信息
            for lineFGSGgff3 in open(x):
                if FGSG in lineFGSGgff3:
                    if 'mRNA\t' in lineFGSGgff3:
                        FGSGsite.append(int(lineFGSGgff3.split('\t')[3]))
                        FGSGsite.append(int(lineFGSGgff3.split('\t')[4]))
            FGSGleft=min(FGSGsite);FGSGright=max(FGSGsite)
            FGSGlen=FGSGright-FGSGleft
            cover=[]
            FGSGrange=list(range(FGSGleft,FGSGright+1))#
            PTrange  =list(range(PTleft,PTright+1))
            for i in PTrange:
                if i in FGSGrange:
                    cover.append(i)#
            percent=float(len(cover))/float(FGSGlen)*100
            if percent<float(Percent):
                #print PTorf,FGSG,percent#,PTleft,PTright,FGSGleft,FGSGright,len(cover),FGSGlen
                myfile2.write('%s\t%s\n'%(linePTFGSG.strip(),int(percent)))
                less50.append(linePTFGSG.split('\t')[0])
            else:
                myfile1.write('%s\t%s\n'%(linePTFGSG.strip(),int(percent)))
    myfile1.close()
    myfile2.close()
    myfile3=open('%sTlesspercent.txt'%out_dir,'a')#
    myfile4=open('%sTMorepercent.txt'%out_dir,'w')#
    for line3 in open('%sOver-percent.txt'%out_dir):
        if line3.split('\t')[0] in less50:
            myfile3.write(line3)
        else:
            myfile4.write(line3)
    myfile3.close()
    myfile4.close()
#longORFcoverFGpercent('filter2FGorf-longORF-gff-FGSG-PT-1-n.txt','gff-FGSG.gff3')
def orfClass(out_dir):#
    PTtran=[]#
    for line1 in open('%sOver-percent.txt'%out_dir):
        if 'PB.' in line1:
            if line1.split('\t')[0] in PTtran:
                pass
            else:
                PTtran.append(line1.split('\t')[0])
    myfile1=open('%sFinal-predict-result.txt'%out_dir,'w')
    myfile1.write('PBtran\tRef\tlongORFPT\tchr\tstrand\tleft\tright\tlen\tCoverage\tORFnum\n')
    for i in PTtran:#
        linePT=[];orflist=[];strand=[];left=[]
        for line2 in open('%sOver-percent.txt'%out_dir):
            if '%s\t'%i in line2:
                orflist.append(line2.split('\t')[2])
                strand.append(line2.split('\t')[4])
                left.append(int(line2.split('\t')[5]))
                linePT.append(line2)
        #ORFclass=sorted(zip(orflist,left),key=lambda x:(-x[1]))#
        if strand[0]=='+':#
            ORFclass=sorted(zip(orflist,left),key=lambda x:(x[1]))
            n=1
            for j in ORFclass:
                for line3 in open('%sOver-percent.txt'%out_dir):
                    if '%s\t'%j[0] in line3:
                        #print line3.strip(),'ORF%s'%n
                        myfile1.write('%s\tORF%s\n'%(line3.strip(),n))
                        n+=1
        if strand[0]=='-':
            ORFclass=sorted(zip(orflist,left),key=lambda x:(-x[1]))
            n=1
            for k in ORFclass:
                for line4 in open('%sOver-percent.txt'%out_dir):
                    if '%s\t'%k[0] in line4:
                        #print line4.strip(),'ORF%s'%n
                        myfile1.write('%s\tORF%s\n'%(line4.strip(),n))
                        n+=1
    myfile1.close()
def PTpredictFromPB2(out_dir,x,y,genome,percent):
    #1 
    print('extrace fasta for gtf')
    os.system('gffread -w %sgtf.fasta -g %s %s'%(out_dir,genome,x))
    # #2 
    print('ORF predict')
    os.system('TransDecoder.LongOrfs -S -t %sgtf.fasta -O %s'%(out_dir,out_dir))
    os.system('gmap_build -D %sdio -d dio %s'%(out_dir,genome))
    #3 
    os.system('gmap -D  %sdio -d dio -f 2 -n 1 -t 15 --no-chimeras --min-intronlength 20 --max-intronlength-middle 4000 --max-intronlength-ends 4000 -z sense_force %slongest_orfs.cds >%sORF.gff3'%(out_dir,out_dir,out_dir))
    #4 
    findPT(out_dir)
    #6 
    PTgtf(out_dir)
    #7 
    CompareToFG(x,y,out_dir)
    #8 
    PTfromGffcompare(x,y,out_dir)
    #9
    LongORFinPT(out_dir)
    #10
    longORFmapSameFG(out_dir)
    #11
    longORFcoverFGpercent(y,out_dir,percent)
    #12
    orfClass(out_dir)

#PTpredictFromPB2('ISO-fusion.collapsed.gtf','gff-FGSG.gff3','gff-FGRAMPH1.gff3')
PTpredictFromPB2(args.out,args.gtf,args.reference,args.genome,args.percent)