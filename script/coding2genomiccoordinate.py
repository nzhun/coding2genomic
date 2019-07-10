#!/usr/bin/env python
import os   ## standard library 
import os.path
import sys  ## standard library
import re
import getopt

## read variant table 
## read gene info table 
## get the positon of a variant 
## extract the sequence using chr:start-end
## output chr start end ref alt
def Indexnumber(word,id):
    i=0
    p=0
    mat=re.compile(r'[0-9]')
    if word[id]=='-':
        id=id+1;
    for i in range(id,len(word)):
        if mat.match(word[i])== None:
             break
        else:
            p=1
    if p==0:
        return -1
    else:
        return i
def searchinterval (cn,gene):
    global dict_cstart 
    starts=sorted(dict_cstart.get(gene))
    ends=sorted(dict_cend.get(gene))
#    print starts
    for j in range(0,len(starts)):
#        print str(starts[j])+'\t'+str(cn)
        if starts[j] > cn:
            break;
    if cn > starts[len(starts)-1]:
        return len(starts)-1
    if j ==0:
        return j
    else:
        return j-1
  
def reverse_seq(word):
    newword=""
    global dict_tab
    for  i in range(len(word)-1,-1):
        newword=newtab+dict_tab.get(word[i])
    return newword
    
def decode_NC(word,gene):
    global dict_cstart
    global dict_gstart
    global dict_tstart
    global dict_chr
    global ref_dict
    refstr=ref_dict.get(gene)
    strand=dict_strand.get(gene)
    starts=sorted(dict_cstart.get(gene))
    gstarts=sorted(dict_gstart.get(gene))
    tstart=dict_tstart.get(gene)
    num=[]
    id1=Indexnumber(word,2)
    n1=int(word[2:id1])
    j=searchinterval(n1,gene)
   # print starts
  #  print str(n1)+'\t'+str(starts[j])
    if n1 >0:
        coor1=gstarts[j]+strand*(n1-starts[j])
    else:
        coor1=gstarts[j]+strand*n1
    alt=""
    ref=""
    i=id1
    while i <len(word):
        id=i+1;
        if word[i]=='+' or word[i]=='-' :
            idx=Indexnumber(word,i+1)
            n=int(word[i+1:idx])
            if word[i]=='+':
                coor1=coor1+int(strand)*n
            elif word[i]=='-':
                coor1=coor1-int(strand)*n
            i=idx
        elif word[i]=="_":
            num.append(coor1)
            id1=Indexnumber(word,i+1)
            n1=int(word[i+1:id1])
            j=searchinterval(n1,gene)
            if n1 >0:
                coor1=gstarts[j]+strand*(n1-starts[j])
            else:
                coor1=gstarts[j]+strand*n1
            i=id1
        elif word[i] in {'d','i','A','C','G','T'}:
            break
    num.append(coor1)
    if len(num)==1:
        num.append(coor1)
    x=word.find("ins")
    y=word.find("dup")
    z=word.find("del")
    u=word.find("delins")
    s=word.find(">")
    if u>-1:
         ref=refstr[num[0]-tstart:num[1]-tstart+1]
         alt=word[u+6:len(word)]
    elif x >-1 or y >-1:
        ref=refstr[num[0]-tstart]
        if y > x:
            x=y
        alt=ref+word[x+3:len(word)]
    elif z >-1:
        num[0]=num[0]-1
        ref=refstr[num[0]-tstart: num[1]-tstart+1]
        alt=refstr[num[0]-tstart]
    elif s>-1:
        alt=word[s+1]
        ref=word[s-1]
        if ref!=refstr[coor1-tstart]:
            print "warning! "+word
    return dict_chr.get(gene)+'\t'+str(num[0])+'\t'+str(num[1])+'\t'+ref+'\t'+alt
          
def loadFasta(ffasta):
    global ref_dict;
    try:
        with open(ffasta, "r") as fp:
            line=fp.readline()
            refstr=""
            tid=""
            while line:
                if line.startswith('>'):
                    if refstr != "":
                        ref_dict[tid]=refstr
                        refstr=""
                        tid=""
                    sets=line.split("|")
                    tid=sets[0][1:len(sets[0])]
                else:
                    
                    refstr=refstr+line.strip().rstrip()
                line=fp.readline()
    except OSError as err:
            print("OS error: {0}".format(err))
    finally:
        fp.close()

def loadGene(gfile):
    global dict_gstart
    global dict_gend
    global dict_cstart
    global dict_cend
    global dict_chr
    global dict_strand
    global dict_tstart
    try:
       # print os.path.isfile(gfile)
        with open(gfile, "r") as fp:
            line=fp.readline()
            #for line in fp:
            while line:
                if line.startswith("Gene"):
                    line=fp.readline()
                    continue;
                sets=line.split('\t')
                strand=sets[9]
                gene=sets[10]
                chrm=sets[4]
                tstart=sets[7]
                dict_tstart[gene]=int(tstart)
                if sets[11] != "":
                    gstart=sets[11]
                    gend=sets[12]
                    if strand == -1 :
                        gstart=sets[12]
                        gend=sets[11]
                    cstart=sets[15]
                    cend=sets[16]
                    dict_chr[gene]=chrm
                    dict_strand[gene]=int(strand)
                    if gene in dict_gstart.keys():
                        dict_gstart.get(gene).add(int(gstart))
                    else:
                        dict_gstart[gene]={int(gstart)}
                    if gene in dict_gend.keys():
                        dict_gend.get(gene).add(int(gend))
                    else:
                        dict_gend[gene]={int(gend)}
                    if gene in dict_cend.keys():
                        dict_cend.get(gene).add(int(cend))
                    else:
                        dict_cend[gene]={int(cend)}
                    if gene in dict_cstart.keys():
                        dict_cstart.get(gene).add(int(cstart))
                    else:
                        dict_cstart[gene]={int(cstart)}
                line=fp.readline()
            fp.close()
    except OSError as err:
        print("OS error: {0}".format(err))

def main(gfile,vfile,ffasta,fout):
    loadGene(gfile)
    loadFasta(ffasta)
    outf=open(fout,'w')
    outf.write("GeneName"+'\t'+"Start"+'\t'+"End"+'\t'+"REF"+'\t'+"ALT"+'\t'+"rawref"+'\t'+"sequence"+'\t')
#    vgene="BMPR2"
 #   var="c.-127936_418+7067del"
 #   print decode_NC(var,vgene)
  #  exit()
    try:
        fv=open(vfile,"r")
        line=fv.readline()
        while line:
            line=line.rstrip()
            if line.startswith("Gene"):
                outf.write(line+'\n')
                line=fv.readline()
                continue
            sets=line.split("\t")
            vgene=sets[0]
            vtype=sets[2]
            var=sets[3]
            coordinate=0;
            coordinate2=0;
            tstart=dict_tstart.get(vgene)
            m=var.find('?')
            if vtype == "Deletion" or m > -1:
                outf.write(line+'\n')
                line=fv.readline()
                next;
            if vtype != "Deletion" and m ==-1:
                snc=decode_NC(var,vgene)
                outf.write(snc+'\t'+line+'\n')
                line=fv.readline()
    finally:
        fv.close()
        
dict_tab={'A':'T','C':'G','G':'C','T':'A'}
dict_gstart={}
dict_gend={}
dict_cstart={}
dict_cend={}
dict_chr={}
dict_strand={}
dict_tstart={}
ref_dict={}
#gfile='/Users/nazhu/Dropbox (CGC)/methods/VT_multiple_dimension/data/GeneInfo.txt'
#vfile='/Users/nazhu/Dropbox (CGC)/methods/VT_multiple_dimension/data/variantinMachado.txt'
#fref='/Users/nazhu/Dropbox (CGC)/methods/VT_multiple_dimension/data/PAHGene.fasta.txt'
myopts,args=getopt.getopt(sys.argv[1:],"g:v:r:o:h")
for s, a in myopts:
    print s+' '+a
    if s=='-g':
        gfile=a
    elif s=='-v':
        vfile=a
    elif s=='-r':
        fref=a
    elif s=='-o':
        fout=a
    elif s=='-h':
        usage()
    else:
        print ("Usage: %s -g geneInfo -v variantInfo -r sequecefile -o output" )
print gfile
print vfile
print "ref "+fref
print "fout "+fout
main(gfile,vfile,fref,fout)
