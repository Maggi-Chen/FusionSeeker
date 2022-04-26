import pysam
import os

def annotate_segment(chrom,start,end):
	allgenes=geneinfo[chrom]
	pin=int(len(allgenes)/2)
	pinstart=0;pinend=len(allgenes)
	ovlpgene=[]
	notfound=True
	while notfound:
		gene=allgenes[pin]
		ovlplen=min(gene[1],end)-max(gene[0],start)

		if ovlplen>=0:
			notfound=False; break
		if start>gene[1]:
			newpin=int((pin+pinend)/2)
			if newpin==pin:
				break
			else:
				pinstart=pin;pin=newpin
				continue
		if end<gene[0]:
			newpin=int((pin+pinstart)/2)
			if newpin==pin:
				break
			else:
				pinend=pin;pin=newpin
				continue

	for gene in allgenes[max(0,pin-10):pin+10]:
		ovlplen=min(gene[1],end)-max(gene[0],start)
		if ovlplen>50 or ovlplen>0.5*(end-start):
			ovlpgene+=[gene[2]]
	return ovlpgene

def getcigarsegment(cigartuple):
	exonsegments=[]
	lastlength=0
	for pair in cigartuple:
		if pair[0] in [3,7,8]:
			exonsegments+=[lastlength]
			lastlength=0
		if pair[0]==0:
			lastlength+=pair[1]
		if pair[0] in [4,5]:
			if lastlength!=0:
				exonsegments+=[lastlength]
				lastlength=0
		if pair[0] in [1,2]:
			continue

	if lastlength!=0:
		exonsegments+=[lastlength]
	return exonsegments


def getsimplecigar(cigartuple):
	leftclip=0
	rightclip=0
	if cigartuple[0][0] in [4,5]:
		leftclip=cigartuple[0][1]
	if cigartuple[-1][0] in [4,5]:
		rightclip=cigartuple[-1][1]
	return [leftclip,0,rightclip]

def simplify_cigar(cigartuple):
	cigarlist=[]
	nummatch=0
	for c in cigartuple:
		if c[0] in [0,7,8]:
			nummatch+=c[1]
		if c[0] in [3]:
			cigarlist+=[str(nummatch)+'M']
			nummatch=0
		if c[0] in [2]:
			continue
		if c[0] in [1]:
			cigarlist+=[str(nummatch)+'M',str(c[1])+'I']
			nummatch=0
		if c[0] in [4,5]:
			if nummatch!=0:
				cigarlist+=[str(nummatch)+'M']
				nummatch=0
	if nummatch!=0:
		cigarlist+=[str(nummatch)+'M']
		nummatch=0
	cigarlist=','.join(cigarlist)
	return cigarlist

def detect_withinread(read,recordread,recordseq):
	candi=[]
	readinfo=[]
	alignpair=read.get_aligned_pairs()
	alignpair=[c for c in alignpair if c[1]!=None and c[0]!=None]
	chrom=read.reference_name
	cigartuple=read.cigartuples
	cigarinfo=getsimplecigar(cigartuple)
	cigarinfo[1]=read.infer_read_length()-cigarinfo[0]-cigarinfo[2]
	exonsegments=getcigarsegment(cigartuple)
	exonsegments=[exonlen for exonlen in exonsegments if exonlen>0]
	nmrate=read.get_tag('NM')*1.0/read.query_alignment_length
	simplecigar=simplify_cigar(cigartuple)

	if read.flag in [0,2048]:
		readstrand='+'
	elif read.flag in [16,2064]:
		readstrand='-'
	else:
		readstrand='.' 


	if recordseq==False or read.is_secondary or read.is_supplementary:
		readsequence=''
		readquality=''
	else:
		readsequence=read.query_sequence
		readqua=read.query_qualities
		readquality=''
		try:
			for c in readqua:
				readquality+=chr(c+33)
		except:
			pass

	start=alignpair[0][1]
	end=alignpair[-1][1]

	if len(exonsegments)==1: # only one exon
		geneboth=annotate_segment(chrom,start,end)
		if recordread:
			readinfo=[chrom,start,end,read.query_name,readstrand,cigarinfo,read.mapping_quality,geneboth,geneboth,[[start,end,exonsegments[0],geneboth]],readsequence,readquality,nmrate,simplecigar,read]

		return [[],readinfo]

	alignsegment=[]
	accumalignlen=0
	for exonlen in exonsegments:
		alignsegment+=[[alignpair[accumalignlen][1],alignpair[accumalignlen+exonlen-1][1],exonlen]]
		accumalignlen+=exonlen
	segmentgenes={}
	for i in range(len(alignsegment)):
		segmentgenes[i]=[]


	for i in range(len(alignsegment)):
		segmentgenes[i]=annotate_segment(chrom,alignsegment[i][0],alignsegment[i][1])

	annotatedalignsegment=[]

	for i in range(len(alignsegment)):
		segment=alignsegment[i]
		annotatedalignsegment+=[segment+[segmentgenes[i]]]


	if recordread:
		readinfo=[chrom,start,end,read.query_name,readstrand,cigarinfo,read.mapping_quality,annotatedalignsegment[0][3],annotatedalignsegment[-1][3],annotatedalignsegment,readsequence,readquality,nmrate,simplecigar,read]

	return [[],readinfo]



def get_raw_signal(bampath,outpath,chrom,recordseq=True):
	bam=pysam.AlignmentFile(bampath,'rb')
	allread=bam.fetch(chrom)
	fusion=[]
	global numberwithin
	numberwithin=0
	split_readinfo=[]
	for read in allread:
		if read.is_secondary:
			continue
		recordread=False
		if read.has_tag('SA'):
			suppalign=read.get_tag('SA').split(';')[:-1]
			suppalign=[c for c in suppalign if int(c.split(',')[3].split('M')[0][::-1].split('S')[0][::-1])>=100]
			if len(suppalign)>0:
				recordread=True
				[withfusion,readinfo]=detect_withinread(read,recordread,recordseq)
				split_readinfo+=[readinfo]
	f=open(outpath+'raw_signal/record_read_'+chrom,'w')
	for c in split_readinfo:
		exoninfo=[]
		for mm in c[9]:
			exoninfo+=[str(mm[0])+','+str(mm[1])+','+str(mm[2])+','+':'.join(mm[3])]
		exoninfo=';'.join(exoninfo)
		cigarinfo=[str(mm) for mm in c[5]]
		f.write(c[0]+'\t'+str(c[1])+'\t'+str(c[2])+'\t'+c[3]+'\t'+c[4]+'\t'+','.join(cigarinfo)+'\t'+str(c[6])+'\t'+str(c[12])+'\t'+','.join(c[7])+'\t'+','.join(c[8])+'\t'+exoninfo+'\t'+c[13]+'\t'+c[10]+'\t'+c[11]+'\n')
	f.close()
	print(chrom)
	return 0


def sortsplitread(a):
	return a.split('\t')[3]


def get_splitgene(read,side):
	genename=''
	maplen=0
	exonlist=read[9]
	if side =='right':
		exonlist=exonlist[::-1]
	allcandigene=exonlist[0][3]
	if len(allcandigene)==1:
		genename=allcandigene[0]
		if genename=='':
			return [genename,0]
		for exon in exonlist:
			if genename in exon[3]:
				maplen+=exon[2]
		return [genename,maplen]

	genemaplen={}
	for candigene in allcandigene:
		genemaplen[candigene]=0
	for exon in exonlist:
		for candigene in allcandigene:
			if candigene in exon[3]:
				genemaplen[candigene]+=exon[2]
	longestgene=''
	longestlen=0
	for candigene in allcandigene:
		if genemaplen[candigene]>longestlen:
			longestgene=candigene
			longestlen=genemaplen[candigene]
	genename=longestgene
	maplen=genemaplen[genename]
	return [genename,maplen]


def remove_ovlp_exon(oldreadinfo,side,removelength):
	readinfo=[]
	exonlist=oldreadinfo[9]
	deletedlen=0

	inspos=oldreadinfo[11]
	if side=='right':
		inspos=inspos[::-1]
	includeins=0
	toremove=removelength
	for c in inspos:
		if toremove<=0:
			break
		if 'M' in c:
			toremove-=int(c[:-1])
			includeins+=int(c[:-1])
		if 'I' in c:
			toremove-=int(c[:-1])
	removelength=includeins
	
	if side=='right':
		while len(exonlist)>0 and removelength>0:
			if removelength>=exonlist[-1][2]-10 :
				removelength-=exonlist[-1][2]
				deletedlen+=exonlist[-1][2]
				exonlist=exonlist[:-1]
				continue
			exonlist[-1]=[exonlist[-1][0],exonlist[-1][1]-removelength,exonlist[-1][2]-removelength,exonlist[-1][3]]
			deletedlen+=removelength
			removelength=0
		if exonlist==[]:
			return []
		readinfo=[oldreadinfo[0],oldreadinfo[1],exonlist[-1][1],oldreadinfo[3],oldreadinfo[4],oldreadinfo[5],oldreadinfo[6],oldreadinfo[7],exonlist[-1][3],exonlist,oldreadinfo[10]]
	else:
		while len(exonlist)>0 and removelength>0:
			if removelength>=exonlist[0][2]-10:
				removelength-=exonlist[0][2]
				deletedlen+=exonlist[0][2]
				exonlist=exonlist[1:]
				continue
			exonlist[0]=[exonlist[0][0]+removelength,exonlist[0][1],exonlist[0][2]-removelength,exonlist[0][3]]
			deletedlen+=removelength
			removelength=0
		if exonlist==[]:
			return []
		readinfo=[oldreadinfo[0],exonlist[0][0],oldreadinfo[2],oldreadinfo[3],oldreadinfo[4],oldreadinfo[5],oldreadinfo[6],exonlist[0][3],oldreadinfo[8],exonlist,oldreadinfo[10]]
	return readinfo



def get_fusion_readpair(read1,read2):
	read1=list(read1)
	read2=list(read2)
	read1[5]=list(read1[5])
	read2[5]=list(read2[5])

	candifusion=[]
	if read1[4]=='-':
		temp=read1[5][0]; read1[5][0]=read1[5][2]; read1[5][2]=temp
	if read2[4]=='-':
		temp=read2[5][0];read2[5][0]=read2[5][2];read2[5][2]=temp

	if read1[5][0]<read2[5][0]:
		leftread=read1
		rightread=read2
	else:
		leftread=read2
		rightread=read1

	gapsize=rightread[5][0]-leftread[5][0]-leftread[5][1]

	
	if 0-gapsize > min(100,leftread[5][1]*0.5,rightread[5][1]*0.5):
		return []
	
	if gapsize<-100:
		
		removelength=0-gapsize
		if leftread[10]>rightread[10]:
			if leftread[4]=='+':
				leftread=remove_ovlp_exon(leftread,'right',removelength)
			else:
				leftread=remove_ovlp_exon(leftread,'left',removelength)
		else:
			if rightread[4]=='+':
				rightread=remove_ovlp_exon(rightread,'left',removelength)
			else:
				rightread=remove_ovlp_exon(rightread,'right',removelength)

	if leftread==[] or rightread==[]:
		return []

	if leftread[4]=='+':
		chrom1=leftread[0]
		bp1=leftread[2]
		[gene1,maplen1]=get_splitgene(leftread,'right')
	else:
		chrom1=leftread[0]
		bp1=leftread[1]
		[gene1,maplen1]=get_splitgene(leftread,'left')


	if rightread[4]=='+':
		chrom2=rightread[0]
		bp2=rightread[1]
		[gene2,maplen2]=get_splitgene(rightread,'left')
	else:
		chrom2=rightread[0]
		bp2=rightread[2]
		[gene2,maplen2]=get_splitgene(rightread,'right')

	if gene1 =='' or gene2=='' or gene1==gene2 or maplen1<100 or maplen2<100:
		return []


	if chrom1<chrom2:
		candifusion=[gene1,gene2,'splitread',chrom1,bp1,chrom2,bp2,read1[3],[leftread[6],rightread[6]],maplen1,maplen2,gapsize]
	elif chrom2<chrom1:
		candifusion=[gene2,gene1,'splitread',chrom2,bp2,chrom1,bp1,read1[3],[rightread[6],leftread[6]],maplen2,maplen1,gapsize]
	else:
		if bp1<bp2:
			candifusion=[gene1,gene2,'splitread',chrom1,bp1,chrom2,bp2,read1[3],[leftread[6],rightread[6]],maplen1,maplen2,gapsize]
		else:
			candifusion=[gene2,gene1,'splitread',chrom2,bp2,chrom1,bp1,read1[3],[rightread[6],leftread[6]],maplen2,maplen1,gapsize]


	return [candifusion]



def get_fusion_sameread(sameread):
	candi=[]
	for i in range(len(sameread)-1):
		read1=sameread[i]
		for read2 in sameread[i+1:]:
			candi+=get_fusion_readpair(read1,read2)
	return candi


def detect_from_split(outpath,good):
	allsplitread=[]
	for chrom in good:
		allsplitread+=open(outpath+'raw_signal/record_read_'+chrom,'r').read().split('\n')[:-1]
	allsplitread.sort(key=sortsplitread)
	fusion=[]
	lastname=''
	sameread=[]

	allsplitreadinfo=[]
	for read in allsplitread:
		read=read.split('\t')
		exoninfo=[]
		for exon in read[10].split(';'):
			exon=exon.split(',')
			exoninfo+=[[int(exon[0]),int(exon[1]),int(exon[2]),exon[3].split(':')]]
		allsplitreadinfo+=[[read[0],int(read[1]),int(read[2]),read[3],read[4],[int(m) for m in read[5].split(',')],int(read[6]),read[8].split(','),read[9].split(','),exoninfo,float(read[7]),read[11].split(',')]]

	allsplitread=''

	for read in allsplitreadinfo:
		if read[3] ==lastname:
			sameread+=[read]
			continue
		if lastname=='':
			lastname=read[3]
			sameread=[read]
			continue
		if len(sameread)>1:
			splitfusion=get_fusion_sameread(sameread)
			fusion+=splitfusion
		sameread=[read]
		lastname=read[3]

	if len(sameread)>1:
		splitfusion=get_fusion_sameread(sameread)
		fusion+=splitfusion


	f=open(outpath+'rawsignal.txt','w')
	for c in fusion:
#		gene1,gene2,'splitread',chrom1,bp1,chrom2,bp2,read1[3],[leftread[6],rightread[6]],maplen1,maplen2,gapsize
		mapq=[str(cc) for cc in c[8]]
		f.write(c[0]+'\t'+c[1]+'\t'+c[2]+'\t'+c[3]+'\t'+str(c[4])+'\t'+c[5]+'\t'+str(c[6])+'\t'+c[7]+'\t'+','.join(mapq)+'\t'+str(c[9])+'\t'+str(c[10])+'\t'+str(c[11])+'\n')
	f.close()
	return 0


