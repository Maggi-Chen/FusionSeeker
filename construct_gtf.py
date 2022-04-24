
def sortgeneinfo(a):
	return a[0]


def write_geneinfo():
	f=open('testlast.ds','w')
	for chrom in geneinfo:
		for c in geneinfo[chrom]:
			f.write(chrom+'\t'+str(c[0])+'\t'+str(c[1])+'\t'+c[2]+'\t'+c[3]+'\t')
			for m in c[4]:
				f.write(str(m[0])+','+str(m[1])+','+str(m[2])+';')
			f.write('\n')
	f.close()
	
def readgeneinfo(geneinfofilename):
	allgene=open(geneinfofilename,'r').read().split('\n')[:-1]
	global geneinfo
	geneinfo={}
	for c in allgene:
		c=c.split('\t')
		eventinfo=[int(c[1]),int(c[2]),c[3],c[4]]
		exoninfo=[]
		for m in c[5].split(';')[:-1]:
			m=m.split(',')
			exoninfo+=[[m[0],int(m[1]),int(m[2])]]
		eventinfo+=[exoninfo]
		if c[0] not in geneinfo:
			geneinfo[c[0]]=[eventinfo]
		else:
			geneinfo[c[0]]+=[eventinfo]

	return geneinfo

def create(gtfinfo,goodchrom,usegeneid,writeds=False):
	global geneinfo
	geneinfo={}
	for chrom in goodchrom:
		geneinfo[chrom]=[]
		allgt=[c for c in gtfinfo if c.split('\t')[0]==chrom and c.split('\t')[2]!='transcript']
		lastgene=''
		genpos=''
		exonpos=[]
		for event in allgt:
			if usegeneid:
				genename=event.split('gene_id "')[1].split('"')[0]
			else:
				try:
					genename=event.split('gene_name "')[1].split('"')[0]
				except:
					continue
			if genename==lastgene:
				exonpos+=[[event.split('\t')[2],int(event.split('\t')[3]),int(event.split('\t')[4])]]
			else:
				if lastgene!='':
					geneinfo[chrom]+=[genepos+[exonpos]]
				lastgene=genename
				genepos=[int(event.split('\t')[3]),int(event.split('\t')[4]),genename,event.split('\t')[6]]
				exonpos=[]
		if exonpos!=[]:
			geneinfo[chrom]+=[genepos+[exonpos]]

	for chrom in geneinfo:
		geneinfo[chrom].sort(key=sortgeneinfo)

	if writeds:
		write_geneinfo()

	return geneinfo


