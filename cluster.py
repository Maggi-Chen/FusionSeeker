
def findmax(neighbor):
	maxlen=0
	maxkey=''
	for c in neighbor:
		if len(neighbor[c])>maxlen:
			maxkey=c
	return c

def sortcount(a):
	return int(a.split('\t')[2])

def merge_pair_same(samepair):
	if len(samepair)==0:
		return []
	if len(samepair)==1:
		return samepair
	samepair.sort(key=sortcount,reverse=True)
	candi={}
	candi[1]=samepair[0]
	numcandi=1
	for c in samepair[1:]:
		ifclose=0
		for i in range(numcandi):
			d=candi[i+1]
			if abs(int(c.split('\t')[4])-int(d.split('\t')[4]))<=2000 and abs(int(c.split('\t')[6])-int(d.split('\t')[6]))<=2000:
				ifclose=1
				d=d.split('\t')
				c=c.split('\t')
				candi[i+1]=d[0]+'\t'+d[1]+'\t'+str(int(d[2])+int(c[2]))+'\t'+d[3]+'\t'+d[4]+'\t'+d[5]+'\t'+d[6]+'\t'+d[7]+'\t'+d[8]
				break
		if ifclose==0:
			candi[numcandi+1]=c
			numcandi+=1
	merged=[]
	for i in range(numcandi):
		merged+=[candi[i+1]]
	return merged

def merge_pair(mergedgf):
	candi=[]
	samepair=[]
	gene1='';gene2=''
	for c in mergedgf:
		cc=c.split('\t')
		if cc[0]==gene1 and cc[1]==gene2:
			samepair+=[c];continue
		candi+=merge_pair_same(samepair)
		samepair=[c]
		gene1=cc[0];gene2=cc[1]
	candi+=merge_pair_same(samepair)
	return candi

def merge_same(samegf,svid):
	samesvinfo=[]
	for c in samegf:
		samesvinfo+=[svid[c]]
	gfchr=samesvinfo[0]
	bp1=[int(c[4]) for c in samesvinfo]
	bp2=[int(c[6]) for c in samesvinfo]
	bp1=int(int(sum(bp1)/len(samegf)))
	bp2=int(int(sum(bp2)/len(samegf)))
	quality=[]
	allq=[c[8] for c in samesvinfo]
	for c in allq:
		quality+=[int(mm) for mm in c.split(',')]
	quality=int(int(sum(quality)/len(samegf)/2))
	reads=','.join([c[7] for c in samesvinfo])
	numsupp=len(set([c[7] for c in samesvinfo]))
	gfinfo=gfchr[0]+'\t'+gfchr[1]+'\t'+str(numsupp)+'\t'+gfchr[3]+'\t'+str(bp1)+'\t'+gfchr[5]+'\t'+str(bp2)+'\t'+str(quality)+'\t'+reads
	return gfinfo

def cluster_same_dbscan(allinfo,maxdistance,outpath):
	candi=[]
	neighbor={}
	svid={}
	for i in range(len(allinfo)):
		svid['gf'+str(i)]=allinfo[i].split('\t')
		neighbor['gf'+str(i)]=[]
	for i in range(len(allinfo)-1):
		gf1=svid['gf'+str(i)]
		for j in range(i+1,len(allinfo)):
			gf2=svid['gf'+str(j)]
			distance=((int(gf1[4])-int(gf2[4]))**2+(int(gf1[6])-int(gf2[6]))**2)**0.5
			if distance<=maxdistance:
				neighbor['gf'+str(i)]+=['gf'+str(j)]
				neighbor['gf'+str(j)]+=['gf'+str(i)]

	done=[]
	while neighbor!={}:
		startkey=findmax(neighbor)
		samegf=[startkey]
		newadd=neighbor.pop(startkey)
		newadd=[c for c in newadd if c in neighbor]
		while newadd!=[]:
			if newadd[0] not in neighbor:
				newadd.remove(newadd[0])
				continue
			samegf+=[newadd[0]]
			newneighbor=neighbor.pop(newadd[0])
			if len(newneighbor)>=3:
				newneighbor=[c for c in newneighbor if c not in samegf+newadd+done]
				newadd+=newneighbor
			newadd.remove(newadd[0])
		candi+=[merge_same(samegf,svid)]
		done+=samegf
	return candi

def cluster_same(allinfo,maxdistance,outpath):
	candi=''

	if len(allinfo)<5:
		pos1=[int(c.split('\t')[4]) for c in allinfo]
		pos2=[int(c.split('\t')[6]) for c in allinfo]

		pos1=1.0*sum(pos1)/len(pos1)
		pos2=1.0*sum(pos2)/len(pos2)
	else:
		pos1=[int(c.split('\t')[4]) for c in allinfo]
		pos2=[int(c.split('\t')[6]) for c in allinfo]
		pos1.sort()
		pos2.sort()
		pos1=pos1[int((len(pos1)-1)/2)]*1.0
		pos2=pos2[int((len(pos2)-1)/2)]*1.0


	distance=[]
	goodinfo=[]
	for c in allinfo:
		distance+=[((int(c.split('\t')[4])-pos1)**2+(int(c.split('\t')[6])-pos2)**2)**0.5]
		if  ((int(c.split('\t')[4])-pos1)**2+(int(c.split('\t')[6])-pos2)**2)**0.5 <=maxdistance:
			goodinfo+=[c]

	if len(goodinfo)==0:
		return ''
	gfinfo=goodinfo[0].split('\t')
	mappingq=[]
	for c in goodinfo:
		mappingq+=c.split('\t')[8].split(',')
	mappingq=round(1.0*sum([int(c) for c in mappingq])/len(mappingq))
	candi=gfinfo[0]+'\t'+gfinfo[1]+'\t'+str(len(goodinfo))+'\t'+gfinfo[3]+'\t'+str(int(round(1.0*sum([int(c.split('\t')[4]) for c in goodinfo])/len(goodinfo))))+'\t'+gfinfo[5]+'\t'+str(int(round(1.0*sum([int(c.split('\t')[6]) for c in goodinfo])/len(goodinfo))))+'\t'+str(mappingq)+'\t'+','.join([c.split('\t')[7] for c in goodinfo])
	return candi


def cluster_bp(outpath,maxdistance,min_supp):
	allsplit=open(outpath+'rawsignal.txt','r').read().split('\n')[:-1]
	allinfo={}
	mergedgf=[]
	global writeid
	writeid=1

	if min_supp==None:
		min_supp=3+int(len(allsplit)/50000)

	for signal in allsplit:
		if signal.split('\t')[0]+'\t'+signal.split('\t')[1] in allinfo:
			allinfo[signal.split('\t')[0]+'\t'+signal.split('\t')[1]]+=[signal]
			continue
		if signal.split('\t')[1]+'\t'+signal.split('\t')[0]  in allinfo:
			signal=signal.split('\t')
			signal=signal[1]+'\t'+signal[0]+'\t'+signal[2]+'\t'+signal[5]+'\t'+signal[6]+'\t'+signal[3]+'\t'+signal[4]+'\t'+signal[7]+'\t'+signal[8]+'\t'+signal[9]+'\t'+signal[10]+'\t'+signal[11]
			allinfo[signal.split('\t')[0]+'\t'+signal.split('\t')[1]]+=[signal]
			continue
		allinfo[signal.split('\t')[0]+'\t'+signal.split('\t')[1]]=[signal]

	for genepair in allinfo:

		if len(allinfo[genepair])==1:
			gfinfp=allinfo[genepair][0].split('\t')
			mergedgf+=[gfinfp[0]+'\t'+gfinfp[1]+'\t1\t'+gfinfp[3]+'\t'+gfinfp[4]+'\t'+gfinfp[5]+'\t'+gfinfp[6]+'\t'+gfinfp[8]+'\t'+gfinfp[7]]
			continue
		candi=cluster_same_dbscan(allinfo[genepair],maxdistance,outpath)

		writeid+=1
		mergedgf+=candi
	f=open(outpath+'clustered_candidate.txt','w')
	for c in mergedgf:
		f.write(c+'\n')
	f.close()

	mergedgf=[c for c in mergedgf if int(c.split('\t')[2])>=min(3,min_supp)]
	mergedgf=merge_pair(mergedgf)
	mergedgf=[c for c in mergedgf if int(c.split('\t')[2])>=min_supp]
	genepairs=[[c.split('\t')[0],c.split('\t')[1]] for c in mergedgf]
	uniquepairs=[]
	for c in genepairs:
		if c not in uniquepairs:
			uniquepairs+=[c]
	genepairs=uniquepairs
	totalgf=len(genepairs)
	genecount={}
	for c in genepairs:
		if c[0] not in genecount:
			genecount[c[0]]=1
		else:
			genecount[c[0]]+=1
		if c[1] not in genecount:
			genecount[c[1]]=1
		else:
			genecount[c[1]]+=1
	repeating=[c for c in genecount if genecount[c]>=max(6,int(totalgf/20))]
	mergedgf=[c for c in mergedgf if c.split('\t')[0] not in repeating and c.split('\t')[1] not in repeating ]
	f=open(outpath+'confident_genefusion.txt','w')
	idnum=1
	for c in mergedgf:
		c=c.split('\t')
		f.write(c[0]+'\t'+c[1]+'\t'+c[2]+'\t'+c[3]+'\t'+c[4]+'\t'+c[5]+'\t'+c[6]+'\tGF0'+str(idnum)+'\t'+c[8]+'\n')
		idnum+=1
	f.close()


	return 0


	

