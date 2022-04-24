import os
import raw_signal
import cluster

def poa_all(outpath,goodchrom):
	try:
		os.mkdir(outpath+'poa_workspace/')
	except:
		pass
	allgf=open(outpath+'confident_genefusion.txt').read().split('\n')[:-1]

	allsupporting=[]
	for c in allgf:
		allsupporting+=c.split('\t')[8].split(',')
	allreadseq={}
	for chrom in goodchrom:
		allinfo=open(outpath+'raw_signal/record_read_'+chrom,'r').read().split('\n')[:-1]
		for c in allinfo:
			if c.split('\t')[3] in allsupporting and c.split('\t')[12]!='':
				allreadseq[c.split('\t')[3]]=c.split('\t')[12]+'\n+\n'+c.split('\t')[13]
	poafile=open(outpath+'poa_workspace/allpoaseq.fa','w')
	for event in allgf:
		event=event.split('\t')
		gfinfo=event[7]+'_'+event[0]+'_'+event[1]+'_'+event[4]+'_'+event[6]
		f=open(outpath+'poa_workspace/suppread_'+gfinfo+'.fastq','w')
		for read in event[8].split(','):
			try:
				if allreadseq[read].split('\n')[1]!='':
					f.write('@'+read+'\n'+allreadseq[read]+'\n')
				else:
					f.write('@'+read+'\n'+allreadseq[read].split('\n')[0]+'\n'+'I'*len(allreadseq[read].split('\n')[0])+'\n')
			except:
				pass
		f.close()
		os.system('bsalign poa '+outpath+'poa_workspace/suppread_'+gfinfo+'.fastq -o '+outpath+'poa_workspace/poa_'+gfinfo+'.fa')
		poactg=open(outpath+'poa_workspace/poa_'+gfinfo+'.fa','r').read().split('>')
		if len(poactg)!=2:
			continue
		poactg=''.join(poactg[1].split('\n')[1:])
		poafile.write('>poa_ctg_'+gfinfo+'\n'+poactg+'\n')
	poafile.close()

	return 0

def polish_bp(outpath,goodchrom,reference,datatype):
	raw_signal.geneinfo=geneinfo
	
	if datatype=='isoseq':	
		os.system('minimap2 -ax splice:hq '+reference+' --secondary=no '+outpath+'poa_workspace/allpoaseq.fa | samtools sort -o '+outpath+'poa_workspace/allpoaseq.bam')
	else:
		os.system('minimap2 -ax splice '+reference+' --secondary=no '+outpath+'poa_workspace/allpoaseq.fa | samtools sort -o '+outpath+'poa_workspace/allpoaseq.bam')

	os.system('samtools index '+outpath+'poa_workspace/allpoaseq.bam')
	os.system('mkdir '+outpath+'poa_workspace/raw_signal/')
	for chrom in goodchrom:
		raw_signal.get_raw_signal(outpath+'poa_workspace/allpoaseq.bam',outpath+'poa_workspace/',chrom,recordseq=False)
	raw_signal.detect_from_split(outpath+'poa_workspace/',goodchrom)
	
	goodpos=open(outpath+'poa_workspace/rawsignal.txt','r').read().split('\n')[:-1]
	allgf=open(outpath+'confident_genefusion.txt').read().split('\n')[:-1]

	gfinfo={}
	for gfeve in allgf:
		event=gfeve.split('\t')
		gfinfo[event[7]+'_'+event[0]+'_'+event[1]+'_'+event[4]+'_'+event[6]]=gfeve

	for goodgf in goodpos:
		goodgf=goodgf.split('\t')
		if goodgf[0]==goodgf[7].split('_')[3] and goodgf[1]==goodgf[7].split('_')[4] :
			oldinfo=gfinfo[goodgf[7][8:]].split('\t')
			newinfo=goodgf[0]+'\t'+goodgf[1]+'\t'+oldinfo[2]+'\t'+goodgf[3]+'\t'+goodgf[4]+'\t'+goodgf[5]+'\t'+goodgf[6]+'\t'+oldinfo[7]+'\t'+oldinfo[8]
			gfinfo[goodgf[7][8:]]=newinfo

	mergedgf=[]
	for gf in gfinfo:
		mergedgf+=[gfinfo[gf]]
	mergedgf=cluster.merge_pair(mergedgf)

	f=open(outpath+'confident_genefusion.txt','w')
	f.write('ID\tGene1\tGene2\tNumSupp\tChrom1\tBreakpoint1\tChrom2\tBreakpoint2\tSupportingReads\n')
	for gf in mergedgf:
		gf=gf.split('\t')
		f.write(gf[7]+'\t'+gf[0]+'\t'+gf[1]+'\t'+gf[2]+'\t'+gf[3]+'\t'+gf[4]+'\t'+gf[5]+'\t'+gf[6]+'\t'+gf[8]+'\n')
	f.close()
	allisoseq=open(outpath+'poa_workspace/allpoaseq.fa','r').read().split('>')[1:]
	poaseq={}
	for c in allisoseq:
		poaseq[c.split('\n')[0].split('_')[2]]=c.split('\n')[1]
	f=open(outpath+'confident_genefusion_transcript_sequence.fa','w')
	for gf in mergedgf:
		gf=gf.split('\t')
		try:
			f.write('>'+gf[7]+'_'+gf[0]+'_'+gf[1]+'_'+gf[3]+'_'+gf[4]+'_'+gf[5]+'_'+gf[6]+'\n'+poaseq[gf[7]]+'\n')
		except:
			pass
	f.close()

	return 0

