try:
	from Bio import SeqIO
	from Bio import Seq
except ImportError:
	SeqIO=None
	Seq=None
import csv
import pickle
import os
import struct
import re
import gzip
import gc
import math
from itertools import zip_longest

AA_MOTIFS=["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
NT_BASES=["A","C","G","T"]
NT_CODONS=["AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG","AGT","ATA","ATC","ATG","ATT","CAA","CAC","CAG","CAT","CCA","CCC","CCG","CCT","CGA","CGC","CGG","CGT","CTA","CTC","CTG","CTT","GAA","GAC","GAG","GAT","GCA","GCC","GCG","GCT","GGA","GGC","GGG","GGT","GTA","GTC","GTG","GTT","TAA","TAC","TAG","TAT","TCA","TCC","TCG","TCT","TGA","TGC","TGG","TGT","TTA","TTC","TTG","TTT"]
SENSE_CODONS=[codon for codon in NT_CODONS if codon not in ["TAA","TAG","TGA"]]


def expand_wildcard_motif(motif,alphabet):
	expanded=[""]
	for character in motif:
		if character=="X":
			expanded=[prefix+base for prefix in expanded for base in alphabet]
		else:
			expanded=[prefix+character for prefix in expanded]
	return expanded


def resolve_posavg_motifs(motif_input,kind,mode):
	kind=int(kind)
	mode=mode.lower()
	motif_input=motif_input.strip()
	compute_scores=False

	# Classic mode preserves legacy behavior exactly.
	if mode=="classic":
		if motif_input=="all":
			compute_scores=True
			if kind==0:
				# Keep legacy codon list behavior.
				motifs=["AAA","AAG","AAC","AAT","AGA","AGG","AGC","AGT","ACA","ACG","ACC","ACT","ATA","ATG","ATC","ATT","GAA","GAG","GAC","GAT","GGA","GGG","GGC","GGT","GCA","GCG","GCC","GCT","GTA","GTG","GTC","GTT","CAA","CAG","CAC","CAT","CGA","CGG","CGC","CGT","CCA","CCG","CCC","CCT","CTA","CTG","CTC","CTT","TAA","TAG","TAC","TAT","TGA","TGG","TGC","TGT","TCA","TCG","TCC","TCT","TTA","TTG","TTC","TTT","TAA","TAG","TGA"]
			elif kind==1:
				motifs=AA_MOTIFS
			else:
				exit()
		else:
			motifs=list(motif_input.split(","))
		return motifs,compute_scores

	# Expanded mode adds wildcard support and score output for expanded sets.
	if motif_input=="all":
		compute_scores=True
		if kind==0:
			return NT_CODONS,compute_scores
		if kind==1:
			return AA_MOTIFS,compute_scores
		exit()

	raw_motifs=[item.strip().upper() for item in motif_input.split(",") if item.strip()!=""]
	motifs=[]
	seen=set()

	for motif in raw_motifs:
		if kind==0 and "X" in motif:
			if len(motif)!=3:
				print("Expanded codon mode requires 3-nt motifs when using X wildcard.")
				exit()
			compute_scores=True
			for expanded in expand_wildcard_motif(motif,NT_BASES):
				if expanded not in seen:
					motifs.append(expanded)
					seen.add(expanded)
		elif kind==1 and "X" in motif:
			compute_scores=True
			for expanded in expand_wildcard_motif(motif,AA_MOTIFS):
				if expanded not in seen:
					motifs.append(expanded)
					seen.add(expanded)
		else:
			if motif not in seen:
				motifs.append(motif)
				seen.add(motif)

	return motifs,compute_scores



# Writegene2 wrapper function
def writegene2_m_wrapper(txtinfile):
	f=open(txtinfile)
	inputparams=[line.strip() for line in f.readlines()] 
	f.close()
	
	variables={}
	nextlineisinput=0
	for inputline in inputparams:
		if inputline=="" or inputline[0]=="#":
			continue
		
		if nextlineisinput!=0:
			variables[nextlineisinput]=inputline
			nextlineisinput=0
		
		if inputline[0:5]=="INPUT":
		# The next line is an input line.
			nextlineisinput=inputline.split("_")[-1]
			
	print("Variables collected:")
	for key in (variables.keys()):
		print(key)
		print(variables[key])
		print("")
	
	
	inputfiles=list(variables["picklenames"].split(","))
	genelist=list(variables["genenames"].split(","))
	
	# Handle Alias conversion:
	aliasdict={}
	f=gzip.open(inputfiles[0],"rb")
	footprints=pickle.load(f)
	f.close()
	for gene in footprints.keys():
		aliasdict[footprints[gene][0]]=gene
	for gnnum in range(len(genelist)):
		if genelist[gnnum][0:3]!="NM_":
			alias=genelist[gnnum]
			if alias in aliasdict:
				genelist[gnnum]=aliasdict[alias]
			else:
				print("Error - alias "+alias+" does not exist. Try the NM_... format name.")
				exit()

	for picklefile in inputfiles:
		writegene2_m(genelist,variables["endmode"],picklefile,variables["outfile"]+"_"+re.split("[./]",picklefile)[-3])


#### writegene2_m function.
### endmode is "all_3" or "all_5"
## This function is used to put the reads that map to a particular gene into a csv file. This can be used for making figures or other analysis.
# genenames is a list of gene names
# endmode is 5' or 3' end alignments
# picklefile is the ribosome density file
# outfile is the name of the output csv file
def writegene2_m(genenames,endmode,picklefile,outfile):
	f=gzip.open(picklefile,"rb")
	footprints=pickle.load(f)
	print("Data loaded.")
	f.close()
	
	writerfile=open(outfile+".csv", "w")	
	writer = csv.writer(writerfile)
		
	for gene in genenames:
		footprintcounts=list(footprints[gene][2][endmode])
		alias=footprints[gene][0]
		writer.writerow([gene]+[alias]+footprintcounts)
	
	writerfile.close()	
	transposecsv(outfile)	
	
def transposecsv(csvfile):
	f=open(csvfile+".csv")
	readercsv=csv.reader(f)
	rows=[row for row in readercsv]
	f.close()

	writerfile=open(csvfile+"_transposed.csv", "w")
	writer = csv.writer(writerfile,delimiter=',')
	for col in zip_longest(*rows, fillvalue=float('nan')):
		writer.writerow(col)
	writerfile.close()
	os.remove(csvfile+".csv")
	os.rename(csvfile+"_transposed.csv",csvfile+".csv")
	
# Extract the given column from a csv file.				
def extractcolumn(rowgen,colnum):
	column=[]
	for row in rowgen:
		if len(row)>colnum:
			column.append(row[colnum])
		else:
			column.append(float('nan'))
	return column
	
	

	

	
# Posavg wrapper function
# This handles the inputs and outputs the metacodon data to a csv file. 
# It allows use of "all" to look at every amino acid or nt combination; in this case pause scores computed from the average data are also written out.
def posavg_m_wrapper(txtinfile):
	f=open(txtinfile)
	inputparams=[line.strip() for line in f.readlines()] 
	f.close()
	
	variables={}
	nextlineisinput=0
	for inputline in inputparams:
		if inputline=="" or inputline[0]=="#":
			continue
		
		if nextlineisinput!=0:
			variables[nextlineisinput]=inputline
			nextlineisinput=0
		
		if inputline[0:5]=="INPUT":
		# The next line is an input line.
			nextlineisinput=inputline.split("_")[-1]
			
	print("Variables collected:")
	for key in (variables.keys()):
		print(key)
		print(variables[key])
		print("")
	
	
	
	inputfiles=list(variables["picklenames"].split(","))
	mode=variables.get("posavgmode","classic")
	motifs,compute_scores=resolve_posavg_motifs(variables["motif"],variables["kind"],mode)
	pause_scores_by_dataset={}
	dataset_labels=[re.split("[./]",picklefile)[-3] for picklefile in inputfiles]
	for dataset_label in dataset_labels:
		pause_scores_by_dataset[dataset_label]={}
		
	
	# Do score and average data.
	writerfile=open(variables["outfile"]+"_avgdata.csv", "w")	# Note no "b" in python 3!
	writer = csv.writer(writerfile)

	batch_size=3
	for batch_start in range(0,len(inputfiles),batch_size):
		batch_files=inputfiles[batch_start:batch_start+batch_size]
		footprints_cache={}
		results_by_pickle={}

		for picklefile in batch_files:
			f=gzip.open(picklefile,"rb")
			footprints_cache[picklefile]=pickle.load(f)
			print("Data loaded for "+picklefile)
			f.close()
			results_by_pickle[picklefile]=posavg_m_multi(
				variables["endmode"],
				picklefile,
				variables["outfile"],
				motifs,
				variables["kind"],
				variables["frame"],
				variables["bkndwindowthresh"],
				variables["bkndwindow"],
				variables["ORFnorm"],
				variables["UTRmode"],
				variables["shift"],
				footprints=footprints_cache[picklefile],
			)

		for motif in motifs:
			for picklefile in batch_files:
				avggene=results_by_pickle[picklefile][motif]
				dataset_label=re.split("[./]",picklefile)[-3]
			
				#Compute pause score
				if compute_scores:	
					bkndwindow=int(variables["bkndwindow"])
					num=sum(avggene[bkndwindow-1:bkndwindow+2])/3		# assume a peak of 3 nt.
					denom=sum(avggene)/len(avggene)	
					if denom==0:
						pause_scores_by_dataset[dataset_label][motif]=float('nan')
					else:
						pause_scores_by_dataset[dataset_label][motif]=num/denom
				
				writer.writerow([dataset_label+"_"+motif]+avggene)

		del results_by_pickle
		del footprints_cache
		gc.collect()
		
	writerfile.close()
	transposecsv(variables["outfile"]+"_avgdata")		
	
	
	# Score data not written out unless doing all. There's not really a reason to do just one score.
	if compute_scores:		
		writerfile=open(variables["outfile"]+"_score.csv", "w")	
		writer = csv.writer(writerfile)
		writer.writerow(["dataset"]+motifs)
		for dataset_label in dataset_labels:
			row=[dataset_label]
			for motif in motifs:
				row.append(pause_scores_by_dataset[dataset_label].get(motif,float('nan')))
			writer.writerow(row)
		writerfile.close()

	
	
	
	

	
	
	
	
#### Codon analysis position average analysis for 64 codons
# endmode is either "all_5" or "all_3".
# picklefile is the ribosome footprints file.
# outfile is the name of the output csv file.
# motif is the nt of aa of the motif that is being averaged. Note you can specify the first or last occurence as noted in the code.
# kind is whether the motif is an amino acid (1) or nucleotide sequence (0)
# frame is frame 0, 1, or 2 to be looking in. Put 3 for all 3 frames.
# bkndwindowthresh this is the minimal rpkm counts within the bkndwindow to be included in the average.
# bkndwindow - This is the size of the half-window around the codon of interest to be included in the average.
# UTRmode is whether doing UTR5,CDS,or UTR3 = 0,1,2, respectively.shift
# ORFnorm will normalize to the ORF for UTRs. 0 does not normalize. A positive value will ORF normalize and specify the rpkm ORF threshold for genes to include. This uses a hard-coded shift of +/-13.
# shift is nt to shift data for site of interest (P site, A site, etc).
# Pause score is computed from the average plot, uses a hardcoded window of +/- 1 nt around the peak.
# There is a hard-coded requirement that for ORFnorm, the UTR region cannot be >10x the ORF. 
def posavg_m(endmode,picklefile,outfile,motif,kind,frame,bkndwindowthresh,bkndwindow,ORFnorm,UTRmode,shift,footprints=None):
	averages=posavg_m_multi(
		endmode,
		picklefile,
		outfile,
		[motif],
		kind,
		frame,
		bkndwindowthresh,
		bkndwindow,
		ORFnorm,
		UTRmode,
		shift,
		footprints=footprints,
	)
	return(averages[motif])


def posavg_m_multi(endmode,picklefile,outfile,motifs,kind,frame,bkndwindowthresh,bkndwindow,ORFnorm,UTRmode,shift,footprints=None):
	if footprints is None:
		f=gzip.open(picklefile,"rb")
		footprints=pickle.load(f)
		print("Data loaded.")
		f.close()

	kind=int(kind)
	frame=int(frame)
	bkndwindowthresh=int(bkndwindowthresh)
	bkndwindow=int(bkndwindow)
	ORFnorm=float(ORFnorm)
	UTRmode=int(UTRmode)
	shift=int(shift)

	if ORFnorm>0 and UTRmode==1:
		print("Invalid entries. Spike checking prevents iORF usage with ORF norm.")
		exit()

	motif_meta=[]
	for motif_raw in motifs:
		firstlast=0
		motif=motif_raw
		if motif[0]=="f":
			motif=motif[1:]
			firstlast=1
		elif motif[0]=="l":
			motif=motif[1:]
			firstlast=2
		if firstlast>0 and kind!=0:
			print("Invalid use of firstlast on an amino acid motif. Nucleotide only.")
			exit()
		motif_meta.append({"raw":motif_raw,"motif":motif,"firstlast":firstlast,"motiflen":len(motif)})

	normal_motif_lookup={}
	firstlast_metas=[]
	for meta in motif_meta:
		if kind==0 and meta["firstlast"]>0:
			firstlast_metas.append(meta)
			continue
		motiflen=meta["motiflen"]
		motif=meta["motif"]
		if motiflen not in normal_motif_lookup:
			normal_motif_lookup[motiflen]={}
		if motif not in normal_motif_lookup[motiflen]:
			normal_motif_lookup[motiflen][motif]=[]
		normal_motif_lookup[motiflen][motif].append(meta)

	averages={meta["raw"]:[0.0 for x in range(2*bkndwindow)] for meta in motif_meta}
	counts_by_motif={meta["raw"]:0 for meta in motif_meta}

	if frame==3:
		framelist=[0,1,2]
	else:
		framelist=[frame]

	for frame in framelist:
		genecount=0
		print("Frame "+str(frame))

		for gene in footprints.keys():
			if gene=="NM_020791.2":
				continue
			if genecount%9000==0:
				print("Genes finished so far="+str(genecount))
			genecount+=1
			ORFstart=int(footprints[gene][3])
			UTR3start=int(footprints[gene][4])

			if shift>=0:
				if UTRmode==0:
					counts=footprints[gene][2][endmode][0:ORFstart-shift]
					genesequence=footprints[gene][1][shift:ORFstart]
					if (ORFstart-shift)<0:
						continue
				elif UTRmode==1:
					counts=footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift]
					genesequence=footprints[gene][1][ORFstart:UTR3start]
					if (ORFstart-shift)<0:
						continue
				elif UTRmode==2:
					if shift==0:
						counts=footprints[gene][2][endmode][UTR3start-shift:]
					else:
						counts=footprints[gene][2][endmode][UTR3start-shift:-shift]
					genesequence=footprints[gene][1][UTR3start:]
					if (UTR3start-shift)<0:
						continue

			if shift<0:
				if UTRmode==0:
					counts=footprints[gene][2][endmode][-shift:ORFstart-shift]
					genesequence=footprints[gene][1][0:ORFstart]
					if (ORFstart-shift)>len(counts):
						continue
				elif UTRmode==1:
					counts=footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift]
					genesequence=footprints[gene][1][ORFstart:UTR3start]
					if (ORFstart-shift)>len(counts) or (UTR3start-shift)>len(counts):
						continue
				elif UTRmode==2:
					counts=footprints[gene][2][endmode][UTR3start-shift:]
					genesequence=footprints[gene][1][UTR3start:shift]
					if (UTR3start-shift)>len(counts):
						continue

			if len(counts)==0:
				continue

			if shift<0:
				shift_loc=-13
			else:
				shift_loc=13
			ORFlevel=(sum(footprints[gene][2][endmode][ORFstart-shift_loc:UTR3start-shift_loc]))/(UTR3start-ORFstart)
			spikecheck=(10*ORFlevel)

			if ORFnorm>0:
				if shift<0:
					if (ORFstart-shift_loc)>len(counts) or (UTR3start-shift_loc)>len(counts):
						continue
				if shift>=0:
					if (ORFstart-shift_loc)<0 or (UTR3start-shift_loc)<0:
						continue
				if ((ORFlevel*1000)<ORFnorm):
					continue

			def get_normfactor_and_check(counts_local,position_nt,genelen_local):
				if not ((position_nt+bkndwindow)<genelen_local and (position_nt-bkndwindow)>0):
					return(None)
				bknd=counts_local[position_nt-bkndwindow:position_nt+bkndwindow]
				if ((sum(bknd))/(len(bknd)/1000))<bkndwindowthresh:
					return(None)
				if ORFnorm==0:
					normfactor=sum(bknd)/(len(bknd))
				else:
					normfactor=ORFlevel
					if (sum(bknd))/(len(bknd))>=spikecheck:
						return(None)
				if normfactor<=0:
					return(None)
				return(normfactor)

			if kind==0:
				genesequence_nt=str(genesequence)
				genelen=len(genesequence_nt)
				for motiflen,motif_lookup in normal_motif_lookup.items():
					i=frame
					while (i+motiflen)<genelen:
						token=genesequence_nt[i:i+motiflen]
						if token in motif_lookup:
							normfactor=get_normfactor_and_check(counts,i,genelen)
							if normfactor is not None:
								for meta in motif_lookup[token]:
									averagegene=averages[meta["raw"]]
									for j in range(len(averagegene)):
										averagegene[j]+=((counts[i-bkndwindow+j])/normfactor)
									counts_by_motif[meta["raw"]]+=1
						i+=motiflen

				for meta in firstlast_metas:
					i=frame
					motif=meta["motif"]
					motiflen=meta["motiflen"]
					firstlast=meta["firstlast"]
					averagegene=averages[meta["raw"]]
					while (i+motiflen)<genelen:
						if motif!=genesequence_nt[i:i+motiflen]:
							i+=motiflen
							continue
						if firstlast==1:
							if motif in genesequence_nt[1:i]:
								i+=motiflen
								continue
						if firstlast==2:
							if motif in genesequence_nt[i+1:]:
								i+=motiflen
								continue
						normfactor=get_normfactor_and_check(counts,i,genelen)
						if normfactor is None:
							i+=motiflen
							continue
						for j in range(len(averagegene)):
							averagegene[j]+=((counts[i-bkndwindow+j])/normfactor)
						counts_by_motif[meta["raw"]]+=1
						i+=motiflen

			elif kind==1:
				counts_frame=counts[frame:]
				genesequence_aa=str(Seq.Seq(genesequence[frame:]).translate())
				genelen_seq=len(genesequence_aa)
				genelen=len(counts_frame)
				for motiflen,motif_lookup in normal_motif_lookup.items():
					i=0
					while (i+motiflen)<genelen_seq:
						token=genesequence_aa[i:i+motiflen]
						if token in motif_lookup:
							position_nt=i*3
							normfactor=get_normfactor_and_check(counts_frame,position_nt,genelen)
							if normfactor is not None:
								for meta in motif_lookup[token]:
									averagegene=averages[meta["raw"]]
									for j in range(len(averagegene)):
										averagegene[j]+=((counts_frame[position_nt-bkndwindow+j])/normfactor)
									counts_by_motif[meta["raw"]]+=1
						i+=motiflen

	for meta in motif_meta:
		motif_raw=meta["raw"]
		count=counts_by_motif[motif_raw]
		if count==0:
			print("Warning: 0 count for motif "+str(motif_raw)+". Leaving average vector at 0.")
			continue
		averagegene=averages[motif_raw]
		for i in range(len(averagegene)):
			averagegene[i]/=count
		print("Number of positions in average")
		print(count)

	return(averages)







# Posstats wrapper function
def posstats_m_wrapper(txtinfile):
	f=open(txtinfile)
	inputparams=[line.strip() for line in f.readlines()] 
	f.close()
	
	variables={}
	nextlineisinput=0
	for inputline in inputparams:
		if inputline=="" or inputline[0]=="#":
			continue
		
		if nextlineisinput!=0:
			variables[nextlineisinput]=inputline
			nextlineisinput=0
		
		if inputline[0:5]=="INPUT":
		# The next line is an input line.
			nextlineisinput=inputline.split("_")[-1]
			
	print("Variables collected:")
	for key in (variables.keys()):
		print(key)
		print(variables[key])
		print("")
	
	inputfiles=list(variables["picklenames"].split(","))
	for picklefile in inputfiles:
		posstats_m(variables["endmode"],picklefile,variables["outfile"]+"_"+re.split("[./]",picklefile)[-3],variables["motif"],variables["kind"],variables["frame"],variables["genethresh"],variables["pkwindow"],variables["shift"],variables["UTRmode"])
	
	

	
### Codon analysis posstats
# shift is a positive value for 5' end aligned, negative for 3' end aligned. It is amount to shift data before computing score (P site, A site, etc).
# UTRmode is whether doing UTR5,CDS,or UTR3 = 0,1,2, respectively.
# endmode is either "all_5" or "all_3".
# picklefile the input footprint file
# outfile is the name of the csv file for output.
# motif - aa or nt sequence to compute scores at
# kind is 0 for nt or 1 for aa.
# frame is 0, 1, or 2 to be looking for motifs in.
# genethresh - this is the minimal counts (rpkm) in the entire gene to be included in the analysis.
# pkwindow - how much on either side of peak to use in numerator of pause score.
# shift is nt to shift data for site of interest (P site, A site, etc).

def posstats_m(endmode,picklefile,outfile,motif,kind,frame,genethresh,pkwindow,shift,UTRmode):
	f=gzip.open(picklefile,"rb")
	footprints=pickle.load(f)
	print("Data loaded.")
	f.close()
	i=0
	motiflen=len(motif)
	
	kind=int(kind)
	frame=int(frame)
	genethresh=int(genethresh)
	pkwindow=int(pkwindow)
	UTRmode=int(UTRmode)
	shift=int(shift)
	
	writerfile=open(outfile+".csv", "w")	# Note no "b" in python 3!
	writer = csv.writer(writerfile)
	
	outputsite=["headers","alias","mRNA_position","localsequence_nt","localsequence_aa","numerator","denominator","pausescore"]
	writer.writerow(outputsite)
	for gene in footprints.keys():
		
		alias=footprints[gene][0]
		ORFstart=int(footprints[gene][3])
		UTR3start=int(footprints[gene][4])
		

		# For positive shifts:			#Note 0 shifts won't work.
		if shift>=0:
			if UTRmode==0:
				counts=footprints[gene][2][endmode][0:ORFstart-shift]
				genesequence=footprints[gene][1][shift:ORFstart]
				if (ORFstart-shift)<0:		# Check for negative indexes: 
					continue
			elif UTRmode==1:
				counts=footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift]
				genesequence=footprints[gene][1][ORFstart:UTR3start]
				if (ORFstart-shift)<0 or (UTR3start-shift)<0:	# Check for negative indexes:
					continue
			elif UTRmode==2:
				counts=footprints[gene][2][endmode][UTR3start-shift:-shift]
				genesequence=footprints[gene][1][UTR3start:]
				if (UTR3start-shift)<0:	# Check for negative indexes:
					continue
		
		# FOR negative SHIFTS (3' end aligned):
		if shift<0:
			if UTRmode==0:
				counts=footprints[gene][2][endmode][-shift:ORFstart-shift]
				genesequence=footprints[gene][1][0:ORFstart]
				if (ORFstart-shift)>len(counts):		# Check for indexes off end:
					continue
			elif UTRmode==1:
				counts=footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift]
				genesequence=footprints[gene][1][ORFstart:UTR3start]
				if (ORFstart-shift)>len(counts) or (UTR3start-shift)>len(counts):	# Check for indexes off end:
					continue
			elif UTRmode==2:
				counts=footprints[gene][2][endmode][UTR3start-shift:]
				genesequence=footprints[gene][1][UTR3start:shift]
				if (UTR3start-shift)>len(counts):	# Check for indexes off end:
					continue
					
					
					
		
		genelen=len(genesequence)
		# Find motifs of interest
		i=frame
		
		if kind==0:
			while (i+motiflen)<genelen: 

				if motif==str(genesequence[i:i+motiflen]):
					if (i+pkwindow)<genelen:
						numerator=sum(counts[i-pkwindow:i+pkwindow+1])/((2*pkwindow+1)/1000)		
						denominator=(sum(counts))/(genelen/1000)

						
						if denominator<genethresh:
							i+=motiflen
							continue
						if denominator>0:
							pausescore=numerator/denominator
						else:
							pausescore=-1

					else:
						i+=motiflen
						continue
					if (i-9)>=0 and (i+motiflen+9)<genelen:
						localsequence=genesequence[i-9:i+motiflen+9]
						localsequenceaa=Seq.Seq(genesequence[i-9:i+motiflen+9]).translate()
					else:
						localsequence="nolocseq"
						localsequenceaa="nolocseq"
						
					if UTRmode==0:
						localposition=i
					elif UTRmode==1:
						localposition=i+ORFstart
					elif UTRmode==2:
						localposition=i+UTR3start
					else:
						exit()
				else:
					i+=motiflen
					continue	
				outputsite=[gene+"_"+str(i),alias,str(localposition),localsequence,localsequenceaa,numerator,denominator,pausescore]
				writer.writerow(outputsite)
				i+=motiflen
		
		elif kind==1:
			i=0
			counts=counts[frame:]
			genesequencent=genesequence[frame:]
			genesequence=Seq.Seq(genesequence[frame:]).translate()
			genelen_seq=len(genesequence)
			genelen=len(counts)
	
			while (i+motiflen)<genelen_seq: 
				if motif==str(genesequence[i:i+motiflen]):
					
					if (i*3+pkwindow)<genelen:
						numerator=sum(counts[i*3-pkwindow:i*3+pkwindow+1])/((2*pkwindow+1)/1000)		
						denominator=(sum(counts))/(genelen/1000)
						if denominator<genethresh:
							i+=motiflen
							continue
						if denominator>0:
							pausescore=numerator/denominator
						else:
							pausescore=-1
					
					else:
						i+=motiflen
						continue
						
					if (i-9)>=0 and (i+motiflen+9)<genelen_seq:
						localsequenceaa=genesequence[i-9:i+motiflen+9]
						localsequence=genesequencent[i*3-18:i*3+motiflen*3+18]
												
						
					else:
						localsequenceaa="nolocseq"
						localsequence="nolocseq"
						
					if UTRmode==0:
						localposition=i
					elif UTRmode==1:
						localposition=i+ORFstart
					elif UTRmode==2:
						localposition=i+UTR3start
					else:
						exit()
				else:
					i+=motiflen
					continue
					
				outputsite=[gene+"_"+str(i),alias,str(localposition),localsequence,localsequenceaa,numerator,denominator,pausescore]
				writer.writerow(outputsite)
				i+=motiflen
	
	
# Metagenewrapper

def metagene_m_wrapper(txtinfile):
	f=open(txtinfile)
	inputparams=[line.strip() for line in f.readlines()] 
	f.close()
	
	variables={}
	nextlineisinput=0
	for inputline in inputparams:
		if inputline=="" or inputline[0]=="#":
			continue
		
		if nextlineisinput!=0:
			variables[nextlineisinput]=inputline
			nextlineisinput=0
		
		if inputline[0:5]=="INPUT":
		# The next line is an input line.
			nextlineisinput=inputline.split("_")[-1]
			
	print("Variables collected:")
	for key in (variables.keys()):
		print(key)
		print(variables[key])
		print("")
	
	outfile=variables["outfile"]
	inputfiles=list(variables["picklenames"].split(","))
	writerfile=open(outfile+".csv", "w")	
	writer = csv.writer(writerfile)
	for picklefile in inputfiles:
		metagene=metagene_m(variables["endmode"],picklefile,variables["kind"],variables["equalweighting"],variables["genethresh"],variables["range5"],variables["range3"])
		writer.writerow([re.split("[./]",picklefile)[-3]]+metagene)
	writerfile.close()	
	transposecsv(outfile)	
	
	
# Metagene
# endmode - must be "all_5". Right now endmode3 is not supported for thresholding since shift 13.
# picklefile - file where footprint density is stored
# kind is whether to do average around start codons (1) or stop codons (2)
# weighting is whether to equally weight all genes (1) or to make an unweighted average based on rpm values (0)
# genethresh is the minimal counts (rpkm) in the entire gene to be included in the analysis.
# range values are the 5' and 3' ranges around the position of interest.

def metagene_m(endmode,picklefile,kind,weighting,genethresh,range5,range3):
	if endmode!="all_5":
		print("ERROR - right now the code doesn't support anything except end5. This is easy to fix though with a little more programming.")
		exit()	
	
	# kind is either 1 or 2 (1 for start codons, 2 for stop codons; other features can come later).
	f=gzip.open(picklefile,"rb")
	footprints=pickle.load(f)
	print("Data loaded.")
	f.close()
	
	kind=int(kind)
	weighting=int(weighting)
	genethresh=int(genethresh)
	range5=int(range5)
	range3=int(range3)
	
	growingmetagene=[0 for x in range(range5+range3)]
	
	genecount=0
	
	for gene in footprints.keys():
		if gene=="NM_020791.2":
			continue	# Skip this gene, bad annotation in Refseq+MANE0.6
		tempmetagene=[0 for x in range(range5+range3)]
		
		ORFstart=int(footprints[gene][3])
		UTR3start=int(footprints[gene][4])
		stopcodon=footprints[gene][1][UTR3start-3:UTR3start]
		
		# For thresholding and equalweighting - assumes no end effects (typically 30 nt) on ends of genes and that shift is 13 roughly.
		if (UTR3start-ORFstart)<=0 or ORFstart<13 or (len(footprints[gene][2][endmode])-UTR3start)<13:
			continue
		else:
			ORFcounts=(sum(footprints[gene][2][endmode][ORFstart+0-13:UTR3start-0-13]))/((UTR3start-ORFstart-0)/1000)
			if ORFcounts<genethresh:
				continue
		
		if kind==1:
			if ORFstart<range5 or (UTR3start-ORFstart)<range3:
				continue
			else:
				tempmetagene=footprints[gene][2][endmode][ORFstart-range5:ORFstart+range3]
				genecount+=1
			
		elif kind==2:
			if (len(footprints[gene][2][endmode])-UTR3start)<range3 or (UTR3start-ORFstart)<range5:
				continue
			else:
				if stopcodon=="TAA":
					tempmetagene=footprints[gene][2][endmode][UTR3start-range5:UTR3start+range3]
					genecount+=1
				else:
					continue

		else:
			print("error")
					
		for position in range(range5+range3):
			if weighting==0:
				growingmetagene[position]+=tempmetagene[position]
			else:
				growingmetagene[position]+=(tempmetagene[position]/(ORFcounts/1000))
			
	for position in range(range5+range3):
		growingmetagene[position]/=genecount
		
	print("Genes included = "+str(genecount))	
	return growingmetagene


def metagene_transcript_m_wrapper(txtinfile):
	f=open(txtinfile)
	inputparams=[line.strip() for line in f.readlines()]
	f.close()

	variables={}
	nextlineisinput=0
	for inputline in inputparams:
		if inputline=="":
			if nextlineisinput!=0:
				variables[nextlineisinput]=""
				nextlineisinput=0
			continue
		if inputline[0]=="#":
			continue
		if inputline[0:5]=="INPUT":
			if nextlineisinput!=0:
				variables[nextlineisinput]=""
			nextlineisinput=inputline[6:]
			continue
		if nextlineisinput!=0:
			variables[nextlineisinput]=inputline
			nextlineisinput=0
	if nextlineisinput!=0:
		variables[nextlineisinput]=""

	print("Variables collected:")
	for key in (variables.keys()):
		print(key)
		print(variables[key])
		print("")

	inputfiles=[x.strip() for x in variables["picklenames"].split(",") if x.strip()!=""]
	outfile=variables["outfile"]

	bins5=int(variables["bins5"])
	binscds=int(variables["binscds"])
	bins3=int(variables["bins3"])
	aggregation_option=int(variables.get("aggregation_option","3"))
	mean_option=int(variables.get("mean_option","1"))
	binning_method=variables.get("binning_method","sample").strip().lower()
	include_gene_ids=_load_gene_list_csv(variables.get("include_genes_csv",""))
	exclude_gene_ids=set([x.strip() for x in variables.get("exclude_genes","").split(",") if x.strip()!=""])
	winsor_percentile=float(variables.get("winsor_percentile","99"))
	trim_fraction=float(variables.get("trim_fraction","0.1"))
	robust_method=variables.get("robust_method","mad").strip().lower()

	row_schema=[]
	for i in range(bins5):
		row_schema.append(("5'UTR",i+1))
	for i in range(binscds):
		row_schema.append(("CDS",i+1))
	for i in range(bins3):
		row_schema.append(("3'UTR",i+1))

	data_by_dataset={}
	dataset_labels=[]
	for picklefile in inputfiles:
		base=os.path.basename(picklefile)
		label=re.sub(r"\.pkl\.gzip$","",base)
		dataset_labels.append(label)
		data_by_dataset[label]=metagene_transcript_m(
			picklefile,
			variables["endmode"],
			variables["equalweighting"],
			variables["genethresh"],
			variables["read_offset5"],
			bins5,
			binscds,
				bins3,
				aggregation_option,
				mean_option,
				include_gene_ids,
				exclude_gene_ids,
				winsor_percentile,
				trim_fraction,
				robust_method,
				binning_method=binning_method,
			)

	writerfile=open(outfile+".csv","w")
	writer=csv.writer(writerfile)
	writer.writerow(["region","region_position","metagene_position"]+dataset_labels)
	start_to_cds_boundary=bins5+1
	cds_to_3utr_boundary=bins5+binscds+1
	print("Boundary 5'UTR->CDS at metagene_position="+str(start_to_cds_boundary))
	print("Boundary CDS->3'UTR at metagene_position="+str(cds_to_3utr_boundary))
	for idx,(site,position) in enumerate(row_schema):
		row=[site,position,idx+1]
		for label in dataset_labels:
			row.append(data_by_dataset[label][idx])
		writer.writerow(row)
	writerfile.close()


def _resample_region(counts,target_bins):
	if target_bins<=0:
		return []
	if len(counts)==0:
		return [0.0 for _ in range(target_bins)]
	if len(counts)==1:
		return [float(counts[0]) for _ in range(target_bins)]
	if target_bins==1:
		return [float(sum(counts))/len(counts)]
	out=[]
	oldmax=len(counts)-1
	newmax=target_bins-1
	for i in range(target_bins):
		pos=(i*oldmax)/float(newmax)
		left=int(pos)
		right=min(left+1,oldmax)
		frac=pos-left
		out.append((1-frac)*float(counts[left])+frac*float(counts[right]))
	return out


def _integrate_region(counts,target_bins):
	if target_bins<=0:
		return []
	if len(counts)==0:
		return [0.0 for _ in range(target_bins)]

	region_len=float(len(counts))
	out=[]
	for i in range(target_bins):
		start=(i*region_len)/float(target_bins)
		end=((i+1)*region_len)/float(target_bins)
		if end<=start:
			out.append(0.0)
			continue

		total=0.0
		left=int(math.floor(start))
		right=int(math.ceil(end))
		for idx in range(left,right):
			if idx<0 or idx>=len(counts):
				continue
			overlap=min(end,idx+1.0)-max(start,float(idx))
			if overlap>0:
				total+=float(counts[idx])*overlap
		out.append(total/(end-start))
	return out


def _shift_counts_5prime(counts,offset):
	shifted=[0.0 for _ in range(len(counts))]
	for idx,val in enumerate(counts):
		newidx=idx+offset
		if newidx>=0 and newidx<len(counts):
			shifted[newidx]+=float(val)
	return shifted


def _load_gene_list_csv(csv_path):
	csv_path=csv_path.strip()
	if csv_path=="":
		return set()

	gene_ids=set()
	with open(csv_path,newline='') as infile:
		reader=csv.reader(infile)
		for row in reader:
			for cell in row:
				token=cell.strip()
				if token!="":
					gene_ids.add(token)
	return gene_ids


def _normalize_gene_token(token):
	return str(token).strip().upper()


def _load_grouped_genes_from_workbook(workbook_path,sheet_name,gene_header,score_header,num_groups,descending):
	scored_rows=[]
	lower_path=workbook_path.lower()

	if lower_path.endswith(".csv"):
		with open(workbook_path,newline='') as infile:
			reader=csv.reader(infile)
			try:
				header_row=next(reader)
			except StopIteration:
				header_row=[]
			header_lookup={}
			for idx,val in enumerate(header_row):
				header_lookup[str(val).strip() if val is not None else ""]=idx

			if gene_header not in header_lookup:
				print("Error: gene column '"+gene_header+"' not found in "+workbook_path)
				exit()
			if score_header not in header_lookup:
				print("Error: score column '"+score_header+"' not found in "+workbook_path)
				exit()

			gene_idx=header_lookup[gene_header]
			score_idx=header_lookup[score_header]
			for row in reader:
				gene=row[gene_idx] if gene_idx<len(row) else None
				score=row[score_idx] if score_idx<len(row) else None
				if gene is None or str(gene).strip()=="":
					continue
				try:
					score_val=float(score)
				except (TypeError,ValueError):
					continue
				if math.isnan(score_val):
					continue
				scored_rows.append((str(gene).strip(),score_val))
	else:
		from openpyxl import load_workbook

		wb=load_workbook(workbook_path,read_only=True,data_only=True)
		ws=wb[sheet_name]
		rows=ws.iter_rows(values_only=True)
		header_row=next(rows)
		header_lookup={}
		for idx,val in enumerate(header_row):
			header_lookup[str(val).strip() if val is not None else ""]=idx

		if gene_header not in header_lookup:
			print("Error: gene column '"+gene_header+"' not found in "+sheet_name)
			exit()
		if score_header not in header_lookup:
			print("Error: score column '"+score_header+"' not found in "+sheet_name)
			exit()

		gene_idx=header_lookup[gene_header]
		score_idx=header_lookup[score_header]
		for row in rows:
			gene=row[gene_idx]
			score=row[score_idx]
			if gene is None or str(gene).strip()=="":
				continue
			try:
				score_val=float(score)
			except (TypeError,ValueError):
				continue
			if math.isnan(score_val):
				continue
			scored_rows.append((str(gene).strip(),score_val))
		wb.close()

	scored_rows.sort(key=lambda x:x[1],reverse=bool(descending))

	groups=[]
	n=len(scored_rows)
	base=n//num_groups
	remainder=n%num_groups
	start=0
	for group_idx in range(num_groups):
		size=base+(1 if group_idx<remainder else 0)
		group_rows=scored_rows[start:start+size]
		start+=size
		if len(group_rows)==0:
			groups.append({"rows":[],"genes":set(),"score_min":None,"score_max":None})
			continue
		group_scores=[row[1] for row in group_rows]
		groups.append({
			"rows":group_rows,
			"genes":set([row[0] for row in group_rows]),
			"score_min":min(group_scores),
			"score_max":max(group_scores),
		})
	return groups


def _resolve_group_genes_for_footprints(footprints,group_genes):
	if len(group_genes)==0:
		return set()
	normalized_group_genes=set([_normalize_gene_token(gene) for gene in group_genes])
	matched=set()
	for gene in footprints.keys():
		alias=footprints[gene][0]
		if _normalize_gene_token(gene) in normalized_group_genes or _normalize_gene_token(alias) in normalized_group_genes:
				matched.add(gene)
				matched.add(alias)
	return matched


def _resolve_group_gene_matches_for_footprints(footprints,group_rows):
	match_records=[]
	include_gene_ids=set()
	group_score_lookup={}
	for gene_name,score_val in group_rows:
		group_score_lookup[_normalize_gene_token(gene_name)]=score_val

	for gene in footprints.keys():
		alias=footprints[gene][0]
		gene_norm=_normalize_gene_token(gene)
		alias_norm=_normalize_gene_token(alias)
		matched_workbook_gene=None
		if gene_norm in group_score_lookup:
			matched_workbook_gene=gene
		elif alias_norm in group_score_lookup:
			matched_workbook_gene=alias
		if matched_workbook_gene is None:
			continue
		include_gene_ids.add(gene)
		include_gene_ids.add(alias)
		match_records.append({
			"gene_id":gene,
			"alias":alias,
			"workbook_gene":matched_workbook_gene,
			"group_score":group_score_lookup[_normalize_gene_token(matched_workbook_gene)],
		})

	return include_gene_ids,match_records


def _percentile_from_sorted(sorted_vals,pct):
	if len(sorted_vals)==0:
		return 0.0
	if pct<=0:
		return float(sorted_vals[0])
	if pct>=100:
		return float(sorted_vals[-1])
	pos=(pct/100.0)*(len(sorted_vals)-1)
	left=int(pos)
	right=min(left+1,len(sorted_vals)-1)
	frac=pos-left
	return (1-frac)*float(sorted_vals[left])+frac*float(sorted_vals[right])


def _aggregate_gene_vectors(gene_vectors,total_len,aggregation_option,mean_option,winsor_percentile,trim_fraction):
	# mean_option is reserved for future mean variants. Option 1 is arithmetic mean.
	if mean_option!=1:
		print("Warning: mean_option currently supports only option 1. Using option 1.")

	out=[0.0 for _ in range(total_len)]
	n=len(gene_vectors)
	if n==0:
		return out

	for pos in range(total_len):
		col=[gene_vectors[g][pos] for g in range(n)]
		col.sort()

		if aggregation_option in [0,1,6]:
			out[pos]=sum(col)/len(col)
		elif aggregation_option==2:
			upper=_percentile_from_sorted(col,winsor_percentile)
			lower=_percentile_from_sorted(col,max(0.0,100.0-winsor_percentile))
			capped=[]
			for val in col:
				if val>upper:
					capped.append(upper)
				elif val<lower:
					capped.append(lower)
				else:
					capped.append(val)
			out[pos]=sum(capped)/len(capped)
		elif aggregation_option==3:
			mid=len(col)//2
			if len(col)%2==1:
				out[pos]=col[mid]
			else:
				out[pos]=(col[mid-1]+col[mid])/2.0
		elif aggregation_option==4:
			k=int(len(col)*trim_fraction)
			if (2*k)>=len(col):
				out[pos]=sum(col)/len(col)
			else:
				trimmed=col[k:len(col)-k]
				out[pos]=sum(trimmed)/len(trimmed)
		elif aggregation_option==5:
			logvals=[math.log1p(max(0.0,v)) for v in col]
			out[pos]=math.expm1(sum(logvals)/len(logvals))
		else:
			print("Warning: invalid aggregation_option "+str(aggregation_option)+". Using option 3 (median).")
			mid=len(col)//2
			if len(col)%2==1:
				out[pos]=col[mid]
			else:
				out[pos]=(col[mid-1]+col[mid])/2.0

	return out


def _robust_scale_vector(vec,robust_method):
	if len(vec)==0:
		return vec
	sorted_vals=sorted(vec)
	mid=len(sorted_vals)//2
	if len(sorted_vals)%2==1:
		median=sorted_vals[mid]
	else:
		median=(sorted_vals[mid-1]+sorted_vals[mid])/2.0

	if robust_method=="zscore":
		mean=sum(vec)/len(vec)
		var=sum([(x-mean)*(x-mean) for x in vec])/len(vec)
		sd=math.sqrt(var)
		if sd==0:
			return [0.0 for _ in vec]
		return [(x-mean)/sd for x in vec]

	madvals=[abs(x-median) for x in vec]
	madvals_sorted=sorted(madvals)
	mid2=len(madvals_sorted)//2
	if len(madvals_sorted)%2==1:
		mad=madvals_sorted[mid2]
	else:
		mad=(madvals_sorted[mid2-1]+madvals_sorted[mid2])/2.0
	scale=1.4826*mad
	if scale==0:
		return [0.0 for _ in vec]
	return [(x-median)/scale for x in vec]


def metagene_transcript_m(picklefile,endmode,equalweighting,genethresh,read_offset5,bins5,binscds,bins3,aggregation_option=3,mean_option=1,include_gene_ids=None,exclude_gene_ids=None,winsor_percentile=99.0,trim_fraction=0.1,robust_method="mad",footprints=None,return_gene_count=False,return_gene_vectors=False,binning_method="sample"):
	if endmode!="all_5":
		print("ERROR - transcript metagene currently expects endmode=all_5 so read_offset5 is well-defined.")
		exit()
	
	if footprints is None:
		f=gzip.open(picklefile,"rb")
		footprints=pickle.load(f)
		print("Data loaded for "+picklefile)
		f.close()

	equalweighting=int(equalweighting)
	genethresh=float(genethresh)
	read_offset5=int(read_offset5)
	aggregation_option=int(aggregation_option)
	mean_option=int(mean_option)
	winsor_percentile=float(winsor_percentile)
	trim_fraction=float(trim_fraction)
	binning_method=binning_method.strip().lower()
	if include_gene_ids is None:
		include_gene_ids=set()
	if exclude_gene_ids is None:
		exclude_gene_ids=set()

	total_len=bins5+binscds+bins3
	gene_vectors=[]
	gene_vector_records=[]
	genecount=0

	for gene in footprints.keys():
		if gene=="NM_020791.2":
			continue
		alias=footprints[gene][0]
		if len(include_gene_ids)>0 and gene not in include_gene_ids and alias not in include_gene_ids:
			continue
		if gene in exclude_gene_ids or alias in exclude_gene_ids:
			continue
		if aggregation_option==1 and gene in exclude_gene_ids:
			continue

		ORFstart=int(footprints[gene][3])
		UTR3start=int(footprints[gene][4])
		counts_raw=footprints[gene][2][endmode]

		if ORFstart<=0 or UTR3start<=ORFstart or UTR3start>len(counts_raw):
			continue

		counts=_shift_counts_5prime(counts_raw,read_offset5)
		five_region=counts[:ORFstart]
		cds_region=counts[ORFstart:UTR3start]
		three_region=counts[UTR3start:]

		if len(cds_region)==0:
			continue
		orf_rpkm=(sum(cds_region))/((len(cds_region))/1000.0)
		if orf_rpkm<genethresh:
			continue

		normfactor=1.0
		if equalweighting==1:
			normfactor=orf_rpkm/1000.0
			if normfactor<=0:
				continue

			if binning_method=="integrate":
				resample_fn=_integrate_region
			else:
				resample_fn=_resample_region

			vec=[]
			vec.extend(resample_fn(five_region,bins5))
			vec.extend(resample_fn(cds_region,binscds))
			vec.extend(resample_fn(three_region,bins3))

			vec=[x/normfactor for x in vec]
			if aggregation_option==6:
				vec=_robust_scale_vector(vec,robust_method)
			gene_vectors.append(vec)
			if return_gene_vectors:
				gene_vector_records.append((gene,alias,vec))
			genecount+=1

	if genecount==0:
		print("No genes passed filters for "+picklefile)
		if return_gene_count and return_gene_vectors:
			return [0.0 for _ in range(total_len)],0,[]
		if return_gene_count:
			return [0.0 for _ in range(total_len)],0
		return [0.0 for _ in range(total_len)]

	grow=_aggregate_gene_vectors(gene_vectors,total_len,aggregation_option,mean_option,winsor_percentile,trim_fraction)

	print("Genes included = "+str(genecount))
	if return_gene_count and return_gene_vectors:
		return grow,genecount,gene_vector_records
	if return_gene_count:
		return grow,genecount
	return grow


def metagene_grouped_m_wrapper(txtinfile):
	f=open(txtinfile)
	inputparams=[line.strip() for line in f.readlines()]
	f.close()

	variables={}
	nextlineisinput=0
	for inputline in inputparams:
		if inputline=="":
			if nextlineisinput!=0:
				variables[nextlineisinput]=""
				nextlineisinput=0
			continue
		if inputline[0]=="#":
			continue
		if inputline[0:5]=="INPUT":
			if nextlineisinput!=0:
				variables[nextlineisinput]=""
			nextlineisinput=inputline[6:]
			continue
		if nextlineisinput!=0:
			variables[nextlineisinput]=inputline
			nextlineisinput=0
	if nextlineisinput!=0:
		variables[nextlineisinput]=""

	print("Variables collected:")
	for key in (variables.keys()):
		print(key)
		print(variables[key])
		print("")

	inputfiles=[x.strip() for x in variables["picklenames"].split(",") if x.strip()!=""]
	outfile=variables["outfile"]
	bins5=int(variables["bins5"])
	binscds=int(variables["binscds"])
	bins3=int(variables["bins3"])
	aggregation_option=int(variables.get("aggregation_option","3"))
	mean_option=int(variables.get("mean_option","1"))
	exclude_gene_ids=set([x.strip() for x in variables.get("exclude_genes","").split(",") if x.strip()!=""])
	winsor_percentile=float(variables.get("winsor_percentile","99"))
	trim_fraction=float(variables.get("trim_fraction","0.1"))
	robust_method=variables.get("robust_method","mad").strip().lower()
	binning_method=variables.get("binning_method","sample").strip().lower()
	num_groups=int(variables.get("num_groups","4"))
	descending=int(variables.get("sort_descending","0"))
	group_prefix=variables.get("group_prefix","Group").strip()
	if group_prefix=="":
		group_prefix="Group"

	group_defs=_load_grouped_genes_from_workbook(
		variables["group_workbook"],
		variables.get("group_sheet","All_gene_pass_rates"),
		variables.get("group_gene_column","Gene"),
		variables.get("group_score_column","average difference"),
		num_groups,
		descending,
	)

	row_schema=[]
	for i in range(bins5):
		row_schema.append(("5'UTR",i+1))
	for i in range(binscds):
		row_schema.append(("CDS",i+1))
	for i in range(bins3):
		row_schema.append(("3'UTR",i+1))

	writerfile=open(outfile+".csv","w")
	writer=csv.writer(writerfile)
	writer.writerow(["group","group_rank","workbook_genes","max_genes_included_across_datasets","region","region_position","metagene_position"]+[re.sub(r"\.pkl\.gzip$","",os.path.basename(picklefile)) for picklefile in inputfiles])

	summaryfile=open(outfile+"_group_summary.csv","w")
	summarywriter=csv.writer(summaryfile)
	summarywriter.writerow(["group","group_rank","dataset","workbook_genes","matched_pickle_genes","genes_included","score_min","score_max"])

	memberfile=open(outfile+"_group_membership.csv","w")
	memberwriter=csv.writer(memberfile)
	memberwriter.writerow(["group","group_rank","workbook_gene","group_score"])

	datasetmemberfile=open(outfile+"_dataset_gene_membership.csv","w")
	datasetmemberwriter=csv.writer(datasetmemberfile)
	datasetmemberwriter.writerow(["group","group_rank","dataset","workbook_gene","group_score","matched_gene_id","matched_alias"])

	vectorfile=open(outfile+"_gene_vectors.csv","w")
	vectorwriter=csv.writer(vectorfile)
	vector_header=["group","group_rank","dataset","workbook_gene","group_score","gene_id","alias"]
	for idx,(site,position) in enumerate(row_schema):
		vector_header.append(site+"_"+str(position)+"_pos"+str(idx+1))
	vectorwriter.writerow(vector_header)

	for group_idx,group_def in enumerate(group_defs,1):
		group_name=group_prefix+str(group_idx)
		for gene_name,score_val in group_def["rows"]:
			memberwriter.writerow([group_name,group_idx,gene_name,score_val])

		group_vectors={}
		group_included_counts={}
		matched_counts={}

		for picklefile in inputfiles:
			f=gzip.open(picklefile,"rb")
			footprints=pickle.load(f)
			print("Data loaded for "+picklefile)
			f.close()

			include_gene_ids,match_records=_resolve_group_gene_matches_for_footprints(footprints,group_def["rows"])
			matched_counts[picklefile]=len(include_gene_ids)
			match_lookup={}
			for record in match_records:
				datasetmemberwriter.writerow([
					group_name,
					group_idx,
					re.sub(r"\.pkl\.gzip$","",os.path.basename(picklefile)),
					record["workbook_gene"],
					record["group_score"],
					record["gene_id"],
					record["alias"],
				])
				match_lookup[record["gene_id"]]=record

			vector,included_count,gene_vector_records=metagene_transcript_m(
				picklefile,
				variables["endmode"],
				variables["equalweighting"],
				variables["genethresh"],
				variables["read_offset5"],
				bins5,
				binscds,
				bins3,
				aggregation_option,
				mean_option,
				include_gene_ids,
				exclude_gene_ids,
				winsor_percentile,
					trim_fraction,
					robust_method,
					footprints=footprints,
					return_gene_count=True,
					return_gene_vectors=True,
					binning_method=binning_method,
				)
			group_vectors[picklefile]=vector
			group_included_counts[picklefile]=included_count
			for gene_id,alias,gene_vec in gene_vector_records:
				record=match_lookup.get(gene_id,{"workbook_gene":"","group_score":""})
				vectorwriter.writerow([
					group_name,
					group_idx,
					re.sub(r"\.pkl\.gzip$","",os.path.basename(picklefile)),
					record["workbook_gene"],
					record["group_score"],
					gene_id,
					alias,
				]+gene_vec)
			summarywriter.writerow([
				group_name,
				group_idx,
				re.sub(r"\.pkl\.gzip$","",os.path.basename(picklefile)),
				len(group_def["genes"]),
				matched_counts[picklefile],
				included_count,
				group_def["score_min"],
				group_def["score_max"],
			])

		genes_included=max(group_included_counts.values()) if len(group_included_counts)>0 else 0
		for idx,(site,position) in enumerate(row_schema):
			row=[group_name,group_idx,len(group_def["genes"]),genes_included,site,position,idx+1]
			for picklefile in inputfiles:
				row.append(group_vectors[picklefile][idx])
			writer.writerow(row)

	writerfile.close()
	summaryfile.close()
	memberfile.close()
	datasetmemberfile.close()
	vectorfile.close()


def _compute_codon_usage_rows(footprints):
	rows=[]
	for gene in footprints.keys():
		if gene=="NM_020791.2":
			continue
		alias=footprints[gene][0]
		ORFstart=int(footprints[gene][3])
		UTR3start=int(footprints[gene][4])
		sequence=str(footprints[gene][1]).upper()
		if ORFstart<0 or UTR3start<=ORFstart or UTR3start>len(sequence):
			continue

		cds=sequence[ORFstart:UTR3start]
		codon_counts={}
		for codon in SENSE_CODONS:
			codon_counts[codon]=0

		for i in range(0,len(cds)-2,3):
			codon=cds[i:i+3]
			if codon in codon_counts:
				codon_counts[codon]+=1

		total_sense_codons=sum(codon_counts.values())
		if total_sense_codons<=0:
			continue

		codon_percents={}
		for codon in SENSE_CODONS:
			codon_percents[codon]=(100.0*codon_counts[codon])/total_sense_codons

		rows.append({
			"gene_id":gene,
			"alias":alias,
			"total_sense_codons":total_sense_codons,
			"codon_counts":codon_counts,
			"codon_percents":codon_percents,
		})
	return rows


def _compute_nt_usage_rows(footprints):
	rows=[]
	for gene in footprints.keys():
		if gene=="NM_020791.2":
			continue
		alias=footprints[gene][0]
		ORFstart=int(footprints[gene][3])
		UTR3start=int(footprints[gene][4])
		sequence=str(footprints[gene][1]).upper()
		if ORFstart<0 or UTR3start<=ORFstart or UTR3start>len(sequence):
			continue

		cds=sequence[ORFstart:UTR3start]
		base_counts={"A":0,"C":0,"G":0,"T":0}
		for base in cds:
			if base in base_counts:
				base_counts[base]+=1

		total_bases=sum(base_counts.values())
		if total_bases<=0:
			continue

		base_percents={}
		for base in ["A","C","G","T"]:
			base_percents[base]=(100.0*base_counts[base])/total_bases

		rows.append({
			"gene_id":gene,
			"alias":alias,
			"total_bases":total_bases,
			"base_counts":base_counts,
			"base_percents":base_percents,
		})
	return rows


def _top_bottom_quartile_rows(scored_rows):
	n=len(scored_rows)
	if n==0:
		return [],[]
	base=n//4
	remainder=n%4
	group_sizes=[]
	for group_idx in range(4):
		group_sizes.append(base+(1 if group_idx<remainder else 0))
	bottom_size=group_sizes[0]
	top_size=group_sizes[3]
	bottom_rows=scored_rows[:bottom_size]
	top_rows=scored_rows[n-top_size:]
	return bottom_rows,top_rows


def metagene_grouped_codonusage_m_wrapper(txtinfile):
	f=open(txtinfile)
	inputparams=[line.strip() for line in f.readlines()]
	f.close()

	variables={}
	nextlineisinput=0
	for inputline in inputparams:
		if inputline=="":
			if nextlineisinput!=0:
				variables[nextlineisinput]=""
				nextlineisinput=0
			continue
		if inputline[0]=="#":
			continue
		if inputline[0:5]=="INPUT":
			if nextlineisinput!=0:
				variables[nextlineisinput]=""
			nextlineisinput=inputline[6:]
			continue
		if nextlineisinput!=0:
			variables[nextlineisinput]=inputline
			nextlineisinput=0
	if nextlineisinput!=0:
		variables[nextlineisinput]=""

	print("Variables collected:")
	for key in (variables.keys()):
		print(key)
		print(variables[key])
		print("")

	inputfiles=[x.strip() for x in variables["picklenames"].split(",") if x.strip()!=""]
	outfile=variables["outfile"]
	bins5=int(variables["bins5"])
	binscds=int(variables["binscds"])
	bins3=int(variables["bins3"])
	aggregation_option=int(variables.get("aggregation_option","3"))
	mean_option=int(variables.get("mean_option","1"))
	exclude_gene_ids=set([x.strip() for x in variables.get("exclude_genes","").split(",") if x.strip()!=""])
	winsor_percentile=float(variables.get("winsor_percentile","99"))
	trim_fraction=float(variables.get("trim_fraction","0.1"))
	robust_method=variables.get("robust_method","mad").strip().lower()
	binning_method=variables.get("binning_method","sample").strip().lower()

	row_schema=[]
	for i in range(bins5):
		row_schema.append(("5'UTR",i+1))
	for i in range(binscds):
		row_schema.append(("CDS",i+1))
	for i in range(bins3):
		row_schema.append(("3'UTR",i+1))

	writerfile=open(outfile+".csv","w")
	writer=csv.writer(writerfile)
	writer.writerow(["dataset","codon","usage_group","quartile_gene_count","genes_included","region","region_position","metagene_position","value"])

	summaryfile=open(outfile+"_group_summary.csv","w")
	summarywriter=csv.writer(summaryfile)
	summarywriter.writerow(["dataset","codon","usage_group","quartile_gene_count","genes_included","usage_min_percent","usage_max_percent"])

	codonusagefile=open(outfile+"_codon_usage.csv","w")
	codonusagewriter=csv.writer(codonusagefile)
	codonusagewriter.writerow(["dataset","gene_id","alias","total_sense_codons"]+SENSE_CODONS)

	memberfile=open(outfile+"_dataset_gene_membership.csv","w")
	memberwriter=csv.writer(memberfile)
	memberwriter.writerow(["dataset","codon","usage_group","gene_id","alias","codon_percent","codon_count","total_sense_codons"])

	for picklefile in inputfiles:
		dataset_label=re.sub(r"\.pkl\.gzip$","",os.path.basename(picklefile))
		f=gzip.open(picklefile,"rb")
		footprints=pickle.load(f)
		print("Data loaded for "+picklefile)
		f.close()

		codon_usage_rows=_compute_codon_usage_rows(footprints)
		row_lookup={}
		for row in codon_usage_rows:
			row_lookup[row["gene_id"]]=row
			codonusagewriter.writerow([dataset_label,row["gene_id"],row["alias"],row["total_sense_codons"]]+[row["codon_percents"][codon] for codon in SENSE_CODONS])

		for codon in SENSE_CODONS:
			scored_rows=[]
			for row in codon_usage_rows:
				scored_rows.append({
					"gene_id":row["gene_id"],
					"alias":row["alias"],
					"percent":row["codon_percents"][codon],
					"count":row["codon_counts"][codon],
					"total_sense_codons":row["total_sense_codons"],
				})
			scored_rows.sort(key=lambda x:(x["percent"],x["gene_id"]))
			bottom_rows,top_rows=_top_bottom_quartile_rows(scored_rows)

			for usage_group,group_rows in [("bottom_quartile",bottom_rows),("top_quartile",top_rows)]:
				include_gene_ids=set([row["gene_id"] for row in group_rows])
				for row in group_rows:
					memberwriter.writerow([dataset_label,codon,usage_group,row["gene_id"],row["alias"],row["percent"],row["count"],row["total_sense_codons"]])

				vector,included_count=metagene_transcript_m(
					picklefile,
					variables["endmode"],
					variables["equalweighting"],
					variables["genethresh"],
					variables["read_offset5"],
					bins5,
					binscds,
					bins3,
					aggregation_option,
					mean_option,
					include_gene_ids,
					exclude_gene_ids,
					winsor_percentile,
					trim_fraction,
					robust_method,
					footprints=footprints,
					return_gene_count=True,
					binning_method=binning_method,
				)

				if len(group_rows)>0:
					usage_min=min([row["percent"] for row in group_rows])
					usage_max=max([row["percent"] for row in group_rows])
				else:
					usage_min=""
					usage_max=""

				summarywriter.writerow([dataset_label,codon,usage_group,len(group_rows),included_count,usage_min,usage_max])
				for idx,(site,position) in enumerate(row_schema):
					writer.writerow([dataset_label,codon,usage_group,len(group_rows),included_count,site,position,idx+1,vector[idx]])

	writerfile.close()
	summaryfile.close()
	codonusagefile.close()
	memberfile.close()


def metagene_grouped_ntusage_m_wrapper(txtinfile):
	f=open(txtinfile)
	inputparams=[line.strip() for line in f.readlines()]
	f.close()

	variables={}
	nextlineisinput=0
	for inputline in inputparams:
		if inputline=="":
			if nextlineisinput!=0:
				variables[nextlineisinput]=""
				nextlineisinput=0
			continue
		if inputline[0]=="#":
			continue
		if inputline[0:5]=="INPUT":
			if nextlineisinput!=0:
				variables[nextlineisinput]=""
			nextlineisinput=inputline[6:]
			continue
		if nextlineisinput!=0:
			variables[nextlineisinput]=inputline
			nextlineisinput=0
	if nextlineisinput!=0:
		variables[nextlineisinput]=""

	print("Variables collected:")
	for key in (variables.keys()):
		print(key)
		print(variables[key])
		print("")

	inputfiles=[x.strip() for x in variables["picklenames"].split(",") if x.strip()!=""]
	outfile=variables["outfile"]
	bins5=int(variables["bins5"])
	binscds=int(variables["binscds"])
	bins3=int(variables["bins3"])
	aggregation_option=int(variables.get("aggregation_option","3"))
	mean_option=int(variables.get("mean_option","1"))
	exclude_gene_ids=set([x.strip() for x in variables.get("exclude_genes","").split(",") if x.strip()!=""])
	winsor_percentile=float(variables.get("winsor_percentile","99"))
	trim_fraction=float(variables.get("trim_fraction","0.1"))
	robust_method=variables.get("robust_method","mad").strip().lower()
	binning_method=variables.get("binning_method","sample").strip().lower()

	row_schema=[]
	for i in range(bins5):
		row_schema.append(("5'UTR",i+1))
	for i in range(binscds):
		row_schema.append(("CDS",i+1))
	for i in range(bins3):
		row_schema.append(("3'UTR",i+1))

	writerfile=open(outfile+".csv","w")
	writer=csv.writer(writerfile)
	writer.writerow(["dataset","base","usage_group","quartile_gene_count","genes_included","region","region_position","metagene_position","value"])

	summaryfile=open(outfile+"_group_summary.csv","w")
	summarywriter=csv.writer(summaryfile)
	summarywriter.writerow(["dataset","base","usage_group","quartile_gene_count","genes_included","usage_min_percent","usage_max_percent"])

	usagefile=open(outfile+"_nt_usage.csv","w")
	usagewriter=csv.writer(usagefile)
	usagewriter.writerow(["dataset","gene_id","alias","total_bases","A","C","G","T"])

	memberfile=open(outfile+"_dataset_gene_membership.csv","w")
	memberwriter=csv.writer(memberfile)
	memberwriter.writerow(["dataset","base","usage_group","gene_id","alias","base_percent","base_count","total_bases"])

	for picklefile in inputfiles:
		dataset_label=re.sub(r"\.pkl\.gzip$","",os.path.basename(picklefile))
		f=gzip.open(picklefile,"rb")
		footprints=pickle.load(f)
		print("Data loaded for "+picklefile)
		f.close()

		nt_usage_rows=_compute_nt_usage_rows(footprints)
		for row in nt_usage_rows:
			usagewriter.writerow([dataset_label,row["gene_id"],row["alias"],row["total_bases"],row["base_percents"]["A"],row["base_percents"]["C"],row["base_percents"]["G"],row["base_percents"]["T"]])

		for base in ["A","C","G","T"]:
			scored_rows=[]
			for row in nt_usage_rows:
				scored_rows.append({
					"gene_id":row["gene_id"],
					"alias":row["alias"],
					"percent":row["base_percents"][base],
					"count":row["base_counts"][base],
					"total_bases":row["total_bases"],
				})
			scored_rows.sort(key=lambda x:(x["percent"],x["gene_id"]))
			bottom_rows,top_rows=_top_bottom_quartile_rows(scored_rows)

			for usage_group,group_rows in [("bottom_quartile",bottom_rows),("top_quartile",top_rows)]:
				include_gene_ids=set([row["gene_id"] for row in group_rows])
				for row in group_rows:
					memberwriter.writerow([dataset_label,base,usage_group,row["gene_id"],row["alias"],row["percent"],row["count"],row["total_bases"]])

				vector,included_count=metagene_transcript_m(
					picklefile,
					variables["endmode"],
					variables["equalweighting"],
					variables["genethresh"],
					variables["read_offset5"],
					bins5,
					binscds,
					bins3,
					aggregation_option,
					mean_option,
					include_gene_ids,
					exclude_gene_ids,
					winsor_percentile,
					trim_fraction,
					robust_method,
					footprints=footprints,
					return_gene_count=True,
					binning_method=binning_method,
				)

				if len(group_rows)>0:
					usage_min=min([row["percent"] for row in group_rows])
					usage_max=max([row["percent"] for row in group_rows])
				else:
					usage_min=""
					usage_max=""

				summarywriter.writerow([dataset_label,base,usage_group,len(group_rows),included_count,usage_min,usage_max])
				for idx,(site,position) in enumerate(row_schema):
					writer.writerow([dataset_label,base,usage_group,len(group_rows),included_count,site,position,idx+1,vector[idx]])

	writerfile.close()
	summaryfile.close()
	usagefile.close()
	memberfile.close()


def metagene_codonusage_slope_summary_m_wrapper(txtinfile):
	f=open(txtinfile)
	inputparams=[line.strip() for line in f.readlines()]
	f.close()

	variables={}
	nextlineisinput=0
	for inputline in inputparams:
		if inputline=="":
			if nextlineisinput!=0:
				variables[nextlineisinput]=""
				nextlineisinput=0
			continue
		if inputline[0]=="#":
			continue
		if inputline[0:5]=="INPUT":
			if nextlineisinput!=0:
				variables[nextlineisinput]=""
			nextlineisinput=inputline[6:]
			continue
		if nextlineisinput!=0:
			variables[nextlineisinput]=inputline
			nextlineisinput=0
	if nextlineisinput!=0:
		variables[nextlineisinput]=""

	print("Variables collected:")
	for key in (variables.keys()):
		print(key)
		print(variables[key])
		print("")

	infile=variables["metagene_csv"]
	outfile=variables["outfile"]
	start_pos=int(variables.get("start_position","106"))
	end_pos=int(variables.get("end_position","395"))

	group_values={}
	feature_field=""
	with open(infile,newline='') as csvfile:
		reader=csv.DictReader(csvfile)
		fieldnames=reader.fieldnames if reader.fieldnames is not None else []
		if "codon" in fieldnames:
			feature_field="codon"
		elif "base" in fieldnames:
			feature_field="base"
		else:
			print("Error: metagene csv must contain either a 'codon' or 'base' column.")
			exit()
		for row in reader:
			dataset=row["dataset"]
			feature_value=row[feature_field]
			usage_group=row["usage_group"]
			metagene_position=int(row["metagene_position"])
			value=float(row["value"])
			key=(dataset,feature_value,usage_group)
			if key not in group_values:
				group_values[key]={}
			group_values[key][metagene_position]=value

	writerfile=open(outfile+".csv","w")
	writer=csv.writer(writerfile)
	writer.writerow(["dataset",feature_field,"first_quartile_slope","fourth_quartile_slope","difference_fourth_minus_first"])

	sorted_keys=sorted(set([(dataset,feature_value) for dataset,feature_value,_usage_group in group_values.keys()]))
	for dataset,feature_value in sorted_keys:
		first_vals=group_values.get((dataset,feature_value,"bottom_quartile"),{})
		fourth_vals=group_values.get((dataset,feature_value,"top_quartile"),{})

		first_slope=_linear_regression_slope_from_positions(first_vals,start_pos,end_pos)
		fourth_slope=_linear_regression_slope_from_positions(fourth_vals,start_pos,end_pos)

		if first_slope=="" or fourth_slope=="":
			diff=""
		else:
			diff=fourth_slope-first_slope

		writer.writerow([dataset,feature_value,first_slope,fourth_slope,diff])

	writerfile.close()


def _linear_regression_slope_from_positions(value_lookup,start_pos,end_pos):
	xvals=[]
	yvals=[]
	for pos in range(start_pos,end_pos+1):
		if pos in value_lookup:
			xvals.append(float(pos))
			yvals.append(float(value_lookup[pos]))

	if len(xvals)<2:
		return ""

	xmean=sum(xvals)/len(xvals)
	ymean=sum(yvals)/len(yvals)
	num=0.0
	den=0.0
	for xval,yval in zip(xvals,yvals):
		num+=((xval-xmean)*(yval-ymean))
		den+=((xval-xmean)*(xval-xmean))
	if den==0:
		return ""
	return num/den



	
	# Wrapper function for counting reads on genes.
# Handles file I/O and calling the counter.
def genelist_m_wrapper(txtinfile):

	f=open(txtinfile)
	inputparams=[line.strip() for line in f.readlines()] 
	f.close()
	
	variables={}
	nextlineisinput=0
	for inputline in inputparams:
		if inputline=="" or inputline[0]=="#":
			continue
		
		if nextlineisinput!=0:
			variables[nextlineisinput]=inputline
			nextlineisinput=0
		
		if inputline[0:5]=="INPUT":
		# The next line is an input line.
			nextlineisinput=inputline.split("_")[-1]
			
	print("Variables collected:")
	for key in (variables.keys()):
		print(key)
		print(variables[key])
		print("")
	
	doextra=int(variables["doextra"])
	
	outfile=variables["outfile"]
	inputfiles=list(variables["picklenames"].split(","))
	
	
	firstwrite=0
	outputtype=int(variables["outputtype"])
	
	if outputtype==0:
		writerfile=open(outfile+".csv", "w")	
		writer = csv.writer(writerfile)
		if doextra==1:
			writerfile2=open(outfile+"_frame.csv", "w")	
			writer2 = csv.writer(writerfile2)

	else:
		print("Currently this functionality (direct counts output - output 1) is disabled until further notice. File sizes are unwieldy - 300 GB each. Transpose fails. Needs thresholding.")	
			#	writerfile=open(outfile+"_"+re.split("[./]",picklefile)[-3]+".csv", "w")	
			#	writer = csv.writer(writerfile)
			#	for gene in genedict.keys():
			#		writer.writerow([gene]+genedict[gene][2][variables["endmode"]])
			#	writerfile.close()
			#	#transposecsv(outfile+"_"+re.split("[./]",picklefile)[-3])		# Can't do big files with this.
			
	for picklefile in inputfiles:
		
		outdictionaries=genelist_m(variables["endmode"],picklefile,variables["shift"],doextra,outputtype)
		
		if outputtype==0:
			
			genedict=outdictionaries[0]
			
			if doextra==1:
				framedict=outdictionaries[1]
			else:
				framedict={}		# Just a placeholder.
			
			if firstwrite==0:
				growinggenedict=genedict
				growingframedict=framedict
				firstwrite=1
				continue
			else:
				for key in genedict.keys():
					growinggenedict[key]+=genedict[key][6:]
				for key in framedict.keys():
					growingframedict[key]+=framedict[key][1:]
		
		
	for key in genedict.keys():
		writer.writerow([key]+growinggenedict[key])					
	for key in framedict.keys():
		writer2.writerow([key]+growingframedict[key])		
			
	
	writerfile.close()		
	if doextra==1:
		writerfile2.close()
		

		
# This function counts reads on genes.
# endmode - whether you want 5' or 3' end aligned reads. "all_5" or "all_3"
# picklefile - input file of ribosome footprints
# shift - amount to shift data before counting. Typically do for P sites.
# doextra - setting this variable to 1 will output an extra file with CDS reads broken down by frame. Other extras could be added in future. 
# outputtype variable sets whether to count reads on genes (0) or just output the reads for every gene (1). 1 is currently disabled.
# NOTE: the dictionaries for gene and frame and pause are called "list" but they are dictionaries.
def genelist_m(endmode,picklefile,shift,doextra,outputtype):
	if endmode!="all_5":
		print("ERROR - right now the genelist code doesn't support anything except end5.")
	
	samplename=re.split("[./]",picklefile)[-3]
	
	# outputtype is either 0 or 1 (0 for traditional counts that are quantitated gene by gene, or 1 for just whole gene list).
	f=gzip.open(picklefile,"rb")
	footprints=pickle.load(f)
	print("Data loaded.")
	f.close()
	
	shift=int(shift)
	
	genecount=0
	outputtype=int(outputtype)
	if outputtype==1:
		
		return footprints
	else:
	
		genelist={}
		genelist["headers"]=["alias","UTR5len","CDSlen","UTR3len","transcriptlen","seq","UTR5_"+samplename,"CDS_"+samplename,"UTR3_"+samplename,"UTR5raw_"+samplename,"CDSraw_"+samplename,"UTR3raw_"+samplename]
		framelist={}
		framelist["headers"]=["alias","CDS0_"+samplename,"CDS1_"+samplename,"CDS2_"+samplename,"CDS0raw_"+samplename,"CDS1raw_"+samplename,"CDS2raw_"+samplename]
		
		# Code to find total mapped reads for picklefile.
		minoverall=1E6
		for gene in footprints.keys():
			tempgene = [i for i in footprints[gene][2][endmode] if i != 0]
			if len(tempgene)!=0:
				mingene=min(tempgene)
				if minoverall>mingene and mingene!=0:
					minoverall=mingene
		
		mappedreads=int(round(1E6/float(minoverall)))
				
		print("Assuming total mapped reads = "+str(mappedreads))
		
		for gene in footprints.keys():
			#if gene=="ENSG00000242485.5":
			#	continue	# Skip MRPL20 - weird stuff from other genome.
			
			if gene=="NM_020791.2": #This is a gene with known bad annotation.
				continue
			
			ORFstart=int(footprints[gene][3])
			UTR3start=int(footprints[gene][4])
			


			if ((ORFstart-shift)<=0) or (len(footprints[gene][2][endmode])-UTR3start)<=0:		# Gene lacks both UTRs.
				CDScount=""
				UTR5count=""
				UTR3count=""
				CDSraw=""
				UTR5raw=""
				UTR3raw=""
				CDScount_0=""
				CDScount_1=""
				CDScount_2=""
				CDScount_0raw=""
				CDScount_1raw=""
				CDScount_2raw=""
				
				continue		#Note this continue may not be ideal to compare transcriptomes.
			else:
				CDScount=sum(footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift])/((UTR3start-ORFstart)/1000)
				UTR5count=sum(footprints[gene][2][endmode][0:ORFstart-shift])/((ORFstart-shift)/1000)
				UTR3count=sum(footprints[gene][2][endmode][UTR3start-shift:-shift])/((len(footprints[gene][2][endmode])-UTR3start)/1000)
				CDSraw=int(round((CDScount*mappedreads/1E6)*((UTR3start-ORFstart)/1000)))
				UTR5raw=int(round((UTR5count*mappedreads/1E6)*((ORFstart-shift)/1000)))
				UTR3raw=int(round((UTR3count*mappedreads/1E6)*((len(footprints[gene][2][endmode])-UTR3start)/1000)))
				
				CDScount_0=sum(footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift][0::3])/((UTR3start-ORFstart)/1000)
				CDScount_1=sum(footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift][1::3])/((UTR3start-ORFstart)/1000)
				CDScount_2=sum(footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift][2::3])/((UTR3start-ORFstart)/1000)
				CDScount_0raw=int(round((CDScount_0*mappedreads/1E6)*((UTR3start-ORFstart)/1000)))
				CDScount_1raw=int(round((CDScount_1*mappedreads/1E6)*((UTR3start-ORFstart)/1000)))
				CDScount_2raw=int(round((CDScount_2*mappedreads/1E6)*((UTR3start-ORFstart)/1000)))
				
		
				
				genecount+=1
		
			genelist[gene]=[footprints[gene][0],ORFstart,UTR3start-ORFstart,len(footprints[gene][2][endmode])-UTR3start,len(footprints[gene][2][endmode]),footprints[gene][1],0,0,0,0,0,0]
			genelist[gene][6]=UTR5count
			genelist[gene][7]=CDScount
			genelist[gene][8]=UTR3count
			genelist[gene][9]=UTR5raw
			genelist[gene][10]=CDSraw
			genelist[gene][11]=UTR3raw
			
			if genelist[gene][4]>32767:
				genelist[gene][5]="toolongforExcel - lookup in fasta instead"
		
			
			# Output frame information		
			framelist[gene]=[footprints[gene][0],CDScount_0,CDScount_1,CDScount_2,CDScount_0raw,CDScount_1raw,CDScount_2raw]
		
			
		print("Genes included = "+str(genecount))	
	
		return [genelist,framelist]


# Function to find small ORFs in UTRs.
# endmode - whether you want 5' or 3' end aligned reads. "all_5" or "all_3"
# picklefile - input mapped ribosome density file
# shift - amount to shift data before analyzing it, typically for P sites.
# lengththresh - this is the threshold of amino acids. Must be longer than this. Negative is only this.
# mismatches - this is how many mismatches to allow in start codon - 1 or 0. For finding near cognates.
# This is whether to do 5'UTR or 3'UTR (5 or 3).

def dorflist_m(endmode,picklefile,shift,lengththresh,mismatches,UTR):

	if int(lengththresh)<0:		# Allow a neg length to be used absolutely.
		abslengththresh=int(lengththresh)*-1
		lengththresh=0
	
	if endmode!="all_5":
		print("ERROR - right now the genelist code doesn't support anything except end5. This is easy to fix though with a little more programming.")
	
	samplename=re.split("[./]",picklefile)[-3]
	
	f=gzip.open(picklefile,"rb")
	footprints=pickle.load(f)
	print("Data loaded "+samplename)
	f.close()
	
	shift=int(shift)
	genelist={}
	######
	genelist["headers"]=["alias","UTR5len","CDSlen","UTR3len","transcriptlen","smallorfseq","smallorfaa","sm_orfstartinUTR","sm_orfstopinUTR","sm_orf_"+samplename,"CDS_"+samplename]

	dorfcounttotal=0
	genecount=0
	for gene in footprints.keys():
		#if gene=="ENSG00000242485.5":
		#	continue	# Skip MRPL20 - weird stuff from other genome.
		
	#	if gene!="ENSG00000184009.11" and gene!="ENSG00000177606.6":		# Actin and Jun
	#		continue
		
		
		dorfcount=0
		ORFstart=int(footprints[gene][3])
		UTR3start=int(footprints[gene][4])
		
		
		if (ORFstart-shift)<0:	# Check for negative indexes:
				continue
		if (UTR3start-shift)<0:	# Check for negative indexes:
				continue
				
		if int(UTR)==3:
			
			counts=footprints[gene][2][endmode][UTR3start-shift:-shift]
			genesequence=footprints[gene][1][UTR3start:]
		else:	
			counts=footprints[gene][2][endmode][0:ORFstart-shift]
			genesequence=footprints[gene][1][shift:ORFstart]
		
		
		
		genelen=len(counts)

		i=0
		dorfcount=0		# Counting number of dorfs per gene.
		while (i+3)<genelen: 
			testmotif=str(genesequence[i:i+3])
			score=0
			if testmotif[0]!="A":
				score+=1
			if testmotif[1]!="T":
				score+=1
			if testmotif[2]!="G":
				score+=1
			if score>int(mismatches):		# How many mismatches to tolerate.
				i+=1
				continue
				
			# Find in frame stop codon.
				
			j=i+3
				
			while (j+3)<genelen:
				testmotif2=str(genesequence[j:j+3])
				if testmotif2=="TAA" or testmotif2=="TAG" or testmotif2=="TGA":
					endpos=j
					break
				j+=3
			else:
				i+=1
				continue # NO in frame stop.
				
			dorfseq=genesequence[i:j+3]
			dorfaa=Seq.Seq(dorfseq[:-3]).translate()
			dorflen=len(dorfaa)
			if dorflen<int(lengththresh) or (lengththresh==0 and dorflen!=abslengththresh):
				i+=1
				continue
				
			if ((ORFstart-shift)<=0):
				CDScount=-10
			else:
				CDScount=sum(footprints[gene][2][endmode][ORFstart-shift:UTR3start-shift][0::3])/((UTR3start-ORFstart)/1000)
				
			# Note we don't shift since this is already shifted. Was error before.
			dorfreads=sum(counts[i:j+3][0::3])/(((dorflen+1)*3)/1000)
						
						
						
			dorfkey=gene+"_smorf_"+str(dorfcount)
		
			
			genelist[dorfkey]=[footprints[gene][0],ORFstart,UTR3start-ORFstart,len(footprints[gene][2][endmode])-UTR3start,len(footprints[gene][2][endmode]),dorfseq,dorfaa,i,j,0,0]
			genelist[dorfkey][9]=dorfreads
			genelist[dorfkey][10]=CDScount
			dorfcount+=1	
			dorfcounttotal+=1
			i+=1
	
		genecount+=1
		if genecount%1000==0:
			print("Genes finished so far="+str(genecount))		
	print("Total dorfs reported = "+str(dorfcounttotal))	
	return genelist


	

# Wrapper function for finding dORFs.	
def dorflist_m_wrapper(txtinfile):
	tossedgenes=0
	f=open(txtinfile)
	inputparams=[line.strip() for line in f.readlines()] 
	f.close()
	
	variables={}
	nextlineisinput=0
	for inputline in inputparams:
		if inputline=="" or inputline[0]=="#":
			continue
		
		if nextlineisinput!=0:
			variables[nextlineisinput]=inputline
			nextlineisinput=0
		
		if inputline[0:5]=="INPUT":
		# The next line is an input line.
			nextlineisinput=inputline.split("_")[-1]
			
	print("Variables collected:")
	for key in (variables.keys()):
		print(key)
		print(variables[key])
		print("")
	
	outfile=variables["outfile"]
	inputfiles=list(variables["picklenames"].split(","))
	
	firstwrite=0
	
	
	writerfile=open(outfile+".csv", "w")	
	writer = csv.writer(writerfile)
	
	
	for picklefile in inputfiles:
		dorfdict=dorflist_m(variables["endmode"],picklefile,variables["shift"],variables["lengththresh"],variables["mismatches"],variables["UTR"])		
		
		if firstwrite==0:
			growinggenedict=dorfdict
			firstwrite=1
			continue
		else:
			for key in dorfdict.keys():
				growinggenedict[key]+=dorfdict[key][9:]
		
	startkey="startkey"
	
	towrite="towrite"
	for key in growinggenedict.keys():
		tossedgenescounter=0
		# 0s for all dorfs are gone.
		for i in range(len(inputfiles)):
			if growinggenedict[key][9+i*2]!=0:
				tossedgenescounter=1
				break
		if tossedgenescounter==0:	
			tossedgenes+=1
			continue
		
		
		# CODE TO TAKE LONGEST ISOFORM ONLY. 
		keyrootnum=key.find("_smorf")
		testkey=key[0:keyrootnum]
		if testkey!=startkey:		# New gene (new key).
			startkey=testkey
			currentstop={}
			if towrite!="towrite":
				writer.writerow(towrite)		# Write out old stop codon.
			
			currentstop[growinggenedict[key][8]]=0		# 0 is just a placeholder.
			towrite=[key]+growinggenedict[key]		
			continue				
		if growinggenedict[key][8] in currentstop:	# Check if we are on a new stop codon. If yes, check if length is better and save it.
			if variables["smallest"]==1:		# Check if we are looking for smaller or longer ORFs. 
				if growinggenedict[key][7] > towrite[8]:		# > signals shortest
					towrite=[key]+growinggenedict[key]
				elif growinggenedict[key][7] == towrite[8]:	#Same length should not occur (not possible).
					print("Error")
					exit()
			else:
				if growinggenedict[key][7] < towrite[8]:		# < signals longest
					towrite=[key]+growinggenedict[key]	
				elif growinggenedict[key][7] == towrite[8]: #Same length should not occur (not possible).
					print("Error")
					exit()
	
		else:	
			writer.writerow(towrite)					# New stop codon, get ready to write out the best option and reset.
			currentstop[growinggenedict[key][8]]=0	
			towrite=[key]+growinggenedict[key]
			
	writer.writerow(towrite)							# Write out very last one.



	writerfile.close()
	print("Genes tossed since no reads in any sample.")
	print(tossedgenes)
	
	
	
