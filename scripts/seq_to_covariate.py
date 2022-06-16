# Object for capturing mutations found in SARS-CoV-2 genomic sequencing data

import sys
import re
from Bio import SeqIO
from collections import OrderedDict, defaultdict
from typing import Dict, List, Tuple
import edlib
import reference_aa

class SeqToCovariate(object):

	DNA_TO_AA = {"TTT": "F","TTC": "F","TTA": "L","TTG": "L","CTT": "L","CTC": "L","CTA": "L","CTG": "L","ATT": "I","ATC": "I",
	"ATA": "I","ATG": "M","GTT": "V","GTC": "V","GTA": "V","GTG": "V","TCT": "S","TCC": "S","TCA": "S","TCG": "S","CCT": "P",
	"CCC": "P","CCA": "P","CCG": "P","ACT": "T","ACC": "T","ACA": "T","ACG": "T","GCT": "A","GCC": "A","GCA": "A","GCG": "A","TAT": "Y",
	"TAC": "Y","TAA": 'STOP',"TAG": 'STOP',"CAT": "H","CAC": "H","CAA": "Q","CAG": "Q","AAT": "N","AAC": "N","AAA": "K","AAG": "K","GAT": "D",
	"GAC": "D","GAA": "E","GAG": "E","TGT": "C","TGC": "C","TGA": 'STOP',"TGG": "W","CGT": "R","CGC": "R","CGA": "R","CGG": "R","AGT": "S",
	"AGC": "S","AGA": "R","AGG": "R","GGT": "G","GGC": "G","GGA": "G","GGG": "G"}

	ORF1a = {"nsp1": (1, 180),"nsp2": (181, 818),"nsp3": (819, 2763), "nsp4": (2764, 3263),"nsp5": (3264, 3569),"nsp6": (3570, 3859),
	"nsp7": (3860, 3942),"nsp8": (3943, 4140),"nsp9": (4141, 4253),"nsp10": (4254, 4392), "nsp11": (4393, 4405)}

	ORF1b = {"nsp12": (1, 932),"nsp13": (933, 1533),"nsp14": (1534, 2060),"nsp15": (2061, 2406),"nsp16": (2407, 2704)}

	# Adapted from https://github.com/nextstrain/ncov/blob/50ceffa/defaults/annotation.gff with slight
	# change to ORF1a/ORF1b based on Vigor4 annotations from the Virus Pathogen Resource
	# Also adapted from Fritz Obermeyer (Broad)
	def _():
		annotation_tsv = """\
		seqname	source	feature	start	end	score	strand	frame	attribute
		.	.	gene	26245	26472	.	+	.	gene_name "E"
		.	.	gene	26523	27191	.	+	.	gene_name "M"
		.	.	gene	28274	29533	.	+	.	gene_name "N"
		.	.	gene	29558	29674	.	+	.	gene_name "ORF10"
		.	.	gene	28734	28955	.	+	.	gene_name "ORF14"
		.	.	gene	266	13483	.	+	.	gene_name "ORF1a"
		.	.	gene	13468	21555	.	+	.	gene_name "ORF1b"
		.	.	gene	25393	26220	.	+	.	gene_name "ORF3a"
		.	.	gene	27202	27387	.	+	.	gene_name "ORF6"
		.	.	gene	27394	27759	.	+	.	gene_name "ORF7a"
		.	.	gene	27756	27887	.	+	.	gene_name "ORF7b"
		.	.	gene	27894	28259	.	+	.	gene_name "ORF8"
		.	.	gene	28284	28577	.	+	.	gene_name "ORF9b"
		.	.	gene	21563	25384	.	+	.	gene_name "Spike"
		"""
		genes = []
		rows = annotation_tsv.split("\n")
		header, rows = rows[0].split("\t"), rows[1:]
		rows = rows[:-1]
		for row in rows:
			if row:
				row = dict(zip(header, row.split("\t")))
				gene_name = row["attribute"].split('"')[1]
				start = int(row["start"])
				end = int(row["end"])
				genes.append(((start, end), gene_name))
		genes.sort()
		return(OrderedDict((gene_name, pos) for pos, gene_name in genes))


	# This maps gene name to the nucleotide position in the genome,
	# as measured in the original Wuhan virus.  It is similar to what Vigor4 does
	GENE_TO_POSITION: Dict[str, Tuple[int, int]] = _()


	# Compute the full covariate sequence (constellation of mutations)
	# Takes in a dictionary with gene names as keys and AA sequences of the query 
	# sequqnce as values, like {'S': 'MFVFLVLLPLVSSQCVNL...'}
	@staticmethod
	def aa_mutations(aa_seqs):

		# Using ORF1a/1b dictionary, get the ORF1a/1b structural protein containing the variant
		def get_orf1ab_structural(pos, orf):
			if orf == 'ORF1a':
				for prot in SeqToCovariate.ORF1a.keys():
					if SeqToCovariate.ORF1a[prot][0] <= pos <= SeqToCovariate.ORF1a[prot][1]:
						new_pos = (pos - SeqToCovariate.ORF1a[prot][0]) + 1
						return(prot, new_pos)
			
			elif orf == 'ORF1b':
				for prot in SeqToCovariate.ORF1b.keys():
					if SeqToCovariate.ORF1b[prot][0] <= pos <= SeqToCovariate.ORF1b[prot][1]:
						new_pos = (pos - SeqToCovariate.ORF1b[prot][0]) + 10 # Account for ORF1a/1b frameshift
						return(prot, new_pos)

		covariate = []
		ambiguous = []
		protein_covariates = {}
		for gene in aa_seqs.keys():
			prot_covariate = defaultdict(list)
			aa_seq = aa_seqs[gene]
			ref_seq = reference_aa.REFERENCE_AA[gene]
			alignment = edlib.align(aa_seq, ref_seq, task = 'path')
			cigar = alignment['cigar']
			if (gene == "ORF1a" or gene == "ORF1b"):
				is_orf = True
				orf = gene
			else:
				is_orf = False

			# Parse cigar only if it has mutations
			if any((c in set('IDX')) for c in cigar):
				alignments = edlib.getNiceAlignment(alignment, aa_seq, ref_seq)
				seq_aligned = alignments['query_aligned']
				insertions = 0
				start = 1
				mutations = []
				for num_muts, idx in re.findall('(\d+)([IDX=])', cigar):
					num = int(num_muts)
					if idx in 'DX':
						if is_orf:
							gene, new_pos = get_orf1ab_structural(start, orf)
						if idx in 'D':
							for i in range(num):
								if is_orf:
									del_mut = gene+"_"+ref_seq[start + i - 1]+str(new_pos + i)+"-"
									prot_covariate[gene].append(del_mut)
									mutations.append(del_mut)
								else:
									del_mut = gene+"_"+ref_seq[start + i - 1]+str(start + i)+"-"
									# Weird glitch in edlib alignment, needs to be manually changed
									if del_mut == "N_G204-":
										del_mut = "N_G204R"
									prot_covariate[gene].append(del_mut)
									mutations.append(del_mut)
						elif idx in 'X':
							for i in range(num):
								if is_orf:
									aa_sub = seq_aligned[(start + i + insertions - 1)]
									if (aa_sub == 'X'): # Ambiguous, collect position then skip
										ambiguous.append((gene, new_pos + i))
										continue
									aa_mut = gene+"_"+ref_seq[start + i - 1]+str(new_pos + i)+aa_sub
									prot_covariate[gene].append(aa_mut)
									mutations.append(aa_mut)
								else:
									aa_sub = seq_aligned[(start + i + insertions - 1)]
									if (aa_sub == 'X'):
										ambiguous.append((gene, start + i))
										continue
									aa_mut = gene+"_"+ref_seq[start + i - 1]+str(start + i)+aa_sub
									prot_covariate[gene].append(aa_mut)
									mutations.append(aa_mut)
					if idx in 'I':
						good_ins = False
						if is_orf:
							gene, new_pos = get_orf1ab_structural(start, orf)
							ins_mut = gene+"_"+"ins"+str(new_pos-1)
						else:
							ins_mut = gene+"_"+"ins"+str(start-1)
						for i in range(num):
							if (seq_aligned[(start + i + insertions - 1)] == "X"):
								continue
							ins_mut = ins_mut+seq_aligned[(start + i + insertions - 1)]
							good_ins = True
						if good_ins:
							# Weird glitch in edlib alignment, needs to be manually changed
							if ins_mut == "N_ins202K":
								ins_mut = "N_R203K"
							prot_covariate[gene].append(ins_mut)
							mutations.append(ins_mut)
						insertions += num
					else:
						start += num

				if mutations != []:
					mutations = ','.join(mutations)
					covariate.append(mutations)
					for prot in prot_covariate.keys():
						protein_covariates[prot] = ','.join(prot_covariate[prot])

		
		# Concatentate the ambiguous chunks together as a representative string
		# denoting protein and positions for which ambiguity is located
		if ambiguous != []:
			ambig_chunks = []
			prev_gene = ambiguous[0][0]
			prev_pos = ambiguous[0][1]
			a = prev_gene+"_"+str(prev_pos)+"-"
			ambiguous = ambiguous[1:]
			for ambig in ambiguous:
				ambig_gene = ambig[0]
				ambig_pos = ambig[1]
				if (ambig_gene == prev_gene and ambig_pos == prev_pos+1):
					prev_pos = ambig_pos
				else:
					a = a+str(prev_pos)
					ambig_chunks.append(a)
					prev_gene = ambig_gene
					prev_pos = ambig_pos
					a = prev_gene+"_"+str(prev_pos)+"-"
			a = a+str(prev_pos)
			ambig_chunks.append(a)
			prev_gene = ambig_gene
			prev_pos = ambig_pos
			a = prev_gene+"_"+str(prev_pos)+"-"
			
			ambig_chunks = ','.join(ambig_chunks)
		else:
			ambig_chunks = "None"


		covariate = ','.join(covariate)
		return(covariate, protein_covariates, ambig_chunks)


	# Adapted from Fritz Obermeyer
	# Takes in a reference and a list of nuc mutations as (1234, 'A'),
	# returns the nonsynonymous and synonymous mutations as disjoint lists
	"""Paramaters:
		ref: SARS-CoV-2 reference sequence, processed fasta, just the string
		muts: Tuple of nucleotide mutations denoted as (position, alt)
	"""
	@staticmethod
	def non_and_syn_mutations(ref, muts):
		ms_by_aa = defaultdict(list)

		for m in muts: 
			position_nuc = m[0]
			mut_nuc = m[1]

			# Need to locate position on codons containing the mut
			for gene, (start, end) in SeqToCovariate.GENE_TO_POSITION.items():
				if start <= position_nuc <= end:
					position_aa = (position_nuc - start) // 3
					position_codon = (position_nuc - start) % 3
					ms_by_aa[gene, position_aa].append((position_codon, mut_nuc, position_nuc))
					break

		synonymous_muts = []
		non_syn_muts = []
		for (gene, position_aa), ms in ms_by_aa.items():
			start, end = SeqToCovariate.GENE_TO_POSITION[gene]

			# Apply mutation to determine new aa.
			pos = start + position_aa * 3
			pos -= 1  # convert from 1-based to 0-based
			old_codon = ref[pos : pos + 3]
			new_codon = list(old_codon)
			for position_codon, mut_nuc, position_nuc in ms:
				new_codon[position_codon] = mut_nuc
			new_codon = "".join(new_codon)

			old_aa = SeqToCovariate.DNA_TO_AA.get(old_codon)
			new_aa = SeqToCovariate.DNA_TO_AA.get(new_codon)
			if new_aa == old_aa:  # synonymous substitution
				for position_codon, mut_nuc, position_nuc in ms:
					syn_mut = ref[position_nuc-1]+str(position_nuc)+mut_nuc
					synonymous_muts.append(syn_mut)
			elif new_aa == None: # Ambiguity in the codon, skip it
				continue
			else:
				for position_codon, mut_nuc, position_nuc in ms:
					non_syn_mut = ref[position_nuc-1]+str(position_nuc)+mut_nuc
					non_syn_muts.append(non_syn_mut)


		synonymous_muts = ','.join(synonymous_muts)
		non_syn_muts = ','.join(non_syn_muts)

		return(synonymous_muts, non_syn_muts)


	# Compute the mutations of a SARS-CoV-2 query sequence based on alignment to the Wuhan reference
	# Takes in a reference nucleotide sequence and the query sequence, and returns non-synonymous mutations, 
	# synonymous mutations, and the covariate (constellation of amino acid substituions)
	"""Parameters:
		ref: SARS-CoV-2 reference sequence, processed fasta, just the string
		seq: SARS-CoV-2 query sequence, processed fasta, just the string
	"""
	@staticmethod
	def mutations(ref, seq):
		# Pairwise alignment to Wuhan reference, get cigar string
		alignment = edlib.align(seq, ref, task = "path")
		alignments = edlib.getNiceAlignment(alignment, seq, ref)
		cigar = alignment['cigar']
		seq_aligned = alignments['query_aligned']

		# Parse the cigar string to get nucleotide substitutions, deletions, and insertions
		insertions = 0
		start = 1
		del_positions = []
		ins_positions = []
		mut_positions = []
		nuc_mutations_all = []
		nuc_mutations_no_indels = []
		for num_muts, idx in re.findall('(\d+)([IDX=])', cigar):
			num = int(num_muts)
			if idx in 'DX':
				if idx in 'D':
					if (start < 266): 
						for i in range(num):
							del_positions.append(start + i)
					elif (start >= 266 and start <= (len(seq)+len(del_positions))):
						for i in range(num):
							del_mut = ref[start + i - 1]+str(start + i)+"del"
							nuc_mutations_all.append(del_mut)
							del_positions.append(start + i)
							mut_positions.append(start + i)
				elif idx in 'X':
					if (start >= 266):
						for i in range(num):
							nuc_mut = seq_aligned[(start + i + insertions - 1)]
							sub = ref[start + i - 1]+str(start + i)+nuc_mut
							nuc_mutations_all.append(sub)
							nuc_mutations_no_indels.append(tuple((start + i, nuc_mut)))
							mut_positions.append(start + i)
			if idx in 'I':
				ins_mut = "ins"+str(start-1)
				for i in range(num):
					ins_mut = ins_mut+seq_aligned[(start + i + insertions - 1)]
					ins_positions.append(start + i)
					mut_positions.append(start + i)
				nuc_mutations_all.append(ins_mut)
				insertions += num
			else:
				start += num

		nuc_mutations_all = ','.join(nuc_mutations_all)

		# Using the computed aligment and insertion and deletion positions, translate the query sequence and store to dict
		# This will only translate the proteins that have mutations to minimize computations.
		bad_sequence = False
		query_seqs = {}
		translated_genes = set()
		for pos in mut_positions:
			for gene, (start, end) in SeqToCovariate.GENE_TO_POSITION.items():
				# Already translated this gene, break inner loop
				if (start <= pos <= end) & (gene in translated_genes):
					break
				elif (start <= pos <= end) & (gene not in translated_genes):
					translated_genes.add(gene)
					# Compute the adjusted start and end of a gene on the query sequence
					# by accounting for insertions and deletions
					seq_start = start + (len([i for i in ins_positions if i < start]) -
						len([j for j in del_positions if j < start]))
					seq_end = end + (len([i for i in ins_positions if i < end]) -
						len([j for j in del_positions if j < end]))
					ref_prot = str(ref[(start-1) : end])
					query = str(seq[(seq_start-1) : seq_end])

					if (abs(len(ref_prot) - len(query)) % 3) == 0: # This will catch frameshifts and terminate the mutation computation
						aa_seq = ""
						for c in range(0, (len(query))-3, 3):
							codon = query[c : c+3]
							aa = SeqToCovariate.DNA_TO_AA.get(codon)
							if aa == None:
								aa_seq += 'X'
							elif aa != 'STOP':
								aa_seq += aa

						query_seqs[gene] = aa_seq
						break
					else:
						bad_sequence = True
						break
			if bad_sequence:
				break

		if bad_sequence:
			return('BAD_SEQUENCE')

		else: 
			synonymous_muts, non_syn_muts = SeqToCovariate.non_and_syn_mutations(ref, nuc_mutations_no_indels)
			covariate, protein_covariates, ambig_chunck = SeqToCovariate.aa_mutations(query_seqs)
			return([non_syn_muts, synonymous_muts, covariate, protein_covariates, ambig_chunck])


if __name__ == "__main__":

	reference = sys.argv[1]
	sequences = sys.argv[2]

	for seq_record in SeqIO.parse(reference, 'fasta'):
		ref = seq_record.seq

	for seq_record in SeqIO.parse(sequences, 'fasta'):
		seq = seq_record.seq
		all_mutations = SeqToCovariate.mutations(ref, seq)
		if all_mutations != 'BAD_SEQUENCE':
			print('\n')
			print("Non-Synonymous:", all_mutations[0])
			print('\n')
			print("Synonymous:", all_mutations[1])
			print('\n')
			print("Covariate:", all_mutations[2])
			print('\n')
			print("Protein Covariates:", all_mutations[3])
			print('\n')
			print("Ambiguous Positions:", all_mutations[4])
			print('\n')
		else:
			print(all_mutations)
			print('\n')