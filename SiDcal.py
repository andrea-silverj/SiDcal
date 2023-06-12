#!/usr/bin/env python

#######################################################################################
#######################--SiDcal - Similarity Index Calculator--########################
#######################################################################################
#######################################################################################

#Author: Andrea Silverj, MSc in Evolutionary Biology - PhD student at C3A|CIBIO, University of Trento, Italy
#Version: 0.9 Beta
#Year: 2020
#Contacts: andrea.silverj@gmail.com; andrea.silverj@unitn.it

#######################################################################################
#######################################################################################

import argparse,sys,os

arguments=sys.argv
help_keys=['-h', '--help']

#########
###CLI###

if len(arguments) > 1:
	if any(hlp in arguments for hlp in help_keys):
		print("########################################################################\n##############################--|SiDcal|--##############################\nVersion: 0.9 Beta\nYear: 2020\n\nAuthor:\nAndrea Silverj, MSc in Evolutionary Biology\nPhD student at the University of Trento, Italy\nFEM-C3A, Rota-Stabelli Lab | CIBIO CM, Segata Lab\nContacts: andrea.silverj@unitn.it\n")
	else:
		pass

parser=argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="Hello dear User, and welcome to SiDcal!\nThis program calculates the similarity between two codon usage tables,\nusing the method described in the paper of Zhou et al., 2013\n\nInterpretation of the values:\nif codon_table_A=codon_table_B, SiD is equal to 0\n########################################################################")

parser.add_argument("-a", help="|e.g., a virus codon table ", metavar='codon_table_A', required=False)
parser.add_argument("-b", help="|e.g., a host codon table", metavar='codon_table_B', required=False)
parser.add_argument("-ex", action='store_true', help="|Print an example of the input format", required=False)
parser.add_argument("-ref", action='store_true', help="|Print the full reference of the Zhou et al.'s paper", required=False)

#parser.add_argument("-rscu", help="Save the 2 resulting RSCU tables in 2 separated files. STILL TO BE IMPLEMENTED", metavar='true|false', required=False, default=False)
#parser.add_argument("-vector", help="Save the 2 tables in 2 separated files in a vector format. STILL TO BE IMPLEMENTED", metavar='true|false', required=False, default=False)

args=parser.parse_args()
arg_parsed=vars(args)

#print(arg_parsed)

codon_tab_input_a=arg_parsed.get("a")
codon_tab_input_b=arg_parsed.get("b")
ref_option=arg_parsed.get("ref")
example_option=arg_parsed.get("ex")

if codon_tab_input_a==None and codon_tab_input_b==None and ref_option==False and example_option==False:
	print('Hello User! Type "SiDcal.py -h" or "SiDcal.py --help" to see how to use this program.')
else:
	pass	

if example_option is True:
	print("This is an example of the input format:\n\nTTT 0.73 (413)  TCT 0.97 (382)  TAT 0.77 (515)  TGT 0.66 (370)\nTTC 1.27 (721)  TCC 1.01 (398)  TAC 1.23 (831)  TGC 1.34 (751)\nTTA 0.60 (295)  TCA 1.33 (523)  TAA 1.50 (5)  TGA 0.00 (0)\nTTG 0.81 (400)  TCG 0.84 (331)  TAG 1.50 (5)  TGG 1.00 (400)\n\nCTT 0.57 (282)  CCT 0.85 (475)  CAT 0.65 (329)  CGT 0.39 (135)\nCTC 0.57 (279)  CCC 0.72 (404)  CAC 1.35 (679)  CGC 1.03 (362)\nCTA 1.39 (681)  CCA 1.18 (659)  CAA 0.88 (600)  CGA 0.63 (220)\nCTG 2.06 (1012)  CCG 1.25 (698)  CAG 1.12 (769)  CGG 0.45 (158)\n\nATT 0.70 (391)  ACT 0.79 (530)  AAT 0.74 (552)  AGT 0.64 (252)\nATC 1.31 (726)  ACC 0.98 (657)  AAC 1.26 (936)  AGC 1.22 (479)\nATA 0.99 (549)  ACA 1.37 (920)  AAA 1.01 (1097)  AGA 2.08 (727)\nATG 1.00 (913)  ACG 0.85 (571)  AAG 0.99 (1069)  AGG 1.42 (497)\n\nGTT 0.52 (389)  GCT 0.65 (501)  GAT 0.59 (538)  GGT 0.69 (364)\nGTC 1.14 (848)  GCC 1.11 (856)  GAC 1.41 (1295)  GGC 0.96 (501)\nGTA 0.94 (702)  GCA 1.45 (1118)  GAA 1.10 (1198)  GGA 1.27 (663)\nGTG 1.40 (1039)  GCG 0.79 (610)  GAG 0.90 (982)  GGG 1.08 (567)\n")
else:
	pass	

if ref_option is True:
	print("This program implements a measure of codon similarity which has been developed in the following paper:\n\nZhou J-h, Zhang J, Sun D-j, Ma Q, Chen H-t, et al. (2013)\nThe Distribution of Synonymous Codon Choice in the Translation Initiation Region of Dengue Virus.\nPLoS ONE 8(10): e77239. doi:10.1371/journal.pone.0077239\n\nPlease, be a gentleman and cite the authors if you use this program.")
else:
	pass	
	
#######################################################
###Import modules and data and calculate RSCU values###

if codon_tab_input_a!=None and codon_tab_input_b!=None:
	import re
	import math
else:
	pass
		
amino_acids ={'Val':'V', 'Ile':'I', 'Leu':'L', 'Glu':'E', 'Gln':'Q', 'Asp':'D', 'Asn':'N', 'His':'H', 'Trp':'W', 'Phe':'F', 'Tyr':'Y', 'Arg':'R', 'Lys':'K', 'Ser':'S', 'Thr':'T', 'Met':'M', 'Ala':'A', 'Gly':'G', 'Pro':'P', 'Cys':'C'}

aa_qcodons={'Phe':2, 'Leu':6, 'Ile':3, 'Met':1, 'Val':4, 'Ser':6, 'Pro':4, 'Thr':4, 'Ala':4, 'Tyr':2, 'Ochre':1, 'Amber':1, 'His':2, 'Gln':2, 'Asn':2, 'Lys':2, 'Asp':2, 'Glu':2, 'Cys':2, 'Opal':1, 'Trp':1, 'Arg':6, 'Gly':4}


def get_counts_codons(input_k):
	cod_tab1=input_k
	cod_tab_o=open(cod_tab1, "r")

	cod_tab1_l=list(cod_tab_o)

	cod_tab1_s=[]
	for ele in cod_tab1_l:
		cod_tab1_s+=ele.strip("\n").split("\t")

	cod_num_freq_no1=list(filter(None, cod_tab1_s))

	cod_num_freq_no2=[]
	for i in cod_num_freq_no1:
		cod_num_freq_no2+=i.split(" ")

	cod_num_freq1=list(filter(None, cod_num_freq_no2))
	cod_num_freq=list(filter(lambda k: 'TER' not in k, cod_num_freq1))

	indexes_num_cod=[]
	for i in range(1, len(cod_num_freq), 3):
		indexes_num_cod+=[i]

	cod_freq=[]
	for cod_data in cod_num_freq:
		if cod_num_freq.index(cod_data) not in indexes_num_cod:
			cod_freq+=[cod_data]

	cod_names=cod_freq[0::2]
	cod_freq_values=cod_freq[1::2]

	cod_w='\t'.join([str(x) for x in cod_names])
	cod_val_w='\t'.join([str(x) for x in cod_freq_values])
	return cod_freq_values

def sum_codon_for_aa(codon_counts_k):
	dict_counts_codon_for_aa={'Phe':codon_counts_k[0]+codon_counts_k[4],'Leu':codon_counts_k[8]+codon_counts_k[12]+codon_counts_k[16]+codon_counts_k[20]+codon_counts_k[24]+codon_counts_k[28], 'Ile':codon_counts_k[32]+codon_counts_k[36]+codon_counts_k[40], 'Met':codon_counts_k[44], 'Val':codon_counts_k[48]+codon_counts_k[52]+codon_counts_k[56]+codon_counts_k[60], 'Ser':codon_counts_k[1]+codon_counts_k[5]+codon_counts_k[9]+codon_counts_k[13]+codon_counts_k[35]+codon_counts_k[39], 'Pro':codon_counts_k[17]+codon_counts_k[21]+codon_counts_k[25]+codon_counts_k[29], 'Thr':codon_counts_k[33]+codon_counts_k[37]+codon_counts_k[41]+codon_counts_k[45], 'Ala':codon_counts_k[49]+codon_counts_k[53]+codon_counts_k[57]+codon_counts_k[61], 'Tyr':+codon_counts_k[2]+codon_counts_k[6], 'Ochre':codon_counts_k[10], 'Amber':codon_counts_k[14], 'His':codon_counts_k[18]+codon_counts_k[22], 'Gln':codon_counts_k[26]+codon_counts_k[30], 'Asn':codon_counts_k[34]+codon_counts_k[38], 'Lys':codon_counts_k[42]+codon_counts_k[46], 'Asp':codon_counts_k[50]+codon_counts_k[54], 'Glu':codon_counts_k[58]+codon_counts_k[62], 'Cys':codon_counts_k[3]+codon_counts_k[7], 'Opal':codon_counts_k[11], 'Trp':codon_counts_k[15], 'Arg':codon_counts_k[19]+codon_counts_k[23]+codon_counts_k[27]+codon_counts_k[31]+codon_counts_k[43]+codon_counts_k[47], 'Gly':codon_counts_k[51]+codon_counts_k[55]+codon_counts_k[59]+codon_counts_k[63]}	
	return dict_counts_codon_for_aa

def rscu_calculator(codon_counts_k, counts_codon_for_aa_k, aa_qcodons):
	if codon_counts_k[10] > 0:
		ochre_count=round(float(codon_counts_k[10]/counts_codon_for_aa_k.get('Ochre')),3)
	else:
		ochre_count=0
	if codon_counts_k[11] > 0:
		opal_count=round(float(codon_counts_k[11]/counts_codon_for_aa_k.get('Opal')),3)
	else:
		opal_count=0
	if codon_counts_k[14] > 0:
		amber_count=round(float(codon_counts_k[14]/counts_codon_for_aa_k.get('Amber')),3)
	else:
		amber_count=0	
	if codon_counts_k[44] > 0:
		met_count=round(float(codon_counts_k[44]/counts_codon_for_aa_k.get('Met')),3)
	else:
		met_count=0	
	rscu_tab={'TTT':round(round((codon_counts_k[0]/counts_codon_for_aa_k.get('Phe')),3)*aa_qcodons.get('Phe'), 3), 'TCT':round(round((codon_counts_k[1]/counts_codon_for_aa_k.get('Ser')),3)*aa_qcodons.get('Ser'), 3), 'TAT':round(round((codon_counts_k[2]/counts_codon_for_aa_k.get('Tyr')),3)*aa_qcodons.get('Tyr'), 3), 'TGT':round(round((codon_counts_k[3]/counts_codon_for_aa_k.get('Cys')),3)*aa_qcodons.get('Cys'), 3), 'TTC':round(round((codon_counts_k[4]/counts_codon_for_aa_k.get('Phe')),3)*aa_qcodons.get('Phe'), 3), 'TCC':round(round((codon_counts_k[5]/counts_codon_for_aa_k.get('Ser')),3)*aa_qcodons.get('Ser'), 3), 'TAC':round(round((codon_counts_k[6]/counts_codon_for_aa_k.get('Tyr')),3)*aa_qcodons.get('Tyr'), 3), 'TGC':round(round((codon_counts_k[7]/counts_codon_for_aa_k.get('Cys')),3)*aa_qcodons.get('Cys'), 3), 'TTA':round(round((codon_counts_k[8]/counts_codon_for_aa_k.get('Leu')),3)*aa_qcodons.get('Leu'), 3), 'TCA':round(round((codon_counts_k[9]/counts_codon_for_aa_k.get('Ser')),3)*aa_qcodons.get('Ser'), 3), 'TAA':round(ochre_count*(aa_qcodons.get('Ochre')+aa_qcodons.get('Amber')+aa_qcodons.get('Opal')), 3), 'TGA':round(opal_count*(aa_qcodons.get('Ochre')+aa_qcodons.get('Amber')+aa_qcodons.get('Opal')), 3), 'TTG':round(round((codon_counts_k[12]/counts_codon_for_aa_k.get('Leu')),3)*aa_qcodons.get('Leu'), 3), 'TCG':round(round((codon_counts_k[13]/counts_codon_for_aa_k.get('Ser')),3)*aa_qcodons.get('Ser'), 3), 'TAG':round(amber_count*(aa_qcodons.get('Ochre')+aa_qcodons.get('Amber')+aa_qcodons.get('Opal')), 3), 'TGG':round(round((codon_counts_k[15]/counts_codon_for_aa_k.get('Trp')),3)*aa_qcodons.get('Trp'), 3), 'CTT':round(round((codon_counts_k[16]/counts_codon_for_aa_k.get('Leu')),3)*aa_qcodons.get('Leu'), 3), 'CCT':round(round((codon_counts_k[17]/counts_codon_for_aa_k.get('Pro')),3)*aa_qcodons.get('Pro'), 3), 'CAT':round(round((codon_counts_k[18]/counts_codon_for_aa_k.get('His')),3)*aa_qcodons.get('His'), 3), 'CGT':round(round((codon_counts_k[19]/counts_codon_for_aa_k.get('Arg')),3)*aa_qcodons.get('Arg'), 3), 'CTC':round(round((codon_counts_k[20]/counts_codon_for_aa_k.get('Leu')),3)*aa_qcodons.get('Leu'), 3), 'CCC':round(round((codon_counts_k[21]/counts_codon_for_aa_k.get('Pro')),3)*aa_qcodons.get('Pro'), 3), 'CAC':round(round((codon_counts_k[22]/counts_codon_for_aa_k.get('His')),3)*aa_qcodons.get('His'), 3), 'CGC':round(round((codon_counts_k[23]/counts_codon_for_aa_k.get('Arg')),3)*aa_qcodons.get('Arg'), 3), 'CTA':round(round((codon_counts_k[24]/counts_codon_for_aa_k.get('Leu')),3)*aa_qcodons.get('Leu'), 3), 'CCA':round(round((codon_counts_k[25]/counts_codon_for_aa_k.get('Pro')),3)*aa_qcodons.get('Pro'), 3), 'CAA':round(round((codon_counts_k[26]/counts_codon_for_aa_k.get('Gln')),3)*aa_qcodons.get('Gln'), 3), 'CGA':round(round((codon_counts_k[27]/counts_codon_for_aa_k.get('Arg')),3)*aa_qcodons.get('Arg'), 3), 'CTG':round(round((codon_counts_k[28]/counts_codon_for_aa_k.get('Leu')),3)*aa_qcodons.get('Leu'), 3), 'CCG':round(round((codon_counts_k[29]/counts_codon_for_aa_k.get('Pro')),3)*aa_qcodons.get('Pro'), 3), 'CAG':round(round((codon_counts_k[30]/counts_codon_for_aa_k.get('Gln')),3)*aa_qcodons.get('Gln'), 3), 'CGG':round(round((codon_counts_k[31]/counts_codon_for_aa_k.get('Arg')),3)*aa_qcodons.get('Arg'), 3), 'ATT':round(round((codon_counts_k[32]/counts_codon_for_aa_k.get('Ile')),3)*aa_qcodons.get('Ile'), 3), 'ACT':round(round((codon_counts_k[33]/counts_codon_for_aa_k.get('Thr')),3)*aa_qcodons.get('Thr'), 3), 'AAT':round(round((codon_counts_k[34]/counts_codon_for_aa_k.get('Asn')),3)*aa_qcodons.get('Asn'), 3), 'AGT':round(round((codon_counts_k[35]/counts_codon_for_aa_k.get('Ser')),3)*aa_qcodons.get('Ser'), 3), 'ATC':round(round((codon_counts_k[36]/counts_codon_for_aa_k.get('Ile')),3)*aa_qcodons.get('Ile'), 3), 'ACC':round(round((codon_counts_k[37]/counts_codon_for_aa_k.get('Thr')),3)*aa_qcodons.get('Thr'), 3), 'AAC':round(round((codon_counts_k[38]/counts_codon_for_aa_k.get('Asn')),3)*aa_qcodons.get('Asn'), 3), 'AGC':round(round((codon_counts_k[39]/counts_codon_for_aa_k.get('Ser')),3)*aa_qcodons.get('Ser'), 3), 'ATA':round(round((codon_counts_k[40]/counts_codon_for_aa_k.get('Ile')),3)*aa_qcodons.get('Ile'), 3), 'ACA':round(round((codon_counts_k[41]/counts_codon_for_aa_k.get('Thr')),3)*aa_qcodons.get('Thr'), 3), 'AAA':round(round((codon_counts_k[42]/counts_codon_for_aa_k.get('Lys')),3)*aa_qcodons.get('Lys'), 3), 'AGA':round(round((codon_counts_k[43]/counts_codon_for_aa_k.get('Arg')),3)*aa_qcodons.get('Arg'), 3), 'ATG':round(met_count*aa_qcodons.get('Met'), 3), 'ACG':round(round((codon_counts_k[45]/counts_codon_for_aa_k.get('Thr')),3)*aa_qcodons.get('Thr'), 3), 'AAG':round(round((codon_counts_k[46]/counts_codon_for_aa_k.get('Lys')),3)*aa_qcodons.get('Lys'), 3), 'AGG':round(round((codon_counts_k[47]/counts_codon_for_aa_k.get('Arg')),3)*aa_qcodons.get('Arg'), 3), 'GTT':round(round((codon_counts_k[48]/counts_codon_for_aa_k.get('Val')),3)*aa_qcodons.get('Val'), 3), 'GCT':round(round((codon_counts_k[49]/counts_codon_for_aa_k.get('Ala')),3)*aa_qcodons.get('Ala'), 3), 'GAT':round(round((codon_counts_k[50]/counts_codon_for_aa_k.get('Asp')),3)*aa_qcodons.get('Asp'), 3), 'GGT':round(round((codon_counts_k[51]/counts_codon_for_aa_k.get('Gly')),3)*aa_qcodons.get('Gly'), 3), 'GTC':round(round((codon_counts_k[52]/counts_codon_for_aa_k.get('Val')),3)*aa_qcodons.get('Val'), 3), 'GCC':round(round((codon_counts_k[53]/counts_codon_for_aa_k.get('Ala')),3)*aa_qcodons.get('Ala'), 3), 'GAC':round(round((codon_counts_k[54]/counts_codon_for_aa_k.get('Asp')),3)*aa_qcodons.get('Asp'), 3), 'GGC':round(round((codon_counts_k[55]/counts_codon_for_aa_k.get('Gly')),3)*aa_qcodons.get('Gly'), 3), 'GTA':round(round((codon_counts_k[56]/counts_codon_for_aa_k.get('Val')),3)*aa_qcodons.get('Val'), 3), 'GCA':round(round((codon_counts_k[57]/counts_codon_for_aa_k.get('Ala')),3)*aa_qcodons.get('Ala'), 3), 'GAA':round(round((codon_counts_k[58]/counts_codon_for_aa_k.get('Glu')),3)*aa_qcodons.get('Glu'), 3), 'GGA':round(round((codon_counts_k[59]/counts_codon_for_aa_k.get('Gly')),3)*aa_qcodons.get('Gly'), 3), 'GTG':round(round((codon_counts_k[60]/counts_codon_for_aa_k.get('Val')),3)*aa_qcodons.get('Val'), 3), 'GCG':round(round((codon_counts_k[61]/counts_codon_for_aa_k.get('Ala')),3)*aa_qcodons.get('Ala'), 3), 'GAG':round(round((codon_counts_k[62]/counts_codon_for_aa_k.get('Glu')),3)*aa_qcodons.get('Glu'), 3), 'GGG':round(round((codon_counts_k[63]/counts_codon_for_aa_k.get('Gly')),3)*aa_qcodons.get('Gly'), 3)}
	return rscu_tab

def power_2(my_list_ai):
	return [round((x**2),3) for x in my_list_ai]

def sid_calculator(rscu_k,rscu_l):
	entriesToRemove = ('TAA', 'TGA', 'TAG', 'TGG', 'ATG')
	for i in entriesToRemove:
		rscu_k.pop(i, None)
		rscu_l.pop(i, None)
	rsk=rscu_k.values()
	rsl=rscu_l.values()
	ai_x_bi=[round(float(k*l),3) for k, l in zip(rsk, rsl)]
	sum_ai_x_bi=round(sum(ai_x_bi),3)
	ai_2=power_2(rsk)
	bi_2=power_2(rsl)
	sum_ai_2=round(sum(ai_2),3)
	sum_bi_2=round(sum(bi_2),3)
	sum_ai_2_x_sum_bi_2=round((sum_ai_2*sum_bi_2),3)
	square_root_sum_ai_2_x_sum_bi_2=round((math.sqrt(sum_ai_2_x_sum_bi_2)),3)
	r_ab=round((sum_ai_x_bi/square_root_sum_ai_2_x_sum_bi_2),3)
	D_AB=round(((1-r_ab)/2),3)
	return D_AB

if codon_tab_input_a != None:
	codon_counts_a=[]
	for i in get_counts_codons(codon_tab_input_a):
		codon_counts_a+=[int(i.strip("(").strip(")"))]

	counts_codon_for_aa_a=sum_codon_for_aa(codon_counts_a)
	rscu_a=rscu_calculator(codon_counts_a,counts_codon_for_aa_a,aa_qcodons)
else:
	pass	

if codon_tab_input_b != None:
	codon_counts_b=[]
	for i in get_counts_codons(codon_tab_input_b):
		codon_counts_b+=[int(i.strip("(").strip(")"))]

	counts_codon_for_aa_b=sum_codon_for_aa(codon_counts_b)
	rscu_b=rscu_calculator(codon_counts_b,counts_codon_for_aa_b,aa_qcodons)
else:
	pass	

if codon_tab_input_a!=None and codon_tab_input_b!=None and rscu_a!=None and rscu_b!=None:
	print("SiD: "+str(sid_calculator(rscu_a,rscu_b)))
else:
	pass	
