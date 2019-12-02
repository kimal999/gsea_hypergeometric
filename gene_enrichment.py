#!/usr/bin/env pyhton

def main():
	
	import re,os
	import numpy as np
	import scipy
	import pandas as pd
	import scipy.stats
	
	parser = argparse.ArgumentParser(description='Calculation of signature score for samples') 

	parser.add_argument('-f1','--file_1', type=str, required=True, 
		help='Gene list for comparsion or the gene signature')
	parser.add_argument('-N','--un_genes', type=int, required=False, default=45956,
                    help='FNo of genes in universe, Default value by GSEA = 45956')
	parser.add_argument('-f2','--file_2', type=str, required=False, 
		help='Gene list for universe')
	parser.add_argument('-f3','--geneset', type=str, required=False, 
		help='Gene sets in gmt format')
	parser.add_argument('-o1','--output1', type=str, required=True, 
		help='Output file with gene enrichment values')
	
	args = parser.parse_args()
	
	gene = {}
	gene_set = {}
	enriched_gene = {}
	enriched_gene_set = {}
	gene_set_no = {}
	all_ant = {}
	
	f1 = open("kegg_enriched_gsea_rppa_TGFB.txt","w")
	f2 = open("kegg_overlapped_gsea_rppa_TGFB.txt","w")
	
	if args.file_2:	
		with open('rppa_all.txt','r') as f:
		    while True:
			text = f.readline()
			if text == "":
			    break
			cols = text.rstrip().split('\t')
			all_ant[cols[0]] = 0
			#enriched_gene[cols[0]] = 0
	
	with open('args.file_1','r') as f:
	    while True:
		text = f.readline()
		if text == "":
		    break
		cols = text.rstrip().split('\t')
		gene[cols[0]] = 0
		enriched_gene[cols[0]] = 0
	
	comparision = len(gene)
	universe = int(args.un_genes)
	
	f1.write("Gene_set\tP-value\t#genes in universe\t#Genes in Gene Set\t#genes in comparison\t#Genes in Overlap\n")
	
	with open('args.geneset','r') as f:
	    while True:
		text = f.readline()
		if text == "":
		    break
		cols = text.rstrip().split('\t')
		gene_set[cols[0]] = {}
		
		gene_set_no[cols[0]] = 0
		for i in range(2,len(cols)):
		    if cols[i] in gene:
			(gene_set[cols[0]])[cols[i]] = cols[i]
			gene_set_no[cols[0]] = gene_set_no[cols[0]] + 1
			gene[cols[i]] = gene[cols[i]] + 1
			
		if args.file_2:	
			ant = 0        
			for i in range(2,len(cols)):
			    if cols[i] in all_ant:
				ant = ant + 1
		esle:
			ant = (len(cols)-2)
			     
		if gene_set_no[cols[0]] > 0:
		    list1 = [universe,ant]
		    list2 = [comparision, gene_set_no[cols[0]]]
		    array = np.array([list1, list2], dtype = np.int32)
		    
		    #print array
		    
		    #results = scipy.stats.fisher_exact(array)
		    results = scipy.stats.chi2_contingency(array)
		    
		    if results[1] <= 1:
			f1.write(cols[0]+"\t"+str(results[1])+"\t"+str(universe)+"\t"+str(ant)+"\t"+str(comparision)+"\t"+str(gene_set_no[cols[0]])+"\n")
			enriched_gene_set[cols[0]] = {}
			for i in range(2,len(cols)):
			    if cols[i] in gene:
				(enriched_gene_set[cols[0]])[cols[i]] = cols[i]
				enriched_gene[cols[i]] = enriched_gene[cols[i]] + 1
			
		    #print cols[0],results[0],results[1]
		    #print  universe,(len(cols)-2),comparision, gene_set[cols[0]],results
	
	"""for items in gene:
	    if gene[items] > 0:
		print items,gene[items]
	f2.write("Gene")
	for j in enriched_gene_set:
	    f2.write("\t"+j)
	f2.write("\n")
	    
	
	for i in enriched_gene:
	    if enriched_gene[i] > 0:
		f2.write(i)
		for j in enriched_gene_set:
		    if i in (enriched_gene_set[j]):
			f2.write("\t1")
		    else:
			f2.write("\t0")
		f2.write("\n")"""
		
	f1.close()
	#f2.close()  
	
if __name__ == '__main__':
    main()
