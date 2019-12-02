#!/usr/bin/env pyhton
import os, sys, argparse, re
import numpy as np
import scipy
import pandas as pd
import scipy.stats

def parse_args():
	parser = argparse.ArgumentParser(description='Calculation of signature score for samples') 
	parser.add_argument('-f1','--file_1', type=str, required=True, 
		help='Gene list for comparsion or the gene signature')
	parser.add_argument('-N','--un_genes', type=int, required=False, default=45956,
                    help='No of genes in universe, Default value by GSEA = 45956')
	parser.add_argument('-f2','--file_2', type=str, required=False, 
		help='Gene list for universe')
	parser.add_argument('-f3','--geneset', type=str, required=True, 
		help='Gene sets in gmt format')
	parser.add_argument('-m','--method', type=str, required=True, 
		help='Method to use for enrichment: fisher or chi2 or hypergeom')
	parser.add_argument('-p','--pvalue', type=float, required=False, default=1,
                    help='P-value for enrichment, Default value = 0.05')
	parser.add_argument('-o1','--output1', type=str, required=True, 
		help='Output file with gene enrichment values')
	
	args = parser.parse_args()
	return args

def universe(total,file_1):
	all_ant = {}
	if file_1:	
		with open(file_1,'r') as f:
		    while True:
			text = f.readline()
			if text == "":
			    break
			cols = text.rstrip().split('\t')
			all_ant[cols[0]] = 0
			#enriched_gene[cols[0]] = 0
		total_genes = len(all_ant)
	else:
		total_genes = total
		all_ant['0'] = 0
		
	return (total_genes,all_ant)
def comparision(file_1):
	gene = {}
	enriched_gene = {}
	with open(file_1,'r') as f:
	    while True:
		text = f.readline()
		if text == "":
		    break
		cols = text.rstrip().split('\t')
		gene[cols[0]] = 0
		enriched_gene[cols[0]] = 0
		
	return (len(gene),gene,enriched_gene)	

def main():
	
	args = parse_args()	
	gene = {}
	gene_set = {}
	enriched_gene = {}
	enriched_gene_set = {}
	gene_set_no = {}
	all_ant = {}
	
	genes_in_universe,all_ant = universe(args.un_genes,args.file_2)
	genes_in_comparision,gene,enriched_gene = comparision(args.file_1)
	
	f1 = open(args.output1,"w")
	f2 = open(args.output1[:-4]+"_overlap_genes.txt","w")
		
	f1.write("Gene_set\tP-value\t#genes in universe\t#Genes in Gene Set\t#genes in comparison\t#Genes in Overlap\tOverlapped genes\n")
	
	with open(args.geneset,'r') as f:
	    while True:
		text = f.readline()
		if text == "":
		    break
		cols = text.rstrip().split('\t')
		gene_set[cols[0]] = {}
		
		gene_set_no[cols[0]] = 0
		for i in range(2,len(cols)):
		    if cols[i] in gene:
		    	    #print gene
		    	    (gene_set[cols[0]])[cols[i]] = cols[i]
		    	    gene_set_no[cols[0]] = gene_set_no[cols[0]] + 1
		    	    gene[cols[i]] = gene[cols[i]] + 1
			
		if args.file_2:	
			ant = 0        
			for i in range(2,len(cols)):
			    if cols[i] in all_ant:
				ant = ant + 1
		else:
			ant = (len(cols)-2)
			     
		if gene_set_no[cols[0]] > 0:
		    list1 = [genes_in_universe,ant]
		    list2 = [genes_in_comparision, gene_set_no[cols[0]]]
		    array = np.array([list1, list2], dtype = np.int32)
		    
		    if args.method == 'fisher':
		    	    results = (scipy.stats.chi2_contingency(array))[1]
		    elif args.method == 'chi2':
		    	    results = (scipy.stats.fisher_exact(array))[1]
		    elif args.method == 'hypergeom':
		    	    results = scipy.stats.hypergeom.sf((gene_set_no[cols[0]]), genes_in_universe, ant, genes_in_comparision, 1)	 
		    else:
		    	    "method is incorrect"
		    
		    if results <= args.pvalue:
						
			enriched_gene_set[cols[0]] = {}
			for i in range(2,len(cols)):
			    if cols[i] in gene:
				(enriched_gene_set[cols[0]])[cols[i]] = cols[i]
				enriched_gene[cols[i]] = enriched_gene[cols[i]] + 1
			
			f1.write(cols[0]+"\t"+str(results)+"\t"+str(genes_in_universe)+"\t"+str(ant)+"\t"+str(genes_in_comparision)+"\t"+str(gene_set_no[cols[0]])+"\t")
			for item in enriched_gene_set[cols[0]]:
				f1.write(item+",")
			f1.write("\n")
		
	for items in gene:
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
		f2.write("\n")
		
	f1.close()
	f2.close()  
		
	
if __name__ == '__main__':
	try:
		main()
	except:
		print "An unknown error occurred.\n"
		raise

