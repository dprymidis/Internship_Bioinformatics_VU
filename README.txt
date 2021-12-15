All Scripts run from inside the APC_Project folder!

The APC_Project folder contains:
	-Data folder:
		-Colorectal
			-Copy_number_variation
				-COAD_CNV.txt
				-READ_CNV.txt
				-CNV_metadata.json
				-CRC_CNV_data_prematched.txt (produced by Parse_CNV_data.R)
				-CRC_CNV_data_matched.txt (produced by Find_Sample_Match.R)
			-Gene_expression
				-Counts
				-Gene_expression_Counts_metadata.json
				-FPKM
				-FPKM-UQ
				-CRC_Gene_expression_all_Counts.txt  (produced by Parse_Gene_exp_data.R)
				-CRC_Gene_expression_Counts.txt, gene expression data for matched samples (produced by Find_Sample_Match.R)
			-Single_nucleotide_variation
				-COAD_SNV.maf
				-READ_SNV.maf

		-CRC_MSS_Samples.txt , all samples are taken into consideration (produced by TMB_and_MS_Status.R)
		-CRC_Gene_expression.txt, gene expression data for matched samples (produced by Find_Sample_Match.R)
		-CRC_Matched_Samples.txt, list of barcodes that match all data (produced by Find_Sample_Match.R) 
		-1_trunc_samples.prefiltered.txt, list of sample barcode and Class (produced by Class_preparation.R)
		-1_trunc_samples.filtered.txt, list of PCA filtered sample barcode and Class (produced by PCA.Rmd)
		-ElNet_genes.Class0_vs_Class1.txt, genes/features selected by the best model of Class0 vs Class 1 (produced by Classification.Rmd)
		-APC_reporter_gene_candidates.txt, list of target genes of Wnt pathway based on Schell et al. 2015
		

	-Scripts (By order of script execution):
		-TMB_and_MS_Status.R -> uses SNV data, calculates tumor mutational burden and colors by MS status, produces TMB plot and CRC_MSS_samples.txt
		-Parse_Gene_exp_data.R -> uses Gene expression data, produces CRC_Gene_expression_all_Counts.txt
		-Parse_CNV_data.R -> uses CNV data, produces CRC_CNV_data_prematched.txt
		-Find_Sample_Match.R -> uses SNV data, CRC_CNV_data_prematched.txt, CRC_Gene_expression_all_Counts.txt and CRC_MSS_samples.txt, produces CRC_Matched_Samples.txt, CRC_CNV_data_matched.txt and CRC_Gene_expression_Counts.txt
		-Class_preparation.R -> uses SNV data and CRC_MSS_Samples.txt, estimates kernel densities for APC truncating mutations and  investigates KRAs mutations, produces KDE plot and 1_trunc_samples.prefiltered.txt
		-PCA.Rmd -> uses 1_trunc_samples.prefiltered.txt and CRC_Gene_expression_all_Counts.txt, performs PCA, produces 1_trunc_samples.filtered.txt
		-Classification.Rmd -> uses 1_trunc_samples.filtered.txt and CRC_Gene_expression_all_Counts.txt, performs classification using 5fold CV repeated 5 times, produces ROC Plots and ElNet_genes.Class0_vs_Class1.txt
		-Wnt_genes.Rmd -> 1_trunc_samples.filtered.txt,  CRC_Gene_expression_all_Counts.txt, and APC_reporter_gene_candidates.txt, performs t-tests on Wnt genes per class combination, produces heatmap of Wnt target gene expression
		-EdgeR.GSEA.Rmd -> uses 1_trunc_samples.filtered.txt, CRC_Gene_expression_all_Counts.txt and ElNet_genes.Class0_vs_Class1.txt, performs Differential expression analysis, Overrepresentation analysis and Gene Set Enrichment Analysis, produces various plots as volcano plot of DE genes and folders containing the ORA and GSEA
		-Differential_expression_Bas_results_reproduction.pos.0v1.top.genes.Rmd -> uses 1_trunc_samples.filtered.txt and CRC_Gene_expression_all_Counts.txt to reproduce various plots of preliminary results
		-EdgeR_genes_int -> uses lists of genes from DE analysis, produces venn diagram of gene intersections
		-survival_analysis.R -> uses 1_trunc_samples.filtered.txt and data from Liu,J. et al. (2018), produces survival analysis plots

	-Other_files:
		-README.txt (this file)
		-GDCdata folder containing legacy data (The file contains downloaded indormation about cohorts, produced by TMB_and_MS_Status.R)
		

Script info: 


Data info:
The Single nucleotide variation data were extracted and moved in the Single_nucleotide_variation folder 
The Copy number variation data were extracted and moved in the Copy_number_variation folder
The Gene expression data were extracted and all the new compressed individual files were moved into the according folder in the Gene_expression folder

