the data analysed in the project was too observe the differences for the development between the healthy t-cell and mutated t-cells(cancer) the data was taken for a puplic data library and both samples verboom and roels both samples was was also analysed using an illumina NGS for the sequencing data
”employed RNA-sequencing, ATAC–sequencing and ChIPmentation on well-defined thymocyte subsets that represent the continuum of human T cell”

Data Type:
Data was taken from the location in the form of txt format with the given accosiation number in them.

Verboom: https://haematologica.org/article/view/8715SRA project: SRP132968SRA samples: SRR6727631, SRR6727635, SRR6727639, SRR6727666T-ALL samples, 4 male TAL subtypes

Roels: https://www.nature.com/articles/s41590-020-0747-9SRA project: SRP263016SRA samples:  SRR11832951, SRR11832952, SRR11832962, SRR11832963Normal single-positive CD4 and CD8 cells, two of each.

Galaxy protocol:
•	Take your data too GalaxyEU https://usegalaxy.eu/ (Works on Chrome Safari Firefox and Microsoft Edge )
•	Register an account  
•	Upload the Data Type (accosiation numbers in this case)
•	Paste data in empty dashboard
•	Name the data types after the source of the samples: Verboom for T-ALL samples and Roels for TCRαβ and TCRγδ T-cell	   development. Now data is saved. 
•	Start work with data through Get Data, and chose the data you have uploaded
	o	Load files in FastQ (Extract reads)
	o	FastQC (Quality Controll)
	o	Fastp (Trimming)
	o	Fasta/FastQ and Fastp (Trimming)
	o	HISAT2 (Mapping to Genome)
	o	HISAT2 log (Splice aware aligners)
	o	MultiQC (Quality Controll)
	o	RNA analysis (Quantification)
	o	edgeR RNA analysis (Differential Expression)
	o	annotateMyIDs (Convert EntrezID to Genename)

•	Download the annotatedID and edgeR files from Galaxy and start R project in Rstudio

The script contain a read in function for the data to the object. it contains a filtering method and the to remove the missing values and the data wrangling to move columns and rownames and re-order based on the mean of more than one column combined , ther is also a staatistical analysis where t.test was preformed on each row through a for-loop a filtering of genes based oon the pvalue of thee genes, optimal pvalue is below 0.05  there after there is some visualizations representing a heatmap, volcano plot and boxplot.
