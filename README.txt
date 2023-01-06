the data analysed in the project was too observe the differences for the development between the healthy t-cell and mutated t-cells(cancer) the data was taken for a puplic data library and both samples verboom and roels both samples was was also analysed using an illumina NGS for the sequencing data
”employed RNA-sequencing, ATAC–sequencing and ChIPmentation on well-defined thymocyte subsets that represent the continuum of human T cell”

Data Type:
Data was taken from the location in the form of txt format with the given accosiation number in them.

Verboom: https://haematologica.org/article/view/8715SRA project: SRP132968SRA samples: SRR6727631, SRR6727635, SRR6727639, SRR6727666T-ALL samples, 4 male TAL subtypes

Roels: https://www.nature.com/articles/s41590-020-0747-9SRA project: SRP263016SRA samples:  SRR11832951, SRR11832952, SRR11832962, SRR11832963Normal single-positive CD4 and CD8 cells, two of each.

Galaxy protocol:

•	Take your data too GalaxyEU https://usegalaxy.eu/ (Works on Chrome Safari Firefox and Microsoft Edge )
•	Register an account  
•	Upload the Data Type (asseccion numers in this case)
•	Paste data in empty dashboard
•	Name the data types after the source of the samples: Verboom for T-ALL samples and Roels for TCRαβ and TCRγδ T-cell development. Now data is saved. 
•	Start work with data through Get Data, and chose the data you have uploaded
	o	Load files in FastQ (Extract reads)
		Load files in FastQ: This step involves extracting reads from raw sequence data in the FastQ format.
	o	FastQC (Quality Controll)
		FastQC: This is a quality control tool that checks for any potential issues with the FastQ data.
	o	Fastp (Trimming)
		Fastp: This step involves trimming the FastQ data to remove low-quality or adapter sequences.
	o	Fasta/FastQ and Fastp (Trimming)
		Fasta/FastQ and Fastp: This step also involves trimming the data, but specifically for Fasta or FastQ formats.
	o	HISAT2 (Mapping to Genome)
		HISAT2: This tool is used to map the trimmed sequence data to a genome reference.
	o	HISAT2 log (Splice aware aligners)
		HISAT2 log: This step involves using splice-aware aligners to analyze the mapping results.
	o	MultiQC (Quality Controll)
		MultiQC: This is another quality control tool that provides an overview of the results from multiple analyses.
	o	RNA analysis (Quantification)
		RNA analysis: This step involves quantifying the RNA in the sequence data.
	o	edgeR RNA analysis (Differential Expression)
		edgeR RNA analysis: This tool is used to identify differential expression of genes in the RNA data.
	o	annotateMyIDs (Convert EntrezID to Genename)
		annotateMyIDs: This step converts EntrezID gene identifiers to gene names for easier interpretation of the results.

•	Download the annotatedID and edgeR files from Galaxy and start R project in Rstudio 

The code performs a two-sided t-test on each row of merged_GeneID, comparing the first four columns (belonging to the "Verboom" sample set) to the last four columns (belonging to the "Roels" sample set). The p-values for each t-test are stored in the Pvalue data frame.

The Pvalue data frame is filtered to only include rows with p-values less than 0.05. This data frame is called signPvalue.

The p-values in signPvalue are adjusted for multiple hypothesis testing using the Benjamini-Hochberg method. The resulting adjusted p-values are stored in a new data frame called AdjustedPvalue.

AdjustedPvalue is filtered to only include rows with adjusted p-values less than 0.05. The resulting data frame is called AdjustedPvalueOFClog.

The mean expression values for the "Verboom" and "Roels" sample sets are calculated for each row of topBottomGenes_tabl and stored in the VerboomMean and RoelsMean columns, respectively.

A scatterplot is created using topBottomGenes_tabl, with the VerboomMean and RoelsMean values plotted on the x-axis and y-axis, respectively. The

A volcano plot is created using AdjustedPvalueOFClog_1, with the FClog and -log10(AdjustedPvalue) values plotted on the x-axis and y-axis, respectively.

A horizontal line is added to the volcano plot at the -log10(0.05) value, and vertical lines are added at x-intercepts of -2 and 2.

A new column called "diffexpressed" is added to AdjustedPvalueOFClog_1 and filled with "NO" values. The values in the "diffexpressed" column are then changed to "UP" if the corresponding FClog value is greater than 2 and the adjusted p-value is less than 0.05, or to "DOWN" if the corresponding FClog value is less than -2 and the adjusted p-value is less than 0.05.

A new volcano plot is created using AdjustedPvalueOFClog_1, with the points colored based on their values in the "diffexpressed" column. A color scale is added to the plot, with the colors corresponding to the levels of "diffexpressed".

Overall, the code performs statistical tests to identify differentially expressed genes between the "Verboom" and "Roels" sample sets, and creates visualizations to display the results.
