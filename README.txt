”employed RNA-sequencing, ATAC–sequencing and ChIPmentation on well-defined thymocyte subsets that represent the continuum of human T cell”

Data Type:
Verboom: https://haematologica.org/article/view/8715SRA project: SRP132968SRA samples: SRR6727631, SRR6727635, SRR6727639, SRR6727666T-ALL samples, 4 male TAL subtypes
Roels: https://www.nature.com/articles/s41590-020-0747-9SRA project: SRP263016SRA samples:  SRR11832951, SRR11832952, SRR11832962, SRR11832963Normal single-positive CD4 and CD8 cells, two of each
Galaxy protocol:

•	Take your data too GalaxyEU https://usegalaxy.eu/ (Works on Chrome Safari Firefox and Microsoft Edge )
•	Register an account  
•	Upload the Data Type (asseccion numers in this case)
•	Paste data in empty dashboard
•	Name the data types after the source of the samples: Verboom for T-ALL samples and Roels for TCRαβ and TCRγδ T-cell development. Now data is saved. 
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