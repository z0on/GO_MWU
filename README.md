#	Rank-based Gene Ontology Analysis with Adaptive Clustering


Mikhail V. Matz
UT Austin, 
matz@utexas.edu

Windows users: make sure perl is installed and is in your system's PATH 
(or you can specify the path to perl executable explicitly).

![alt tag](https://github.com/z0on/GO_MWU/blob/master/heats_cc_gomwu.png)

Short Guide 
-----------

1. Put all this into the same directory:
	- scripts: GO_MWU.R, gomwu_a.pl, gomwu_b.pl, gomwu.functions.R
	- GO hierarchy file 
		(go.obo, http://www.geneontology.org/GO.downloads.ontology.shtml)
	- table of GO annotations for your sequences: two-column (gene id - GO terms), 
		tab-delimited, one line per gene, multiple GO terms separated by semicolon. 
		If you have multiple lines per gene, use nrify_GOtable.pl to merge them.
	- table of measure of interest for your sequences: two columns of comma-separated 
		values: gene id, continuous measure of significance such as log(fold-change) or
		-log(p-value). To perform standard GO enrichment analysis based on Fisher's 
		exact test, use binary measure (1 or 0, i.e., either sgnificant or not). To analyze modules derived from WGCNA, specify 0 for genes not included in the module and the kME value (number between 0 and 1, module membership score) for genes included in the module.
	
	It is important to have the latter two tables representing the whole 
	genome (or transcriptome) - at least the portion that was measured -
	rather than some select group of genes since the test relies on comparing
	the behavior of individual GO categories to the whole.

2. Make sure you have perl and R. Windows-based people: install perl from [here](http://strawberryperl.com/) and specify your perl.exe file with full path to it as the perlPath argument to function gomwuStats (inside the GO_MWU.R script). The R part requires package "ape", which 
you might need to install prior to running this method.

3. Open GO_MWU.R script; edit the input file names, mark and execute bits of code
separated by blank lines one by one. Follow instructions given as comments in the script.

4. Drag corner of the plot to rescale and match text and tree. After this, to
achieve better "word map" effect, rerun gomwuPlot with modified "txtsize" parameter.

5. Save the plot as pdf file.
   
###	The Output:

The plot consists of three parts: 

-	Hierarchical clustering tree of significant GO categories based on shared genes 
	in the current dataset. Categories with no branch length between them are subsets 
	of each other and their significance is most likely driven by the same genes.

-	Category names, plotted in different colors and fonts. Fonts indicate the level of 
	statistical significance, colors indicate enrichment of GO categories with either 
	up- (red) or down- (blue) regulated genes. The category names are preceded by the 
	fraction indicating the number of "good candidates" relative to the total number of 
	genes belonging to this category. The "good candidates" are the genes exceeding an 
	arbitrary 'absValue' cutoff in their significance measure. Adjust 'absValue' parameter
	according to what your measure is. By default it is set to -log(0.05,10), assuming 
	that the measure is a signed log p-value (so, the "good candidates" would be the ones
	with raw p-value < 0.05). Ideally we would like to see more than one such gene per 
	displayed GO category. With 'level1=1' the script will display all the categories 
	containing "good candidates", which is a good way to summarize the whole GO content 
	of the experiment. Note that 'absValue' parameter does not affect statistics and 
	serves just the illustrative purpose. In the Fisher-test mode (binary significance
	measure) and signed WGCNA module analysis the colors are not used; in that case specify 
	absValue=0.001 to make the script display the fraction of genes with non-zero measure
	within a GO category. 
	
-	The legend giving the correspondence of the fonts to significance thresholds. The 
	method corrects the p-values using Benjamini-Hochberg false discovery rate procedure except when analyzing WGCNA modules; in that case the false discovery rate is determined from ten permutations where significance measures are randomly shuffled among genes. 
	To set different thresholds for plotting, change parameters 'level1', 'level2' and 
	'level3' in gomwuPlot.

In addition, the script prints out the number of GO categories displayed and the fraction
of "good candidates" that these categories account for. This is useful to evaluate whether
the generated GO summary really accounts for a substantial portion of what was going on.


Suggested citation
------------------
In its present form GO_MWU method was first used in:
Wright, R. M., Aglyamova, G. V., Meyer, E.  and Matz, M. V. Gene expression associated with white syndromes in a reef-building coral, Acropora hyacinthus. BMC Genomics 2015, 16: 371. 
( http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1540-2 )


How it works
------------

In contrast to most other "GO enrichment analysis" methods (e.g., GeneMerge or DAVID), this one does not look for GO categories enriched 
among "significant" genes. 

Instead, it measures whether each GO category is significantly enriched by either up or 
down-regulated genes. Basically, the method tests whether the genes belonging 
to a certain GO category are significantly bunched up near the top or the bottom 
of the global ranked list of genes, instead of being spread evenly all over it. The
test used is called the Mann-Whitney U (MWU) test.

The major advantage of this approach is that the experimenter does not have to 
impose an arbitrary threshold for initial selection of "significant genes", and thus 
the whole dataset can be used to gain information. 

In fact, no preliminary statistical test is required prior to the analysis; 
the method is best suited to analyze the distribution of raw measures, 
such as dN/dS values, log-fold-changes of gene expression, raw p-values, 
or kME (correlation) values from WGCNA.

The method can also be run in a traditional mode, looking for GO categories
significantly over-represented among "significant genes" (based on Fisher's exact test).
To make the method work in this mode, the measure of significance should be binary
(1 or 0, i.e., significant or not).

The method also considerably more powerful than standard Fisher's exact test for analysis of WGCNA modules, since it derives additional power from the module membership scores (kME values).

The method automatically retrieves all the missing parental terms for the lower-level 
GO categories. Then, fully redundant GO categories (i.e., containing exactly the same 
genes) are collapsed under name of the lower-level (more specific) term. Then,
highly similar categories are merged according to complete linkage clustering based on
the fraction of shared genes. The distance measure for clustering, introduced in 
Kosiol et al 2008, is the number of genes shared among the two GO categories within the analyzed dataset divided by the size of the smaller of the two categories. The resulting hierarchical tree is then “cut” at the adjustable 
“height” ('cutTreeHeight' parameter in the call to gomwuStats) to merge clustered 
categories. The default for cutTreeHeight is 0.25, implying that a group of categories 
will be merged if the most dissimilar two of them share >75% of genes included in the 
smaller of the two. The merged categories inherit the name of the largest one. This 
simplifies the GO hierarchy, generates biologically meaningful groups of categories 
tailored for the particular dataset, and improves the multiple testing situation. 

In the final plot, the method shows hierarchical clustering of GO categories based on 
the number of genes shared between them, to indicate which categories might be 
significant because of the same genes. 

Where does it come from
-----------------------

The MWU-based method of GO analysis was first introduced in Nielsen et al PLoS Biol 2005, 
3:e170. Its was used together with the hierarchical clustering of displayed GO categories 
in Kosiol et al PLoS Genet 2008, 4:e1000144 and Voolstra et al PLoS ONE 2011, 
6(5): e20392. A related rank-based method of GO analysis is GSEA: doi:pnas.0506580102.

Details on the input format
---------------------------

The GO annotations table should have two tab-delimited columns: gene name, and a string 
of concatenated GO terms separates by semicolons, like this: 
```
isogroup0	GO:0016301;GO:0005515;GO:0007507;GO:0030239;GO:0065007
isogroup10	GO:0044424
isogroup100	GO:0006810;GO:0080090;GO:0023033;GO:0065008;GO:0044237;GO:0051649
isogroup10001	GO:0009987;GO:0000323;GO:0016787
isogroup10002	GO:0005488
isogroup10004	unknown
....
```
(the genes without annotation should be called "unknown", if you want to analyze these too)

NB: The table must contain just a single line per gene. If you have a table in which
a gene appears on several lines, use the included perl script nrify_GOtable.pl 
to compile a non-redundant table. 
nrify_GOtable.pl takes a single argument, which is the filename 
to process, and prints the result to STDOUT, so use it like this:
nrify_GOtable.pl [filename of the redundant table] > [filename of the non-redundant table]

The GO table provided with the scripts is called amil_defog_iso2go.tab

The table of significance measures: it is the comma-separated (CSV) table of continuous
measures that must be associated with GO enrichment (for example, kME value, p-value, 
or dN/dS value). The table should have a header line, but what is in it is irrelevant. 
The first column should contain the  gene name, and the second - the measure of interest:

```sh
gene,logP
isogroup0,8.3
isogroup1,2
isogroup10,9.9
isogroup100,2.6
isogroup1000,2.4
isogroup10000,1
isogroup10001,0
isogroup10002,0.9
isogroup10003,-0.3
isogroup10004,1.5
isogroup10006,-6.9
...
```
The test file provided here is called "heats.csv". It contains the results of
gene expression profiling of coral response to long-term heat stress, in the 
form of "signed negative log p-values". These measures are negative decimal
logarithms of the raw (uncorrected) p-value for each gene, multiplied
by -1 if the gene was down-regulated. As a result, highly significant 
up-regulated genes get highly positive values, and highly significant
down-regulated genes get highly negative values.

Analyzing WGCNA modules
------------ 

In this case the method does two layers of testing: first, global Fisher's exact test for presence-absence of functional categories in the module, and second, within-module MWU test for association of the included functional categories with higher kME values (module membership scores). The product of the two p-values becomes the new test statistic. The false discovery rate is then determined from ten permutations where significance measures are randomly shuffled among genes.

To perform this analysis, the input data file should list all genes that were used in WGCNA; the genes that are not included in the module should receive a significance measure of 0, and genes within the module - a kME value. If you are analysing signed WGCNA modules, add 'Module=TRUE,Alternative="g"' option when running gomwuStats(); for unsigned modules, just 'Module=TRUE' (see comments in the GO_MWU.R script).
 
Output Files
------------ 

The script generates three tables:

(GO division)_(input filename) : main data table containing reformatted and augmented 
GO terms for each gene (in addition to the originally listed terms, the script finds 
all their parental terms if any were missing), and measures of interest. 

dissim_(GO division)_(go-to-gene table filename) : dissimilarity matrix of GO categories 
based on the number of genes shared between them in the dataset. 

MWU_(GO division)_(input filename) : the results of MWU test
