# RNAvis

## Description

This R/Shiny application has been created to share the results of the FLEPseq2 experiments we performed in our lab and with the scientific community. The version of the FLEPseq2 pipeline has been modified to improve performance (no fastq to fasta conversion), and contains an additional step that allows the detection of other untemolated nucleotide additions after the poly(A) tail, such as uridylation.



## Quick start

In order to use this application, follow these steps:

## 1. “Run selection” menu:

Choose a FLEPseq2 experiment and click on the “Get the data” Button. The corresponding files will be loaded and some general information about the run will appear on screen.

**Download Sample Data:** This tab allows the user to see the correspondence between the name of the file and the sample name that is biologically relevant. It is also possible to download the raw tables containing the FLEPseq2 results by sample. To do so, simply select on the dropdown menu the file that you want to download and click the “Download” button.

**Run Mapping Stats:** This tab will contain general information about the run. First, a table containing the total number of mapped reads and the number of mapped reads by chromosome. The following plots contain the percentage of mapped reads by chromosome and the number of detected genes for each sample.

**Poly(A) Distribution:** In this tab are displayed transcriptome-wide bulk and intergenic distributions of the poly(A) tail length. Those graphs are not generated as the user is using the application, because it is a computationally intensive step that can crash the server. Should you need other representations or to explore the tables yourself, the data is available for download.

## 2. “Transcript Specific” menu:

If you want to have more information about specific transcripts of one gene, just write the AGI in the box entitled “Select your transcript of interest” and the hit the “Get/update data”. 

**Transcript overview:** This tab will show you a schema of the gene with the introns and the exons. The table under this representation is the table extracted from the GTF annotation file and will give you the same informations.

**Intronic Profiles:** This plot will display a plot showing the transcripts with its different parts such as intron, exon, poly(A) tail and additional tail by a colored rectangular shape. An arrow indicates the direction of the transcript. When introns are spliced, the rectangular shape is not displayed, but they might be retained and will therefore be kept in the transcript. The user can choose the intronic profile, and the default loaded profile is the one with no introns retained (regular splicing).

**FLEPseq Results:** In this tab, the user can download results of the FLEPseq2 runs for only the gene of interest, in opposition to downloading the whole dataset. The user can also select a subset of columns to be displayed before generating and downloading the resulting table.

**Poly(A) Tail:** This tab shows the distribution of the poly(A) tail length of the transcripts for the selected gene. The user can either see the result as a smoothed distribution or as a histogram.

**Additional Tail:** in this section, we can see the transcript split into three categories: the ones that don’t have any additional tail, the ones that have one but not specifically considered as a U-tail, and the U-tails. The user can use the slider to modify the threshold of  the percentage of U residues that are needed to consider the transcript as being uridylated.

## 3. “Transcripts List” menu:

This section allows the user to select multiple genes either from a list (as in the previous section), but also to upload a list containing one AGI by line. The resulting plots will be the same as for one transcript, except that the intronic profiles cannot be displayed with several genes.


## Notes:

**Zooming**: Some plots of the application can be zoomed in. To do so, click and select an area of the plot and double-click on it. It is possible to zoom several times. A double-click with no selection will reset the plot to its original size.


## References 
Jia, J., Lu, W., Liu, B., Fang, H., Yu, Y., Mo, W., Zhang, H., Jin, X., Shu, Y., Long, Y., Pei, Y., Zhai, J., 2022. An atlas of plant full-length RNA reveals tissue-specific and monocots–dicots conserved regulation of poly(A) tail length. Nat. Plants 8, 1118–1126. https://doi.org/10.1038/s41477-022-01224-9

You can find the code for the pipeline on the following links:

- The original FLEPseq2 pipeline: https://github.com/ZhaiLab-SUSTech/FLEPSeq
- The version used here (with the additional tail detection): https://github.com/jackson-peter/FLEPseq2


