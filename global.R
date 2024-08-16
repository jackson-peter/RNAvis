# Global.R is loaded once at App start (not at browser refresh!)
# FLEP seq results organisation
mapping_dir="2_Mapping"
tail_dir="4_Tail"

tail_ext=".sorted.csv.gz"
polya_bulk_plot_ext="_polyA_bulk.png"
polya_ig_plot_ext="_polyA_ig.png"
bulk_ig_plot_ext="_polyA_bulk_ig.png"
cumul_polyA_plot_ext="_polyA_cumul.png"
mapping_ext=".bam.cov.txt"
mapping_cols=c("rname","startpos","endpos", "numreads", "covbases", "coverage", "meandepth","meanbaseq","meanmapq")

intronic_prfl_ftrs=c("exon", "intron", "polyA_tail", "add_tail")
intronic_prfl_clrs= carto_pal(length(intronic_prfl_ftrs), "Vivid")
names(intronic_prfl_clrs)= intronic_prfl_ftrs

index_ext=".sorted.csv.gz"
tabix_l_ext=".list.tsv"
sample_table="barcode_correspondance.*$"
transcript_col=4
start_col=6
end_col=7
DOC_DIR <- file.path("www","doc_markdowns")
