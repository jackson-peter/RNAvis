HTSLIB_PATH = "/shared/biotools/htslib/1.18/bin/"
README="/shared/home/jpeter/Scripts/RNAvis/README.md"
datasets_path <- "/shared/home/jpeter/ShinyApps/RNAvis/data/FLEPseq_runs_new/"

runs_infos <- read.table("/shared/home/jpeter/ShinyApps/RNAvis/data/Runs_infos.tsv", col.names =c("Run_name", "Run_infos"), sep = '\t', fill=T) 

flep_runs_df <- get_flepruns(datasets_path)
flep_runs <- setNames(flep_runs_df$flep_run, flep_runs_df$run_desc)

ref_gtf <- "/shared/home/jpeter/Data/ReferenceGenomes/A_thaliana/Araport11/Araport11_GTF_genes_transposons.current.gtf"
gtf_colnames <- c("seqnames", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")