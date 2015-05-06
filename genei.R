# Rscript to check gene symbols
# for previous symbols / synonymous
# extract approved symbols and entrez Ids

args <- commandArgs(TRUE)
print (args[1])

#genelist <- read.table("/Volumes/research/IIHG/Genetic Renal Platform/GRD_v2_genelist.txt", sep="\t")
genelist <- read.table(args[1], sep="\t")

allgenes <- paste("^", as.character(genelist[,1]), sep="", collapse="$|") 
hgnc <- read.delim(url("http://www.genenames.org/cgi-bin/download?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=gd_pub_eg_id&col=gd_pub_refseq_ids&status=Approved&status=Entry+Withdrawn&status_opt=2&where=%28%28gd_pub_chrom_map+not+like+%27%25patch%25%27+and+gd_pub_chrom_map+not+like+%27%25alternate+reference+locus%25%27%29+or+gd_pub_chrom_map+IS+NULL%29+and+gd_locus_type+%3D+%27gene+with+protein+product%27&order_by=gd_app_sym_sort&format=text&limit=&hgnc_dbtag=on&submit=submit"))

mappedgenes <- hgnc[c(grep(allgenes,hgnc[,c(2)]), grep(allgenes,hgnc[,c(6)]), grep(allgenes,hgnc[,c(5)])),c(1,2,5,6)]

write.table(unique(mappedgenes), paste(args[1], "mapped.txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE)