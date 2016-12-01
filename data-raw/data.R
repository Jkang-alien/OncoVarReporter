library(devtools)
library(dplyr)
library(stringr)
library(KEGGREST)

### PFAM

pfam_domains <- read.table(gzfile("data-raw/pfam.domains.tsv.gz"), sep="\t",quote="",stringsAsFactors=F, header=TRUE)
pfam_domains$pfam_id <- stringr::str_replace(pfam_domains$pfam_id,"\\.[0-9]{1,}$","")
devtools::use_data(pfam_domains,overwrite=T)

### UNIPROT/SWISSPROT PROTEIN FEATURES

swissprot_features <- read.table(gzfile("data-raw/uniprot.features.tsv.gz"),sep="\t",na.strings = c("NA"),stringsAsFactors = F,header=T,comment.char="",quote="")
devtools::use_data(swissprot_features,overwrite=T)

### CLINVAR

clinvar <- read.table(gzfile("data-raw/clinvar.tsv.gz"), header=T, stringsAsFactors = F, quote="", sep="\t",na.strings=c("","-",NA),comment.char="")
clinvar$var_id <- paste(paste0('chr',clinvar$chrom),clinvar$pos,clinvar$ref,clinvar$alt,sep="_")
clinvar$build <- 'GRCh37'
clinvar$CLINVAR_TRAITS_ALL <- paste(clinvar$variant_origin,clinvar$trait,sep=" - ")
clinvar$CLINVAR_MSID <- as.character(clinvar$measureset_id)
devtools::use_data(clinvar,overwrite=T)

civic_biomarkers <- read.table(file="data-raw/civic.biomarkers.tsv",sep="\t",na.strings = c("NA"),stringsAsFactors = F,header=T,comment.char="#",quote="")
devtools::use_data(civic_biomarkers,overwrite=T)

### GENE/TRANSCRIPT ANNOTATIONS AND X-REFS

gene_xref <- read.table(gzfile("data-raw/gene.transcript.onco_xref.GRCh37.tsv.gz"),header=T,stringsAsFactors = F,quote="",sep="\t",na.strings=c("","NA"),comment.char="#")
devtools::use_data(gene_xref, overwrite = T)

### MUTATIONAL SIGNATURES

signatures_aetiologies <- read.table(file="data-raw/signatures_aetiologies.tsv",sep="\t",quote="",header=T,na.strings=c("NA"),stringsAsFactors = F)
devtools::use_data(signatures_aetiologies, overwrite = T)

### DOCM

docm_processed <- read.table(file="data-raw/docm_processed.tsv",sep="\t",quote="",header=T,na.strings=c("NA"),stringsAsFactors = F)
docm_processed <- dplyr::select(docm_processed, docm_id,disease,pmid,citation)
docm_processed <- dplyr::rename(docm_processed, OTHER_DISEASE_DOCM = disease, VAR_ID = docm_id)

tmp1 <- dplyr::select(docm_processed,VAR_ID,pmid)
tmp2 <- tidyr::separate_rows(tmp1,pmid,sep=",")

tmp3 <- dplyr::select(docm_processed,VAR_ID,citation)
tmp4 <- tidyr::separate_rows(tmp3,citation,sep="\\|")

tmp8 <- data.frame('citation'=tmp4$citation,'pmid'=tmp2$pmid, stringsAsFactors = F) %>% dplyr::distinct()
tmp8$html_citation <- paste0('<a href=\'http://www.ncbi.nlm.nih.gov/pubmed/',tmp8$pmid,'\' target=\'_blank\'>',tmp8$citation,'</a>')

tmp9 <- dplyr::left_join(tmp2,tmp8)

docm_literature <- as.data.frame(dplyr::select(tmp9,VAR_ID,html_citation) %>% dplyr::group_by(VAR_ID) %>% dplyr::summarise(OTHER_LITERATURE_DOCM = paste(html_citation, collapse=", ")))
docm_literature <- dplyr::left_join(docm_literature, dplyr::select(docm_processed,VAR_ID,OTHER_DISEASE_DOCM))
docm_literature$OTHER_DISEASE_DOCM <- stringr::str_replace_all(docm_literature$OTHER_DISEASE_DOCM,",",", ")

devtools::use_data(docm_literature, overwrite = T)

### KEGG

kegg_section_subsection <- read.table('data-raw/kegg.section.subsection.tsv',sep='\t',header=T,stringsAsFactors=F)

kegg_pathway <- as.data.frame(keggList('pathway','hsa'),stringsAsFactors=FALSE)
kegg_pathway$name <- str_replace(kegg_pathway[,1],' - Homo sapiens \\(human\\)','')
kegg_pathway$pathway_id <- str_replace(rownames(kegg_pathway),'path:','')
kegg_pathway <- dplyr::select(kegg_pathway, pathway_id, name)
kegg_pathway <- merge(kegg_pathway, kegg_section_subsection,by.x="pathway_id",by.y="pathway_id",all.x=TRUE)
row.names(kegg_pathway) <- NULL

kegg_disease <- as.data.frame(keggList('disease'),stringsAsFactors=FALSE)
kegg_disease$name <- str_replace(kegg_disease[,1],' - Homo sapiens (human)','')
kegg_disease$disease_id <- str_replace(rownames(kegg_disease),'ds:','')
kegg_disease <- dplyr::select(kegg_disease, disease_id, name)
row.names(kegg_disease) <- NULL

i <- 1
kegg_pathway_gene <- data.frame('pathway_id'=character(),'gene_id'=integer())
kegg_pathway_disease <- data.frame('pathway_id'=character(),'disease_id'=character())
while(i <= nrow(kegg_pathway)){
  pathway_id <- kegg_pathway[i,]$pathway_id
  kegg_ref <- keggGet(pathway_id)
  if (!is.null(kegg_ref[[1]]$DISEASE)){
    disease_ids <- names(kegg_ref[[1]]$DISEASE)
    kegg_pathway_disease <- rbind(kegg_pathway_disease, data.frame('pathway_id'=rep(pathway_id,length(names(kegg_ref[[1]]$DISEASE))),'disease_id'=names(kegg_ref[[1]]$DISEASE)))
  }
  if (!is.null(kegg_ref[[1]]$GENE)){
    genelist <- kegg_ref[[1]]$GENE
    odd_indices <- seq(1,length(genelist),2)
    gene_indices <- as.integer(genelist[odd_indices])
    kegg_pathway_gene <- rbind(kegg_pathway_gene, data.frame('pathway_id'=rep(pathway_id,length(gene_indices)),'gene_id'=gene_indices,stringsAsFactors=FALSE))
  }
  cat("Pathway number ",i," - ",kegg_pathway[i,]$name,"\n")

  i <- i + 1
}

devtools::use_data(kegg_pathway, overwrite=T)
devtools::use_data(kegg_disease, overwrite=T)
devtools::use_data(kegg_pathway_gene, overwrite=T)
devtools::use_data(kegg_pathway_disease, overwrite=T)

