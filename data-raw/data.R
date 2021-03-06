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
civic_biomarkers <- dplyr::select(civic_biomarkers, -c(pubmed_id,entrezgene,variant_groups,disease_url,disease_ontology_id,drug_interaction_type,variant_groups,api_url))
civic_biomarkers <- dplyr::rename(civic_biomarkers, biomarker_description = variant_description)

devtools::use_data(civic_biomarkers,overwrite=T)

cbmdb_biomarkers <- read.table(file="data-raw/cbmdb.biomarkers.tsv",sep="\t",na.strings = c("NA"),stringsAsFactors = F,header=T,comment.char="#",quote="")
cbmdb_biomarkers <- dplyr::select(cbmdb_biomarkers, -c(sources,variant_description))
cbmdb_biomarkers$mapping_category <- NA
cbmdb_biomarkers$mapping_rank <- NA
cbmdb_biomarkers$biomarker_description <- NA

devtools::use_data(cbmdb_biomarkers,overwrite=T)


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

pcgr_all_annotation_columns <- c('VARIANT_CLASS',
                                 'CHROM',
                                 'REF',
                                 'ALT',
                                 'POS',
                                 'CANCER_CENSUS_GERMLINE',
                                 'CANCER_CENSUS_SOMATIC',
                                 'CANCER_MUTATION_HOTSPOT',
                                 'INTOGEN_DRIVER',
                                 'INTOGEN_DRIVER_MUT',
                                 'CBMDB_ID',
                                 'CIVIC_ID',
                                 'CIVIC_ID_2',
                                 'CCDS',
                                 'ENTREZ_ID',
                                 'UNIPROT_ID',
                                 'ONCOGENE',
                                 'TUMOR_SUPPRESSOR',
                                 'CLINVAR',
                                 'CLINVAR_MSID',
                                 'CLINVAR_VARIANT_ORIGIN',
                                 'CLINVAR_SIG',
                                 'CLINVAR_TRAITS_ALL',
                                 'COSMIC',
                                 'COSMIC_SITE_HISTOLOGY',
                                 'COSMIC_DRUG_RESISTANCE',
                                 'COSMIC_MUTATION_ID',
                                 'OTHER_LITERATURE_DOCM',
                                 'OTHER_DISEASE_DOCM',
                                 'DBSNP',
                                 'DBSNPRSID',
                                 'DBSNPBUILDID',
                                 'DBSNP_VALIDATION',
                                 'GLOBAL_AF_EXAC',
                                 'GLOBAL_AF_1KG',
                                 'ANTINEOPLASTIC_DRUG_INTERACTION',
                                 'ANTINEOPLASTIC_DRUG_INTERACTIONS',
                                 'GENOMIC_CHANGE',
                                 'GENOME_VERSION',
                                 'VCF_SAMPLE_ID',
                                 'SYMBOL',
                                 'ONCOSCORE',
                                 'CONSEQUENCE',
                                 'PROTEIN_CHANGE',
                                 'PROTEIN_FEATURE',
                                 'GENE_NAME',
                                 'PROTEIN_DOMAIN',
                                 'EXON',
                                 'CDS_CHANGE',
                                 'EFFECT_PREDICTIONS',
                                 'VEP_ALL_CONSEQUENCE',
                                 'CALL_CONFIDENCE',
                                 'OXIDATION_ARTEFACT',
                                 'DP_TUMOR',
                                 'AF_TUMOR',
                                 'DP_NORMAL',
                                 'AF_NORMAL')
devtools::use_data(pcgr_all_annotation_columns, overwrite=T)


pcgr_tsv_tiered_columns <- c('GENOMIC_CHANGE',
                             'GENOME_VERSION',
                             'VCF_SAMPLE_ID',
                             'SYMBOL',
                             'GENE_NAME',
                             'CCDS',
                             'ENTREZ_ID',
                             'UNIPROT_ID',
                             'ONCOSCORE',
                             'ONCOGENE',
                             'TUMOR_SUPPRESSOR',
                             'INTOGEN_DRIVER',
                             'CANCER_CENSUS_SOMATIC',
                             'CANCER_CENSUS_GERMLINE',
                             'CONSEQUENCE',
                             'PROTEIN_CHANGE',
                             'PROTEIN_DOMAIN',
                             'CDS_CHANGE',
                             'EFFECT_PREDICTIONS',
                             'CANCER_MUTATION_HOTSPOT',
                             'INTOGEN_DRIVER_MUT',
                             'VEP_ALL_CONSEQUENCE',
                             'DBSNP',
                             'COSMIC',
                             'COSMIC_SITE_HISTOLOGY',
                             'COSMIC_DRUG_RESISTANCE',
                             'CLINVAR',
                             'CLINVAR_SIG',
                             'GLOBAL_AF_EXAC',
                             'GLOBAL_AF_1KG',
                             'CALL_CONFIDENCE',
                             'DP_TUMOR',
                             'AF_TUMOR',
                             'DP_NORMAL',
                             'AF_NORMAL',
                             'TIER',
                             'TIER_DESCRIPTION')

devtools::use_data(pcgr_tsv_tiered_columns, overwrite=T)

tier2_tags_display <- c("SYMBOL",
                        "PROTEIN_CHANGE",
                        "GENE_NAME",
                        "PROTEIN_DOMAIN",
                        "PROTEIN_FEATURE",
                        "CDS_CHANGE",
                        "CANCER_MUTATION_HOTSPOT",
                        "CANCER_CENSUS_SOMATIC",
                        "INTOGEN_DRIVER_MUT",
                        "CONSEQUENCE",
                        "EFFECT_PREDICTIONS",
                        "DBSNP",
                        "COSMIC",
                        "COSMIC_SITE_HISTOLOGY",
                        "CLINVAR",
                        "OTHER_LITERATURE_DOCM",
                        "OTHER_DISEASE_DOCM",
                        "ANTINEOPLASTIC_DRUG_INTERACTIONS",
                        "VEP_ALL_CONSEQUENCE",
                        "GENOME_VERSION",
                        "GENOMIC_CHANGE",
                        "ONCOGENE",
                        "TUMOR_SUPPRESSOR",
                        "GLOBAL_AF_EXAC",
                        "GLOBAL_AF_1KG",
                        "CLINVAR_TRAITS_ALL",
                        "CLINVAR_MSID",
                        "CALL_CONFIDENCE",
                        "DP_TUMOR",
                        "AF_TUMOR",
                        "DP_NORMAL",
                        "AF_NORMAL",
                        "OXIDATION_ARTEFACT")

devtools::use_data(tier2_tags_display, overwrite=T)

tier3_tags_display <- c("SYMBOL",
                        "ONCOSCORE",
                        "PROTEIN_CHANGE",
                        "GENE_NAME",
                        "PROTEIN_DOMAIN",
                        "PROTEIN_FEATURE",
                        "CDS_CHANGE",
                        "CANCER_MUTATION_HOTSPOT",
                        "CANCER_CENSUS_SOMATIC",
                        "INTOGEN_DRIVER_MUT",
                        "CONSEQUENCE",
                        "EFFECT_PREDICTIONS",
                        "DBSNP",
                        "COSMIC",
                        "COSMIC_SITE_HISTOLOGY",
                        "CLINVAR",
                        "ANTINEOPLASTIC_DRUG_INTERACTIONS",
                        "VEP_ALL_CONSEQUENCE",
                        "GENOME_VERSION",
                        "GENOMIC_CHANGE",
                        "ONCOGENE",
                        "TUMOR_SUPPRESSOR",
                        "CALL_CONFIDENCE",
                        "DP_TUMOR",
                        "AF_TUMOR",
                        "DP_NORMAL",
                        "AF_NORMAL",
                        "OXIDATION_ARTEFACT")

devtools::use_data(tier3_tags_display, overwrite=T)

tier4_tags_display <- c("SYMBOL",
                        "ONCOSCORE",
                        "PROTEIN_CHANGE",
                        "GENE_NAME",
                        "PROTEIN_DOMAIN",
                        "PROTEIN_FEATURE",
                        "CDS_CHANGE",
                        "CONSEQUENCE",
                        "EFFECT_PREDICTIONS",
                        "DBSNP",
                        "COSMIC",
                        "COSMIC_SITE_HISTOLOGY",
                        "CLINVAR",
                        "ANTINEOPLASTIC_DRUG_INTERACTIONS",
                        "GENOME_VERSION",
                        "VEP_ALL_CONSEQUENCE",
                        "GENOMIC_CHANGE",
                        "CALL_CONFIDENCE",
                        "DP_TUMOR",
                        "AF_TUMOR",
                        "DP_NORMAL",
                        "AF_NORMAL",
                        "OXIDATION_ARTEFACT")

devtools::use_data(tier4_tags_display, overwrite=T)

tier5_tags_display <- c("SYMBOL",
                        "ONCOSCORE",
                        "CONSEQUENCE",
                        "GENE_NAME",
                        "DBSNP",
                        "COSMIC",
                        "CDS_CHANGE",
                        "COSMIC_SITE_HISTOLOGY",
                        "CLINVAR",
                        "GENOME_VERSION",
                        "VEP_ALL_CONSEQUENCE",
                        "GENOMIC_CHANGE",
                        "CALL_CONFIDENCE",
                        "DP_TUMOR",
                        "AF_TUMOR",
                        "DP_NORMAL",
                        "AF_NORMAL",
                        "OXIDATION_ARTEFACT")

devtools::use_data(tier5_tags_display, overwrite=T)

tier1_tags_display <- c("SYMBOL",
                       "PROTEIN_CHANGE",
                       "CLINICAL_SIGNIFICANCE",
                       "EVIDENCE_LEVEL",
                       "EVIDENCE_TYPE",
                       "EVIDENCE_DIRECTION",
                       "DISEASE_NAME",
                       "DESCRIPTION",
                       "GENE_NAME",
                       "CITATION",
                       "DRUG_NAMES",
                       "RATING",
                       "PROTEIN_DOMAIN",
                       "PROTEIN_FEATURE",
                       "CDS_CHANGE",
                       "CANCER_MUTATION_HOTSPOT",
                       "OTHER_LITERATURE_DOCM",
                       "OTHER_DISEASE_DOCM",
                       "INTOGEN_DRIVER_MUT",
                       "CONSEQUENCE",
                       "EFFECT_PREDICTIONS",
                       "VEP_ALL_CONSEQUENCE",
                       "DBSNP",
                       "COSMIC",
                       "COSMIC_SITE_HISTOLOGY",
                       "CLINVAR",
                       "GENOME_VERSION",
                       "GENOMIC_CHANGE",
                       "CALL_CONFIDENCE",
                       "DP_TUMOR",
                       "AF_TUMOR",
                       "DP_NORMAL",
                       "AF_NORMAL")

devtools::use_data(tier1_tags_display, overwrite=T)
