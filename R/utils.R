library(BSgenome.Hsapiens.UCSC.hg19)
library(magrittr)
library(deconstructSigs)
library(data.table)
library(ggplot2)



#' Function that plots call numbers
#'
#' @param project_directory name of roject directory
#' @param strelka_snv_vcf_gz name of bgzipped VCF file with annotated Strelka SNVs
#' @param mutect_snv_vcf_gz name of bgzipped VCF file with annotated MuTect SNVs
#' @param strelka_indel_vcf_gz name of bgzipped VCF file with annotated Strelka indels
#' @param project_sample_stats sample_stats data frame
#' @param signatures.limit Number of signatures to limit mutational signature analysis
#' @param strong_calls_only Logical indicating if report is based on strong calls only or all calls
#'
#' @return p
#' @export
#'
generate_project_reports <- function(project_directory, strelka_snv_vcf_gz, mutect_snv_vcf_gz, strelka_indel_vcf_gz = NULL, project_sample_stats = NULL, signatures.limit = 10, strong_calls_only = FALSE, project_name = 'NCGC_tumor_exome'){
  vcf_data_df_mutect <- OncoVarReporter::get_calls(mutect_snv_vcf_gz)
  vcf_data_df_strelka <- OncoVarReporter::get_calls(strelka_snv_vcf_gz)
  vcf_data_df_strelka_indel <- OncoVarReporter::get_calls(strelka_indel_vcf_gz)

  vcf_data_df_mutect <- dplyr::select(vcf_data_df_mutect, -c(end,width,strand,totalDepth,refDepth,altDepth,sampleNames))
  vcf_data_df_strelka <- dplyr::select(vcf_data_df_strelka, -c(end,width,strand,totalDepth,refDepth,altDepth,sampleNames))
  vcf_data_df_strelka_indel <- dplyr::select(vcf_data_df_strelka_indel, -c(end,width,strand,totalDepth,refDepth,altDepth,sampleNames))

  common_columns <- intersect(intersect(colnames(vcf_data_df_mutect), colnames(vcf_data_df_strelka)),colnames(vcf_data_df_strelka_indel))
  vcf_data_df_mutect <- dplyr::select(vcf_data_df_mutect, dplyr::one_of(common_columns))
  vcf_data_df_strelka <- dplyr::select(vcf_data_df_strelka, dplyr::one_of(common_columns))
  vcf_data_df_strelka_indel <- dplyr::select(vcf_data_df_strelka_indel,dplyr::one_of(common_columns))

  all_tier_tsv_variants <- NULL
  all_tier_tsv_fname <- paste0(project_name,'.snvs_indels.weak_strong.tiers.tsv')
  if(strong_calls_only == TRUE){
    all_tier_tsv_fname <- paste0(project_name,'.snvs_indels.strong.tiers.tsv')
  }

  for(sample_id in unique(vcf_data_df_strelka$VCF_SAMPLE_ID)){
    cat(sample_id,"\n")
    tumor_id <- unlist(stringr::str_split(sample_id,"_",2))[2]
    control_id <- unlist(stringr::str_split(sample_id,"_",2))[1]

    sample_coverage_statistics <- NULL
    show_coverage_stats <- FALSE
    if(!('T_ID' %in% colnames(project_sample_stats) & 'C_ID' %in% colnames(project_sample_stats))){
      cat('WARNING: Did not find \'T_ID\' and \'C_ID\' in sample_stats',sep="\n")
    }
    else{
      sample_coverage_statistics <- project_sample_stats %>% dplyr::filter(T_ID == tumor_id & C_ID == control_id)
      if(nrow(sample_coverage_statistics) > 0 & 'T_Mean_coverage' %in% colnames(sample_coverage_statistics) & 'T_Median_coverage' %in% colnames(sample_coverage_statistics) & 'C_Mean_coverage' %in% colnames(sample_coverage_statistics) & 'C_Median_coverage' %in% colnames(sample_coverage_statistics)){
        show_coverage_stats <- TRUE
      }
    }

    sample_calls_mutect <- dplyr::filter(vcf_data_df_mutect, VCF_SAMPLE_ID == sample_id)
    sample_calls_strelka <- dplyr::filter(vcf_data_df_strelka, VCF_SAMPLE_ID == sample_id)
    sample_calls_strelka_indel <- dplyr::filter(vcf_data_df_strelka_indel, VCF_SAMPLE_ID == sample_id)
    if(nrow(sample_calls_strelka_indel) > 0){
      sample_calls_strelka_indel$CALL_CONFIDENCE <- 'STRONG'
    }

    var_id_STRONG <- dplyr::inner_join(dplyr::select(sample_calls_mutect,VAR_ID),dplyr::select(sample_calls_strelka,VAR_ID),by=c("VAR_ID"))
    var_id_WEAK_MUTECT <- dplyr::anti_join(dplyr::select(sample_calls_mutect,VAR_ID),var_id_STRONG,by=c("VAR_ID"))
    var_id_WEAK_STRELKA <- dplyr::anti_join(dplyr::select(sample_calls_strelka,VAR_ID),var_id_STRONG,by=c("VAR_ID"))

    strong_calls <- dplyr::inner_join(sample_calls_mutect, var_id_STRONG,by=c("VAR_ID"))
    if(nrow(strong_calls) > 0){
      strong_calls$CALL_CONFIDENCE <- 'STRONG'
    }
    weak_mutect_calls <- dplyr::anti_join(sample_calls_mutect, var_id_STRONG,by=c("VAR_ID"))
    if(nrow(weak_mutect_calls) > 0){
      weak_mutect_calls$CALL_CONFIDENCE <- 'WEAK_MUTECT'
    }
    weak_strelka_calls <- dplyr::anti_join(sample_calls_strelka, var_id_STRONG,by=c("VAR_ID"))
    if(nrow(weak_strelka_calls) > 0){
      weak_strelka_calls$CALL_CONFIDENCE <- 'WEAK_STRELKA'
    }

    sample_calls <- rbind(strong_calls,weak_strelka_calls,weak_mutect_calls,sample_calls_strelka_indel)
    if(strong_calls_only == TRUE){
      sample_calls <- rbind(strong_calls, sample_calls_strelka_indel)
    }

    sample_calls$OXIDATION_ARTEFACT <- NA

    report_data <- OncoVarReporter::generate_report_data(sample_calls,
                                                         sample_id = sample_id,
                                                         minimum_n_signature_analysis = 50,
                                                         signatures.limit = signatures.limit)
    if(!is.null(report_data$tsv_variants)){
      all_tier_tsv_variants <- rbind(all_tier_tsv_variants, report_data$tsv_variants)
    }
    report_data$sequencing_approach <- sample_coverage_statistics$sequencing_approach
    report_data$sample_id <- sample_id

    cnv_plot <- FALSE
    cnv_report_tsgene_loss <- FALSE
    cnv_report_oncogene_gain <- FALSE
    cnv_report_biomarkers <- FALSE
    cnv_report_segments <- FALSE

    cnv_plot_facets_png <- paste0(project_directory,'/cnv/facets/',sample_id,'_plots.png')
    if(file.exists(cnv_plot_facets_png)){
      cnv_plot <- TRUE
    }
    cnv_segments_file <- paste0(project_directory,'/cnv/facets/',sample_id,'_segment_values.tsv')
    sample_cnv_tsv <- paste0(project_directory,'/',sample_id,'.facets.cnv.annotated.tsv')
    if(file.exists(cnv_segments_file)){
      cnv_data <- OncoVarReporter::cnv_segment_annotation(cnv_segments_file, format='facets')
      if(nrow(cnv_data$ranked_segments) > 0){
        cnv_report_segments <- TRUE
      }
      if(nrow(cnv_data$tsgene_homozygous_deletion) > 0){
        cnv_report_tsgene_loss <- TRUE
      }
      if(nrow(cnv_data$oncogene_amplified) > 0){
        cnv_report_oncogene_gain <- TRUE
      }
      if(nrow(cnv_data$cnv_df_for_print) > 0){
        write.table(cnv_data$cnv_df_for_print,file=sample_cnv_tsv,col.names = T,row.names = F,quote=F,sep="\t")
        gzip_command <- paste0('gzip -f ',sample_cnv_tsv)
        system(gzip_command, intern=F)
      }
      if(!is.null(cnv_data$cna_biomarkers)){
        if(nrow(cnv_data$cna_biomarkers) > 0){
          cnv_report_biomarkers <- TRUE
        }
      }
    }

    tier1_report <- FALSE
    tier2_report <- FALSE
    tier3_report <- FALSE
    tier4_report <- FALSE
    tier5_report <- FALSE
    signature_report <- FALSE
    missing_signature_data <- FALSE

    tier1_report <- report_data$tier1_report
    tier2_report <- report_data$tier2_report
    tier3_report <- report_data$tier3_report
    tier4_report <- report_data$tier4_report
    tier5_report <- report_data$tier5_report
    signature_report <- report_data$signature_report
    if(signature_report == FALSE){
      missing_signature_data <- TRUE
    }
    show_data_sources <- TRUE

    if(strong_calls_only == TRUE){
      rmarkdown::render(system.file("templates","report_strong.Rmd", package="OncoVarReporter"), output_file = paste0(sample_id,'.strong.tumor_report.html'), output_dir = project_directory, params = list(signature_report = signature_report, tier1_report = tier1_report, tier2_report = tier2_report, tier3_report = tier3_report, tier4_report = tier4_report, tier5_report = tier5_report, cnv_report_tsgene_loss = cnv_report_tsgene_loss, cnv_report_oncogene_gain = cnv_report_oncogene_gain, cnv_plot = cnv_plot, cnv_report_segments = cnv_report_segments, cnv_report_biomarkers = cnv_report_biomarkers, show_coverage_stats = show_coverage_stats, show_data_sources = show_data_sources),quiet=T)
    }
    else{
      rmarkdown::render(system.file("templates","report.Rmd", package="OncoVarReporter"), output_file = paste0(sample_id,'.strong_weak.tumor_report.html'), output_dir = project_directory, params = list(signature_report = signature_report, tier1_report = tier1_report, tier2_report = tier2_report, tier3_report = tier3_report, tier4_report = tier4_report, tier5_report = tier5_report, cnv_report_tsgene_loss = cnv_report_tsgene_loss, cnv_report_oncogene_gain = cnv_report_oncogene_gain, cnv_plot = cnv_plot, cnv_report_segments = cnv_report_segments, cnv_report_biomarkers = cnv_report_biomarkers, show_coverage_stats = show_coverage_stats, show_data_sources = show_data_sources),quiet=T)
    }
  }
  if(!is.null(all_tier_tsv_variants)){
    write.table(all_tier_tsv_variants,file=all_tier_tsv_fname, sep="\t",col.names = T,row.names = F,quote = F)
  }

}

#' Function that plots call numbers
#'
#' @param sample_calls data frame with sample variants
#' @param title Title for plot
#'
#' @return p
#' @export
#'
plot_call_statistics <- function(sample_calls, title = 'Title'){
  plot_df <- sample_calls %>% dplyr::select(VARIANT_CLASS, CALL_CONFIDENCE)
  plot_df$VARIANT_CLASS <- toupper(plot_df$VARIANT_CLASS)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = VARIANT_CLASS, fill=CALL_CONFIDENCE)) + ggplot2::geom_bar(stat="count",position=ggplot2::position_dodge()) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(title) +
    ggplot2::scale_color_brewer(palette='Dark2') +
    ggplot2::theme(
      title = ggplot2::element_text(size=10, family="Helvetica"),
      axis.title.x = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size=8, family="Helvetica"),
      axis.text.x = ggplot2::element_text(family="Helvetica",size=8),
      axis.text.y = ggplot2::element_text(family="Helvetica",size=8),
      axis.title.y = ggplot2::element_text(family="Helvetica",size=8,vjust=1.5),
      plot.margin = (ggplot2::unit(c(1, 0, 1, 1), "cm")),
      plot.title = ggplot2::element_text(family="Helvetica",size=8,vjust=2)) + ggplot2::ylab('Count')
  return(p)
}


#' Function that annotates CNV segment files (FACETS)
#'
#' @param cnv_file from FACETS
#' @param format
#'
#' @return cnv_data
#' @export
#'

cnv_segment_annotation <- function(cnv_file, format = 'facets'){

  cnv_df <- NULL
  if(format == 'facets'){
    cnv_df <- read.table(file=cnv_file,header = T,stringsAsFactors = F,comment.char="")
    cnv_df <- dplyr::mutate(cnv_df, chromosome = X.chromosome, LogR = CN_logR_median, cf = cellular_fraction_MM_EM_optimized, tmp_cnTotal = total_copy_number_MM_EM_optimized, tmp_cnMinor = minor_copy_number_MM_EM_optimized)
    cnv_df <- as.data.frame(cnv_df %>% dplyr::rowwise() %>% dplyr::mutate(cnTotal = max(tmp_cnTotal,tmp_cnMinor,na.rm=T)))
    cnv_df <- as.data.frame(cnv_df %>% dplyr::rowwise() %>% dplyr::mutate(cnMinor = min(tmp_cnTotal,tmp_cnMinor)))
    cnv_df <- dplyr::select(cnv_df, -c(tmp_cnTotal,tmp_cnMinor))
    cnv_df$LogR <- round(cnv_df$LogR, digits = 3)
    cnv_df$cf <- round(as.numeric(cnv_df$cf), digits = 3)
    cnv_df$chromosome <- paste0("chr",cnv_df$chromosome)
  }
  else{
    if(format == 'tcga'){
      cnv_df <- read.table(file=cnv_file,header = T,stringsAsFactors = F,comment.char="")
      cnv_df <- dplyr::rename(cnv_df, chromosome = Chromosome, LogR = Segment_Mean, segment_start = Start, segment_end = End)
      cnv_df <- dplyr::select(cnv_df, -c(Sample,Num_Probes)) %>% dplyr::distinct()
      cnv_df$chromosome <- paste0("chr",cnv_df$chromosome)
      cnv_df$cnTotal <- NA
      cnv_df$cnTotal <- as.integer(round((2^cnv_df$LogR) * 2, digits = 0))
      cnv_df$cnMinor <- NA
      cnv_df$cf <- NA
    }
  }

  if(nrow(cnv_df[cnv_df$chromosome == 'chr23',])){
    cnv_df[cnv_df$chromosome == 'chr23',]$chromosome <- 'chrX'
  }
  if(nrow(cnv_df[cnv_df$chromosome == 'chr24',])){
    cnv_df[cnv_df$chromosome == 'chr24',]$chromosome <- 'chrY'
  }
  seqinfo_hg19 <- GenomeInfoDb::Seqinfo(seqnames = GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19)), seqlengths = GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19)), genome = 'hg19')

  cnv_gr <- GenomicRanges::makeGRangesFromDataFrame(cnv_df, keep.extra.columns = T, seqinfo = seqinfo_hg19, seqnames.field = 'chromosome',start.field = 'segment_start', end.field = 'segment_end', ignore.strand = T, starts.in.df.are.0based = T)

  ensembl_genes_xref <- dplyr::filter(gene_xref, !is.na(chromosome_name) & !is.na(transcript_start) & !is.na(transcript_end))
  ensembl_genes_xref <- dplyr::select(ensembl_genes_xref, chromosome_name,transcript_start,transcript_end,ensembl_gene_id,ensembl_transcript_id,entrezgene,gencode_v19,gene_biotype,symbol,name,cancer_census_somatic,cancer_census_germline,tsgene,ts_oncogene,intogen_drivers,antineoplastic_drugs_dgidb)
  ensembl_genes_gr <- GenomicRanges::makeGRangesFromDataFrame(ensembl_genes_xref, keep.extra.columns = T, seqinfo = seqinfo_hg19, seqnames.field = 'chromosome_name', start.field = 'transcript_start', end.field = 'transcript_end', ignore.strand = T, starts.in.df.are.0based = T)

  hits <- GenomicRanges::findOverlaps(ensembl_genes_gr, cnv_gr, type="within", select="all")
  ranges <- ensembl_genes_gr[queryHits(hits)]
  mcols(ranges) <- c(mcols(ranges),mcols(cnv_gr[subjectHits(hits)]))

  df <- as.data.frame(mcols(ranges))
  df$segment_start <- start(ranges(cnv_gr[subjectHits(hits)]))
  df$segment_end <- end(ranges(cnv_gr[subjectHits(hits)]))
  df$segment_length <- paste(round((as.numeric((df$segment_end - df$segment_start)/1000000)),digits = 2),"Mb")
  df$transcript_start <- start(ranges)
  df$transcript_end <- end(ranges)
  df$chrom <- as.character(seqnames(ranges))
  df$segment_link <- paste0("<a href='",paste0('http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=',paste0(df$chrom,':',df$segment_start,'-',df$segment_end)),"' target=\"_blank\">",paste0(df$chrom,':',df$segment_start,'-',df$segment_end),"</a>")
  df_print <- df
  df_print <- dplyr::select(df_print,chrom,segment_start,segment_end,cf,LogR,cnTotal,cnMinor,ensembl_gene_id,symbol,ensembl_transcript_id,transcript_start,transcript_end,name,gene_biotype,cancer_census_germline,cancer_census_somatic,tsgene,ts_oncogene,intogen_drivers,antineoplastic_drugs_dgidb,gencode_v19)


  chrOrder <- c(as.character(paste0('chr',c(1:22))),"chrX","chrY")
  df_print$chrom <- factor(df_print$chrom, levels=chrOrder)
  df_print <- df_print[order(df_print$chrom),]
  df_print$segment_start <- as.integer(df_print$segment_start)
  df_print$segment_end <- as.integer(df_print$segment_end)

  df_print_sorted <- NULL
  for(chrom in chrOrder){
    if(nrow(df_print[df_print$chrom == chrom,]) > 0){
      chrom_regions <- df_print[df_print$chrom == chrom,]
      chrom_regions_sorted <- chrom_regions[with(chrom_regions, order(segment_start, segment_end)),]
      df_print_sorted <- rbind(df_print_sorted, chrom_regions_sorted)
    }
  }

  df <- dplyr::select(df, -ensembl_transcript_id) %>% dplyr::filter(gene_biotype == 'protein_coding') %>% dplyr::distinct()
  df$cancer_census_somatic <- stringr::str_replace_all(df$cancer_census_somatic,"&",", ")

  df <- dplyr::rename(df, ANTINEOPLASTIC_DRUG_INTERACTION = antineoplastic_drugs_dgidb)
  df$VAR_ID <- rep(1:nrow(df))
  df <- OncoVarReporter::annotate_variant_link(df, vardb = 'DGIDB')
  df <- dplyr::rename(df, ONCOGENE = ts_oncogene, TUMOR_SUPPRESSOR = tsgene, ENTREZ_ID = entrezgene, CANCER_CENSUS_SOMATIC = cancer_census_somatic, GENE = symbol, CHROMOSOME = chrom, GENENAME = name, ANTINEOPLASTIC_DRUG_INTERACTIONS = DGIDBLINK, SEGMENT_LENGTH = segment_length, SEGMENT = segment_link)
  df <- OncoVarReporter::annotate_variant_link(df, vardb = 'NCBI_GENE')
  df <- dplyr::rename(df, GENE_NAME = NCBI_GENE_LINK)

  df <- dplyr::select(df, CHROMOSOME, GENE, GENE_NAME, CANCER_CENSUS_SOMATIC, TUMOR_SUPPRESSOR, ONCOGENE, ANTINEOPLASTIC_DRUG_INTERACTIONS,SEGMENT_LENGTH, SEGMENT,cnTotal,cnMinor,LogR,cf) %>% dplyr::distinct()
  df <- df %>% dplyr::distinct()

  segments <- NULL
  segments <- dplyr::select(df, SEGMENT, SEGMENT_LENGTH, LogR, cnTotal, cnMinor,cf) %>% dplyr::distinct()
  segments <- dplyr::filter(segments, cnTotal != 2)
  segments <- segments %>% dplyr::arrange(desc(LogR))

  oncogene_amplified <- NULL
  oncogene_amplified <- dplyr::filter(df, !is.na(ONCOGENE) & (cnTotal >= 5))
  oncogene_amplified <- dplyr::select(oncogene_amplified, -c(TUMOR_SUPPRESSOR, ONCOGENE))
  oncogene_amplified <- oncogene_amplified %>% dplyr::arrange(ANTINEOPLASTIC_DRUG_INTERACTIONS)
  if(nrow(oncogene_amplified) > 0){
    oncogene_amplified$CNA_TYPE <- 'gain'
  }
  tsgene_homozygous_deletion <- NULL
  tsgene_homozygous_deletion <- dplyr::filter(df, !is.na(TUMOR_SUPPRESSOR) & (cnTotal == 0))
  tsgene_homozygous_deletion <- dplyr::select(tsgene_homozygous_deletion, -c(TUMOR_SUPPRESSOR, ONCOGENE))
  tsgene_homozygous_deletion <- tsgene_homozygous_deletion %>% dplyr::arrange(CANCER_CENSUS_SOMATIC)
  if(nrow(tsgene_homozygous_deletion) > 0){
    tsgene_homozygous_deletion$CNA_TYPE <- 'loss'
  }
  civic_cna_biomarkers <- dplyr::filter(civic_biomarkers, alteration_type == 'CNA') %>% dplyr::select(genesymbol,evidence_type,evidence_level,evidence_description,disease_name,evidence_direction,pubmed_html_link,drug_names,rating,clinical_significance,civic_consequence)
  names(civic_cna_biomarkers) <- toupper(names(civic_cna_biomarkers))
  civic_cna_biomarkers <- dplyr::rename(civic_cna_biomarkers, GENE = GENESYMBOL, CNA_TYPE = CIVIC_CONSEQUENCE, DESCRIPTION = EVIDENCE_DESCRIPTION, CITATION = PUBMED_HTML_LINK)

  cna_biomarkers <- NULL
  if(!is.null(tsgene_homozygous_deletion)){
    if(nrow(tsgene_homozygous_deletion) > 0){
      civic_biomarker_hits1 <- dplyr::inner_join(tsgene_homozygous_deletion, civic_cna_biomarkers)
      cna_biomarkers <- rbind(cna_biomarkers,civic_biomarker_hits1)
    }
  }
  if(!is.null(oncogene_amplified)){
    if(nrow(oncogene_amplified) > 0){
      civic_biomarker_hits2 <- dplyr::inner_join(oncogene_amplified, civic_cna_biomarkers)
      cna_biomarkers <- rbind(cna_biomarkers,civic_biomarker_hits2)
    }
  }

  if(!is.null(cna_biomarkers)){
    cna_biomarkers <- cna_biomarkers[c("CHROMOSOME","GENE","CNA_TYPE","EVIDENCE_LEVEL","CLINICAL_SIGNIFICANCE","EVIDENCE_TYPE","DESCRIPTION","DISEASE_NAME","EVIDENCE_DIRECTION","DRUG_NAMES","CITATION","RATING","GENE_NAME","CANCER_CENSUS_SOMATIC","ANTINEOPLASTIC_DRUG_INTERACTIONS","SEGMENT_LENGTH", "SEGMENT","LogR","cnTotal","cf")]
    cna_biomarkers <- cna_biomarkers %>% dplyr::arrange(EVIDENCE_LEVEL,RATING)
  }

  cnv_data <- list('ranked_segments' = segments, 'oncogene_amplified' = oncogene_amplified, 'tsgene_homozygous_deletion' = tsgene_homozygous_deletion,'cnv_df_for_print' = df_print_sorted, 'cna_biomarkers' = cna_biomarkers)
  return(cnv_data)
}

#' Function that generates dense and tiered annotated variant datasets
#'
#' @param report_data report data structure generated from pcgr
#' @param tier1_variants
#' @param tier2_variants
#' @param tier3_variants
#' @param tier4_variants
#' @param tier5_variants
#' @param sample_id Sample identifier
#'
#' @return tsv_variants data frame with tier-annotated list of variants for tab-separated output
#' @export
#'
generate_tier_tsv <- function(tier1_variants, tier2_variants, tier3_variants, tier4_variants, tier5_variants, sample_id = 'test'){

  tier1_tsv <- tier1_variants
  tsv_variants <- NULL
  if(nrow(tier1_tsv) > 0){
    tier1_tsv$TIER <- 'TIER 1'
    tier1_tsv$TIER_DESCRIPTION <- 'Clinical biomarker - prognostic/diagnostic/drug sensitivity/resistance'
    tier1_tsv$VCF_SAMPLE_ID <- sample_id
    tsv_variants <- rbind(tsv_variants, dplyr::select(tier1_tsv, dplyr::one_of(pcgr_tsv_tiered_columns))) %>% dplyr::distinct()
  }
  tier2_tsv <- tier2_variants
  if(nrow(tier2_tsv) > 0){
    tier2_tsv$TIER <- 'TIER 2'
    tier2_tsv$TIER_DESCRIPTION <- 'Other cancer mutation hotspot/predicted driver mutation/curated cancer-associated mutation'
    tier2_tsv$VCF_SAMPLE_ID <- sample_id
    tsv_variants <- rbind(tsv_variants, dplyr::select(tier2_tsv, dplyr::one_of(pcgr_tsv_tiered_columns)))
  }
  tier3_tsv <- tier3_variants
  if(nrow(tier3_tsv) > 0){
    tier3_tsv$TIER <- 'TIER 3'
    tier3_tsv$TIER_DESCRIPTION <- 'Other cancer census gene/proto-oncogene/tumor suppressor mutation'
    tier3_tsv$VCF_SAMPLE_ID <- sample_id
    tsv_variants <- rbind(tsv_variants, dplyr::select(tier3_tsv, dplyr::one_of(pcgr_tsv_tiered_columns)))
  }
  tier4_tsv <- tier4_variants
  if(nrow(tier4_tsv) > 0){
    tier4_tsv$TIER <- 'TIER 4'
    tier4_tsv$VCF_SAMPLE_ID <- sample_id
    tier4_tsv$TIER_DESCRIPTION <- 'Other coding mutation'
    tsv_variants <- rbind(tsv_variants, dplyr::select(tier4_tsv, dplyr::one_of(pcgr_tsv_tiered_columns)))
  }
  tier5_tsv <- tier5_variants
  if(nrow(tier5_tsv) > 0){
    tier5_tsv$TIER <- 'TIER 4'
    tier5_tsv$VCF_SAMPLE_ID <- sample_id
    tier5_tsv$TIER_DESCRIPTION <- 'Non-coding mutation'
    tsv_variants <- rbind(tsv_variants, dplyr::select(tier5_tsv, dplyr::one_of(pcgr_tsv_tiered_columns)))
  }
  tsv_variants$COSMIC <- unlist(lapply(stringr::str_match_all(tsv_variants$COSMIC,"COSM[0-9]{1,}"),paste,collapse=","))
  tsv_variants$DBSNP <- unlist(lapply(stringr::str_match_all(tsv_variants$DBSNP,">rs[0-9]{1,}<"),paste,collapse=","))
  tsv_variants$DBSNP <- stringr::str_replace_all(tsv_variants$DBSNP,">|<","")
  tsv_variants$GENE_NAME <- unlist(lapply(stringr::str_match_all(tsv_variants$GENE_NAME,">.+<"),paste,collapse=","))
  tsv_variants$GENE_NAME <- stringr::str_replace_all(tsv_variants$GENE_NAME,">|<","")
  tsv_variants$CLINVAR <- unlist(lapply(stringr::str_match_all(tsv_variants$CLINVAR,">.+<"),paste,collapse=","))
  tsv_variants$CLINVAR <- stringr::str_replace_all(tsv_variants$CLINVAR,">|<","")
  tsv_variants$PROTEIN_DOMAIN <- unlist(lapply(stringr::str_match_all(tsv_variants$PROTEIN_DOMAIN,">.+<"),paste,collapse=","))
  tsv_variants$PROTEIN_DOMAIN <- stringr::str_replace_all(tsv_variants$PROTEIN_DOMAIN,">|<","")

  return(tsv_variants)
}


#' Function that generates multiple reports
#'
#' @param sample_calls data frame with list of variant calls
#' @param sample_id sample identifier
#' @param minimum_n_signature_analysis minimum number of mutations for signature analysis
#' @param signatures.limit limit the number of possible mutational signatures
#'
#' @return report_data data frame with all report elements
#' @export
#'
generate_report_data <- function(sample_calls, sample_id = NULL, minimum_n_signature_analysis = 50, signatures.limit = 6){

  tier1_report <- FALSE
  tier2_report <- FALSE
  tier3_report <- FALSE
  tier4_report <- FALSE
  tier5_report <- FALSE
  clinical_evidence_items_tier1A <- data.frame()
  clinical_evidence_items_tier1B <- data.frame()
  clinical_evidence_items_tier1C <- data.frame()
  variants_tier1 <- data.frame()
  variants_tier2 <- data.frame()
  variants_tier3 <- data.frame()
  variants_tier4 <- data.frame()
  variants_tier5 <- data.frame()

  signature_report <- FALSE
  signature_call_set <- data.frame()

  sample_calls_coding <- sample_calls %>% dplyr::filter(stringr::str_detect(CONSEQUENCE,"stop_gained|stop_lost|start_lost|frameshift_variant|missense_variant|splice_donor|splice_acceptor|inframe_deletion|inframe_insertion"))
  sample_calls_noncoding <- sample_calls %>% dplyr::filter(!stringr::str_detect(CONSEQUENCE,"stop_gained|stop_lost|start_lost|frameshift_variant|missense_variant|splice_donor|splice_acceptor|inframe_deletion|inframe_insertion"))

  sample_stats_plot_all <- OncoVarReporter::plot_call_statistics(sample_calls,"Somatic calls - all")
  sample_stats_plot_coding <- OncoVarReporter::plot_call_statistics(sample_calls_coding,"Somatic calls - coding")

  min_variants_for_signature <- minimum_n_signature_analysis
  if(any(grepl(paste0("VARIANT_CLASS$"),names(sample_calls)))){
    if(nrow(sample_calls[sample_calls$VARIANT_CLASS == 'SNV',]) >= min_variants_for_signature){
      signature_call_set <- sample_calls[sample_calls$VARIANT_CLASS == 'SNV',]
      signature_call_set <- dplyr::filter(signature_call_set, CHROM != 'MT')
      signature_call_set$VCF_SAMPLE_ID <- sample_id
      signature_report <- TRUE
    }
  }

  clinical_evidence_items_tier1A <- OncoVarReporter::get_clinical_associations_civic_cbmdb(sample_calls_coding)
  clinical_evidence_items_tier1B <- OncoVarReporter::get_clinical_associations_civic_cbmdb(sample_calls_coding, mapping = 'codon')
  clinical_evidence_items_tier1C <- OncoVarReporter::get_clinical_associations_civic_cbmdb(sample_calls_coding, mapping = 'exon')
  variants_tier1 <- rbind(clinical_evidence_items_tier1A,clinical_evidence_items_tier1B,clinical_evidence_items_tier1C)
  variants_tier1_tsv <- variants_tier1
  if(nrow(clinical_evidence_items_tier1B) > 0){
    clinical_evidence_items_tier1B <- dplyr::select(clinical_evidence_items_tier1B, -ONCOSCORE)
  }
  if(nrow(clinical_evidence_items_tier1A) > 0){
    clinical_evidence_items_tier1A <- dplyr::select(clinical_evidence_items_tier1A, -ONCOSCORE)
  }
  if(nrow(clinical_evidence_items_tier1C) > 0){
    clinical_evidence_items_tier1C <- dplyr::select(clinical_evidence_items_tier1C, -ONCOSCORE)
  }
  if(nrow(variants_tier1) > 0){
    variants_tier1 <- variants_tier1 %>% dplyr::select(GENOMIC_CHANGE) %>% dplyr::distinct()
    tier1_report <- TRUE
  }

  ## Analyze Tier 2: curated mutations, cancer mutation hotspots and predicted driver mutations
  tier2_tags <- c("SYMBOL","ONCOSCORE","PROTEIN_CHANGE","GENE_NAME","PROTEIN_DOMAIN","PROTEIN_FEATURE","CDS_CHANGE","CANCER_MUTATION_HOTSPOT","CANCER_CENSUS_SOMATIC","INTOGEN_DRIVER_MUT","CONSEQUENCE","EFFECT_PREDICTIONS","DBSNP","COSMIC","COSMIC_SITE_HISTOLOGY","CLINVAR","DOCM_DISEASE","ANTINEOPLASTIC_DRUG_INTERACTIONS","VEP_ALL_CONSEQUENCE","GENOME_VERSION","GENOMIC_CHANGE","CALL_CONFIDENCE","DP_TUMOR","AF_TUMOR","DP_NORMAL","AF_NORMAL","OXIDATION_ARTEFACT")
  variants_tier2 <- dplyr::select(sample_calls_coding, dplyr::one_of(tier2_tags))
  variants_tier2 <- variants_tier2 %>% dplyr::filter(!is.na(INTOGEN_DRIVER_MUT) | !is.na(CANCER_MUTATION_HOTSPOT) | !is.na(DOCM_DISEASE))
  if(nrow(variants_tier1) > 0){
    variants_tier2 <- dplyr::anti_join(variants_tier2, variants_tier1, by=c("GENOMIC_CHANGE"))
  }
  tier12 <- variants_tier1
  if(nrow(variants_tier2) > 0){
    tier2_report <- TRUE
    if(nrow(variants_tier2[is.na(variants_tier2$ONCOSCORE),]) > 0){
      variants_tier2[is.na(variants_tier2$ONCOSCORE),]$ONCOSCORE <- 0
    }
    variants_tier2 <- variants_tier2 %>% dplyr::arrange(desc(ONCOSCORE))
    tier12 <- rbind(variants_tier1,dplyr::select(variants_tier2,GENOMIC_CHANGE)) %>% dplyr::distinct()
  }

  ## Analyze Tier 3: coding mutations in oncogenes/tumor suppressors/cancer census genes
  tier3_tags <- c("SYMBOL","ONCOSCORE","PROTEIN_CHANGE","GENE_NAME","PROTEIN_DOMAIN","PROTEIN_FEATURE","CDS_CHANGE","CANCER_MUTATION_HOTSPOT","INTOGEN_DRIVER_MUT","CONSEQUENCE","EFFECT_PREDICTIONS","DBSNP","COSMIC","COSMIC_SITE_HISTOLOGY","CLINVAR","ANTINEOPLASTIC_DRUG_INTERACTIONS","VEP_ALL_CONSEQUENCE","GENOME_VERSION","GENOMIC_CHANGE","ONCOGENE","TUMOR_SUPPRESSOR","CANCER_CENSUS_SOMATIC","CALL_CONFIDENCE","DP_TUMOR","AF_TUMOR","DP_NORMAL","AF_NORMAL","OXIDATION_ARTEFACT")
  variants_tier3 <- dplyr::select(sample_calls_coding, dplyr::one_of(tier3_tags))
  variants_tier3 <- variants_tier3 %>% dplyr::filter(!is.na(CANCER_CENSUS_SOMATIC) | ONCOGENE == TRUE | TUMOR_SUPPRESSOR == TRUE)
  if(nrow(tier12) > 0){
    variants_tier3 <- dplyr::anti_join(variants_tier3,tier12, by=c("GENOMIC_CHANGE"))
  }
  tier123 <- tier12
  if(nrow(variants_tier3) > 0){
    tier3_report <- TRUE
    if(nrow(variants_tier3[is.na(variants_tier3$ONCOSCORE),]) > 0){
      variants_tier3[is.na(variants_tier3$ONCOSCORE),]$ONCOSCORE <- 0
    }
    variants_tier3 <- variants_tier3 %>% dplyr::arrange(desc(ONCOSCORE))
    tier123 <- rbind(tier12,dplyr::select(variants_tier3,GENOMIC_CHANGE)) %>% dplyr::distinct()
  }

  ## Analyze Tier 4: Other coding mutations
  tier4_tags <- c("SYMBOL","ONCOSCORE","PROTEIN_CHANGE","GENE_NAME","PROTEIN_DOMAIN","PROTEIN_FEATURE","CDS_CHANGE","CONSEQUENCE","EFFECT_PREDICTIONS","DBSNP","COSMIC","COSMIC_SITE_HISTOLOGY","CLINVAR","ANTINEOPLASTIC_DRUG_INTERACTIONS","GENOME_VERSION","VEP_ALL_CONSEQUENCE","GENOMIC_CHANGE","CALL_CONFIDENCE","DP_TUMOR","AF_TUMOR","DP_NORMAL","AF_NORMAL","OXIDATION_ARTEFACT")
  variants_tier4 <- dplyr::select(sample_calls_coding, dplyr::one_of(tier4_tags))
  if(nrow(tier123) > 0){
    variants_tier4 <- dplyr::anti_join(variants_tier4,tier123, by=c("GENOMIC_CHANGE"))
  }
  if(nrow(variants_tier4) > 0){
    if(nrow(variants_tier4[is.na(variants_tier4$ONCOSCORE),]) > 0){
      variants_tier4[is.na(variants_tier4$ONCOSCORE),]$ONCOSCORE <- 0
    }
    variants_tier4 <- variants_tier4 %>% dplyr::arrange(desc(ONCOSCORE))
    tier4_report <- TRUE
  }


  ## Analyze Tier 5: Non-coding mutations
  tier5_tags <- c("SYMBOL","ONCOSCORE","CONSEQUENCE","PROTEIN_CHANGE","GENE_NAME","PROTEIN_DOMAIN","PROTEIN_FEATURE","DBSNP","COSMIC","CDS_CHANGE","EFFECT_PREDICTIONS","COSMIC_SITE_HISTOLOGY","CLINVAR","ANTINEOPLASTIC_DRUG_INTERACTIONS","GENOME_VERSION","VEP_ALL_CONSEQUENCE","GENOMIC_CHANGE","CALL_CONFIDENCE","DP_TUMOR","AF_TUMOR","DP_NORMAL","AF_NORMAL","OXIDATION_ARTEFACT")
  variants_tier5 <- dplyr::select(sample_calls_noncoding, dplyr::one_of(tier5_tags))
  #variants_tier5 <- dplyr::filter(variants_tier5, !(stringr::str_detect(CONSEQUENCE,"intron_variant|intergenic_variant") & !stringr::str_detect(CONSEQUENCE,"splice_")))

  if(nrow(variants_tier5) > 0){
    if(nrow(variants_tier5[is.na(variants_tier5$ONCOSCORE),]) > 0){
      variants_tier5[is.na(variants_tier5$ONCOSCORE),]$ONCOSCORE <- 0
    }
    variants_tier5 <- variants_tier5 %>% dplyr::arrange(desc(ONCOSCORE))
    tier5_report <- TRUE
  }

  tsv_variants <- OncoVarReporter::generate_tier_tsv(variants_tier1_tsv,
                                                     variants_tier2,
                                                     variants_tier3,
                                                     variants_tier4,
                                                     variants_tier5,
                                                     sample_id = sample_id)

  variants_tier5 <- dplyr::select(variants_tier5,-c(PROTEIN_CHANGE,PROTEIN_DOMAIN,PROTEIN_FEATURE,COSMIC_SITE_HISTOLOGY,ANTINEOPLASTIC_DRUG_INTERACTIONS))
  variants_tier2 <- dplyr::select(variants_tier5, -ONCOSCORE)
  report_data <- list('sample_stats_plot_all' = sample_stats_plot_all, 'sample_stats_plot_coding' = sample_stats_plot_coding,'tier1_report' = tier1_report, 'tier2_report' = tier2_report, 'tier3_report' = tier3_report, 'tier4_report' = tier4_report, 'tier5_report' = tier5_report, 'clinical_evidence_items_tier1A' = clinical_evidence_items_tier1A, 'clinical_evidence_items_tier1B' = clinical_evidence_items_tier1B, 'clinical_evidence_items_tier1C' = clinical_evidence_items_tier1C, 'tsv_variants' = tsv_variants, 'variants_tier1' = variants_tier1, 'variants_tier2' = variants_tier2, 'variants_tier3' = variants_tier3, 'variants_tier4' = variants_tier4,'variants_tier5' = variants_tier5, 'signature_report' = signature_report, 'signature_call_set' = signature_call_set, 'signatures.limit' = signatures.limit, 'sample_name' = sample_id)


  return(report_data)

}

#' Function that adds HTML links to different genetic variant identifiers
#'
#' @param var_df data frame with variants
#' @param vardb type of variant database for which HTML links is to be provided
#' @param linktype type of link to be generated
#'
#' @return var_df
#' @export
#'
annotate_variant_link <- function(var_df, vardb = "DBSNP", linktype = "dbsource"){
  if(vardb == 'DBSNP'){

    if(any(grepl(paste0("^DBSNPRSID$"),names(var_df))) & any(grepl(paste0("^VAR_ID$"),names(var_df)))){
      var_df_unique_slim <- dplyr::select(var_df, DBSNPRSID, VAR_ID) %>% dplyr::filter(!is.na(DBSNPRSID)) %>% dplyr::distinct()
      if(nrow(var_df_unique_slim) > 0){
        var_df_unique_slim_melted <- var_df_unique_slim %>% tidyr::separate_rows(DBSNPRSID,sep=",")
        if(linktype == "dbsource"){
          var_df_unique_slim_melted <- var_df_unique_slim_melted %>% dplyr::mutate(tmp_dbsnp = paste0("<a href='http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=",DBSNPRSID,"' target=\"_blank\">rs",DBSNPRSID,"</a>"))
        }

        var_df_links <- dplyr::group_by(var_df_unique_slim_melted, VAR_ID) %>% dplyr::summarise(DBSNPLINK = unlist(paste(tmp_dbsnp, collapse = ", ")))
        var_df_links <- dplyr::select(var_df_links, VAR_ID, DBSNPLINK)
        var_df <- dplyr::left_join(var_df, var_df_links,by=c("VAR_ID" = "VAR_ID"))
      }
      else{
        var_df$DBSNPLINK <- NA
      }
    }
    else{
      cat("WARNING: Could not generate DBSNP links - no DBSNPRSID provided in annotated VCF",sep="\n")
      var_df$DBSNPLINK <- NA
    }
  }
  if(vardb == 'CLINVAR'){

    if(any(grepl(paste0("^CLINVAR_MSID$"),names(var_df))) & any(grepl(paste0("^VAR_ID$"),names(var_df))) & any(grepl(paste0("^CLINVAR_TRAITS_ALL$"),names(var_df)))){
      var_df_unique_slim <- dplyr::select(var_df, VAR_ID, CLINVAR_MSID, CLINVAR_TRAITS_ALL) %>% dplyr::filter(!is.na(CLINVAR_MSID)) %>% dplyr::distinct()
      if(nrow(var_df_unique_slim) > 0){
        if(linktype == "dbsource"){
          var_df_unique_slim <- var_df_unique_slim %>% dplyr::mutate(CLINVARLINK = paste0("<a href='http://www.ncbi.nlm.nih.gov/clinvar/variation/",CLINVAR_MSID,"' target=\"_blank\">",CLINVAR_TRAITS_ALL,"</a>"))
        }
        var_df_links <- var_df_unique_slim %>% dplyr::select(VAR_ID, CLINVARLINK)
        var_df <- dplyr::left_join(var_df, var_df_links,by=c("VAR_ID" = "VAR_ID"))
      }
      else{
        var_df$CLINVARLINK <- NA
      }
    }
    else{
      cat("WARNING: Could not generate CLINVAR links - no CLINVAR_MSID provided in annotated VCF",sep="\n")
      var_df$CLINVARLINK <- NA
    }
  }

  if(vardb == 'NCBI_GENE'){

    if(any(grepl(paste0("^GENENAME$"),names(var_df))) & any(grepl(paste0("^VAR_ID$"),names(var_df))) & any(grepl(paste0("^ENTREZ_ID$"),names(var_df)))){
      var_df_unique_slim <- dplyr::select(var_df, VAR_ID, GENENAME, ENTREZ_ID) %>% dplyr::filter(!is.na(ENTREZ_ID)) %>% dplyr::distinct()
      if(nrow(var_df_unique_slim) > 0){
        if(linktype == "dbsource"){
          var_df_unique_slim <- var_df_unique_slim %>% dplyr::mutate(NCBI_GENE_LINK = paste0("<a href='http://www.ncbi.nlm.nih.gov/gene/",ENTREZ_ID,"' target=\"_blank\">",GENENAME,"</a>"))
        }
        var_df_links <- var_df_unique_slim %>% dplyr::select(VAR_ID, NCBI_GENE_LINK)
        var_df <- dplyr::left_join(var_df, var_df_links,by=c("VAR_ID" = "VAR_ID"))
      }
      else{
        var_df$NCBI_GENE_LINK <- NA
      }
    }
    else{
      cat("WARNING: Could not generate NCBI_GENE links - no ENTREZ_ID provided in annotated VCF",sep="\n")
      var_df$NCBI_GENE_LINK <- NA
    }
  }

  if(vardb == 'COSMIC'){

    if(any(grepl(paste0("^COSMIC_MUTATION_ID$"),names(var_df))) & any(grepl(paste0("^VAR_ID$"),names(var_df)))){
      var_df_unique_slim <- dplyr::select(var_df, VAR_ID, COSMIC_MUTATION_ID) %>% dplyr::filter(!is.na(COSMIC_MUTATION_ID)) %>% dplyr::distinct()
      if(nrow(var_df_unique_slim) > 0){
        var_df_unique_slim_melted <- var_df_unique_slim %>% tidyr::separate_rows(COSMIC_MUTATION_ID,sep="&")
        if(linktype == "dbsource"){
          var_df_unique_slim_melted <- var_df_unique_slim_melted %>% dplyr::mutate(tmp_cosmic = paste0("<a href='http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=",stringr::str_replace(COSMIC_MUTATION_ID,"COSM",""),"' target=\"_blank\">",COSMIC_MUTATION_ID,"</a>"))
        }
        var_df_links <- dplyr::group_by(var_df_unique_slim_melted, VAR_ID) %>% dplyr::summarise(COSMICLINK = unlist(paste(tmp_cosmic, collapse = ", ")))
        var_df_links <- dplyr::select(var_df_links, VAR_ID, COSMICLINK)
        var_df <- dplyr::left_join(var_df, var_df_links,by=c("VAR_ID" = "VAR_ID"))
      }
      else{
        var_df$COSMICLINK <- NA
      }
    }
    else{
      cat("WARNING: Could not add COSMIC links - no COSMIC_MUTATION_ID provided in annotated VCF",sep="\n")
      var_df$COSMICLINK <- NA
    }
  }
  if(vardb == 'UCSC'){
    if(any(grepl(paste0("^GENOMIC_CHANGE$"),names(var_df))) & any(grepl(paste0("^VAR_ID$"),names(var_df)))){
      var_df_unique <- dplyr::select(var_df, VAR_ID, GENOMIC_CHANGE) %>% dplyr::filter(!is.na(GENOMIC_CHANGE)) %>% dplyr::distinct()
      var_df_unique$CHROM <- stringr::str_split_fixed(var_df_unique$VAR_ID,"_",4)[,1]
      var_df_unique$POS <- stringr::str_split_fixed(var_df_unique$VAR_ID,"_",4)[,2]
      var_df_unique$UCSC_LINK <- paste0("<a href='",paste0('http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=',paste0('chr',var_df_unique$CHROM,':',var_df_unique$POS,'-',var_df_unique$POS)),"' target=\"_blank\">",var_df_unique$GENOMIC_CHANGE,"</a>")

      var_df_links <- dplyr::select(var_df_unique, VAR_ID, UCSC_LINK)
      var_df <- dplyr::left_join(var_df, var_df_links,by=c("VAR_ID" = "VAR_ID"))
    }
  }


  if(vardb == 'DGIDB'){
    if(any(grepl(paste0("^ANTINEOPLASTIC_DRUG_INTERACTION$"),names(var_df))) & any(grepl(paste0("^VAR_ID$"),names(var_df)))){
      var_df_unique_slim <- dplyr::select(var_df, VAR_ID, ANTINEOPLASTIC_DRUG_INTERACTION) %>% dplyr::filter(!is.na(ANTINEOPLASTIC_DRUG_INTERACTION)) %>% dplyr::distinct()
      if(nrow(var_df_unique_slim) > 0){
        var_df_unique_slim_melted <- var_df_unique_slim %>% tidyr::separate_rows(ANTINEOPLASTIC_DRUG_INTERACTION,sep="&")
        if(linktype == "dbsource"){
          var_df_unique_slim_melted <- var_df_unique_slim_melted %>% dplyr::mutate(tmp_dgidb = paste0("<a href='http://dgidb.genome.wustl.edu/drugs/",stringr::str_replace(ANTINEOPLASTIC_DRUG_INTERACTION,"\\S{1,}:",""),"' target=\"_blank\">",ANTINEOPLASTIC_DRUG_INTERACTION,"</a>"))
        }

        var_df_links <- dplyr::group_by(var_df_unique_slim_melted, VAR_ID) %>% dplyr::summarise(DGIDBLINK = unlist(paste(tmp_dgidb, collapse = ", ")))
        var_df_links <- dplyr::select(var_df_links, VAR_ID, DGIDBLINK)
        var_df <- dplyr::left_join(var_df, var_df_links,by=c("VAR_ID" = "VAR_ID"))
      }
      else{
        var_df$DGIDBLINK <- NA
      }
    }
    else{
      cat("WARNING: Could not generate DGIdb links - no DGIDB info provided in annotated VCF",sep="\n")
      var_df$DGIDBLINK <- NA
    }
  }

  return(var_df)

}

#' A function that converts the INFO tags in a VRanges object into a basic types
#'
#' @param vr VRanges object
#' @return vr Vranges object with INFO tags more simply formatted
#'
#' @export

postprocess_vranges_info <- function(vr){
  ## Convert IntegerLists and CharacterLists to basic character lists
  vcf_annotations_df <- NULL
  for(tag in colnames(GenomicRanges::mcols(vr))){
    mcol_class <- class(GenomicRanges::mcols(vr)[,c(tag)])[1]

    if(mcol_class != "character" & mcol_class != "integer" & mcol_class != "logical" & mcol_class != "numeric"){
      annotation_track <- NULL
      #cat("TAG: ",tag, ', type:',mcol_class,'\n')
      if(mcol_class == "CompressedCharacterList"){
        annotation_track <- data.frame(val = as.character(Biostrings::unstrsplit(GenomicRanges::mcols(vr)[[tag]], sep=',')))
      }else{
        annotation_track <- data.frame(val = as.character(sapply(GenomicRanges::mcols(vr)[[tag]], paste, collapse=",")))
      }
      if(is.null(vcf_annotations_df)){
        vcf_annotations_df <- data.frame(annotation_track$val)
        names(vcf_annotations_df) <- c(tag)
        vcf_annotations_df[,c(tag)] <- as.character(vcf_annotations_df[,c(tag)])
      }
      else{
        vcf_annotations_df[,c(tag)] <- as.character(annotation_track$val)
      }
      ## add NA to empty values
      if(nrow(as.data.frame(vcf_annotations_df[nchar(vcf_annotations_df[,c(tag)]) == 0,])) != 0){
        if(dim(vcf_annotations_df)[2] == 1){
          vcf_annotations_df[nchar(vcf_annotations_df[,c(tag)]) == 0,] <- NA
        }
        else{
          vcf_annotations_df[nchar(vcf_annotations_df[,c(tag)]) == 0,][,c(tag)] <- NA
        }
      }
    }
    else{
      #cat("TAG: ",tag, ', type:',mcol_class,'\n')
      if(is.null(vcf_annotations_df)){
        vcf_annotations_df <- data.frame(GenomicRanges::mcols(vr)[,c(tag)])
        names(vcf_annotations_df) <- c(tag)
      }
      else{
        vcf_annotations_df[,c(tag)] <- GenomicRanges::mcols(vr)[,c(tag)]
      }
    }
  }

  ## add variant_id and sample_id
  position <- GenomicRanges::start(GenomicRanges::ranges(vr))
  vcf_annotations_df['VAR_ID'] <- paste(as.character(GenomeInfoDb::seqnames(vr)),position,VariantAnnotation::ref(vr),VariantAnnotation::alt(vr),sep="_")
  vcf_annotations_df['VCF_SAMPLE_ID'] <- as.character(VariantAnnotation::sampleNames(vr))
  GenomicRanges::mcols(vr) <- S4Vectors::DataFrame(vcf_annotations_df)

  return(vr)

}

#' Function that adds PFAM name descriptions to PFAM identifiers
#'
#' @param vcf_data_df
#'
#' @return vcf_data_df_pfam
#' @export
#'
add_pfam_domain_links <- function(vcf_data_df){

  if("DOMAINS" %in% colnames(vcf_data_df) & "VAR_ID" %in% colnames(vcf_data_df)){
    pfam_df <- dplyr::select(vcf_data_df,DOMAINS,VAR_ID) %>% dplyr::filter(!is.na(DOMAINS))
    if(nrow(pfam_df) == 0){
      return(vcf_data_df)
    }
    pfam_df <- pfam_df %>% dplyr::distinct() %>% tidyr::separate_rows(DOMAINS,sep="&") %>% dplyr::filter(stringr::str_detect(DOMAINS,"Pfam_domain"))
    pfam_df$DOMAINS <- stringr::str_replace(pfam_df$DOMAINS,"Pfam_domain:","")
    pfam_df <- dplyr::left_join(pfam_df,pfam_domains,by=c("DOMAINS" = "pfam_id")) %>% dplyr::select(VAR_ID,url)
    pfam_df <- dplyr::rename(pfam_df, PD = url)
    pfam_ret <- as.data.frame(dplyr::group_by(pfam_df, VAR_ID) %>% dplyr::summarise(PROTEIN_DOMAIN = paste(PD, collapse=", ")))

    vcf_data_df <- dplyr::left_join(vcf_data_df,pfam_ret,by=c("VAR_ID" = "VAR_ID"))
  }

  return(vcf_data_df)
}

#' Function that appends clinical annotations for somatic cancer variants
#'
#' @param vcf_data_df data frame with variants
#' @param mapping - one of 'exact' (allele-specific) or 'approximate' (codon or exon-level biomarkers)
#' @param ncgc - logical indicating whether NCGC-specific tags are to be appended
#'
#' @return vcf_data_df
#' @export
#'
get_clinical_associations_civic_cbmdb <- function(vcf_data_df, mapping = 'exact', variant_origin = 'Somatic Mutation', ncgc = FALSE){

  tags <- c('SYMBOL','ONCOSCORE','CONSEQUENCE','PROTEIN_CHANGE','GENE_NAME','PROTEIN_DOMAIN','PROTEIN_FEATURE','EFFECT_PREDICTIONS','COSMIC','COSMIC_SITE_HISTOLOGY','DBSNP','CLINVAR','EXON','CDS_CHANGE','CANCER_MUTATION_HOTSPOT','OTHER_LITERATURE_DOCM','OTHER_DISEASE_DOCM','INTOGEN_DRIVER_MUT','ANTINEOPLASTIC_DRUG_INTERACTIONS','VEP_ALL_CONSEQUENCE','CIVIC_ID','CIVIC_ID_2','CBMDB_ID','GENOME_VERSION','GENOMIC_CHANGE','CALL_CONFIDENCE','DP_TUMOR','AF_TUMOR','DP_NORMAL','AF_NORMAL')
  if("pubmed_html_link" %in% colnames(civic_biomarkers)){
    civic_biomarkers <- dplyr::rename(civic_biomarkers, citation = pubmed_html_link)
  }
  if("evidence_description" %in% colnames(civic_biomarkers)){
    civic_biomarkers <- dplyr::rename(civic_biomarkers, description = evidence_description)
  }
  if("pubmed_html_link" %in% colnames(cbmdb_biomarkers)){
    cbmdb_biomarkers <- dplyr::rename(cbmdb_biomarkers, citation = pubmed_html_link)
  }
  if("evidence_description" %in% colnames(cbmdb_biomarkers)){
    cbmdb_biomarkers <- dplyr::rename(cbmdb_biomarkers, description = evidence_description)
  }
  clinical_evidence_items <- data.frame()
  if(mapping == 'exact'){
    vcf_data_df_civic <- vcf_data_df %>% dplyr::filter(!is.na(CIVIC_ID))
    if(nrow(vcf_data_df_civic) > 0){
      tmp <- dplyr::select(vcf_data_df_civic,CIVIC_ID,VAR_ID)
      tmp <- tmp %>% tidyr::separate_rows(CIVIC_ID,sep=",")
      vcf_data_df_civic <- merge(tmp,dplyr::select(vcf_data_df_civic,-c(CIVIC_ID)),by.x = "VAR_ID",by.y = "VAR_ID")
      civic_calls <- dplyr::select(vcf_data_df_civic,dplyr::one_of(tags))
      eitems <- dplyr::left_join(civic_calls,dplyr::filter(dplyr::select(civic_biomarkers,-c(civic_exon,civic_consequence)),alteration_type == 'MUT'),by=c("CIVIC_ID" = "civic_id"))
      names(eitems) <- toupper(names(eitems))
      eitems$BIOMARKER_MAPPING <- 'exact'
      clinical_evidence_items <- rbind(clinical_evidence_items, eitems)
    }
    vcf_data_df_cbmdb <- vcf_data_df %>% dplyr::filter(is.na(CIVIC_ID) & !is.na(CBMDB_ID))
    if(nrow(vcf_data_df_cbmdb) > 0){
      tmp <- dplyr::select(vcf_data_df_cbmdb,CBMDB_ID,VAR_ID)
      tmp <- tmp %>% tidyr::separate_rows(CBMDB_ID,sep=",")
      tmp$CBMDB_ID <- as.integer(tmp$CBMDB_ID)
      vcf_data_df_cbmdb <- merge(tmp,dplyr::select(vcf_data_df_cbmdb,-c(CBMDB_ID)),by.x = "VAR_ID",by.y = "VAR_ID")
      cbmdb_calls <- dplyr::select(vcf_data_df_cbmdb,dplyr::one_of(tags))
      eitems <- dplyr::left_join(cbmdb_calls,dplyr::filter(dplyr::select(cbmdb_biomarkers,-c(drug_family,transvar_id)),alteration_type == 'MUT'),by=c("CBMDB_ID" = "CBMDB_ID"))
      names(eitems) <- toupper(names(eitems))
      eitems$BIOMARKER_MAPPING <- 'exact'
      clinical_evidence_items <- rbind(clinical_evidence_items, eitems)
    }
  }
  else{
    vcf_data_df_civic <- vcf_data_df %>% dplyr::filter(!is.na(CIVIC_ID_2))
    if(nrow(vcf_data_df_civic) > 0){
      tmp <- dplyr::select(vcf_data_df_civic,CIVIC_ID_2,VAR_ID)
      tmp <- tmp %>% tidyr::separate_rows(CIVIC_ID_2,sep=",")
      vcf_data_df_civic <- merge(tmp,dplyr::select(vcf_data_df_civic,-c(CIVIC_ID_2)),by.x = "VAR_ID",by.y = "VAR_ID")
      civic_calls <- dplyr::select(vcf_data_df_civic,dplyr::one_of(tags))
      clinical_evidence_items <- dplyr::left_join(civic_calls,dplyr::filter(civic_biomarkers,alteration_type == 'MUT'),by=c("CIVIC_ID_2" = "civic_id"))
      clinical_evidence_items <- clinical_evidence_items %>% dplyr::filter(mapping_category != 'gene')
      names(clinical_evidence_items) <- toupper(names(clinical_evidence_items))
      if(nrow(clinical_evidence_items) == 0){
        return(clinical_evidence_items)
      }
      else{
        if(mapping == 'codon' & 'CIVIC_CODON' %in% colnames(clinical_evidence_items)){
          if(nrow(clinical_evidence_items[!is.na(clinical_evidence_items$CIVIC_CODON),]) > 0){
            clinical_evidence_items$CIVIC_CODON <- as.character(clinical_evidence_items$CIVIC_CODON)
            clinical_evidence_items <- dplyr::filter(clinical_evidence_items, !is.na(CIVIC_CODON))
            clinical_evidence_items$CODON <- as.character(stringr::str_replace_all(clinical_evidence_items$PROTEIN_CHANGE,"p\\.[A-Z]{1,}|[A-Z]|[A-Z]{1,}$|fs$|del$|dup$|delins[A-Z]{1,}$",""))
            clinical_evidence_items$CIVIC_CODON <- as.character(clinical_evidence_items$CIVIC_CODON)
            clinical_evidence_items <- clinical_evidence_items %>% dplyr::filter(startsWith(CODON,CIVIC_CODON))
            if(nrow(clinical_evidence_items) > 0){
              clinical_evidence_items <- clinical_evidence_items %>% dplyr::filter(!is.na(CIVIC_CONSEQUENCE) & startsWith(CONSEQUENCE,CIVIC_CONSEQUENCE) | is.na(CIVIC_CONSEQUENCE))
              clinical_evidence_items <- clinical_evidence_items %>% dplyr::select(-c(EXON,CODON,CIVIC_CONSEQUENCE,MAPPING_CATEGORY,CIVIC_CODON,CIVIC_EXON))
              names(clinical_evidence_items) <- toupper(names(clinical_evidence_items))
            }
            else{
              clinical_evidence_items <- data.frame()
            }
          }
          else{
            clinical_evidence_items <- data.frame()
          }
        }
        if(mapping == 'exon' & 'CIVIC_EXON' %in% colnames(clinical_evidence_items)){
          if(nrow(clinical_evidence_items[!is.na(clinical_evidence_items$CIVIC_EXON),]) > 0){
            clinical_evidence_items <- dplyr::filter(clinical_evidence_items, !is.na(CIVIC_EXON))
            clinical_evidence_items$EXON <- as.integer(stringr::str_split_fixed(clinical_evidence_items$EXON,"/",2)[,1])
            clinical_evidence_items <- clinical_evidence_items %>% dplyr::filter(EXON == CIVIC_EXON)
            if(nrow(clinical_evidence_items) > 0){
              clinical_evidence_items <- clinical_evidence_items %>% dplyr::filter(!is.na(CIVIC_CONSEQUENCE) & startsWith(CONSEQUENCE,CIVIC_CONSEQUENCE) | is.na(CIVIC_CONSEQUENCE))
              clinical_evidence_items <- clinical_evidence_items %>% dplyr::select(-c(EXON,CIVIC_CONSEQUENCE,MAPPING_CATEGORY,CIVIC_CODON,CIVIC_EXON))
              names(clinical_evidence_items) <- toupper(names(clinical_evidence_items))
            }
            else{
              clinical_evidence_items <- data.frame()
            }
          }
          else{
            clinical_evidence_items <- data.frame()
          }
        }
      }
    }

  }

  if(nrow(clinical_evidence_items) > 0){
    clinical_evidence_items <- clinical_evidence_items[,c("SYMBOL","PROTEIN_CHANGE","CLINICAL_SIGNIFICANCE","EVIDENCE_LEVEL","EVIDENCE_TYPE","EVIDENCE_DIRECTION","DISEASE_NAME","DESCRIPTION","GENE_NAME","CITATION","DRUG_NAMES","RATING","PROTEIN_DOMAIN","PROTEIN_FEATURE","CDS_CHANGE","CANCER_MUTATION_HOTSPOT",'OTHER_LITERATURE_DOCM','OTHER_DISEASE_DOCM',"INTOGEN_DRIVER_MUT","CONSEQUENCE","EFFECT_PREDICTIONS","VEP_ALL_CONSEQUENCE","DBSNP","COSMIC","COSMIC_SITE_HISTOLOGY","CLINVAR","GENOME_VERSION","GENOMIC_CHANGE","ONCOSCORE","CALL_CONFIDENCE","DP_TUMOR","AF_TUMOR","DP_NORMAL","AF_NORMAL")]
    unique_variants <- clinical_evidence_items %>% dplyr::select(SYMBOL,CONSEQUENCE,PROTEIN_CHANGE,CDS_CHANGE) %>% dplyr::distinct()
    clinical_evidence_items <- clinical_evidence_items %>% dplyr::arrange(EVIDENCE_LEVEL,RATING)
    cat(nrow(clinical_evidence_items),' clinical evidence items found .. (',nrow(unique_variants),' unique variants), mapping = ',mapping,'\n',sep="")
    cat('Underlying variants:','\n')
    for(i in 1:nrow(unique_variants)){
      cat(paste(unique_variants[i,],collapse=" "),'\n')
    }
  }
  else{
    cat(nrow(clinical_evidence_items),' clinical evidence items found .. mapping = ',mapping,'\n')
  }
  return(clinical_evidence_items)

}


#' Function that adds read support (depth, allelic fraction) for tumor and normal
#'
#' @param vcf_data_df
#'
#' @return vcf_data_df
#' @export
#'
add_read_support <- function(vcf_data_df){
  for(v in c('DP_TUMOR','AF_TUMOR','DP_NORMAL','AF_NORMAL')){
    vcf_data_df[v] <- NA
  }
  if(!is.null(vcf_data_df$ADT) && !is.null(vcf_data_df$ADC)){
    tmp_tumor <- dplyr::select(tidyr::separate(vcf_data_df,ADT, c('refn','altn'),convert=T), refn, altn)
    tmp_normal <- dplyr::select(tidyr::separate(vcf_data_df,ADC, c('refn','altn'),convert=T), refn, altn)

    ## for variants called in which alternate and reference allele make up a fraction of total read depth
    tmp_tumor_dp <- as.integer(vcf_data_df$DPT)
    vcf_data_df$DP_TUMOR <- tmp_tumor$refn + tmp_tumor$altn
    tmp2 <- as.data.frame(vcf_data_df %>% dplyr::rowwise() %>% dplyr::mutate(DP_TUMOR_REVISED = max(DP_TUMOR,DPT,na.rm=T)))
    vcf_data_df$DP_TUMOR <- tmp2$DP_TUMOR_REVISED

    vcf_data_df$AF_TUMOR <- round(tmp_tumor$altn / vcf_data_df$DP_TUMOR, digits=4)
    vcf_data_df$DP_NORMAL <- tmp_normal$refn + tmp_normal$altn
    vcf_data_df$AF_NORMAL <- round(tmp_normal$altn / vcf_data_df$DP_NORMAL, digits=4)
  }
  return(vcf_data_df)
}


#' Function that adds SwissProt feature descriptions based on keys coming from OncoVarExplorer
#'
#' @param vcf_data_df
#'
#' @return vcf_data_df
#' @export
#'
add_swissprot_feature_descriptions <- function(vcf_data_df){

  swissprot_features$UNIPROT_FEATURE <- paste(paste(swissprot_features$uniprot_id,swissprot_features$feature_type,sep=":"),paste(swissprot_features$aa_start,swissprot_features$aa_stop,sep="-"),sep=":")
  swissprot_features$PF <- paste(paste(swissprot_features$type_description,paste(swissprot_features$aa_start,swissprot_features$aa_stop,sep="-"),sep=":"),swissprot_features$description,sep=":")

  if("UNIPROT_FEATURE" %in% colnames(vcf_data_df) & "VAR_ID" %in% colnames(vcf_data_df)){
    feature_df <- dplyr::select(vcf_data_df,UNIPROT_FEATURE,VAR_ID) %>% dplyr::distinct()
    if(nrow(feature_df) == 0){
      return(vcf_data_df)
    }
    feature_df <- feature_df %>% tidyr::separate_rows(UNIPROT_FEATURE,sep="&")
    feature_df <- as.data.frame(dplyr::left_join(feature_df,dplyr::select(swissprot_features,UNIPROT_FEATURE,PF),by=c("UNIPROT_FEATURE")))
    feature_df <- as.data.frame(dplyr::group_by(feature_df, VAR_ID) %>% dplyr::summarise(PROTEIN_FEATURE = paste(PF, collapse=", ")))
    feature_df[feature_df$PROTEIN_FEATURE == "NA",]$PROTEIN_FEATURE <- NA

    vcf_data_df <- dplyr::left_join(dplyr::select(vcf_data_df,-UNIPROT_FEATURE),feature_df, by=c("VAR_ID" = "VAR_ID"))
  }

  return(vcf_data_df)

}

#' Function that reads a VCF from OncoVarExplorer pipeline
#'
#' @param vcf_gz_file
#'
#' @return vcf_data_df
#' @export
#'
get_calls <- function(vcf_gz_file){
  vcf_data_vr <- VariantAnnotation::readVcfAsVRanges(vcf_gz_file,genome = "hg19")
  vcf_data_vr <- vcf_data_vr[!is.na(vcf_data_vr$GT) & !(vcf_data_vr$GT == '.'),]
  vcf_data_vr <- OncoVarReporter::postprocess_vranges_info(vcf_data_vr)

  vcf_data_df <- as.data.frame(vcf_data_vr)
  vcf_data_df$GENOME_VERSION <- 'GRCh37'
  vcf_data_df <- dplyr::rename(vcf_data_df, CHROM = seqnames, POS = start, REF = ref, ALT = alt, CONSEQUENCE = Consequence, PROTEIN_CHANGE = HGVSp_short)
  vcf_data_df$GENOMIC_CHANGE <- paste(paste(paste(paste0("g.chr",vcf_data_df$CHROM),vcf_data_df$POS,sep=":"),vcf_data_df$REF,sep=":"),vcf_data_df$ALT,sep=">")

  vcf_data_df <- OncoVarReporter::add_pfam_domain_links(vcf_data_df)
  vcf_data_df <- OncoVarReporter::add_swissprot_feature_descriptions(vcf_data_df)
  vcf_data_df <- OncoVarReporter::add_read_support(vcf_data_df)
  vcf_data_df <- dplyr::left_join(vcf_data_df, docm_literature, by=c("VAR_ID"))

  gencode_xref <- dplyr::rename(gene_xref, Gene = ensembl_gene_id, GENENAME = name, ENTREZ_ID = entrezgene)
  gencode_xref <- gencode_xref %>% dplyr::filter(!is.na(Gene)) %>% dplyr::select(Gene,GENENAME,ENTREZ_ID) %>% dplyr::distinct()
  gencode_xref$GENENAME <- stringr::str_replace(gencode_xref$GENENAME," \\[.{1,}$","")
  gencode_xref$ENTREZ_ID <- as.character(gencode_xref$ENTREZ_ID)
  gencode_xref <- dplyr::filter(gencode_xref, !is.na(GENENAME) & !is.na(ENTREZ_ID))

  vcf_data_df <- dplyr::left_join(vcf_data_df,gencode_xref,by=c("ENTREZ_ID","Gene"))

  tmp <- dplyr::select(clinvar, CLINVAR_TRAITS_ALL, CLINVAR_MSID)
  if ("CLINVAR_MSID" %in% colnames(vcf_data_df)){
    vcf_data_df <- dplyr::left_join(vcf_data_df,tmp,by=c("CLINVAR_MSID"))
  }
  if("COSMIC_SITE_HISTOLOGY" %in% colnames(vcf_data_df)){
    vcf_data_df$COSMIC_SITE_HISTOLOGY <- stringr::str_replace_all(vcf_data_df$COSMIC_SITE_HISTOLOGY,"&",", ")
  }
  if("EFFECT_PREDICTIONS" %in% colnames(vcf_data_df)){
    vcf_data_df$EFFECT_PREDICTIONS <- stringr::str_replace_all(vcf_data_df$EFFECT_PREDICTIONS,"\\.&|\\.$","NA&")
    vcf_data_df$EFFECT_PREDICTIONS <- stringr::str_replace_all(vcf_data_df$EFFECT_PREDICTIONS,"&$","")
    vcf_data_df$EFFECT_PREDICTIONS <- stringr::str_replace_all(vcf_data_df$EFFECT_PREDICTIONS,"&",", ")
  }
  if("INTOGEN_DRIVER_MUT" %in% colnames(vcf_data_df)){
    vcf_data_df$INTOGEN_DRIVER_MUT <- stringr::str_replace_all(vcf_data_df$INTOGEN_DRIVER_MUT,"&",", ")
  }
  if("CANCER_CENSUS_SOMATIC" %in% colnames(vcf_data_df)){
    vcf_data_df$CANCER_CENSUS_SOMATIC <- stringr::str_replace_all(vcf_data_df$CANCER_CENSUS_SOMATIC,"&",", ")
  }
  if("COSMIC_DRUG_RESISTANCE" %in% colnames(vcf_data_df)){
    vcf_data_df$COSMIC_DRUG_RESISTANCE <- stringr::str_replace_all(vcf_data_df$COSMIC_DRUG_RESISTANCE,"&",", ")
  }
  if("VEP_ALL_CONSEQUENCE" %in% colnames(vcf_data_df)){
    vcf_data_df$VEP_ALL_CONSEQUENCE <- stringr::str_replace_all(vcf_data_df$VEP_ALL_CONSEQUENCE,",",", ")
  }
  if("DOCM_DISEASE" %in% colnames(vcf_data_df)){
    vcf_data_df$DOCM_DISEASE <- stringr::str_replace_all(vcf_data_df$DOCM_DISEASE,",",", ")
  }

  ## Add HTML links for COSMIC, DBSNP, DGIDB and CLINVAR entries
  if(!("COSMIC" %in% colnames(vcf_data_df))){
    vcf_data_df <- OncoVarReporter::annotate_variant_link(vcf_data_df, vardb = 'COSMIC')
    vcf_data_df <- dplyr::rename(vcf_data_df, COSMIC = COSMICLINK)
  }
  if(!("DBSNP" %in% colnames(vcf_data_df))){
    vcf_data_df <- OncoVarReporter::annotate_variant_link(vcf_data_df, vardb = 'DBSNP')
    vcf_data_df <- dplyr::rename(vcf_data_df, DBSNP = DBSNPLINK)
  }
  #if(!("GENOMIC_CHANGE" %in% colnames(vcf_data_df))){
  #vcf_data_df <- OncoVarReporter::annotate_variant_link(vcf_data_df, vardb = 'UCSC')
  #vcf_data_df <- dplyr::rename(vcf_data_df, GENOMIC_CHANGE = UCSC_LINK)

  if(!("CLINVAR" %in% colnames(vcf_data_df))){
    vcf_data_df <- OncoVarReporter::annotate_variant_link(vcf_data_df, vardb = 'CLINVAR')
    vcf_data_df <- dplyr::rename(vcf_data_df, CLINVAR = CLINVARLINK)
  }
  if(!("ANTINEOPLASTIC_DRUG_INTERACTIONS" %in% colnames(vcf_data_df))){
    vcf_data_df <- OncoVarReporter::annotate_variant_link(vcf_data_df, vardb = 'DGIDB')
    vcf_data_df <- dplyr::rename(vcf_data_df, ANTINEOPLASTIC_DRUG_INTERACTIONS = DGIDBLINK)
  }
  if(!("GENE_NAME" %in% colnames(vcf_data_df))){
    vcf_data_df <- OncoVarReporter::annotate_variant_link(vcf_data_df, vardb = 'NCBI_GENE')
    vcf_data_df <- dplyr::rename(vcf_data_df, GENE_NAME = NCBI_GENE_LINK)
  }

  return(vcf_data_df)

}

#' Function that filters a data frame with variants according to population-specific germline frequencies
#'
#' @param var_df data frame with variants
#' @param pop population ('european' or 'global')
#' @param dbquery '1KG' or 'ExAC'
#' @param min_af minimum allele frequency required for variant to be filtered
#'
#' @return var_df
#' @export
#'

filter_db_germline_variants <- function(var_df, pop='european',dbquery = '1KG', min_af = 0.05){

  if(pop == 'nor'){
    if(any(grepl(paste0("^AF_NOR$"),names(var_df)))){
      var_df <- var_df %>% dplyr::filter(is.na(AF_NOR) | AF_NOR < min_af)
    }
  }
  else{
    if(pop == 'european' | pop == 'global'){
      if(dbquery == '1KG'){
        if(any(grepl(paste0("^EUR_AF_1KG$"),names(var_df)))){
          var_df <- var_df %>% dplyr::filter(is.na(EUR_AF_1KG) | EUR_AF_1KG < min_af)
        }
      }
      else{
        if(any(grepl(paste0("^NFE_AF_EXAC$"),names(var_df)))){
          var_df <- var_df %>% dplyr::filter(is.na(NFE_AF_EXAC) | NFE_AF_EXAC < min_af)
        }
      }
    }
    if(pop == 'global'){
      if(dbquery == '1KG'){
        pop_tags <- c('EAS_AF_1KG','AMR_AF_1KG','AFR_AF_1KG','SAS_AF_1KG')
        for(poptag in pop_tags){
          if(any(grepl(paste0("^",poptag,"$"),names(var_df)))){
            var_df <- var_df[is.na(var_df[poptag]) | var_df[poptag] < min_af,]
          }
        }
      }
      else{
        pop_tags <- c('SAS_AF_EXAC','EAS_AF_EXAC','AMR_AF_EXAC','AFR_AF_EXAC','FIN_AF_EXAC','OTH_AF_EXAC')
        for(poptag in pop_tags){
          if(any(grepl(paste0("^",poptag,"$"),names(var_df)))){
            var_df <- var_df[is.na(var_df[poptag]) | var_df[poptag] < min_af,]
          }
        }
      }
    }
  }
  return(var_df)
}


#' Function that assigns a category ('Rare','Common' etc) to population-specific germline frequencies
#'
#' @param var_df data frame with variants
#' @param pop_af_column population_column
#'
#' @return var_df
#' @export
#'
assign_poplevel_frequency_class <- function(var_df, pop_af_column){

  if(any(grepl(paste0("^",pop_af_column,"$"),names(var_df)))){
    if(nrow(var_df[!is.na(var_df[pop_af_column]) & var_df[pop_af_column] >= 0.05,]) > 0){
      var_df[(!is.na(var_df[pop_af_column]) & var_df[pop_af_column] >= 0.05),]$pop_common <- 'Common'
    }
    if(nrow(var_df[(!is.na(var_df[pop_af_column]) & var_df[pop_af_column] >= 0.01 & var_df[pop_af_column] < 0.05),]) > 0){
      var_df[(!is.na(var_df[pop_af_column]) & var_df[pop_af_column] >= 0.01 & var_df[pop_af_column] < 0.05),]$pop_lowfreq <- 'LowFreq'
    }
    if(nrow(var_df[(!is.na(var_df[pop_af_column]) & var_df[pop_af_column] >= 0.001 & var_df[pop_af_column] < 0.01),]) > 0){
      var_df[(!is.na(var_df[pop_af_column]) & var_df[pop_af_column] >= 0.001 & var_df[pop_af_column] < 0.01),]$pop_rare <- 'Rare'
    }
    if(nrow(var_df[!is.na(var_df[pop_af_column]) & var_df[pop_af_column] < 0.001 & var_df[pop_af_column] > 0,]) > 0){
      var_df[(!is.na(var_df[pop_af_column]) & var_df[pop_af_column] < 0.001 & var_df[pop_af_column] > 0),]$pop_veryrare <- 'VeryRare'
    }
    if(nrow(var_df[!is.na(var_df[pop_af_column]) & var_df[pop_af_column] == 0.00,]) > 0){
      var_df[(!is.na(var_df[pop_af_column]) & var_df[pop_af_column] == 0.00),]$pop_monomorphic <- 'Monomorphic'
    }
  }
  return(var_df)

}

#' Function that assigns a category ('Rare','Common' etc) to population-specific germline frequencies
#'
#' @param var_df data frame with variants
#' @param dbquery 1KG or ExAC
#' @param pop population
#' @param result_tag name of result column
#'
#' @return var_df
#' @export
#'
assign_poplevel_frequency <- function(var_df, dbquery='1KG', pop='european', result_tag = 'FREQ_EXAC_EUROPEAN'){

  var_df$pop_monomorphic <- rep(NA,nrow(var_df))
  var_df$pop_common <- rep(NA,nrow(var_df))
  var_df$pop_rare <- rep(NA,nrow(var_df))
  var_df$pop_veryrare <- rep(NA,nrow(var_df))
  var_df$pop_lowfreq <- rep(NA,nrow(var_df))

  pop_db <- data.frame('population' = 'american','db' = '1KG', 'tag' = 'AMR_AF_1KG', stringsAsFactors = F)
  pop_db <- rbind(pop_db, data.frame('population' = 'african', 'db' = '1KG', 'tag' = 'AFR_AF_1KG'))
  pop_db <- rbind(pop_db, data.frame('population' = 'european', 'db' = '1KG', 'tag' = 'EUR_AF_1KG'))
  pop_db <- rbind(pop_db, data.frame('population' = 'east_asian', 'db' = '1KG', 'tag' = 'EAS_AF_1KG'))
  pop_db <- rbind(pop_db, data.frame('population' = 'south_asian', 'db' = '1KG','tag' = 'SAS_AF_1KG'))
  pop_db <- rbind(pop_db, data.frame('population' = 'african', 'db' = 'ExAC', 'tag' = 'AFR_AF_EXAC'))
  pop_db <- rbind(pop_db, data.frame('population' = 'american', 'db' = 'ExAC', 'tag' = 'AMR_AF_EXAC'))
  pop_db <- rbind(pop_db, data.frame('population' = 'european', 'db' = 'ExAC', 'tag' = 'NFE_AF_EXAC'))
  pop_db <- rbind(pop_db, data.frame('population' = 'east_asian', 'db' = 'ExAC', 'tag' = 'EAS_AF_EXAC'))
  pop_db <- rbind(pop_db, data.frame('population' = 'south_asian', 'db' = 'ExAC','tag' = 'SAS_AF_EXAC'))
  pop_db <- rbind(pop_db, data.frame('population' = 'global', 'db' = 'ExAC','tag' = 'GLOBAL_AF_EXAC'))

  tags <- character()
  tags <- c(dplyr::filter(pop_db, population == pop & db == dbquery)$tag)
  for(i in 1:length(tags)){
    var_df <- assign_poplevel_frequency_class(var_df, tags[i])
  }
  var_df[result_tag] <- stringr::str_replace_all(stringr::str_replace_all(paste(var_df$pop_monomorphic,var_df$pop_rare,var_df$pop_veryrare,var_df$pop_lowfreq,var_df$pop_common,sep=","), "(,{0,}NA(,){0,}){1,}",","),"(^,)|(,$)","")
  var_df[nchar(var_df[,result_tag]) == 0,result_tag] <- NA

  var_df <- dplyr::select(var_df, -pop_monomorphic, -pop_common, -pop_rare, -pop_veryrare, -pop_lowfreq)
  return(var_df)
}
