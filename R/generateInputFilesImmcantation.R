#' @title Generate input files for the Immcantation Pipeline
#' @description Takes an input directory of 10x VDJ results and generates a combined output annotation CSV and FASTA file of all VDJ contigs (BCR or TCR)
#' @param vdj.results.dir Directory containing VDJ results across all samples
#' @param contig.annotation.file String of the file name of the contig annotation results files <default: 'filtered_contig_annotations.csv'>
#' @param contig.fasta.file String of the file name of the contig FASTA results files <default: 'filtered_contig.fasta'>
#' @param output.file.prefix String of a custom prefix that would be added to the beginning of all files outputted from this function <default: 'combined'>
#' @param contig.id.col Column name within the contig annotations containing all contig IDs <default: 'contig_id'>
#' @param output.dir Output directory for output files
#' @return A combined output annotation CSV and FASTA file of all VDJ contigs.
#' @examples generateInputFilesImmcantation(vdj.results.dir,output.dir,output.file.prefix = "all_BCR")
#' @import Biostrings
#' @export
generateInputFilesImmcantation <- function(vdj.results.dir,
                                           contig.annotation.file = "filtered_contig_annotations.csv",
                                           contig.fasta.file = "filtered_contig.fasta",
                                           output.file.prefix = "combined",
                                           contig.id.col = "contig_id",
                                           output.dir) {
  sample_list <- as.list(list.files(vdj.results.dir))
  annot_file_list <- file.path(list.files(vdj.results.dir,full.names = T),contig.annotation.file)
  fasta_file_list <- file.path(list.files(vdj.results.dir,full.names = T),contig.fasta.file)
  
  # Create combined annotation CSV
  annot_df_list <- lapply(sample_list, formatVDJannotation,
                          vdj.results.dir=vdj.results.dir,
                          contig.annotation.file=contig.annotation.file,
                          contig.id.col=contig.id.col)
  
  annot_df <- do.call("rbind",annot_df_list)
  
  # Create combined FASTA
  fasta_list <- lapply(sample_list, formatVDJfasta,
                       vdj.results.dir=vdj.results.dir,
                       contig.fasta.file=contig.fasta.file)
  fasta_obj <- Reduce(c, fasta_list)
  
  # Output files
  fasta_file_name <- paste(output.file.prefix,contig.fasta.file,sep = "_")
  annot_file_name <- paste(output.file.prefix,contig.annotation.file,sep = "_")
  output_fasta_fh <- file.path(output.dir,fasta_file_name)
  output_annot_fh <- file.path(output.dir,annot_file_name)
  
  writeXStringSet(fasta_obj, output_fasta_fh)
  write.csv(annot_df,output_annot_fh,quote = F,row.names = F)
}


formatVDJannotation <- function(sample,vdj.results.dir,contig.annotation.file,contig.id.col) {
  annot_fh <- file.path(vdj.results.dir,sample,contig.annotation.file)
  annot_df <- read.csv(annot_fh,stringsAsFactors = F)
  annot_df[,contig.id.col] <- paste(sample,annot_df[,contig.id.col],sep = ":")
  
  return(annot_df)
}


formatVDJfasta <- function(sample,vdj.results.dir,contig.fasta.file) {
  fasta_fh <- file.path(vdj.results.dir,sample,contig.fasta.file)
  fasta <- readDNAStringSet(fasta_fh)
  names(fasta) <- paste(sample,names(fasta),sep = ":")
  
  return(fasta)
}
