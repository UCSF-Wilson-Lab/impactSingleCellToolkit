#' @title Generate input files for GLIPH2
#' @description Takes an input TCR VDJ data and generates all files needed to run GLIPH2
#' @param tcr.df TCR VDJ data outputted from the Immcantation pipeline
#' @param patient.col Column name which contains the patient ID <default: 'patient'>
#' @param clone.col Column name which contains the clone ID <default: 'clone_id'>
#' @param locus.col Column name which contains the type of BCR or TCR chain <default: 'locus'>
#' @param unique.id.col Column name which contains unique cell IDs <default: 'unique_cell_id'>
#' @param vcall.col Column name which contains V gene calls <default: 'v_call'>
#' @param jcall.col Column name which contains J gene calls <default: 'j_call'>
#' @param cdr3.aa.col Column name which contains the AA sequences of the junction including the CDR3 region <default: 'junction_aa'>
#' @param results.dir Path to results directory <default: './FilesForGLIPH2'>
#' @return No variables returned, but a set of files will be outputted to the results directory provided.
#' @references 
#' Huang H., Wang C., et al. \emph{Nat Biotechnol} 2020 Apr 27; https://doi.org/10.1038/s41587-020-0505-4
#' @examples generateGLIPH2inputFiles(tcr.df = tcr_df,results.dir = "~/FilesForGLIPH2")
#' @import data.table
#' @import stringr
#' @import dplyr
#' @export

generateGLIPH2inputFiles <- function(tcr.df,
                                     patient.col = "patient",
                                     clone.col   = "clone_id",
                                     locus.col   = "locus",
                                     unique.id.col = "unique_cell_id",
                                     vcall.col     = "v_call",
                                     jcall.col     = "j_call",
                                     cdr3.aa.col   = "junction_aa",
                                     results.dir = "./FilesForGLIPH2") {
  # List output files
  dir.create(results.dir,recursive = T)
  output_tsv   <- file.path(results.dir,"TCR_input_table.GLIPH2.tsv")
  output_cdr3b <- file.path(results.dir,"ref.txt")
  output_trb   <- file.path(results.dir,"ref_V.txt")
  output_len   <- file.path(results.dir,"ref_L.txt")
  
  # Add a column for patient:clonotype
  tcr.df$pt_clonotype <- paste(tcr.df[,patient.col],tcr.df[,clone.col],sep = ":")
  
  ### 2. Create empty GLIPH data frame
  pt_clonotype_list <- unique(tcr.df$pt_clonotype)

  # Column names: CDR3b		TRBV	TRBJ	CDR3a		TRAV		TRAJ	PatientCounts
  input_df <- data.frame(matrix(nrow = length(pt_clonotype_list),ncol = 6))
  names(input_df) <- c("CDR3b","TRBV","TRBJ","CDR3a","subject:condition","count")
  
  row.names(input_df) <- pt_clonotype_list
  
  ### 3. Fill Input data frame with Alpha and Beta info from Immcantation Results 
  for (clone in pt_clonotype_list) {
    clone_df   <- tcr.df[tcr.df$pt_clonotype %in% clone,]
    chains     <- unique(clone_df[,locus.col])
    cell_count <- length(unique(clone_df[,unique.id.col]))
    
    if(length(chains) != 2){next} # At this stage every TCR should have an alpha and beta chain
    
    tra_row <- clone_df[clone_df[,locus.col] %in% c("TRA"),]
    trb_row <- clone_df[clone_df[,locus.col] %in% c("TRB"),]
    
    # Populate input DF
    input_row <- input_df[clone,]
    
    # Beta
    input_row$CDR3b <- unique(trb_row[,cdr3.aa.col])
    vgene           <- unique(tstrsplit(as.character(trb_row[,vcall.col]),"\\*")[[1]])
    jgene           <- unique(tstrsplit(as.character(trb_row[,jcall.col]),"\\*")[[1]])
    input_row$TRBV  <- vgene
    input_row$TRBJ  <- jgene
    
    # Alpha
    input_row$CDR3a <- unique(tra_row[,cdr3.aa.col])
    
    # Patient ID ('Patient' column)
    patient   <- unique(clone_df[,patient.col])
    input_row$`subject:condition` <- patient
    
    # Cell count
    input_row$count <- cell_count
    
    # Add row back to original input_df
    input_df[row.names(input_df) %in% clone,] <- input_row
  }
  
  
  ### 4. Write TSV for GLIPH
  write.table(input_df,file = output_tsv,row.names = F,col.names = F,quote = F,sep = "\t")
  
  ### 5. Create CDR3b files for GLIPH2
  cdr3b_df <- input_df[,c("CDR3b","TRBV","TRBJ")]
  write.table(cdr3b_df,file = output_cdr3b,row.names = F,col.names = F,quote = F,sep = "\t")
  
  ### 6. TRBV freq
  trb_df        <- as.data.frame(table(cdr3b_df$TRBV))
  names(trb_df) <- c("vgene","Freq")
  total_v_usage <- sum(trb_df$Freq)
  trb_df$Freq   <- trb_df$Freq /total_v_usage
  write.table(trb_df,file = output_trb,row.names = F,col.names = F,quote = F,sep = "\t")
  
  ### 7. CDR3b lengths freq
  cdr3b_len_df <- as.data.frame(cdr3b_df[,c("CDR3b")])
  colnames(cdr3b_len_df) <- c("CDR3b")
  cdr3b_len_df$AA_len <- nchar(as.character(cdr3b_len_df$CDR3b))
  
  len_freq_df <- as.data.frame(table(cdr3b_len_df$AA_len))
  names(len_freq_df) <- c("AA_length","Freq")
  total       <- sum(len_freq_df$Freq)
  len_freq_df$Freq   <- len_freq_df$Freq /total
  write.table(len_freq_df,file = output_len,row.names = F,col.names = F,quote = F,sep = "\t")
  
}