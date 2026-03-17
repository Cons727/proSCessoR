#' Run the VDJ Pipeline
#'
#' This function acts as a wrapper to execute the bash-based pairing pipeline, handling module loads, fasta merging, and IgBLAST alignment.
#'
#' @param outname Character. The output prefix for the run (e.g., "HVTN301_BCS2040").
#' @param ref Character. Absolute path to the reference directory.
#' @param db_name Character. Database name to use. Defaults to "ogrdb".
#' @param vdj Character. VDJ type. Defaults to "vdj_b".
#' @param fasta Character. Path to the fasta file. Defaults to Processing_QC/${outname}_singlets_prod_contig.fasta.
#' @param igm_igd Character. Flag for IGM/IGD contig selection. Defaults to "FALSE".
#' @param work_dir Character. The working directory where the raw data folders (e.g., P01) are located. Defaults to current working directory.
#' @param python_site_package Character. Path to python libraries where all python dependancies are installed for example "/home/mikey/.local/lib/python3.8/site-packages".
#' 
#' @return Invisible integer of the system call exit status (0 means success).
#' @export

pair_and_merge <- function(outname,
                        ref,
                        db_name = "ogrdb",
                        vdj = "vdj_b",
                        fasta = NULL,
                        igm_igd = "FALSE",
                        work_dir = getwd(),
                        python_site_packages) {
  
  # Locate the bash script in the installed package's bin/ directory
  pkg_bin_dir <- system.file("bin", package = "proSCessoR")
  
  if (pkg_bin_dir == "") {
    stop("Could not find bin/ directory.")
  }
  
  # Point to the specific wrapper script
  script_path <- file.path(pkg_bin_dir, "sbatch_pair_merge.sh")
  
  old_dir <- getwd()
  on.exit(setwd(old_dir))
  
  setwd(work_dir)
  
  # Get the current PATH and prepend your package's bin directory
  current_path <- Sys.getenv("PATH")
  new_path <- paste0(pkg_bin_dir, ":", current_path)
  
  # Get default fasta path
  if(is.null(fasta)){
    fasta <- paste0("Processing_QC/", outname, "_singlets_prod_contig.fasta")
  }

  # Export all variables, including the updated PATH
  Sys.setenv(OUTNAME = outname,
             REF = ref,
             DB_NAME = db_name,
             VDJ = vdj,
             FASTA = fasta,
             IGM_IGD = igm_igd,
             USER_PYTHON_PATH = python_site_packages,
             PATH = new_path) # The scripts should now find each other natively
  
  message("Starting pairings pipeline for: ", outname)
  
  status <- system2("bash", args = c("-l", script_path))
  
  if (status != 0) {
    stop("Pipeline failed with exit status: ", status)
  }
  
  message("Pipeline completed successfully.")
  return(invisible(status))
}