#' OmicSelector_download_tissue_miRNA_data_from_TCGA
#'
#' Function allowing to dowload all miRNA-seq data from TCGA.
#'
#' @param data_folder Path. Where to dowload the data?
#'
#' @export
OmicSelector_download_tissue_miRNA_data_from_TCGA = function(data_folder = getwd())
{
  oldwd = getwd()
  setwd(data_folder)
  suppressMessages(library(GDCRNATools))
  suppressMessages(library(data.table))

  zacznij_od_projektu = 1

  projects = TCGAbiolinks::getGDCprojects()
  str(projects$project_id)
  fwrite(projects, "projects.csv")

  for (i in zacznij_od_projektu:length(projects$project_id)) {
    project = projects$project_id[i]
    print(project)

    #rnadir <- paste('RNAseq',project, sep='/')
    mirdir <- paste('miRNAs',project, sep='/')
    clinicaldir <- paste('Clinical',project, sep='/')

    result = tryCatch({
      ####### Download mature miRNA data #######
      gdcRNADownload(project.id     = project,
                     data.type      = 'miRNAs',
                     write.manifest = FALSE,
                     method         = 'gdc-client',
                     directory      = mirdir)
      ####### Download clinical data #######

      gdcClinicalDownload(project.id     = project,
                          write.manifest = FALSE,
                          method         = 'gdc-client',
                          directory      = clinicaldir)

      ####### Parse miRNAs metadata #######
      metaMatrix.MIR <- gdcParseMetadata(project.id = project,
                                         data.type  = 'miRNAs',
                                         write.meta = FALSE)

      ####### Filter duplicated samples in miRNAs metadata #######
      metaMatrix.MIR <- gdcFilterDuplicate(metaMatrix.MIR)

      ####### Filter non-Primary Tumor and non-Solid Tissue Normal samples in miRNAs metadata #######
      metaMatrix.MIR <- gdcFilterSampleType(metaMatrix.MIR)

      mirCounts <- gdcRNAMerge(metadata  = metaMatrix.MIR,
                               path      = mirdir, # the folder in which the data stored
                               organized = FALSE, # if the data are in separate folders
                               data.type = 'miRNAs')
      clinicalDa <- gdcClinicalMerge(path = clinicaldir, key.info = TRUE)
      clinicalDa$patient = row.names(clinicalDa)

      temp = t(mirCounts)
      temp2 = cbind(metaMatrix.MIR,temp)
      temp3 = merge(x=temp2,y=clinicalDa,by="patient", all.x=TRUE)
      fwrite(temp3, paste0("miRNA_",project,".csv"))

    }, warning = function(w) {
      cat("\n\n\nDownloaded with warnings:\n")
      print(w)
    }, error = function(e) {
      cat("\n\n\nError!!!")
      print(e)
    }, finally = {
      cat("\n\n\nGoing to the next project...")
    })
  }

  setwd(data_folder)
}
