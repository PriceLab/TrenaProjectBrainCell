#---------------------------------------------------------------------------------------------------
#' @import methods
#' @import TrenaProject
#' @importFrom AnnotationDbi select
#' @import org.Hs.eg.db
#'
#' @title TrenaProjectLiver-class
#'
#' @name TrenaProjectLiver-class
#' @rdname TrenaProjectLiver-class
#' @aliases TrenaProjectLiver
#' @exportClass TrenaProjectLiver
#'

.TrenaProjectBrainCell <- setClass("TrenaProjectBrainCell",
                                  contains="TrenaProjectHG38")

#----------------------------------------------------------------------------------------------------
#' Define an object of class TrenaProjectLiver
#'
#' @description
#' Expression, variant and covariate data for the genes of interest (perhaps unbounded) for pre-term birth studies
#'
#' @rdname TrenaProjectLiver-class
#'
#' @export
#'
#' @return An object of the TrenaProjectLiver class
#'

TrenaProjectBrainCell <- function(quiet=TRUE)

{
   genomeName <- "hg38"

   directory <- system.file(package="TrenaProjectBrainCell", "extdata", "geneSets")
   geneSet.files <- list.files(directory)
   geneSets <- list()

   for(file in geneSet.files){
      full.path <- file.path(directory, file)
      genes <- scan(full.path, sep="\t", what=character(0), quiet=TRUE, comment.char="#")
      geneSet.name <- sub(".txt", "", file)
      geneSets[[geneSet.name]] <- genes
      }

   footprintDatabaseNames <- c("brain_hint_16",  "brain_hint_20", "brain_wellington_16", "brain_wellington_20")
   expressionDirectory <- system.file(package="TrenaProjectBrainCell", "extdata", "expression")
   variantsDirectory <- system.file(package="TrenaProjectBrainCell", "extdata", "variants")
   footprintDatabaseHost <- "khaleesi.systemsbiology.net"
   footprintDatabasePort <- 5432

   covariatesFile <- NA_character_;

   stopifnot(file.exists(expressionDirectory))

   .TrenaProjectBrainCell(TrenaProjectHG38(projectName="TrenaProjectBrainCell",
                                       supportedGenes=geneSets[[1]],
                                       footprintDatabaseHost=footprintDatabaseHost,
                                       footprintDatabasePort=footprintDatabasePort,
                                       footprintDatabaseNames=footprintDatabaseNames,
                                       expressionDirectory=expressionDirectory,
                                       variantsDirectory=variantsDirectory,
                                       covariatesFile=covariatesFile,
                                       quiet=quiet
                                       ))

} # TrenaProjectBrainCell, the constructor
#----------------------------------------------------------------------------------------------------
