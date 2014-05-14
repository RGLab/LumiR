#' blum class
#' 
#' The class for bead level information for a Luminex xMap experiment. It 
#' contains the fluorescence intensity for each bead as well as design 
#' information regarding the experiment.
#' 
#' @slot phenoData An \code{AnnotatedDataFrame}. Contains the information 
#'  regarding the samples (e.g: sample_type, sample_name, well, filename, ...).
#' @slot featureData An \code{AnnotatedDataFrame}. Contains the information 
#'  regarding the analytes: ID and name.
#' @slot exprs A \code{data.table}. Contains the intensities measured for each 
#'  bead.
#' 
#' 
#' @seealso \code{read.experiment}, \code{\link{ExpressionSet}}
#' @author Renan Sauteraud
#'
#' @aliases
#' blum
#' blum-class
#' @importClassesFrom Biobase AnnotatedDataFrame
#' @import data.table
#' @export
#' 
setClass("blum",
         representation=representation( #cannot be ExpressionSet because exprs is not a matrix
           ## Contains information about samples           
           phenoData="AnnotatedDataFrame",
           ## contains information about analytes           
           featureData="AnnotatedDataFrame",
           ## list of bead level data
           ## Stored as samples -> analytes
           exprs="data.table")
         )

#' slum class
#' 
#' The class contains the summarized information for a Luminex xMap experiment.
#' It contains the median fluorescence intensity for each well and analytes as
#' well as information regarding the desing of the experiment.
#' 
#' @slot phenoData A \code{AnnotatedDataFrame}.Contains the information 
#'  regarding the samples (e.g: sample_type, sample_name, well, filename, ...).
#' @slot featureData An \code{AnnotatedDataFrame}. Contains the information 
#'  regarding the analytes: ID and name.
#' @slot assayData An \code{environment} with two identically dimensioned 
#'  matrices with analytes as rownames and sample_id as colnames.
#'  \itemize{
#'    \item \code{exprs}: The \code{matrix} of the MFI.
#'    \item \code{concentration}: The \code{matrix} of the concentrations 
#'      calculated for the MFIs in \code{exprs}.
#'  }
#' @slot fit A \code{data.frame}. Contains the information regarding the 
#'  standard curve fitting: The values of the parameters to use in \code{inv} 
#'  for each couple analyte / sample_id.
#' @slot formula A \code{formula}. The formula used for the standard curve 
#'  fitting. By default a 5-Parameters Logistic.
#' @slot inv A \code{function}. The inverse function of the formula.
#' @slot unit A \code{character}, "MFI".
#' @slot protocolData An optional \code{AnnotatedDataFrame}, inherited from 
#'  \code{ExpressionSet}.
#' @slot experimentData An optional \code{MIAME}, inherited from 
#'  \code{ExpressionSet}.
#' 
#' @seealso \code{slummarize}, \code{\link{ExpressionSet}}, \code{blum}
#' @author Renan Sauteraud
#' 
#' @aliases
#' slum
#' slum-class
#' @importClassesFrom Biobase ExpressionSet
#' @export
#' 
setClass("slum", 
  contains="ExpressionSet", 
  representation(unit="character", formula="formula", inv="function", fit="data.frame"))


setClassUnion("blumORslum", members=c("blum", "slum"))
