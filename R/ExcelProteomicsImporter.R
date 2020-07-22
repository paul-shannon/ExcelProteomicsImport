#' @import org.Hs.eg.db
#' @importFrom methods new
#'
#' @title ExcelProteomicsImporter
#------------------------------------------------------------------------------------------------------------------------
#' @name ExcelProteomicsImporter-class
#' @rdname ExcelProteomicsImporter-class
#' @aliases ExcelProteomicsImporter
#'
#' @import methods

.ExcelProteomicsImporter <- setClass("ExcelProteomicsImporter",
                          representation = representation(
                             rawTable="data.frame",
                             state="environment"
                             )
                          )
#------------------------------------------------------------------------------------------------------------------------
setGeneric('getRawTable', signature='obj', function(obj) standardGeneric ('getRawTable'))
#------------------------------------------------------------------------------------------------------------------------
#' Create a ExcelProteomicsImporter connection
#'
#' @rdname ExcelProteomicsImporter-class
#'
#' @param tabDelimitedFilename character
#'
#' @return An object of the ExcelProteomicsImporter class
#'
#' @export
#'
ExcelProteomicsImporter <- function(tabDelimitedFilename)
{
   state <- new.env(parent=emptyenv())
   stopifnot(file.exists(tabDelimitedFilename))
   tbl.raw <- read.table(tabDelimitedFilename, sep="\t", header=TRUE, as.is=TRUE)
   .ExcelProteomicsImporter(rawTable=tbl.raw, state=state)

} # ctor
#------------------------------------------------------------------------------------------------------------------------
#' get the original, raw data.frame
#'
#' @rdname getRawTable
#' @aliases getRawTable
#'
#' @export
#'
setMethod('getRawTable', 'ExcelProteomicsImporter',

     function(obj){
        obj@rawTable
        })

#------------------------------------------------------------------------------------------------------------------------


