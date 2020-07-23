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
setGeneric('validateColumnNames', signature='obj', function(obj) standardGeneric ('validateColumnNames'))
setGeneric('getGroupCount', signature='obj', function(obj) standardGeneric ('getGroupCount'))
setGeneric('identifyGroupBoundaries', signature='obj', function(obj) standardGeneric ('identifyGroupBoundaries'))
setGeneric('getGroupBoundaries', signature='obj', function(obj) standardGeneric ('getGroupBoundaries'))
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
#' learn experiment conditions and timepoints from the column names
#'
#' @rdname validateColumnNames
#' @aliases validateColumnNames
#'
#' @export
#'
setMethod('validateColumnNames', 'ExcelProteomicsImporter',

     function(obj){
        tbl <- obj@rawTable
        colnames <- colnames(tbl)
        stopifnot("PeptideLabel" %in% colnames | "Peptide.Label" %in% colnames)
        if("PeptideLabel" %in% colnames)
            analytes <- tbl$PeptideLabel
        if("Peptide.Label" %in% colnames)
            analytes <- tbl[, "Peptide.Label"]

        stopifnot(length(unique(analytes)) == nrow(tbl))

        stopifnot(any(grepl("avg", colnames)))
        stopifnot(any(grepl("stdev", colnames)))
        avg.cols <- grep("avg", colnames)
        sd.cols  <- grep("stdev", colnames)
        stopifnot(length(avg.cols) == length(sd.cols))

        groupCount <- length(avg.cols)   # expect 1 or 2
        stopifnot(groupCount %in% c(1,2))
        obj@state$groupCount <- groupCount
        return(TRUE)
        })

#------------------------------------------------------------------------------------------------------------------------
#' identifyGroupBoundaries: where avg and sd columns begin and end
#'
#' @rdname identifyGroupBoundaries
#' @aliases identifyGroupBoundaries
#'
#' @export
#'
setMethod('identifyGroupBoundaries', 'ExcelProteomicsImporter',

     function(obj){
        colnames <- colnames(obj@state$rawTable)
        avg.cols <- grep("avg", colnames)
        sd.cols  <- grep("stdev", colnames)
        if(obj@state$groupCount == 1){
           obj@stage$group.1.avg <- c(avg.cols[1] + 1, sd.cols[1] - 1)
           obj@stage$group.1.sd  <- c(sd.cols[1] + 1, length(colnames))
           obj@stage$group.2.avg <- NA
           obj@stage$group.2.sd <- NA
           }
        if(obj@state$groupCount == 2){
           obj@stage$group.1.avg <- c(avg.cols[1] + 1, avg.cols[2] - 1)
           obj@stage$group.1.sd  <- c(sd.cols[1] + 1,  sd.cols[2] - 1)
           obj@stage$group.2.avg <- c(avg.cols[2] + 1, sd.cols[1] - 1)
           obj@stage$group.2.sd <-  c(sd.cols[2] + 1,  length(colnames))
           }
        }) # identifyGroupBoundaries

#------------------------------------------------------------------------------------------------------------------------
#' get groupCount, the number of logically separate condition groups in the spreadsheet
#'
#' @rdname getGroupCount
#' @aliases getGroupCount
#'
#' @export
#'
setMethod('getGroupCount', 'ExcelProteomicsImporter',

     function(obj){
        obj@stage$groupCount
        })

#------------------------------------------------------------------------------------------------------------------------
#' get groupGroupBoundaries, the start and stop of columns for avg group/s, sd group/s
#'
#' @rdname getGroupBoundaries
#' @aliases getGroupBoundaries
#'
#' @export
#'
setMethod('getGroupBoundaries', 'ExcelProteomicsImporter',

     function(obj){
         list(avg.1=obj@stage$group.1.avg,
              avg.2=obj@stage$group.2.avg,
              sd.1 =obj@stage$group.1.sd,
              sd.2 =obj@stage$group.2.sd)
        })

#------------------------------------------------------------------------------------------------------------------------


