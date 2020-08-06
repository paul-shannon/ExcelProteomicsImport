#' @import org.Hs.eg.db
#' @import yaml
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
                             filename="character",
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
setGeneric('getConditionColnames', signature='obj', function(obj) standardGeneric ('getConditionColnames'))
setGeneric('conditionToYAMLreadyList', signature='obj', function(obj, colname) standardGeneric ('conditionToYAMLreadyList'))
setGeneric('toYAML', signature='obj', function(obj) standardGeneric ('toYAML'))
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
   tbl.raw <- read.table(tabDelimitedFilename, sep="\t", header=TRUE, as.is=TRUE, check.names=FALSE)
   obj <- .ExcelProteomicsImporter(filename=tabDelimitedFilename, rawTable=tbl.raw, state=state)
   validateColumnNames(obj)
   identifyGroupBoundaries(obj)
   obj

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
        stopifnot("PeptideLabel" %in% colnames | "Peptide Label" %in% colnames)
        if("PeptideLabel" %in% colnames)
            analytes <- tbl$PeptideLabel
        if("Peptide Label" %in% colnames)
            analytes <- tbl[, "Peptide Label"]

        stopifnot(length(unique(analytes)) == nrow(tbl))

             # just one "avg" and one "stdev" column allowed
        avg.cols <- grep("avg", colnames)
        sd.cols  <- grep("stdev", colnames)
        stopifnot(length(avg.cols) == 1)
        stopifnot(length(sd.cols) == 1)
        obj@state$groupCount <- 1         # remnant of abandoned support for multiple groups
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
        colnames <- colnames(obj@rawTable)
        avg.cols <- grep("avg", colnames)
        sd.cols  <- grep("stdev", colnames)
        if(obj@state$groupCount == 1){
           obj@state$group.1.avg <- c(avg.cols[1] + 1, sd.cols[1] - 1)
           obj@state$group.1.sd  <- c(sd.cols[1] + 1, length(colnames))
           obj@state$group.2.avg <- NA
           obj@state$group.2.sd <- NA
           }
        if(obj@state$groupCount == 2){
           obj@state$group.1.avg <- c(avg.cols[1] + 1, avg.cols[2] - 1)
           obj@state$group.1.sd  <- c(sd.cols[1] + 1,  sd.cols[2] - 1)
           obj@state$group.2.avg <- c(avg.cols[2] + 1, sd.cols[1] - 1)
           obj@state$group.2.sd <-  c(sd.cols[2] + 1,  length(colnames))
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
        obj@state$groupCount
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
         list(avg.1=obj@state$group.1.avg,
              avg.2=obj@state$group.2.avg,
              sd.1 =obj@state$group.1.sd,
              sd.2 =obj@state$group.2.sd)
        })

#------------------------------------------------------------------------------------------------------------------------
#' getConditionColnames
#'
#' @rdname getConditionColnames
#' @aliases getConditionColnames
#'
#' @export
#'
setMethod('getConditionColnames', 'ExcelProteomicsImporter',

     function(obj){
        col.start <-  obj@state$group.1.avg[1]
        col.end <-  obj@state$group.1.avg[2]
        colnames(obj@rawTable)[col.start:col.end]
        })

#------------------------------------------------------------------------------------------------------------------------
#' conditionToYAMLreadyList
#'
#' @rdname conditionToYAMLreadyList
#' @aliases conditionToYAMLreadyList
#'
#' @export
#'
setMethod('conditionToYAMLreadyList', 'ExcelProteomicsImporter',

    function(obj, colname){
       string.is.numeric <- function(string){
          tryCatch(expr = is.numeric(as.numeric(string)), warning = function(w) {return(FALSE)})
          }
       sample.name = colname
       tokens <-  strsplit(colname, split="\\+|_")[[1]]
       sample.type <- tokens[1]
       elements <- tokens[-1]
       elementCount <- 0
       treatments <- list()
       last.element <- elements[length(elements)]
       for(el in elements[-(length(elements))]){
           elementCount <- elementCount + 1
           if(grepl("Gy", el)){
               name <- "radiation"
               level <- as.numeric(sub("Gy", "", el))
               units <- "Gy"
           } else {
               name <- el
               level <- ''
               units <- ''
           }
           treatments[[elementCount]] <- list(name=name, level=level, units=units)
           } # for el
          # now work on the last element
       elementCount <- elementCount + 1
       if(string.is.numeric(last.element)){ # assume this is time, and units are hours?
          timepoint <- as.numeric(last.element)
          treatments[[elementCount]] <- list(name="time", level=timepoint, units="hours")
       } else {
          treatments[[elementCount]] <- list(name=last.element, level="", units="")
          }
       list(name=sample.name, sampleType=sample.type, treatments=treatments)
       }) # conditionToYAMLreadyList

#------------------------------------------------------------------------------------------------------------------------
#' get the matrix of averages
#'
#' @rdname toYAML
#' @aliases toYAML
#'
#' @export
#'
setMethod('toYAML', 'ExcelProteomicsImporter',

     function(obj){
         colnames <- getConditionColnames(obj)
         conditionLists <- lapply(colnames, function(colname) conditionToYAMLreadyList(obj, colname))
         sample.name <- strsplit(colnames[1], "_")[[1]][1]  # first token relaibly the sample name
         list(filename=obj@filename,
              sample=sample.name,
              mutation="",
              conditions=conditionLists)
        })

#------------------------------------------------------------------------------------------------------------------------


