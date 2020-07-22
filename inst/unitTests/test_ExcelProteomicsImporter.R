# test_ExcelProteomicsImporter.R
#------------------------------------------------------------------------------------------------------------------------
library(ExcelProteomicsImporter)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_ctor()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_ctor <- function()
{
   message(sprintf("--- test_ctor"))

   filename <- "test_data_LCL57cell_line_experiment.txt"
   file.path <- system.file(package="ExcelProteomicsImporter", "extdata", filename)
   checkTrue(file.exists(file.path))

   importer <- ExcelProteomicsImporter(file.path)
   tbl.raw <- getRawTable(importer)
   checkEquals(dim(tbl.raw), c(37, 27))

} # test_ctor
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
