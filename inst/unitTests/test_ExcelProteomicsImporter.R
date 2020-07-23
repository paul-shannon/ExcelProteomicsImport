# test_ExcelProteomicsImporter.R
#------------------------------------------------------------------------------------------------------------------------
library(ExcelProteomicsImporter)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
test.files <- file.path(system.file(package="ExcelProteomicsImporter", "extdata"),
                        c("test_data2_LCL57cell_line_experiment2.txt",
                          "test_data3_PBMC_experiment.txt",
                          "test_data4_Glioma_cell_line_experiment.txt",
                          "test_data_LCL57cell_line_experiment.txt"))
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_ctor()
   test_validateColumnNames()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_ctor <- function()
{
   message(sprintf("--- test_ctor"))

   for(file in test.files){
      checkTrue(file.exists(file))
      importer <- ExcelProteomicsImporter(file)
      tbl.raw <- getRawTable(importer)
      checkTrue(nrow(tbl.raw) %in% c(39, 37, 33, 37))
      checkTrue(ncol(tbl.raw) %in% c(15, 37, 29, 27))
      }

} # test_ctor
#------------------------------------------------------------------------------------------------------------------------
test_validateColumnNames <- function()
{
   message(sprintf("--- test_validateColumnNames"))

   for(file in test.files){
      printf("--- %s", file)
      importer <- ExcelProteomicsImporter(file)
      checkTrue(validateColumnNames(importer))
      }

} # test_validateColumnNames
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
