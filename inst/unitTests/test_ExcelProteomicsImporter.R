# test_ExcelProteomicsImporter.R
#------------------------------------------------------------------------------------------------------------------------
Sys.setlocale("LC_ALL", "C")
#------------------------------------------------------------------------------------------------------------------------
library(ExcelProteomicsImporter)
library(RUnit)
library(yaml)
#------------------------------------------------------------------------------------------------------------------------
test.files <- file.path(system.file(package="ExcelProteomicsImporter", "extdata"),
                        c("test_data2_LCL57cell_line_experiment2.txt",
                          "test_data3_PBMC_experiment.txt",
                          "test_data4_Glioma_cell_line_experiment.txt",
                          "test_data_LCL57cell_line_experiment_clean.txt"))
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_ctor()
   test_validateColumnNames()
   test_groupCount()

   test_groupBoundaries_file1()
   test_groupBoundaries_file2()
   test_groupBoundaries_file3()
   test_groupBoundaries_file4()

   test_getConditionFilenames()
   test_toYAMLreadyList_file1()
   test_toYAMLreadyList_file2()
   test_toYAMLreadyList_file3()
   test_toYAMLreadyList_file4()

   test_toYAML_file1()
   test_toYAML_file2()
   test_toYAML_file3()
   test_toYAML_file4()

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
      checkTrue(ncol(tbl.raw) %in% c(15, 37, 29, 25))
      }

} # test_ctor
#------------------------------------------------------------------------------------------------------------------------
test_validateColumnNames <- function()
{
   message(sprintf("--- test_validateColumnNames"))

   for(file in test.files){
      #printf("--- %s", file)
      importer <- ExcelProteomicsImporter(file)
      checkTrue(validateColumnNames(importer))
      }

} # test_validateColumnNames
#------------------------------------------------------------------------------------------------------------------------
test_groupCount <- function()
{
   message(sprintf("--- test_groupCount"))

     # first file has just one group, so just one "avg" colname, one "stdev"
   x <- ExcelProteomicsImporter(test.files[1])
   checkEquals(getGroupCount(x), 1)
   checkEquals(length(grep("avg", colnames(getRawTable(x)))), getGroupCount(x))
   checkEquals(length(grep("stdev", colnames(getRawTable(x)))), getGroupCount(x))

     # second file has just one group, so just one "avg" colname, one "stdev"
   x <- ExcelProteomicsImporter(test.files[2])
   checkEquals(getGroupCount(x), 1)
   checkEquals(length(grep("avg", colnames(getRawTable(x)))), getGroupCount(x))
   checkEquals(length(grep("stdev", colnames(getRawTable(x)))), getGroupCount(x))

     # third file has just one group, so just one "avg" colname, one "stdev"
   x <- ExcelProteomicsImporter(test.files[3])
   checkEquals(getGroupCount(x), 1)
   checkEquals(length(grep("avg", colnames(getRawTable(x)))), getGroupCount(x))
   checkEquals(length(grep("stdev", colnames(getRawTable(x)))), getGroupCount(x))

     # fourth file has two groups, so two "avg" colnames, two "stdev"
   x <- ExcelProteomicsImporter(test.files[4])
   checkEquals(getGroupCount(x), 1)
   checkEquals(length(grep("avg", colnames(getRawTable(x)))), getGroupCount(x))
   checkEquals(length(grep("stdev", colnames(getRawTable(x)))), getGroupCount(x))

} # test_groupCount
#------------------------------------------------------------------------------------------------------------------------
test_groupBoundaries_file1 <- function()
{
   message(sprintf("--- test_groupBoundaries_file1"))
   expected.cols <- c("LCL57_IR_00Gy_1", "LCL57_IR_01Gy_1", "LCL57_IR_02Gy_1", "LCL57_IR_05Gy_1", "LCL57_IR_10Gy_1")

     #-------------------------------------------------------------------------
     # first file has just one group, so just one "avg" colname, one "stdev"
     #-------------------------------------------------------------------------

   x <- ExcelProteomicsImporter(test.files[1])
   tbl <- getRawTable(x)

   checkEquals(dim(tbl), c(39, 15))
   checkEquals(getGroupCount(x), 1)
   checkEquals(length(grep("avg", colnames(getRawTable(x)))), getGroupCount(x))
   checkEquals(length(grep("stdev", colnames(getRawTable(x)))), getGroupCount(x))

     #--------------------------------------------------------------------------------
     # one good test: avg colnames match stdev colnames, with latter having .1 suffix
     #--------------------------------------------------------------------------------

   cols <- getGroupBoundaries(x)
   avg.cols <- colnames(tbl)[cols$avg.1[1]: cols$avg.1[2]]
   checkEquals(length(avg.cols), 5)
   checkEquals(avg.cols, expected.cols)
   stdev.cols <- colnames(tbl)[cols$sd.1[1]: cols$sd.1[2]]
   checkEquals(length(stdev.cols), 5)
   checkEquals(stdev.cols, expected.cols)

} # test_groupBoundaries_file1
#------------------------------------------------------------------------------------------------------------------------
test_groupBoundaries_file2 <- function()
{
   message(sprintf("--- test_groupBoundaries_file2"))
   expected.cols <- c("PBMC+SEB_DMSO_0Gy_1",   "PBMC+SEB_DMSO_5Gy_1",   "PBMC+SEB_ATMi_0Gy_1",
                      "PBMC+SEB_ATMi_5Gy_1",   "PBMC+SEB_DNAPKi_0Gy_1", "PBMC+SEB_DNAPKi_5Gy_1",
                      "PBMC+SEB_ATRi_0Gy_1",   "PBMC+SEB_ATRi_5Gy_1",   "PBMC-SEB_DMSO_0Gy_1",
                      "PBMC-SEB_DMSO_5Gy_1",   "PBMC-SEB_ATMi_0Gy_1",   "PBMC-SEB_ATMi_5Gy_1",
                      "PBMC-SEB_DNAPKi_0Gy_1", "PBMC-SEB_DNAPKi_5Gy_1", "PBMC-SEB_ATRi_0Gy_1",
                      "PBMC-SEB_ATRi_5Gy_1")

     #-------------------------------------------------------------------------
     # first file has just one group, so just one "avg" colname, one "stdev"
     #-------------------------------------------------------------------------

   x <- ExcelProteomicsImporter(test.files[2])
   tbl <- getRawTable(x)
   checkEquals(dim(tbl), c(37, 37))
   checkEquals(getGroupCount(x), 1)
   checkEquals(length(grep("avg", colnames(getRawTable(x)))), getGroupCount(x))
   checkEquals(length(grep("stdev", colnames(getRawTable(x)))), getGroupCount(x))

     #--------------------------------------------------------------------------------
     # one good test: avg colnames match stdev colnames, with latter having .1 suffix
     #--------------------------------------------------------------------------------

   cols <- getGroupBoundaries(x)
   avg.cols <- colnames(tbl)[cols$avg.1[1]: cols$avg.1[2]]
   checkEquals(length(avg.cols), 16)
   checkEquals(avg.cols, expected.cols)
   stdev.cols <- colnames(tbl)[cols$sd.1[1]: cols$sd.1[2]]
   checkEquals(length(stdev.cols), 16)
   checkEquals(stdev.cols, expected.cols)

} # test_groupBoundaries_file2
#------------------------------------------------------------------------------------------------------------------------
test_groupBoundaries_file3 <- function()
{
   message(sprintf("--- test_groupBoundaries_file3"))
   expected.cols <- c("G7B_none_0Gy_1",  "G7B_none_5Gy_1",  "G7S_none_0Gy_1",
                      "G7S_none_5Gy_1",  "G7S_ATMi_5Gy_1",  "G7S_Olaparib_5Gy_1",
                      "R10B_none_0Gy_1", "R10B_none_5Gy_1", "R10S_none_0Gy_1",
                      "R10S_none_5Gy_1", "R10S_ATMi_5Gy_1", "R10S_Olaparib_5Gy_1")

     #-------------------------------------------------------------------------
     # first file has just one group, so just one "avg" colname, one "stdev"
     #-------------------------------------------------------------------------

   x <- ExcelProteomicsImporter(test.files[3])
   tbl <- getRawTable(x)
   checkEquals(dim(tbl), c(33, 29))
   checkEquals(getGroupCount(x), 1)
   checkEquals(length(grep("avg", colnames(tbl))), getGroupCount(x))
   checkEquals(length(grep("stdev", colnames(tbl))), getGroupCount(x))

     #--------------------------------------------------------------------------------
     # one good test: avg colnames match stdev colnames, with latter having .1 suffix
     #--------------------------------------------------------------------------------

   cols <- getGroupBoundaries(x)
   avg.cols <- colnames(tbl)[cols$avg.1[1]: cols$avg.1[2]]
   checkEquals(length(avg.cols), 12)
   checkEquals(avg.cols, expected.cols)
   stdev.cols <- colnames(tbl)[cols$sd.1[1]: cols$sd.1[2]]
   checkEquals(length(stdev.cols), 12)
   checkEquals(stdev.cols, expected.cols)

} # test_groupBoundaries_file3
#------------------------------------------------------------------------------------------------------------------------
test_groupBoundaries_file4 <- function()
{
   message(sprintf("--- test_groupBoundaries_file4"))

   expected.cols <- c("LCL57_ATMi_0Gy_1", "LCL57_ATMi_5Gy_0.25", "LCL57_ATMi_5Gy_1",
                      "LCL57_ATMi_5Gy_6", "LCL57_ATMi_5Gy_24",
                      "LCL57_DMSO_0Gy_1", "LCL57_DMSO_5Gy_0.25", "LCL57_DMSO_5Gy_1",
                      "LCL57_DMSO_5Gy_6", "LCL57_DMSO_5Gy_24")

     #-------------------------------------------------------------------------
     # first file has just one group, so just one "avg" colname, one "stdev"
     #-------------------------------------------------------------------------

   x <- ExcelProteomicsImporter(test.files[4])
   tbl <- getRawTable(x)
   checkEquals(dim(tbl), c(37, 25))
   checkEquals(getGroupCount(x), 1)
   checkEquals(length(grep("avg", colnames(tbl))), getGroupCount(x))
   checkEquals(length(grep("stdev", colnames(tbl))), getGroupCount(x))

     #--------------------------------------------------------------------------------
     # one good test: avg colnames match stdev colnames, with latter having .1 suffix
     #--------------------------------------------------------------------------------

   cols <- getGroupBoundaries(x)
   avg.cols <- colnames(tbl)[cols$avg.1[1]: cols$avg.1[2]]
   checkEquals(length(avg.cols), 10)
   checkEquals(avg.cols, expected.cols)

   stdev.cols <- colnames(tbl)[cols$sd.1[1]: cols$sd.1[2]]
   checkEquals(length(stdev.cols), 10)
   checkEquals(stdev.cols, expected.cols)

} # test_groupBoundaries_file4
#----------------------------------------------------------------------------------------------------
test_getConditionFilenames <- function()
{
   message(sprintf("--- test_getConditionFilenames"))
   x <- ExcelProteomicsImporter(test.files[1])
   colnames <- getConditionColnames(x)

   checkEquals(getConditionColnames(ExcelProteomicsImporter(test.files[1])),
               c("LCL57_IR_00Gy_1","LCL57_IR_01Gy_1","LCL57_IR_02Gy_1","LCL57_IR_05Gy_1","LCL57_IR_10Gy_1"))
   checkEquals(getConditionColnames(ExcelProteomicsImporter(test.files[2])),
               c("PBMC+SEB_DMSO_0Gy_1","PBMC+SEB_DMSO_5Gy_1","PBMC+SEB_ATMi_0Gy_1","PBMC+SEB_ATMi_5Gy_1",
                 "PBMC+SEB_DNAPKi_0Gy_1","PBMC+SEB_DNAPKi_5Gy_1","PBMC+SEB_ATRi_0Gy_1","PBMC+SEB_ATRi_5Gy_1",
                 "PBMC-SEB_DMSO_0Gy_1","PBMC-SEB_DMSO_5Gy_1","PBMC-SEB_ATMi_0Gy_1","PBMC-SEB_ATMi_5Gy_1",
                 "PBMC-SEB_DNAPKi_0Gy_1","PBMC-SEB_DNAPKi_5Gy_1","PBMC-SEB_ATRi_0Gy_1","PBMC-SEB_ATRi_5Gy_1"))

   checkEquals(getConditionColnames(ExcelProteomicsImporter(test.files[3])),
               c("G7B_none_0Gy_1","G7B_none_5Gy_1","G7S_none_0Gy_1","G7S_none_5Gy_1","G7S_ATMi_5Gy_1",
                 "G7S_Olaparib_5Gy_1","R10B_none_0Gy_1","R10B_none_5Gy_1","R10S_none_0Gy_1",
                 "R10S_none_5Gy_1","R10S_ATMi_5Gy_1","R10S_Olaparib_5Gy_1"))

   checkEquals(getConditionColnames(ExcelProteomicsImporter(test.files[4])),
               c("LCL57_ATMi_0Gy_1","LCL57_ATMi_5Gy_0.25","LCL57_ATMi_5Gy_1","LCL57_ATMi_5Gy_6",
                 "LCL57_ATMi_5Gy_24","LCL57_DMSO_0Gy_1","LCL57_DMSO_5Gy_0.25","LCL57_DMSO_5Gy_1",
                 "LCL57_DMSO_5Gy_6","LCL57_DMSO_5Gy_24"))

} # test_getConditionFilenames
#----------------------------------------------------------------------------------------------------
test_toYAMLreadyList_file1 <- function()
{
   message(sprintf("--- test_toYAMLreadyList_file1"))

   x <- ExcelProteomicsImporter(test.files[1])
   colnames <- getConditionColnames(x)
   yaml.string <- as.yaml(conditionToYAMLreadyList(x, colnames[1]))
   yaml.lines <- strsplit(yaml.string, "\n")[[1]]

   checkEquals(yaml.lines[1], "name: LCL57_IR_00Gy_1")
   checkEquals(yaml.lines[2], "sampleType: LCL57")
   checkEquals(yaml.lines[3], "treatments:")
   checkEquals(yaml.lines[4], "- name: IR")
   checkEquals(yaml.lines[5], "  level: ''")
   checkEquals(yaml.lines[6], "  units: ''")
   checkEquals(yaml.lines[7], "- name: radiation")
   checkEquals(yaml.lines[8], "  level: 0.0")
   checkEquals(yaml.lines[9], "  units: Gy")
   checkEquals(yaml.lines[10], "- name: time")
   checkEquals(yaml.lines[11], "  level: 1.0")
   checkEquals(yaml.lines[12], "  units: hours")

} # test_toYAMLreadyList_file1
#----------------------------------------------------------------------------------------------------
test_toYAMLreadyList_file2 <- function()
{
   message(sprintf("--- test_toYAMLreadyList_file2"))

   x <- ExcelProteomicsImporter(test.files[2])
   colnames <- getConditionColnames(x)
   yaml.string <- as.yaml(conditionToYAMLreadyList(x, colnames[1]))
   yaml.lines <- strsplit(yaml.string, "\n")[[1]]

   checkEquals(yaml.lines[1], "name: PBMC+SEB_DMSO_0Gy_1")
   checkEquals(yaml.lines[2], "sampleType: PBMC")
   checkEquals(yaml.lines[3], "treatments:")
   checkEquals(yaml.lines[4], "- name: SEB")
   checkEquals(yaml.lines[5], "  level: ''")
   checkEquals(yaml.lines[6], "  units: ''")
   checkEquals(yaml.lines[7], "- name: DMSO")
   checkEquals(yaml.lines[8], "  level: ''")
   checkEquals(yaml.lines[9], "  units: ''")
   checkEquals(yaml.lines[10], "- name: radiation")
   checkEquals(yaml.lines[11], "  level: 0.0")
   checkEquals(yaml.lines[12], "  units: Gy")
   checkEquals(yaml.lines[13], "- name: time")
   checkEquals(yaml.lines[14], "  level: 1.0")
   checkEquals(yaml.lines[15], "  units: hours")

} # test_toYAMLreadyList_file2
#----------------------------------------------------------------------------------------------------
test_toYAMLreadyList_file3 <- function()
{
   message(sprintf("--- test_toYAMLreadyList_file3"))

   x <- ExcelProteomicsImporter(test.files[3])
   colnames <- getConditionColnames(x)
   yaml.string <- as.yaml(conditionToYAMLreadyList(x, colnames[1]))
   yaml.lines <- strsplit(yaml.string, "\n")[[1]]

   checkEquals(yaml.lines[1],  "name: G7B_none_0Gy_1")
   checkEquals(yaml.lines[2],  "sampleType: G7B")
   checkEquals(yaml.lines[3],  "treatments:")
   checkEquals(yaml.lines[4],  "- name: none")
   checkEquals(yaml.lines[5],  "  level: ''")
   checkEquals(yaml.lines[6],  "  units: ''")
   checkEquals(yaml.lines[7],  "- name: radiation")
   checkEquals(yaml.lines[8],  "  level: 0.0")
   checkEquals(yaml.lines[9],  "  units: Gy")
   checkEquals(yaml.lines[10], "- name: time")
   checkEquals(yaml.lines[11], "  level: 1.0")
   checkEquals(yaml.lines[12], "  units: hours")

} # test_toYAMLreadyList_file3
#----------------------------------------------------------------------------------------------------
test_toYAMLreadyList_file4 <- function()
{
   message(sprintf("--- test_toYAMLreadyList_file4"))

   x <- ExcelProteomicsImporter(test.files[4])
   colnames <- getConditionColnames(x)
   yaml.string <- as.yaml(conditionToYAMLreadyList(x, colnames[1]))
   yaml.lines <- strsplit(yaml.string, "\n")[[1]]

   checkEquals(yaml.lines[1], "name: LCL57_ATMi_0Gy_1")
   checkEquals(yaml.lines[2], "sampleType: LCL57")
   checkEquals(yaml.lines[3], "treatments:")
   checkEquals(yaml.lines[4], "- name: ATMi")
   checkEquals(yaml.lines[5], "  level: ''")
   checkEquals(yaml.lines[6], "  units: ''")
   checkEquals(yaml.lines[7], "- name: radiation")
   checkEquals(yaml.lines[8], "  level: 0.0")
   checkEquals(yaml.lines[9], "  units: Gy")
   checkEquals(yaml.lines[10], "- name: time")
   checkEquals(yaml.lines[11], "  level: 1.0")
   checkEquals(yaml.lines[12], "  units: hours")

} # test_toYAMLreadyList_file4
#----------------------------------------------------------------------------------------------------
test_toYAML_file1 <- function()
{
   message(sprintf("--- test_toYAML_file1"))

   x <- ExcelProteomicsImporter(test.files[1])
   list <- toYAML(x)
   checkEquals(length(list$conditions), 5)
   checkEquals(unlist(lapply(list$conditions, function(cond) cond$name)),
      c("LfCL57_IR_00Gy_1", "LCL57_IR_01Gy_1", "LCL57_IR_02Gy_1", "LCL57_IR_05Gy_1", "LCL57_IR_10Gy_1"))

   yaml.string <- as.yaml(list)
   yaml.lines <- strsplit(yaml.string, "\n")[[1]]
   yaml.out <- sprintf("%s.yaml", sub(".txt", "", basename(test.files[1])))
   writeLines(yaml.lines, yaml.out)

} # test_toYAML_file1
#----------------------------------------------------------------------------------------------------
test_toYAML_file2 <- function()
{
   message(sprintf("--- test_toYAML_file2"))

   x <- ExcelProteomicsImporter(test.files[2])
   list <- toYAML(x)
   checkEquals(length(list$conditions), 16)
   yaml.out <- sprintf("%s.yaml", sub(".txt", "", basename(test.files[2])))
   writeLines(yaml.lines, yaml.out)

} # test_toYAML_file2
#----------------------------------------------------------------------------------------------------
test_toYAML_file3 <- function()
{
   message(sprintf("--- test_toYAML_file3"))

   x <- ExcelProteomicsImporter(test.files[3])
   list <- toYAML(x)
   checkEquals(length(list$conditions), 12)
   yaml.out <- sprintf("%s.yaml", sub(".txt", "", basename(test.files[3])))
   writeLines(yaml.lines, yaml.out)

} # test_toYAML_file3
#----------------------------------------------------------------------------------------------------
test_toYAML_file4 <- function()
{
   message(sprintf("--- test_toYAML_file4"))

   x <- ExcelProteomicsImporter(test.files[4])
   list <- toYAML(x)
   checkEquals(length(list$conditions), 10)

   yaml.out <- sprintf("%s.yaml", sub(".txt", "", basename(test.files[4])))
   writeLines(yaml.lines, yaml.out)


} # test_toYAML_file4
#----------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
