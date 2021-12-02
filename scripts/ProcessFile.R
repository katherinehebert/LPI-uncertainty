ProcessFile <- function (DatasetName, ref_year, MODEL_SELECTION_FLAG, GAM_GLOBAL_FLAG, 
          DATA_LENGTH_MIN, AVG_TIME_BETWEEN_PTS_MAX, GLOBAL_GAM_FLAG_SHORT_DATA_FLAG, 
          AUTO_DIAGNOSTIC_FLAG, LAMBDA_MIN, LAMBDA_MAX, ZERO_REPLACE_FLAG, 
          OFFSET_ALL, OFFSET_NONE, OFFSET_DIFF, LINEAR_MODEL_SHORT_FLAG, 
          CAP_LAMBDAS, SHOW_PROGRESS, basedir) 
{
  md5val <- tools::md5sum(DatasetName)
  Data = read.table(DatasetName, header = TRUE)
  SpeciesSSet = Data[1]
  IDSSet = Data[2]
  YearSSet = Data[3]
  PopvalueSSet = Data[4]
  rm(Data)
  FinalYear = max(YearSSet)
  if (min(YearSSet) < ref_year) {
    InitialYear = ref_year
  }
  else {
    InitialYear = min(YearSSet)
  }
  InitialYear = ref_year
  cat("Calculating LPI for Species\n")
  pop_lambda_filename <- file.path(basedir, gsub(".txt", "_PopLambda.txt", 
                                                 DatasetName))
  Pop_Headers <- t(c("population_id", as.vector(InitialYear:FinalYear)))
  write.table(Pop_Headers, file = pop_lambda_filename, sep = ",", 
              eol = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
  SpeciesLambda = CalcLPI_saveGAM(Species = SpeciesSSet, ID = IDSSet, 
                          Year = YearSSet, Popvalue = PopvalueSSet, InitialYear = InitialYear, 
                          FinalYear = FinalYear, DatasetName = DatasetName, MODEL_SELECTION_FLAG = MODEL_SELECTION_FLAG, 
                          GAM_GLOBAL_FLAG = GAM_GLOBAL_FLAG, DATA_LENGTH_MIN = DATA_LENGTH_MIN, 
                          AVG_TIME_BETWEEN_PTS_MAX = AVG_TIME_BETWEEN_PTS_MAX, 
                          GLOBAL_GAM_FLAG_SHORT_DATA_FLAG = GLOBAL_GAM_FLAG_SHORT_DATA_FLAG, 
                          AUTO_DIAGNOSTIC_FLAG = AUTO_DIAGNOSTIC_FLAG, LAMBDA_MIN = LAMBDA_MIN, 
                          LAMBDA_MAX = LAMBDA_MAX, ZERO_REPLACE_FLAG = ZERO_REPLACE_FLAG, 
                          OFFSET_ALL = OFFSET_ALL, OFFSET_NONE = OFFSET_NONE, 
                          OFFSET_DIFF = OFFSET_DIFF, LINEAR_MODEL_SHORT_FLAG = LINEAR_MODEL_SHORT_FLAG, 
                          CAP_LAMBDAS = CAP_LAMBDAS, show_progress = SHOW_PROGRESS, 
                          basedir = basedir)
  DataFileName = file.path(basedir, "lpi_temp", paste0(md5val, 
                                                       "_splambda.csv"))
  cat(sprintf("Saving species lambda to file: %s\n", DataFileName))
  write.table(SpeciesLambda, DataFileName, sep = ",", col.names = FALSE, 
              row.names = FALSE)
  sp.count <- as.data.frame(table(SpeciesSSet))
  DataFileName = file.path(basedir, gsub(".txt", "_lambda.csv", 
                                         DatasetName))
  rownames(SpeciesLambda) <- t(unique(SpeciesSSet))
  colnames(SpeciesLambda) <- InitialYear:FinalYear
  sorted_lambdas <- SpeciesLambda[order(rownames(SpeciesLambda)), 
  ]
  sorted_lambdas_count = cbind(sp.count, sorted_lambdas)
  cat(sprintf("Saving species lambda to file: %s\n", DataFileName))
  write.table(sorted_lambdas_count, DataFileName, sep = ",", 
              col.names = NA)
  cat("Calculating DTemp\n")
  DTemp = matrix(0, 1, dim(SpeciesLambda)[2])
  for (I in 1:dim(SpeciesLambda)[2]) {
    YearData = SpeciesLambda[, I]
    if (!CAP_LAMBDAS) {
      Index = which(YearData != -1)
    }
    else {
      Index = which(!is.na(YearData))
    }
    if (length(Index) > 0) {
      DTemp[I] = mean(YearData[Index])
    }
    else DTemp[I] = -99
  }
  DataFileName = file.path(basedir, "lpi_temp", paste0(md5val, 
                                                       "_dtemp.csv"))
  cat("Saving DTemp to file: ", DataFileName, "\n")
  colnames(DTemp) <- InitialYear:FinalYear
  write.table(DTemp, DataFileName, sep = ",", row.names = FALSE)
  DataFileName = file.path(basedir, gsub(".txt", "_dtemp.csv", 
                                         DatasetName))
  cat("Saving DTemp to file: ", DataFileName, "\n")
  write.table(DTemp, DataFileName, sep = ",", row.names = FALSE)
  return(dim(SpeciesLambda)[2])
}
