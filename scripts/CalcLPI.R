CalcLPI_saveGAM <- function (Species, ID, Year, Popvalue, InitialYear, FinalYear, 
          DatasetName, MODEL_SELECTION_FLAG, GAM_GLOBAL_FLAG, DATA_LENGTH_MIN, 
          AVG_TIME_BETWEEN_PTS_MAX, GLOBAL_GAM_FLAG_SHORT_DATA_FLAG, 
          AUTO_DIAGNOSTIC_FLAG, LAMBDA_MIN, LAMBDA_MAX, ZERO_REPLACE_FLAG, 
          OFFSET_ALL, OFFSET_NONE, OFFSET_DIFF, LINEAR_MODEL_SHORT_FLAG, 
          CAP_LAMBDAS = FALSE, show_progress = FALSE, basedir = ".") 
{
  noRecs = max(dim(Popvalue))
  sNames = unique(Species)
  sID = unique(ID)
  noSpecies = max(dim(sNames))
  noPop = max(dim(unique(ID)))
  PopNotProcessed = matrix(0, 1, noPop)
  MethodFlag = matrix(0, 1, noPop)
  PopNotProcessedCounter = 0
  PopProcessedGAMCounter = 0
  sNamesCounter = 0
  sNamesArray = sNames
  sIDArray = sID
  PopProcessedGAM = matrix(0, 1, noPop)
  cat(sprintf("Number of species: %s (in %s populations)\n", 
              noSpecies, noPop))
  SpeciesLambda = matrix(0, noSpecies, FinalYear - InitialYear + 
                           1)
  if (show_progress) {
    prog <- txtProgressBar(min = 0, max = noSpecies, char = "*", 
                           style = 3)
  }
  MethodFlagLoop = 0
  for (I in 1:noSpecies) {
    Index = 1
    if (length(which(objects() == "sIndex")) != 0) 
      rm(sIndex)
    sIndex = which(Species == toString(sNames[I, 1]))
    PopID = unique(ID[sIndex, 1])
    if (length(which(objects() == "PopLambda")) != 0) 
      rm(PopLambda)
    PopIDSize = length(PopID)
    PopLambda = matrix(-1, PopIDSize, FinalYear - InitialYear + 
                         1)
    JIndex = 1
    for (J in 1:PopIDSize) {
      IndexPop = which(ID == PopID[J])
      YearPop = Year[IndexPop, 1]
      PopN = Popvalue[IndexPop, 1]
      DataTimeLength = max(YearPop) - min(YearPop)
      AvgTimeBetweenPts = DataTimeLength/length(PopN)
      if ((length(PopN) >= DATA_LENGTH_MIN) & (AvgTimeBetweenPts < 
                                               AVG_TIME_BETWEEN_PTS_MAX)) {
        if (OFFSET_ALL) {
          cat(sprintf("Offsetting all time-series by 1 to avoid log(0)\n"))
          PopN <- PopN + 1
        }
        else if (OFFSET_NONE) {
          IndexZero = which(PopN == 0)
          if (length(IndexZero) > 0) {
            OffsetVal = 1e-17
            PopN = PopN + OffsetVal
          }
          else {
            PopN <- PopN
          }
        }
        else if (OFFSET_DIFF) {
          IndexZero = which(PopN == 0)
          if (length(IndexZero) > 0) {
            if (mean(PopN) == 0) {
              OffsetVal = 1e-17
            }
            else {
              if (max(PopN) >= 1) {
                PopN = PopN + 1
              }
              else {
                IndexNonZero = which(PopN != 0)
                OffsetVal = mean(PopN[IndexNonZero]) * 
                  0.01
                PopN = PopN + OffsetVal
              }
            }
          }
        }
        else {
          IndexZero = which(PopN == 0)
          if (ZERO_REPLACE_FLAG == 1) {
            if (mean(PopN) == 0) {
              OffsetVal = 1e-17
            }
            else {
              IndexNonZero = which(PopN != 0)
              OffsetVal = mean(PopN[IndexNonZero]) * 
                0.01
            }
          }
          else {
            IndexNonZero = which(PopN != 0)
            OffsetVal = min(PopN[IndexNonZero])
          }
          if (ZERO_REPLACE_FLAG == 2) {
            if (mean(PopN) == 0) {
              OffsetVal = 1e-17
            }
            else {
              IndexNonZero = which(PopN != 0)
              OffsetVal = 1
            }
          }
          if (length(IndexZero) > 0) 
            PopN = PopN + OffsetVal
        }
        SortResults = sort(YearPop, index.return = TRUE)
        YearPop = SortResults$x
        TempI = SortResults$ix
        PopN = PopN[TempI]
        if (length(which(objects() == "YearPopInt")) != 
            0) 
          rm(YearPopInt)
        if (length(which(objects() == "PopNInt")) != 
            0) 
          rm(PopNInt)
        YearPopInt = YearPop[1]:YearPop[length(YearPop)]
        PopNLog = log(PopN)
        Flag = 0
        GAMFlag = GAM_GLOBAL_FLAG
        if (mean(PopN) == PopN[1]) {
          GAMFlag = 0
        }
        if (GAMFlag == 1) {
          if (MODEL_SELECTION_FLAG == 0) {
            SmoothParm = round(length(PopN)/2)
            if (SmoothParm >= 3) {
              s <- mgcv::s
              model <- mgcv::gam(PopNLog ~ s(YearPop, 
                                             k = SmoothParm), fx = TRUE)
              # save the GAM
              saveRDS(model, paste0(here::here("models/rlpi/"), gsub("_pops.txt", "", DatasetName), ".rds"))
              if (AUTO_DIAGNOSTIC_FLAG == 1) {
                rsd <- residuals(model)
                s <- mgcv::s
                modelres <- mgcv::gam(rsd ~ s(YearPop, 
                                              k = length(PopN), bs = "cs"), gamma = 1.4)
                if ((abs(sum(modelres$edf) - 1)) < 0.01) {
                  PopNInt <- predict(model, data.frame(YearPop = YearPopInt))
                  PopNInt = exp(PopNInt)
                  Flag = 1
                  PopProcessedGAMCounter = PopProcessedGAMCounter + 
                    1
                  PopProcessedGAM[PopProcessedGAMCounter] = PopID[J]
                }
              }
              else {
                summary(model)
                readline(prompt = "Press any key to continue")
                plot(model, pages = 1, residuals = TRUE, 
                     all.terms = TRUE, shade = TRUE, shade.col = 2)
                readline(prompt = "Press any key to continue")
                mgcv::gam.check(model)
                Char = readline(prompt = "Press 'Y' to accept model, 'N' to reject GAM model and use default method")
                while ((Char != "Y") & (Char != "N")) {
                  Char = readline(prompt = "Press 'Y' to accept model, 'N' to reject GAM model and use default method")
                }
                if (Char == "Y") {
                  PopNInt <- predict(model, data.frame(YearPop = YearPopInt))
                  PopNInt = exp(PopNInt)
                  Flag = 1
                  PopProcessedGAMCounter = PopProcessedGAMCounter + 
                    1
                  PopProcessedGAM[PopProcessedGAMCounter] = PopID[J]
                }
              }
            }
          }
          else {
            if (length(PopN) >= 6) {
              SmoothParm = 3
              if (AUTO_DIAGNOSTIC_FLAG == 1) {
                while ((length(PopN) >= SmoothParm) & 
                       (Flag == 0)) {
                  s <- mgcv::s
                  model <- mgcv::gam(PopNLog ~ s(YearPop, 
                                                 k = SmoothParm), fx = TRUE)
                  rsd <- residuals(model)
                  modelres <- mgcv::gam(rsd ~ s(YearPop, 
                                                k = length(PopN), bs = "cs"), gamma = 1.4)
                  if ((abs(sum(modelres$edf) - 1)) < 
                      0.01) {
                    Flag = 1
                    PopNInt <- predict(model, data.frame(YearPop = YearPopInt))
                    PopNInt = exp(PopNInt)
                    PopProcessedGAMCounter = PopProcessedGAMCounter + 
                      1
                    PopProcessedGAM[PopProcessedGAMCounter] = PopID[J]
                  }
                  else {
                    SmoothParm = SmoothParm + 1
                  }
                }
              }
              else {
                while ((length(PopN) >= SmoothParm) & 
                       (Flag == 0)) {
                  s <- mgcv::s
                  model <- mgcv::gam(PopNLog ~ s(YearPop, 
                                                 k = SmoothParm), fx = TRUE)
                  summary(model)
                  readline(prompt = "Press any key to continue")
                  plot(model, pages = 1, residuals = TRUE, 
                       all.terms = TRUE, shade = TRUE, 
                       shade.col = 2)
                  readline(prompt = "Press any key to continue")
                  mgcv::gam.check(model)
                  Char = readline(prompt = "Press 'Y' to accept model, 'N' to reject model")
                  while ((Char != "Y") & (Char != "N")) {
                    Char = readline(prompt = "Press 'Y' to accept model, 'N' to reject model")
                  }
                  if (Char == "Y") {
                    PopNInt <- predict(model, data.frame(YearPop = YearPopInt))
                    PopNInt = exp(PopNInt)
                    Flag = 1
                    PopProcessedGAMCounter = PopProcessedGAMCounter + 
                      1
                    PopProcessedGAM[PopProcessedGAMCounter] = PopID[J]
                  }
                  else {
                    SmoothParm = SmoothParm + 1
                  }
                }
              }
            }
          }
        }
        if (Flag == 0) {
          if (GLOBAL_GAM_FLAG_SHORT_DATA_FLAG == 1) {
            SmoothParm = length(PopN)
            s <- mgcv::s
            model <- mgcv::gam(PopNLog ~ s(YearPop, 
                                           k = SmoothParm), fx = TRUE)
            PopNInt <- predict(model, data.frame(YearPop = YearPopInt))
            PopNInt = exp(PopNInt)
            PopProcessedGAMCounter = PopProcessedGAMCounter + 
              1
            PopProcessedGAM[PopProcessedGAMCounter] = PopID[J]
          }
          else {
            if (LINEAR_MODEL_SHORT_FLAG == TRUE) {
              MethodFlagLoop = MethodFlagLoop + 1
              MethodFlag[MethodFlagLoop] = PopID[J]
              model <- lm(PopNLog ~ YearPop)
              PopNInt <- predict(model, data.frame(YearPop = YearPopInt))
              PopNInt = exp(PopNInt)
            }
            else {
              MethodFlagLoop = MethodFlagLoop + 1
              MethodFlag[MethodFlagLoop] = PopID[J]
              PopNInt = matrix(-1, 1, length(YearPopInt))
              for (K in 1:length(YearPopInt)) {
                k = which(YearPop == YearPopInt[K])
                if (length(k) > 0) {
                  PopNInt[K] = PopN[k]
                }
                else {
                  YearStart = YearPopInt[K]
                  YearStart = YearStart - 1
                  k = which(YearPop == YearStart)
                  while (length(k) == 0) {
                    YearStart = YearStart - 1
                    k = which(YearPop == YearStart)
                  }
                  PopNStart = PopN[k]
                  YearEnd = YearPopInt[K]
                  YearEnd = YearEnd + 1
                  k = which(YearPop == YearEnd)
                  while (length(k) == 0) {
                    YearEnd = YearEnd + 1
                    k = which(YearPop == YearEnd)
                  }
                  PopNEnd = PopN[k]
                  PopNInt[K] = PopNStart * ((PopNEnd/PopNStart)^((YearPopInt[K] - 
                                                                    YearStart)/(YearEnd - YearStart)))
                }
              }
            }
          }
        }
        YearPop = InitialYear:FinalYear
        PopN = matrix(0, 1, length(YearPop))
        k = which(PopNInt == 0)
        k1 = which(PopNInt > 0)
        TempVal = 0
        if (length(k) > 0) {
          if (length(k1) > 0) {
            if (ZERO_REPLACE_FLAG == 1) {
              TempVal = mean(PopNInt[k1]) * 0.01
            }
            else {
              TempVal = min(PopNInt[k1])
            }
            PopNInt = PopNInt + TempVal
          }
        }
        for (K in InitialYear:FinalYear) {
          k = which(YearPopInt == K)
          if (length(k) > 0) {
            if (PopNInt[k] == 0) {
              PopN[K - InitialYear + 1] = -1
            }
            else PopN[K - InitialYear + 1] = log10(PopNInt[k])
          }
          else PopN[K - InitialYear + 1] = -1
        }
        PopLambda[JIndex, 1] = 1
        StartYear = InitialYear + 1
        for (K in StartYear:FinalYear) {
          if ((PopN[K - InitialYear + 1] != -1) & (PopN[K - 
                                                        InitialYear] != -1)) 
            PopLambda[JIndex, K - InitialYear + 1] = PopN[K - 
                                                            InitialYear + 1] - PopN[K - InitialYear]
          else PopLambda[JIndex, K - InitialYear + 1] = -1
        }
        JIndex = JIndex + 1
      }
      else {
        PopNotProcessedCounter = PopNotProcessedCounter + 
          1
        PopNotProcessed[PopNotProcessedCounter] = PopID[J]
      }
    }
    PopData <- cbind(as.vector(PopID), PopLambda)
    pop_lambda_filename <- file.path(basedir, gsub(".txt", 
                                                   "_PopLambda.txt", DatasetName))
    write.table(PopData, sep = ",", eol = "\n", file = pop_lambda_filename, 
                quote = FALSE, col.names = FALSE, row.names = FALSE, 
                append = TRUE)
    EndYear = FinalYear - InitialYear + 1
    for (K in 1:EndYear) {
      k = which(PopLambda[, K] != -1)
      if (length(k) > 0) {
        PopLambdaTemp = PopLambda[k, K]
        IndexTemp = which(PopLambdaTemp < LAMBDA_MAX)
        IndexTempBad_max = which(PopLambdaTemp > LAMBDA_MAX)
        if (length(IndexTemp) > 0) {
          PopLambdaTemp1 = PopLambdaTemp[IndexTemp]
          IndexTemp = which(PopLambdaTemp1 > LAMBDA_MIN)
          IndexTempBad_min = which(PopLambdaTemp1 < 
                                     LAMBDA_MIN)
          if (length(IndexTemp) > 0) {
            if (CAP_LAMBDAS) {
              SpeciesLambda[I, K] = mean(c(PopLambdaTemp1[IndexTemp], 
                                           rep(LAMBDA_MAX, length(IndexTempBad_max)), 
                                           rep(LAMBDA_MIN, length(IndexTempBad_min))))
            }
            else {
              SpeciesLambda[I, K] = mean(PopLambdaTemp1[IndexTemp])
            }
          }
          else {
            if (CAP_LAMBDAS) {
              SpeciesLambda[I, K] = LAMBDA_MIN
            }
            else {
              SpeciesLambda[I, K] = NA
            }
          }
        }
        else {
          if (CAP_LAMBDAS) {
            SpeciesLambda[I, K] = LAMBDA_MAX
          }
          else {
            SpeciesLambda[I, K] = NA
          }
        }
      }
      else {
        SpeciesLambda[I, K] = NA
      }
    }
    sNamesCounter = sNamesCounter + 1
    sNamesArray[sNamesCounter, 1] = sNames[I, 1]
    sIDArray[sNamesCounter] = ID[I, 1]
    if (show_progress) 
      setTxtProgressBar(prog, I)
  }
  if (show_progress) 
    close(prog)
  cat("\n")
  PopNotProcessed1 = PopNotProcessed[1, 1:PopNotProcessedCounter]
  write.table(PopNotProcessed1, file = file.path(basedir, 
                                                 "lpi_temp", "PopNotProcessed.txt"))
  MethodFlag1 = MethodFlag[1, 1:MethodFlagLoop]
  if (LINEAR_MODEL_SHORT_FLAG == 1) {
    write.table(MethodFlag1, file = file.path(basedir, "lpi_temp", 
                                              "PopProcessed_LM.txt"))
  }
  else {
    write.table(MethodFlag1, file = file.path(basedir, "lpi_temp", 
                                              "PopProcessed_Chain.txt"))
  }
  PopProcessedGAM1 = PopProcessedGAM[1, 1:PopProcessedGAMCounter]
  write.table(PopProcessedGAM1, file = file.path(basedir, 
                                                 "lpi_temp", "PopProcessedGAM.txt"))
  sNamesArray1 = sNamesArray[1:sNamesCounter, 1]
  write.table(sNamesArray1, file = file.path(basedir, "lpi_temp", 
                                             "SpeciesName.txt"), quote = FALSE)
  Headers <- t(c("Species", as.vector(InitialYear:FinalYear)))
  SpeciesData <- cbind(as.vector(sNamesArray), SpeciesLambda)
  lambda_filename <- file.path(basedir, gsub(".txt", "_Lambda.txt", 
                                             DatasetName))
  write.table(Headers, file = lambda_filename, sep = ",", 
              eol = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
  write.table(SpeciesData, sep = ",", eol = "\n", file = lambda_filename, 
              quote = FALSE, col.names = FALSE, row.names = FALSE, 
              append = TRUE)
  return(SpeciesLambda)
}
