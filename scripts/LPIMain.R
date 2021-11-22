# customised to save the bootstrap samples
LPIMain_custom <- function (infile = "Infile.txt", basedir = ".", REF_YEAR = 1970, 
          PLOT_MAX = 2017, force_recalculation = 0, use_weightings = 0, 
          use_weightings_B = 0, title = "", CI_FLAG = 1, LEV_FLAG = 0, 
          SWITCH_PT_FLAG = 0, BOOT_STRAP_SIZE = 100, save_plots = 1, 
          plot_lpi = 1, goParallel = FALSE, MODEL_SELECTION_FLAG = 0, 
          GAM_GLOBAL_FLAG = 1, DATA_LENGTH_MIN = 2, AVG_TIME_BETWEEN_PTS_MAX = 100, 
          GLOBAL_GAM_FLAG_SHORT_DATA_FLAG = 0, AUTO_DIAGNOSTIC_FLAG = 1, 
          LAMBDA_MIN = -1, LAMBDA_MAX = 1, ZERO_REPLACE_FLAG = 1, 
          OFFSET_ALL = 0, OFFSET_NONE = FALSE, OFFSET_DIFF = FALSE, 
          LINEAR_MODEL_SHORT_FLAG = FALSE, CAP_LAMBDAS = TRUE, VERBOSE = TRUE, 
          SHOW_PROGRESS = TRUE, scenario_name) 
{
  ptm <- proc.time()
  `%op%` <- if (goParallel) 
    foreach::`%dopar%`
  else foreach::`%do%`
  doParallel::registerDoParallel()
  success = dir.create(basedir, showWarnings = FALSE)
  if (success) {
    print(sprintf("** Created folder: %s", basedir))
  }
  dir.create(file.path(basedir, "lpi_temp"), showWarnings = FALSE)
  if (success) {
    print(sprintf("** Created folder: %s", file.path(basedir, 
                                                     "lpi_temp")))
  }
  FileTable = read.table(infile, header = TRUE)
  FileNames = FileTable$FileName
  Group = FileTable[2]
  GroupList = unique(Group[[1]])
  Weightings = FileTable[3]
  if (use_weightings == 1) {
    WeightingsA = FileTable[3]
    cat(sprintf("Weightings...\n"))
    for (i in 1:length(GroupList)) {
      print(paste("Group:", GroupList[i]))
      cat("\t")
      print(Weightings[Group == GroupList[i]])
      cat("\t")
      print("Normalised weights (sum to 1)")
      Weightings[Group == GroupList[i]] = Weightings[Group == 
                                                       GroupList[i]]/sum(Weightings[Group == GroupList[i]])
      cat("\t")
      print(Weightings[Group == GroupList[i]])
    }
    cat("\n")
  }
  if (use_weightings_B == 1) {
    FileWeightingsB = FileTable[4]
    WeightingsB = unique(cbind(Group, FileWeightingsB))$WeightingB
    cat(sprintf("WeightingsB...\n"))
    print(WeightingsB)
    cat("\n")
  }
  NoGroups = length(unique(Group[[1]]))
  cat("Number of groups: ", NoGroups, "\n")
  NoFiles = max(dim(Group))
  DSizes <- foreach::foreach(FileNo = 1:NoFiles, .combine = cbind) %op% 
    {
      md5val <- tools::md5sum(toString(FileNames[FileNo]))
      if ((force_recalculation == 1) || (!file.exists(file.path(basedir, 
                                                                "lpi_temp", paste0(md5val, "_dtemp.csv")))) || 
          (!file.exists(file.path(basedir, "lpi_temp", 
                                  paste0(md5val, "_splambda.csv"))))) {
        cat(sprintf("processing file: %s\n", toString(FileNames[FileNo])))
        ProcessFile(DatasetName = toString(FileNames[FileNo]), 
                    ref_year = REF_YEAR, MODEL_SELECTION_FLAG = MODEL_SELECTION_FLAG, 
                    GAM_GLOBAL_FLAG = GAM_GLOBAL_FLAG, DATA_LENGTH_MIN = DATA_LENGTH_MIN, 
                    AVG_TIME_BETWEEN_PTS_MAX = AVG_TIME_BETWEEN_PTS_MAX, 
                    GLOBAL_GAM_FLAG_SHORT_DATA_FLAG = GLOBAL_GAM_FLAG_SHORT_DATA_FLAG, 
                    AUTO_DIAGNOSTIC_FLAG = AUTO_DIAGNOSTIC_FLAG, 
                    LAMBDA_MIN = LAMBDA_MIN, LAMBDA_MAX = LAMBDA_MAX, 
                    ZERO_REPLACE_FLAG = ZERO_REPLACE_FLAG, OFFSET_ALL = OFFSET_ALL, 
                    OFFSET_NONE = OFFSET_NONE, OFFSET_DIFF = OFFSET_DIFF, 
                    LINEAR_MODEL_SHORT_FLAG = LINEAR_MODEL_SHORT_FLAG, 
                    CAP_LAMBDAS = CAP_LAMBDAS, SHOW_PROGRESS = SHOW_PROGRESS, 
                    basedir = basedir)
      }
    }
  DSize = PLOT_MAX - REF_YEAR + 2
  SpeciesLambdaArray = data.frame(NULL)
  SpeciesNamesArray = data.frame(NULL)
  DTempArrayTemp <- matrix(data = NA, nrow = NoFiles, ncol = DSize)
  DTempArray = data.frame(DTempArrayTemp)
  DataSizeArray = matrix(0, NoFiles, 2)
  fileindex = NULL
  for (FileNo in 1:NoFiles) {
    md5val <- tools::md5sum(toString(as.character(FileNames[FileNo])))
    FileName = file.path(basedir, "lpi_temp", paste0(md5val, 
                                                     "_splambda.csv"))
    SpeciesLambda = read.table(FileName, header = FALSE, 
                               sep = ",")
    debug_print(VERBOSE, sprintf("Loading previously analysed species lambda file for '%s' from MD5 hash: %s\n", 
                                 as.character(FileNames[FileNo]), FileName))
    species_names = read.table(file.path(basedir, "lpi_temp/SpeciesName.txt"))
    cat(sprintf("%s, Number of species: %s\n", as.character(FileNames[FileNo]), 
                dim(SpeciesLambda)[1]))
    SpeciesLambdaArray <- plyr::rbind.fill(SpeciesLambdaArray, 
                                           SpeciesLambda)
    SpeciesNamesArray <- plyr::rbind.fill(SpeciesNamesArray, 
                                          species_names)
    fileindex = c(fileindex, rep(FileNo, dim(SpeciesLambda)[1]))
    FileName = file.path(basedir, "lpi_temp", paste0(md5val, 
                                                     "_dtemp.csv"))
    debug_print(VERBOSE, sprintf("Loading previously analysed dtemp file from MD5 hash: %s\n", 
                                 FileName))
    DTemp = read.table(FileName, header = T, sep = ",")
    DTempArray[FileNo, 1:dim(DTemp)[2]] = t(DTemp)
  }
  f_name = file = file.path(basedir, gsub(".txt", "_dtemp_array.txt", 
                                          infile))
  cat("Saving DTemp Array to file: ", f_name, "\n")
  write.table(DTempArray, f_name)
  dtemp_df <- data.frame(filenames = FileNames, dtemps = DTempArray)
  colnames(dtemp_df) <- c("filename", seq(REF_YEAR, REF_YEAR + 
                                            DSize - 1))
  if (save_plots) {
    variable <- value <- filename <- NULL
    width = 20
    height = 8
    pdf(file.path(basedir, gsub(".txt", "_dtemp_array_plot.pdf", 
                                infile)), width = width, height = height)
    df.m <- reshape2::melt(dtemp_df, id.vars = "filename")
    df.m$value[df.m$value == -99] = NA
    p_line <- ggplot2::ggplot(df.m, ggplot2::aes(variable, 
                                                 value, group = filename, col = filename)) + ggplot2::geom_line() + 
      ggplot2::theme(text = ggplot2::element_text(size = 16), 
                     axis.text.x = ggplot2::element_text(size = 8, 
                                                         angle = 90, hjust = 1), legend.position = "bottom") + 
      ggplot2::guides(col = ggplot2::guide_legend(nrow = 6))
    print(p_line)
    dev.off()
  }
  f_name = file = file.path(basedir, gsub(".txt", "_dtemp_array_named.csv", 
                                          infile))
  cat("Saving DTemp Array with filesnames to file: ", f_name, 
      "\n")
  write.csv(dtemp_df, f_name, row.names = FALSE)
  t1 <- proc.time() - ptm
  cat(sprintf("[Calculating LPI...] System: %f, User: %f, Elapsed: %f\n", 
              t1[1], t1[2], t1[3]))
  I = calculate_index(DTempArray, fileindex, DSize, Group, 
                      Weightings, use_weightings, use_weightings_B, WeightingsB)
  Ifinal <- I
  valid_index_years = ((!is.na(Ifinal)) & (Ifinal != -99))
  cat(sprintf("Number of valid index years: %d (of possible %d)\n", 
              sum(valid_index_years), length(valid_index_years)))
  if (CI_FLAG == 1) {
    t1 <- proc.time() - ptm
    cat(sprintf("[Calculating CIs...] System: %f, User: %f, Elapsed: %f\n", 
                t1[1], t1[2], t1[3]))
    BootI = matrix(0, BOOT_STRAP_SIZE, DSize)
    BootIFlag = matrix(0, 1, BOOT_STRAP_SIZE)
    BootI <- foreach::foreach(Loop = 1:BOOT_STRAP_SIZE) %op% 
      {
        bootstrap_lpi(SpeciesLambdaArray, fileindex, 
                      DSize, Group, Weightings, use_weightings, 
                      use_weightings_B, WeightingsB, CAP_LAMBDAS)
      }
    cat("\n")
    t1 <- proc.time() - ptm
    cat(sprintf("[CIs calculated] System: %f, User: %f, Elapsed: %f\n", 
                t1[1], t1[2], t1[3]))
    BootI <- do.call(cbind, BootI)
    BootI <- t(BootI)
    saveRDS(BootI, here::here(paste0("outputs/rlpi/", scenario_name, "_bootstrap_rlpi.rds")))
    CIx = matrix(0, DSize, 2)
    CIx[1, 1] = 1
    CIx[1, 2] = 1
    for (J in 2:DSize) {
      if (valid_index_years[J]) {
        BootIVal = BootI[, J]
        CIx[J, 1] = quantile(BootIVal, 0.025, names = FALSE)
        CIx[J, 2] = quantile(BootIVal, 0.975, names = FALSE)
      }
      else {
        CIx[J, 1] = NA
        CIx[J, 2] = NA
      }
    }
  }
  if (LEV_FLAG == 1) {
    t1 <- proc.time() - ptm
    cat(sprintf("[Calculating species leverages...] System: %f, User: %f, Elapsed: %f\n", 
                t1[1], t1[2], t1[3]))
    leverage_results = list()
    leverage_diff = list()
    leverage_species = list()
    overall_lambdas = calc_lambdas(Ifinal)
    for (i in 1:nrow(SpeciesLambdaArray)) {
      lev_I = calc_leverage_lpi(SpeciesLambdaArray[-i, 
      ], fileindex, DSize, Group, Weightings, use_weightings, 
      use_weightings_B, WeightingsB)
      leverage_results[[i]] = lev_I
      leverage_diff[[i]] = calc_lambdas(lev_I) - overall_lambdas
      leverage_species[[i]] = SpeciesNamesArray[i, ]
    }
    cat("\n")
    t1 <- proc.time() - ptm
    cat(sprintf("[Species leverages calculated] System: %f, User: %f, Elapsed: %f\n", 
                t1[1], t1[2], t1[3]))
    leverage_results <- do.call(cbind, leverage_results)
    leverage_results <- t(leverage_results)
    leverage_diff <- do.call(cbind, leverage_diff)
    leverage_diff <- t(leverage_diff)
    leverage_results_table <- data.frame(leverage_results)
    colnames(leverage_results_table) <- seq(REF_YEAR, REF_YEAR + 
                                              DSize - 1)
    leverage_results_table$id <- unlist(leverage_species)
    write.csv(leverage_results_table, file = file.path(basedir, 
                                                       "species_leverage_lpi_results.csv"))
    leverage_diff_table <- data.frame(leverage_diff)
    colnames(leverage_diff_table) <- seq(REF_YEAR, REF_YEAR + 
                                           DSize - 1)
    leverage_diff_table$total <- rowSums(leverage_diff_table)
    leverage_diff_table$id <- unlist(leverage_species)
    write.csv(leverage_diff_table, file = file.path(basedir, 
                                                    "species_leverage_diff_lambdas_results.csv"))
  }
  if (SWITCH_PT_FLAG == 1) {
    t1 <- proc.time() - ptm
    cat(sprintf("[Calculating Switch Points...] System: %f, User: %f, Elapsed: %f\n", 
                t1[1], t1[2], t1[3]))
    sp_prog <- txtProgressBar(min = 0, max = BOOT_STRAP_SIZE, 
                              char = "*", style = 3)
    SecondDerivBoot = matrix(0, BOOT_STRAP_SIZE, DSize)
    for (Loop in 1:BOOT_STRAP_SIZE) {
      DTempArrayTemp <- matrix(data = NA, nrow = NoFiles, 
                               ncol = DSize)
      DTempArray = data.frame(DTempArrayTemp)
      for (FileNo in 1:NoFiles) {
        SpeciesLambda = SpeciesLambdaArray[fileindex == 
                                             FileNo, ]
        n = length(SpeciesLambda[, 1])
        BootIndex = 1:n
        BootSam <- sample(BootIndex, replace = T)
        DTemp = matrix(0, 1, dim(SpeciesLambda)[2])
        for (LoopI in 1:dim(SpeciesLambda)[2]) {
          SpeciesLambdaVal = SpeciesLambda[, LoopI]
          BootVal = SpeciesLambdaVal[BootSam]
          Index = which(BootVal != -1)
          if (length(Index) > 0) {
            DTemp[LoopI] = mean(BootVal[Index])
          }
          else DTemp[LoopI] = -99
        }
        if (dim(SpeciesLambda)[2] > DSize) 
          DSize = dim(SpeciesLambda)[2]
        DTempArray[FileNo, 1:dim(t(DTemp))[1]] = t(DTemp)
      }
      I = calculate_index(DTempArray, fileindex, DSize, 
                          Group, Weightings)
      h = 1
      d = 6
      interval = 1
      SecDeriv = CalcSDev(I, h, d, interval)
      SecondDerivBoot[Loop, ] = SecDeriv
      setTxtProgressBar(sp_prog, Loop)
    }
    close(sp_prog)
    CI = matrix(0, DSize, 2)
    for (J in 1:DSize) {
      SecondDerivBootVal = SecondDerivBoot[, J]
      Index = which(SecondDerivBootVal != -1)
      SecondDerivBootVal = SecondDerivBoot[Index, J]
      CI[J, 1] = quantile(SecondDerivBootVal, 0.025, names = FALSE)
      CI[J, 2] = quantile(SecondDerivBootVal, 0.975, names = FALSE)
    }
    SwitchingPt = matrix(0, 1, DSize)
    for (J in 1:DSize) {
      if ((CI[J, 1] > 0) & (CI[J, 2] > 0)) 
        SwitchingPt[J] = 1
      if ((CI[J, 1] < 0) & (CI[J, 2] < 0)) 
        SwitchingPt[J] = -1
    }
  }
  CI2 <- data.frame(CIx)
  lowerCI <- t(CI2$X1)
  upperCI <- t(CI2$X2)
  if (plot_lpi) {
    if (CI_FLAG == 1) {
      plot_lpi(Ifinal, REF_YEAR, PLOT_MAX, CI_FLAG, lowerCI, 
               upperCI)
    }
    else {
      plot_lpi(Ifinal, REF_YEAR, PLOT_MAX)
    }
    if (nchar(title) > 0) {
      title <- paste("[", title, "] ", sep = "")
    }
    if (CI_FLAG == 1) {
      title(paste(title, "Calculated Index; Bootstraps = ", 
                  BOOT_STRAP_SIZE, sep = ""))
    }
    else {
      title(paste(title, "Calculated Index"))
    }
  }
  if (SWITCH_PT_FLAG) {
    LPIdaplta <- cbind(Ifinal, CI2, t(SwitchingPt))
    colnames(LPIdata) <- c("LPI_final", "CI_low", "CI_high", 
                           "SwitchPoint")
  }
  else if (CI_FLAG) {
    LPIdata <- cbind(Ifinal, CI2)
    colnames(LPIdata) <- c("LPI_final", "CI_low", "CI_high")
  }
  else {
    LPIdata <- Ifinal
    colnames(LPIdata) <- c("LPI_final")
  }
  rownames(LPIdata) <- seq(REF_YEAR, REF_YEAR + DSize - 1)
  f_name = file = file.path(basedir, gsub(".txt", "_Results.txt", 
                                          infile))
  cat("Saving final output to file: ", f_name, "\n")
  write.table(LPIdata, f_name)
  FileTable = read.table(infile, header = TRUE)
  FileNames = FileTable$FileName
  Group = FileTable[2]
  NoGroups = length(unique(Group))
  NoFiles = max(dim(Group))
  for (FileNo in 1:NoFiles) {
    Dataset <- toString(FileNames[FileNo])
    Data <- read.table(Dataset, header = TRUE)
    colnames(Data) <- c("Binomial", "ID", "year", "popvalue")
    year <- NULL
    minmax <- plyr::ddply(Data, "ID", plyr::summarise, min_year = min(year), 
                          max_year = max(year))
    f_name = file = file.path(basedir, gsub(".txt", "_Minmax.txt", 
                                            Dataset))
    cat("Saving Min/Max file to: ", f_name, "\n")
    write.table(minmax, sep = ",", eol = "\n", f_name, quote = FALSE, 
                append = FALSE, row.names = F, col.names = T)
  }
  if (save_plots) {
    output_file <- file.path(basedir, gsub(".txt", ".pdf", 
                                           infile))
    cat("Saving Plot to PDF: ", output_file, "\n")
    dev.copy(pdf, output_file)
    dev.off()
  }
  t1 <- proc.time() - ptm
  cat(sprintf("[END] System: %f, User: %f, Elapsed: %f\n", 
              t1[1], t1[2], t1[3]))
  return(LPIdata)
}
