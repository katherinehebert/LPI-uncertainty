## use Hmsc to model lambda vs. year with latent factor per site 

# MODELING growth rate vs. lag-1 growth rate ----

## prepare data for model construction

# make one data frame (which includes a column of growth rates)
ex_l_all <- bind_rows(ex_l, .id = "ID")

# create year x species matrix  
Y = ex_l_all %>% 
  subset(select = c(ID, year, Location, dt)) %>%
  pivot_wider(names_from = ID, values_from = dt)
# remove columns with at least one NA
Y <- Y[ , colSums(is.na(Y)) == 0]

# create XData dataframe of previous year growth rates
temp <- subset(Y, select = -c(year,Location)) %>% apply(2, lag)
temp[1,] <- 0 # assign 0 to first year as a baseline
XData = data.frame(x = temp)


# RANDOM EFFECTS ####

# extract study design (i.e. siteID)
studyDesign = data.frame(Year = factor(sort(Y$year)))

# remove unneeded columns from Y
Y = subset(Y, select = -c(Location, year))
# convert to numeric 
Y = apply(Y, 2, as.numeric)

# settings for posterior distribution sampling
nChains = 2          # run 2 independent MCMC chains
thin = 1             # record every step of the iterations
samples = 100#0       # obtain 1000 samples 
transient = 5*thin # 500*thin usually - ignore iterations before the 2500th
verbose = 5*thin # 500*thin usually

# create random effect at level of UNIQUE SITE
rL = HmscRandomLevel(units = studyDesign$Year)

# create model object with random effect
m = Hmsc(Y = Y, 
         XData = XData,
         studyDesign = studyDesign, 
         ranLevels = list(Year = rL), XScale = FALSE)
# estimate parameters
m = sampleMcmc(m, thin = thin, samples = samples, transient = transient,
               nChains = nChains, verbose = verbose)
# save model object
save(m, file = "~/Documents/GitHub/LPI-sensitivity/DFOexample/hmsc.rds")


# evaluate species associations - which ones have significant associations?

# look at the estimated species-to-species associations
OmegaCor = computeAssociations(m) # extract estimated associations from model object
heatmap(OmegaCor[[1]]$mean, Rowv = NA, Colv = NA)
# note: computeAssociations() also converts covariances to -1 to 1

# choose to plot associations for which the posterior probability for being negative or positive is at least 0.95
supportLevel = 0.95
toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean
# plot all
corrplot::corrplot(OmegaCor[[1]]$mean, 
                   method = "color",
                   col = colorRampPalette(c("blue","white","red"))(200),
                   title = paste("random effect level:", m$rLNames[1]),
                   mar = c(0,0,1,0))

# plot chord diagram of species associations
mat = OmegaCor[[1]]$mean
diag(mat) = NA
#pal = RColorBrewer::brewer.pal(ncol(Y), "Set2")
pal = rainbow(ncol(Y))

# only the significant ones!
chordDiagram(toPlot, annotationTrack = "grid", preAllocateTracks = 1, 
             grid.col = pal)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", 
              niceFacing = TRUE, adj = c(0, 0.5), col = "grey")}, 
  bg.border = NA)    

# predict values from estimated model paramets
preds = Hmsc::computePredictedValues(m)
# evaluate model fit
MF = evaluateModelFit(hM=m, predY=preds)

# get estimate of beta parameter
postBeta = getPostEstimate(m, parName = "Beta")
# plot significant parameters
par(mar = c(6,6,3,2))
plotBeta(m, postBeta, param = "Support", supportLevel = .05)

# get parameter estimates
mpost = convertToCodaObject(m)
summary(mpost$Beta[[1]])

# biplot to visualize the responses of the species to the latent variables alone
rL$nfMin=2
rL$nfMax=2
# null model
m.0 = Hmsc(Y = Y, XData = XData, XFormula = ~1, 
           studyDesign = studyDesign, ranLevels = list(Year=rL))
m.0 = sampleMcmc(m.0, thin = thin, samples = samples, transient = transient,
                 nChains = nChains, verbose = verbose)
# extract posterior means of site and species loadings
etaPost = getPostEstimate(m.0, "Eta") # site loadings (Eta - fancy n)
lambdaPost = getPostEstimate(m.0, "Lambda") # species loadings (lambda)

# biplot
colVar = "Year" # variable to colour point by
par(mar = c(5,4,4,2) + 0.1)
biPlot(m.0, etaPost = etaPost, 
       lambdaPost = lambdaPost, 
       factors = c(1,2), 
       #colors = colorRampPalette(c("blue","white","red"))(10), 
       cex = 3)
# circles : site loadings
# triangles: species loadings


# MODELING growth rate vs. time ------------------------------------------------

# create year x species matrix  
Y2 = ex_l_all %>% 
  subset(select = c(ID, year, Location, dt)) %>%
  pivot_wider(names_from = ID, values_from = dt)
# remove columns with at least one NA
Y2 <- Y2[ , colSums(is.na(Y2)) == 0]

# create XData dataframe of time
XData2 = data.frame(x = Y2$year)

# RANDOM EFFECTS ####

# extract study design (i.e. siteID)
studyDesign2 = data.frame(Year = factor(sort(Y2$year)))

# remove unneeded columns from Y
Y2 = subset(Y2, select = -c(Location, year))
# convert to numeric 
Y2 = apply(Y2, 2, as.numeric)

# create random effect at level of UNIQUE SITE
rL2 = HmscRandomLevel(units = studyDesign2$Year)

# create model object with random effect
m2 = Hmsc(Y = Y2, XFormula = ~ x + 1,
          XData = XData2,
          studyDesign = studyDesign2, 
          ranLevels = list(Year = rL2), XScale = FALSE)
# estimate parameters
m2 = sampleMcmc(m2, thin = thin, samples = samples, transient = transient,
                nChains = nChains, verbose = verbose)
# save model object
save(m2, file = "~/Documents/GitHub/LPI-sensitivity/DFOexample/hmsc2.rds")

# evaluate species associations - which ones have significant associations?

# look at the estimated species-to-species associations
OmegaCor2 = computeAssociations(m2) # extract estimated associations from model object
heatmap(OmegaCor2[[1]]$mean, Rowv = NA, Colv = NA)
# note: computeAssociations() also converts covariances to -1 to 1

# choose to plot associations for which the posterior probability for being negative or positive is at least 0.95
toPlot2 = ((OmegaCor2[[1]]$support>supportLevel)
           + (OmegaCor2[[1]]$support<(1-supportLevel))>0)*OmegaCor2[[1]]$mean
# plot all
corrplot::corrplot(OmegaCor2[[1]]$mean, 
                   method = "color",
                   col = colorRampPalette(c("blue","white","red"))(200),
                   title = paste("random effect level:", m2$rLNames[1]),
                   mar = c(0,0,1,0))

# plot chord diagram of species associations
mat2 = OmegaCor2[[1]]$mean
diag(mat2) = NA

# plot the significant ones!
chordDiagram(toPlot2, annotationTrack = "grid", preAllocateTracks = 1, 
             grid.col = pal)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", 
              niceFacing = TRUE, adj = c(0, 0.5), col = "grey")}, 
  bg.border = NA)    

# predict values from estimated model paramets
preds2 = Hmsc::computePredictedValues(m2)
# evaluate model fit
MF2 = evaluateModelFit(hM=m2, predY=preds2)

# get estimate of beta parameter
postBeta2 = getPostEstimate(m2, parName = "Beta")
# plot significant parameters
par(mar = c(6,6,3,2))
plotBeta(m2, postBeta2, param = "Support", supportLevel = .05)

# biplot to visualize the responses of the species to the latent variables alone
rL2$nfMin=2
rL2$nfMax=2
# null model
m.02 = Hmsc(Y = Y2, XData = XData2, XFormula = ~1, 
            studyDesign = studyDesign2, ranLevels = list(Year=rL2))
m.02 = sampleMcmc(m.02, thin = thin, samples = samples, transient = transient,
                  nChains = nChains, verbose = verbose)
# extract posterior means of site and species loadings
etaPost2 = getPostEstimate(m.02, "Eta") # site loadings (Eta - fancy n)
lambdaPost2 = getPostEstimate(m.02, "Lambda") # species loadings (lambda)

