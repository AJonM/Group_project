library(ape)
library(picante)
library(letsR)
library(geiger)
library(phytools)
library(dismo)
library(XML)
library(maptools)
library(sp)
library(raster)
library(foreign)
library(rgdal)
library(rgeos)
library(readxl)
library(geosphere)

# Functions and projections
'%!in%' <- function(x,y)!('%in%'(x,y))
behrmann.proj <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +ellps=WGS84 +units=km")
azimuthal.proj <-CRS("+proj=aeqd +lat_0=0 +lon_0=-0")
latlong = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")


#############################################################
################# List species per island ###################
#############################################################
# Load IUCN shapefile (2019)
mammals_IUCN <- readOGR("TERRESTRIAL_MAMMALS.shp","TERRESTRIAL_MAMMALS", dropNULLGeometries=T)

# Remove species that were introduced, not mapped or with presence uncertain
remove <- c("Extant & Introduced (resident)", "Extinct & Introduced", "Presence Uncertain", 
            "Introduced", "Probably Extant & Introduced (resident)", "Not Mapped", "Presence Uncertain & Vagrant")
to_remove <- which(mammals_IUCN$legend %in% remove)
if(length(to_remove) != 0){ mammals_IUCN <- mammals_IUCN[-which(mammals_IUCN$legend %in% remove),]}

# List bat species
bat_species <- as.character(unique(mammals_IUCN[which(mammals_IUCN$order_name %in% "CHIROPTERA"),]$binomial))
bat_species <- c(bat_species, "Pteropus brunneus", "Pteropus pilosus", "Pteropus subniger", "Pteropus tokudae") # for some reason some bat species are not listed. Add them manually

# Load island shapefile (GADM 3.6 edited in QGIS to remove continental areas)
island_polygons_edited <- readOGR("gadm36_edited.shp", layer = "gadm36_edited", dropNULLGeometries=T)
island_polygons_edited <- spTransform(island_polygons_edited, crs(mammals_IUCN)) 

# Load island data from Weigelt et al. 2013 - PNAS. Spatialized the points by using the coordinates from the centroids informed in Weigelt's spreadsheet from the sup material
PNAS_island_edited <- readOGR("Weigelt_PNAS_spatial_points.shp", layer = "Weigelt_PNAS_spatial_points", dropNULLGeometries=T)
PNAS_island_edited <- spTransform(PNAS_island_edited, crs(mammals_IUCN)) 

# Overlap island polygons with Weigelt's data points
island_matched_PNAS <- over(PNAS_island_edited, island_polygons_edited, returnList  = F)

# Add polygons ID to PNAS data (ID_gadm)
PNAS_island_edited$ID_gadm <- island_matched_PNAS$id

# Subset island polygons to remain only those with data
island_polygons_edited <- subset(island_polygons_edited, island_polygons_edited$id %in% island_matched_PNAS$id)

# Calculate species presence per island
# !!! This step was only done after visually inspecting in QGIS if island shape and species polygon matched!!!
{
  island_poly <- list()
  for(i in 1:length(island_polygons_edited)){
    island_poly[[i]] <- island_polygons_edited[i,]
    print(i)
  }
  
  get_spp_per_island <- function(island_poly, IUCN_poly){
    library(sp)
    resu <- as.character(unique(IUCN_poly[island_poly,]$binomial))
    return(resu)
  }
  
  # Paralelize and run in the cluster
  source("https://raw.githubusercontent.com/csdambros/R-functions/master/passSSHid.R")
  passSSHid() 
  library (snow) 
  cl <- makeMPIcluster(count = 20)
  {
    stops <- c(1, 10, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11024) # add breaks in case something goes wrong
    for(i in 1:(length(stops)-1)){
      resu <- clusterApply(cl = cl, x = island_poly[c(stops[i]:stops[i+1])], fun = get_spp_per_island,  IUCN_poly = mammals_IUCN)
      save(resu, file = paste("resu_overlap_island_IUCN", stops[i], stops[i+1], sep = "_"))
    }
    stopCluster(cl)
  }
  
  island_over_IUCN <- list()
  stops <- c(1, 10, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11024)
  for(i in 1:(length(stops)-1)){
    load(paste("resu_overlap_island_IUCN", stops[i], stops[i+1], sep = "_"))
    island_over_IUCN[c(stops[i]:stops[i+1])] <- resu
  }
  
  names(island_over_IUCN) <- island_polygons_edited$id
  save(island_over_IUCN, file = "island_over_IUCN_full_2019")
}

# Build presence and absence matrix
spp_per_island <- matrix(0, nrow = length(unique(na.omit(PNAS_island_edited$id_spatial_island))), ncol = length(unique(mammals_IUCN@data$binomial)))
rownames(spp_per_island) = unique(na.omit(PNAS_island_edited$id_spatial_island))
colnames(spp_per_island) = unique(mammals_IUCN@data$binomial)
# Add presence to species that occur in each island
for(i in 1:length(island_over_IUCN)){
  if(length(island_over_IUCN[[i]]) != 0){
    spp_per_island[names(island_over_IUCN[i]), island_over_IUCN[[i]]] <- 1      
  }
}

## Add extinct species
  # Phylacine
endemicity <- read.csv("Phylacine_trait_data.csv", sep = ",", na.strings = "NA", stringsAsFactors = F)
endemicity$Binomial.1.2 <- gsub("_", " ", endemicity$Binomial.1.2, fixed=TRUE)
endemicity <- endemicity[endemicity$IUCN.Status.1.2 == "EX",]
endemicity <- endemicity[endemicity$Marine == 0,]
endemicity <- endemicity[endemicity$Island.Endemicity %in% c("Occurs only on isolated islands", "Occurs on small land bridge islands", "Occurs on large land bridge islands"),]
  # Get information on which island they lived in - Faurby & Svenning (2016)
extinct_faurby <- read.csv("faurby_with_island.csv", sep = ";", na.strings = "NA", stringsAsFactors = F)
extinct_faurby <- extinct_faurby[!is.na(extinct_faurby$PNAS.ID.NA),] # get only data to which we know the island
for(i in 1:nrow(extinct_faurby)){
  lines <- !is.na(extinct_faurby[i,c(1:ncol(extinct_faurby))])
  temporary <- extinct_faurby[i,lines]
 if(i==1){list_extinct_spp <- temporary[1,c(5:ncol(temporary))]}
 if(i==2){list_extinct_spp <- unlist(c(list_extinct_spp,temporary[1,c(5:ncol(temporary))]), use.names=FALSE)}
 if(sum(i!=c(1,2))>=2){
     list_extinct_spp <- c(list_extinct_spp,temporary[1,c(5:ncol(temporary))])
     if(is.list(list_extinct_spp)){
       list_extinct_spp <- unlist(list_extinct_spp, use.names=FALSE)
     }
   }
 }
  # Retrieve which island they belong to
for(i in 1:length(extinct_faurby$PNAS.ID)){
 if(length(which(PNAS_island_subset$ID %in% extinct_faurby$PNAS.ID[i]))!=0){
   extinct_faurby$PNAS_id_spatial[i] <- PNAS_island_subset$id_spatial_island[which(PNAS_island_subset$ID %in% extinct_faurby$PNAS.ID[i])]
   extinct_faurby$PNAS_id_spatial_island[i] <- as.character(PNAS_island_subset$Island[which(PNAS_island_subset$ID %in% extinct_faurby$PNAS.ID[i])])
  }
  else{
    extinct_faurby$PNAS_id_spatial[i] <- NA
    extinct_faurby$PNAS_id_spatial_island[i] <- NA
   }
 }
for(i in 1:length(endemicity$Binomial.1.2)){ # for each extinct species on island...
  location <- which(extinct_faurby == endemicity$Binomial.1.2[i], arr.ind = T) # get data on which island it lived in
  loc_row <- location[1,1]
  loc_col <- location[1,2]
  if(!is.na(extinct_faurby$PNAS_id_spatial[loc_row])){
    row_position <- which(rownames(spp_per_island) %in% extinct_faurby$PNAS_id_spatial[loc_row]) # which row in presab matrix
    col_position <- which(colnames(spp_per_island) %in% endemicity$Binomial.1.2[i])
    if(length(col_position) != 0){ # if we already have this species on the presab matrix...
      spp_per_island[row_position,endemicity$Binomial.1.2[i]] <- 1
    }
    if(length(col_position) == 0){ # if not
      spp_per_island[,endemicity$Binomial.1.2[i]] <- 0
      spp_per_island[row_position,endemicity$Binomial.1.2[i]] <- 1
    }
  }
}

## Summarize richness and SIE per island
# Put island data in the same order as in the presab object
PNAS_island_subset <- PNAS_island_subset[order(PNAS_island_subset$id_spatial_island),]
spp_per_island <- spp_per_island[order(as.numeric(rownames(spp_per_island))),]

# Calculate SR and save in island data
PNAS_island_subset$Richness_mammal <- rowSums(spp_per_island)
PNAS_island_subset$Richness_bat <- rowSums(spp_per_island[,which(colnames(spp_per_island) %in% bat_species)])
PNAS_island_subset$Richness_nonVol <- PNAS_island_subset$Richness_mammal  - PNAS_island_subset$Richness_bat

# List 100 richest islands
richest_islands <- rownames(spp_per_island)[order(rowSums(spp_per_island), decreasing = T)][c(1:100)]
rowSums(spp_per_island)[order(rowSums(spp_per_island), decreasing = T)[c(1:100)]] # check richness
PNAS_island_edited[which(PNAS_island_edited$id_spatial_island %in% richest_islands),]@data[,c("Island", "Area", "id_spatial_island")]

# List 100 most frequent species on islands
frequency_on_island <- colnames(spp_per_island)[order(colSums(spp_per_island), decreasing = T)][c(1:100)]
colSums(spp_per_island)[order(colSums(spp_per_island), decreasing = T)[c(1:100)]]

# Calculate SIE
# Which species are present in only one island?
SIE <- names(which(colSums(spp_per_island) == 1))
length(SIE)
# Exclude those that occur in the continent
mammals_IUCN <- readOGR("D:/Onedrive/Documentos/Academia/Doutorado/Tese/Banco_de_dados/Filogenias_e_distribuicao/Dados_distribuicao/IUCN_july_2018/TERRESTRIAL_MAMMALS.shp","TERRESTRIAL_MAMMALS", dropNULLGeometries=T) # Load shapefile
# Remove introduced species
continents <- readOGR("IUCN_basemap_continents.shp", layer = "IUCN_basemap_continents", dropNULLGeometries=T)
continents <- spTransform(continents, crs(mammals_IUCN)) 
spp_in_continents <- vector()
for(i in 1:length(continents)){
   spp_in_continents <- c(spp_in_continents, as.character(unique(mammals_IUCN[continents[i,],]$binomial)))
   print(i)
}

SIE <- SIE[which(SIE %!in% spp_in_continents)]

# Add SIE data to spatial data
for(i in 1:nrow(spp_per_island)){
  temp <- which(colnames(spp_per_island)[which(spp_per_island[i,]>0)] %in% SIE) 
  if(length(temp) != 0){
    PNAS_island_subset@data[which(PNAS_island_subset@data$id_spatial_island %in% rownames(spp_per_island)[i]),"SIE_total"] <- length(temp)
  }
  else{
    PNAS_island_subset@data[which(PNAS_island_subset@data$id_spatial_island %in% rownames(spp_per_island)[i]),"SIE_total"] <- 0
  }
  
  temp2 <- which(colnames(spp_per_island[,which(colnames(spp_per_island) %in% bat_species)])[which(spp_per_island[,which(colnames(spp_per_island) %in% bat_species)][i,]>0)] %in% SIE) # how many bat species in that island are SIE?
  if(length(temp2) != 0){
    PNAS_island_subset@data[which(PNAS_island_subset@data$id_spatial_island %in% rownames(spp_per_island[,which(colnames(spp_per_island) %in% bat_species)])[i]),"SIE_bats"] <- length(temp2)
  }
  else{
    PNAS_island_subset@data[which(PNAS_island_subset@data$id_spatial_island %in% rownames(spp_per_island)[i]),"SIE_bats"] <- 0
  }
}

# Proportion of SIE
PNAS_island_subset$pSIE <- PNAS_island_subset$SIE_total/PNAS_island_subset$Richness_mammal
PNAS_island_subset$pSIE_bats <- PNAS_island_subset$SIE_bats/PNAS_island_subset$Richness_bat
PNAS_island_subset$pSIE_nonVol <- PNAS_island_subset$SIE_nonVol/PNAS_island_subset$Richness_nonVol
PNAS_island_subset$SIE_nonVol <- PNAS_island_subset$SIE_total - PNAS_island_subset$SIE_bats

PNAS_island_subset@data[is.na(PNAS_island_subset@data$pSIE),"pSIE"] <- 0
PNAS_island_subset@data[is.na(PNAS_island_subset@data$pSIE_bats),"pSIE_bats"] <- 0
PNAS_island_subset@data[is.na(PNAS_island_subset@data$pSIE_nonVol),"pSIE_nonVol"] <- 0


###################################################################
#################### ANALYSIS! ####################################
###################################################################
library(ape)
library(picante)
library(letsR)
library(geiger)
library(phytools)
library(dismo)
library(XML)
library(maptools)
library(sp)
library(raster)
library(foreign)
library(rgdal)
library(rgeos)
library(readxl)
library(geosphere)
library(dplyr)
library(sjPlot)
library(ggplot2)
library(corrplot)
library(tidyverse)
library(stargazer)
library(ncf)
library(glmmTMB)
library(TMB)
library(effects)
library(arm)
library(ggpubr)
library(DHARMa)
library(ggpubr)
library(lme4)
library(visreg)
library(MASS)

island_data = read.table("Appendix_1.xlsx")
island_data_endemism = subset(island_data, island_data$SIE_mammal > 0)

## Check multicolinearity and relationship among predictors
predictor = c("Area", "Current_isolation",	"Past_isolation",	"Climate_velocity",	"Temperature_mean",	"Temperature_sd",	"Precipitation_mean",
              "Precipitation_sd",	"Elevation_sd",	"bioregion",	"bioregion_SIE")
# All islands
island_data[,predictor] = apply(island_data[,predictor], 2, as.numeric)
print(corrplot::corrplot(cor(island_data[,predictor]),
                        method = "color", type = "upper", order = "alphabet", number.cex = .7, addCoef.col = "black", tl.col = "black", 
                        tl.srt = 90, diag = FALSE,  title = "all islands", mar=c(0,0,1,0)))
dev.off()
# Only endemics
png(height=1200, width=1500, pointsize=15, file="multicollinearity_endemism.png")
island_data_endemism[,predictor] = apply(island_data_endemism[,predictor], 2, as.numeric)
print(corrplot::corrplot(cor(island_data_endemism[,predictor]),
                         method = "color", type = "upper", order = "alphabet", number.cex = .7, addCoef.col = "black", tl.col = "black", 
                         tl.srt = 90, diag = FALSE,  title = "islands with endemics", mar=c(0,0,1,0)))
dev.off()


### Fit models
# Scale the predictors
island_data[,predictor] = scale(island_data[,predictor])

# Species richness
island_data$Richness_mammal_with0 <- island_data$Richness_mammal-1
final_model_glmmTMB_Richness_mammal = glmmTMB(Richness_mammal~Area*Current_isolation + Temperature_mean + Temperature_sd + Past_isolation + Climate_velocity +
                                       Precipitation_sd + Precipitation_mean + Elevation_sd + (1|bioregion) +
                                       (0+Area*Current_isolation|bioregion) + (0+Temperature_mean|bioregion) + (0+Temperature_sd|bioregion) +
                                       (0+Past_isolation|bioregion) + (0+Climate_velocity|bioregion) + (0+Precipitation_sd|bioregion) + (0+Precipitation_mean|bioregion) +
                                       (0+Elevation_sd|bioregion), data = island_data, family = truncated_nbinom2(link = "log"))
final_model_glmmTMB_Richness_bat = glmmTMB(Richness_bat~Area*Current_isolation + Temperature_mean + Temperature_sd + Past_isolation + Climate_velocity +
                                        Precipitation_sd + Precipitation_mean + Elevation_sd + (1|bioregion) +
                                        (0+Area*Current_isolation|bioregion) + (0+Temperature_mean|bioregion) + (0+Temperature_sd|bioregion) +
                                        (0+Past_isolation|bioregion) + (0+Climate_velocity|bioregion) + (0+Precipitation_sd|bioregion) + (0+Precipitation_mean|bioregion) +
                                        (0+Elevation_sd|bioregion), data = island_data, family = nbinom2)
final_model_glmmTMB_Richness_nonVol = glmmTMB(Richness_nonVol~Area*Current_isolation + Temperature_mean + Temperature_sd + Past_isolation + Climate_velocity +
                                          Precipitation_sd + Precipitation_mean + Elevation_sd + (1|bioregion) +
                                          (0+Area*Current_isolation|bioregion) + (0+Temperature_mean|bioregion) + (0+Temperature_sd|bioregion) +
                                          (0+Past_isolation|bioregion) + (0+Climate_velocity|bioregion) + (0+Precipitation_sd|bioregion) + (0+Precipitation_mean|bioregion) +
                                          (0+Elevation_sd|bioregion), data = island_data, family = nbinom2)

# Single Island Endemics (SIE)
island_data_endemism$SIE_total_with0 <- island_data_endemism$SIE_total-1
island_data_endemism[,predictor] = scale(island_data_endemism[,predictor])

final_model_SIE_total = glm.nb(SIE_total_with0~Area+Current_isolation + Temperature_mean + Temperature_sd + Past_isolation + Climate_velocity +
                                 Precipitation_sd + Precipitation_mean + Elevation_sd, data = island_data_endemism)
final_model_SIE_total = glm.nb(SIE_total_with0~Area+Current_isolation + Temperature_mean + Temperature_sd + Past_isolation + Climate_velocity +
                                 Precipitation_sd + Precipitation_mean + Elevation_sd, data = island_data_endemism, init.theta = final_model_SIE_total$theta)

final_model_SIE_bats = glm.nb(SIE_bats~Area+Current_isolation + Temperature_mean + Temperature_sd + Past_isolation + Climate_velocity +
                                Precipitation_sd + Precipitation_mean + Elevation_sd, data = island_data_endemism)

final_model_SIE_nonVol = glm.nb(SIE_nonVol~Area+Current_isolation + Temperature_mean + Temperature_sd + Past_isolation + Climate_velocity +
                                  Precipitation_sd + Precipitation_mean + Elevation_sd, data = island_data_endemism)
final_model_SIE_nonVol = glm.nb(SIE_total_with0~Area+Current_isolation + Temperature_mean + Temperature_sd + Past_isolation + Climate_velocity +
                                 Precipitation_sd + Precipitation_mean + Elevation_sd, 
                               data = island_data_endemism, init.theta = final_model_SIE_nonVol$theta)

plot_models(final_model_SIE_total, final_model_SIE_bats,final_model_SIE_nonVol)


# pSIE
final_model_pSIE_total = glmer(pSIE_total~Area+Current_isolation + Temperature_mean + Temperature_sd + Past_isolation + Climate_velocity +
                                         Precipitation_sd + Precipitation_mean + Elevation_sd + (1|bioregion_SIE) + (0+Current_isolation|bioregion_SIE) +
                                         (0+Area|bioregion_SIE) + (0+Temperature_mean|bioregion_SIE) + (0+Temperature_sd|bioregion_SIE) +
                                         (0+Past_isolation|bioregion_SIE) + (0+Climate_velocity|bioregion_SIE) + (0+Precipitation_sd|bioregion_SIE) + (0+Precipitation_mean|bioregion_SIE) +
                                         (0+Elevation_sd|bioregion_SIE), data = island_data_endemism, 
                                          family = binomial(link = "logit"), weights= Richness_mammal)
ss <- getME(final_model_pSIE_total,c("theta","fixef"))
final_model_pSIE_total <- update(final_model_pSIE_total,start=ss)
summary(final_model_pSIE_total)

final_model_pSIE_bats = glmer(pSIE_bats~Area+Current_isolation + Temperature_mean + Temperature_sd + Past_isolation + Climate_velocity +
                                          Precipitation_sd + Precipitation_mean + Elevation_sd + (1|bioregion_SIE) + (0+Current_isolation|bioregion_SIE) +
                                          (0+Area|bioregion_SIE) + (0+Temperature_mean|bioregion_SIE) + (0+Temperature_sd|bioregion_SIE) +
                                          (0+Past_isolation|bioregion_SIE) + (0+Climate_velocity|bioregion_SIE) + (0+Precipitation_sd|bioregion_SIE) + (0+Precipitation_mean|bioregion_SIE) +
                                          (0+Elevation_sd|bioregion_SIE), data = island_data_endemism, 
                                          family = binomial(link = "logit"), weights= Richness_bat)
summary(final_model_pSIE_bats)

final_model_pSIE_nonVol = glmer(pSIE_nonVol~Area+Current_isolation + Temperature_mean + Temperature_sd + Past_isolation + Climate_velocity +
                                            Precipitation_sd + Precipitation_mean + Elevation_sd + (1|bioregion_SIE) + (0+Current_isolation|bioregion_SIE) +
                                            (0+Area|bioregion_SIE) + (0+Temperature_mean|bioregion_SIE) + (0+Temperature_sd|bioregion_SIE) +
                                            (0+Past_isolation|bioregion_SIE) + (0+Climate_velocity|bioregion_SIE) + (0+Precipitation_sd|bioregion_SIE) + (0+Precipitation_mean|bioregion_SIE) +
                                            (0+Elevation_sd|bioregion_SIE), data = island_data_endemism, 
                                            family = binomial(link = "logit"), weights= Richness_nonVol)
ss <- getME(final_model_pSIE_nonVol,c("theta","fixef"))
final_model_pSIE_nonVol <- update(final_model_pSIE_nonVol,start=ss)
summary(final_model_pSIE_nonVol)

# Check for zero-inflation and for overdispersion
final_model = list(final_model_glmmTMB_Richness_mammal, final_model_glmmTMB_Richness_bat, final_model_glmmTMB_Richness_nonVol,
                    final_model_SIE_total, final_model_SIE_bats, final_model_SIE_nonVol,
                    final_model_pSIE_total, final_model_pSIE_bats, final_model_pSIE_nonVol)

library(DHARMa)
for(i in 1:length(final_model)){
  simulated <- simulateResiduals(fittedModel = final_model[[i]], re.form = NULL)
  hist(simulated$scaledResiduals)
  hist(simulated$scaledResidualsNormal) 
  zeroinfl <- testZeroInflation(simulated)
  zeroinfl
  zeroinflation = F
  if(zeroinfl$p.value <= 0.05){
    original_formula <-  formula(final_model[[i]])
    final_model[[i]] <- glmmTMB(original_formula, data = island_data, family="poisson", ziformula = ~1)    
    simulated <- simulateResiduals(fittedModel = final_model[[i]])
    zeroinflation = T
  }
  overdisp_fun <- function(model) {
    rdf <- df.residual(model)
    rp <- residuals(model)
    Pearson.chisq <- sum(rp^2)
    prat <- Pearson.chisq/rdf
    pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
    c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
  }
  overdisp_fun(final_model[[i]])[4]
  simulated <- simulateResiduals(fittedModel = final_model[[i]], re.form = NULL)
  testDispersion(simulated, alternative = "greater")
  testDispersion(simulated, alternative = "two.sided")
  # Check the residuals
  plot(simulated) # qq-plot to detect overall deviations from the expected distribution
  plot(final_model[[i]])
  # Plot residuals against predictor using simulated residuals
  vars = colnames(final_model[[i]]@frame)[-c(length(colnames(final_model[[i]]@frame)))] # we would expect niformity in y direction if we plot against any predictor
  par(mfrow = c(3, ceiling(length(vars)/3)))
  for(i in 1:length(vars)){
    plotResiduals(island_data[,vars[i]], simulated$scaledResiduals, xlab = vars[i])
  }
  # Plot residuals against predictor using true residuals
  par(mfrow = c(3, ceiling(length(vars)/3)))
  for(i in 1:length(vars)){
    plot(island_data[,vars[i]], residuals(final_model[[i]], type = "pearson"), xlab = vars[i])
  }
  dev.off()
  # Check for outliers
  testOutliers(simulated, alternative = "two.sided")
  check_outliers(final_model[[i]]) 
}

# Calculate R square
#library(MuMIn)
for(i in 1:length(final_model)){
  MuMIn::r.squaredGLMM(final_model[[i]])
}

# Plot all models coefficients
# By SR, SIE, pSIE
SR = plot_models(final_model_glmer_Richness_mammal,
                final_model_glmer_Richness_bat,
                final_model_glmmTMB_Richness_nonVol,
            transform = NULL,
            colors = c("tomato1", "red2", "red4"),
            title = "Species richness", dot.size = 3.5, spacing = 0.5,
            #show.values = T,
            show.p = T, 
            axis.title = "Stantardized beta coefficient",
            axis.labels = c("Area", "Area x SLMP", "CCVT",
                            "Elevation sd", "Past_isolation", "Mean precipitation",
                            "Precipitation sd", "SLMP", "Mean temperature", "Temperature sd"),
            m.labels = c("All species", 
                         "Bats",
                         "Non-volants"))
SR = SR + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = 0, linetype="dotted", colour = "grey", size = 1)  + theme(legend.position=c(0.02,0.15))
pdf("coefficients_SR.pdf", onefile = T, useDingbats=FALSE)
print(SR)
dev.off()

SIE = plot_models(final_model_SIE_total,
                  final_model_SIE_bats,
                  final_model_SIE_nonVol,
                transform = NULL,
                colors = c("lightblue", "blue2", "blue4"),
                title = "Single Island Endemics (SIE)", dot.size = 3.5, spacing = 0.5,
                #show.values = T,
                show.p = T, 
                axis.title = "Stantardized beta coefficient",
                axis.labels = c("Area", "CCVT",
                                "Elevation sd", "Past_isolation", "Mean precipitation",
                                "Precipitation sd", "SLMP", "Mean temperature", "Temperature sd"),
                m.labels = c("All species", 
                             "Bats",
                             "Non-volants"))
SIE = SIE + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = 0, linetype="dotted", colour = "grey", size = 1)  + theme(legend.position=c(0.02,0.15))
pdf("coefficients_SIE.pdf", onefile = T, useDingbats=FALSE)
print(SIE)
dev.off()

pSIE = plot_models(final_model_pSIE_total,
                   final_model_pSIE_bats,
                   final_model_pSIE_nonVol,
                transform = NULL,
                colors = c("green1", "green3", "green4"),
                title = "Proportion of Single Island Endemics (pSIE)", dot.size = 3.5, spacing = 0.5,
                #show.values = T,
                show.p = T, 
                axis.title = "Stantardized beta coefficient",
                axis.labels = c("Area", "CCVT",
                                "Elevation sd", "Past_isolation", "Mean precipitation",
                                "Precipitation sd", "SLMP", "Mean temperature", "Temperature sd"),
                m.labels = c("All species", 
                             "Bats",
                             "Non-volants"))
pSIE = pSIE + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = 0, linetype="dotted", colour = "grey", size = 1)  + theme(legend.position=c(0.02,0.15))
pdf("coefficients_pSIE.pdf", onefile = T, useDingbats=FALSE)
print(pSIE)
dev.off()

pdf("coefficients_SR_SIE_pSIE.pdf", onefile = T, width = 16, useDingbats=FALSE)
ggarrange(SR, SIE, pSIE, nrow = 1, ncol = 3)
dev.off()

## Compare islands with bats and non-flying, only bats and only non-flying
## Species richness
island_data$occupancy <- ifelse(island_data$Richness_bat > 0 & island_data$Richness_nonVol == 0, "Only bats", NA)
island_data$occupancy <- ifelse(island_data$Richness_bat > 0 & island_data$Richness_nonVol > 0, "Both", island_data$occupancy)
island_data$occupancy <- ifelse(island_data$Richness_bat == 0 & island_data$Richness_nonVol > 0, "Only non-flying", island_data$occupancy)
island_data$occupancy <- as.factor(island_data$occupancy)
summary(island_data$occupancy)

# Area
diff_area = lm(Area~occupancy, data = island_data) 
summary(diff_area)
anova(diff_area)
library(multcomp)
PostHocTeste=glht(diff_area,linfct = mcp(occupancy = "Tukey"))
summary(PostHocTeste)
my_comparisons <- list(c("Both", "Only bats"), c("Both", "Only non-flying"), c("Only bats", "Only non-flying"))
p1 <- ggboxplot(island_data, x = "occupancy", y = "Area", ylab = "Area (log)", xlab = "",
              color = "occupancy", palette =c("red3", "blue3", "green3")) +
              stat_compare_means(comparisons = my_comparisons,label = "p.signif") + 
              ggtitle("Area") + theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

# SLMP
diff_slmp = lm(Current_isolation ~occupancy, data = island_data) 
summary(diff_slmp)
anova(diff_slmp)
library(multcomp)
PostHocTeste=glht(diff_slmp,linfct = mcp(occupancy = "Tukey"))
summary(PostHocTeste)
p2 = ggboxplot(island_data, x = "occupancy", y = "Current_isolation", ylab = "SLMP", xlab = "",
               color = "occupancy", palette =c("red3", "blue3", "green3")) +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") + 
  ggtitle("SLMP") + theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

# Temperature
diff_temp = lm(Temperature_mean ~occupancy, data = island_data) 
summary(diff_temp)
anova(diff_temp)
library(multcomp)
PostHocTeste=glht(diff_temp,linfct = mcp(occupancy = "Tukey"))
summary(PostHocTeste)
p3= ggboxplot(island_data, x = "occupancy", y = "Temperature_mean", ylab = "Mean temperature", xlab = "",
              color = "occupancy", palette =c("red3", "blue3", "green3")) +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") + 
  ggtitle("Mean temperature") + theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

# Temp_sd
diff_temp_sd = lm(Temperature_sd ~occupancy, data = island_data) 
summary(diff_temp_sd)
anova(diff_temp_sd)
library(multcomp)
PostHocTeste=glht(diff_temp_sd,linfct = mcp(occupancy = "Tukey"))
summary(PostHocTeste)
p4 = ggboxplot(island_data, x = "occupancy", y = "Temperature_sd", ylab = "Temperature sd", xlab = "",
               color = "occupancy", palette =c("red3", "blue3", "green3")) +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") + 
  ggtitle("Temperature sd") + theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

# Past_isolation 
diff_Past_isolation = lm(Past_isolation ~occupancy, data = island_data) 
summary(diff_Past_isolation)
anova(diff_Past_isolation)
library(multcomp)
PostHocTeste=glht(diff_Past_isolation,linfct = mcp(occupancy = "Tukey"))
summary(PostHocTeste)
island_data$Past_isolation_cat <- ifelse(island_data$Past_isolation == 1, "Connected", island_data$Past_isolation)
island_data$Past_isolation_cat <- ifelse(island_data$Past_isolation_cat == 0, "Disconnected", island_data$Past_isolation_cat)
#ggplot(island_data, aes(occupancy, ..count..)) + geom_bar(aes(fill = Past_isolation_cat), position = "dodge")
p10 = ggplot(island_data, aes(Past_isolation_cat, ..count..)) + geom_bar(aes(fill = occupancy), position = "dodge") + 
  scale_fill_manual(values = c("red3", "blue3", "green3")) + labs(x= "", y="Number of islands") + theme_classic()

# Climate_velocity 
diff_Climate_velocity  = lm(Climate_velocity ~occupancy, data = island_data) 
summary(diff_Climate_velocity )
anova(diff_Climate_velocity )
library(multcomp)
PostHocTeste=glht(diff_Climate_velocity ,linfct = mcp(occupancy = "Tukey"))
summary(PostHocTeste)
p5 = ggboxplot(island_data, x = "occupancy", y = "Climate_velocity", ylab = "CCVT (log)", xlab = "",
               color = "occupancy", palette =c("red3", "blue3", "green3")) +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") + 
  ggtitle("Climate-change velocity") + theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

# Precipitation_sd
diff_prec_sd = lm(Precipitation_sd~occupancy, data = island_data) 
summary(diff_prec_sd)
anova(diff_prec_sd)
library(multcomp)
PostHocTeste=glht(diff_prec_sd,linfct = mcp(occupancy = "Tukey"))
summary(PostHocTeste)
p6 = ggboxplot(island_data, x = "occupancy", y = "Precipitation_sd", ylab = "Precipitation sd (log)", xlab = "",
               color = "occupancy", palette =c("red3", "blue3", "green3")) +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") + 
  ggtitle("Precipitation sd") + theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

# Prec mean
diff_prec_mean = lm(Precipitation_mean ~occupancy, data = island_data) 
summary(diff_prec_mean)
anova(diff_prec_mean)
library(multcomp)
PostHocTeste=glht(diff_prec_mean,linfct = mcp(occupancy = "Tukey"))
summary(PostHocTeste)
p7 = ggboxplot(island_data, x = "occupancy", y = "Precipitation_mean", ylab = "Mean precipitation (sqrt)", xlab = "",
               color = "occupancy", palette =c("red3", "blue3", "green3")) +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") + 
  ggtitle("Mean precipitation") + theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

# Elevation sd
diff_elevation = lm(Elevation_sd ~occupancy, data = island_data) 
summary(diff_elevation)
anova(diff_elevation)
library(multcomp)
PostHocTeste=glht(diff_elevation,linfct = mcp(occupancy = "Tukey"))
summary(PostHocTeste)
p8 = ggboxplot(island_data, x = "occupancy", y = "Elevation_sd", ylab = "Elevation sd (log)", xlab = "",
               color = "occupancy", palette =c("red3", "blue3", "green3")) +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") + 
  ggtitle("Elevation sd") + theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

# Latitude
diff_latitude = lm(Lat_centroid ~occupancy, data = island_data) 
summary(diff_latitude)
anova(diff_latitude)
library(multcomp)
PostHocTeste=glht(diff_latitude,linfct = mcp(occupancy = "Tukey"))
summary(PostHocTeste)
p9= ggboxplot(island_data, x = "occupancy", y = "Lat_centroid", ylab = "Latitude", xlab = "",
              color = "occupancy", palette =c("red3", "blue3", "green3")) +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") + 
  ggtitle("Lattitude") + theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9, p10, nrow = 3, ncol = 4, common.legend = T)


## Plot coefficients per realm (bubble plot)
## SR total
library(reshape)
coefs = coefficients(final_model_glmer_Richness_mammal)[[1]][,-1]
colnames(coefs) = c("Area", "SLMP", "Temperature mean", "Temperature sd", "Past_isolation",
                    "CCVT", "Precipitation sd", "Precipitation mean", "Elevation sd", "Area*SLMP")
coefs = cbind(coefs, rownames(coefs))
colnames(coefs)[ncol(coefs)] = "realm"
coefs <- melt(coefs, id.vars = c('realm'))
coefs$valueDiverging = abs(coefs$value)
coefs$valueDiscrete = ifelse(coefs$value > 0, "Positive", "Negative")
coefs$variable <- as.character(coefs$variable)
coefs$variable <- factor(coefs$variable, levels=c("Area", "SLMP", "Past_isolation", "CCVT", "Temperature mean", "Temperature sd", "Precipitation mean",
                                               "Precipitation sd", "Elevation sd", "Area*SLMP"))

# Discrete colour
pdf("bubble_plot_coefficients.pdf", onefile = T, width = 8, height = 8, useDingbats=FALSE)
ggplot(coefs, aes(variable, realm)) + ggtitle("SR total") +
  scale_y_discrete(limits = rev(levels(coefs$realm))) + xlab("") + ylab("") + scale_x_discrete(position = "top") +
  geom_exec(geom_point, data = coefs, size = "valueDiverging", fill = "valueDiscrete", shape = 21, color = "valueDiscrete") +
  scale_size(range = c(3, 13), breaks = c(0,0.1,0.2,0.3,0.5,0.75,1,1.25,1.5,1.75,2)) + theme_minimal()+
  scale_fill_manual(values=c("blue3", "red3")) + scale_color_manual(values=c("blue3", "red3"))+ 
  theme(axis.text.x = element_text(angle = 45, colour = "black", size = 12), axis.text.y = element_text(colour = "black", size = 12))


## SR bats
library(reshape)
coefs = coefficients(final_model_glmer_Richness_bat)[[1]][,-1]
colnames(coefs) = c("Area", "SLMP", "Temperature mean", "Temperature sd", "Past_isolation",
                    "CCVT", "Precipitation sd", "Precipitation mean", "Elevation sd", "Area*SLMP")
coefs = cbind(coefs, rownames(coefs))
colnames(coefs)[ncol(coefs)] = "realm"
coefs <- melt(coefs, id.vars = c('realm'))
coefs$valueDiverging = abs(coefs$value)
coefs$valueDiscrete = ifelse(coefs$value > 0, "Positive", "Negative")
coefs$variable <- as.character(coefs$variable)
coefs$variable <- factor(coefs$variable, levels=c("Area", "SLMP", "Past_isolation", "CCVT", "Temperature mean", "Temperature sd", "Precipitation mean",
                                                  "Precipitation sd", "Elevation sd", "Area*SLMP"))

# Discrete colour
ggplot(coefs, aes(variable, realm)) + ggtitle("SR bats") +
  scale_y_discrete(limits = rev(levels(coefs$realm))) + xlab("") + ylab("") + scale_x_discrete(position = "top") +
  geom_exec(geom_point, data = coefs, size = "valueDiverging", fill = "valueDiscrete", shape = 21, color = "valueDiscrete") +
  scale_size(range = c(3, 15), breaks = c(0,0.1,0.2,0.3,0.5,0.75,1,1.25,1.5,1.75,2,2.5,2.8)) + theme_minimal()+
  scale_fill_manual(values=c("blue3", "red3")) + scale_color_manual(values=c("blue3", "red3"))+ 
  theme(axis.text.x = element_text(angle = 45, colour = "black", size = 12), axis.text.y = element_text(colour = "black", size = 12))


## SR non-flying
library(reshape)
coefs = coefficients(final_model_glmmTMB_Richness_nonVol)[[1]][[1]][,-1]
colnames(coefs) = c("Area", "SLMP", "Temperature mean", "Temperature sd", "Past_isolation",
                    "CCVT", "Precipitation sd", "Precipitation mean", "Elevation sd", "Area*SLMP")
coefs = cbind(coefs, rownames(coefs))
colnames(coefs)[ncol(coefs)] = "realm"
coefs <- melt(coefs, id.vars = c('realm'))
coefs$valueDiverging = abs(coefs$value)
coefs$valueDiscrete = ifelse(coefs$value > 0, "Positive", "Negative")
coefs$variable <- as.character(coefs$variable)
coefs$variable <- factor(coefs$variable, levels=c("Area", "SLMP", "Past_isolation", "CCVT", "Temperature mean", "Temperature sd", "Precipitation mean",
                                                  "Precipitation sd", "Elevation sd", "Area*SLMP"))

# Discrete colour
ggplot(coefs, aes(variable, realm)) + ggtitle("SR non-flying") +
  scale_y_discrete(limits = rev(levels(coefs$realm))) + xlab("") + ylab("") + scale_x_discrete(position = "top") +
  geom_exec(geom_point, data = coefs, size = "valueDiverging", fill = "valueDiscrete", shape = 21, color = "valueDiscrete") +
  scale_size(range = c(3, 15), breaks = c(0,0.1,0.2,0.3,0.5,0.75,1,1.25,1.5,1.75,2,2.5,2.8)) + theme_minimal()+
  scale_fill_manual(values=c("blue3", "red3")) + scale_color_manual(values=c("blue3", "red3"))+ 
  theme(axis.text.x = element_text(angle = 45, colour = "black", size = 12), axis.text.y = element_text(colour = "black", size = 12))

## pSIE all
library(reshape)
coefs = coefficients(final_model_pSIE_total)[[1]][,-1]
colnames(coefs) = c("Area", "SLMP", "Temperature mean", "Temperature sd", "Past_isolation",
                    "CCVT", "Precipitation sd", "Precipitation mean", "Elevation sd")
coefs = cbind(coefs, rownames(coefs))
colnames(coefs)[ncol(coefs)] = "realm"
coefs <- melt(coefs, id.vars = c('realm'))
coefs$valueDiverging = abs(coefs$value)
coefs$valueDiscrete = ifelse(coefs$value > 0, "Positive", "Negative")
coefs$variable <- as.character(coefs$variable)
coefs$variable <- factor(coefs$variable, levels=c("Area", "SLMP", "Past_isolation", "CCVT", "Temperature mean", "Temperature sd", "Precipitation mean",
                                                  "Precipitation sd", "Elevation sd"))

# Discrete colour
ggplot(coefs, aes(variable, realm)) + ggtitle("pSIE all") +
  scale_y_discrete(limits = rev(levels(coefs$realm))) + xlab("") + ylab("") + scale_x_discrete(position = "top") +
  geom_exec(geom_point, data = coefs, size = "valueDiverging", fill = "valueDiscrete", shape = 21, color = "valueDiscrete") +
  scale_size(range = c(3, 11), breaks = c(0,0.1,0.2,0.3,0.5,0.75,1,1.25,1.45)) + theme_minimal()+
  scale_fill_manual(values=c("blue3", "red3")) + scale_color_manual(values=c("blue3", "red3"))+ 
  theme(axis.text.x = element_text(angle = 45, colour = "black", size = 12), axis.text.y = element_text(colour = "black", size = 12))

## pSIE bats
library(reshape)
coefs = coefficients(final_model_pSIE_bats)[[1]][,-1]
colnames(coefs) = c("Area", "SLMP", "Temperature mean", "Temperature sd", "Past_isolation",
                    "CCVT", "Precipitation sd", "Precipitation mean", "Elevation sd")
coefs = cbind(coefs, rownames(coefs))
colnames(coefs)[ncol(coefs)] = "realm"
coefs <- melt(coefs, id.vars = c('realm'))
coefs$valueDiverging = abs(coefs$value)
coefs$valueDiscrete = ifelse(coefs$value > 0, "Positive", "Negative")
coefs$variable <- as.character(coefs$variable)
coefs$variable <- factor(coefs$variable, levels=c("Area", "SLMP", "Past_isolation", "CCVT", "Temperature mean", "Temperature sd", "Precipitation mean",
                                                  "Precipitation sd", "Elevation sd"))

# Discrete colour
ggplot(coefs, aes(variable, realm)) + ggtitle("pSIE bats") +
  scale_y_discrete(limits = rev(levels(coefs$realm))) + xlab("") + ylab("") + scale_x_discrete(position = "top") +
  geom_exec(geom_point, data = coefs, size = "valueDiverging", fill = "valueDiscrete", shape = 21, color = "valueDiscrete") +
  scale_size(range = c(3, 12), breaks = c(0,0.1,0.2,0.3,0.5,0.75,1,1.25,1.5,1.65)) + theme_minimal()+
  scale_fill_manual(values=c("blue3", "red3")) + scale_color_manual(values=c("blue3", "red3"))+ 
  theme(axis.text.x = element_text(angle = 45, colour = "black", size = 12), axis.text.y = element_text(colour = "black", size = 12))

## pSIE non-flying
library(reshape)
coefs = coefficients(final_model_pSIE_nonVol)[[1]][,-1]
colnames(coefs) = c("Area", "SLMP", "Temperature mean", "Temperature sd", "Past_isolation",
                    "CCVT", "Precipitation sd", "Precipitation mean", "Elevation sd")
coefs = cbind(coefs, rownames(coefs))
colnames(coefs)[ncol(coefs)] = "realm"
coefs <- melt(coefs, id.vars = c('realm'))
coefs$valueDiverging = abs(coefs$value)
coefs$valueDiscrete = ifelse(coefs$value > 0, "Positive", "Negative")
coefs$variable <- as.character(coefs$variable)
coefs$variable <- factor(coefs$variable, levels=c("Area", "SLMP", "Past_isolation", "CCVT", "Temperature mean", "Temperature sd", "Precipitation mean",
                                                  "Precipitation sd", "Elevation sd"))

# Discrete colour
ggplot(coefs, aes(variable, realm)) + ggtitle("pSIE non-flying") +
  scale_y_discrete(limits = rev(levels(coefs$realm))) + xlab("") + ylab("") + scale_x_discrete(position = "top") +
  geom_exec(geom_point, data = coefs, size = "valueDiverging", fill = "valueDiscrete", shape = 21, color = "valueDiscrete") +
  scale_size(range = c(3, 13), breaks = c(0,0.1,0.2,0.3,0.5,0.75,1,1.25,1.5,1.75)) + theme_minimal()+
  scale_fill_manual(values=c("blue3", "red3")) + scale_color_manual(values=c("blue3", "red3"))+ 
  theme(axis.text.x = element_text(angle = 45, colour = "black", size = 12), axis.text.y = element_text(colour = "black", size = 12))
dev.off()


## Plot map with diversity (richness, endemism)
# Bubble chart
library(rworldmap)
continents <- readOGR("IUCN_basemap_continents.shp")
continents <- spTransform(continents, crs(PNAS_island_subset))
plot(continents)

library(rworldmap)
world <- getMap(resolution = "less islands")
world <- spTransform(world, crs(PNAS_island_subset))
world <- world[-which(world$REGION %in% "Antarctica"),]
world <- aggregate(world, "REGION")
plot(world)

PNAS_island_subset$Latitude <- round(coordinates(PNAS_island_subset)[,2],0)
PNAS_island_subset$Longitude <- round(coordinates(PNAS_island_subset)[,1],0)

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
library(classInt)
# SR total
classes = classIntervals(PNAS_island_subset@data$Richness_mammal, n = 10, style = "jenks")
breaks = range01(classes$brks)
pdf("Richness_mammal_plot_size_variation.pdf", onefile = T, width = 15, height = 7, useDingbats=FALSE)
ggplot(PNAS_island_subset@data, aes(Longitude, Latitude)) + theme_void() + 
  #geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "gray90", colour = "gray20", size = 0.005) +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "gray90", colour = "gray20", size = 0.005) +
  #geom_polygon(data = PNAS_island_subset, aes(x = long, y = lat, group = group), fill = "gray90", colour = "gray20", size = 0.005) +
  geom_exec(geom_point, data = PNAS_island_subset@data, size = "Richness_mammal", fill = "Richness_mammal", alpha = 0.5, shape = 21, color = "black") +
  scale_size(range = c(3, 13), breaks = classes$brks) +
  scale_fill_viridis_c(breaks = classes$brks, 
                       values = rev(breaks), guide = "legend")
dev.off()

# SR bats
classes = classIntervals(PNAS_island_subset@data$Richness_bat, n = 10, style = "jenks")
breaks = range01(classes$brks)
pdf("Richness_bat_plot_size_variation.pdf", onefile = T, width = 15, height = 7, useDingbats=FALSE)
ggplot(PNAS_island_subset@data, aes(Longitude, Latitude)) + theme_void() + 
  #geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "gray90", colour = "gray20", size = 0.005) +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "gray90", colour = "gray20", size = 0.005) +
  #geom_polygon(data = PNAS_island_subset, aes(x = long, y = lat, group = group), fill = "gray90", colour = "gray20", size = 0.005) +
  geom_exec(geom_point, data = PNAS_island_subset@data, size = "Richness_bat", fill = "Richness_bat", alpha = 0.5, shape = 21, color = "black") +
  scale_size(range = c(3, 13), breaks = classes$brks) +
  scale_fill_viridis_c(breaks = classes$brks, 
                       values = rev(breaks), guide = "legend")
dev.off()

# SR nonVol
classes = classIntervals(PNAS_island_subset@data$Richness_nonVol, n = 10, style = "jenks")
breaks = range01(classes$brks)
pdf("Richness_nonVol_plot_size_variation.pdf", onefile = T, width = 15, height = 7, useDingbats=FALSE)
ggplot(PNAS_island_subset@data, aes(Longitude, Latitude)) + theme_void() + 
  #geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "gray90", colour = "gray20", size = 0.005) +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "gray90", colour = "gray20", size = 0.005) +
  #geom_polygon(data = PNAS_island_subset, aes(x = long, y = lat, group = group), fill = "gray90", colour = "gray20", size = 0.005) +
  geom_exec(geom_point, data = PNAS_island_subset@data, size = "Richness_nonVol", fill = "Richness_nonVol", alpha = 0.5, shape = 21, color = "black") +
  scale_size(range = c(3, 13), breaks = classes$brks) +
  scale_fill_viridis_c(breaks = classes$brks, 
                       values = rev(breaks), guide = "legend")
dev.off()

# Endemism
island_data_endemism$Latitude <- round(island_data_endemism$Lat_centroid,0)
island_data_endemism$Longitude <- round(island_data_endemism$Long_centroid,0)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
library(classInt)

# SIE all
classes = classIntervals(island_data_endemism$SIE_total, n = 7, style = "jenks")
breaks = range01(classes$brks)
pdf("SIE_total_plot_size_variation.pdf", onefile = T, width = 15, height = 7, useDingbats=FALSE)
ggplot(PNAS_island_subset@data, aes(Longitude, Latitude)) + theme_void() + 
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "gray90", colour = "gray20", size = 0.005) +
  geom_polygon(data = PNAS_island_subset, aes(x = long, y = lat, group = group), fill = "gray90", colour = "gray20", size = 0.005) +
  geom_exec(geom_point, data = island_data_endemism, size = "SIE_total", fill = "SIE_total", alpha = 0.8, shape = 21, color = "black") +
  scale_size(range = c(3, 13), breaks = classes$brks) +
  scale_fill_viridis_c(breaks = classes$brks, 
                       values = rev(breaks), guide = "legend")
dev.off()

# SIE bats
classes = classIntervals(island_data_endemism$SIE_bats, n = 7, style = "jenks")
breaks = range01(classes$brks)
pdf("SIE_bats_plot_size_variation.pdf", onefile = T, width = 15, height = 7, useDingbats=FALSE)
ggplot(PNAS_island_subset@data, aes(Longitude, Latitude)) + theme_void() + 
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "gray90", colour = "gray20", size = 0.005) +
  geom_polygon(data = PNAS_island_subset, aes(x = long, y = lat, group = group), fill = "gray90", colour = "gray20", size = 0.005) +
  geom_exec(geom_point, data = island_data_endemism, size = "SIE_bats", fill = "SIE_bats", alpha = 0.8, shape = 21, color = "black") +
  scale_size(range = c(3, 13), breaks = classes$brks) +
  scale_fill_viridis_c(breaks = classes$brks, 
                       values = rev(breaks), guide = "legend")
dev.off()

# SIE nonVol
classes = classIntervals(island_data_endemism$SIE_nonVol, n = 7, style = "jenks")
breaks = range01(classes$brks)
pdf("SIE_nonVol_plot_size_variation.pdf", onefile = T, width = 15, height = 7, useDingbats=FALSE)
ggplot(PNAS_island_subset@data, aes(Longitude, Latitude)) + theme_void() + 
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "gray90", colour = "gray20", size = 0.005) +
  geom_polygon(data = PNAS_island_subset, aes(x = long, y = lat, group = group), fill = "gray90", colour = "gray20", size = 0.005) +
  geom_exec(geom_point, data = island_data_endemism, size = "SIE_nonVol", fill = "SIE_nonVol", alpha = 0.8, shape = 21, color = "black") +
  scale_size(range = c(3, 13), breaks = classes$brks) +
  scale_fill_viridis_c(breaks = classes$brks, 
                       values = rev(breaks), guide = "legend")
dev.off()

{
# Equal sizes
pdf("Richness_mammal_plot_one_size.pdf", onefile = T, width = 15, height = 7, useDingbats=FALSE)
ggplot(PNAS_island_subset@data, aes(Longitude, Latitude)) + theme_void() +
  #geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "gray90", colour = "gray20", size = 0.005) +
  #geom_polygon(data = PNAS_island_subset, aes(x = long, y = lat, group = group), fill = "gray90", colour = "gray20", size = 0.005) +
  geom_exec(geom_point, data = PNAS_island_subset@data, size = 4, fill = "Richness_mammal", alpha = 0.8, shape = 21, color = "black") +
  scale_fill_viridis_c(breaks = classes$brks, values = rev(breaks), guide = "legend")
dev.off()

# Different sizes ring
pdf("Richness_mammal_plot_size_variation_ring2.pdf", onefile = T, width = 15, height = 7, useDingbats=FALSE)
ggplot(PNAS_island_subset@data, aes(Longitude, Latitude)) + theme_void() + 
  #geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "gray90", colour = "gray20", size = 0.005) +
  #geom_polygon(data = PNAS_island_subset, aes(x = long, y = lat, group = group), fill = "gray90", colour = "gray20", size = 0.005) +
  geom_point(data = PNAS_island_subset@data, aes(Longitude, Latitude, color = Richness_mammal, size = Richness_mammal), 
             alpha = 1, shape = 21, stroke = 1.6) +
  scale_size(range = c(2, 6), breaks = classes$brks) +
  scale_color_viridis_c(option = "D", breaks = classes$brks, values = rev(breaks), guide = "legend")
dev.off()

# Equal sizes ring
pdf("Richness_mammal_plot_one_size_ring.pdf", onefile = T, width = 15, height = 7, useDingbats=FALSE)
ggplot(PNAS_island_subset@data, aes(Longitude, Latitude)) + theme_void() + 
  #geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "gray90", colour = "gray20", size = 0.005) +
  #geom_polygon(data = PNAS_island_subset, aes(x = long, y = lat, group = group), fill = "gray90", colour = "gray20", size = 0.005) +
  geom_point(data = PNAS_island_subset@data, aes(Longitude, Latitude, color = Richness_mammal), 
             alpha = 0.6, shape = 21, stroke = 1.7, size = 4) +
  scale_color_viridis_c(breaks = classes$brks, 
                        values = rev(breaks), guide = "legend")
dev.off()
}

## Histogram of species richness per bioregion
# Species richness
bioreg <- levels(PNAS_island_subset$bioregion)
pdf("Histogram_SR.pdf", onefile = T, width = 9, height = 8, useDingbats=FALSE)
par(mfrow = c(3,3))
for(i in 1:length(bioreg)){
  temp <- subset(PNAS_island_subset, PNAS_island_subset$bioregion %in% bioreg[i])
  hist(log(temp@data$Richness_mammal), xlab = "log Species richness", ylab = "Frequency", main = bioreg[i])
  hist(log(temp@data$Richness_bat), xlab = "log Species richness", ylab = "Frequency", main = bioreg[i])
  hist(log(temp@data$Richness_nonVol), xlab = "log Species richness", ylab = "Frequency", main = bioreg[i])
}
dev.off()


## Histogram of variables after and before transformation
pdf("Variables_histogram.pdf", onefile = T, width = 9, height = 8, useDingbats=FALSE)
par(mfrow = c(3,3))
variables <- c("Area", "Area", "Current_isolation", "Temperature_mean", "Temperature_sd", "Past_isolation", 
               "CCVT", "Climate_velocity", "Prec_sd_Chelsa", "Precipitation_sd", "Prec_mean_Chelsa", "Precipitation_mean", 
               "Elev_sd_GTOPO", "Elevation_sd")
for(i in 1:length(variables)){
  hist(island_data[,variables[i]], xlab = variables[i], ylab = "Frequency", main = variables[i])
}
dev.off()
