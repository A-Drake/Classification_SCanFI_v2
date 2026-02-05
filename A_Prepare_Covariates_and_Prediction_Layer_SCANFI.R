###################################################
### Species Classification 
### A. Prepare Covariate and Prediction Raster layers - SCANFI V2
### by Anna Drake
### February 3 2026 ####
##################################################

require(terra)
require(dplyr)
require(plyr)
require(sf)
require(data.table)
require(ggplot2)
require(gdata)
require(tidyterra)
require(exactextractr)
require(tidyr)
require(R.utils)

`%unlike%`<-Negate(`%like%`)

# Set up path for data on google drive and output folder ----
root<-"G:/Shared drives/BAM_NationalModels5" #PC link

#Get shortcuts to SCANFI 2 files ----
SCANFIage <- filePath(file.path(root, "SppHabitatAssoc_byBCR", "SCANFI2","SCANFI_age_v2_20260119.lnk"), expandLinks="any")
SCANFIatt <- filePath(file.path(root, "SppHabitatAssoc_byBCR", "SCANFI2","SCANFI_attributes_v2_20260119.lnk"), expandLinks="any")
out<-"C:/Users/andrake/OneDrive - NRCan RNCan/Desktop/HabitatAssoications_SCANFIage"
setwd(out)

# Get file lists for attributes and age  ------------
files.at <- data.frame(path = list.files(file.path(SCANFIatt), recursive = TRUE, full.names = TRUE))
files.ag <- data.frame(path = list.files(file.path(SCANFIage), recursive = TRUE, full.names = TRUE))

files.at2 <- files.at |>
  separate(path, into=c("drive","f1", "f2","f3","file"), remove=FALSE, sep="/") |>
  separate(file, into=c("derive","type","value","year","version","date","filetype1", "filetype2","filetype3")) |>
  mutate(year = as.numeric(year))

files.ag3 <- files.ag |>
  separate(path, into=c("drive", "f1", "f2","f3","file"), remove=FALSE, sep="/") |>
  separate(file, into=c("derive","type","value","year","version","date","filetype1", "filetype2","filetype3")) |>
  mutate(year = as.numeric(year))

# Import shapefiles ----
# Canada ----
can <- read_sf(file.path(root, "Regions", "CAN_adm", "CAN_adm0.shp")) %>% 
  st_transform(crs=crs("EPSG:3348")) # Canadian boundary

# BCRs ---- # get 61
bcr<- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")) %>%
  st_transform(crs=crs("EPSG:3348"))|> st_crop(can)%>% st_intersection(can) #%>%
# Remove nothern portion ----
bcrS<-bcr[-c(1,11),]
plot(bcrS$geometry)

########################################
# Import species
########################################
# Species list ---
bootfolder<-"G:/Shared drives/BAM_NationalModels5/output/07_predictions/"
specieslist<-list.files(bootfolder)

wetland<-c("ALFL","COYE","WISN","SPSA","NESP","MAWR","YHBL","RWBL","SEWR","SWSP","BEKI","SOSA","GRYE","KILL","LEYE","UPSA","RUBL","DUNL")
grassland<-c("BOBO","GRSP","HOLA","LCSP","SAVS","UPSA","VESP")
irruptive<-c("WWCR","RECR","PISI","EVGR","PIGR")
anthr<-c("HOSP","ROPI","EUST")# keep anthropological expansion spp
open<-c("WIPT","ROPT","WIPT","RKPT","RLHA","GYRF","SNOW","BLUE","NOWH","EYWG","HASP","LALO","SMLO","SNBU","PEFA","AMPI")  # adding AMPI keeping CHCH, HORE and CORE

# retain these where not in the above categories ----
nonforest<-c("AMCR","AMGO","AMPI","ATSP","BAOR","BANS","BARS","BEKI","BBMA","BLJA",
             "BGGN","BOBO","BRBL","BRTH","BHCO","CHSP","CCSP","CLSW","COGR","COYE","EABL","EAKI","EAPH","EATO","EUST","FISP","FOSP",
             "GCSP","GWWA","GRSP","GRCA","GRYE","HOLA","HOSP","HOWR","INBU","KILL","LALO","LCSP","LEYE","LISP","MAWR","MOBL","MODO","NESP",
             "NOCA","PAWA","RWBL","ROPI","RUBL","SAVS","SEWR","SOSA","SOSP","SPSA","SWSP","TRES","UPSA","VESP","WCSP","WIPT","WISN",
             "WIWA","YHBL","ALFL","AMRO","BBCU","YEWA","CORA","NOFL","CEDW")

marginal<-c("BGGN","BWWA","GWWA","RBWO","RHWO","ROSA","WITU","YBCU","YTVI") # species barely in Canada

exclude<-c(wetland,grassland,open,marginal)%>%unique()

specieslist<-specieslist[specieslist %notin% exclude]


###############################################
# 1. PREDICTION LAYERS 
###############################################

# Landscape scale extraction using circle matrix on coarsened 30m layers

# Number of cells to cover radius
n_cells <- 6 #(2000/360); 360=12*30
radius<-2000

# Create a matrix of distances from the center
focal_mat <- matrix(0, nrow=2*n_cells + 1, ncol=2*n_cells + 1)
center <- n_cells + 1

for (i in 1:nrow(focal_mat)) {
  for (j in 1:ncol(focal_mat)) {
    # Euclidean distance from center in meters
    dist <- sqrt((i - center)^2 + (j - center)^2) * 360
    if (dist <= radius) {
      focal_mat[i, j] <- 1
    }
  }
}

## Import tree and age data -----

Template<-rast(file.path(root,"PredictionRasters/Wetland/1km/SurfaceWater_1km.tif")) %>%
  project("EPSG:3348", method="near") %>% crop(bcrS)%>%mask(bcrS)

# get conifer layers ---
Conlist<-files.at2$path %>% 
  .[.%like% "2025"]%>%
  .[.%unlike% ".tif."]%>%
  .[.%unlike% "prcD"]%>%
  .[.%unlike% "biomass"]%>%
    .[.%unlike% "closure"]%>%
      .[.%unlike% "height"]%>%
  .[.%unlike% "broadleaf"]%>%
  .[. %unlike% "nfiLandcover"]

prcC<-rast(file.path(Conlist))

#sum and aggregate values ~200m  
prcC<-prcC%>% aggregate(.,12,"mean", na.rm=T)%>% round(., digits=0)
prcCs<-app(prcC,fun="sum",na.rm=T)

#Get deciduous layer ---

Dlist<-files.at2$path %>% 
  .[.%like% "2025"]%>%
  .[.%unlike% ".tif."]%>%
  .[.%like% "broadleaf"]

prcD<-rast(file.path(Dlist)) %>% aggregate(.,12,"mean", na.rm=T)%>% round(., digits=0)

## Calc Conifer/Deciduous dominance @ 200m and reproject ----
Conifer<-prcCs/(prcCs+prcDs) %>% round(2)
plot(Conifer)
Con200m<-Conifer%>% project(., Template, method="near") %>% crop(bcrS)%>%mask(bcrS) 

## Percent forested ----
pForest<-prcC+prcD
For200m<-round(pForest,0)%>% project(., Template, method="near") %>% crop(bcrS)%>%mask(bcrS) 

## Set predict year set to 2025 (to match landcover) and method to PC ----
method<-YR<-Con200m
values(YR)[values(!is.na(YR))]<-2025
values(method)[values(!is.na(method))]<-0

### Age files --------------
Alist<-files.ag3$path %>% 
  .[.%like% "2025"]%>%
  .[.%unlike% ".tif."]
ForestAge<-rast(file.path(Alist))

#Now aggregate age at same dimensions ----
AgeAgg<-aggregate(ForestAge,12,"median",ties="random", na.rm=T)  # dominant age class at ~210 m radius
#AgeAge<-mask*AgeAgg # keep if at least 14 estimates are within radius
Age200m<-AgeAgg%>%project(., Template, method="near") %>% crop(bcrS)%>%mask(bcrS) 

# Get landscape values  -------------
Con2k<-focal(Conifer,w=focal_mat, "mean", na.rm=T)%>% 
  crop(bcrS)%>%mask(bcrS) # 2km mean conifer mix
For2k<-focal(pForest,w=focal_mat, "mean", na.rm=T)%>% 
  crop(bcrS)%>%mask(bcrS)  # 2km mean coverage
Age2k<-focal(AgeAgg, w=focal_mat, "median", ties="random", na.rm=T)%>% 
  crop(bcrS)%>%mask(bcrS)  # 2km mean coverage

####################################################
# Write out prediction rasters and prediction stack
####################################################

writeRaster(Age200m,"Loc_Age1km.tif", overwrite=TRUE)
writeRaster(Age2k,"Land_Age1km.tif", overwrite=TRUE)
writeRaster(Con200m,"Loc_Con1km.tif", overwrite=TRUE)
writeRaster(Con2k,"Land_Con1km.tif", overwrite=TRUE)
writeRaster(For200m,"Loc_For1km.tif", overwrite=TRUE)
writeRaster(For2k,"Land_For1km.tif", overwrite=TRUE)

pred<-c(Con200m,Con2k,For200m,For2k,Age200m,Age2k, YR, method)
names(pred)<-c("LocCon","LandCon","LocFor","LandFor","LocAge","LandAge","year","method")

# Remove lake areas -------------
Lakes<-read_sf(file.path("C:/Users/andrake/OneDrive - NRCan RNCan/Desktop/Partitioning Change","gis","lhy_000c16a_e.shp"))%>%st_transform(crs=crs(pred))
Lakes$area<-sf::st_area(Lakes$geometry)
summary(Lakes$area)
Lakes<-Lakes%>%filter(area>units::set_units(27840000, m^2))
Lakes<-rasterize(Lakes,pred)
pred<-pred%>%mask(Lakes, inverse=T)
writeRaster(pred, filename="predictionlayer.tif", overwrite=TRUE)

###########################################################
# 2. MODEL COVARIATES
###########################################################

#6. Load data package ----
load(file.path(root, "data/04_NM5.0_data_stratify.Rdata"))
rm('bird','offsets') # remove heavy items
gc()

# Get Canadian points and buffer to 200m -----  # confirm projection lines up
loc.n <- visit |> 
  dplyr::select(id, project, location, lat, lon, year) |> 
  unique() |> 
  st_as_sf(coords=c("lon", "lat"), crs=4326, remove=FALSE) |> 
  st_transform(crs=crs("EPSG:3348"))|> st_crop(bcr) |>
  st_transform(crs=crs(ForestAge))

remove("visit")
gc()

#### Buffer for extraction ----
loc.buff <- st_buffer(loc.n, 200)
reg.buff <- st_buffer(loc.n, 2000)

# Extract modal forest age in 200m buffer --------
Age<-rast(file.path(Alist))

Age200mcov <- loc.buff |> 
  exact_extract(x=Age, 'median', force_df=TRUE) |> 
  data.table::setnames("LocAge")|> 
  cbind(loc.n)

gc()

# Extract all the tree data and process afterwards
Trees<-files.at2$path %>% 
  .[.%like% "2025"]%>%
  .[.%unlike% ".tif."]%>%
  .[.%unlike% "prcD"]%>%
  .[.%unlike% "biomass"]%>%
  .[.%unlike% "closure"]%>%
  .[.%unlike% "height"]%>%
  .[. %unlike% "nfiLandcover"]

Trees200mcov <- loc.buff |> 
  exact_extract(x=Trees, 'mean', force_df=TRUE) |> 
  cbind(loc.n)

rm("loc.buff")

# Extract modal forest age in 2000m buffer - v. slow ----
AgeReg <- reg.buff |> 
  exact_extract(x=Age, 'median', force_df=TRUE) |> 
  data.table::setnames("LandAge")|> 
  cbind(loc.n)

Trees2kcov <- reg.buff |> 
  exact_extract(x=Trees, 'mean', force_df=TRUE) |> 
  cbind(loc.n)

rm("reg.buff")

gc()

#drop geometry and round to whole years ----
AgeReg<-AgeReg[-8]


cov2<-cov[,c(1,18,22,26,28,30,36,38,40,42,44,103)]


# There are 10 forest metrics - 9 conifer and 1 deciduous ---

# Conifer ----
cov2$C<-round((cov2$SCANFIBalsamFir_1km+
                      cov2$SCANFIBlackSpruce_1km+
                      cov2$SCANFIprcC_1km+
                      cov2$SCANFIDouglasFir_1km +
                      cov2$SCANFIJackPine_1km+
                      cov2$SCANFILodgepolePine_1km+
                      cov2$SCANFIPonderosaPine_1km +
                       cov2$SCANFITamarack_1km +
                       cov2$SCANFIWhiteRedPine_1km),0)
#total forest ----
cov2$Forest<-round((cov2$SCANFIBalsamFir_1km+
                      cov2$SCANFIBlackSpruce_1km+
                      cov2$SCANFIprcC_1km+
                      cov2$SCANFIDouglasFir_1km +
                      cov2$SCANFIJackPine_1km+
                      cov2$SCANFILodgepolePine_1km+
                      cov2$SCANFIPonderosaPine_1km +
                      cov2$SCANFITamarack_1km +
                      cov2$SCANFIWhiteRedPine_1km +
                      cov2$SCANFIprcD_1km),0)

cov2$Forest<-ifelse(cov2$Forest>100,100,cov2$Forest)

# Proportion conifer -----
cov2$Conifer<-round((cov2$C/cov2$Forest),2)

names(AgeReg)
names(cov2)
# Put together into a single covariate set ----
Covar<-merge(cov2[,c(1,12:13,15:16)],Agecov[,c(2,8)], by="id", all.x=T)
Covar<-merge(Covar,AgeReg[,c(2,8)], by="id",all.x=T)  ###TO HERE
View(Covar)
write.csv(Covar,"ClassificationCovariates.csv")

######### End of input prep, survey data comes from V5 objects  #############