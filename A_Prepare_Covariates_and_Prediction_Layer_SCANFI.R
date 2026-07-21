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
require(sp)
require(data.table)
require(ggplot2)
require(gdata)
require(tidyterra)
require(exactextractr)
require(tidyr)
require(R.utils)

`%unlike%`<-Negate(`%like%`)

terra::tmpFiles(current=T,orphan=T,old=T,remove=T)

# Set up path for data on google drive and output folder ----
root<-"G:/XXXXXXXXXXX" #PC link

#Get shortcuts to SCANFI 2 files ----
SCANFIage <- filePath(file.path(root, "SppHabitatAssoc_byBCR", "SCANFI2","SCANFI_age_v2_20260119.lnk"), expandLinks="any")
SCANFIatt <- filePath(file.path(root, "SppHabitatAssoc_byBCR", "SCANFI2","SCANFI_attributes_v2_20260119.lnk"), expandLinks="any")

# Select out folder ---
out<-"C:/Users/andrake/OneDrive - NRCan RNCan/Desktop/HabitatAssociation_SCANFIage"
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
  st_transform(crs=crs("EPSG:5072")) # Canadian boundary

# BCRs ----
bcr<- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")) %>%
  st_transform(crs=crs("EPSG:5072"))|> st_crop(can)%>% st_intersection(can) 

# Remove nothern portion ----
bcrS<-bcr[-c(1,11),]

###############################################
# 1. PREDICTION LAYERS 
###############################################

# Set up a template for prediction projection
Template<-rast(file.path(root,"PredictionRasters/ClimateNormal/FFP_1km.tif"))%>% crop(bcrS)%>%mask(bcrS)

# Note on aggregate method
# 30m to 200m radius area is 3% over
# 30m to 2000m radius area is 0.2% over

## Import tree and age data -----
# Conifer layers ---
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

# Deciduous layer ---
Dlist<-files.at2$path %>% 
  .[.%like% "2025"]%>%
  .[.%unlike% ".tif."]%>%
  .[.%like% "broadleaf"]

prcD<-rast(file.path(Dlist)) %>% aggregate(.,12,"mean", na.rm=T)%>% round(., digits=0)

## Calc Conifer/Deciduous dominance @ 200m and reproject ----
Conifer<-prcCs/(prcCs+prcD) %>% round(2)
Con200m<-Conifer%>% project(., Template, method="near") %>% crop(bcrS)%>%mask(bcrS) 

## Percent forested ----
pForest<-prcCs+prcD
For200m<-round(pForest,0)%>% project(., Template, method="near") %>% crop(bcrS)%>%mask(bcrS) 

### Age files --------------
Alist<-files.ag3$path %>% 
  .[.%like% "2025"]%>%
  .[.%unlike% ".tif."]
ForestAge<-rast(file.path(Alist))

#Now aggregate age at same dimensions using median ----
AgeAgg<-aggregate(ForestAge,12,"median",ties="random", na.rm=T) 
Age200m<-AgeAgg%>%project(., Template, method="near") %>% crop(bcrS)%>%mask(bcrS) 

### Save aggregated layers ----
writeRaster(Conifer,"ConAggregate.tif", overwrite=TRUE)
writeRaster(pForest,"ForAggregate.tif", overwrite=TRUE)

#Get landscape values  -------------
# Mean of means works for these so aggregate to speed up ----

For2k1<-pForest%>% aggregate(.,2,"mean", na.rm=T) %>% focal(., w=5, fun="mean", na.rm=T) # get a focal average (could use focal_mat but it removes a lot of area at 720 res)
For2k<-For2k1%>% terra::project("EPSG:3348") %>% terra::project(., Template, method="near") %>% # reproject to 1 km through 3348 to prevent errors
  crop(bcrS)%>%mask(bcrS) # 2km mean forest %

Con2k1<-Conifer%>% aggregate(.,2,"mean", na.rm=T)%>% focal(., w=5, fun="mean", na.rm=T)
Con2k<-Con2k1%>% terra::project("EPSG:3348") %>% terra::project(., Template, method="near") %>% 
  crop(bcrS)%>%mask(bcrS)

# Age: Can't use median of median so have to work from 30m (slow) --------

# Create 1km raster template in the same projection as ForestAge
Template2<-Template%>% terra::project("EPSG:3348")%>% terra::project(crs(ForestAge)) # for projection to work need to run through 3348
points_df <- as.data.frame(Template2, xy = TRUE)
points_sf <- st_as_sf(points_df, coords = c("x", "y"), crs = crs(Template2))
reg.buff <- st_buffer(points_sf, 2000)  # 8 million points at 1km across Canada
rm(points_sf)
gc()

# Use exact extract to get values for each point and then convert to raster (slow) ----
Age2kcov <-exact_extract(x = Age, y = reg.buff, 'median', progress = TRUE)

# back up ----
write.csv(Age2kcov,"2kAge.csv")

#transform back to raster ----------
points_df<-cbind(points_df,LandAge=round(Age2kcov,0))
Age2kr<-terra::rast(points_df, crs=crs(Template2))
Age2kFinal<-Age2kr%>% terra::project("EPSG:3348") %>% terra::project(Template)
Age2kFinal<-Age2kFinal%>% crop(bcrS)%>%mask(bcrS)

####################################################
# Write out prediction rasters and prediction stack
####################################################
# Local 200m
writeRaster(Age200m,"Loc_Age1km.tif", overwrite=TRUE)
writeRaster(Con200m,"Loc_Con1km.tif", overwrite=TRUE)
writeRaster(For200m,"Loc_For1km.tif", overwrite=TRUE)

# Landscape 2000m
writeRaster(For2k,"Land_For1km.tif", overwrite=TRUE)
writeRaster(Age2kFinal,"Land_Age1km.tif", overwrite=TRUE)
writeRaster(Con2k,"Land_Con1km.tif", overwrite=TRUE)

## Set predict year set to 2025 (to match landcover) and method to "PC" ----
method<-YR<-Con200m
values(YR)[values(!is.na(YR))]<-2025
values(method)[values(!is.na(method))]<-0

# Stack all --------
pred<-c(Con200m,Con2k,For200m,For2k,Age200m,Age2k, YR, method)
names(pred)<-c("Loc_pConifer","Land_pConifer","Loc_closure","Land_closure","LocAge","LandAge","year","method")

# Remove lake areas -------------
Lakes<-read_sf(file.path("C:/Users/andrake/OneDrive - NRCan RNCan/Desktop/Partitioning Change","gis","lhy_000c16a_e.shp"))%>%st_transform(crs=crs(pred))
Lakes$area<-sf::st_area(Lakes$geometry)
Lakes<-Lakes%>%filter(area>units::set_units(2500000, m^2))  #the mean of the shapefiles
Lakes<-rasterize(Lakes,pred)

pred<-pred%>%mask(Lakes, inverse=T)
writeRaster(pred, filename="predictionlayer2025.tif", overwrite=TRUE)

###########################################################
# 2. EXTRACT MODEL COVARIATES
###########################################################

#6. Load data package ----
load(file.path(root, "data/04_NM5.0_data_stratify.Rdata"))
rm('bird','offsets','cov','bootlist','bcrlist') # remove heavy items
gc()

# Import SCANFI template
Template<-rast(files.ag3[1,1])
bcrS2<-bcrS%>%st_transform(crs=crs(Template))

# Get Canadian points and buffer to 200m ----
loc.n <- visit |> 
  dplyr::select(id, project, location, lat, lon, year) |> 
  unique() |> 
  st_as_sf(coords=c("lon", "lat"), crs=4326, remove=FALSE) |> 
  st_transform(crs=crs(Template))|> st_crop(bcrS2) 

# Remove ID and get unique location x year values (cuts total by ~60,000) 
loc.u<-visit |> 
  dplyr::select(project, location, lat, lon, year) |> 
  unique() |> 
  st_as_sf(coords=c("lon", "lat"), crs=4326, remove=FALSE) |> 
  st_transform(crs=crs(Template))|> st_crop(bcrS2) 
remove("visit")
gc()

####################################################
# Extract time-matched values at each survey point
####################################################

#Split surveys into years
bin<-c(0,1985,1990,1995,2000,2005,2010,2015,2020,2025,2030)

for (i in c(1:9)){ # run through each 5-year period

# Median  forest age in 200m buffer --------
Alist<-files.ag3$path %>% 
    .[.%like% bin[i+1]]%>% #Get rasters for period i+1 
    .[.%unlike% ".tif."]

Age <- rast(Alist)

# Extract all the tree data and process afterwards
Tlist<-files.at2$path %>% 
    .[.%like%  bin[i+1]]%>% #Get rasters for period i+1 
    .[.%unlike% ".tif."]%>%
    .[.%unlike% "biomass"]%>%
    .[.%unlike% "height"]%>%
    .[. %unlike% "nfiLandcover"]
  
Trees<-rast(Tlist)
print(bin[i+1])

# extract surveys that occurred in this period --------------
loc.y<-loc.u|>filter(year>bin[i]&year<=bin[i+1])

# still too large so split into 20 chunks for processing -----------------
count<-round(nrow(loc.y)/20,0)  
  
for (j in c(0:19)){
start<-(j*count)+1
end<-(j+1)*count
  
if(end>nrow(loc.y)){
  end<-nrow(loc.y)
  }

loc.y.j<-loc.y[c(start:end),] # chunk j

# Buffer for extraction ----
loc.buff <- st_buffer(loc.y.j, 200)
reg.buff <- st_buffer(loc.y.j, 2000)
print("done buffering")

Age_values <- exact_extract(x = Age, y = loc.buff, 'median', progress = FALSE)

Trees200mcov <- exact_extract(x=Trees,y=loc.buff, 'mean', progress= FALSE)
colnames(Trees200mcov) <- gsub("mean.SCANFI_spsCC_", "Loc_", colnames(Trees200mcov))
print("done local")
rm("loc.buff")

# Extract modal forest age in 2000m buffer (v. slow) ----
Age_2k <- exact_extract(x = Age, y = reg.buff, 'median', progress = FALSE)

Trees2kcov <- exact_extract(x=Trees,y=reg.buff, 'mean', progress= FALSE)
colnames(Trees2kcov) <- gsub("mean.SCANFI_spsCC_", "Land_", colnames(Trees2kcov))

AllOut <- cbind(
  sf::st_drop_geometry(loc.y.j), 
  LocAge = round(Age_values, 0), 
  Trees200mcov,
  LandAge = round(Age_2k, 0), 
  Trees2kcov
)

write.csv(AllOut,paste0("Classification_covariates_",i,"_",j,".csv"))
rm("reg.buff","AllOut")
gc()
  
} # end of chunk
} # end of period

#### Done getting buffered values ######################

###################################################
## Roll values into single dataset 
###################################################
bin<-c(0,1985,1990,1995,2000,2005,2010,2015,2020,2025,2030)

outputsFinal<-list()
for (i in c(1:9)){
  
yr<-bin[i+1] # year of age raster used
PointVal <- list.files(file.path(out, "RawCovariate"), recursive = TRUE, full.names = TRUE)
PointVal.i<-PointVal%>%.[.%like% paste0("_",i,"_")]
outputs<-lapply(PointVal.i, read.csv)
outJ<-do.call(rbind,outputs)

#Simplify and unify column names
#Land and local
outJ <- outJ %>%
  rename_with(~ gsub("_\\d{4}.*$", "", .), starts_with("Land_"))
outJ <- outJ %>%
  rename_with(~ gsub("_\\d{4}.*$", "", .), starts_with("Loc_"))
# Mean closure
outJ <- outJ %>%
  rename_with(~ c("Loc_closure","Land_closure"), starts_with("mean."))

### Correct Age class to survey year (within the 5-year window of the layer) -----
  
outJ$LocAge<-outJ$LocAge-(yr-outJ$year) # correct age by survey year
outJ$LocAge<-ifelse(outJ$LocAge<0,NA,outJ$LocAge) # remove any stands that were disturbed in that window (negative age)

outJ$LandAge<-outJ$LandAge-(yr-outJ$year) # correct age by survey year
outJ$LandAge<-ifelse(outJ$LandAge<0,NA,outJ$LandAge) # remove any stands that were disturbed in that window 
outJ<-outJ[-1]
outputsFinal[[i]]<-outJ

}
outFinal<-do.call(rbind,outputsFinal)
outFinal<-unique(outFinal) # some duplicates - skip these
locationID <- st_drop_geometry(loc.n)
Dataset<-merge(locationID,outFinal, by=c("project", "location", "lat","lon", "year"))
write.csv(Dataset,"SurveyCovariates.csv") # back-up
# confirm: should be 1,329,023 records - 9 missing

# Limit to final covariates ---
# Calc Conifer/Deciduous dominance----

Dataset$Loc_pConifer<-(Dataset$Loc_balsamFir+
                      Dataset$Loc_blackSpruce + 
                      Dataset$Loc_douglasFir+ 
                      Dataset$Loc_jackPine+ 
                      Dataset$Loc_lodgepolePine + 
                      Dataset$Loc_otherConiferous +
                      Dataset$Loc_ponderosaPine + 
                      Dataset$Loc_tamarack+        
                      Dataset$Loc_whiteRedPine)/
    (Dataset$Loc_balsamFir+ 
     Dataset$Loc_blackSpruce + 
     Dataset$Loc_broadleaf + # broadleaf = deciduous
     Dataset$Loc_douglasFir+ 
     Dataset$Loc_jackPine+ 
     Dataset$Loc_lodgepolePine + 
     Dataset$Loc_otherConiferous +
     Dataset$Loc_ponderosaPine + 
     Dataset$Loc_tamarack+        
     Dataset$Loc_whiteRedPine)
Dataset$Loc_pConifer<-round(Dataset$Loc_pConifer,digits=2)

# landscape values

Dataset$Land_pConifer<-(Dataset$Land_balsamFir+
                         Dataset$Land_blackSpruce + 
                         Dataset$Land_douglasFir+ 
                         Dataset$Land_jackPine+ 
                         Dataset$Land_lodgepolePine + 
                         Dataset$Land_otherConiferous +
                         Dataset$Land_ponderosaPine + 
                         Dataset$Land_tamarack+        
                         Dataset$Land_whiteRedPine)/
  (Dataset$Land_balsamFir+ 
     Dataset$Land_blackSpruce + 
     Dataset$Land_broadleaf + # broadleaf = deciduous
     Dataset$Land_douglasFir+ 
     Dataset$Land_jackPine+ 
     Dataset$Land_lodgepolePine + 
     Dataset$Land_otherConiferous +
     Dataset$Land_ponderosaPine + 
     Dataset$Land_tamarack+        
     Dataset$Land_whiteRedPine)
Dataset$Land_pConifer<-round(Dataset$Land_pConifer,digits=2)

DatasetS<-Dataset[c("project","location","lat","lon","year",
                    "id","LocAge","Loc_closure","LandAge",             
                    "Land_closure","Loc_pConifer","Land_pConifer")]  
DatasetS$Loc_closure<-round(DatasetS$Loc_closure,2)
DatasetS$Land_closure<-round(DatasetS$Land_closure,2)
write.csv(DatasetS,"ClassificationCovariates.csv")

############ end of dataset ##############
