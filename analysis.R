#
#
#

# load libraries
library(tidyverse); library(sp); library(sf); library(raster); library(data.table); library(mapdata)
library(maptools); library(gridExtra); library(nngeo); library(stringr); library(rgdal); library(scales)

# basemap
crs.1 <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
land <- st_read("ne_50m_admin_0_countries.shp")
land <- st_transform(land, st_crs("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"))

rcp8.5_2080 <- raster("covariates/fwvel_ensemble_rcp85_2085.tif")
rcp8.5_2080 <- projectRaster(rcp8.5_2080, 
                             crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m",
                             res=10000, method="bilinear")

# read in base grid data
grd_100 <- st_read("100km_gridClean.shp") %>% dplyr::select(X,Y)
grd_100 <- unique(grd_100) %>% mutate(FID=row_number())
st_crs(grd_100) <- crs.1

grd_200 <- st_read("200km_gridClean.shp") %>% dplyr::select(X,Y)
grd_200 <- unique(grd_200) %>% mutate(FID=row_number())
st_crs(grd_100) <- crs.1

grd_400 <- st_read("400km_gridClean.shp") %>% dplyr::select(X,Y)
grd_400 <- unique(grd_400) %>% mutate(FID=row_number())
st_crs(grd_100) <- crs.1

# grd_800 <- st_read("800km_gridClean.shp") %>% dplyr::select(X,Y)
# grd_800 <- unique(grd_800) %>% mutate(FID=row_number())
# st_crs(grd_100) <- 102008

# human footprint map
footprint <- raster("covariates/footprint_clipped.tif")
footprint <- projectRaster(footprint, 
                           crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m",
                           res=10000, method="bilinear")
footprint <- crop(footprint, grd_100)
footprint[footprint > 50] <- NA

# WWF biomes
biomes <- raster("covariates/WWFBiomes.tif")
crs(biomes) <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
biomes <- crop(biomes, grd_100)
biomes <- projectRaster(biomes,
                        crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m",
                           res=10000, method="bilinear")
biomes[biomes > 50] <- NA

# read in fishnet data and count overall richness by overlap
# 100 km resolution
fsh_100_Buff <- fread("100km_Buffer.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
fsh_100_Buff <- st_as_sf(fsh_100_Buff, coords=c("X", "Y"))
st_crs(fsh_100_Buff) <- crs.1

grd_100$Buff100 <- lengths(st_intersects(grd_100, fsh_100_Buff))

# 200 km resolution 
fsh_200_noBuff <- fread("200km_NoBuffer.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
fsh_200_noBuff <- st_as_sf(fsh_200_noBuff, coords=c("X", "Y"))
st_crs(fsh_200_noBuff) <- crs.1

grd_200$noBuff200 <- lengths(st_intersects(grd_200, fsh_200_noBuff))


# 400 km resolution
fsh_400_noBuff <- fread("400km_NoBuffer.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
fsh_400_noBuff <- st_as_sf(fsh_400_noBuff, coords=c("X", "Y"))
st_crs(fsh_400_noBuff) <- crs.1

grd_400$noBuff400 <- lengths(st_intersects(grd_400, fsh_400_noBuff))


# 800 km resolution
# fsh_800_noBuff <- fread("800km_NoBuffer.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
# fsh_800_noBuff <- st_as_sf(fsh_800_noBuff, coords=c("X", "Y"))
# st_crs(fsh_800_noBuff) <- 102008
# 
# grd_800$noBuff800 <- lengths(st_intersects(grd_800, fsh_800_noBuff))
# 
# fsh_800_Buff <- fread("800km_Buffer.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
# fsh_800_Buff <- st_as_sf(fsh_800_Buff, coords=c("X", "Y"))
# st_crs(fsh_800_Buff) <- 102008
# 
# grd_800$Buff800 <- lengths(st_intersects(grd_800, fsh_800_Buff))

############################################  
# read in and merge occurrence information #
############################################
ebut <- as_tibble(fread("ebut_intersections.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)) 
ebut <- ebut %>% 
  mutate(year=as.numeric(str_extract(Date.Observed, "^\\d{4}"))) %>%
  dplyr::select(Family, species, X, Y, year, inRM) %>% 
  mutate(basis="HUMAN_OBSERVATION", family=str_to_sentence(Family)) %>%
  dplyr::select(-Family)

idig <- as_tibble(fread("idig_intersections.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)) 
idig <- idig %>%
  mutate(year=as.numeric(str_extract(dwc.eventDate, "^\\d{4}"))) %>%
  dplyr::select(dwc.family, species, X, Y, inRM, year, basis=dwc.basisOfRecord) %>%
  mutate(family=str_to_sentence(dwc.family)) %>%
  dplyr::select(-dwc.family)

gbif <- as_tibble(fread("total_gbif_intersections.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)) 
gbif <- gbif %>% 
  dplyr::select(family, species, X, Y, inRM, year, basis=basisOfRecord) %>%
  mutate(family=str_to_sentence(family))

occur <- rbind(ebut, idig, gbif)
taxa <- unique(dplyr::select(occur, family, species))

rm(list=c("ebut", "idig", "gbif"))

##################################
# update families in fishnetting #
##################################
# 100 km resolution
fsh_100_Buff <- fsh_100_Buff %>%
  left_join(taxa, by=c("scientificName"="species")) %>%
  dplyr::select(-scientificName) %>%
  mutate(family=str_to_sentence(family))

# 200 km resolution
fsh_200_noBuff <- fsh_200_noBuff %>%
  left_join(taxa, by=c("scientificName"="species")) %>%
  dplyr::select(-scientificName) %>%
  mutate(family=str_to_sentence(family))

# 400 km resolution
fsh_400_noBuff <- fsh_400_noBuff %>%
  left_join(taxa, by=c("scientificName"="species")) %>%
  dplyr::select(-scientificName) %>%
  mutate(family=str_to_sentence(family))

# 800 km resolution
# fsh_800_noBuff <- fsh_800_noBuff %>%
#   left_join(taxa, by=c("scientificName"="species")) %>%
#   dplyr::select(-scientificName) %>%
#   mutate(family=str_to_sentence(family))

###################################
# filter fishnets to family level #
###################################

# 100 meter resolution 
fsh_100_Buff_nym <- filter(fsh_100_Buff, family=="Nymphalidae")
grd_100$Buff100_nym <- lengths(st_intersects(grd_100, fsh_100_Buff_nym))

fsh_100_Buff_pap <- filter(fsh_100_Buff, family=="Papilionidae")
grd_100$Buff100_pap <- lengths(st_intersects(grd_100, fsh_100_Buff_pap))

fsh_100_Buff_lyc <- filter(fsh_100_Buff, family=="Lycaenidae")
grd_100$Buff100_lyc <- lengths(st_intersects(grd_100, fsh_100_Buff_lyc))

fsh_100_Buff_hes <- filter(fsh_100_Buff, family=="Hesperiidae")
grd_100$Buff100_hes <- lengths(st_intersects(grd_100, fsh_100_Buff_hes))

fsh_100_Buff_pie <- filter(fsh_100_Buff, family=="Pieridae")
grd_100$Buff100_pie <- lengths(st_intersects(grd_100, fsh_100_Buff_pie))

fsh_100_Buff_rio <- filter(fsh_100_Buff, family=="Riodinidae")
grd_100$Buff100_rio <- lengths(st_intersects(grd_100, fsh_100_Buff_rio))

# 200 meter resolution 
fsh_200_noBuff_nym <- filter(fsh_200_noBuff, family=="Nymphalidae")
grd_200$noBuff200_nym <- lengths(st_intersects(grd_200, fsh_200_noBuff_nym))

fsh_200_noBuff_pap <- filter(fsh_200_noBuff, family=="Papilionidae")
grd_200$noBuff200_pap <- lengths(st_intersects(grd_200, fsh_200_noBuff_pap))

fsh_200_noBuff_lyc <- filter(fsh_200_noBuff, family=="Lycaenidae")
grd_200$noBuff200_lyc <- lengths(st_intersects(grd_200, fsh_200_noBuff_lyc))

fsh_200_noBuff_hes <- filter(fsh_200_noBuff, family=="Hesperiidae")
grd_200$noBuff200_hes <- lengths(st_intersects(grd_200, fsh_200_noBuff_hes))

fsh_200_noBuff_pie <- filter(fsh_200_noBuff, family=="Pieridae")
grd_200$noBuff200_pie <- lengths(st_intersects(grd_200, fsh_200_noBuff_pie))

fsh_200_noBuff_rio <- filter(fsh_200_noBuff, family=="Riodinidae")
grd_200$noBuff200_rio <- lengths(st_intersects(grd_200, fsh_200_noBuff_rio))

# 400 meter resolution 
fsh_400_noBuff_nym <- filter(fsh_400_noBuff, family=="Nymphalidae")
grd_400$noBuff400_nym <- lengths(st_intersects(grd_400, fsh_400_noBuff_nym))

fsh_400_noBuff_pap <- filter(fsh_400_noBuff, family=="Papilionidae")
grd_400$noBuff400_pap <- lengths(st_intersects(grd_400, fsh_400_noBuff_pap))

fsh_400_noBuff_lyc <- filter(fsh_400_noBuff, family=="Lycaenidae")
grd_400$noBuff400_lyc <- lengths(st_intersects(grd_400, fsh_400_noBuff_lyc))

fsh_400_noBuff_hes <- filter(fsh_400_noBuff, family=="Hesperiidae")
grd_400$noBuff400_hes <- lengths(st_intersects(grd_400, fsh_400_noBuff_hes))

fsh_400_noBuff_pie <- filter(fsh_400_noBuff, family=="Pieridae")
grd_400$noBuff400_pie <- lengths(st_intersects(grd_400, fsh_400_noBuff_pie))

fsh_400_noBuff_rio <- filter(fsh_400_noBuff, family=="Riodinidae")
grd_400$noBuff400_rio <- lengths(st_intersects(grd_400, fsh_400_noBuff_rio))

# 800 meter resolution 
# fsh_800_Buff_nym <- filter(fsh_800_noBuff, family=="Nymphalidae")
# grd_800$Buff800_nym <- lengths(st_intersects(grd_800, fsh_800_Buff_nym))
# 
# fsh_800_Buff_pap <- filter(fsh_800_noBuff, family=="Papilionidae")
# grd_800$Buff800_pap <- lengths(st_intersects(grd_800, fsh_800_Buff_pap))
# 
# fsh_800_Buff_lyc <- filter(fsh_800_noBuff, family=="Lycaenidae")
# grd_800$Buff800_lyc <- lengths(st_intersects(grd_800, fsh_800_Buff_lyc))
# 
# fsh_800_Buff_hes <- filter(fsh_800_noBuff, family=="Hesperiidae")
# grd_800$Buff800_hes <- lengths(st_intersects(grd_800, fsh_800_Buff_hes))
# 
# fsh_800_Buff_pie <- filter(fsh_800_noBuff, family=="Pieridae")
# grd_800$Buff800_pie <- lengths(st_intersects(grd_800, fsh_800_Buff_pie))
# 
# fsh_800_Buff_rio <- filter(fsh_800_noBuff, family=="Riodinidae")
# grd_800$Buff800_rio <- lengths(st_intersects(grd_800, fsh_800_Buff_rio))

#############################################
# convert occurrence data to spatial object #
#############################################
# filter for records between 1950-2019
occur <- occur %>% filter(between(year, 1950, 2019))
occur <- st_as_sf(occur, coords=c("X", "Y"))
st_crs(occur) <- crs.1

###############################################
# apply filters for basis of record attribute #
###############################################
occur <- occur
occur_inRM <- filter(occur, inRM=="Yes")
occur_spec <- filter(occur_inRM, basis!="HUMAN_OBSERVATION", basis!="UNKNOWN", basis!="machineobservation", basis!="MACHINE_OBSERVATION")
occur_obse <- filter(occur_inRM, basis=="HUMAN_OBSERVATION")

nrow(filter(occur, inRM=="Yes"))/nrow(occur) # number in range overall

inRMYears <- as.data.frame(occur) %>% filter(inRM=="Yes") %>%
  group_by(year) %>% summarise(n=n())

outRMYears <- as.data.frame(occur) %>% filter(inRM=="No") %>%
  group_by(year) %>% summarise(n=n())

totalRMYears <- merge(outRMYears, inRMYears, by="year") %>%
  filter(between(year, 1950, 2019)) %>% mutate(perc.n.in = n.y/(n.x+n.y))
plot(totalRMYears$year, totalRMYears$perc.n.in)
mean(totalRMYears$perc.n.in)
1.96*sd(totalRMYears$perc.n.in)/sqrt(nrow(totalRMYears))
mean(filter(totalRMYears, between(year, 2010, 2019))$perc.n.in)
1.96*sd(filter(totalRMYears, between(year, 2010, 2019))$perc.n.in)/sqrt(nrow(filter(totalRMYears, between(year, 2010, 2019))))

occur_nym <- filter(occur_inRM, family=="Nymphalidae")
occur_pap <- filter(occur_inRM, family=="Papilionidae")
occur_lyc <- filter(occur_inRM, family=="Lycaenidae")
occur_hes <- filter(occur_inRM, family=="Hesperiidae")
occur_pie <- filter(occur_inRM, family=="Pieridae")

occur_spec_nym <- filter(occur_spec, family=="Nymphalidae")
occur_spec_pap <- filter(occur_spec, family=="Papilionidae")
occur_spec_lyc <- filter(occur_spec, family=="Lycaenidae")
occur_spec_hes <- filter(occur_spec, family=="Hesperiidae")
occur_spec_pie <- filter(occur_spec, family=="Pieridae")

occur_obse_nym <- filter(occur_obse, family=="Nymphalidae")
occur_obse_pap <- filter(occur_obse, family=="Papilionidae")
occur_obse_lyc <- filter(occur_obse, family=="Lycaenidae")
occur_obse_hes <- filter(occur_obse, family=="Hesperiidae")
occur_obse_pie <- filter(occur_obse, family=="Pieridae")
# 
# occur_t1 <- filter(occur, between(year, 1950, 1969))
# occur_t2 <- filter(occur, between(year, 1970, 1989))
# occur_t3 <- filter(occur, between(year, 1990, 2009))
# occur_t4 <- filter(occur, between(year, 2010, 2019))

#################################
# count unique species in grids #
#################################

# 100 meter resolution
grd_100_allRich <- occur %>% st_join(grd_100) %>% group_by(FID) %>% summarise(n_speciesAll=n_distinct(species))
grd_100_inRMRich <- occur_inRM %>% st_join(grd_100) %>% group_by(FID) %>% summarise(n_speciesInRM=n_distinct(species))
grd_100_inRMspecRich <- occur_spec %>% st_join(grd_100) %>% group_by(FID) %>% summarise(n_speciesInRMSpec=n_distinct(species))
grd_100_inRMobseRich <- occur_obse %>% st_join(grd_100) %>% group_by(FID) %>% summarise(n_speciesInRMObse=n_distinct(species))

grd_100_inRMnym <- occur_nym %>% st_join(grd_100) %>% group_by(FID) %>% summarise(n_speciesInRMNym=n_distinct(species))
grd_100_inRMpap <- occur_pap %>% st_join(grd_100) %>% group_by(FID) %>% summarise(n_speciesInRMPap=n_distinct(species))
grd_100_inRMlyc <- occur_lyc %>% st_join(grd_100) %>% group_by(FID) %>% summarise(n_speciesInRMlyc=n_distinct(species))
grd_100_inRMhes <- occur_hes %>% st_join(grd_100) %>% group_by(FID) %>% summarise(n_speciesInRMhes=n_distinct(species))
grd_100_inRMpie <- occur_pie %>% st_join(grd_100) %>% group_by(FID) %>% summarise(n_speciesInRMpie=n_distinct(species))

grd_100_inRMspecnym <- occur_spec_nym %>% st_join(grd_100) %>% group_by(FID) %>% summarise(n_speciesInRMspecNym=n_distinct(species))
grd_100_inRMspecpap <- occur_spec_pap %>% st_join(grd_100) %>% group_by(FID) %>% summarise(n_speciesInRMspecPap=n_distinct(species))
grd_100_inRMspeclyc <- occur_spec_lyc %>% st_join(grd_100) %>% group_by(FID) %>% summarise(n_speciesInRMspeclyc=n_distinct(species))
grd_100_inRMspeches <- occur_spec_hes %>% st_join(grd_100) %>% group_by(FID) %>% summarise(n_speciesInRMspeches=n_distinct(species))
grd_100_inRMspecpie <- occur_spec_pie %>% st_join(grd_100) %>% group_by(FID) %>% summarise(n_speciesInRMspecpie=n_distinct(species))

grd_100_inRMobsenym <- occur_obse_nym %>% st_join(grd_100) %>% group_by(FID) %>% summarise(n_speciesInRMobseNym=n_distinct(species))
grd_100_inRMobsepap <- occur_obse_pap %>% st_join(grd_100) %>% group_by(FID) %>% summarise(n_speciesInRMobsePap=n_distinct(species))
grd_100_inRMobselyc <- occur_obse_lyc %>% st_join(grd_100) %>% group_by(FID) %>% summarise(n_speciesInRMobselyc=n_distinct(species))
grd_100_inRMobsehes <- occur_obse_hes %>% st_join(grd_100) %>% group_by(FID) %>% summarise(n_speciesInRMobsehes=n_distinct(species))
grd_100_inRMobsepie <- occur_obse_pie %>% st_join(grd_100) %>% group_by(FID) %>% summarise(n_speciesInRMobsepie=n_distinct(species))

# grd_100_t1 <- occur_t1 %>% st_join(grd_100) %>% group_by(FID) %>% summarise(n_speciesT1=n_distinct(species))
# grd_100_t2 <- occur_t2 %>% st_join(grd_100) %>% group_by(FID) %>% summarise(n_speciest2=n_distinct(species))
# grd_100_t3 <- occur_t3 %>% st_join(grd_100) %>% group_by(FID) %>% summarise(n_speciest3=n_distinct(species))
# grd_100_t4 <- occur_t4 %>% st_join(grd_100) %>% group_by(FID) %>% summarise(n_speciest4=n_distinct(species))

# merge back with 100 grid
grd_100 <- grd_100 %>% left_join(as.data.frame(grd_100_allRich), by="FID")
grd_100 <- grd_100 %>% left_join(as.data.frame(grd_100_inRMRich), by="FID")
grd_100 <- grd_100 %>% left_join(as.data.frame(grd_100_inRMspecRich), by="FID")
grd_100 <- grd_100 %>% left_join(as.data.frame(grd_100_inRMobseRich), by="FID")

grd_100 <- grd_100 %>% left_join(as.data.frame(grd_100_inRMnym), by="FID")
grd_100 <- grd_100 %>% left_join(as.data.frame(grd_100_inRMpap), by="FID")
grd_100 <- grd_100 %>% left_join(as.data.frame(grd_100_inRMlyc), by="FID")
grd_100 <- grd_100 %>% left_join(as.data.frame(grd_100_inRMhes), by="FID")
grd_100 <- grd_100 %>% left_join(as.data.frame(grd_100_inRMpie), by="FID")

grd_100 <- grd_100 %>% left_join(as.data.frame(grd_100_inRMspecnym), by="FID")
grd_100 <- grd_100 %>% left_join(as.data.frame(grd_100_inRMspecpap), by="FID")
grd_100 <- grd_100 %>% left_join(as.data.frame(grd_100_inRMspeclyc), by="FID")
grd_100 <- grd_100 %>% left_join(as.data.frame(grd_100_inRMspeches), by="FID")
grd_100 <- grd_100 %>% left_join(as.data.frame(grd_100_inRMspecpie), by="FID")

grd_100 <- grd_100 %>% left_join(as.data.frame(grd_100_inRMobsenym), by="FID")
grd_100 <- grd_100 %>% left_join(as.data.frame(grd_100_inRMobsepap), by="FID")
grd_100 <- grd_100 %>% left_join(as.data.frame(grd_100_inRMobselyc), by="FID")
grd_100 <- grd_100 %>% left_join(as.data.frame(grd_100_inRMobsehes), by="FID")
grd_100 <- grd_100 %>% left_join(as.data.frame(grd_100_inRMobsepie), by="FID")

# grd_100 <- grd_100 %>% left_join(as.data.frame(grd_100_t1), by="FID")
# grd_100 <- grd_100 %>% left_join(as.data.frame(grd_100_t2), by="FID")
# grd_100 <- grd_100 %>% left_join(as.data.frame(grd_100_t3), by="FID")
# grd_100 <- grd_100 %>% left_join(as.data.frame(grd_100_t4), by="FID")

# 200 meter resolution
grd_200_allRich <- occur %>% st_join(grd_200) %>% group_by(FID) %>% summarise(n_speciesAll=n_distinct(species))
grd_200_inRMRich <- occur_inRM %>% st_join(grd_200) %>% group_by(FID) %>% summarise(n_speciesInRM=n_distinct(species))
grd_200_inRMSpecRich <- occur_spec %>% st_join(grd_200) %>% group_by(FID) %>% summarise(n_speciesInRMSpec=n_distinct(species))
grd_200_inRMObseRich <- occur_obse %>% st_join(grd_200) %>% group_by(FID) %>% summarise(n_speciesInRMObse=n_distinct(species))

grd_200_inRMObsenym <- occur_obse_nym %>% st_join(grd_200) %>% group_by(FID) %>% summarise(n_speciesInRMObseNym=n_distinct(species))
grd_200_inRMObsepap <- occur_obse_pap %>% st_join(grd_200) %>% group_by(FID) %>% summarise(n_speciesInRMObsePap=n_distinct(species))
grd_200_inRMObselyc <- occur_obse_lyc %>% st_join(grd_200) %>% group_by(FID) %>% summarise(n_speciesInRMObselyc=n_distinct(species))
grd_200_inRMObsehes <- occur_obse_hes %>% st_join(grd_200) %>% group_by(FID) %>% summarise(n_speciesInRMObsehes=n_distinct(species))
grd_200_inRMObsepie <- occur_obse_pie %>% st_join(grd_200) %>% group_by(FID) %>% summarise(n_speciesInRMObsepie=n_distinct(species))

# grd_200_t1 <- occur_t1 %>% st_join(grd_200) %>% group_by(FID) %>% summarise(n_speciesT1=n_distinct(species))
# grd_200_t2 <- occur_t2 %>% st_join(grd_200) %>% group_by(FID) %>% summarise(n_speciest2=n_distinct(species))
# grd_200_t3 <- occur_t3 %>% st_join(grd_200) %>% group_by(FID) %>% summarise(n_speciest3=n_distinct(species))
# grd_200_t4 <- occur_t4 %>% st_join(grd_200) %>% group_by(FID) %>% summarise(n_speciest4=n_distinct(species))

# merge back with 200 grid
grd_200 <- grd_200 %>% left_join(as.data.frame(grd_200_allRich), by="FID")
grd_200 <- grd_200 %>% left_join(as.data.frame(grd_200_inRMRich), by="FID")
grd_200 <- grd_200 %>% left_join(as.data.frame(grd_200_inRMSpecRich), by="FID")
grd_200 <- grd_200 %>% left_join(as.data.frame(grd_200_inRMObseRich), by="FID")

grd_200 <- grd_200 %>% left_join(as.data.frame(grd_200_inRMObsenym), by="FID")
grd_200 <- grd_200 %>% left_join(as.data.frame(grd_200_inRMObsepap), by="FID")
grd_200 <- grd_200 %>% left_join(as.data.frame(grd_200_inRMObselyc), by="FID")
grd_200 <- grd_200 %>% left_join(as.data.frame(grd_200_inRMObsehes), by="FID")
grd_200 <- grd_200 %>% left_join(as.data.frame(grd_200_inRMObsepie), by="FID")

# grd_200 <- grd_200 %>% left_join(as.data.frame(grd_200_t1), by="FID")
# grd_200 <- grd_200 %>% left_join(as.data.frame(grd_200_t2), by="FID")
# grd_200 <- grd_200 %>% left_join(as.data.frame(grd_200_t3), by="FID")
# grd_200 <- grd_200 %>% left_join(as.data.frame(grd_200_t4), by="FID")

# 400 meter resolution
grd_400_allRich <- occur %>% st_join(grd_400) %>% group_by(FID) %>% summarise(n_speciesAll=n_distinct(species))
grd_400_inRMRich <- occur_inRM %>% st_join(grd_400) %>% group_by(FID) %>% summarise(n_speciesInRM=n_distinct(species))
grd_400_inRMSpecRich <- occur_spec %>% st_join(grd_400) %>% group_by(FID) %>% summarise(n_speciesInRMSpec=n_distinct(species))
grd_400_inRMObseRich <- occur_obse %>% st_join(grd_400) %>% group_by(FID) %>% summarise(n_speciesInRMObse=n_distinct(species))

grd_400_inRMObsenym <- occur_obse_nym %>% st_join(grd_400) %>% group_by(FID) %>% summarise(n_speciesInRMObseNym=n_distinct(species))
grd_400_inRMObsepap <- occur_obse_pap %>% st_join(grd_400) %>% group_by(FID) %>% summarise(n_speciesInRMObsePap=n_distinct(species))
grd_400_inRMObselyc <- occur_obse_lyc %>% st_join(grd_400) %>% group_by(FID) %>% summarise(n_speciesInRMObselyc=n_distinct(species))
grd_400_inRMObsehes <- occur_obse_hes %>% st_join(grd_400) %>% group_by(FID) %>% summarise(n_speciesInRMObsehes=n_distinct(species))
grd_400_inRMObsepie <- occur_obse_pie %>% st_join(grd_400) %>% group_by(FID) %>% summarise(n_speciesInRMObsepie=n_distinct(species))

# grd_400_t1 <- occur_t1 %>% st_join(grd_400) %>% group_by(FID) %>% summarise(n_speciesT1=n_distinct(species))
# grd_400_t2 <- occur_t2 %>% st_join(grd_400) %>% group_by(FID) %>% summarise(n_speciest2=n_distinct(species))
# grd_400_t3 <- occur_t3 %>% st_join(grd_400) %>% group_by(FID) %>% summarise(n_speciest3=n_distinct(species))
# grd_400_t4 <- occur_t4 %>% st_join(grd_400) %>% group_by(FID) %>% summarise(n_speciest4=n_distinct(species))

# merge back with 400 grid
grd_400 <- grd_400 %>% left_join(as.data.frame(grd_400_allRich), by="FID")
grd_400 <- grd_400 %>% left_join(as.data.frame(grd_400_inRMRich), by="FID")
grd_400 <- grd_400 %>% left_join(as.data.frame(grd_400_inRMSpecRich), by="FID")
grd_400 <- grd_400 %>% left_join(as.data.frame(grd_400_inRMObseRich), by="FID")

grd_400 <- grd_400 %>% left_join(as.data.frame(grd_400_inRMObsenym), by="FID")
grd_400 <- grd_400 %>% left_join(as.data.frame(grd_400_inRMObsepap), by="FID")
grd_400 <- grd_400 %>% left_join(as.data.frame(grd_400_inRMObselyc), by="FID")
grd_400 <- grd_400 %>% left_join(as.data.frame(grd_400_inRMObsehes), by="FID")
grd_400 <- grd_400 %>% left_join(as.data.frame(grd_400_inRMObsepie), by="FID")

# grd_400 <- grd_400 %>% left_join(as.data.frame(grd_400_t1), by="FID")
# grd_400 <- grd_400 %>% left_join(as.data.frame(grd_400_t2), by="FID")
# grd_400 <- grd_400 %>% left_join(as.data.frame(grd_400_t3), by="FID")
# grd_400 <- grd_400 %>% left_join(as.data.frame(grd_400_t4), by="FID")

# 800 meter resolution
# grd_800_allRich <- occur %>% st_join(grd_800) %>% group_by(FID) %>% summarise(n_speciesAll=n_distinct(species))
# grd_800_inRMRich <- occur_inRM %>% st_join(grd_800) %>% group_by(FID) %>% summarise(n_speciesInRM=n_distinct(species))
# grd_800_inRMSpecRich <- occur_spec %>% st_join(grd_800) %>% group_by(FID) %>% summarise(n_speciesInRMSpec=n_distinct(species))
# grd_800_inRMObseRich <- occur_obse %>% st_join(grd_800) %>% group_by(FID) %>% summarise(n_speciesInRMObse=n_distinct(species))
# 
# grd_800_inRMObsenym <- occur_obse_nym %>% st_join(grd_800) %>% group_by(FID) %>% summarise(n_speciesInRMObseNym=n_distinct(species))
# grd_800_inRMObsepap <- occur_obse_pap %>% st_join(grd_800) %>% group_by(FID) %>% summarise(n_speciesInRMObsePap=n_distinct(species))
# grd_800_inRMObselyc <- occur_obse_lyc %>% st_join(grd_800) %>% group_by(FID) %>% summarise(n_speciesInRMObselyc=n_distinct(species))
# grd_800_inRMObsehes <- occur_obse_hes %>% st_join(grd_800) %>% group_by(FID) %>% summarise(n_speciesInRMObsehes=n_distinct(species))
# grd_800_inRMObsepie <- occur_obse_pie %>% st_join(grd_800) %>% group_by(FID) %>% summarise(n_speciesInRMObsepie=n_distinct(species))
# 
# grd_800_t1 <- occur_t1 %>% st_join(grd_800) %>% group_by(FID) %>% summarise(n_speciesT1=n_distinct(species))
# grd_800_t2 <- occur_t2 %>% st_join(grd_800) %>% group_by(FID) %>% summarise(n_speciest2=n_distinct(species))
# grd_800_t3 <- occur_t3 %>% st_join(grd_800) %>% group_by(FID) %>% summarise(n_speciest3=n_distinct(species))
# grd_800_t4 <- occur_t4 %>% st_join(grd_800) %>% group_by(FID) %>% summarise(n_speciest4=n_distinct(species))
# 
# # merge back with 800 grid
# grd_800 <- grd_800 %>% left_join(as.data.frame(grd_800_allRich), by="FID")
# grd_800 <- grd_800 %>% left_join(as.data.frame(grd_800_inRMRich), by="FID")
# grd_800 <- grd_800 %>% left_join(as.data.frame(grd_800_inRMSpecRich), by="FID")
# grd_800 <- grd_800 %>% left_join(as.data.frame(grd_800_inRMObseRich), by="FID")
# 
# grd_800 <- grd_800 %>% left_join(as.data.frame(grd_800_inRMObsenym), by="FID")
# grd_800 <- grd_800 %>% left_join(as.data.frame(grd_800_inRMObsepap), by="FID")
# grd_800 <- grd_800 %>% left_join(as.data.frame(grd_800_inRMObselyc), by="FID")
# grd_800 <- grd_800 %>% left_join(as.data.frame(grd_800_inRMObsehes), by="FID")
# grd_800 <- grd_800 %>% left_join(as.data.frame(grd_800_inRMObsepie), by="FID")
# 
# grd_800 <- grd_800 %>% left_join(as.data.frame(grd_800_t1), by="FID")
# grd_800 <- grd_800 %>% left_join(as.data.frame(grd_800_t2), by="FID")
# grd_800 <- grd_800 %>% left_join(as.data.frame(grd_800_t3), by="FID")
# grd_800 <- grd_800 %>% left_join(as.data.frame(grd_800_t4), by="FID")

####################
# calculate ratios #
####################

# 100 km resolution
grd_100 <- grd_100 %>% mutate(logInRMOccur_Buff = n_speciesInRM/Buff100,
                              logInRMSpec_Buff = n_speciesInRMSpec/Buff100,
                              logInRMObse_Buff = n_speciesInRMObse/Buff100,
                              nymobsRatio = n_speciesInRMobseNym/Buff100_nym,
                              papobsRatio = n_speciesInRMobsePap/Buff100_pap,
                              lycobsRatio = n_speciesInRMobselyc/Buff100_lyc,
                              hesobsRatio = n_speciesInRMobsehes/Buff100_hes,
                              pieobsRatio = n_speciesInRMobsepie/Buff100_pie,
                              nymspecRatio = n_speciesInRMspecNym/Buff100_nym,
                              papspecRatio = n_speciesInRMspecPap/Buff100_pap,
                              lycspecRatio = n_speciesInRMspeclyc/Buff100_lyc,
                              hesspecRatio = n_speciesInRMspeches/Buff100_hes,
                              piespecRatio = n_speciesInRMspecpie/Buff100_pie,
                              nymRatio = n_speciesInRMNym/Buff100_nym,
                              papRatio = n_speciesInRMPap/Buff100_pap,
                              lycRatio = n_speciesInRMlyc/Buff100_lyc,
                              hesRatio = n_speciesInRMhes/Buff100_hes,
                              pieRatio = n_speciesInRMpie/Buff100_pie) %>%
                              # t1Ratio = n_speciesT1/Buff100,
                              # t2Ratio = n_speciest2/Buff100,
                              # t3Ratio = n_speciest3/Buff100,
                              # t4Ratio = n_speciest4/Buff100) %>%
  mutate(logInRMOccur_Buff = ifelse(logInRMOccur_Buff > 1, 1, logInRMOccur_Buff),
         logInRMSpec_Buff = ifelse(logInRMSpec_Buff > 1, 1, logInRMSpec_Buff),
         logInRMObse_Buff = ifelse(logInRMObse_Buff > 1, 1, logInRMObse_Buff),
         nymobsRatio = ifelse(nymobsRatio > 1, 1, nymobsRatio),
         papobsRatio = ifelse(papobsRatio > 1, 1, papobsRatio),
         lycobsRatio = ifelse(lycobsRatio > 1, 1, lycobsRatio),
         hesobsRatio = ifelse(hesobsRatio > 1, 1, hesobsRatio),
         pieobsRatio = ifelse(pieobsRatio > 1, 1, pieobsRatio),
         nymspecRatio = ifelse(nymspecRatio > 1, 1, nymspecRatio),
         papspecRatio = ifelse(papspecRatio > 1, 1, papspecRatio),
         lycspecRatio = ifelse(lycspecRatio > 1, 1, lycspecRatio),
         hesspecRatio = ifelse(hesspecRatio > 1, 1, hesspecRatio),
         piespecRatio = ifelse(piespecRatio > 1, 1, piespecRatio),
         nymRatio = ifelse(nymRatio > 1, 1, nymRatio),
         papRatio = ifelse(papRatio > 1, 1, papRatio),
         lycRatio = ifelse(lycRatio > 1, 1, lycRatio),
         hesRatio = ifelse(hesRatio > 1, 1, hesRatio),
         pieRatio = ifelse(pieRatio > 1, 1, pieRatio))
         # t1Ratio = ifelse(t1Ratio > 1, 1, t1Ratio),
         # t2Ratio = ifelse(t2Ratio > 1, 1, t2Ratio),
         # t3Ratio = ifelse(t3Ratio > 1, 1, t3Ratio),
         # t4Ratio = ifelse(t4Ratio > 1, 1, t4Ratio))

# 200 km resolution
grd_200 <- grd_200 %>% mutate(logAllOccur_NoBuff = n_speciesAll/noBuff200,
                              logInRMOccur_NoBuff = n_speciesInRM/noBuff200,
                              logInRMOccur_Buff = n_speciesInRM/noBuff200,
                              logInRMSpec_Buff = n_speciesInRMSpec/noBuff200,
                              logInRMObse_Buff = n_speciesInRMObse/noBuff200,
                              nymObsRatio = n_speciesInRMObseNym/noBuff200_nym,
                              papObsRatio = n_speciesInRMObsePap/noBuff200_pap,
                              lycObsRatio = n_speciesInRMObselyc/noBuff200_lyc,
                              hesObsRatio = n_speciesInRMObsehes/noBuff200_hes,
                              pieObsRatio = n_speciesInRMObsepie/noBuff200_pie) %>%
                              # t1Ratio = n_speciesT1/Buff200,
                              # t2Ratio = n_speciest2/Buff200,
                              # t3Ratio = n_speciest3/Buff200,
                              # t4Ratio = n_speciest4/Buff200) %>%
  mutate(logInRMOccur_Buff = ifelse(logInRMOccur_Buff > 1, 1, logInRMOccur_Buff),
         logInRMSpec_Buff = ifelse(logInRMSpec_Buff > 1, 1, logInRMSpec_Buff),
         logInRMObse_Buff = ifelse(logInRMObse_Buff > 1, 1, logInRMObse_Buff),
         nymObsRatio = ifelse(nymObsRatio > 1, 1, nymObsRatio),
         papObsRatio = ifelse(papObsRatio > 1, 1, papObsRatio),
         lycObsRatio = ifelse(lycObsRatio > 1, 1, lycObsRatio),
         hesObsRatio = ifelse(hesObsRatio > 1, 1, hesObsRatio),
         pieObsRatio = ifelse(pieObsRatio > 1, 1, pieObsRatio))
         # t1Ratio = ifelse(t1Ratio > 1, 1, t1Ratio),
         # t2Ratio = ifelse(t2Ratio > 1, 1, t2Ratio),
         # t3Ratio = ifelse(t3Ratio > 1, 1, t3Ratio),
         # t4Ratio = ifelse(t4Ratio > 1, 1, t4Ratio))

# 400 km resolution
grd_400 <- grd_400 %>% mutate(logAllOccur_NoBuff = n_speciesAll/noBuff400,
                              logInRMOccur_NoBuff = n_speciesInRM/noBuff400,
                              logInRMOccur_Buff = n_speciesInRM/noBuff400,
                              logInRMSpec_Buff = n_speciesInRMSpec/noBuff400,
                              logInRMObse_Buff = n_speciesInRMObse/noBuff400,
                              nymObsRatio = n_speciesInRMObseNym/noBuff400_nym,
                              papObsRatio = n_speciesInRMObsePap/noBuff400_pap,
                              lycObsRatio = n_speciesInRMObselyc/noBuff400_lyc,
                              hesObsRatio = n_speciesInRMObsehes/noBuff400_hes,
                              pieObsRatio = n_speciesInRMObsepie/noBuff400_pie) %>%
                              # t1Ratio = n_speciesT1/Buff400,
                              # t2Ratio = n_speciest2/Buff400,
                              # t3Ratio = n_speciest3/Buff400,
                              # t4Ratio = n_speciest4/Buff400) %>%
  mutate(logInRMOccur_Buff = ifelse(logInRMOccur_Buff > 1, 1, logInRMOccur_Buff),
         logInRMSpec_Buff = ifelse(logInRMSpec_Buff > 1, 1, logInRMSpec_Buff),
         logInRMObse_Buff = ifelse(logInRMObse_Buff > 1, 1, logInRMObse_Buff),
         nymObsRatio = ifelse(nymObsRatio > 1, 1, nymObsRatio),
         papObsRatio = ifelse(papObsRatio > 1, 1, papObsRatio),
         lycObsRatio = ifelse(lycObsRatio > 1, 1, lycObsRatio),
         hesObsRatio = ifelse(hesObsRatio > 1, 1, hesObsRatio),
         pieObsRatio = ifelse(pieObsRatio > 1, 1, pieObsRatio))
         # t1Ratio = ifelse(t1Ratio > 1, 1, t1Ratio),
         # t2Ratio = ifelse(t2Ratio > 1, 1, t2Ratio),
         # t3Ratio = ifelse(t3Ratio > 1, 1, t3Ratio),
         # t4Ratio = ifelse(t4Ratio > 1, 1, t4Ratio))

# # 800 km resolution
# grd_800 <- grd_800 %>% mutate(logAllOccur_NoBuff = n_speciesAll/noBuff800,
#                               logInRMOccur_NoBuff = n_speciesInRM/noBuff800,
#                               logInRMOccur_Buff = n_speciesInRM/Buff800,
#                               logInRMSpec_Buff = n_speciesInRMSpec/Buff800,
#                               logInRMObse_Buff = n_speciesInRMObse/Buff800,
#                               nymObsRatio = n_speciesInRMObseNym/Buff800_nym,
#                               papObsRatio = n_speciesInRMObsePap/Buff800_pap,
#                               lycObsRatio = n_speciesInRMObselyc/Buff800_lyc,
#                               hesObsRatio = n_speciesInRMObsehes/Buff800_hes,
#                               pieObsRatio = n_speciesInRMObsepie/Buff800_pie,
#                               t1Ratio = n_speciesT1/Buff800,
#                               t2Ratio = n_speciest2/Buff800,
#                               t3Ratio = n_speciest3/Buff800,
#                               t4Ratio = n_speciest4/Buff800) %>%
#   mutate(logInRMOccur_Buff = ifelse(logInRMOccur_Buff > 1, 1, logInRMOccur_Buff),
#          logInRMSpec_Buff = ifelse(logInRMSpec_Buff > 1, 1, logInRMSpec_Buff),
#          logInRMObse_Buff = ifelse(logInRMObse_Buff > 1, 1, logInRMObse_Buff),
#          nymObsRatio = ifelse(nymObsRatio > 1, 1, nymObsRatio),
#          papObsRatio = ifelse(papObsRatio > 1, 1, papObsRatio),
#          lycObsRatio = ifelse(lycObsRatio > 1, 1, lycObsRatio),
#          hesObsRatio = ifelse(hesObsRatio > 1, 1, hesObsRatio),
#          pieObsRatio = ifelse(pieObsRatio > 1, 1, pieObsRatio),
#          t1Ratio = ifelse(t1Ratio > 1, 1, t1Ratio),
#          t2Ratio = ifelse(t2Ratio > 1, 1, t2Ratio),
#          t3Ratio = ifelse(t3Ratio > 1, 1, t3Ratio),
#          t4Ratio = ifelse(t4Ratio > 1, 1, t4Ratio))

#######################
## Sample Covariates ##
#######################

grd_100 <- unique(grd_100)

# RCP 8.5 Climate Velocity into 2085
grd_100$rcp85_2080 <- extract(rcp8.5_2080, grd_100, fun=mean, na.rm=TRUE)

grd_100 <- grd_100 %>% mutate(clim.rank = percent_rank(rcp85_2080)) %>% 
  mutate(clim.cat = ifelse(clim.rank >= 0.95, "2", ifelse(clim.rank >= 0.80, "1", "0")))
grd_100$clim.cat <- as.character(grd_100$clim.cat)

grd_100pt <- st_centroid(grd_100) 
grd_100pt <- grd_100pt %>% mutate(x.crd = st_coordinates(grd_100pt)[,1], y.crd = st_coordinates(grd_100pt)[,2])

climate <- dplyr::select(grd_100pt, clim.cat, logInRMOccur_Buff)


grd_100pt <- grd_100pt %>% mutate(logInRMOccur_Buff=replace_na(logInRMOccur_Buff, 0))

# plot climate velocities on grid, highlight 90th and 95th percentiles
tiff("climateMap.tiff", units="cm", width=12.5, height=12.5, res=350)
ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  #geom_sf(grd_100, mapping=aes(fill=clim.cat), color=NA, alpha=0.7)+
  geom_point(grd_100pt, mapping=aes(x=x.crd, y=y.crd,
                                    color=clim.cat,
                                    size=1-logInRMOccur_Buff), 
             alpha=0.5, shape=15)+
  scale_color_manual(values=c("#3c9ab2", "#e8c927", "#f22300"))+
  #scale_color_manual(values=c("black", "white"))+
  scale_size_continuous(range=c(0.01, 1.99999))+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void() + theme(legend.position = "none")
dev.off()

climate <- filter(climate, !is.na(rcp85_2080))

# get percent sampled over 80% for each climate percentile.
1-nrow(filter(climate, clim.cat==1, logInRMOccur_Buff >= 0.80))/nrow(filter(climate, clim.cat==1))
1-nrow(filter(climate, clim.cat==2, logInRMOccur_Buff >= 0.80))/nrow(filter(climate, clim.cat==2))

# human footprint and majority biome
library(exactextractr)
getMode <- function(x,...){
  uniqx <- unique(x)
  uniqx <- uniqx[!is.na(uniqx)]
  uniqx[which.max(tabulate(match(x, uniqx)))]
}


nature <- fread("covariates/PA_grid/PA_grid_area.csv")
grd_100$ID <- paste(grd_100$X, grd_100$Y)
grd_100 <- grd_100 %>% left_join(nature, by="ID")

grd_100$footprint <- exact_extract(footprint, grd_100, fun="mean")
grd_100$biome <- exact_extract(biomes, grd_100, fun="majority")

fit.all <- lm(logInRMOccur_Buff~footprint+grid_area, data=grd_100)
fit.spec <- lm(logInRMSpec_Buff~footprint+grid_area, data=grd_100)
fit.obse <- lm(logInRMObse_Buff~footprint+grid_area, data=grd_100)

fit.all.1 <- lm(logInRMOccur_Buff~footprint, data=grd_100)
fit.spec.1 <- lm(logInRMSpec_Buff~footprint, data=grd_100)
fit.obse.1 <- lm(logInRMObse_Buff~footprint, data=grd_100)

extractAIC(fit.all)-extractAIC(fit.all.1)

extractAIC(fit.spec)-extractAIC(fit.spec.1)

extractAIC(fit.obse)-extractAIC(fit.obse.1)

t.test(grd_100$logInRMSpec_Buff, grd_100$logInRMObse_Buff)
sd(na.omit(grd_100$logInRMObse_Buff))/sqrt(length(na.omit(grd_100$logInRMObse_Buff)))
sd(na.omit(grd_100$logInRMSpec_Buff))/sqrt(length(na.omit(grd_100$logInRMSpec_Buff)))

# get summary stats for biomes
grd_100$biome <- round(grd_100$biome, digits=0)
biom.st <- as.data.frame(dplyr::select(grd_100, biome, logInRMOccur_Buff, logInRMObse_Buff, logInRMSpec_Buff,
                                       nymobsRatio, papobsRatio, lycobsRatio, hesobsRatio, pieobsRatio)) %>%
  mutate(biome=floor(biome))

biom.fl <- biom.st %>% group_by(biome) %>%
  summarise(meanAll=mean(logInRMOccur_Buff, na.rm=TRUE), meanobse=mean(logInRMObse_Buff, na.rm=TRUE), meanSpec=mean(logInRMSpec_Buff, na.rm=TRUE),
            meanNym=mean(nymobsRatio, na.rm=TRUE), meanPap=mean(papobsRatio, na.rm=TRUE), meanLyc=mean(lycobsRatio, na.rm=TRUE), meanHes=mean(hesobsRatio, na.rm=TRUE),
            meanPie=mean(pieobsRatio, na.rm=TRUE),
            sdAll=sd(logInRMOccur_Buff, na.rm=TRUE), sdobse=sd(logInRMObse_Buff, na.rm=TRUE), sdSpec=sd(logInRMSpec_Buff, na.rm=TRUE),
            sdNym=sd(nymobsRatio, na.rm=TRUE), sdPap=sd(papobsRatio, na.rm=TRUE), sdLyc=sd(lycobsRatio, na.rm=TRUE), sdHes=sd(hesobsRatio, na.rm=TRUE),
            sdPie=sd(pieobsRatio, na.rm=TRUE), count=n()) %>%
  mutate(colorAll = ifelse(meanAll >= 0.8, "> 80%", ifelse(meanAll >= 0.5, "50% <= x < 80%", "< 50%")),
         colorobs = ifelse(meanobse >= 0.8, "> 80%", ifelse(meanobse >= 0.5, "50% <= x < 80%", "< 50%")),
         colorSpe = ifelse(meanSpec >= 0.8, "> 80%", ifelse(meanSpec >= 0.5, "50% <= x < 80%", "< 50%")))

biom.st <- biom.st %>% group_by(biome) %>%
  summarise(meanAll=mean(logInRMOccur_Buff, na.rm=TRUE), meanobse=mean(logInRMObse_Buff, na.rm=TRUE), meanSpec=mean(logInRMSpec_Buff, na.rm=TRUE),
            meanNym=mean(nymobsRatio, na.rm=TRUE), meanPap=mean(papobsRatio, na.rm=TRUE), meanLyc=mean(lycobsRatio, na.rm=TRUE), meanHes=mean(hesobsRatio, na.rm=TRUE),
            meanPie=mean(pieobsRatio, na.rm=TRUE),
            sdAll=sd(logInRMOccur_Buff, na.rm=TRUE), sdobse=sd(logInRMObse_Buff, na.rm=TRUE), sdSpec=sd(logInRMSpec_Buff, na.rm=TRUE),
            sdNym=sd(nymobsRatio, na.rm=TRUE), sdPap=sd(papobsRatio, na.rm=TRUE), sdLyc=sd(lycobsRatio, na.rm=TRUE), sdHes=sd(hesobsRatio, na.rm=TRUE),
            sdPie=sd(pieobsRatio, na.rm=TRUE), count=n()) %>%
  mutate(colorAll = ifelse(meanAll >= 0.8, ">= 80%", ifelse(meanAll >= 0.5, "50% <= x < 80%", "< 50%")),
         colorobs = ifelse(meanobse >= 0.8, ">= 80%", ifelse(meanobse >= 0.5, "50% <= x < 80%", "< 50%")),
         colorSpe = ifelse(meanSpec >= 0.8, ">= 80%", ifelse(meanSpec >= 0.5, "50% <= x < 80%", "< 50%"))) %>%
  filter(count > 10) %>%
  mutate(coord=row_number(), coordAll=1, coordobse=3, coordSpec=5, coordNym=8, coordPap=10, coordLyc=12, coordHes=14, coordPie=16)
  

library(ggforce)

tiff("circles.tiff", units="cm", width=8.3, height=8.3, res=350)
ggplot(biom.st)+
  geom_circle(aes(x0=coordAll, y0=coord, r=meanAll/2+sdAll/2, color=colorAll, fill=colorAll), show.legend=FALSE)+
  geom_circle(aes(x0=coordAll, y0=coord, r=meanAll/2-sdAll/2, color=colorAll), fill="white")+
  geom_circle(aes(x0=coordAll, y0=coord, r=meanAll/2), color="black", fill=NA)+
  geom_circle(aes(x0=coordobse, y0=coord, r=meanobse/2+sdobse/2, color=colorobs, fill=colorobs), show.legend=FALSE)+
  geom_circle(aes(x0=coordobse, y0=coord, r=meanobse/2-sdobse/2, color=colorobs), fill="white")+
  geom_circle(aes(x0=coordobse, y0=coord, r=meanobse/2), color="black", fill=NA)+
  geom_circle(aes(x0=coordSpec, y0=coord, r=meanSpec/2+sdSpec/2, color=colorSpe, fill=colorSpe), show.legend=FALSE)+
  geom_circle(aes(x0=coordSpec, y0=coord, r=meanSpec/2-sdSpec/2, color=colorSpe), fill="white")+
  geom_circle(aes(x0=coordSpec, y0=coord, r=meanSpec/2), color="black", fill=NA)+
  # geom_circle(aes(x0=coordHes, y0=coord, r=meanHes/2+sdHes/2), color="grey", fill="grey")+
  # geom_circle(aes(x0=coordHes, y0=coord, r=meanHes/2-sdHes/2), color="grey", fill="white")+
  # geom_circle(aes(x0=coordHes, y0=coord, r=meanHes/2), bolor="black", fill=NA)+
  # geom_circle(aes(x0=coordLyc, y0=coord, r=meanLyc/2+sdLyc/2), color="grey", fill="grey")+
  # geom_circle(aes(x0=coordLyc, y0=coord, r=meanLyc/2-sdLyc/2), color="grey", fill="white")+
  # geom_circle(aes(x0=coordLyc, y0=coord, r=meanLyc/2), bolor="black", fill=NA)+
  # geom_circle(aes(x0=coordNym, y0=coord, r=meanNym/2+sdNym/2), color="grey", fill="grey")+
  # geom_circle(aes(x0=coordNym, y0=coord, r=meanNym/2-sdNym/2), color="grey", fill="white")+
  # geom_circle(aes(x0=coordNym, y0=coord, r=meanNym/2), bolor="black", fill=NA)+
  # geom_circle(aes(x0=coordPap, y0=coord, r=meanPap/2+sdPap/2), color="grey", fill="grey")+
  # geom_circle(aes(x0=coordPap, y0=coord, r=meanPap/2-sdPap/2), color="grey", fill="white")+
  # geom_circle(aes(x0=coordPap, y0=coord, r=meanPap/2), bolor="black", fill=NA)+
  # geom_circle(aes(x0=coordPie, y0=coord, r=meanPie/2+sdPie/2), color="grey", fill="grey")+
  # geom_circle(aes(x0=coordPie, y0=coord, r=meanPie/2-sdPie/2), color="grey", fill="white")+
  # geom_circle(aes(x0=coordPie, y0=coord, r=meanPie/2), color="black", fill=NA)+
  #geom_circle(aes(x0=c(0,1,2,3), y0=c(0,0,0,0), r=c(0.25/2, 0.5/2, 0.75/2, 1/2)))+
scale_fill_manual(values=c("firebrick2", "seagreen3", "goldenrod1"))+
scale_color_manual(values=c("firebrick2", "seagreen3", "goldenrod1"))+
xlim(0, 11) + ylim(0, 11)+
  theme_void()+
  theme(legend.position = "none") + labs(tag="(a)")
dev.off()

b.ply <- st_read("covariates/WWFBiomes.shp")
b.ply <- st_transform(b.ply, crs.1)
b.ply <- b.ply %>% left_join(biom.st, by=c("BIOME" = "biome"))

tiff("circles_legend.tiff", units="cm", width=8.3, height=8.3, res=350)
ggplot()+
  geom_circle(aes(x0=c(1,2,3,4), y0=c(1,1,1,1), r=c(0.25/2, 0.5/2, 0.75/2, 1/2)))+
  xlim(0, 11) + ylim(0, 11)+
  theme_void() + labs(tag="   ")
dev.off()


# ggplot()+
#   geom_sf(data=land, fill="grey", color=NA)+
#   geom_sf(b.ply, mapping=aes(fill=as.factor(BIOME)), color="grey", size=0.05)+
#   scale_fill_manual(values=c("#839791", "#839791", "#839791",
#                              "#839791", "#839791", "#839791",
#                              "#F28F3B", "#F28F3B", "#F28F3B",
#                              "#F28F3B", "#575761", "#F4D6CC",
#                              "#F4D6CC", "#F4D6CC", NA, NA))+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()

tiff("biome_map.tiff", units="in", width=8.3, height=8.3, res=350)
ggplot()+
  geom_sf(data=land, fill="grey", color=NA)+
  geom_sf(b.ply, mapping=aes(fill=colorAll), color=NA, size=0.05)+
  scale_fill_manual(values=c("firebrick2", "seagreen3", "goldenrod1"))+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void()+labs(tag="(b)") 
dev.off()

###################
# yearly analysis #
###################
yearly <- data.frame(startYear <- integer(),
                     endYear <- integer(),
                     specCom <- double(),
                     obsCom <- double(),
                     specAve <- double(),
                     obsAve <- double())

grd_100c <- grd_100

# biyearly
i = 1950
while(i < 2019){
  print(paste("Processing years: ", i, "-", i+1))
  grd_100c <- grd_100
  
  yearly[i-1949,1] <- i
  yearly[i-1949,2] <- i+1
  
  occur_spec <- filter(occur_inRM, between(year, i, i+1), basis!="HUMAN_OBSERVATION")
  occur_obse <- filter(occur_inRM, between(year, i, i+1), basis=="HUMAN_OBSERVATION")
  
  grd_100_biyearlyspec <- occur_spec %>% st_join(grd_100c) %>% group_by(FID) %>% summarise(n_speciesInRMSpec=n_distinct(species))
  grd_100_biyearlyobs <- occur_obse %>% st_join(grd_100c) %>% group_by(FID) %>% summarise(n_speciesInRMObse=n_distinct(species))
  
  grd_100c <- grd_100c %>% left_join(as.data.frame(grd_100_biyearlyspec), by="FID")
  grd_100c <- grd_100c %>% left_join(as.data.frame(grd_100_biyearlyobs), by="FID")
  
  grd_100c <- grd_100c %>% mutate(logInRMSpec_Buff = n_speciesInRMSpec.y/Buff100,
                                logInRMObse_Buff = n_speciesInRMObse.y/Buff100) %>%
    mutate(logInRMSpec_Buff = ifelse(logInRMSpec_Buff > 1, 1, logInRMSpec_Buff),
           logInRMObse_Buff = ifelse(logInRMObse_Buff > 1, 1, logInRMObse_Buff))
  
  yearly[i-1949,3] <- nrow(filter(grd_100c, logInRMSpec_Buff >= 0.8))
  yearly[i-1949,4] <- nrow(filter(grd_100c, logInRMObse_Buff >= 0.8))
  yearly[i-1949,5] <- mean(na.omit(grd_100c$logInRMSpec_Buff))
  yearly[i-1949,6] <- mean(na.omit(grd_100c$logInRMObse_Buff))

  i = i+2
}

yearly <- na.omit(yearly)

yearly <- yearly %>% mutate(total = specCom....double..+obsCom....double..+
                              (nrow(grd_100)-specCom....double..+obsCom....double..),
                            percSpec = specCom....double../total,
                            obsCom = obsCom....double../total,
                            unsamp = (nrow(grd_100)-specCom....double..+obsCom....double..)/total,
                            ratioS = specCom....double../(specCom....double..+obsCom....double..),
                            ratioO = obsCom....double../(specCom....double..+obsCom....double..))

# ggplot()+
#   geom_line(aes(x=yearly$startYear....integer.., y=yearly$specCom....double..), color="red4", lwd=1.1)+
#   geom_point(aes(x=yearly$startYear....integer.., y=yearly$specCom....double..), color="red4", fill="white", pch=21, size=2)+
#   geom_line(aes(x=yearly$startYear....integer.., y=yearly$obsCom....double..), color="royalblue4", lwd=1.1)+
#   geom_point(aes(x=yearly$startYear....integer.., y=yearly$obsCom....double..), color="royalblue4", fill="white", pch=21, size=2)+
#   labs(y="Number of 100km Cells Over 80% Complete", x="Years of Sampling")+
#   scale_fill_manual(values=c("red4", "royalblue4"), labels=c("Specimen Data", "Observation Data"))+
#   theme_minimal()

# ggplot()+
#   geom_line(aes(x=yearly$startYear....integer.., y=yearly$specAve....double..), color="red4", lwd=1.1)+
#   geom_point(aes(x=yearly$startYear....integer.., y=yearly$specAve....double..), color="red4", fill="white", pch=21, size=2)+
#   geom_line(aes(x=yearly$startYear....integer.., y=yearly$obsAve....double..), color="royalblue4", lwd=1.1)+
#   geom_point(aes(x=yearly$startYear....integer.., y=yearly$obsAve....double..), color="royalblue4", fill="white", pch=21, size=2)+
#   labs(y="Average Sampling Completeness Across All Cells", x="Years of Sampling")+
#   theme_minimal()

yearly.long <- yearly %>% dplyr::select(startYear....integer.., specCom....double.., obsCom....double..) %>% reshape2::melt(id.vars="startYear....integer..")

# Supplemental 1
p1 <- ggplot()+
  geom_bar(aes(y=value, x=startYear....integer.., fill=variable), data=yearly.long, stat="identity", width = 1.8)+
  theme_minimal()+labs(x="Sampling Start Year", y="Count of Cells over 80% Complete")+
  scale_x_continuous(breaks=seq(1950, 2019, 2))+
  scale_fill_manual(values=c("red3", "royalblue3"), labels = c("Museum Specimens", "Community Observations"))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2), panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), legend.title = element_blank(),
        legend.position = "none", axis.text=element_text(size=12)) + labs(tag="(b)")

sup <- as.data.frame(occur_inRM) %>% mutate(basis=recode(basis, preservedspecimen="PRESERVED_SPECIMEN",
                                                         OBSERVATION="HUMAN_OBSERVATION",
                                                         MATERIAL_SAMPLE="PRESERVED_SPECIMEN"))
sup <-  sup %>% filter(basis %in% c("HUMAN_OBSERVATION",
                                    "PRESERVED_SPECIMEN")) %>% group_by(year, basis) %>% tally()

p2 <- ggplot()+
  geom_bar(aes(y=n, x=year, fill=basis), data=sup, stat="identity", width=1)+
  theme_minimal()+labs(x="Sampling Year", y="Number of Records")+
  scale_fill_manual(values=c("royalblue3", "red3"), labels = c("Museum Specimens", "Community Observations"))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.2), panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), legend.title = element_blank(),
        legend.position = "none", axis.text=element_text(size=12))+
  scale_y_continuous(label=comma) + labs(tag="(a)")
  
tiff("Supp1.tiff", units="in", width=7, height=6, res=350)
grid.arrange(p2, p1, nrow=2)
dev.off()

library(ggsignif)
# Between family ANOVAs
fams <- as.data.frame(grd_100) %>% dplyr::select(FID, nymRatio, papRatio, lycRatio, hesRatio, pieRatio)
fams <- fams %>% reshape2::melt(id.vars="FID")
aov.all <- aov(value~variable, data=fams)
TukeyHSD(aov.all)

fams.obs <- as.data.frame(grd_100) %>% dplyr::select(FID, nymobsRatio, papobsRatio, lycobsRatio, hesobsRatio, pieobsRatio)
fams.obs <- fams.obs %>% reshape2::melt(id.vars="FID")
aov.obs <- aov(value~variable, data=fams.obs)
TukeyHSD(aov.obs)

fams.spec <- as.data.frame(grd_100) %>% dplyr::select(FID, nymspecRatio, papspecRatio, lycspecRatio, hesspecRatio, piespecRatio)
fams.spec <- fams.spec %>% reshape2::melt(id.vars="FID")
aov.spec <- aov(value~variable, data=fams.spec)
TukeyHSD(aov.spec)

fams.chi <- as.data.frame(grd_100) %>% filter(logInRMOccur_Buff >= 0.5) %>% dplyr::select(FID, nymspecRatio, papspecRatio, lycspecRatio, hesspecRatio, piespecRatio,
                                                 nymobsRatio, papobsRatio, lycobsRatio, hesobsRatio, pieobsRatio)

# report counts for chi-square
nrow(filter(fams.chi, nymspecRatio > 0.5))
nrow(filter(fams.chi, nymobsRatio > 0.5))

nrow(filter(fams.chi, papspecRatio > 0.5))
nrow(filter(fams.chi, papobsRatio > 0.5))

nrow(filter(fams.chi, lycspecRatio > 0.5))
nrow(filter(fams.chi, lycobsRatio > 0.5))

nrow(filter(fams.chi, hesspecRatio > 0.5))
nrow(filter(fams.chi, hesobsRatio > 0.5))

nrow(filter(fams.chi, piespecRatio > 0.5))
nrow(filter(fams.chi, pieobsRatio > 0.5))


#write.csv(fams.chi, "chi_raw.csv")
library(chisq.posthoc.test)
fams.chi <- read.csv("chi_raw.csv")
rownames(fams.chi) <- fams.chi[,1]
fams.chi <- fams.chi[,-1]
chisq.test(fams.chi)
chisq.posthoc.test(fams.chi, method="bonferroni", round=6)

fams.chi <- read.csv("chi_plots.csv")
chi.nym <- ggplot()+
  geom_bar(data=filter(fams.chi, Family=="Nym"), aes(x=Type, y=Value, col=Type, fill=Type), stat="identity", 
           show.legend=FALSE)+xlab("")+ylab("")+scale_fill_manual(values=c("royalblue3", "red3"))+
  scale_color_manual(values=c("royalblue3", "red3"))+ylim(0,600)+
  theme_minimal() + labs(tag="   ")
chi.pap <- ggplot()+
  geom_bar(data=filter(fams.chi, Family=="Pap"), aes(x=Type, y=Value), stat="identity", 
           show.legend=FALSE)+xlab("")+ylab("")+ylim(0,600)+
  theme_minimal()+ labs(tag="   ")
chi.lyc <- ggplot()+
  geom_bar(data=filter(fams.chi, Family=="Lyc"), aes(x=Type, y=Value), stat="identity", 
           show.legend=FALSE)+xlab("")+ylab("")+ylim(0,600)+
  theme_minimal()+ labs(tag="   ")
chi.hes <- ggplot()+
  geom_bar(data=filter(fams.chi, Family=="Hes"), aes(x=Type, y=Value), stat="identity", 
           show.legend=FALSE)+xlab("")+ylab("")+ylim(0,600)+
  theme_minimal() + labs(tag="(d)")
chi.pie <- ggplot()+
  geom_bar(data=filter(fams.chi, Family=="Pie"), aes(x=Type, y=Value, col=Type, fill=Type), stat="identity", 
           show.legend=FALSE)+xlab("")+ylab("")+scale_fill_manual(values=c("royalblue3", "red3"))+
  scale_color_manual(values=c("royalblue3", "red3"))+ylim(0,600)+
  theme_minimal() + labs(tag="   ")

tiff("tukeys.tiff", height = 8, width = 8, res=400)
ggplot(fams, aes(x=variable, y=value))+
  geom_boxplot(fill="grey80", color="black")+
  scale_x_discrete() + xlab("Family")+
  ylab("Completeness")+
  geom_signif(comparisons=list(c("nymRatio", "papRatio"),
                               c("nymRatio", "lycRatio"),
                               c("nymRatio", "hesRatio"),
                               c("nymRatio", "pieRatio"),
                               c("lycRatio", "papRatio"),
                               c("hesRatio", "papRatio"),
                               c("pieRatio", "papRatio"),
                               c("hesRatio", "lycRatio"),
                               c("pieRatio", "lycRatio"),
                               c("pieRatio", "hesRatio")),
              map_signif_level = TRUE,
              tip_length=0,
              y_position = c(-0.1, -0.2, -0.3, -0.4, 1.1, 1.2, 1.3, -0.6, -0.7, 1.5)) +
  theme_minimal()
dev.off()

p2 <- ggplot(fams.obs, aes(x=variable, y=value))+
  geom_boxplot(fill="grey80", color="black")+
  scale_x_discrete() + xlab("Family")+
  ylab("Completeness")+
  geom_signif(comparisons=list(c("nymobsRatio", "papobsRatio"),
                               c("nymobsRatio", "lycobsRatio"),
                               c("nymobsRatio", "hesobsRatio"),
                               c("nymobsRatio", "pieobsRatio"),
                               c("lycobsRatio", "papobsRatio"),
                               c("hesobsRatio", "papobsRatio"),
                               c("pieobsRatio", "papobsRatio"),
                               c("hesobsRatio", "lycobsRatio"),
                               c("pieobsRatio", "lycobsRatio"),
                               c("pieobsRatio", "hesobsRatio")),
              map_signif_level = TRUE,
              tip_length=0,
              y_position = c(-0.1, -0.2, -0.3, -0.4, 1.1, 1.2, 1.3, -0.6, -0.7, 1.5)) +
  theme_minimal()

p3 <- ggplot(fams.spec, aes(x=variable, y=value))+
  geom_boxplot(fill="grey80", color="black")+
  scale_x_discrete() + xlab("Family")+
  ylab("Completeness")+
  geom_signif(comparisons=list(c("nymspecRatio", "papspecRatio"),
                               c("nymspecRatio", "lycspecRatio"),
                               c("nymspecRatio", "hesspecRatio"),
                               c("nymspecRatio", "piespecRatio"),
                               c("lycspecRatio", "papspecRatio"),
                               c("hesspecRatio", "papspecRatio"),
                               c("piespecRatio", "papspecRatio"),
                               c("hesspecRatio", "lycspecRatio"),
                               c("piespecRatio", "lycspecRatio"),
                               c("piespecRatio", "hesspecRatio")),
              map_signif_level = TRUE,
              tip_length=0,
              y_position = c(-0.1, -0.2, -0.3, -0.4, 1.1, 1.2, 1.3, -0.6, -0.7, 1.5)) +
  theme_minimal()

grid.arrange(p1, p2, p3, ncol=3)
# decades
# i = 1950
# while(i < 2020){
#   print(paste("Processing years: ", i, "-", i+9))
#   
#   yearly[i-1949,]$startYear <- i
#   yearly[i-1949,]$endYear <- i+1
#   
#   grd_100_biyearlyspec <- occur_spec %>% st_join(grd_100) %>% group_by(FID) %>% summarise(n_speciesInRMSpec=n_distinct(species))
#   grd_100_biyearlyobs <- occur_obse %>% st_join(grd_100) %>% group_by(FID) %>% summarise(n_speciesInRMObse=n_distinct(species))
#   
#   yearly[i-1949,]$specCom <- mean(grd_100_biyearlyspec$n_speciesInRMSpec)
#   yearly[i-1949,]$obsCom <- mean(grd_100_biyearlyspec$n_speciesInRMObse)
#   
#   i = i+10
# }

##################
# Visualizations #
##################

# 100 km fishnet no buffer
# km100_noBuffer <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_100, mapping=aes(fill=noBuff100), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, midpoint=300)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# ggsave("km100_noBuffer.pdf", plot=km100_noBuffer)
# 
# # 100 km fishnet buffer
# km100_Buffer <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_100, mapping=aes(fill=log(Buff100)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, midpoint=2.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# ggsave("km100_Buffer.pdf", plot=km100_Buffer)
# 
# # visualize the occurrence richness values
# # 100 km all occurrences
# km100_allOccur <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_100, mapping=aes(fill=log(n_speciesAll)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, midpoint=2.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# ggsave("km100_allOccur.pdf", plot=km100_allOccur)
# 
# # 100 km in range occurrences
# km100_allOccurInRM <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_100, mapping=aes(fill=log(n_speciesInRM)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, midpoint=2.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# ggsave("km100_allOccurInRM.pdf", plot=km100_allOccurInRM)
# 
# # 100 km in range occurrences from specimens
# km100_specimensInRM <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_100, mapping=aes(fill=log(n_speciesInRMSpec)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, midpoint=2.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# ggsave("km100_specimensInRM.pdf", plot=km100_specimensInRM)
# 
# # 100 km in range occurrences from human observations
# km100_observationsInRM <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_100, mapping=aes(fill=log(n_speciesInRMObse)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, midpoint=2.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# ggsave("km100_observationsInRM.pdf", plot=km100_observationsInRM)
# 
# # 100 km all occurrences no buffer
# ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_100, mapping=aes(fill=log(logAllOccur_NoBuff)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # 100 km range occurrences no buffer
# ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_100, mapping=aes(fill=logInRMOccur_NoBuff), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()

###########################################################
##### ALL RECORDS IN RANGE MAPS WITH BUFFERING FIGURE #####
###########################################################
# 100 km resolution
# composite data
grd_100_p1 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(grd_100, mapping=aes(fill=logInRMOccur_Buff), color=NA, alpha=0.7, show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void() + labs(tag="(a)")

# museum specimens
grd_100_p2 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(grd_100, mapping=aes(fill=logInRMSpec_Buff), color=NA, alpha=0.7, show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void() + labs(tag="(b)")

# human observations
grd_100_p3 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(grd_100, mapping=aes(fill=logInRMObse_Buff), color=NA, alpha=0.7, show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void() + labs(tag="(c)")


# 200 km resolution
# composite data
grd_200_p1 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(grd_200, mapping=aes(fill=logInRMOccur_Buff), color=NA, alpha=0.7, show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void()+ labs(tag="   ")

# museum specimens
grd_200_p2 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(grd_200, mapping=aes(fill=logInRMSpec_Buff), color=NA, alpha=0.7, show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void()+ labs(tag="   ")

# human observations
grd_200_p3 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(grd_200, mapping=aes(fill=logInRMObse_Buff), color=NA, alpha=0.7, show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void()+ labs(tag="   ")


# 400 km resolution
# composite data
grd_400_p1 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(grd_400, mapping=aes(fill=logInRMOccur_Buff), color=NA, alpha=0.7, show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void()+ labs(tag="   ")

# museum specimens
grd_400_p2 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(grd_400, mapping=aes(fill=logInRMSpec_Buff), color=NA, alpha=0.7, show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void()+ labs(tag="   ")

# human observations
grd_400_p3 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(grd_400, mapping=aes(fill=logInRMObse_Buff), color=NA, alpha=0.7, show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void()+ labs(tag="   ")

tiff("Basis.tiff", units="cm", width=16.6, height=18, res=400)
grid.arrange(grd_100_p1, grd_100_p2, grd_100_p3,
             grd_200_p1, grd_200_p2, grd_200_p3,
             grd_400_p1, grd_400_p2, grd_400_p3, nrow=3, ncol=3)
dev.off()

# # 800 km resolution
# # composite data
# grd_800_p1 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_800, mapping=aes(fill=logInRMOccur_Buff), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # museum specimens
# grd_800_p2 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_800, mapping=aes(fill=logInRMSpec_Buff), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # human observations
# grd_800_p3 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_800, mapping=aes(fill=logInRMObse_Buff), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# g4 <- grid.arrange(grd_800_p1, grd_800_p2, grd_800_p3, nrow=3)

########################################################
##### FAMILY SPECIFIC SPECIES RICHNESS FROM RANGES #####
########################################################
# 100 km resolution #
# Hesperiidae
# grd_100_p1 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_100, mapping=aes(fill=log(Buff100_hes)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,6), midpoint=2.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Lycaenidae
# grd_100_p2 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_100, mapping=aes(fill=log(Buff100_lyc)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,6), midpoint=2.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Nymphalidae
# grd_100_p3 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_100, mapping=aes(fill=log(Buff100_nym)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,6), midpoint=2.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Papilionidae
# grd_100_p4 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_100, mapping=aes(fill=log(Buff100_pap)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,6), midpoint=2.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Pieridae
# grd_100_p5 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_100, mapping=aes(fill=log(Buff100_pie)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,6), midpoint=2.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# grid.arrange(grd_100_p1, grd_100_p2, grd_100_p3, grd_100_p4, grd_100_p5, nrow=5, ncol=1)

# # 200 km resolution #
# # Hesperiidae
# grd_200_p1 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_200, mapping=aes(fill=log(Buff200_hes)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,6), midpoint=2.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Lycaenidae
# grd_200_p2 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_200, mapping=aes(fill=log(Buff200_lyc)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,6), midpoint=2.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Nymphalidae
# grd_200_p3 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_200, mapping=aes(fill=log(Buff200_nym)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,6), midpoint=2.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Papilionidae
# grd_200_p4 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_200, mapping=aes(fill=log(Buff200_pap)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,6), midpoint=2.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Pieridae
# grd_200_p5 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_200, mapping=aes(fill=log(Buff200_pie)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,6), midpoint=2.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# grid.arrange(grd_200_p1, grd_200_p2, grd_200_p3, grd_200_p4, grd_200_p5, nrow=3, ncol=2)
# 
# # 400 km resolution #
# # Hesperiidae
# grd_400_p1 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_400, mapping=aes(fill=log(Buff400_hes)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,6), midpoint=2.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Lycaenidae
# grd_400_p2 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_400, mapping=aes(fill=log(Buff400_lyc)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,6), midpoint=2.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Nymphalidae
# grd_400_p3 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_400, mapping=aes(fill=log(Buff400_nym)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,6), midpoint=2.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Papilionidae
# grd_400_p4 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_400, mapping=aes(fill=log(Buff400_pap)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,6), midpoint=2.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Pieridae
# grd_400_p5 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_400, mapping=aes(fill=log(Buff400_pie)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,6), midpoint=2.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# grid.arrange(grd_400_p1, grd_400_p2, grd_400_p3, grd_400_p4, grd_400_p5, nrow=3, ncol=2)
# 
# # 800 km resolution #
# # Hesperiidae
# grd_800_p1 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_800, mapping=aes(fill=log(Buff800_hes)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,6), midpoint=2.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Lycaenidae
# grd_800_p2 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_800, mapping=aes(fill=log(Buff800_lyc)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,6), midpoint=2.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Nymphalidae
# grd_800_p3 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_800, mapping=aes(fill=log(Buff800_nym)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,6), midpoint=2.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Papilionidae
# grd_800_p4 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_800, mapping=aes(fill=log(Buff800_pap)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,6), midpoint=2.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Pieridae
# grd_800_p5 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_800, mapping=aes(fill=log(Buff800_pie)), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,6), midpoint=2.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# grid.arrange(grd_800_p1, grd_800_p2, grd_800_p3, grd_800_p4, grd_800_p5, nrow=3, ncol=2)

##### FAMILIES LEVEL COMPLETENESS #####
# 100 km resolution observations only
# Hesperiidae
grd_100_p1_obse <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(grd_100, mapping=aes(fill=hesobsRatio), color=NA, alpha=0.7, show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void() + labs(tag="(c)")

# Lycaenidae
grd_100_p2_obse <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(grd_100, mapping=aes(fill=lycobsRatio), color=NA, alpha=0.7, show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void()+ labs(tag="   ")

# Nymphalidae
grd_100_p3_obse <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(grd_100, mapping=aes(fill=nymobsRatio), color=NA, alpha=0.7, show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void()+ labs(tag="   ")

# Papilionidae
grd_100_p4_obse <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(grd_100, mapping=aes(fill=papobsRatio), color=NA, alpha=0.7, show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void()+ labs(tag="   ")

# Pieridae
grd_100_p5_obse <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(grd_100, mapping=aes(fill=pieobsRatio), color=NA, alpha=0.7, show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void()+ labs(tag="   ")

# tiff("100kmFamilyObse.tiff", units="in", width=7, height=10, res=350)
# grid.arrange(grd_100_p1, grd_100_p2, grd_100_p3, grd_100_p4, grd_100_p5, nrow=5, ncol=1)
# dev.off()

# 100 km resolution specimens only
# Hesperiidae
grd_100_p1_spec <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(grd_100, mapping=aes(fill=hesspecRatio), color=NA, alpha=0.7, show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void() + labs(tag="(b)")

# Lycaenidae
grd_100_p2_spec <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(grd_100, mapping=aes(fill=lycspecRatio), color=NA, alpha=0.7, show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void()+ labs(tag="   ")

# Nymphalidae
grd_100_p3_spec <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(grd_100, mapping=aes(fill=nymspecRatio), color=NA, alpha=0.7, show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void()+ labs(tag="   ")

# Papilionidae
grd_100_p4_spec <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(grd_100, mapping=aes(fill=papspecRatio), color=NA, alpha=0.7, show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void()+ labs(tag="   ")

# Pieridae
grd_100_p5_spec <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(grd_100, mapping=aes(fill=piespecRatio), color=NA, alpha=0.7, show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void()+ labs(tag="   ")

# tiff("100kmFamilySpec.tiff", units="in", width=7, height=10, res=350)
# grid.arrange(grd_100_p1, grd_100_p2, grd_100_p3, grd_100_p4, grd_100_p5, nrow=5, ncol=1)
# dev.off()

# 100 km all records
# Hesperiidae
grd_100_p1 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(grd_100, mapping=aes(fill=hesRatio), color=NA, alpha=0.7, show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void()  + labs(tag="(a)")

# Lycaenidae
grd_100_p2 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(grd_100, mapping=aes(fill=lycRatio), color=NA, alpha=0.7, show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void()+ labs(tag="   ")

# Nymphalidae
grd_100_p3 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(grd_100, mapping=aes(fill=nymRatio), color=NA, alpha=0.7, show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void()+ labs(tag="   ")

# Papilionidae
grd_100_p4 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(grd_100, mapping=aes(fill=papRatio), color=NA, alpha=0.7, show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void()+ labs(tag="   ")

# Pieridae
grd_100_p5 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(grd_100, mapping=aes(fill=pieRatio), color=NA, alpha=0.7, show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
  theme_void()+ labs(tag="   ")

tiff("100kmFamilyAll_withChi.tiff", units="cm", width=19.5, height=20, res=350)
grid.arrange(grd_100_p1, grd_100_p1_spec, grd_100_p1_obse, chi.hes,
             grd_100_p2, grd_100_p2_spec, grd_100_p2_obse, chi.lyc,
             grd_100_p3, grd_100_p3_spec, grd_100_p3_obse, chi.nym,
             grd_100_p4, grd_100_p4_spec, grd_100_p4_obse, chi.pap,
             grd_100_p5, grd_100_p5_spec, grd_100_p5_obse, chi.pie,
             nrow=5, ncol=4)
dev.off()

# # 200 km resolution
# # Hesperiidae
# grd_200_p1 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_200, mapping=aes(fill=hesObsRatio), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Lycaenidae
# grd_200_p2 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_200, mapping=aes(fill=lycObsRatio), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Nymphalidae
# grd_200_p3 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_200, mapping=aes(fill=nymObsRatio), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Papilionidae
# grd_200_p4 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_200, mapping=aes(fill=papObsRatio), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Pieridae
# grd_200_p5 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_200, mapping=aes(fill=pieObsRatio), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# grid.arrange(grd_200_p1, grd_200_p2, grd_200_p3, grd_200_p4, grd_200_p5, nrow=5, ncol=1)
# 
# # 400 km resolution
# # Hesperiidae
# grd_400_p1 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_400, mapping=aes(fill=hesObsRatio), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Lycaenidae
# grd_400_p2 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_400, mapping=aes(fill=lycObsRatio), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Nymphalidae
# grd_400_p3 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_400, mapping=aes(fill=nymObsRatio), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Papilionidae
# grd_400_p4 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_400, mapping=aes(fill=papObsRatio), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Pieridae
# grd_400_p5 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_400, mapping=aes(fill=pieObsRatio), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# grid.arrange(grd_400_p1, grd_400_p2, grd_400_p3, grd_400_p4, grd_400_p5, nrow=5, ncol=1)
# 
# # 800 km resolution
# # Hesperiidae
# grd_800_p1 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_800, mapping=aes(fill=hesObsRatio), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Lycaenidae
# grd_800_p2 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_800, mapping=aes(fill=lycObsRatio), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Nymphalidae
# grd_800_p3 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_800, mapping=aes(fill=nymObsRatio), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Papilionidae
# grd_800_p4 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_800, mapping=aes(fill=papObsRatio), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# # Pieridae
# grd_800_p5 <- ggplot()+
#   geom_sf(data=land, fill="white", color="black")+
#   geom_sf(grd_800, mapping=aes(fill=pieObsRatio), color=NA, alpha=0.7)+
#   scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
#   coord_sf(xlim=c(-3935000, 2957500), ylim=c(-2838545,4563455))+
#   theme_void()
# 
# grid.arrange(grd_800_p1, grd_800_p2, grd_800_p3, grd_800_p4, grd_800_p5, nrow=5, ncol=1)

















################
# BOREAL PLOTS #
################

boreal <- st_read("borealtundra.shp")
boreal <- st_transform(boreal, st_crs(102008))
boreal <- st_crop(boreal, grd_100)

# 100km scale over time #
p1 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(data=boreal, fill="grey", color=NA)+
  geom_sf(grd_100, mapping=aes(fill=t1Ratio, color=(t1Ratio>=0.25)), alpha=0.7,
          show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  scale_color_manual(name="T1 Sufficiently Sampled", values=setNames(c("black", "white"), c(T, F)))+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-260000,4563455))+
  theme_void()+
  theme(legend.title = element_blank())

p2 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(data=boreal, fill="grey", color=NA)+
  geom_sf(grd_100, mapping=aes(fill=t2Ratio, color=(t1Ratio>=0.25 & t2Ratio >=0.25)), alpha=0.7,
          show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  scale_color_manual(name="T1 Sufficiently Sampled", values=setNames(c("black", "white"), c(T, F)))+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-260000,4563455))+
  theme_void()+
  theme(legend.title = element_blank())

p3 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(data=boreal, fill="grey", color=NA)+
  geom_sf(grd_100, mapping=aes(fill=t3Ratio, color=(t1Ratio>=0.25 & t2Ratio >=0.25 & t3Ratio >=0.25)), alpha=0.7,
          show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5,)+
  scale_color_manual(name="T1 Sufficiently Sampled", values=setNames(c("black", "white"), c(T, F)))+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-260000,4563455))+
  theme_void()+
  theme(legend.title = element_blank())

p4 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(data=boreal, fill="grey", color=NA)+
  geom_sf(grd_100, mapping=aes(fill=t4Ratio, color=(t1Ratio>=0.25 & t2Ratio>=0.25 & t4Ratio>=0.25)), alpha=0.7,
          show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  scale_color_manual(name="T1 Sufficiently Sampled", values=setNames(c("black", "white"), c(T, F)))+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-260000,4563455))+
  theme_void()

grid.arrange(p1, p2, p3, p4, nrow=1, ncol=4)

# 200 km scale over time #
p1 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(data=boreal, fill="grey", color=NA)+
  geom_sf(grd_200, mapping=aes(fill=t1Ratio, color=(t1Ratio>=0.25)), alpha=0.7,
          show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  scale_color_manual(name="T1 Sufficiently Sampled", values=setNames(c("black", "white"), c(T, F)))+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-260000,4563455))+
  theme_void()+
  theme(legend.title = element_blank())

p2 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(data=boreal, fill="grey", color=NA)+
  geom_sf(grd_200, mapping=aes(fill=t2Ratio, color=(t1Ratio>=0.25 & t2Ratio >=0.25)), alpha=0.7,
          show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  scale_color_manual(name="T1 Sufficiently Sampled", values=setNames(c("black", "white"), c(T, F)))+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-260000,4563455))+
  theme_void()+
  theme(legend.title = element_blank())

p3 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(data=boreal, fill="grey", color=NA)+
  geom_sf(grd_200, mapping=aes(fill=t3Ratio, color=(t1Ratio>=0.25 & t2Ratio >=0.25 & t3Ratio >=0.25)), alpha=0.7,
          show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5,)+
  scale_color_manual(name="T1 Sufficiently Sampled", values=setNames(c("black", "white"), c(T, F)))+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-260000,4563455))+
  theme_void()+
  theme(legend.title = element_blank())

p4 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(data=boreal, fill="grey", color=NA)+
  geom_sf(grd_200, mapping=aes(fill=t4Ratio, color=(t1Ratio>=0.25 & t2Ratio>=0.25 & t4Ratio>=0.25)), alpha=0.7,
          show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  scale_color_manual(name="T1 Sufficiently Sampled", values=setNames(c("black", "white"), c(T, F)))+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-260000,4563455))+
  theme_void()

grid.arrange(p1, p2, p3, p4, nrow=1, ncol=4)

# 400 km scale #
p1 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(data=boreal, fill="grey", color=NA)+
  geom_sf(grd_400, mapping=aes(fill=t1Ratio, color=(t1Ratio>=0.25)), alpha=0.7,
          show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  scale_color_manual(name="T1 Sufficiently Sampled", values=setNames(c("black", "white"), c(T, F)))+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-260000,4563455))+
  theme_void()+
  theme(legend.title = element_blank())

p2 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(data=boreal, fill="grey", color=NA)+
  geom_sf(grd_400, mapping=aes(fill=t2Ratio, color=(t1Ratio>=0.25 & t2Ratio >=0.25)), alpha=0.7,
          show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  scale_color_manual(name="T1 Sufficiently Sampled", values=setNames(c("black", "white"), c(T, F)))+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-260000,4563455))+
  theme_void()+
  theme(legend.title = element_blank())

p3 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(data=boreal, fill="grey", color=NA)+
  geom_sf(grd_400, mapping=aes(fill=t3Ratio, color=(t1Ratio>=0.25 & t2Ratio >=0.25 & t3Ratio >=0.25)), alpha=0.7,
          show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5,)+
  scale_color_manual(name="T1 Sufficiently Sampled", values=setNames(c("black", "white"), c(T, F)))+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-260000,4563455))+
  theme_void()+
  theme(legend.title = element_blank())

p4 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(data=boreal, fill="grey", color=NA)+
  geom_sf(grd_400, mapping=aes(fill=t4Ratio, color=(t1Ratio>=0.25 & t2Ratio>=0.25 & t4Ratio>=0.25)), alpha=0.7,
          show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  scale_color_manual(name="T1 Sufficiently Sampled", values=setNames(c("black", "white"), c(T, F)))+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-260000,4563455))+
  theme_void()

grid.arrange(p1, p2, p3, p4, nrow=1, ncol=4)

# 800 km scale #
p1 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(data=boreal, fill="grey", color=NA)+
  geom_sf(grd_800, mapping=aes(fill=t1Ratio, color=(t1Ratio>=0.25)), alpha=0.7,
          show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  scale_color_manual(name="T1 Sufficiently Sampled", values=setNames(c("black", "white"), c(T, F)))+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-260000,4563455))+
  theme_void()+
  theme(legend.title = element_blank())

p2 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(data=boreal, fill="grey", color=NA)+
  geom_sf(grd_800, mapping=aes(fill=t2Ratio, color=(t1Ratio>=0.25 & t2Ratio >=0.25)), alpha=0.7,
          show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  scale_color_manual(name="T1 Sufficiently Sampled", values=setNames(c("black", "white"), c(T, F)))+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-260000,4563455))+
  theme_void()+
  theme(legend.title = element_blank())

p3 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(data=boreal, fill="grey", color=NA)+
  geom_sf(grd_800, mapping=aes(fill=t3Ratio, color=(t1Ratio>=0.25 & t2Ratio >=0.25 & t3Ratio >=0.25)), alpha=0.7,
          show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5,)+
  scale_color_manual(name="T1 Sufficiently Sampled", values=setNames(c("black", "white"), c(T, F)))+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-260000,4563455))+
  theme_void()+
  theme(legend.title = element_blank())

p4 <- ggplot()+
  geom_sf(data=land, fill="white", color="black")+
  geom_sf(data=boreal, fill="grey", color=NA)+
  geom_sf(grd_800, mapping=aes(fill=t4Ratio, color=(t1Ratio>=0.25 & t2Ratio>=0.25 & t4Ratio>=0.25)), alpha=0.7,
          show.legend=FALSE)+
  scale_fill_gradient2(low="#3c9ab2", mid="#e8c927", high="#f22300", na.value=NA, limits=c(0,1), midpoint=0.5)+
  scale_color_manual(name="T1 Sufficiently Sampled", values=setNames(c("black", "white"), c(T, F)))+
  coord_sf(xlim=c(-3935000, 2957500), ylim=c(-260000,4563455))+
  theme_void()

grid.arrange(p1, p2, p3, p4, nrow=1, ncol=4)


