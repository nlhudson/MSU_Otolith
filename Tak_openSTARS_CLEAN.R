#make sure data.table is not loaded before running 
#remove.packages("data.table")
#install.packages("data.table", type = "source",
#                 repos = "https://Rdatatable.gitlab.io/data.table")
library(openSTARS)
#library(data.table)i
# give paths to GRASS and where to store the GRASS data base
# Mac e.g. /Library/Frameworks/R.framework/Versions/4.2/Resources/library/openSTARS/extdata/taku_riv THIS IS OLD LOCATION 
# changed because redownload of R: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/openSTARS/extdata/taku_riv NEW CURRENT LOCATION!
grass_program_path <- "/Applications/GRASS-8.2.app/Contents/Resources"
working_dir <- file.path(tempdir(), "taku_workflow")
grass_db_path <- file.path(working_dir, "grassDB")
dir.create(working_dir)
setwd(tempdir())

# OPEN GRASS and set up session with grass_db_path
grass_db_path
# specify the path to the digital elevation model
dem_path <- system.file("extdata", "taku_riv", "tdem.tif", package = "openSTARS")


setup_grass_environment(dem = dem_path, 
                        gisBase = grass_program_path,
                        gisDbase = grass_db_path,
                        location = "t_location",
                        remove_GISRC = TRUE,
                        override = TRUE
)
preds_r_path <- c(system.file("extdata", "taku_riv", "ppt_100m.tif", package = "openSTARS"),
                  system.file("extdata", "taku_riv", "tdem.tif", package = "openSTARS"),
                  system.file("extdata", "taku_riv", "pfi.tif", package = "openSTARS"),
                  system.file("extdata", "taku_riv", "takglac_R.tif", package = "openSTARS"),
                  system.file("extdata", "taku_riv", "slope_tak.tif", package = "openSTARS"),
                  system.file("extdata", "taku_riv", "rf_rm1.tif", package = "openSTARS"),
                  system.file("extdata", "taku_riv", "rf_srsrq1.tif", package = "openSTARS"),
                  system.file("extdata", "taku_riv", "rf_srsrq3.tif", package = "openSTARS"),
                  system.file("extdata", "taku_riv", "rf_agemu.tif", package = "openSTARS"),
                  system.file("extdata", "taku_riv", "xx_glimcode.tif", package = "openSTARS")
)
preds_v_path <- c(system.file("extdata", "taku_riv", "takglac_V.shp", package = "openSTARS"),
                  system.file("extdata", "taku_riv", "glimall.shp", package = "openSTARS"),
                  system.file("extdata", "taku_riv", "geol_glim.shp", package = "openSTARS"),
                  system.file("extdata", "taku_riv", "Grxclass.shp", package = "openSTARS"), #only rock class from geologic maps
                  system.file("extdata", "taku_riv", "Gglimsimp.shp", package = "openSTARS"), #simplified glim classes B, I, A
                  system.file("extdata", "taku_riv", "gs_tats.shp", package = "openSTARS"), #gen stock points
                  system.file("extdata", "taku_riv", "gs_tsms.shp", package = "openSTARS"), #gen stock points
                  system.file("extdata", "taku_riv", "gs_ksl.shp", package = "openSTARS"), #gen stock points
                  system.file("extdata", "taku_riv", "gs_trap.shp", package = "openSTARS"), #gen stock points
                  system.file("extdata", "taku_riv", "gs_kuth.shp", package = "openSTARS"),
                  system.file("extdata", "taku_riv", "gs_tmeni.shp", package = "openSTARS") #gen stock points
                  )

sites_path <- system.file("extdata", "taku_riv", "tak_sites1.shp", package = "openSTARS")
#path that uses original location of HACK01
sites_path <- system.file("extdata", "taku_riv", "tak_sites_hack_og.shp", package = "openSTARS")

#Make df with original HACK01 cords----
hydrositesH <- st_read("/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/openSTARS/extdata/taku_riv/tak_sites1.shp")
      # New latitude and longitude
      new_lat <- 58.189
      new_lon <- -131.651
      # Create an sf object with the new point in WGS84 CRS (EPSG:4326)
      new_point <- st_point(c(new_lon, new_lat))
      new_point_sf <- st_sfc(new_point, crs = 4326)
      # Transform the new point to the same CRS as your original data frame (NAD83 / UTM zone 10N)
      original_crs <- st_crs(hydrositesH)
      new_point_transformed <- st_transform(new_point_sf, original_crs)
      # Find the index of the NAKI03 point in your data frame
      hack01_index <- which(hydrositesH$Site == "HACK01")
# Update the geometry of NAKI03
hydrositesH$geometry[hack01_index] <- new_point_transformed
#plot to check location
library(tmap)
tmap_mode("view")
tm_shape(hydrositesH)+tm_dots()
st_write(hydrositesH, "/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/openSTARS/extdata/taku_riv/tak_sites_hack_og.shp", append = F)
#    sites_path <- system.file("extdata", "taku_riv", "genstock.shp", package = "openSTARS") #gen stock obs point locations



# existing stream network----
#streams_path <- system.file("extdata", "taku_riv", "stream750.shp", package = "openSTARS")
# removed tuls tributary arcid!=451 & arcid!=464 & arcid!=511
streams_path <- system.file("extdata", "taku_riv", "streams_tuls.shp", package = "openSTARS")

# centered pred sites
predsites_path <- system.file("extdata", "taku_riv", "cent_preds.shp", package = "openSTARS")

#import all layer
import_data(dem = dem_path, sites = sites_path, streams = streams_path, pred_sites = predsites_path,
            predictor_vector = preds_v_path, predictor_v_names = c("takglac_V", "glimall", "geol_glim", "Grxclass",
                                                                    "Gglimsimp", 'gs_tats','gs_tsms', 'gs_ksl', 
                                                                   'gs_trap', 'gs_kuth', 'gs_tmeni'), 
            predictor_raster = preds_r_path, predictor_r_names = c("ppt_100m","dem","pfi","takglac_R", "slope_tak", 
                                                                   "rf_rm1","rf_srsrq1","rf_srsrq3",  "rf_agemu", 
                                                                   "xx_glimcodes"))

#getAnywhere('derive_streams')

#FUNCTION to derive.streams() with stream order----

derive.streamsworder <- function (burn = 0, accum_threshold = 700, condition = TRUE, 
                                  min_stream_length = 0, dem_name = NULL, clean = TRUE, mem = FALSE) 
{
  if (condition == TRUE & (ifelse(is.null(dem_name), FALSE, 
                                  dem_name != "dem"))) 
    stop("Only an unmodified DEM should be used for conditioning.")
  if (burn != 0 & (ifelse(is.null(dem_name), FALSE, grepl("burn", 
                                                          dem_name)))) 
    stop("Only an unburnt DEM should be used for burn in.")
  if (is.null(dem_name)) 
    dem_name <- "dem"
  dem_name_out <- dem_name
  vect <- execGRASS("g.list", parameters = list(type = "vect"), 
                    intern = TRUE)
  rast <- execGRASS("g.list", parameters = list(type = "rast"), 
                    intern = TRUE)
  if (!dem_name %in% rast) 
    stop("DEM not found. Did you run import_data()?")
  if (condition) {
    message("Conditioning DEM ...")
    execGRASS("r.hydrodem", flags = c("overwrite"), parameters = list(input = dem_name, 
                                                                      output = "dem_cond"))
    dem_name <- "dem_cond"
    dem_name_out <- "dem_cond"
  }
  if ("streams_o" %in% vect & burn > 0) {
    message("Burning streams into DEM ...")
    execGRASS("v.to.rast", flags = c("overwrite", "quiet"), 
              parameters = list(input = "streams_o", type = "line", 
                                output = "streams_or", use = "val", value = 1))
    dem_name_out <- paste0(dem_name, "_burn", burn)
    if (.Platform$OS.type == "windows") {
      execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), 
                parameters = list(expression = paste0("\"dem2 = if(isnull(streams_or),  ", 
                                                      dem_name, ", ", dem_name, "-", burn, ")\"")))
      execGRASS("g.copy", flags = c("overwrite", "quiet"), 
                parameters = list(raster = paste0("dem2,", dem_name_out)))
      execGRASS("g.remove", flags = c("quiet", "f"), parameters = list(type = "raster", 
                                                                       name = "dem2"))
    }
    else {
      execGRASS("r.mapcalc", flags = c("overwrite", "quiet"), 
                parameters = list(expression = paste0("\"", dem_name_out, 
                                                      " = if(isnull(streams_or),  ", dem_name, ", ", 
                                                      dem_name, "-", burn, ")\"")))
    }
    if (clean) {
      execGRASS("g.remove", flags = c("quiet", "f"), parameters = list(type = "raster", 
                                                                       name = "streams_or"))
    }
  }
  if (mem) {
    fl <- "m"
  }
  else {
    fl <- NULL
  }
  execGRASS("r.watershed", flags = c("overwrite", "quiet", 
                                     fl), parameters = list(elevation = dem_name_out, accumulation = "accums"))
  message("Deriving streams from DEM ...")
  ncell <- execGRASS("g.region", flags = "p", intern = T)
  ncell <- as.numeric(unlist(strsplit(ncell[grep("cells", ncell)], 
                                      split = ":"))[2])
  execGRASS("r.stream.extract", flags = c("overwrite", "quiet"), 
            parameters = list(elevation = dem_name_out, accumulation = "accums", 
                              threshold = accum_threshold, d8cut = ncell, stream_length = min_stream_length, 
                              stream_raster = "streams_r", direction = "dirs"))
  message("Calculating stream topology ...")
  execGRASS("r.stream.order", flags = c("overwrite", "quiet", 
                                        "z", "m"), parameters = list(stream_rast = "streams_r", 
                                                                     direction = "dirs", elevation = dem_name_out, accumulation = "accums", 
                                                                     stream_vect = "streams_v"), ignore.stderr = T)
  execGRASS("g.remove", flags = c("f", "quiet"), type = "raster", 
            name = "streams_r")
  execGRASS("v.db.renamecolumn", flags = "quiet", parameters = list(map = "streams_v", 
                                                                    column = "next_stream,next_str"))
  execGRASS("v.db.renamecolumn", flags = "quiet", parameters = list(map = "streams_v", 
                                                                    column = "flow_accum, flow_accu"))
  execGRASS("v.db.dropcolumn", flags = c("quiet"), parameters = list(map = "streams_v", 
                                                                     columns = c("horton", "shreve", "hack", "topo_dim", 
                                                                                 "scheidegger", "drwal_old", "stright", "sinosoid", 
                                                                                 "source_elev", "outlet_elev", "elev_drop", "out_drop", 
                                                                                 "gradient")))
  a <- execGRASS("v.extract", flags = c("overwrite", "quiet"), 
                 parameters = list(input = "streams_v", output = "streams_v1", 
                                   type = "line", where = paste0("length > 0")), intern = TRUE, 
                 ignore.stderr = TRUE)
  execGRASS("g.copy", flags = c("overwrite", "quiet"), parameters = list(vector = "streams_v1,streams_v"), 
            intern = TRUE, ignore.stderr = TRUE)
  execGRASS("g.remove", flags = c("quiet", "f"), parameters = list(type = "vector", 
                                                                   name = "streams_v1"))
  message("Derived streams saved as 'streams_v'.")
}


# deriver streams including stream orders 
derive.streamsworder(accum_threshold = 500, condition = T, clean = T, burn = 10)
    
    streamso <- readVECT('streams_v',ignore.stderr = TRUE)
    # tmap
    map_streams <- tm_shape(streamso) + 
      tm_lines()
    tmap_mode("view") #or view
    map_streams
#derive_streams() #runs with default arguments ----
#derive_streams(accum_threshold = 500, condition = T, clean = TRUE, burn = 10)


#check complex confluences----
cp <- check_compl_confluences() 

if (cp) +
      correct_compl_confluences()

#optional
#delete lakes #ISSUE breaks lakes and doesnt repair network, because it is not suppost to repair network----
#    lakes_path <- system.file("extdata", "taku_riv", "tlakes.shp", package = "openSTARS")
#    delete_lakes(lakes = lakes_path, keep = F)


# plot -----
streams <- readVECT('streams_v',ignore.stderr = TRUE)
streamsf <- st_as_sf(streams)
str1 <- streamsf %>% select(prev_str01,geometry) 
str1 <- str1[str1$prev_str01 != 0,]
plot(str1)

#NOT needed if use derive.streamworder()#####************* 
#str2 <- streamsf %>% select(prev_str02,geometry) %>% st_drop_geometry()
#str2 <- str2[str2$prev_str02 != 0,] 
#plot(str2)
#stream order index by rid (Used in fish Assignment)
#str2_id <- str2 %>% select(prev_str02) %>% rename(rid=prev_str02) %>% st_drop_geometry()

#streams_with_lakes <- readVECT('streams_v_prev_lakes', ignore.stderr = TRUE)
#plot(streams_with_lakes)
#lakes <- readVECT('lakes', ignore.stderr = TRUE)
#streams$col <- 1
#streams_with_lakes$col <- 1
#lakes$col <- 1
#tmap_mode("plot")
## tmap mode set to plotting
tm_shape(lakes, bbox = streams ) +
  tm_polygons(col = "col", palette = "grey",
              title = "", n = 1, legend.format = list(fun=function(x) "lakes")) +
  tm_shape(streams_with_lakes) +
  tm_lines(col = "col", lwd = 4, legend.col.show = TRUE, palette = "dodgerblue",
           title.col = "", legend.format = list(fun=function(x) "original streams")) +
  tm_shape(streams) +
  tm_lines(col = "col", lwd = 1, lty = 1, palette = "darkblue", legend.col.show = TRUE,
           title.col = "", legend.format = list(fun=function(x) "streams, lakes deleted")) +
  tm_scale_bar() +
  tm_layout(scale = 1, legend.bg.color = T, legend.position = c("left","bottom"))


dem <- read_RAST("dem", ignore.stderr = T)



#be derived for the streams and stored in a new vector map edges.----
calc_edges()

#Prepare sites#----
calc_sites()

# use pred sites genearted for every edge----
calc_sites(predictions = 'cent_preds')
      
# plot-----
      pred_sites <- readVECT("cent_preds", ignore.stderr = T)
      sites <- read_VECT("sites_o", ignore.stderr = TRUE)
      plot(streams)
      points(pred_sites)
      points(sites)

#
#plot -----
library(sp)
streams <- readVECT('streams_v', ignore.stderr = TRUE)
sites <- read_VECT("sites_o", ignore.stderr = TRUE)
ppt <- readRAST('ppt', ignore.stderr = TRUE)

sites$Site


# Prep pred sites----- no need to run if run above with predictions = ""
# calc_prediction_sites(predictions = "preds", dist = 50,  netIDs = 3)



#
# plot-----
pred_sites <- readVECT("cent_preds", ignore.stderr = T)
plot(streams)
points(pred_sites)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### alc edges attributes NOTE
#one giving the attribute for the rca of the edge ("attribute_name_e") and one
#for the attribute of the total catchment of the edge ("attribute_name_c")


#SITES calc edges attributes----
calc_attributes_edges(input_raster = c("dem","dem","dem", 
                                       "slope_tak", "pfi", "ppt_100m", "rf_agemu", 
                                       "rf_rm1","rf_srsrq1", "rf_srsrq3",
                                       "takglac_R"),
                      stat_rast = c("max","min", rep("mean", 3),"sum",c(rep("mean", 4)), "percent"),
                      attr_name_rast = c("maxelev","minelev","avelev", "avslope","avpfi", "acc_ppt","avrfage",
                                         "avrfrm1", "avrfq1", "avrfq3",
                                         "glac"),
                      input_vector = c("geol_glim",  "Gglimsimp", "Grxclass"), 
                      stat_vect = c(rep("percent",3)), 
                      attr_name_vect = c("Glim", "simp_glim", 'rx_class'),
                      round_dig = 6
                      
)

edges <- readVECT("edges", ignore.stderr = T)
head(edges@data, n=5)
names(edges@data)

# OPTIONAL gen stock edge calc #gs_tats, gs_tsms, gs_ksl, gs_kuth, gs_tmeni, gs_tats------------
calc_attributes_edges(input_vector = c('gs_tats','gs_tsms', 'gs_ksl', 'gs_trap', 'gs_kuth', 'gs_tmeni'),
                      stat_vect = c(rep("count",6)),
                      attr_name_vect = c("tats", 'tsms', 'ksl', 'trap', 'kuth', 'tmeni')
)
edges <- readVECT("edges", ignore.stderr = T)
edgesf <- st_as_sf(edges)
print(unique(edgesf$tmenic_c))

#gen stock pred site predictions
calc_attributes_sites_approx(sites_map = "cent_preds", 
                             input_attr_name = c("tats", 'tsms', 'ksl', 'trap', 'kuth', 'tmeni'),
                             output_attr_name = c("tatsC", 'tsmsC', 'kslC', 'trapC', 'kuthC', 'tmeniC'),
                             stat = c(rep("count",6)), 
                             overwrite=TRUE)
preds <- readVECT("cent_preds", ignore.stderr = T)
print(unique(preds@data$tmenC))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#SITES calculate approx. catchment area, e.g., average slope per catchment and percentage of land use classes at each site----
calc_attributes_sites_approx(sites_map = "sites", 
                             input_attr_name = c("maxelev", "minelev",
                                                 "INTRUp","METAMOp", "SEDZp", "UMAFICp", "VOLCp",
                                                 "a_ignp","i_ignp","b_ignp", "metp", "pyrop", "carbp","sedsp",
                                                 "avelev","avslope","avpfi", "acc_ppt", 
                                                 "avrfage","avrfrm1", "avrfq1", "avrfq3",
                                                 "glacp", "mtp", "pap", "pbp", "pip",
                                                 "pyp", "scp", "ssp", "vap", "vbp",
                                                 "vip"),
                             output_attr_name = c("maxelevA", "minelevA",
                                                  "INTRU_P","METAMO_P", "SEDZ_P", "UMAFIC_P", "VOLC_P",
                                                  "a_ignP","i_ignP","b_ignP", "metP", "pyroP", "carbP","sedsP",
                                                  "elevMu","slopeA", "pfiA", "pptAcc", "rf_agemuA", 
                                                  "rf_rm1A","rf_srsrq1A", "rf_srsrq3A",
                                                  "glacP", "mtP", "paP", "pbP", "piP",
                                                  "pyP", "scP", "ssP", "vaP", "vbP",
                                                  "viP"),
                             stat = c(rep("mean",6),rep("mean",7),rep("mean",2), rep("mean", 7), rep("percent",11)), 
                             overwrite=TRUE)

sites <- readVECT("sites", ignore.stderr = TRUE)
head(sites@data, n=26)

#PRED_SITES calculate approx. catchment area, e.g., average slope per catchment and percentage of land use classes at each site----
calc_attributes_sites_approx(sites_map = "cent_preds", 
                             input_attr_name = c("maxelev", "minelev",
                                                 "INTRUp","METAMOp", "SEDZp", "UMAFICp", "VOLCp",
                                                 "a_ignp","i_ignp","b_ignp", "metp", "pyrop", "carbp","sedsp",
                                                 "avelev", "avslope","avpfi", "acc_ppt", 
                                                 "avrfage","avrfrm1", "avrfq1", "avrfq3",
                                                 "glacp", "mtp", "pap", "pbp", "pip",
                                                 "pyp", "scp", "ssp", "vap", "vbp",
                                                 "vip"),
                             output_attr_name = c("maxelevA", "minelevA",
                                                  "INTRU_P","METAMO_P", "SEDZ_P", "UMAFIC_P", "VOLC_P",
                                                  "a_ignP","i_ignP","b_ignP", "metP", "pyroP", "carbP","sedsP",
                                                  "elevMu",
                                                  "slopeA", "pfiA", "pptAcc", "rf_agemuA", 
                                                  "rf_rm1A","rf_srsrq1A", "rf_srsrq3A",
                                                  "glacP", "mtP", "paP", "pbP", "piP",
                                                  "pyP", "scP", "ssP", "vaP", "vbP",
                                                  "viP"),
                             stat = c(rep("mean",6),rep("mean",7),rep("mean",2), rep("mean", 7), rep("percent",11)), 
                             overwrite=TRUE)
preds <- readVECT("cent_preds", ignore.stderr = TRUE)
head(preds@data, n=5)

#OPTION: write all files to an wd ssn folder----
simpath <- "~/bigdata/ssn/takuriv_ssn_file"
simp_ssn_dir <- file.path(simpath, 'tak102822.ssn')
export_ssn(simp_ssn_dir, predictions="preds", delete_directory = T)
list.files(simp_ssn_dir)


#
#best OPTION: write all files to local directory ssn folder----
ssnlocpath <- "/Users/kylebrennan/Documents/ssn_object_file"
locssn_dir <- file.path(ssnlocpath, 'takssn_hack01_ogloc.ssn') #change DATE on ssn
export_ssn(locssn_dir, predictions = 'cent_preds', delete_directory = TRUE)
#
#OPTION: write to a temp directory----
ssn_dir <- file.path(tempdir(), 'tak102822.ssn')
export_ssn(ssn_dir, predictions = 'preds', delete_directory = T)
list.files(ssn_dir)
