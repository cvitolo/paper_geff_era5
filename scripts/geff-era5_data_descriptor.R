# INSTALL AND LOAD PACKAGES ####################################################
packs <- c("cffdrs", "leaflet", "ggplot2", "dplyr", "raster", "rgeos",
           "sp", "lubridate", "httr", "rgdal", "ggmap", "gridExtra", "viridis",
           "lutz", "xtable")
new.packages <- packs[!(packs %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packs, require, character.only = TRUE)
rm(packs, new.packages)

#### ALGORITHM VALIDATION ######################################################

# Compare CFFDRS and GEFF algorithms
y <- structure(list(long = -91.82, lat = 45.98, hour = 18,
                    yr = 2017L, mon = 1L, day = 1L, temp = 17,
                    rh = 42, ws = 25, prec = 0),
               .Names = c("long", "lat", "hour", "yr", "mon", "day", "temp",
                          "rh", "ws", "prec"),
               row.names = 1L, class = "data.frame")
init <- data.frame(ffmc = 85, dmc = 6, dc = 15, lat = y$lat)
round(cffdrs::fwi(input = y, init = init, out = "fwi"), 2)

# ALGORITHM        & FFMC  & DMC  & DC    & ISI   & BUI  & FWI  & DSR\\
# R package cffdrs & 87.69 & 7.29 & 17.76 & 10.85 & 7.28 & 9.46 & 1.45\\

# GEFF - manual run of GEFF (Fortran code) from the above initial conditions
# Here are the results:
# ALGORITHM  & FFMC  & DMC  & DC    & ISI   & BUI  & FWI   & DSR\\
# GEFF-ERAI & 87.70 & 8.54 & 19.01 & 10.80 & 8.49 & 10.10 & 1.63\\
# GEFF-ERA5 & 87.69 & 8.54 & 19.01 & 10.85 & 8.49 & 10.10 & 1.63\\

rm(list = ls())

# GET DATA FROM GFWED for 2017 #################################################
myDates <- seq.Date(from = as.Date("2017-01-01"), to = as.Date("2017-12-31"),
                    by = "day")
for (myday in seq_along(myDates)){
  print(myday)
  tmpdir <- "data/GFWED"
  reformatted_date <- gsub(pattern = "-", replacement = "", myDates[myday])
  tmpfilename <- paste0("FWI.GEOS-5.Daily.Default.", reformatted_date, ".nc")
  download.file(url = paste0("https://portal.nccs.nasa.gov/datashare/GlobalFWI/",
                              "v2.0/fwiCalcs.GEOS-5/Default/GEOS-5/2017/",
                              tmpfilename),
                destfile = file.path(tmpdir, tmpfilename))
}

# EXTRACT YEAR 2017 ONLY FROM ERAI AND ERA5 ####################################

# FWI based on ERAI
# Get fwi.nc (v3.0) from Zenodo: https://doi.org/10.5281/zenodo.3251000
# Rename the file "fwi_erai.nc", then load it as a brick
erai <- brick("/hugetmp/Downloads/fwi_erai.nc")
# Get the indices of the year 2017
idx <- which(substr(names(erai), 2, 5) == 2017)
# Subset the datacube (2017 only)
erai <- erai[[idx]]
# Save it as "data/fwi2017erai.nc"
writeRaster(erai, filename = "data/fwi2017erai.nc", format = "CDF",
            overwrite = TRUE)

rm(list = ls())

# FWI based on ERA5
# Fire danger indices are available from the Copernicus Climate Data Store,
# however the index FWI is also available from Zenodo.
# Get data from CDS or Zenodo: https://doi.org/10.5281/zenodo.3269270
# Rename the file "fwi_era5.nc", then load it as a brick
era5 <- brick("/hugetmp/Downloads/fwi_era5.nc")
# Get the indices of the year 2017
idx <- which(substr(names(era5), 2, 5) == 2017)
# Subset the datacube (2017 only)
era5 <- era5[[idx]]
# Save it as "data/fwi2017erai.nc"
writeRaster(era5, filename = "data/fwi2017era5.nc", format = "CDF",
            overwrite = TRUE)

rm(list = ls())

# GET DATA FROM SYNOP STATIONS #################################################
# The code in this section can only run within ECMWF internal network because it
# uses internal resources and tools (stvl)

# Set suitable dates
myDates <- seq.Date(from = as.Date("2017-01-01"), to = as.Date("2017-12-31"),
                    by = "day")
myParams <- c("tp", "2t", "10dd", "10ff", "2d")
log_error <- NULL
# Loop over the hours
for (myHour in seq(0, 23, 1)){
  df <- matrix(NA, nrow = 0, ncol = 8)
  # Loop over the parameters
  for (j in seq_along(myParams)){
    myParam <- myParams[j]
    for (i in seq_along(myDates)){
      print(paste0(myHour, "UTC ", myParam, " ", myDates[[i]]))

      # Define request parameters
      params <- list(table = "observation",
                     dataset = "synop",
                     year = year(myDates[[i]]),
                     month = month(myDates[[i]]),
                     day = day(myDates[[i]]),
                     hour = myHour,
                     parameter = myParam,
                     get_data = "ascii")

      # Total precipitation should be accumulated over period = 24h*3600s
      if (myParam == "tp"){params[["period"]] <- 86400}

      # Make POST request
      y <- httr::POST("http://ecverify.ecmwf.int/cgi-bin/stvl/stvl_browse.py",
                      body = params)

      # Extract table (if available)
      x <- try(expr = read.table(text = content(y, "text", encoding = "utf-8"),
                                 sep = " "),
               silent = TRUE)
      if (class(x) != "try-error"){
        x$time <- myDates[[i]]
        x$param <- myParam
        x$hour <- myHour
        df <- rbind(df, x)
      }else{
        print("Data unavailable")
        log_error <- c(log_error,
                       paste0(myHour, "UTC ", myParam, " ", myDates[[i]]))
      }
    }
  }
  saveRDS(df, paste0("data/synop/tempdf_", myHour, ".rds"))
}

# Check log_error before proceeding.

rm(list = ls())

# Add local time to the attributes, using lutz::tz_lookup_coords().

for (myHour in 0:23){

  print(paste(myHour, "UTC"))

  # Read data.frame
  df <- readRDS(paste0("data/synop/tempdf_", myHour, ".rds"))
  names(df) <- c("lat", "long", "id", "z", "value", "date", "param", "UTC_time")
  # Remove the column with z coordinate and re-arrange
  df <- df[, c(3, 1, 2, 6, 8, 7, 5)]

  # Collate info on synop stations
  temp <- unique(df[, c("id", "lat", "long")])
  if (exists("stations_all")){
    stations_all <- rbind(stations_all, temp)
    stations_all <- unique(stations_all)
  }else{
    stations_all <- temp
  }

  # Convert from long to wide format
  spreadDF <- tidyr::spread(data = df, param, value); rm(df, temp)

  if (all(c("10dd", "10ff", "2d", "2t", "tp") %in% names(spreadDF))) {
    # Remove rows with incomplete records
    df <- spreadDF %>% filter(complete.cases(`10dd`, `10ff`, `2d`, `2t`, `tp`))
    rm(spreadDF)

    # Add time zone
    df$tzid <- lutz::tz_lookup_coords(lat = df$lat, lon = df$long,
                                      method = "fast") # accurate

    # Remove records with NA tzid or timezone = "uninhabited"
    temp <- which(is.na(df$tzid) | df$tzid == "uninhabited")
    if (length(temp) > 0) {
      df <- df[-temp, ]
    }

    # UTC time as needed by with_tz()
    df$timestamp_utc <- as.POSIXct(paste0(as.character(df$date), " ",
                                          myHour, ":00"),
                                   format = "%Y-%m-%d %H:%M", tz = "UTC")

    df <- df %>% group_by(tzid) %>%
      mutate(local_time = hour(with_tz(timestamp_utc, tzid)))

    # Remove records with NA local time
    temp <- which(is.na(df$local_time))
    if (length(temp) > 0) {
      df <- df[-temp, ]
    }

    # Keep only records corresponding to 11-13
    # This is to capture noon local time even with Daylight Saving Time
    temp <- which(df$local_time %in% 11:13)
    if (length(temp) > 0) {
      df <- df[temp, ]

      # Calculate relative humidity
      df$`2d` <- df$`2d` - 273.15
      df$`2t` <- df$`2t` - 273.15
      df$rh <- (6.11 * 10 ^ (7.5 * (df$`2d`) / (237.7 + (df$`2d`)))) /
        (6.11 * 10 ^ (7.5 * (df$`2t`) / (237.7 + (df$`2t`)))) * 100

      # Convert to match cffdrs::fwi() requirements
      df$ws <- df$`10ff` * 3.6
      df$yr <- lubridate::year(df$date)
      df$mon <- lubridate::month(df$date)
      df$day <- lubridate::day(df$date)

      # Add info about hemisphere and season
      df$hemisphere <- ifelse(df$lat >= 0, "North", "South")
      df$season <- NA
      # North Hemisphere
      df$season[which(df$lat >= 0 & between(df$mon, 4, 9))] <- "Dry"
      df$season[which(df$lat >= 0 &
                        (between(df$mon, 1, 3) |
                           between(df$mon, 10, 12)))] <- "Wet"
      # South Hemisphere
      df$season[which(df$lat < 0 & between(df$mon, 4, 9))] <- "Wet"
      df$season[which(df$lat < 0 &
                        (between(df$mon, 1, 3) |
                           between(df$mon, 10, 12)))] <- "Dry"

      if (exists("df_all")){
        df_all <- rbind(df_all, df)
      }else{
        df_all <- df
      }

      print(paste("records:", dim(df_all)[1],
                  "- stations:", dim(stations_all)[1]))

    } else {
      message("No data at local noon")
    }

  }

}

# AVERAGE OVER THE 11-13 TIME WINDOW
df <- df_all %>%
  group_by(id, lat, long, tzid, season, yr, mon, day) %>%
  summarise(temp = last(`2t`), prec = last(tp), rh = last(rh), ws = last(ws))

# Save all the stations as they are
saveRDS(df, "data/df_all.rds")
# This dataset is saved in the data folder of the repo as it can only be obtained using tools internal to ECMWF.

rm(list = ls())

# PUT TOGETHER DATA FROM SYNOP STATIONS, ERAI AND ERA5 #########################

erai <- raster::brick("data/fwi2017erai.nc")
# The dataset below needs to be rotated because the original longitude range is
# [0,360] while erai is in [-180, +180]
era5 <- raster::rotate(raster::brick("data/fwi2017era5.nc"))

# If the layers have no dates, we add them.
names(erai) <- names(era5) <- seq.Date(from = as.Date("2017-01-01"),
                                       to = as.Date("2017-12-31"),
                                       by = "day")

# Load unique stations and data
df <- readRDS("data/df_all.rds")
# There are stations with same id but different lat/lon, let's separate them
station_unique <- unique(df[, c("id", "lat", "long")])

df_used <- data.frame(matrix(NA, nrow = 0, ncol = ncol(df) + 3))
names(df_used) <- c(names(df), "OBS", "ERAI", "ERA5")
for (i in seq_along(station_unique$id)){

  # We filter over the station
  dfx <- df %>% filter(id == station_unique$id[i],
                       lat == station_unique$lat[i],
                       long == station_unique$long[i])

  dim_dfx <- dim(dfx)[1]

  # Discard stations with less than 30 days recording
  if (dim_dfx >= 30){

    pt <- data.frame(long = dfx$long[1], lat = dfx$lat[1])
    initial_condition <- data.frame(ffmc = 85, dmc = 6, dc = 15,
                                    lat = dfx$lat[1])
    dfx$OBS <- cffdrs::fwi(input = dfx,
                              init = initial_condition,
                              out = "fwi")$FWI

    # Extract the modelled FWI from ERAI
    temp <- t(raster::extract(x = erai, y = SpatialPoints(pt)))
    mytimestamps <- which(substr(names(erai), 2, 11) %in%
                            paste0(dfx$yr, ".", sprintf("%02d", dfx$mon), ".",
                                   sprintf("%02d", dfx$day)))
    dfx$ERAI <- temp[mytimestamps]

    # Extract the modelled FWI from ERA5
    temp <- t(raster::extract(x = era5, y = SpatialPoints(pt)))
    mytimestamps <- which(substr(names(era5), 2, 11) %in%
                            paste0(dfx$yr, ".", sprintf("%02d", dfx$mon), ".",
                                   sprintf("%02d", dfx$day)))
    dfx$ERA5 <- temp[mytimestamps]

    df_used <- dplyr::bind_rows(df_used, dfx)
    rm(mytimestamps, temp, initial_condition, dfx)
  }else{
    print(paste("Station n.", i, "discarded because contains < 30 records"))
  }

}

saveRDS(df_used, "data/df_geff_erai_era5.rds")
rm(list = ls())

#### Data validation: comparison with observed FWI and old reanalysis (ERAI) ###

df <- readRDS("data/df_geff_erai_era5.rds")

# Split tzid into region and subregion
dfx <- df %>%
  filter(OBS <= 250) %>% # Remove oddly high value in observations
  na.omit %>%
  mutate(region = sapply(strsplit(tzid, "/"), `[`, 1),
         subregion = sapply(strsplit(tzid, "/"), `[`, 2)) %>%
  filter(region != "Etc") %>% # remove undefined zones
  # dplyr::select(id, lat, long, region, season, OBS, ERAI, ERA5) %>%
  dplyr::select(id, lat, long, region, OBS, ERAI, ERA5) %>%
  # group_by(id, lat, long, region, season) %>%
  group_by(id, lat, long, region) %>%
  #add_tally() %>% # add station count
  #filter(n >= 30) %>% # remove stations with less than 30 records
  summarise( # Bias and Anomaly correlation by time, id and season
            bias_erai = mean(OBS - ERAI, na.rm = TRUE),
            bias_era5 = mean(OBS - ERA5, na.rm = TRUE),
            anomaly_correlation_erai = cor(OBS - mean(OBS, na.rm = T),
                                           ERAI - mean(ERAI, na.rm = T),
                                           use = "complete.obs"),
            anomaly_correlation_era5 = cor(OBS - mean(OBS, na.rm = T),
                                           ERA5 - mean(ERA5, na.rm = T),
                                           use = "complete.obs")) %>%
  na.omit

# As the standard deviation for ERAI generates NAs, some records are removed and
# this affects the skill of ERA5 as well

# Summarise results in a table
dfx_summary <- dfx %>%
  # group_by(region, season) %>%
  group_by(region) %>%
  summarise(count = n(),
            bias_erai = round(median(bias_erai, na.rm = T), 2),
            bias_era5 = round(median(bias_era5, na.rm = T), 2),
            ac_erai = round(median(anomaly_correlation_erai, na.rm = T), 2),
            ac_era5 = round(median(anomaly_correlation_era5, na.rm = T), 2))
# copy-paste this into latex main.tex
xtable::xtable(x = dfx_summary, caption = "Validation")

range(dfx_summary$bias_erai)
range(dfx_summary$bias_era5)

# Is ERA5 bias statistically significantly different from ERAI bias? NO!
t.test(dfx_summary$bias_erai, dfx_summary$bias_era5, alternative = "two.sided",
       var.equal = TRUE)

stations_used <- dfx %>% group_by(id) %>%
  summarize(lat = as.numeric(names(which.max(table(lat)))),
            long = as.numeric(names(which.max(table(long))))) # 3077 stations

# Check locations on interactive map
leaflet(data = stations_used) %>%
  addTiles() %>%
  addMarkers(~long, ~lat, label = ~as.character(id))

# How many stations are in the North hemisphere?
round(prop.table(table(stations_used$lat >= 0)), 2)

############################# FIGURE 1 #########################################

# Screenshot of the CDS web interface

############################# FIGURE 2 #########################################

rm(list = ls())

# REANALYSIS 2017 only
fwi2017 <- raster::rotate(raster::brick("data/fwi2017era5.nc")[[1]])

# These are all the synop stations
df <- readRDS("data/df_geff_erai_era5.rds")

# Generate figure

# Convert to sp objects
sp::coordinates(dfx) <- ~long+lat

# Use GEFF-RE grid to plot the world
world <- fwi2017[[1]]
world[world > 0] <- 0

pdf(file = "Synops.pdf", width = 10, height = 6.7)
raster::plot(world, col = "gray95", legend = FALSE)
raster::plot(dfx, col = "#2A69A2", pch = 19, cex = 0.1, add = TRUE)
legend(x = "bottom", legend = "Stations used for validation",
       col = "#2A69A2", pch = 20, bty = "n", horiz = TRUE)
dev.off()

### This figure was prettified in QGIS

############################# FIGURE 3 #########################################

dfx <- df %>%
  filter(OBS <= 250) %>% # Remove oddly high value in observations
  # filter(OBS > 0) %>%
  na.omit %>%
  mutate(region = sapply(strsplit(tzid, "/"), `[`, 1),
         subregion = sapply(strsplit(tzid, "/"), `[`, 2)) %>%
  filter(region != "Etc") # remove undefined zones

# Explore global distributions
x <- reshape2::melt(data = dfx[, c("region", "OBS", "ERAI", "ERA5")],
                    id = "region")

# Boxplots by regions
ggplot(x, aes(x=variable, y=value, fill = region)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~region, scales = "free_y") +
  xlab("") + ylab("FWI") + theme_bw() +
  theme(text = element_text(size=20)) +
  coord_cartesian(ylim = c(0, 150)) +
  scale_fill_discrete(name = "Region", h = c(0, 360), c = 80, l = 60)

# Custom palette:
#E16A86 Africa
#C7821C America
#909800 Arctic
#00A846 Asia
#00AD9A Atlantic
#00A2D3 Australia
#9183E6 Europe
#D766C9 Indian
#E16A86 Pacific

rm(list = ls())

### Comparison with ENSO #######################################################

# Crop reanalysis over SE Asia
system("cdo sellonlatbox,90,132,-14,21 /scratch/rd/nen/perClaudia/era5/fwi.nc data/fwi_era5_seasia.nc")
r <- brick("data/fwi_era5_seasia.nc")
seasia_clima_98 <- caliver::daily_clima(r, probs = 0.98)
# writeRaster(seasia_clima_98, filename = "data/clima_98_seasia.nc",
#             format = "CDF", overwrite = TRUE)
days_above_98 <- stack()
for (myyear in 1980:2018){ # myyear <- 1980
  print(myyear)
  if (lubridate::leap_year(myyear)){
    day_indices <- 1:366
  }else{
    day_indices <- c(1:59, 61:366) # skip 29th Feb
  }
  idx <- which(substr(names(r), 2, 5) == myyear)
  r_year <- r[[idx]]
  days_above_98thp <- calc(r_year > seasia_clima_98[[day_indices]], sum)
  days_above_98 <- stack(days_above_98, days_above_98thp)
}
# writeRaster(days_above_98, filename = "data/days_above_98_seasia.nc",
#             format = "CDF", overwrite = TRUE)

p98 <- na.omit(setNames(as.data.frame(days_above_98), 1980:2018)); p98$P <- 98
df <- reshape2::melt(p98, id.vars = "P")
df$P <- as.factor(df$P)
df$col <- ifelse(df$variable %in% c(1982, 1983, 1997, 1998, 2015, 2016,
                                    1987, 1988, 1991, 1992), "red",
                 ifelse(df$variable %in% c(1988, 1989, 1999, 2000, 2007, 2008,
                                           2010, 2011), "blue", "gray"))
df$alph <- ifelse(df$variable %in% c(1982, 1983, 1997, 1998, 2015, 2016, 1988,
                                     1989, 1999, 2000, 2007, 2008, 2010, 2011),
                  0.5, 0.1)
# https://ggweather.com/enso/oni.htm
vstrong_nino <- c(1982, 1983, 1997, 1998, 2015, 2016)
strong_nino <- c(1987, 1988, 1991, 1992)
mod_nino <- c(1986, 1994, 2002)
weak_nino <- c(2004, 2014, 2018)
strong_nina <- c(1989, 1999, 2000, 2007, 2008, 2010, 2011)
mod_nina <- c(1995, 1996, 2012)
weak_nina <- c(1984, 1985, 2001, 2005, 2006, 2009)
df$Intensities <- ifelse(df$variable %in% vstrong_nino, "Very Strong El Nino",
                         ifelse(df$variable %in% strong_nino,
                                "Strong El Nino",
                                ifelse(df$variable %in% mod_nino,
                                       "Moderate El Nino",
                                       ifelse(df$variable %in% strong_nina,
                                              "Strong La Nina",
                                              ifelse(df$variable %in% mod_nina,
                                                     "Moderate La Nina",
                                                     "Weak (either)")))))
df$Intensities <- as.factor(x = df$Intensities)
levels(df$Intensities) <- c("Very Strong El Nino", "Strong El Nino",
                            "Moderate El Nino", "Strong La Nina",
                            "Moderate La Nina", "Weak (either)")

# ggplot(df, aes(x = value, y = P, color = variable, alpha = alph), size = 0.1) +
#   geom_point() + geom_jitter() +
#   scale_color_manual(name = "Year",
#                      values=c(rep("lightgray", 2), rep("red", 2),
#                               rep("lightgray", 4), rep("blue", 2),
#                               rep("lightgray", 7), rep("red", 2),
#                               rep("blue", 2), rep("lightgray", 6),
#                               rep("blue", 2), rep("lightgray", 1),
#                               rep("blue", 2), rep("lightgray", 3),
#                               rep("red", 2), rep("lightgray", 2))) +
#   scale_alpha(guide = "none") + xlab("") + ylab("Percentile") + theme_bw() +
#   ggtitle("SE Asia")
#
# ggplot(df, aes(y = value, x = P, fill = Intensities)) +
#   geom_boxplot() + xlab("Percentile") + ylab("") + theme_bw() +
#   ggtitle("SE Asia") + coord_flip()

############################# FIGURE 4 #########################################

ggplot(df, aes(x=variable, y=value, fill = col)) +
  geom_boxplot(outlier.alpha = 0) +
  theme(text = element_text(size=15),
        axis.text.x = element_text(angle = 90)) +
  xlab("") + ylab("Number of days above threshold") + ylim(0, 60) +
  scale_fill_manual(name = "",
                    labels = c("La Nina", "Either", "El Nino"),
                    breaks = c("blue", "gray", "red"),
                    values = c("blue", "gray", "red"))

############################# FIGURE 5 #########################################

# Map of days above threshold
myMap <- get_stamenmap(bbox = c(left = bbox(days_above_98)[[1]],
                                bottom = bbox(days_above_98)[[2]],
                                right = bbox(days_above_98)[[3]],
                                top = bbox(days_above_98)[[4]]),
                       maptype="toner-lite", color="bw", zoom=5, crop = T)

myplot <- ggmap(myMap) + xlab("Longitude") + ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))

# For year 1997
rtp <- rasterToPolygons(days_above_98[[18]])
names(rtp) <- "layer"
rtp$layer <- cut(rtp$layer, seq(0,100,10), include.lowest=T)

myplot + ggtitle("1997") +
  geom_polygon(data = rtp,
               aes(x = long, y = lat, group = group,
                   fill = rep(rtp$layer, each = 5)),
               size = 0, alpha = 0.5) +
  scale_fill_manual(name = "Number of days\nabove threshold",
                    values = c("gray20", "gray40", "gray60", "gray80", "yellow",
                               "green", "cyan", "blue", "red", "darkred"),
                    drop = FALSE)

# For year 2010
rtp <- rasterToPolygons(days_above_98[[31]])
names(rtp) <- "layer"
rtp$layer <- cut(rtp$layer, seq(0,100,10), include.lowest=T)

myplot + ggtitle("2010") +
  geom_polygon(data = rtp,
               aes(x = long, y = lat, group = group,
                   fill = rep(rtp$layer, each = 5)),
               size = 0, alpha = 0.5) +
  scale_fill_manual(name = "Number of days\nabove threshold",
                    values = c("gray20", "gray40", "gray60", "gray80", "yellow",
                               "green", "cyan", "blue", "red", "darkred"),
                    drop = FALSE)
