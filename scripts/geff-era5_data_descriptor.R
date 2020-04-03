# LOAD PACKAGES ################################################################

library("raster")
library("dplyr")
library("leaflet")
library("ggplot2")
library("httr")
library("lubridate")
library("ggmap")
library("colorspace")
library("htmlwidgets")

# GET DATA FROM GFWED, ERAI AND ERA5 FOR 2017 ##################################

# FWI based on GFWED
myDates <- seq.Date(from = as.Date("2017-01-01"),
                    to = as.Date("2017-12-31"),
                    by = "day")

for (myday in seq_along(myDates)){
  print(myday)
  reformatted_date <- gsub(pattern = "-", replacement = "", myDates[myday])
  tmpfilename <- paste0("FWI.GEOS-5.Daily.Default.", reformatted_date, ".nc")
  download.file(url = paste0("https://portal.nccs.nasa.gov/datashare/",
                             "GlobalFWI/v2.0/fwiCalcs.GEOS-5/Default/",
                             "GEOS-5/2017/",
                             tmpfilename),
                destfile = file.path("/hugetmp/reanalysis/GFWED/all",
                                     tmpfilename))
}
system("cdo select,name=GEOS-5_FWI /hugetmp/reanalysis/GFWED/all/* data/fwi2017gfwed.nc")
system("cdo remapnn,data/fwi2017era5.nc data/fwi2017gfwed.nc data/fwi2017gfwed_remapped.nc")

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
# Fire danger indices are available from the Copernicus Climate Data Store
system("cdo cat /hugetmp/reanalysis/GEFF-ERA5/hres/fwi/ECMWF_FWI_2017* data/fwi2017era5_temp.nc")
# The dataset needs to be rotated because the original longitude range is
# [0,360] while erai and gfwed are in [-180, +180]
system("cdo sellonlatbox,-180,+180,-90,+90 data/fwi2017era5_temp.nc data/fwi2017era5.nc")
rm(list = ls())

# Get observations from SYNOP stations using STVL and R (no longer used!) ######

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
      x <- try(expr = read.table(text = httr::content(y, "text", encoding = "utf-8"),
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

df <- df_all %>%
  group_by(id, lat, long, tzid, season, yr, mon, day) %>%
  summarise(temp = last(`2t`), prec = last(tp), rh = last(rh), ws = last(ws))

# Save all the stations as they are
saveRDS(df, "data/df_all_stvl_R.rds")

# Get observations from SYNOP stations using STVL and Python3 ##################

# The following Python3 script can only run within ECMWF internal network
# because it uses internal resources and tools (STVL)
library("reticulate")
reticulate::py_run_file("scripts/STVL_get_observations.py")

# Add local time to the attributes, using lutz::tz_lookup_coords()

for (i in 1:12){

  print(i)

  # Read data.frame
  df <- read.csv(paste0("/hugetmp/SYNOP/Observations_2017",
                        sprintf("%02d", i), ".csv"))
  # Remove unecessary columns
  df <- df[, -which(names(df) %in% c("X", "level", "elevation"))]

  # Convert from long to wide format
  spreadDF <- df %>% tidyr::pivot_wider(id_cols = names(df)[1:4],
                                        names_from = param,
                                        values_from = value_0)

  if (all(c("10dd", "10ff", "2d", "2t", "tp") %in% names(spreadDF))) {
    # Remove rows with incomplete records
    df <- spreadDF %>% filter(complete.cases(`10dd`, `10ff`, `2d`, `2t`, `tp`))
    rm(spreadDF)

    # Add time zone
    df$tzid <- lutz::tz_lookup_coords(lat = df$latitude,
                                      lon = df$longitude,
                                      method = "accurate") # or method = "fast"

    # Remove records with NA tzid or timezone = "uninhabited"
    temp <- which(is.na(df$tzid) | df$tzid == "uninhabited")
    if (length(temp) > 0) {
      df <- df[-temp, ]
    }

    # UTC time as needed by with_tz()
    df$timestamp_utc <- as.POSIXct(df$date,
                                   format = "%Y-%m-%d %H:%M", tz = "UTC")

    df <- df %>%
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
      df$rh <- caliver:::relative_humidity(t2m = df$`2t`, d2m = df$`2d`)

      # Convert to match cffdrs::fwi() requirements
      df$`2d` <- df$`2d` - 273.15
      df$`2t` <- df$`2t` - 273.15
      df$ws <- df$`10ff` * 3.6
      df$yr <- lubridate::year(df$date)
      df$mon <- lubridate::month(df$date)
      df$day <- lubridate::day(df$date)

      names(df)[1:3] <- c("id", "lat", "long")

      if (exists("df_all")){
        df_all <- rbind(df_all, df)
      }else{
        df_all <- df
      }

      print(paste("records:", dim(df_all)[1]))

    } else {
      message("No data for this month")
    }

  }

}

# # Average over the 11-13 time window
df <- df_all %>%
  group_by(id, lat, long, tzid, yr, mon, day) %>%
  summarise(temp = last(`2t`), prec = last(tp), rh = last(rh), ws = last(ws))

# Cleanup

summary(df$lat)
summary(df$long)

summary(df$temp) # Temperatures exceed the highest value ever recorded!
plot(density(df$temp))
df$temp[which(df$temp > 56.7)] <- NA # Remove
df <- df[!is.na(df$temp), ]
plot(density(df$temp))

summary(df$prec)
plot(density(df$prec))

summary(df$rh) # RH exceed 100%!
plot(density(df$rh))
df$rh[which(df$rh > 100.0001)] <- NA
df <- df[!is.na(df$rh), ]
plot(density(df$rh))

summary(df$ws) # WS exceed strongest wind ever recorded ~ 115m/s (no tornado)
plot(density(df$ws))
df$ws[which(df$ws > 115)] <- NA
df <- df[!is.na(df$ws), ]
plot(density(df$ws))

# Save all the stations as they are
saveRDS(df, "data/df_all_stvl_python_clean.rds")
rm(list = ls())

# PUT TOGETHER DATA FROM SYNOP STATIONS, ERAI, ERA5 AND GFWED ##################

erai <- raster::brick("data/fwi2017erai.nc")
gfwed <- raster::brick("data/fwi2017gfwed_remapped.nc")
era5 <- raster::brick("data/fwi2017era5.nc")

# If the layers have no dates, we add them.
names(erai) <- names(era5) <- names(gfwed) <-
  seq.Date(from = as.Date("2017-01-01"), to = as.Date("2017-12-31"), by = "day")

# Load unique stations and data
df <- readRDS("data/df_all_stvl_python_clean.rds")
ids <- unique(df$id)

df_used <- data.frame(matrix(NA, nrow = 0, ncol = ncol(df) + 4))
names(df_used) <- c(names(df), "OBS", "ERAI", "GFWED", "ERA5")

for (i in seq_along(ids)){

  print(i)

  # We filter over the station id
  dfx <- df %>% filter(id == ids[i]) %>% select(id, lat, long, tzid,
                                                yr, mon, day,
                                                temp, prec, rh, ws)

  # Discard stations with less than 30 days recording
  if (dim(dfx)[1] >= 30){

    initial_condition <- data.frame(ffmc = 85, dmc = 6, dc = 15,
                                    lat = dfx$lat[1])

    dfx$OBS <- cffdrs::fwi(input = dfx,
                           init = initial_condition,
                           out = "fwi")$FWI

    # Generate spatial points
    pt <- SpatialPoints(data.frame(long = dfx$long[1], lat = dfx$lat[1]))

    # Extract the modelled FWI from ERAI
    temp <- t(raster::extract(x = erai, y = pt))
    mytimestamps <- which(substr(names(erai), 2, 11) %in%
                            paste0(dfx$yr, ".", sprintf("%02d", dfx$mon), ".",
                                   sprintf("%02d", dfx$day)))
    dfx$ERAI <- temp[mytimestamps]

    # Extract the modelled FWI from GFWED
    temp <- t(raster::extract(x = gfwed, y = pt))
    dfx$GFWED <- temp[mytimestamps]

    # Extract the modelled FWI from ERA5
    temp <- t(raster::extract(x = era5, y = pt))
    dfx$ERA5 <- temp[mytimestamps]

    df_used <- dplyr::bind_rows(df_used, dfx)
    rm(mytimestamps, temp, initial_condition, dfx)
  }else{
    print(paste("Station n.", i, "discarded because contains < 30 records"))
  }

}

df_used$GFWED[df_used$GFWED == "NaN"] <- NA

# Clean
summary(df_used$OBS) # Obs FWI is not expected to exceed 200!
summary(df_used$ERAI)
summary(df_used$GFWED)
summary(df_used$ERA5)

plot(density(df_used$OBS))
df_used$OBS[which(df_used$OBS > 200)] <- NA # Remove
df_used <- df_used[!is.na(df_used$OBS), ]
plot(density(df_used$OBS))

plot(density(na.omit(df_used$ERAI)))
plot(density(na.omit(df_used$GFWED)))
plot(density(na.omit(df_used$ERA5)))

# Save
saveRDS(df_used, "data/df_geff_erai_gfwed_era5_stvl_python_clean.rds")
rm(list = ls())

############################# FIGURE 1 #########################################

# Screenshot of the CDS web interface

############################# FIGURE 2 #########################################

# Data validation: comparison with observed FWI, ERAI and GFWED

df_used <- readRDS("data/df_geff_erai_gfwed_era5_stvl_python_clean.rds")

minimum_size <- 1
maximum_bias <- 200
scaling_function <- function(x) {sqrt(x)}

df_to_map_era5 <- df_used %>%
  # Split tzid into region and subregion
  mutate(region = sapply(strsplit(tzid, "/"), `[`, 1),
         subregion = sapply(strsplit(tzid, "/"), `[`, 2)) %>%
  filter(region != "Etc") %>% # remove undefined zones
  filter(!is.na(OBS), !is.na(ERA5)) %>% # make sure you have data to compare
  group_by(id) %>%
  add_tally() %>% # add station count
  filter(n >= 30) %>% # remove stations with less than 30 days in a year
  summarise(lat = mean(lat),
            long = mean(long),
            region = names(which.max(table(region))),
            obs = round(mean(OBS, na.rm = TRUE), 2),
            erai = round(mean(ERAI, na.rm = TRUE), 2),
            gfwed = round(mean(GFWED, na.rm = TRUE), 2),
            era5 = round(mean(ERA5, na.rm = TRUE), 2),
            # Bias and Anomaly correlation (and p-values)
            bias_era5 = round(mean(ERA5 - OBS, na.rm = TRUE), 2),
            ac_era5 = round(cor(OBS - mean(OBS, na.rm = TRUE),
                                ERA5 - mean(ERA5, na.rm = TRUE)), 2),
            p_value = cor.test(OBS - mean(OBS, na.rm = TRUE),
                               ERA5 - mean(ERA5, na.rm = TRUE))$p.value < 0.05) %>%
  mutate(color = ifelse(abs(bias_era5) >= maximum_bias, 1,
                        ifelse(p_value == TRUE, 2, 3)),
         radius = ifelse(abs(bias_era5) <= 1, minimum_size,
                         ifelse(abs(bias_era5) >= maximum_bias,
                                minimum_size + scaling_function(maximum_bias),
                                minimum_size +
                                  scaling_function(abs(bias_era5)))),
         opacity = abs(1 - ac_era5))

idx <- as.numeric(names(table(df_to_map_era5$color)))
pal <- colorFactor(palette = c("grey", "navy", "red")[idx],
                   domain = df_to_map_era5$color)
labels <- c("Unreliable observations",
            "Cor. is stat. significant",
            "Cor. is not stat. significant")[idx]

# Check locations on interactive map
m <- leaflet(data = df_to_map_era5) %>%
  # Base groups
  # addTiles(group = "OpenStreetMap") %>%
  addProviderTiles(providers$CartoDB.Positron, group = "CartoDB (default)") %>%
  addProviderTiles(providers$OpenTopoMap, group = "OpenTopoMap") %>%
  addCircleMarkers(~long, ~lat,
                   radius = ~radius,
                   color = ~pal(color),
                   fillOpacity = ~opacity,
                   stroke = FALSE,
                   popup = ~paste("<strong>", "Bias:", "</strong>",
                                  bias_era5, "<br>",
                                  "<strong>", "Anomaly correlation:",
                                  "</strong>", ac_era5, "<br>",
                                  "<strong>", "Statistically significant:",
                                  "</strong>", p_value),
                   group = "Validation points") %>%
  setMaxBounds(lng1 = -180, lat1 = -90, lng2 = 180, lat2 = 90) %>%
  addLegend(position = "topright",
            pal = pal,
            title = "GEFF-ERA5 vs OBS",
            values = ~color,
            opacity = 1,
            labFormat = function(type, cuts, p) {  # Here's the trick
              paste0(labels)
            }
  ) %>%
  # Layers control
  addLayersControl(
    baseGroups = c("CartoDB (default)", "OpenTopoMap"),
    overlayGroups = "Validation points",
    options = layersControlOptions(collapsed = FALSE)
  )

# PUBLISH ON RPUBS THE INTERACTIVE MAP, THEN TAKE A SCREENSHOT FOR THE PAPER
saveWidget(m, file = "GEFF-ERA5_2017_diagnostic_map.html", selfcontained = TRUE)

############################# TABLE 1 ##########################################

df_to_compare <- df_used %>%
  # Split tzid into region and subregion
  mutate(region = sapply(strsplit(tzid, "/"), `[`, 1),
         subregion = sapply(strsplit(tzid, "/"), `[`, 2)) %>%
  filter(region != "Etc") %>% # remove undefined zones
  filter(!is.na(OBS), !is.na(ERAI), !is.na(GFWED), !is.na(ERA5)) %>%
  group_by(id) %>%
  add_tally() %>% # add station count
  filter(n >= 30) %>% # remove stations with less than 30 days in a year
  summarise(lat = as.numeric(names(which.max(table(lat)))),
            long = as.numeric(names(which.max(table(long)))),
            region = names(which.max(table(region))),
            obs = round(mean(OBS, na.rm = TRUE), 2),
            erai = round(mean(ERAI, na.rm = TRUE), 2),
            gfwed = round(mean(GFWED, na.rm = TRUE), 2),
            era5 = round(mean(ERA5, na.rm = TRUE), 2),
            # Bias and Anomaly correlation by time, id
            bias_erai = round(mean(ERAI - OBS, na.rm = TRUE), 2),
            bias_gfwed = round(mean(GFWED - OBS, na.rm = TRUE), 2),
            bias_era5 = round(mean(ERA5 - OBS, na.rm = TRUE), 2),
            ac_erai = round(cor(OBS - mean(OBS, na.rm = TRUE),
                                ERAI - mean(ERAI, na.rm = TRUE)), 2),
            ac_gfwed = round(cor(OBS - mean(OBS, na.rm = TRUE),
                                 GFWED - mean(GFWED, na.rm = TRUE)), 2),
            ac_era5 = round(cor(OBS - mean(OBS, na.rm = TRUE),
                                ERA5 - mean(ERA5, na.rm = TRUE)), 2),
            p_value_erai = cor.test(OBS - mean(OBS, na.rm = TRUE),
                                    ERAI - mean(ERAI, na.rm = TRUE))$p.value < 0.05,
            p_value_gfwed = cor.test(OBS - mean(OBS, na.rm = TRUE),
                                     GFWED - mean(GFWED, na.rm = TRUE))$p.value < 0.05,
            p_value_era5 = cor.test(OBS - mean(OBS, na.rm = TRUE),
                                    ERA5 - mean(ERA5, na.rm = TRUE))$p.value < 0.05) %>%
  filter(abs(bias_era5) < maximum_bias,
         p_value_erai == TRUE, p_value_gfwed == TRUE, p_value_era5 == TRUE)

# Summary table for large regions - copy-paste this into latex main.tex
dfx_to_table <- df_to_compare %>%
  group_by(region) %>%
  summarise(n = n(), # Bias and Anomaly correlation by time, id and
            bias_erai = round(mean(bias_erai), 2),
            bias_gfwed = round(mean(bias_gfwed), 2),
            bias_era5 = round(mean(bias_era5), 2),
            ac_erai = round(mean(ac_erai), 2),
            ac_gfwed = round(mean(ac_gfwed), 2),
            ac_era5 = round(mean(ac_era5), 2))

print(xtable::xtable(x = dfx_to_table, caption = "Validation"),
      include.rownames = FALSE)

############################# FIGURE 3 #########################################

# Explore global distributions
x <- reshape2::melt(data = df_to_compare[, c("region", "obs",
                                             "erai", "gfwed", "era5")],
                    id.vars = "region")

means <- aggregate(value ~  variable + region, x, mean)

p <- ggplot(x, aes(x = variable, y = value)) +
  geom_boxplot(outlier.shape = NA, width = 0.8) +
  geom_point(data = means, aes(x = variable, y = value),
             shape = 21, colour = "black", fill = "white") + #col = "red"
  facet_wrap(~region, scales = "free_y") +
  xlab("") + ylab("FWI") + theme_bw() +
  scale_x_discrete(labels = c("OBS", "ERAI", "GFWED", "ERA5")) +
  theme(text = element_text(size = 20), legend.position = "none")
ggsave(filename = "images/boxplots.eps",
       plot = p,
       device = "eps",
       width = 300,
       height = 250,
       units = "mm",
       dpi = 300)

# Calculate median for Atlantic
x %>% filter(region == "Atlantic") %>%
  group_by(variable) %>%
  summarise(meanFWI = median(value))

# Calculate mean for Atlantic
x %>% filter(region == "Atlantic") %>%
  group_by(variable) %>%
  summarise(meanFWI = mean(value))

rm(list = ls())

############################# FIGURE 4 #########################################

df_ens <- data.frame(matrix(NA, nrow = 0, ncol = 16))
dates_2017 <- seq.Date(from = as.Date("2017-01-01"),
                       to = as.Date("2017-12-31"),
                       by = "day")
for (i in seq_along(dates_2017)){

  one_date <- dates_2017[i]
  print(one_date)

  # Get ERA5 ENS
  era5 <- raster::stack()
  for (ens_member in 0:9){
    era5 <- raster::stack(era5,
                          file.path("/hugetmp/reanalysis/GEFF-ERA5/ens",
                                    paste0("ECMWF_FWI_",
                                           gsub("-", "", one_date),
                                           "_1200_0", ens_member, "_fwi.nc")))
  }
  era5 <- raster::rotate(era5)

  # Filter points with OBS on a given date
  dfx_date <- dfx_points[dfx_points$date == one_date, ]
  # Add placeholders columns
  dfx_date[, 7:16] <- NA
  # Convert to sp objects
  sp::coordinates(dfx_date) <- ~long+lat
  # Point inspection
  dfx_date@data[, 5:14] <- raster::extract(x = era5, y = dfx_date)

  x <- as.data.frame(dfx_date)
  x <- x[complete.cases(x),]
  df_ens <- rbind(df_ens, x)

}

saveRDS(df_ens, "data/df_ens.rds")

df_melt <- df_ens %>%
  reshape2::melt(id.vars = c("region", "lat", "long", "date", "OBS", "ERA5"))

ggplot(df_melt) +
  facet_wrap(~region, ncol = 4, scale = "free") +
  # OBS
  stat_ecdf(aes(x = OBS, linetype = "Observations"), lwd = 0.5) +
  # ERA5
  stat_ecdf(aes(x = ERA5, linetype = "HRES Reanalysis"), lwd = 0.5) +
  # ENS
  stat_ecdf(aes(x = value, colour = region, group = variable)) +
  scale_linetype_discrete(name = "") +
  theme(text = element_text(size=20)) +
  xlab("") + ylab("") + xlim(0, 100) + ylim(0.5, 1) +
  scale_colour_discrete(name = "ENS Reanalysis", h = c(0, 360), c = 80, l = 60)

# Comparison with ENSO #########################################################

# Crop reanalysis over SE Asia
system("cdo sellonlatbox,90,132,-14,21 /scratch/rd/nen/perClaudia/era5/fwi_1980_2019.nc /perm/mo/moc0/repos/GEFF-ERA5/data/fwi_era5_seasia.nc")

r <- brick("data/fwi_era5_seasia.nc")
seasia_clima_98 <- caliver::daily_clima(r, probs = 0.98)
writeRaster(seasia_clima_98, filename = "data/clima_98_seasia.nc",
            format = "CDF", overwrite = TRUE)

days_above_98 <- stack()
for (myyear in 1980:2018){
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
writeRaster(days_above_98, filename = "data/days_above_98_seasia.nc",
            format = "CDF", overwrite = TRUE)

# https://ggweather.com/enso/oni.htm
df_ens <- na.omit(setNames(as.data.frame(days_above_98), 1980:2018)) %>%
  mutate(P = 98) %>%
  reshape2::melt(id.vars = "P") %>%
  mutate(P = as.factor(P),
         col = ifelse(variable %in% c(1982, 1983, 1997, 1998, 2015),
                      "red",
                      ifelse(variable %in% c(1988, 1989, 1999, 2000,
                                             2007, 2008, 2010, 2011), "blue",
                             "gray")),
         alph = ifelse(variable %in% c(1982, 1983, 1988, 1989, 1997, 1998, 1999,
                                       2000, 2007, 2008, 2010, 2011, 2015),
                       0.5, 0.1))

############################# FIGURE 5 #########################################

ggplot(df_ens, aes(x = variable, y = value, fill = col)) +
  geom_boxplot(outlier.alpha = 0) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab("") + ylab("Number of days above 98th percentile") + ylim(0, 60) +
  scale_fill_manual(name = "",
                    labels = c("La Nina", "Either", "El Nino"),
                    breaks = c("blue", "gray", "red"),
                    values = c("blue", "gray", "red"))

############################# FIGURE 6 #########################################

# Map of days above threshold
myMap <- get_stamenmap(bbox = c(left = bbox(days_above_98)[[1]],
                                bottom = bbox(days_above_98)[[2]],
                                right = bbox(days_above_98)[[3]],
                                top = bbox(days_above_98)[[4]]),
                       maptype="toner-lite", color="bw", zoom=5, crop = T)

myplot <- ggmap(myMap) + xlab("Longitude") + ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))

# For year 1997
rtp_1997 <- rasterToPolygons(days_above_98[[18]])
names(rtp_1997) <- "layer"
common_max <- max(rtp_1997$layer) + 10
rtp_1997$layer <- cut(rtp_1997$layer, seq(0, common_max, 10),
                      include.lowest = TRUE)
# For year 2010
rtp_2010 <- rasterToPolygons(days_above_98[[31]])
names(rtp_2010) <- "layer"
rtp_2010$layer <- cut(rtp_2010$layer, seq(0, common_max, 10),
                      include.lowest = TRUE)

# Choose palette
common_palette <- colorspace::sequential_hcl(12, palette = "RedYellow",
                                             rev = TRUE)

myplot + ggtitle("1997") +
  geom_polygon(data = rtp_1997,
               aes(x = long, y = lat, group = group,
                   fill = rep(rtp_1997$layer, each = 5)),
               size = 0, alpha = 0.7) +
  scale_fill_manual(name = "Number of days\nabove 98th\npercentile",
                    values = common_palette,
                    drop = FALSE) +
  theme(text = element_text(size=20))

myplot + ggtitle("2010") +
  geom_polygon(data = rtp_2010,
               aes(x = long, y = lat, group = group,
                   fill = rep(rtp_2010$layer, each = 5)),
               size = 0, alpha = 0.7) +
  scale_fill_manual(name = "Number of days\nabove 98th\npercentile",
                    values = common_palette,
                    drop = FALSE) +
  theme(text = element_text(size=20))
