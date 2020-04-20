# LOAD PACKAGES ################################################################

library("raster")
# library("caliver")
# library("dplyr")
# library("leaflet")
# library("ggplot2")
# library("httr")
# library("lubridate")
# library("ggmap")
# library("colorspace")
# library("htmlwidgets")
# library("lutz")

# Get data from GFWED, ERAI AND ERA5 for year 2017 #############################

# FWI based on ERA5
# Fire danger indices are available from the Copernicus Climate Data Store
system(paste0("cdo cat /hugetmp/reanalysis/GEFF-ERA5/hres/fwi/ECMWF_FWI_2017* ",
              "/perm/mo/moc0/repos/GEFF-ERA5/data/fwi2017era5_temp.nc"))
# The dataset needs to be rotated because the original longitude range is
# [0,360] while erai and gfwed are in [-180, +180]
system(paste0("cdo sellonlatbox,-180,+180,-90,+90 ",
              "/perm/mo/moc0/repos/GEFF-ERA5/data/fwi2017era5_temp.nc ",
              "/perm/mo/moc0/repos/GEFF-ERA5/data/fwi2017era5.nc"))

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
system("cdo remapnn,data/fwi2017era5.nc data/fwi2017erai.nc data/fwi2017erai_remapped.nc")

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

rm(list = ls())

# TECHNICAL VALIDATION #########################################################

# GET FWI from reanalysis datasets
era5 <- raster::brick("data/fwi2017era5.nc")
erai <- raster::brick("data/fwi2017erai_remapped.nc")
gfwed <- raster::brick("data/fwi2017gfwed_remapped.nc")

# GFWED is only available in certain periods in some regions,
# for a fair comparison ERA5 and ERAI will be masked using GFWED
era5_masked <- raster::mask(x = era5, mask = gfwed, progress = "text"); rm(era5)
erai_masked <- raster::mask(x = erai, mask = gfwed, progress = "text"); rm(erai)

# Calculate means
mean_era5 <- raster::calc(x = era5_masked, fun = mean,
                          na.rm = TRUE, progress = "text")
mean_erai <- raster::calc(x = erai_masked, fun = mean,
                          na.rm = TRUE, progress = "text")
mean_gfwed <- raster::calc(x = gfwed, fun = mean,
                           na.rm = TRUE, progress = "text")

# Calculate anomalies
anom_era5 <- era5_masked - mean_era5
anom_erai <- erai_masked - mean_erai
anom_gfwed <- gfwed - mean_gfwed

# Calculate anomaly correlations and related p.values,
# mask-out cells with p-value > 0.05
ac_era5_erai <- raster::corLocal(anom_era5, anom_erai, test = TRUE)
ac_era5_erai_masked <- mask(x = ac_era5_erai[[1]],
                            mask = ac_era5_erai[[2]] > 0.05,
                            maskvalue = TRUE)
ac_era5_gfwed <- corLocal(anom_era5, anom_gfwed, test = TRUE)
ac_era5_gfwed_masked <- mask(x = ac_era5_gfwed[[1]],
                             mask = ac_era5_gfwed[[2]] > 0.05,
                             maskvalue = TRUE)

# Calculate biases
bias_era5_erai <- era5_masked - erai_masked
bias_era5_gfwed <- era5_masked - gfwed
mean_bias_era5_erai <- calc(x = bias_era5_erai, fun = mean, na.rm = TRUE)
mean_bias_era5_gfwed <- calc(x = bias_era5_gfwed, fun = mean, na.rm = TRUE)
mean_bias_era5_erai_masked <- mask(x = mean_bias_era5_erai,
                                    mask = ac_era5_erai[[2]] > 0.05,
                                    maskvalue = TRUE)
mean_bias_era5_gfwed_masked <- mask(x = mean_bias_era5_gfwed,
                                    mask = ac_era5_gfwed[[2]] > 0.05,
                                    maskvalue = TRUE)

writeRaster(mean_bias_era5_erai_masked,
            filename = "data/mean_bias_era5_erai.nc",
            format = "CDF", overwrite = TRUE)
writeRaster(mean_bias_era5_gfwed_masked,
            filename = "data/mean_bias_era5_gfwed.nc",
            format = "CDF", overwrite = TRUE)
writeRaster(ac_era5_erai_masked,
            filename = "data/ac_era5_erai_masked.nc",
            format = "CDF", overwrite = TRUE)
writeRaster(ac_era5_gfwed_masked,
            filename = "data/ac_era5_gfwed_masked.nc",
            format = "CDF", overwrite = TRUE)

############################# FIGURE 1: CDS screenshot #########################

# Screenshot of the CDS web interface

############################# FIGURE 2 #########################################

library("raster")
library("rnaturalearth")
library("rnaturalearthdata")
library("colorspace")

bias_era5_gfwed <- raster("data/mean_bias_era5_gfwed.nc")
bias_era5_erai <- raster("data/mean_bias_era5_erai.nc")
ac_era5_gfwed <- raster("data/ac_era5_gfwed_masked.nc")
ac_era5_erai <- raster("data/ac_era5_erai_masked.nc")

climate_map <- raster("/perm/mo/moc0/Beck_KG_V1/Beck_KG_V1_present_0p5.tif")
climate_map_remapped <- resample(x = climate_map, y = bias_era5_gfwed, method = "ngb")
code = c("Af   Tropical, rainforest",
               "Am   Tropical, monsoon",
               "Aw   Tropical, savannah",
               "BWh  Arid, desert, hot",
               "BWk  Arid, desert, cold",
               "BSh  Arid, steppe, hot",
               "BSk  Arid, steppe, cold",
               "Csa  Temperate, dry summer, hot summer",
               "Csb  Temperate, dry summer, warm summer",
               "Csc  Temperate, dry summer, cold summer",
               "Cwa  Temperate, dry winter, hot summer",
               "Cwb  Temperate, dry winter, warm summer",
               "Cwc  Temperate, dry winter, cold summer",
               "Cfa  Temperate, no dry season, hot summer",
               "Cfb  Temperate, no dry season, warm summer",
               "Cfc  Temperate, no dry season, cold summer",
               "Dsa  Cold, dry summer, hot summer",
               "Dsb  Cold, dry summer, warm summer",
               "Dsc  Cold, dry summer, cold summer",
               "Dsd  Cold, dry summer, very cold winter",
               "Dwa  Cold, dry winter, hot summer",
               "Dwb  Cold, dry winter, warm summer",
               "Dwc  Cold, dry winter, cold summer",
               "Dwd  Cold, dry winter, very cold winter",
               "Dfa  Cold, no dry season, hot summer",
               "Dfb  Cold, no dry season, warm summer",
               "Dfc  Cold, no dry season, cold summer",
               "Dfd  Cold, no dry season, very cold winter",
               "ET   Polar, tundra",
               "EF   Polar, frost")
climate_legend <- data.frame(zone = 1:30, code)

# MEAN
func <- "mean"
df_bias_era5_gfwed <- data.frame(zonal(x = bias_era5_gfwed,
                                       z = climate_map_remapped, func))
names(df_bias_era5_gfwed)[2] <- "bias_era5_gfwed"
df_bias_era5_erai <- data.frame(zonal(x = bias_era5_erai,
                                      z = climate_map_remapped, func))
names(df_bias_era5_erai)[2] <- "bias_era5_erai"
df_ac_era5_gfwed <- data.frame(zonal(x = ac_era5_gfwed,
                                     z = climate_map_remapped, func))
names(df_ac_era5_gfwed)[2] <- "ac_era5_gfwed"
df_ac_era5_erai <- data.frame(zonal(x = ac_era5_erai,
                                    z = climate_map_remapped, func))
names(df_ac_era5_erai)[2] <- "ac_era5_erai"
x <- merge(df_bias_era5_gfwed, climate_legend, by = "zone", all = TRUE)
x <- merge(x, df_bias_era5_erai, by = "zone", all = TRUE)
x <- merge(x, df_ac_era5_gfwed, by = "zone", all = TRUE)
x <- merge(x, df_ac_era5_erai, by = "zone", all = TRUE)
x <- x[complete.cases(x), ]
print(xtable::xtable(x = x[, c(1, 3, 2, 4, 5, 6)], caption = "Validation"),
      include.rownames = FALSE)

# Set up breaks
breaks_bias <- c(-60, -30, -15, -10, -5, 5, 10, 15, 30, 60)
breaks_ac <- c(-1, -0.6, -0.2, 0, 0.2, 0.6, 1)
# Set up palettes (Bias = Broc, AC = Oslo)
pal_bias <- rev(colorspace::diverge_hcl(length(breaks_bias) - 1, palette = "Blue-Red"))
pal_ac <- colorspace::sequential_hcl(length(breaks_ac) - 1, palette = "Inferno")

# STATIC FIGURES
ratioWH <- 1.77
W <- 10
H <- W/ratioWH

cairo_ps("images/bias_era5_erai.eps", width = W, height = H)
plot(bias_era5_erai, breaks = breaks_bias,
     main = "(a) ERA5 vs ERAI",
     col = pal_bias, legend = FALSE)
plot(coastline110, add = TRUE)
legend("bottomleft", horiz = FALSE, inset = 0.01, title = "Mean bias",
       fill = pal_bias, border = "gray", box.col = NA,
       legend = c("[-60, -30[", "[-30, -15[", "[-15, -10[", "[-10, -5[",
                  "[-5, +5[",
                  "[+5, +10[", "[+10, +15[", "[+15, +30[", "[+30, +60]"))
dev.off()

cairo_ps("images/bias_era5_gfwed.eps", width = W, height = H)
plot(bias_era5_gfwed, breaks = breaks_bias,
     main = "(b) ERA5 vs GFWED",
     col = pal_bias, legend = FALSE)
plot(coastline110, add = TRUE)
legend("bottomleft", horiz = FALSE, inset = 0.01, title = "Mean bias",
       fill = pal_bias, border = "gray", box.col = NA,
       legend = c("[-60, -30[", "[-30, -15[", "[-15, -10[", "[-10, -5[",
                  "[-5, +5[",
                  "[+5, +10[", "[+10, +15[", "[+15, +30[", "[+30, +60]"))
dev.off()

cairo_ps("images/ac_era5_erai.eps", width = W, height = H)
plot(ac_era5_erai, breaks = breaks_ac,
     main = "(a) ERA5 vs ERAI",
     col = pal_ac, legend = FALSE)
plot(coastline110, add = TRUE)
legend("bottomleft", horiz = FALSE, inset = 0.01, title = "Anomaly correlation",
       fill = pal_ac, border = "gray", box.col = NA,
       legend = c("[-1.0, -0.6[", "[-0.6, -0.2[",
                  "[-0.2, 0[", "[0, +0.2[",
                  "[+0.2, +0.6[", "[+0.6, +1.0]"))
dev.off()

cairo_ps("images/ac_era5_gfwed.eps", width = W, height = H)
plot(ac_era5_gfwed, breaks = breaks_ac,
     main = "(b) ERA5 vs GFWED",
     col = pal_ac, legend = FALSE)
plot(coastline110, add = TRUE)
legend("bottomleft", horiz = FALSE, inset = 0.01, title = "Anomaly correlation",
       fill = pal_ac, border = "gray", box.col = NA,
       legend = c("[-1.0, -0.6[", "[-0.6, -0.2[",
                  "[-0.2, 0[", "[0, +0.2[",
                  "[+0.2, +0.6[", "[+0.6, +1.0]"))
dev.off()

# INTERACTIVE FIGURES
library(leaflet)

# Bias GFWED
library("mapview")
mapview::mapview(c(bias_era5_gfwed, climate_map))

mapview::mapview(bias_era5_gfwed, at = breaks_bias, col.regions = pal_bias,
                 alpha.regions = 0.5, map.types = mapviewGetOption("OpenTopoMap"))

leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron, group = "CartoDB (default)") %>%
  addProviderTiles(providers$OpenTopoMap, group = "OpenTopoMap") %>%
  addRasterImage(bias_era5_gfwed,
                 colors = colorBin(palette = pal_bias,
                                   bins = breaks_bias,
                                   na.color = "transparent"),
                 opacity = 0.5, group = "mean_bias_era5_gfwed") %>%
  addRasterImage(climate_map,
                 group = "climate_map") %>%
  addLayersControl(
    baseGroups = c("CartoDB (default)", "OpenTopoMap"),
    overlayGroups = c("mean_bias_era5_gfwed", "climate_map"),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  addMiniMap()

# Bias ERAI
leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron, group = "CartoDB (default)") %>%
  addProviderTiles(providers$OpenTopoMap, group = "OpenTopoMap") %>%
  addRasterImage(bias_era5_erai,
                 colors = colorBin(palette = pal_bias,
                                   bins = breaks_bias,
                                   na.color = "transparent"),
                 opacity = 0.5, group = "mean_bias_era5_erai") %>%
  addRasterImage(climate_map,
                 group = "climate_map") %>%
  addLayersControl(
    baseGroups = c("CartoDB (default)", "OpenTopoMap"),
    overlayGroups = c("mean_bias_era5_erai", "climate_map"),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  addMiniMap() %>%
  addLegend(pal = pal_bias,
            values = values(bias_era5_erai),
            title = "Mean bias ERA5 vs ERAI")

# Check locations on interactive map
m <- leaflet(data = df_era5_gfwed) %>%
  # Base groups
  addTiles() %>%
  #addProviderTiles(providers$CartoDB.Positron, group = "CartoDB (default)") %>%
  #addProviderTiles(providers$OpenTopoMap, group = "OpenTopoMap") %>%
  addCircleMarkers(~x, ~y,
                   color = ~pal(bias),
                   fillOpacity = ~abs(1 - ac),
                   stroke = FALSE,
                   popup = ~paste("<strong>", "Bias:", "</strong>",
                                  bias, "<br>",
                                  "<strong>", "Anomaly correlation:",
                                  "</strong>", ac, "<br>"),
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

############################# FIGURE 3 boxplots of regional distributions ######

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

############################# FIGURE 4: ecdf ENS ###############################

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

############################# FIGURE 5: comparison with ENSO (boxplot) #########

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

ggplot(df_ens, aes(x = variable, y = value, fill = col)) +
  geom_boxplot(outlier.alpha = 0) +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab("") + ylab("Number of days above 98th percentile") + ylim(0, 60) +
  scale_fill_manual(name = "",
                    labels = c("La Nina", "Either", "El Nino"),
                    breaks = c("blue", "gray", "red"),
                    values = c("blue", "gray", "red"))

############################# FIGURE 6: comparison with ENSO (maps) ############

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
