# LOAD PACKAGES ################################################################

library("raster")
library("rnaturalearthdata")
library("colorspace")
library("dplyr")
library("lubridate")
library("xts")
library("ggplot2")
library("ggmap")
library("dygraphs")
library("forecast")

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
dates_2017 <- seq.Date(from = as.Date("2017-01-01"),
                       to = as.Date("2017-12-31"),
                       by = "day")

for (myday in seq_along(dates_2017)){
  print(myday)
  reformatted_date <- gsub(pattern = "-", replacement = "", dates_2017[myday])
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

# TECHNICAL VALIDATION #########################################################

# GET FWI from reanalysis hres datasets
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

# Calculate biases (masked as the AC plots, for visual consistency)
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

errors <- stack(mean_bias_era5_gfwed_masked, mean_bias_era5_erai_masked,
                ac_era5_gfwed_masked, ac_era5_erai_masked)
names(errors) <- c("bias_gfwed", "bias_erai", "ac_gfwed", "ac_erai")

writeRaster(errors, filename = "data/errors.nc",
            format = "CDF", overwrite = TRUE)

# Assess ensemble
era5_ens_mean_stack <- era5_ens_sd_stack <- stack()
for (i in seq_along(dates_2017)){

  one_date <- dates_2017[i]
  print(one_date)

  era5_hr <- raster::raster(file.path("/hugetmp/reanalysis/GEFF-ERA5/hres/fwi/",
                                      paste0("ECMWF_FWI_",
                                             gsub("-", "", one_date),
                                             "_1200_hr_fwi.nc")))

  # Get ERA5 ENS
  era5_ens <- raster::stack()
  for (ens_member in 0:9){
    era5_ens <- raster::stack(era5_ens,
                              file.path("/hugetmp/reanalysis/GEFF-ERA5/ens",
                                        paste0("ECMWF_FWI_",
                                               gsub("-", "", one_date),
                                               "_1200_0", ens_member, "_fwi.nc")))
  }
  # Mean
  era5_ens_mean <- calc(era5_ens, mean)
  era5_ens_mean_stack <- stack(era5_ens_mean_stack, era5_ens_mean)
  # Spread
  era5_ens_sd <- calc(era5_ens, sd)
  era5_ens_sd_stack <- stack(era5_ens_sd_stack, era5_ens_sd)

}
era5_ens_mean_stack_mean <- rotate(calc(era5_ens_mean_stack, mean))
era5_ens_sd_stack_mean <- rotate(calc(era5_ens_sd_stack, mean))

writeRaster(era5_ens_mean_stack_mean,
            filename = "data/era5_ens_mean_stack_mean.nc",
            format = "CDF", overwrite = TRUE)
writeRaster(era5_ens_sd_stack_mean, filename = "data/era5_ens_sd_stack_mean.nc",
            format = "CDF", overwrite = TRUE)

rm(list = ls())

# FIGURE 1: CDS screenshot #####################################################

# Screenshot of the CDS web interface

# FIGURE 2: HRES MEAN ##########################################################

era5_hres <- raster::brick("data/fwi2017era5.nc")
era5_hres_mean <- raster::calc(era5_hres, mean)

# Define size of static maps
ratioWH <- 1.77
W <- 10
H <- W/ratioWH

# Set up breaks
breaks_fwi <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, +Inf)

# Set up palettes
pal_fwi <- rev(colorspace::sequential_hcl(length(breaks_fwi) - 1,
                                          palette = "Red-Yellow"))

cairo_ps("images/hres_mean.eps", width = W, height = H)
plot(era5_hres_mean, breaks = breaks_fwi, main = "",
     col = pal_fwi, legend = FALSE)
plot(coastline110, add = TRUE)
legend("bottomleft", horiz = FALSE, inset = 0.01, title = "Mean of FWI HRES",
       fill = pal_fwi, border = "gray", box.col = NA,
       legend = c("[0, 10[", "[10, 20[", "[20, 30[", "[30, 40[", "[40, 50[",
                  "[50, 60[", "[60, 70[", "[70, 80[", "[80, 90[", "[90, 100+["))
dev.off()

# FIGURE 3: ENS MEAN AND SPREAD ################################################
era5_ens_mean <- raster::raster("data/era5_ens_mean_stack_mean.nc")
era5_ens_sd <- raster::raster("data/era5_ens_sd_stack_mean.nc")

cairo_ps("images/ens_mean.eps", width = W, height = H)
plot(era5_ens_mean, breaks = breaks_fwi, main = "",
     col = pal_fwi, legend = FALSE)
plot(coastline110, add = TRUE)
legend("bottomleft", horiz = FALSE, inset = 0.01,
       title = "Mean of FWI\n ENS mean",
       fill = pal_fwi, border = "gray", box.col = NA,
       legend = c("[0, 10[", "[10, 20[", "[20, 30[", "[30, 40[", "[40, 50[",
                  "[50, 60[", "[60, 70[", "[70, 80[", "[80, 90[", "[90, 100+["))
dev.off()

breaks_sd <- c(0, 1, 2, 3, 4, 5, 6, 7, 8)
pal_sd <- colorspace::sequential_hcl(length(breaks_sd) - 1, palette = "Oslo")

cairo_ps("images/ens_spread.eps", width = W, height = H)
plot(era5_ens_sd, breaks = breaks_sd, main = "",
     col = pal_sd, legend = FALSE)
plot(coastline110, add = TRUE)
legend("bottomleft", horiz = FALSE, inset = 0.01,
       title = "Mean of FWI\n ENS spread",
       fill = pal_sd, border = "gray", box.col = NA,
       legend = c("[0, 1[", "[1, 2[", "[2, 3[", "[3, 4[", "[4, 5[", "[5, 6[",
                  "[6, 7[", "[7, 8]"))
dev.off()

# FIGURES 4&5: BIAS and AC #####################################################

errors <- raster::brick("data/errors.nc")
names(errors) <- c("bias_gfwed", "bias_erai", "ac_gfwed", "ac_erai")

gfed <- readRDS("data/GFED4_BasisRegions.rds")

df <- data.frame(zonal(x = errors, z = rasterize(x = gfed, y = errors[[1]])))

df$Region <- c("Boreal North America",
               "Temperate North America",
               "Central America",
               "North Hemisphere South America",
               "South Hemisphere South America",
               "Europe",
               "Middle East",
               "North Hemisphere Africa",
               "South Hemisphere Africa",
               "Boreal Asia",
               "Central Asia",
               "Southeast Asia",
               "Equatorial Asia",
               "Australia and New Zealand")
df$Region_short <- gfed@data$Region
print(xtable::xtable(x = df[, c(6, 2:5)], caption = "Validation"),
      include.rownames = FALSE)

# Set up breaks
breaks_bias <- c(-60, -30, -15, -10, -5, 5, 10, 15, 30, 60)
breaks_ac <- c(-1, -0.6, -0.2, 0, 0.2, 0.6, 1)
# Set up palettes
pal_bias <- rev(colorspace::diverge_hcl(length(breaks_bias) - 1,
                                        palette = "Blue-Red"))
pal_ac <- colorspace::sequential_hcl(length(breaks_ac) - 1,
                                     palette = "Inferno")

cairo_ps("images/bias_erai.eps", width = W, height = H)
plot(bias_erai, breaks = breaks_bias, main = "(a) ERA5 vs ERAI",
     col = pal_bias, legend = FALSE)
plot(coastline110, add = TRUE)
legend("bottomleft", horiz = FALSE, inset = 0.01, title = "Mean bias",
       fill = pal_bias, border = "gray", box.col = NA,
       legend = c("[-60, -30[", "[-30, -15[", "[-15, -10[", "[-10, -5[",
                  "[-5, +5[",
                  "[+5, +10[", "[+10, +15[", "[+15, +30[", "[+30, +60]"))
dev.off()

cairo_ps("images/bias_era5_gfwed.eps", width = W, height = H)
plot(bias_gfwed, breaks = breaks_bias,
     main = "(b) ERA5 vs GFWED",
     col = pal_bias, legend = FALSE)
plot(coastline110, add = TRUE)
legend("bottomleft", horiz = FALSE, inset = 0.01, title = "Mean bias",
       fill = pal_bias, border = "gray", box.col = NA,
       legend = c("[-60, -30[", "[-30, -15[", "[-15, -10[", "[-10, -5[",
                  "[-5, +5[",
                  "[+5, +10[", "[+10, +15[", "[+15, +30[", "[+30, +60]"))
dev.off()

cairo_ps("images/ac_erai.eps", width = W, height = H)
plot(ac_erai, breaks = breaks_ac,
     main = "(a) ERA5 vs ERAI",
     col = pal_ac, legend = FALSE)
plot(coastline110, add = TRUE)
legend("bottomleft", horiz = FALSE, inset = 0.01, title = "Anomaly correlation",
       fill = pal_ac, border = "gray", box.col = NA,
       legend = c("[-1.0, -0.6[", "[-0.6, -0.2[",
                  "[-0.2, 0[", "[0, +0.2[",
                  "[+0.2, +0.6[", "[+0.6, +1.0]"))
dev.off()

cairo_ps("images/ac_gfwed.eps", width = W, height = H)
plot(ac_gfwed, breaks = breaks_ac,
     main = "(b) ERA5 vs GFWED",
     col = pal_ac, legend = FALSE)
plot(coastline110, add = TRUE)
legend("bottomleft", horiz = FALSE, inset = 0.01, title = "Anomaly correlation",
       fill = pal_ac, border = "gray", box.col = NA,
       legend = c("[-1.0, -0.6[", "[-0.6, -0.2[",
                  "[-0.2, 0[", "[0, +0.2[",
                  "[+0.2, +0.6[", "[+0.6, +1.0]"))
dev.off()

# FIGURE 6: ENS use case #######################################################

# Pedrogao Grande
spoint <- sf::st_as_sf(data.frame(long = -8.23, lat = 39.95),
                       coords = c("long", "lat"), crs = 4326)
plot(coastline110); plot(spoint, add = TRUE, col = "red")

for (i in seq_along(dates_2017)){

  one_date <- dates_2017[i]
  print(one_date)

  era5_hr <- raster::raster(file.path("/hugetmp/reanalysis/GEFF-ERA5/hres/fwi/",
                                      paste0("ECMWF_FWI_",
                                             gsub("-", "", one_date),
                                             "_1200_hr_fwi.nc")))
  era5_hr <- raster::rotate(era5_hr)
  era5_hr <- raster::extract(era5_hr, spoint)

  # Get ERA5 ENS
  era5_ens <- raster::stack()
  for (ens_member in 0:9){
    era5_ens <- raster::stack(era5_ens,
                              file.path("/hugetmp/reanalysis/GEFF-ERA5/ens",
                                        paste0("ECMWF_FWI_",
                                               gsub("-", "", one_date),
                                               "_1200_0", ens_member, "_fwi.nc")))
  }

  # Get data at point
  era5_ens <- as.numeric(raster::extract(rotate(era5_ens), spoint))
  x <- data.frame(hr = era5_hr,
                  ens_00 = era5_ens[1], ens_01 = era5_ens[2],
                  ens_02 = era5_ens[3], ens_03 = era5_ens[4],
                  ens_04 = era5_ens[5], ens_05 = era5_ens[6],
                  ens_06 = era5_ens[7], ens_07 = era5_ens[8],
                  ens_08 = era5_ens[9], ens_09 = era5_ens[10])

  if (exists("dfx")){
    dfx <- rbind(dfx, x)
  }else{
    dfx <- x
  }

}

dfx$date <- dates_2017
dfx$ens_min <- apply(dfx[, 2:11], 1, min)
dfx$ens_max <- apply(dfx[, 2:11], 1, max)
dfx$ens_mean <- apply(dfx[, 2:11], 1, mean)

# saveRDS(dfx, "data/dfx_ens.rds")
# dfx <- readRDS("data/dfx_ens.rds")

ggplot(dfx, aes(date)) +
  geom_ribbon(aes(ymin = ens_min, ymax = ens_max), fill = "grey80") +
  geom_line(aes(y = ens_mean, color = "ENS mean"), size = 0.5) +
  # geom_line(aes(y = hr, color = "HRES"), size = 0.5) +
  lims(x = c(as.Date("2017-06-01"), as.Date("2017-06-30")), y = c(0, 45)) +
  geom_vline(xintercept = as.Date("2017-06-17"),
             linetype = "dashed", color = "red", size = 1) +
  annotate(geom = "text", x = as.Date("2017-06-21"),
           y = 5, label = "Event in Pedrógão Grande", color = "red") +
  scale_color_discrete(name = "") + xlab("") + ylab("FWI") + theme_bw()

ggsave(filename = "images/PG_ens_2017.eps", plot = last_plot(),
       device = "eps", width = W, height = H, units = "in")

# Full time series and trends at Pedrogao Grande
dates <- seq.Date(from = as.Date("1980-01-01"),
                  to = as.Date("2019-12-31"),
                  by = "day")

ncfname <- "data/cds_fwi_ens_PG.nc"
x <- raster::raster(ncfname)
df <- as.data.frame(x, xy = TRUE) %>%
  tidyr::pivot_wider(names_from = "x", values_from = "Index") %>%
  dplyr::arrange(y)  %>%
  mutate(dates = dates)
names(df)[2:11] <- paste0("em", 0:9)
df$mean <- apply(df[, 2:11], 1, mean)
df$min <- apply(df[, 2:11], 1, min)
df$max <- apply(df[, 2:11], 1, max)
df$spread <- df$max - df$min

rm(x, dates, ncfname)

# Aggregate mean/min/max to weekly values
dfts <- xts::as.xts(df[, -12], order.by = as.Date(df$dates))
dfm <- xts::apply.monthly(dfts, mean)
# Keep June only
dfmjune <- dfm[lubridate::month(index(dfm)) == 6, ]

rm(dfts, dfm)

# Decomposition - 40 years - JUNE ONLY - force linear trend ####################
df_june <- df[lubridate::month(df$dates) == 6, ]
# Frequency is 1 month of daily data per year
emean <- ts(data = df_june$mean, start = c(1980, 01), frequency = 365.25/12)
decompose_df <- tslm(emean ~ trend + fourier(emean, 2))
trend <- coef(decompose_df)[1] + coef(decompose_df)['trend']*seq_along(emean)
components <- cbind(
  data = emean,
  trend = trend,
  season = emean - trend - residuals(decompose_df),
  remainder = residuals(decompose_df)
)
autoplot(components, facet=TRUE)

cairo_ps("images/trend_ens.eps", width = W, height = H)
barplot(dfmjune$spread)
barplot(height = c(rep(NA, 37),
                   dfmjune$spread[lubridate::year(index(dfmjune)) == 2017],
                   NA, NA),
        add = TRUE, col = "red")
abline(a = 0, b = coef(decompose_df)[2]*length(emean)/length(dfmjune$spread),
       col = "red", lwd = 2)
dev.off()

# More advanced plots
df_ts <- xts::xts(x = df[, c(2:11, 13)], order.by = df$dates)
dygraph(df_ts, main = "Daily FWI from GEFF-ERA5 ENS at Pedrógão Grande") %>%
  dyGroup(paste0("em", 0:9), color = rep("gray", 10)) %>%
  dySeries("mean", color = "red")
df_mon <- xts::apply.monthly(df_ts, FUN = mean)
dygraph(df_mon, main = "Monthly FWI from GEFF-ERA5 ENS at Pedrógão Grande") %>%
  dyGroup(paste0("em", 0:9), color = rep("gray", 10)) %>%
  dySeries("mean", color = "red")

# FIGURE 7: comparison with ENSO (boxplot) #####################################

# Crop reanalysis over SE Asia
system("cdo sellonlatbox,90,132,-14,21 /scratch/rd/nen/perClaudia/era5/fwi_1980_2019.nc /perm/mo/moc0/repos/GEFF-ERA5/data/fwi_seasia.nc")

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
# days_above_98 <- raster::brick("data/days_above_98_seasia.nc")

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

df_ens_plot <- ggplot(df_ens, aes(x = variable, y = value, fill = col)) +
  geom_boxplot(outlier.alpha = 0) +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab("") + ylab("Number of days above 98th percentile") + ylim(0, 60) +
  scale_fill_manual(name = "",
                    labels = c("La Nina", "Either", "El Nino"),
                    breaks = c("blue", "gray", "red"),
                    values = c("blue", "gray", "red"))

ggsave(filename = "images/enso_timeseries_sea.eps", plot = df_ens_plot,
       device = "eps", width = W, height = H, units = "in")

# FIGURE 8: comparison with ENSO (maps) ########################################

# Map of days above threshold
myMap <- ggmap::get_stamenmap(bbox = c(left = bbox(days_above_98)[[1]],
                                       bottom = bbox(days_above_98)[[2]],
                                       right = bbox(days_above_98)[[3]],
                                       top = bbox(days_above_98)[[4]]),
                              maptype = "toner-lite",
                              color = "bw", zoom = 5, crop = T)

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

myplot + ggtitle("(a) Very strong positive ENSO in 1997") +
  geom_polygon(data = rtp_1997,
               aes(x = long, y = lat, group = group,
                   fill = rep(rtp_1997$layer, each = 5)),
               size = 0, alpha = 0.7) +
  scale_fill_manual(name = "Number of days\nabove 98th\npercentile",
                    values = common_palette,
                    drop = FALSE) +
  theme(text = element_text(size = 12))

ggsave(filename = "images/Y1997.pdf", plot = last_plot(),
       device = "pdf", width = W, height = H, units = "in")

myplot + ggtitle("(b) Very strong negative ENSO in 2010") +
  geom_polygon(data = rtp_2010,
               aes(x = long, y = lat, group = group,
                   fill = rep(rtp_2010$layer, each = 5)),
               size = 0, alpha = 0.7) +
  scale_fill_manual(name = "Number of days\nabove 98th\npercentile",
                    values = common_palette,
                    drop = FALSE) +
  theme(text = element_text(size = 12))

ggsave(filename = "images/Y2010.pdf", plot = last_plot(),
       device = "pdf", width = W, height = H, units = "in")
