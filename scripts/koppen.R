climate_map <- raster("data/Beck_KG_V1/Beck_KG_V1_present_0p5.tif")
climate_map_remapped <- resample(x = climate_map, y = bias_gfwed, method = "ngb")
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

# Tests
x <- bias_erai
x[x <= 10] <- NA

plot(x, breaks = breaks_bias, col = pal_bias)
plot(gfed, add = TRUE)
i <- 6
region <- gfed[i, ]; plot(region, add = TRUE, col = NA, border = "red")

errors_masked <- mask(x, region)
climate_map_masked <- mask(climate_map_remapped, region)
df_mean <- data.frame(zonal(x = errors_masked, z = climate_map_masked, "mean"))
df_sum <- data.frame(zonal(x = errors_masked, z = climate_map_masked, "sum"))
df_mean$counts <- df_sum$sum/df_mean$mean
df_mean[order(df_mean$counts, decreasing = TRUE), ]

for (i in gfed@data$ID){
  region <- gfed[i, ]
  errors_masked <- mask(errors, region)
  climate_map_masked <- mask(climate_map_remapped, region)
  df <- data.frame(zonal(x = errors_masked, z = climate_map_masked, "mean"))
  df$region <- region@data$Region
  if (exists("df_all")){
    df_all <- rbind(df_all, df)
  }else{
    df_all <- df
  }
}

print(xtable::xtable(x = x[, c(1, 3, 2, 4, 5, 6)], caption = "Validation"),
      include.rownames = FALSE)