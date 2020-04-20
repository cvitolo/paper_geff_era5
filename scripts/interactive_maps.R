# INTERACTIVE FIGURES ##########################################################
library(leaflet)

# Bias GFWED
library("mapview")
mapview::mapview(c(bias_gfwed, climate_map))

mapview::mapview(bias_gfwed, at = breaks_bias, col.regions = pal_bias,
                 alpha.regions = 0.5, map.types = mapviewGetOption("OpenTopoMap"))

leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron, group = "CartoDB (default)") %>%
  addProviderTiles(providers$OpenTopoMap, group = "OpenTopoMap") %>%
  addRasterImage(bias_gfwed,
                 colors = colorBin(palette = pal_bias,
                                   bins = breaks_bias,
                                   na.color = "transparent"),
                 opacity = 0.5, group = "mean_bias_gfwed") %>%
  addRasterImage(climate_map,
                 group = "climate_map") %>%
  addLayersControl(
    baseGroups = c("CartoDB (default)", "OpenTopoMap"),
    overlayGroups = c("mean_bias_gfwed", "climate_map"),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  addMiniMap()

# Bias ERAI
leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron, group = "CartoDB (default)") %>%
  addProviderTiles(providers$OpenTopoMap, group = "OpenTopoMap") %>%
  addRasterImage(bias_erai,
                 colors = colorBin(palette = pal_bias,
                                   bins = breaks_bias,
                                   na.color = "transparent"),
                 opacity = 0.5, group = "mean_bias_erai") %>%
  addRasterImage(climate_map,
                 group = "climate_map") %>%
  addLayersControl(
    baseGroups = c("CartoDB (default)", "OpenTopoMap"),
    overlayGroups = c("mean_bias_erai", "climate_map"),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  addMiniMap() %>%
  addLegend(pal = pal_bias,
            values = values(bias_erai),
            title = "Mean bias ERA5 vs ERAI")

# Check locations on interactive map
m <- leaflet(data = df_gfwed) %>%
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