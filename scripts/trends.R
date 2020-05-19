# Time series decomposition: trend + seasonality + random
# USES 1 ENSEMBLE MEMBER
library(dplyr)
library(ggplot2)
library(forecast)

df <- readRDS("df_PG.rds")
# Decomposition - 40 years #####################################################
df %>%
  # filter(lubridate::month(dates) == 6) %>% # only June
  select(em0) %>%
  msts(seasonal.periods = c(365.25)) %>%
  mstl() %>% autoplot()

df %>%
  # filter(lubridate::month(dates) == 6) %>% # only June
  select(em0) %>%
  msts(seasonal.periods = c(365.25, 365.25 * 10)) %>%
  mstl() %>% autoplot()

# Decomposition - 40 years - JUNE ONLY #########################################
df %>%
  filter(lubridate::month(dates) == 6) %>% # only June
  select(em0) %>%
  msts(seasonal.periods = c(30)) %>%
  mstl() %>% autoplot()

# Decomposition - 40 years - force linear trend ################################
em0 <- ts(data = df$em0, start = c(1980, 01), frequency = 365.25)
decompose_df <- tslm(em0 ~ trend + fourier(em0, 2))
trend <- coef(decompose_df)[1] + coef(decompose_df)['trend']*seq_along(em0)
components <- cbind(
  data = em0,
  trend = trend,
  season = em0 - trend - residuals(decompose_df),
  remainder = residuals(decompose_df)
)
autoplot(components, facet=TRUE)

# Decomposition - 40 years - JUNE ONLY - force linear trend ####################
df_june <- df[lubridate::month(df$dates) == 6, ]
emean <- ts(data = df_june$mean, start = c(1980, 01), frequency = 365.25)
decompose_df <- tslm(emean ~ trend + fourier(emean, 2))
trend <- coef(decompose_df)[1] + coef(decompose_df)['trend']*seq_along(emean)
components <- cbind(
  data = emean,
  trend = trend,
  season = emean - trend - residuals(decompose_df),
  remainder = residuals(decompose_df)
)
autoplot(components, facet=TRUE)
