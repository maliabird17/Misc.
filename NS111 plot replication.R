#### Loading packages and data and processing from Nature paper ####
library("mgcViz")
library("tidyverse")
library("mgcv")
library("gratia")
library("patchwork")
library("scales")
library("RColorBrewer")
# Code and data from: https://github.com/dm807cam/nucella2021belgium

### Temperature data ------
# Import data
belgium_air_temperature <- read.csv("https://raw.githubusercontent.com/dm807cam/nucella2021belgium/main/data/env/air_temperature_comp.csv", sep = ",", header = T)
glimpse(belgium_air_temperature)

belgium_sea_temperature <- read.csv("https://raw.githubusercontent.com/dm807cam/nucella2021belgium/main/data/env/sea_surface_temperature.csv", sep = ",", header = T)
glimpse(belgium_sea_temperature)

belgium_wind_speed <- read.csv("https://raw.githubusercontent.com/dm807cam/nucella2021belgium/main/data/env/wind_speed_comp.csv", sep = ",", header = T)
glimpse(belgium_wind_speed)

# Functions -------
signifD <- function(x, d, upper, lower, eval = 0) {
  miss <- upper > eval & lower < eval
  incr <- decr <- x
  want <- d > eval
  incr[!want | miss] <- NA
  want <- d < eval
  decr[!want | miss] <- NA
  list(incr = incr, decr = decr)
}
# GAMMs for temperature and wind speed time series -------------------------
# Extend controls
ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")

# Fit Air temperature GAMM
gamm_air_temp <- gamm(air_temperature_celsius ~ s(year, bs = "tp") + 
                        s(month, bs = "cr", k = 12) + 
                        te(year,month, bs = c("tp", "cr")), 
                      data = belgium_air_temperature, 
                      method="REML", 
                      family=gaussian(link = "identity"),
                      correlation = corARMA(form = ~ 1|year, p=2),
                      control = ctrl)

# Check for autocorrelation
layout(matrix(1:2, ncol = 2))
res <- resid(gamm_air_temp$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(2) errors")
pacf(res, lag.max = 36, main = "pACF- AR(2) errors")
layout(1)

# Check model diagnostics
appraise(gamm_air_temp$gam)
draw(gamm_air_temp$gam)
summary(gamm_air_temp$gam)

# ---

# Fit Sea surface temperature GAMM
gamm_sea_temp <- gamm(sea_temperature_celsius ~ 
                        s(year, bs = "tp") + 
                        s(month, bs = "cr", k = 12) + 
                        te(year,month, bs = c("tp", "cr")), 
                      data = belgium_sea_temperature, 
                      method="REML", 
                      family=gaussian(link = "identity"),
                      correlation = corARMA(form = ~ 1|year, p=1),
                      control = ctrl)

# Check for autocorrelation
layout(matrix(1:2, ncol = 2))
res <- resid(gamm_sea_temp$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(2) errors")
pacf(res, lag.max = 36, main = "pACF- AR(2) errors")
layout(1)

# Check model diagnostics
appraise(gamm_sea_temp$gam)
draw(gamm_sea_temp$gam)
summary(gamm_sea_temp$gam)

# Export all models -------------------------------------------------------
# Predict air temperature -------------------------------------------------

want <- seq(1, nrow(belgium_air_temperature),
            length.out = 200)

pdatAT <- with(belgium_air_temperature,
               data.frame(year = year[want],
                          month = month[want]))

p2 <- data.frame(predict(
  gamm_air_temp$gam,
  newdata = pdatAT,
  type = "terms",
  se.fit = TRUE
))

pdatAT <- transform(pdatAT,
                    p2 = p2[, 1],
                    se2 = p2[, 4])

df.res <- df.residual(gamm_air_temp$gam)

crit.t <- qt(0.025, df.res, lower.tail = FALSE)

pdatAT <- transform(pdatAT,
                    upper = p2 + (crit.t * se2),
                    lower = p2 - (crit.t * se2))

gamm_air_temp.d <- fderiv(gamm_air_temp)

gamm_air_temp.dci <- confint(gamm_air_temp.d,
                             parm = "year",
                             type = "confidence")

gamm_air_temp.dsig <- signifD(
  pdatAT$p2,
  d = gamm_air_temp.d$derivatives$`s(year)`$deriv,
  gamm_air_temp.dci$upper,
  gamm_air_temp.dci$lower
)

# Predict sea surface temperature -------------------------------------------------

want <- seq(1, nrow(belgium_air_temperature),
            length.out = 200)

pdatST <- with(belgium_air_temperature,
               data.frame(year = year[want],
                          month = month[want]))

p2 <- data.frame(predict(
  gamm_sea_temp$gam,
  newdata = pdatST,
  type = "terms",
  se.fit = TRUE
))

pdatST <- transform(pdatST,
                    p2 = p2[, 1],
                    se2 = p2[, 4])

df.res <- df.residual(gamm_sea_temp$gam)

crit.t <- qt(0.025, df.res, lower.tail = FALSE)

pdatST <- transform(pdatST,
                    upper = p2 + (crit.t * se2),
                    lower = p2 - (crit.t * se2))

gamm_sea_temp.d <- fderiv(gamm_sea_temp)

gamm_sea_temp.dci <-
  confint(gamm_sea_temp.d, parm = "year", type = "confidence")

gamm_sea_temp.dsigST <- signifD(
  pdatST$p2,
  d = gamm_sea_temp.d$derivatives$`s(year)`$deriv,
  gamm_sea_temp.dci$upper,
  gamm_sea_temp.dci$lower
)

# Plot air temperature -------------------------------------------------
ATC <- ggplot(pdatAT, aes(year,
                          p2)) +
  geom_line(linetype = 2) +
  geom_ribbon(aes(ymin = lower,
                  ymax = upper),
              fill = "black",
              alpha = .15) +
  geom_line(
    pdatAT,
    mapping = aes(year,
                  unlist(gamm_air_temp.dsig$incr)),
    colour = brewer.pal(n = 4, name = "Spectral")[1]
  ) +
  geom_line(
    pdatAT,
    mapping = aes(year,
                  unlist(gamm_air_temp.dsig$decr)),
    colour = brewer.pal(n = 4, name = "Spectral")[4]
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "(a1)",
    vjust = 1.5,
    hjust = -.25,
    parse = TRUE,
    size = 3,
    family = "Times"
  ) +
  scale_fill_grey(start = 0.2,
                  end = 0.8) +
  theme(text = element_text(size = 10, family="Times"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank()
  ) +
  ylab(expression(AT ~ (degree * C * ":" ~ centred))) +
  scale_x_continuous(expand = c(0.015, 0.015),
                     breaks = pretty_breaks(n = 6))

ATC

# Plot sea surface temperature -------------------------------------------------
STC <- ggplot(pdatST, aes(year,
                          p2)) +
  geom_line(linetype = 2) +
  geom_ribbon(aes(ymin = lower,
                  ymax = upper),
              fill = "black",
              alpha = .15) +
  geom_line(
    pdatST,
    mapping = aes(year,
                  unlist(gamm_sea_temp.dsigST$incr)),
    colour = brewer.pal(n = 4, name = "Spectral")[1]
  ) +
  geom_line(
    pdatST,
    mapping = aes(year, unlist(gamm_sea_temp.dsigST$decr)),
    colour = brewer.pal(n = 4, name = "Spectral")[4]
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "(b1)",
    vjust = 1.5,
    hjust = -.25,
    parse = TRUE,
    size = 3,
    family = "Times"
  ) +
  scale_fill_grey(start = 0.2,
                  end = 0.8) +
  theme(text = element_text(size = 10, family="Times"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank()
  ) +
  ylab(expression(SST ~ (degree * C * ":" ~ centred))) +
  scale_x_continuous(expand = c(0.015, 0.015), breaks = pretty_breaks(n = 6))
STC


# CALCITE LAYER 
nucella_thickness_o1 <- read.csv("https://raw.githubusercontent.com/dm807cam/nucella2021belgium/main/data/shell/thickness/nucella_thickness_o1a.csv", sep = ",", header = T)
glimpse(nucella_thickness_o1)

nucella_thickness_o2 <- read.csv("https://raw.githubusercontent.com/dm807cam/nucella2021belgium/main/data/shell/thickness/nucella_thickness_o2a.csv", sep = ",", header = T)
glimpse(nucella_thickness_o2)

nucella_thickness_o3 <- read.csv("https://raw.githubusercontent.com/dm807cam/nucella2021belgium/main/data/shell/thickness/nucella_thickness_o3a.csv", sep = ",", header = T)
glimpse(nucella_thickness_o3)


nucella_storage <- read_delim("https://raw.githubusercontent.com/dm807cam/nucella2021belgium/main/data/shell/sample_storage.csv")


# Assign session name to shell thickness data
nucella_thickness_o1$session <- "o1"
nucella_thickness_o2$session <- "o2"
nucella_thickness_o3$session <- "o3"

# Combine shell layer thickness into single data frame
nucella_thickness <- na.omit(rbind(
  nucella_thickness_o1,
  nucella_thickness_o2,
  nucella_thickness_o3))

# Convert shell layer thickness from um to mm
nucella_thickness$aragonite_thickness_mm <- nucella_thickness$aragonite_thickness_um / 1000
nucella_thickness$calcite_thickness_mm <- nucella_thickness$calcite_thickness_um / 1000

# Combine storage conditions data frame with shell layer thickness data frame
nucella_thickness <- full_join(nucella_thickness, nucella_storage, by="sample")


nucella_shape <- read.csv("https://raw.githubusercontent.com/dm807cam/nucella2021belgium/main/data/shell/shell_shape.csv", sep= ",", header = T)
nucella_shape$location <- as.factor(nucella_shape$location)

# GAM: Calcite Layer Thickness --------------------------------------------
m_ct_I <- gam(calcite_thickness_mm ~ 
                s(year, k=3, bs="tp") + 
                s(shell_height, k=3, bs="tp") + 
                s(location, bs="re"), 
              family = gaussian(link = "identity"), 
              data = nucella_shape)





par(mfrow=c(2,2), mai=c(.5,.5,.5,.5))
DHARMa::plotResiduals(m_ct_I)
DHARMa::plotQQunif(m_ct_I)
acf(resid(m_ct_I), lag.max = 36, main = "ACF")
pacf(resid(m_ct_I), lag.max = 36, main = "pACF")
dev.off()

gratia::appraise(m_ct_I)
# Model diagnostic looks unsatisfying. Changing the model family link
# does not improve the diagnostic plots. log-transform dependent variable.
nucella_shape$calcite_thickness_mm_log <- log(nucella_shape$calcite_thickness_mm)

m_ct_null <- gam(calcite_thickness_mm_log ~ 1,
                 family = gaussian(link = "identity"), 
                 data = nucella_shape)

m_ct_II <- gam(calcite_thickness_mm_log ~ 
                 s(year, k=4, bs="tp") + 
                 s(shell_height, k=4, bs="tp") + 
                 s(location, bs="re"), 
               family = gaussian(link = "identity"), 
               data = nucella_shape)

par(mfrow=c(2,2), mai=c(.5,.5,.5,.5))
DHARMa::plotResiduals(m_ct_II)
DHARMa::plotQQunif(m_ct_II)
acf(resid(m_ct_II), lag.max = 36, main = "ACF")
pacf(resid(m_ct_II), lag.max = 36, main = "pACF")
dev.off()

gratia::appraise(m_ct_II)
gratia::draw(m_ct_II)
summary(m_ct_II)



# Predict calcite layer thickness -----------------------------------------
pdat_calcite <- with(nucella_shape,
                     data.frame(
                       year = seq(min(year), max(year), length = 1000),
                       shell_height = median(nucella_shape$shell_height),
                       location = "ZEE"
                     ))

pred_calcite <- predict(
  m_ct_II,
  newdata = pdat_calcite,
  exclude = "location",
  type = "response",
  se.fit = T
)

predframe_calcite <- data.frame(
  year = pdat_calcite$year,
  pred_calcite = pred_calcite$fit,
  se_calcite = pred_calcite$se.fit
)


#### Edits made to their code and my plotting #####
pdatAT <- pdatAT%>% filter(year >= 1888)
pdatST <- pdatST%>% filter(year >= 1888)


### Plotting ###
plot_malia_original <- ggplot() + 
  geom_line(predframe_calcite, mapping = (aes(year, exp(pred_calcite))),
            colour = "black", linetype = 2) +
  geom_ribbon(predframe_calcite, mapping = aes(year, ymin = (exp(pred_calcite) - se_calcite), 
                                               ymax = (exp(pred_calcite) + se_calcite)),
              fill = "black",
              alpha = .15) +
  geom_line(pdatST, mapping = (aes(year, p2)), linetype = 2, colour = "blue") +
  geom_ribbon(pdatST, mapping = (aes(year, ymin = lower, ymax = upper)),
                alpha = .15, 
              fill = "blue") +
  geom_line(pdatAT, mapping = (aes(year, p2)), linetype = 2, colour = "red") +
  geom_ribbon(pdatAT, mapping = (aes(year, ymin = lower, ymax = upper)),
              fill = "red",
                alpha = .15) +
  ggtitle("Calcite Thickness and Environmental Temperature Trends") +
  theme(text = element_text(size = 10, family="Arial"),
        axis.title.y.left = element_text(colour = "black"),
        axis.title.y.right = element_text(colour = "blue"),
        legend.position = "none",
        plot.margin = unit(c(0.3, 0.3, 0, 0.3), "cm")) +
  scale_y_continuous(name = expression(paste(Calcite ~ layer ~ "(" * mm * ")")),
                     sec.axis = sec_axis(trans = ~., 
                                         name = expression("Temperature (celcius)"))) + 
  scale_x_continuous(name = expression("Year"), expand = c(0.015, 0.015), breaks = pretty_breaks(n = 6))


plot_malia_original
