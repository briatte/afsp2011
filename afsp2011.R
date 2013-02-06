# What: Replication code for "Quand la biopolitique..." (AFSP 2011) 
# Who:  François Briatte
# When: 2011-08-15

# Dependencies.
require(ggplot2)
require(MASS)
require(mgcv)
require(nlme)
require(plyr)
require(RColorBrewer)

setwd("~/Documents/Research/Cancer/Papers/2011-08-AFSP/afsp2011/")

# ========
# = DATA =
# ========

# Health expenditure data
# Source: OECD Health 2011
data <- read.csv("data/oecd_hexp.csv",sep="\t")
data <- subset(data, year < 2010)

eums <- c("Austria", "Belgium", "Denmark", "Finland", "France", "Germany", "Greece", "Ireland", "Italy", "Luxembourg", "Netherlands", "Portugal", "Spain", "Sweden")
eu14 <- subset(data, cty %in% eums)
eu15 <- subset(data, cty %in% c(eums, "United Kingdom"))
uk <- subset(data, cty == "United Kingdom")

# Collapsed data
eu14.summary <- ddply(eu14, 
                      .(year), 
                      summarise, 
                      sum.gdp  = sum(hexp_gdp,na.rm = TRUE),
                      mean.gdp = mean(hexp_gdp,na.rm = TRUE), 
                      sum.ppp  = sum(hexp_ppp,na.rm = TRUE), 
                      mean.ppp = mean(hexp_ppp,na.rm = TRUE)
                      )

eu15.summary <- ddply(eu15,
                      .(year),
                      summarise, 
                      sum.gdp  = sum(hexp_gdp,na.rm = TRUE),
                      mean.gdp = mean(hexp_gdp,na.rm = TRUE), 
                      sum.ppp  = sum(hexp_ppp,na.rm = TRUE), 
                      mean.ppp = mean(hexp_ppp,na.rm = TRUE)
                      )

yrange <- range(uk$year)

# Definition of UK government party
gov <- read.csv("data/uk_gov.csv")

# =========
# = MODEL =
# =========

# Custom functions for gamm() diagnostics plots and derivatives of fitted splines
source("functions/tsa.R")
source("functions/gamm.R")

# Lag plots
lag.plot1(ts(uk$hexp_gdp),12)

# Model 1
# GAMM with 15 years fitted smoother
#
m1 <- gamm(hexp_gdp ~ s(year, k = 15), data = uk); summary(m1$gam)
m1.resid <- resid(m1$lme, type = "normalized")

# Diagnostics
acf2(m1.resid)

# Model 2
# AR(1)
#
m2 <- gamm(hexp_gdp ~ s(year, k = 15), data = uk, correlation = corARMA(form = ~ year, p = 1))

# Model 3
# AR(2)
#
m3 <- gamm(hexp_gdp ~ s(year, k = 15), data = uk, correlation = corARMA(form = ~ year, p = 2))

# Diagnostics
# Notes:
# - AIC and BIC both select Model 3 as the best overall fit.
# - The distribution of residuals for Model 3 is leptokurtic.
anova(m1$lme, m2$lme, m3$lme)
plot(m3$gam, residuals = TRUE, pch = 19, cex = 0.75)
with(eu14.summary, tsDiagGamm(m3, timevar = year, observed = uk$hexp_gdp))

# Predicted values
# Notes:
# - Models 1 and 3 respectively provide the GAMM and AR(2) fitted smoothers.
# - See Fig. 10 for a comparison of error terms in both models.
pred <- with(eu14.summary, data.frame(year = seq(min(year), max(year), length = 200)))
pred.m1 <- predict(m1$gam, newdata = pred)
pred.m3 <- predict(m3$gam, newdata = pred)

# First derivatives
# Notes:
# - Computations are based on the method of finite differences.
# - See Fig. 11 for the 95% and 99% confidence intervals.
m3.d <- Deriv(m3, n = 200)

# Confidence intervals
# Notes:
# - 99% confidence intervals are used to identify significant changes.
# - See Fig. 12 for a summary of model results over the actual data.
m3.d.95ci <- confint(m3.d, alpha = 0.05)
m3.d.99ci <- confint(m3.d, alpha = 0.01)
m3.d.inc95 <- ifelse(m3.d.95ci$year$lower > 0, m3.d$year$deriv,NA)
m3.d.inc99 <- ifelse(m3.d.99ci$year$lower > 0, m3.d$year$deriv,NA)
sig <- signifD(pred.m3, m3.d$year$deriv, m3.d.99ci$year$upper, m3.d.99ci$year$lower, eval = 0)

# Simulated trends

# Graphical options
layout(matrix(c(1,2), 2, 1), widths = c(1,1), heights = c(1.25,1))

# Simulation from posterior distribution of beta
Rbeta <- mvrnorm(n = 1000, coef(m3$gam), vcov(m3$gam))
Xp <- predict(m3$gam, newdata = pred, type = "lpmatrix")
sim1 <- Xp %*% t(Rbeta)

# Simulation plots
par(mar = c(3,3,3,2) + 0.1)
set.seed(321)
want <- sample(1000, 25)
ylim <- range(sim1[,want], uk$hexp_gdp)
plot(hexp_gdp ~ year, data = uk, ylim = ylim, type = "n", xlab = "", ylab="")
matlines(pred$year, sim1[,want], col = "grey", lty = 1, pch = NA)
points(hexp_gdp ~ year, data = uk, col = "orange", bg = "yellow", pch = 21, cex = 1)
abline(v = 2000, lty = 2,col = "grey")
title("Simulated trends for Model 3")

# Post-1999 observations
par(mar = c(3,3,1,2) + 0.1)
set.seed(321)
want <- sample(1000, 50)
rwant <- with(pred, which(year >= 1999))
twant <- with(uk, which(year >= 1999))
ylim <- range(sim1[rwant,want], uk$hexp_gdp[twant])
plot(hexp_gdp ~ year, data = uk, ylim = ylim, xlim = c(1999, 2009), type = "n", xlab = "", ylab = "")
matlines(pred$year, sim1[,want], col = "grey", lty = 1, pch = NA)
points(hexp_gdp ~ year, data = uk, col = "orange", bg = "yellow", pch = 21, cex = 1.5)
abline(v = 2000,lty = 2, col = "grey")

# ==========
# = GRAPHS =
# ==========

# Graphical options
mark <- geom_vline(xintercept = 2000, linetype = "dashed", colour = "#999999")
line <- geom_hline(yintercept = 1, colour = "black", linetype = "solid")
cols <- brewer.pal(8, "Set1")

# Fig. 1
# Total health expenditure as % of GDP, 1960-2009
#
fig1.hexp_gdp <- ggplot() + 
	geom_boxplot(
    data = eu15, aes(x = year, y = hexp_gdp,group=year), 
    colour = "darkgrey", outlier.colour = "darkgrey"
    ) +
	geom_point(data = eu14.summary, aes(x = year, y = mean.gdp, size=mean.ppp, colour = "EU")) +
	geom_line(data = eu14.summary, aes(x = year, y = mean.gdp, colour = "EU")) +
	geom_point(data = uk, aes(x = year, y = hexp_gdp, size = hexp_ppp, colour = "UK")) +
	geom_line(data = uk, aes(x = year, y = hexp_gdp, colour = "UK")) +
	scale_size_area(name = "PPP") + scale_colour_manual("GDP", values = c("EU" = cols[2], "UK" = cols[1])) +
	xlab("") + ylab(""); fig1.hexp_gdp

# Fig. 2
# Total health expenditure as US$ PPP, 1960-2009
#
fig2.hexp_ppp <- ggplot() + 
	geom_boxplot(data = eu15, aes(x = year, y = hexp_ppp,group=year), colour = "darkgrey", outlier.colour = "darkgrey") +
	geom_point(data = eu14.summary, aes(x = year, y = mean.ppp, size=mean.gdp, colour = "EU")) +
	geom_line(data = eu14.summary, aes(x = year, y = mean.ppp, colour = "EU")) +
	geom_point(data = uk, aes(x = year, y = hexp_ppp, size = hexp_gdp, colour = "UK")) +
	geom_line(data = uk, aes(x = year, y = hexp_ppp, colour = "UK")) +
  scale_size_area("% PPP") +
  scale_colour_manual("Countries", values = c("EU" = cols[2], "UK" = cols[1])) +
	xlab("") + ylab(""); fig2.hexp_ppp

# Contribution to EU mean (as % of GDP)
eu15$contrib <- 0
for (y in yrange[1]:yrange[2]) {
	yeartotal <- eu15.summary$sum.gdp[eu15.summary$year == y]
	eu15$contrib[eu15$year == y] <- 100 * eu15$hexp_gdp[eu15$year == y] / yeartotal
}

# Fig. 3
# Distribution of mean(EU) and UK contributions
#
eu15$isuk <- ifelse(eu15$cty=="United Kingdom", "UK", "EU")
fig3.contrib_density <- ggplot(eu15) + 
	geom_density(
    aes(x = contrib, group = isuk, colour = isuk, fill= isuk), 
    alpha = .2, binwidth = .5, na.rm = TRUE
    ) + 
  scale_fill_manual("", values = c("EU" = cols[2], "UK" = cols[1])) + 
  scale_colour_manual("", values = c("EU" = cols[2], "UK" = cols[1])) +
	xlab("") + ylab(""); fig3.contrib_density

#
# (G1) RAW UK AND OECD SERIES
#

# Bounds for GDP areas
gdp.bounds <- range(eu14.summary$mean.gdp, na.rm = TRUE)
gdp.min <- min(gdp.bounds)
gdp.max <- max(gdp.bounds)

# Bounds for PPP areas
ppp.bounds <- range(eu14.summary$mean.ppp, na.rm = TRUE)
ppp.min <- min(ppp.bounds)
ppp.max <- max(ppp.bounds)

# Fig. 4
# Total health expenditure as % of GDP, 1960-2009
#
fig4.comp_gdp <- ggplot() +
	geom_line(data = eu14.summary, aes(x = year, y = mean.gdp, colour = "EU")) +
	geom_line(data = uk, aes(x = year, y = hexp_gdp, colour = "UK")) +
	scale_colour_manual("", values = c("EU" = "black", "UK" = "darkgrey")) +
	geom_rect(
    data = gov, 
    aes(NULL, NULL, xmin = start, xmax = end, fill = party), 
    ymin = gdp.min, 
    ymax = gdp.max, 
    alpha = .2
    ) +
	scale_fill_manual("", values = c("blue", "red")) +
	mark + xlab("") + ylab(""); fig4.comp_gdp

# Fig. 5
# Total health expenditure as US$ PPP, 1960-2009
#
fig5.comp_ppp <- ggplot() +
  geom_line(data = eu14.summary, aes(x = year, y = mean.ppp, colour = "EU")) +
  geom_line(data = uk, aes(x = year, y = hexp_ppp, colour = "UK")) +
  scale_colour_manual("", values = c("EU" = "black", "UK" = "darkgrey")) +
  geom_rect(
    data = gov, 
    aes(NULL, NULL, xmin = start, xmax = end, fill = party), 
    ymin = ppp.min, 
    ymax = ppp.max, 
    alpha = .2
  ) +
  scale_fill_manual("", values = c("blue", "red")) +
  mark + xlab("") + ylab(""); fig5.comp_ppp

#
# (G2) UK/OECD RATIOS
#

# Definition of ratios
ratio.gdp <- uk$hexp_gdp/eu14.summary$mean.gdp
ratio.ppp <- uk$hexp_ppp/eu14.summary$mean.ppp
ratio.gdp.a <- uk$hexp_gdp/eu15.summary$mean.gdp
ratio.ppp.a <- uk$hexp_ppp/eu15.summary$mean.ppp

# Bounds
rbounds <- c(range(ratio.gdp), range(ratio.ppp))
rmin <- min(rbounds)
rmax <- max(rbounds)
rbounds.a <- c(range(ratio.gdp.a), range(ratio.ppp.a))
rmin.a <- min(rbounds.a)
rmax.a <- max(rbounds.a)

# Bounds for GDP ratios
rbounds.gdp <- c(range(ratio.gdp), range(ratio.gdp.a))
rmin.gdp <- min(rbounds.gdp)
rmax.gdp <- max(rbounds.gdp)

# Bounds for PPP ratios
rbounds.ppp <- c(range(ratio.ppp), range(ratio.ppp.a))
rmin.ppp <- min(rbounds.ppp)
rmax.ppp <- max(rbounds.ppp)

# Fig. 6
# UK/EU health expenditure ratio as % of GDP
# 
fig6.ratios_gdp <- ggplot(uk) +
	geom_line(aes(x = year, y = ratio.gdp, colour = "EU-14")) +
	geom_line(aes(x = year, y = ratio.gdp.a, colour = "EU-15")) +
  scale_colour_manual("", values = c("EU-14" = "black","EU-15" = "darkgrey")) +
  geom_rect(
    data = gov, 
    aes(NULL, NULL, xmin = start, xmax = end, fill = party), 
    ymin = rmin.gdp, 
    ymax = rmax.gdp, 
    alpha = .2
  ) +
  scale_fill_manual("", values = c("blue", "red")) +
	line + mark + xlab("") + ylab(""); fig6.ratios_gdp

# Fig. 7
# OECD/UK health expenditure ratio as US$ PPP
# 
fig7.ratios_ppp <- ggplot(uk) +
	geom_line(aes(x = year, y = ratio.ppp, colour = "EU-14")) +
	geom_line(aes(x = year, y = ratio.ppp.a, colour = "EU-15")) +
  scale_colour_manual("", values = c("EU-14" = "black","EU-15" = "darkgrey")) +
  geom_rect(
    data = gov, 
    aes(NULL, NULL, xmin = start, xmax = end, fill = party), 
    ymin = rmin.ppp, 
    ymax = rmax.ppp, 
    alpha = .2
  ) +
  scale_fill_manual("", values = c("blue", "red")) +
	line + mark + xlab("") + ylab(""); fig7.ratios_ppp

# Fig. 8
# OECD/UK EU-14 ratios (both measures)
#
fig8.ratios <- ggplot(uk) +
	geom_line(aes(x = year, y = ratio.gdp, colour = "GDP")) +
	geom_line(aes(x = year, y = ratio.ppp, colour = "PPP")) +
	scale_colour_manual("", values = c("GDP" = "black", "PPP" = "darkgrey")) +
	geom_rect(data = gov, aes(NULL, NULL, xmin = start, xmax = end, fill = party), ymin = rmin, ymax = rmax, alpha = .2) +
	scale_fill_manual("", values = c("blue", "red")) +
	line + mark + xlab("") + ylab(""); fig8.ratios

# Fig. 9
# OECD/UK EU-15 ratios (both measures)
#
fig9.ratios_a <- ggplot(uk) +
	geom_line(aes(x = year, y = ratio.gdp.a, colour = "GDP")) +
	geom_line(aes(x = year, y = ratio.ppp.a, colour = "PPP")) +
	scale_colour_manual("", values = c("GDP" = "black", "PPP" = "darkgrey")) +
	geom_rect(data = gov, aes(NULL, NULL, xmin = start, xmax = end, fill = party), ymin = rmin.a, ymax = rmax.a, alpha = .2) +
  scale_fill_manual("", values = c("blue", "red")) +
	line + mark + xlab("") + ylab(""); fig9.ratios_a

#
# (G3) MODEL RESULTS
#

# Fig. 10
# Uncorrelated vs. AR(2) errors
#
fig10.fitted.smoothers <- ggplot() + 
	 geom_point(data = uk, aes(x = year, y = hexp_gdp)) +
	 geom_line(data = pred, aes(x = year, y = pred.m1, colour = "GAMM")) +
	 geom_line(data = pred, aes(x = year, y = pred.m3, colour = "AR(2)")) +
	 scale_colour_manual("", values = c("GAMM" = cols[1], "AR(2)" = cols[2])) +
	 xlab("") + ylab(""); fig10.fitted.smoothers

# Fig. 11
# First derivatives (finite differences)
#
fig11.fitted.derivatives <- ggplot() +
	geom_ribbon(aes(x = pred$year, ymin = m3.d.99ci$year$lower, ymax = m3.d.95ci$year$lower), fill = "red",alpha=0.2) +
	geom_ribbon(aes(x = pred$year, ymin = m3.d.95ci$year$upper, ymax = m3.d.99ci$year$upper), fill = "red",alpha=0.2) +
	geom_ribbon(aes(x = pred$year, ymin = m3.d.95ci$year$lower, ymax = m3.d.95ci$year$upper), fill = "orange",alpha=0.2) +
	geom_line(aes(x = pred$year, y = m3.d$year$deriv)) +
	geom_line(aes(x = pred$year, y = m3.d.99ci$year$upper), linetype = "dashed", colour = "red") +
	geom_line(aes(x = pred$year, y = m3.d.99ci$year$lower), linetype = "dashed", colour = "red") + 
	geom_line(aes(x = pred$year, y = m3.d.95ci$year$upper), linetype = "dashed", colour = "orange") +
	geom_line(aes(x = pred$year, y = m3.d.95ci$year$lower), linetype = "dashed", colour = "orange") +
	geom_line(aes(x = pred$year, y = m3.d.inc95), size = 2, colour = "orange") +
	geom_line(aes(x = pred$year, y = m3.d.inc99), size = 2, colour = "red") +
	geom_hline(yintercept = 0) + xlab("") + ylab(""); fig11.fitted.derivatives

# Fig. 12
# Fitted additive model with AR(2) errors and superimposed periods of significant increase in expenditure
# 
fig12.fitted.changes <- ggplot() +
	geom_line(data = pred, aes(x = year, y = pred.m3, colour = "Fitted")) +
	geom_line(data = pred, aes(x = year, y = sig$incr, colour = "Increase")) +
	geom_point(data = uk, aes(x = year, y = hexp_gdp, colour = "Actual")) +
	scale_colour_manual("", values = c("Actual" = "black", "Fitted" = "darkgrey", "Increase"="red")) +
	geom_rect(data = gov, aes(NULL, NULL, xmin = start, xmax = end, fill = party), ymin = min(uk$hexp_gdp), ymax = max(uk$hexp_gdp), alpha = .2) +
	scale_fill_manual("", values = c("blue", "red")) +
	mark + xlab("") + ylab(""); fig12.fitted.changes

# List figures
print(graphs <- ls(pattern="fig"))