QE <- 1.602e-19
ME <- 9.11e-31
KB <- 1.38e-23
EPSILON0 <- 8.85419e-12

#descrip_file <- "output/2014\\ to\\ 2015\\ CX1/np-with-sph/A/pot-run-descrip.txt"
#macro_file <- "output/2014\\ to\\ 2015\\ CX1/np-with-sph/A/pot-macroscopic.dat"
#descrip_file <- "output/2014\\ to\\ 2015\\ CX1/B-5/A/pot-run-descrip.txt"
#macro_file <- "output/2014\\ to\\ 2015\\ CX1/B-5/pot-macro-all.dat"
#descrip_file <- "output/2014-06\\ to\\ 2014-08/CX1\\ pot\\ B0/J/pot-run-descrip.txt"
#macro_file <- "output/2014-06\\ to\\ 2014-08/CX1\\ pot\\ B0/pot-macro-all.dat"
descrip_file <- "output/pot-run-descrip.txt"
macro_file <- "output/pot-macroscopic.dat"
#descrip_file <- "~/Downloads/pot-run-descrip.txt"
#macro_file <- "~/Downloads/pot-macroscopic.dat"

do_omlchafitting <- TRUE

num_from_grepped_cmd <- function(string, cut_arg="-c28-")
{
	cmd <- paste("grep", string, descrip_file, "| tail -n 1 | cut", cut_arg)
	return(as.numeric(system(cmd, intern=TRUE))[1])
}

dust_r <- num_from_grepped_cmd("'n radius (m)'", "-c23-")
pixel_size <- num_from_grepped_cmd("pixel")
size_in_pixels <- num_from_grepped_cmd("'n radius (p'")
drift <- num_from_grepped_cmd("'ift vel'")
if (is.na(size_in_pixels)) {
	# domain's cubical, not radial
	size_in_pixels <- num_from_grepped_cmd("'n length (p'") / 2
}
n <- num_from_grepped_cmd("density")
settle_time <- num_from_grepped_cmd("settling")
if (is.na(settle_time)) {
	settle_time <- num_from_grepped_cmd("'warm\\-*up'")
}

d <- read.table(gsub("\\\\", "", macro_file), header=TRUE, nrows=1e6)

if (settle_time < min(d$TIME)) {
	settle_time <- min(d$TIME)
}

if (do_omlchafitting) {
	post_pipe <- "| tail -n 2 | head -n 1 | sed 's/^ *//'"
	if (drift < 1e-9) {
		omlchafit <- system(paste("chafitoml/chafitoml", dust_r, n/2, settle_time, macro_file, post_pipe), intern=TRUE)
	} else {
		omlchafit <- system(paste("chafitoml/chafitsoml", dust_r, n/2, drift, settle_time, macro_file, post_pipe), intern=TRUE)
	}
	omlchafit <- c(sapply(strsplit(omlchafit, " +"), as.numeric))
}

# decide how much to smooth the really noisy time series
if (length(d[,1]) < 210) {
	smoo <- 1
} else if (length(d[,1]) < 400) {
	smoo <- 3
} else {
	smoo <- 7
}

library(LambertW)  # no longer helpfully defines `erf`

erf <- function(x)
{
	return(2 * pnorm(sqrt(2) * x) - 1)
}

s1 <- function(v)
{
	v[v==0] <- 1e-9
	return(exp(-v^2) / 2 + sqrt(pi) * erf(v) * (1 + 2*v^2) / (4*v))
}

s2 <- function(v)
{
	v[v==0] <- 1e-9
	return(sqrt(pi) * erf(v) / (2*v))
}

draw_expected_soml <- function(Q_not_eta=FALSE)
{
	mu <- sqrt(1836.15)
	dsub <- d[d$TIME >= settle_time,]
	upsilon <- QE^2 * dust_r * (n/2) / (EPSILON0 * sqrt(2*pi*ME*KB * dsub$TE[1]))
	y <- rep(NA, length(dsub$TIME))
	y[1] <- (dsub$TIME[1] - settle_time) * upsilon *
	        (sqrt(dsub$TI[1] / dsub$TE[1]) / mu - 1)
	for (i in 2:length(dsub$TIME)) {
		upsilon <- QE^2 * dust_r * (n/2) /
		           (EPSILON0 * sqrt(2*pi*ME*KB * dsub$TE[i]))
		beta <- dsub$TI[i] / dsub$TE[i]
		v <- drift / (sqrt(2 * KB * dsub$TI[i] / ME) / mu)
		y[i] <- y[i-1] + (dsub$TIME[i] - dsub$TIME[i-1]) * upsilon *
		                 (sqrt(beta) * (s1(v) - s2(v) * y[i-1] / beta) / mu - exp(y[i-1]))
	}
	if (Q_not_eta) {
		lines(dsub$TIME,
		      runmed((4*pi*EPSILON0*KB/QE) * dsub$TE * dust_r * y / QE, smoo),
		      col="green")
	} else {
		lines(dsub$TIME, y, col="green")
	}
}

draw_omlchafit <- function(Q_not_eta=FALSE)
{
	old_dust_r <- dust_r
	old_n <- n
	dust_r <- omlchafit[1]
	n <- 2 * omlchafit[2]
	mu <- sqrt(omlchafit[3])
	dsub <- d[d$TIME >= settle_time,]
	upsilon <- QE^2 * dust_r * (n/2) / (EPSILON0 * sqrt(2*pi*ME*KB * dsub$TE[1]))
	y <- rep(NA, length(dsub$TIME))
	y[1] <- (dsub$TIME[1] - settle_time) * upsilon *
	        (sqrt(dsub$TI[1] / dsub$TE[1]) / mu - 1)
	for (i in 2:length(dsub$TIME)) {
		upsilon <- QE^2 * dust_r * (n/2) /
		           (EPSILON0 * sqrt(2*pi*ME*KB * dsub$TE[i]))
		beta <- dsub$TI[i] / dsub$TE[i]
		if (drift < 1e-9) {
			v <- 0
		} else {
			v <- omlchafit[4] / sqrt(2 * KB * dsub$TI[i] / (omlchafit[3] * ME))
		}
		y[i] <- y[i-1] + (dsub$TIME[i] - dsub$TIME[i-1]) * upsilon *
		                 (sqrt(beta) * (s1(v) - s2(v) * y[i-1] / beta) / mu - exp(y[i-1]))
	}
	if (Q_not_eta) {
		lines(dsub$TIME,
		      runmed((4*pi*EPSILON0*KB/QE) * dsub$TE * dust_r * y / QE, smoo),
		      col="orange")
	} else {
		lines(dsub$TIME, y * dust_r / old_dust_r, col="orange")
	}
	dust_r <- old_dust_r
	n <- old_n
}

latest_TI <- median(tail(d$TI[!is.na(d$TI)], 0.1 * length(d$TI)), na.rm=TRUE)
latest_TE <- median(tail(d$TE[!is.na(d$TE)], 0.1 * length(d$TE)), na.rm=TRUE)
beta <- latest_TI / latest_TE
expec_v <- drift / sqrt(2 * KB * latest_TI / (1836.153 * ME))
srat <- s1(expec_v) / s2(expec_v)

expected_eta <- (srat * beta) - W(sqrt(1836.153 * beta) * exp(srat * beta) / s2(expec_v))

dust_chgs_to_potl <- function(chgs)
{
	return(chgs * QE / (4 * pi * EPSILON0 * dust_r))
}

# floating potential, normalized
d$SURFPOT <- dust_chgs_to_potl(d$DGQ)
d$PFN <- QE * d$SURFPOT / (KB * d$TE)

# try fitting an exponential decay curve to grain charge time evolution
expected_phif <- expected_eta * KB * median(tail(d$TE, 10)) / QE
expected_qf <- 4 * pi * EPSILON0 * dust_r * expected_phif
fit <- NA
fit_pars <- NA
try({
	fit <- nls(DGQ ~ (TIME >= START) * ASYM * (1 - exp(-(TIME-START)/TAU)), d,
	           start=list(ASYM=expected_qf / QE,
	                      TAU=median(d$TIME),
	                      START=min(c(0, d$TIME[d$DGQ < 0]))))
	fit_pars <- fit$m$getPars()
	fit_ses <- summary(fit)[["parameters"]][, "Std. Error"]
	# apply quick & dirty autocorrel. correction based on
	# Foster & Rahmstorf (2011): compute the sample lag-1 residual
	# autocorrelation and compute a std. err. correction factor
	autocorr_q <- acf(residuals(fit)[d$TIME >= settle_time],
	                  plot=FALSE)$acf[2]
	nu_q <- (1 + autocorr_q) / (1 - autocorr_q)
	fit_ses <- fit_ses * sqrt(nu_q)
}, silent=TRUE)

# also compute a quick & dirty autocorrel. correction factor for eta
# using the lag-1 autocorrelation
autocorr_pfn <- acf(d$PFN[d$TIME >= settle_time],
                    plot=FALSE, na.action=na.omit)$acf[2]
nu_pfn <- (1 + autocorr_pfn) / (1 - autocorr_pfn)

omlchafit_eqm <- function(cutoff)
{
	r <- omlchafit[1]
	n_e <- omlchafit[2]
	ds <- d[d$TIME > cutoff,]
	if (drift < 1e-9) {
		v <- 0
	} else {
		v <- omlchafit[4] / sqrt(2 * KB * ds$TI / (omlchafit[3] * ME))
	}
	beta <- median(ds$TI / ds$TE)
	srbe <- median((s1(v) / s2(v)) * ds$TI / ds$TE)
	return((srbe - W(sqrt(omlchafit[3] * beta) * exp(srbe) / median(s2(v))))
	       * r / dust_r)
}

# optional extra, useful for estimating equilibrium eta from OML fit
summary_for_comp <- function(cutoff)
{
	ds <- d[d$TIME > cutoff,]
	print(round(quantile(ds$TI / ds$TE, c(0.05, 0.5, 0.95)), 4))
	print(round(quantile(ds$PFN, c(0.05, 0.5, 0.95)), 3))
	print(round(omlchafit_eqm(cutoff), 3))
	return(c(median(pixel_size * size_in_pixels / ds$DL),
	         median(ds$DGQ),
	         min(d$TIME[d$DGQ < (1 - exp(-1)) * median(ds$DGQ)]) - settle_time,
	         2 * pi * sqrt(ME * EPSILON0 / (n/2)) / QE,
	         drift / median(sqrt(KB * ds$TE / (1836.153 * ME)))))
}

par(mfrow=c(2,3), mar=c(5,5,3,2), cex.axis=1.5, cex.lab=1.5, cex.main=1.5)

smoothed_TE <- runmed(d$TE[!is.na(d$TE)], smoo)
plot(d$TIME[!is.na(d$TE)], smoothed_TE,
     type="l", col="blue",
     ylim=quantile(c(runmed(d$TI, smoo), smoothed_TE), c(0, 0.999)),
     xlab="time (s)", ylab="temperature (K)",
     main=expression(paste(T[i], " and ", T[e], " (smoothed)")))
abline(h=latest_TE, col="blue", lty="dashed")
abline(h=latest_TI, col="red", lty="dashed")
grid()
lines(d$TIME, runmed(d$TI, smoo), col="red")
#text(5e-6, 219, "B = 20 T", cex=2)

plot(d$TIME[!is.na(d$TE)], runmed((d$TI / d$TE)[!is.na(d$TE)], smoo),
     type="l", xlab="time (s)", ylab=expression(T[i]/T[e]),
     main=expression(paste(T[i] / T[e], " (smoothed)")))
abline(h=1, lty="dotted", col="grey")
abline(h=latest_TI/latest_TE, lty="dashed")

plot(d$TIME, d$KE/QE, type="l",
     ylim=quantile(c(d$KE, d$KI, d$PE, d$PI)/QE, c(1e-4, 0.9999), na.rm=TRUE),
     col="blue", lwd=2, lty="dashed",
     main="mean electronic & ionic KE & PE",
     xlab="time (s)", ylab="energy (eV)")
grid()
abline(h=0, lty="dotted")
lines(d$TIME, d$PE/QE, col="#0000ff9a")
lines(d$TIME, d$KI/QE, col="red", lwd=2, lty="dashed")
lines(d$TIME, d$PI/QE, col="#ff00007a")
#if (length(d$PE) > 500) {
#	lines(d$TIME, predict(loess(PE/QE ~ TIME, d, span=0.04)),
#	      col="#0000ff", lwd=2)
#	lines(d$TIME, predict(loess(PI/QE ~ TIME, d, span=0.04)),
#	      col="#990000", lwd=2)
#}

alpha <- QE^2 / (4*pi*EPSILON0*KB * latest_TE * dust_r)
stoc_moms <- system(paste("../OML\\ statistics/moments",
                          alpha, "0 42.85", latest_TI/latest_TE, expec_v),
                    intern=TRUE)
stoc_moms <- as.numeric(unlist(strsplit(stoc_moms, " ")))
smsd <- sqrt(stoc_moms[2])

plot(d$TIME, d$DGQ, type="n", main="charges on grain",
     xlab="time (s)", ylab="elementary charges")
rect(settle_time, stoc_moms[1] - (2 * smsd),
     max(d$TIME), stoc_moms[1] + (2 * smsd),
     col="#00000010", border=NA)
rect(settle_time, stoc_moms[1] - smsd,
     max(d$TIME), stoc_moms[1] + smsd,
     col="#00000010", border=NA)
text(settle_time, stoc_moms[1] + (0.5 * smsd),
     expression(phantom(0) %+-% 1 * sigma))
text(settle_time, stoc_moms[1] + (1.5 * smsd),
     expression(phantom(0) %+-% 2 * sigma))
if (do_omlchafitting) {
	draw_expected_soml(TRUE)
	draw_omlchafit(TRUE)
}
grid()
if ((length(fit) > 1) && (length(fit_pars) == 3)) {
	curve((x >= fit_pars[3]) * fit_pars[1]
	      * (1 - exp(-(x-fit_pars[3]) / fit_pars[2])),
	      add=TRUE, col="red", lty="dashed", n=501)
	text(0.5 * max(d$TIME), max(d$DGQ) + min(d$DGQ) / 9,
	     bquote(paste(tau == (.(signif(fit_pars[2], 3))
	                  %+-% .(signif(fit_ses[2], 2))), "s")),
	     cex=1.5, col="red")
	text(0.5 * max(d$TIME), max(d$DGQ) + min(d$DGQ) / 5,
	     bquote(paste(Q == (.(signif(fit_pars[1], 4))
	                  %+-% .(signif(fit_ses[1], 2))), "e")),
	     cex=1.5, col="red")
}
lines(d$TIME, d$DGQ)  # plot the actual data last for clarity!

if (max(abs(d$PFN))) {
	y_limits <- range(d$PFN, na.rm=TRUE)
	y_limits[1] <- min(c(expected_eta, y_limits[1]))
	plot(d$TIME, runmed(d$PFN, smoo), type="l", ylim=y_limits,
	     main=expression(paste("grain surface ", eta, " (smoothed)")),
	     xlab="time (s)")
	if (do_omlchafitting) {
		draw_expected_soml()
		draw_omlchafit()
	}
	grid()
	abline(h=expected_eta, lty="dashed", col="blue")
	text(max(d$TIME) / 6, 0.12 + expected_eta, signif(expected_eta, 4),
	     cex=1.4, col="blue")
} else {
	acf(d$TE, lag.max=min(c(100, length(d$TE))),
	    main=expression(paste(T[e], " autocorrelation")))
}

plot(d$TIME, runmed((pixel_size * size_in_pixels) / d$DL, smoo), type="l",
     main=expression(paste(R / lambda[D], " (smoothed)")),
     xlab="time (s)", ylab=expression(paste(R / lambda[D])))
grid()

#if (sd(d$TE) > 0) {
#	acf(d$TE, lag.max=150, main="electron temperature autocorrelation",
#	    na.action=na.pass)
#	grid()
#}

par(mfrow=c(1,1))
