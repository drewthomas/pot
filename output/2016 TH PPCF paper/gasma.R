library(nortest)

KB <- 1.38065e-23
ME <- 9.109383e-31
QE <- 1.6021766e-19

mime <- 1836.153
Temp <- 220
a <- 2.5e-6
R <- 4e-4
N_e <- 75000

rea <- function(partial_path)
{
	cat("Reading", partial_path, "...\n")
	d <- read.table(paste("./", partial_path, sep=""), header=TRUE)

	if (grepl("macro", partial_path)) {
		# Downsample macroscopic data time series by 2/3 after the first
		# in-simulation microsecond so the plot PDF files are smaller.
		n <- dim(d)[1]
		indices <- (1:n)[d$TIME < 1e6]
		indices <- unique(c(indices, seq(1, n, 3)))
		d <- d[indices,]
	}
	
	return(d)
}

his <- function(x, ...)
{
	co <- "#696969"
	hist(x, freq=FALSE, breaks=49, col=co, border=co, main="", ylab="",
	     cex.lab=1.6, cex.axis=1.4, ...)
	grid(col="grey")
}

his_vel_comp <- function(mass_ratio=1, v, x_la, ...)
{
	vth <- sqrt(KB * Temp / (mass_ratio * ME))
	his(v, xlab=x_la)
	curve(dnorm(x, 0, vth / 1e3), lwd=3, lty="dashed", add=TRUE)
	if (mass_ratio == 1) {
		main_label <- expression(e^-phantom(0))
	} else {
		main_label <- expression(p^+phantom(0))
	}
	text(quantile(v, 0.995), 0.27 / (vth / 1e3), main_label, cex=2)
}

maxwell_speed_pdf <- function(mass_ratio=1, x)
{
	a <- sqrt(KB * Temp / (mass_ratio * ME))
	return(sqrt(2/pi) * (x^2) * exp(-(x^2) / (2*a*a)) / (a^3))
}

plot_v <- function(d, species="e", path="")
{
	d <- d[d$S == species,]
	if (species == "e") {
		mass_ratio <- 1
	} else {
		mass_ratio <- mime
	}

	if (path != "") {
		pdf(path, width=6, height=7)
	}

	par(las=1, mar=c(4, 4.4, 0.4, 0.4), mfrow=c(2,2))

	his_vel_comp(mass_ratio, d$VX / 1e3, expression(v[x] ~ (km %.% s^-1)))
	his_vel_comp(mass_ratio, d$VY / 1e3, expression(v[y] ~ (km %.% s^-1)))
	his_vel_comp(mass_ratio, d$VZ / 1e3, expression(v[z] ~ (km %.% s^-1)))

	speeds <- sqrt(d$VX^2 + d$VY^2 + d$VZ^2) / 1e3
	his(speeds, xlab=expression("speed " ~ (km %.% s^-1)))
	curve(1e3 * maxwell_speed_pdf(mass_ratio, 1e3*x),
	      add=TRUE, lwd=3, lty="dashed")
	if (mass_ratio == 1) {
		main_label <- expression(e^-phantom(0))
	} else {
		main_label <- expression(p^+phantom(0))
	}
	text(quantile(speeds, 0.997), 0.6 / mean(speeds), main_label, cex=2)

	if (path != "") {
		dev.off()
	} else {
		par(mfrow=c(1,1))
	}

	print(signif(c(ad.test(d$VX)$p.value,
	               ad.test(d$VY)$p.value,
	               ad.test(d$VZ)$p.value), 2))
}

oml_q_curve <- function(mac, n_e=2.797645e14, settle=5e-7)
{
	aleph <- 4 * pi * (a^2) * n_e * sqrt(KB * Temp / (2 * pi * ME))
	mac <- mac[mac$TIME >= settle,]
	q <- data.frame(TIME=mac$TIME, Q=rep(0, length(mac$TIME)))

	for (i in 2:length(mac$TIME)) {
		ro <- mac[i,]
		eta <- q$Q[i-1] * (QE^2) / (4 * pi * 8.85e-12 * KB * Temp * a)
		dq <- aleph * (((sqrt(1 / mime)) * (1 - eta)) - exp(eta))
		q$Q[i] <- q$Q[i-1] + ((mac$TIME[i] - mac$TIME[i-1]) * dq)
	}

	lines(1e6 * q$TIME, q$Q, lty="dashed", lwd=3)
}

dnpns <- rea("np-no-sph.dat")
plot_v(dnpns, "e", "v-e-pot-np-ns.pdf")
plot_v(dnpns, "i", "v-i-pot-np-ns.pdf")

dnpws_mac <- rea("np-with-sph-macro.dat")
dnpws_mac_sub <- dnpws_mac[dnpws_mac$TIME > 1e-6,]

pdf("q-pot-np-ws.pdf", width=5, height=5)
par(las=1, mar=c(4.2, 4.2, 0.1, 0.1))
plot(1e6 * dnpws_mac$TIME[dnpws_mac$TIME >= 5e-7],
     dnpws_mac$DGQ[dnpws_mac$TIME >= 5e-7], type="n",
     xlab=expression("time " * (mu * s)), ylab="", cex.axis=1.3, cex.lab=1.3)
mtext("q", 2, 3, cex=1.4)
mu_sig <- c(-82.6761, sqrt(25.6239))
rect(0.5, mu_sig[1] - mu_sig[2], 6, mu_sig[1] + mu_sig[2], col="lightgrey",
     border=NA)
grid(col="grey")
lines(1e6 * dnpws_mac$TIME[dnpws_mac$TIME >= 5e-7],
      dnpws_mac$DGQ[dnpws_mac$TIME >= 5e-7])
oml_q_curve(dnpws_mac)
legend(2.05, -22, c("q from simulation", "OML theory prediction"),
       col="black", lty=c("solid", "dashed"), lwd=c(1,3), bg="white",
       cex=1.05)
dev.off()

dns_mac <- rea("no-sph-1-macro.dat")

dns_mac_mini <- dns_mac[dns_mac$TIME < 5e-7,]
dns_mac_mini$TIME <- 1e9 * dns_mac_mini$TIME
dns_mac_mini$PE <- dns_mac_mini$PE / (KB * Temp)
dns_mac_mini$PI <- dns_mac_mini$PI / (KB * Temp)
dns_mac_mini_2 <- dns_mac_mini[dns_mac_mini$TIME <= 30,]

fit_osc <- TRUE

if (fit_osc) {

	omega_pe <- 2 * pi / 6.65874  # period of 6.65874 in nanoseconds, notice
	omega_pi <- omega_pe / sqrt(1836.15)
	nl_mod_1 <- nls(PE ~ offset + (TIME / tau_s)^2
	                     + (V_a * exp(-TIME / tau_d)
	                        * cos((sqrt((omega_pe^2) + (tau_d^-2)) * TIME)
    	                          + phase)),
	                dns_mac_mini_2,
    	            start=list(offset=-0.2, tau_s=99, V_a=0.2,
        	                   tau_d=2, phase=0))
	nl_mod_2 <- nls(PE ~ offset
	                     + (V_a * exp(-TIME / tau_d)
	                        * cos((sqrt((omega_pi^2) + (tau_d^-2)) * TIME)
	                              + phase)),
    	            dns_mac_mini[dns_mac_mini$TIME > 8,],
	                start=list(offset=0, V_a=-0.2, tau_d=50, phase=0))
	#print(summary(nl_mod_1))
	#print(summary(nl_mod_2))

	tau_ds <- as.vector(c(nl_mod_1$m$getPars()["tau_d"],
	                      nl_mod_2$m$getPars()["tau_d"]))
	omegas <- as.vector(c(sqrt((omega_pe^2) + (tau_ds[1]^-2)),
	                      sqrt((omega_pi^2) + (tau_ds[2]^-2))))

	## Print both the oscillation frequencies (`omegas`) and their standard
	## errors, as estimated by the delta method.
	#print(signif(omegas, 3))
	#print(1 / ((tau_ds^3) * omegas))

	nl_mod_3 <- nls(PE ~ offset + (TIME / tau_s)^2
	                     + (V_a * exp(-TIME / tau_d)
	                        * cos((omega * TIME) + phase)),
	                dns_mac_mini_2,
	                start=list(offset=-0.2, tau_s=99, V_a=0.2,
	                           tau_d=2, omega=1, phase=0))
	nl_mod_4 <- nls(PE ~ offset
	                     + (V_a * exp(-TIME / tau_d)
	                        * cos((omega * TIME) + phase)),
	                dns_mac_mini[dns_mac_mini$TIME > 8,],
	                start=list(offset=0, V_a=0.2, tau_d=50,
	                           omega=omega_pi, phase=pi))
	print(summary(nl_mod_3))
	print(summary(nl_mod_4))

}

pdf("wobble-pot-ns.pdf", width=8, height=5)
par(las=1, mar=c(4.4, 4.5, 0.1, 0.2), mfrow=c(1,2))
plot(dns_mac_mini_2$TIME, dns_mac_mini_2$PE, type="n", ylim=c(-0.26, -0.02),
     xlab="time (ns)", ylab=expression("mean electron PE " * (k[B] * T)))
grid(col="grey")
if (fit_osc) {
	lines(dns_mac_mini_2$TIME, predict(nl_mod_3, dns_mac_mini_2), lwd=4,
	      col="darkgrey")
}
lines(dns_mac_mini_2$TIME, dns_mac_mini_2$PE)
plot(dns_mac_mini$TIME / 1e3, dns_mac_mini$PE, type="n",
     xlab=expression("time " * (mu * s)),
     ylab=expression("mean electron PE " * (k[B] * T)))
grid(col="grey")
if (fit_osc) {
	lines(dns_mac_mini$TIME / 1e3, predict(nl_mod_4, dns_mac_mini), lwd=4,
	      col="darkgrey")
}
lines(dns_mac_mini$TIME / 1e3, dns_mac_mini$PE)
dev.off()
