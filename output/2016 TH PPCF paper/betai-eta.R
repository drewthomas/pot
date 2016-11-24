phl <- read.table("phl-2007-fig-11a.dat", header=TRUE)
mine <- read.table("cx1-betai-eta.dat", header=TRUE)

eta_of_beta_i <- function(beta_i, mu=sqrt(1836.153), theta=1)
{
	f <- function(eta)
	{
		iota_star <- function(beta)
		{
			coeffs <- c(1, -0.0946, -0.305, 0.950, -2.200, 1.150)
			z <- beta / (1 + beta)
			return(sum(coeffs * (z^(0:5))))
		}

		A_star <- function(w)
		{
			return((0.678 * w) + (1.543 * w^2) - (1.212 * w^3))
		}

		beta_e <- mu * sqrt(theta) * beta_i
		beta_e <- 1e-10 + beta_e  # prevent divisions by zero below
		w <- (-eta / beta_e) / (1 - (eta / beta_e))
		A_st <- A_star(w)
		e_denom <- (A_st + ((1 - A_st) * iota_star(beta_e))) * exp(eta)
		return((((sqrt(theta) / mu) * (1 - (eta / theta)) / e_denom) - 1)^2)
	}

	op_result <- optimize(f, c(-2.6, -1.6))
	return(op_result$minimum)
}

phl_blue <- phl[phl$TI == 1,]

beta_i_samp <- seq(0, 0.5, 2e-3)
oml_ion_eta <- sapply(beta_i_samp, eta_of_beta_i)

plo <- function(x_lim=NA, y_lim=NA,...)
{
	x_lab <- expression(beta[i] %==% symbol("\341") * a / r[gi]
	                    * symbol("\361"))

	if ((length(x_lim) == 1) && is.na(x_lim)) {
		x_lim <- range(c(phl_blue$BETAI, mine$BETAI))
	}
	if ((length(y_lim) == 1) && is.na(y_lim)) {
		y_lim=range(c(phl_blue$PHIF, mine$ETAMED))
	}

	plot(phl_blue$BETAI, phl_blue$PHIF, pch=5,
	     xlim=x_lim, ylim=y_lim,
	     cex=1.1, cex.axis=1.3, cex.lab=1.5, xlab=x_lab, ylab="")
	grid(col="grey")
	points(mine$BETAI, mine$ETAMED, pch=16, cex=1.1)
}

pdf("betai-eta-1.pdf", width=4.4, height=4.4)
par(las=1, mar=c(4, 3.8, 0.2, 0.2))

plo(range(phl_blue$BETAI), c(-2.7, -2.1))
mtext(expression(atop(eta[a], atop(phantom(1), phantom(1)))), 2, 2, cex=1.8)
lines(beta_i_samp, oml_ion_eta, lwd=2)

legend(0.17, -2.25, pch=c(5, 16), pt.cex=1.1, cex=1.1,
       c("SCEPTIC", "pot, sample median"), bg="white")
text(0.16, -2.15, "unmagnetized-ion\nprediction", cex=1.2)

dev.off()

pdf("betai-eta-2.pdf", width=4.4, height=4.4)
par(las=1, mar=c(4, 3.8, 0.2, 0.2))

plo()
mtext(expression(eta[a]), 2, 2, cex=1.8)

legend(0.4 * max(mine$BETAI), -2.61, pch=c(5, 16), pt.cex=1.1, cex=1.1,
       c("SCEPTIC", "pot, sample median"), bg="white")

dev.off()
