library(LambertW)

# SOML floating potential in singly-ionized plasma
soml <- function(theta, mime, u)
{
	s1 <- function(u)
	{
		s1_core <- function(u)
		{
			return((sqrt(pi) * (1 + (2 * u^2)) * erf(u) / (4 * u))
			       + (exp(-u^2) / 2))
		}
		result <- u
		result[result > 1e-14] <- s1_core(result[result > 1e-14])
		result[result < 1] <- 1
		return(result)
	}

	s2 <- function(u)
	{
		result <- u
		result[result > 1e-14] <-
			sqrt(pi) * erf(result[result > 1e-14]) /
			(2 * result[result > 1e-14])
		result[result <= 1e-14] <- 1
		return(result)
	}

	s1s2 <- s1(u) / s2(u)
	ts1s2 <- theta * s1s2
	return(ts1s2 - W(sqrt(mime * theta) * exp(ts1s2) / s2(u)))
}

pot <- read.table("pot-soml-eta.dat", header=TRUE)

# Turn the date string in `pot` into a real date.
pot$DAT <- as.Date(pot$DAT)

# Eliminate runs with the T_i / T_e ratio of about 1.485, characteristic
# of the naive, pre-SCEPTIC-algorithm reinjection algorithm.
pot <- pot[pot$TITEMED < 1.48,]

# The remaining runs will have T_i / T_e ratios of essentially 0.1 or
# essentially 1. Categorize them nicely according to their ratio.
pot$COOL <- pot$TITEMED < 0.3

# Renormalize `U` so that it's normalized by the special ion thermal speed
# rather than the Bohm speed.
pot$URE <- pot$U / sqrt(2 * signif(pot$TITEMED, 1))

# Pick out the subset of pot runs with adequately short time steps.
posh <- pot[pot$DT < 1.1,]

begin_soml_plot <- function(partial_path, d, x_name, x_la, over_soml=FALSE, y_li=-c(2.9, 1.9), y_name="ETAMED")
{
	pdf(paste("SOML-pot-", partial_path, ".pdf", sep=""),
	    width=5.7, height=5.7)
	if (over_soml) {
		par(las=1, mar=c(4.4, 6.8, 0.1, 0.5))
		y <- d[[y_name]] / soml(d$TITEMED, 1836.15, d$U / sqrt(2 * d$TITEMED))
		plot(d[[x_name]], y, type="n", ylim=y_li,
		     cex=1.2, cex.axis=1.4, cex.lab=1.4, xaxt="n", xlab=x_la, ylab="")
		abline(h=seq(0.95, 1.2, 0.05), col="grey", lty="dashed")
		abline(v=unique(d[[x_name]]), col="grey", lty="dotted")
		x_labels <- sort(unique(d[[x_name]]))
		if (x_name == "DT") {
			# Wicked-sick formatting abuse here, but it's the easiest way
			# to get nice tick labels (which fit) at both 0.5 ps and 1 ps.
			x_labels[1] <- expression(atop(phantom(0),
			                               atop(phantom(0),
			                                    atop(phantom(0), frac(1, 2)))))
		}
		axis(1, sort(unique(d[[x_name]])), x_labels, cex.axis=1.5)
		points(d[[x_name]], y, pch=c(21 - (20 * d$COOL)), bg="black")
		mtext(expression(frac(eta[a], eta[SOML])), 2, 3, cex=1.6)
	} else {
		par(las=1, mar=c(4.2, 4.6, 0.2, 0.2))
		plot(d[[x_name]], d[[y_name]], pch=c(21 - (20 * d$COOL)), bg="black",
		     ylim=y_li, cex=1.2, cex.axis=1.4, cex.lab=1.4,
		     xlab=x_la, ylab="")
		mtext(expression(atop(phantom(0), eta[a])), 2, 3, cex=1.6)
		grid(col="grey")
	}
	if (grepl("OMLFIT", x_name)) {
		legend_y <- -2.71
	} else {
		legend_y <- sum(c(0.412, 0.588) * y_li)
	}
	legend(sum(c(0.4, 0.6) * range(d[[x_name]])), legend_y,
	       c(expression(Theta == 0.1), expression(Theta == 1)),
	       pch=c(1, 21), pt.bg="black", bg="white", cex=1.4)
}

#begin_soml_plot(1, pot, "U", expression("drift speed " * (c[s])))
#curve(soml(1, 1836.15, x / sqrt(2 * 1)), add=TRUE, lwd=2)
#curve(soml(0.1, 1836.15, x / sqrt(2 * 0.1)), add=TRUE, lwd=2, lty="dashed")
#dev.off()
#
#begin_soml_plot(2, pot, "URE", 
#                expression("drift speed " * (sqrt(2 * k[B] * T[i] / m[i]))))
#curve(soml(1, 1836.15, x), add=TRUE, lwd=2)
#curve(soml(0.1, 1836.15, x), add=TRUE, lwd=2, lty="dashed")
#dev.off()
#
#begin_soml_plot(3, pot, "U", expression("drift speed " * (c[s])),
#                TRUE, c(0.955, 1.197))
#abline(h=1, lwd=2, lty="twodash")
#dev.off()
#
#begin_soml_plot(4, pot, "DT", expression(Delta * t * " (ps)"),
#                TRUE, c(0.955, 1.197))
#abline(h=1, lwd=2, lty="twodash")
#dev.off()

begin_soml_plot(5, posh, "U", expression("drift speed " * (c[s])))
curve(soml(1, 1836.15, x / sqrt(2 * 1)), add=TRUE, lwd=2)
curve(soml(0.1, 1836.15, x / sqrt(2 * 0.1)), add=TRUE, lwd=2, lty="dashed")
dev.off()
