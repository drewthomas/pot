deb <- read.table("debug.dat")
colnames(deb) <- c("A", "N", "MIME", "NLL")
N <- length(deb[,1])
par(mar=c(4,4,1,1))
plot(deb$A, deb$N, pch=21, col=hsv(0, 0, (N:1)/N, 0.8), bg=hsv(0, 0, 0, 0.04),
#plot(deb$A, deb$N, pch=4,
     cex=5*(N:1)/N, xlab="a", ylab=expression(n[e]), log="y")
grid()
lines(deb$A, deb$N, lty="dotted")

#library(scatterplot3d)
#scatterplot3d(deb$A, deb$N, deb$MIME, type="b")

library(rgl)
plot3d(deb$A, deb$N, deb$MIME)
