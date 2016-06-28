# ECS 289
# Project
source("OE.R")
library(png)
true_matrix = readPNG('image.png')[,,1]
svd = svd(true_matrix)
plot(svd$d, main="singular values of the true matrix")
plot(cumsum(svd$d))
par(mfrow=c(1,1))
par(mar=rep(1,4))
myplot(true_matrix)

visible_matrix = OE_visible(true_matrix = true_matrix, visibility = 0.3)
query_budge = round(length(true_matrix)*0.1)
r=2
system.time(PI <- OE_order(visible_matrix = visible_matrix, r, adjustment = FALSE))
system.time(r2t1 <- OE_extend(PI,visible_matrix, true_matrix, r, query_budge, condition_limit=1.1))
system.time(r2t1.5 <- OE_extend(PI,visible_matrix, true_matrix, r, query_budge, condition_limit=1.5))
system.time(r2t2 <- OE_extend(PI,visible_matrix, true_matrix, r, query_budge, condition_limit=2))
system.time(r2t2.5 <- OE_extend(PI,visible_matrix, true_matrix, r, query_budge, condition_limit=2.5))
system.time(r2t3 <- OE_extend(PI,visible_matrix, true_matrix, r, query_budge, condition_limit=3))
system.time(r2t3.5 <- OE_extend(PI,visible_matrix, true_matrix, r, query_budge, condition_limit=3.5))
system.time(r2t4 <- OE_extend(PI,visible_matrix, true_matrix, r, query_budge, condition_limit=4))
system.time(r2t1.01 <- OE_extend(PI,visible_matrix, true_matrix, r, query_budge, condition_limit=1.01))
par(mfrow=c(2,4))
par(mar=rep(1,4))
myplot(visible_matrix)
myplot(r2t1$recover)
myplot(r2t1.5$recover)
myplot(r2t2$recover)
myplot(r2t2.5$recover)
myplot(r2t3$recover)
myplot(r2t3.5$recover)
myplot(r2t4$recover)

r=5
system.time(PI <- OE_order(visible_matrix = visible_matrix, r, adjustment = FALSE))
system.time(r5t1 <- OE_extend(PI,visible_matrix, true_matrix, r, query_budge, condition_limit=1.1))
system.time(r5t1.5 <- OE_extend(PI,visible_matrix, true_matrix, r, query_budge, condition_limit=1.5))
system.time(r5t2 <- OE_extend(PI,visible_matrix, true_matrix, r, query_budge, condition_limit=2))
system.time(r5t2.5 <- OE_extend(PI,visible_matrix, true_matrix, r, query_budge, condition_limit=2.5))
system.time(r5t3 <- OE_extend(PI,visible_matrix, true_matrix, r, query_budge, condition_limit=3))
system.time(r5t3.5 <- OE_extend(PI,visible_matrix, true_matrix, r, query_budge, condition_limit=3.5))
system.time(r5t4 <- OE_extend(PI,visible_matrix, true_matrix, r, query_budge, condition_limit=4))
par(mfrow=c(2,4))
par(mar=rep(1,4))
myplot(visible_matrix)
myplot(r5t1$recover)
myplot(r5t1.5$recover)
myplot(r5t2$recover)
myplot(r5t2.5$recover)
myplot(r5t3$recover)
myplot(r5t3.5$recover)
myplot(r5t4$recover)

mylist = list(r2t1, r2t1.5,
              r2t2, r2t2.5,
              r2t3, r2t3.5,
              r2t4)
mylist2 = list(r5t1, r5t1.5,
              r5t2, r5t2.5,
              r5t3, r5t3.5,
              r5t4)

par(mfrow=c(1,2))
par(mar=rep(4,4))
plot(x=c(1,1.5,2,2.5,3,3.5,4), y=sapply(mylist, function(i) length(which(!is.na(i$query)))), lwd=2, type="l",lty=1, main="number of query v.s. threshold", ylab="number of query", xlab="threshold")
lines(x=c(1,1.5,2,2.5,3,3.5,4), y=sapply(mylist2, function(i) length(which(!is.na(i$query)))), lwd=2, lty=2, main="number of query v.s. threshold", ylab="number of query", xlab="threshold")
legend("topright", lty=c(1,2), lwd=2, c("rank 2","rank 5"), cex=0.7)
plot(x=c(1,1.5,2,2.5,3,3.5,4), y=sapply(mylist, function(i) i$error), type="l",lwd=2,lty=1, main="relative error v.s. threshold", ylab="relative error", xlab="threshold", ylim=c(0,1))
lines(x=c(1,1.5,2,2.5,3,3.5,4), y=sapply(mylist2, function(i) i$error), lty=2, lwd=2,main="number of query v.s. threshold", ylab="number of query", xlab="threshold")
legend("bottomleft", lty=c(1,2), lwd=2, c("rank 2","rank 5"), cex=0.7)