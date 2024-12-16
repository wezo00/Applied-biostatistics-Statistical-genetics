library(readr)
beerhall <- read.table("beerhall.txt")
View(beerhall)
attach(beerhall)

# numerical univariate
summary(Crime)
summary(Ale)
summary(School)
summary(Worship)

# graphical univariate
layout(matrix(1:2,ncol=2))
plot(Crime)
plot(Ale)
plot(School)
plot(Worship)
par(mfrow=c(1,2))
#par(mfrow=c(2,2))
par(pty="s")
par(pty="m")

barplot(Crime)
# numerical bivariate
beerhall[c(5,6,7)]
data.frame(beerhall[c(5,6,7)])
cor(data.frame(beerhall[c(4,5,6,7)]))
matrix(beerhall[c(5,6,7)])
# graphical bivariate
pairs( ~ Ale + School + Worship, data=beerhall)
pairs( ~ Crime + Ale + School + Worship, data=beerhall, gap = 1/5)

names(beerhall)
beerhall$Ale <- factor(beerhall$Ale)
beerhall$School<- factor(beerhall$School)
beerhall$Worship <- factor(beerhall$Worship)
plot(beerhall)
crimemod <- lm(Crime ~ Ale + School + Worship, data=beerhall)
summary(crimemod)

pairs( ~ Ale + School + Worship, data=beerhall)

lm.null<-lm(Crime ~ 1)
summary(lm.null)
lm.ale<-lm(Crime ~ Ale)
summary(lm.ale)
lm.worship<-lm(Crime ~ Worship)
summary(lm.worship)
lm.school<-lm(Crime ~ School)
summary(lm.school)
lm.ale.school<-lm(Crime ~ Ale + School)
summary(lm.ale.school)
lm.ale.worship<-lm(Crime ~ Ale + Worship)
summary(lm.ale.worship)

#qq normal plot of residuals and resididuals vs fitted
plot(lm.null, which = 2)
plot(lm.ale.school, which = 1)
plot(lm.ale.school, which = 2)

plot(lm.ale.worship, which = 1)
plot(lm.ale.worship, which = 2)

lm.full<-lm(Crime ~ Ale + Worship + School)
summary(lm.full)

AIC(lm.null); AIC(lm.ale); AIC(lm.school); AIC(lm.ale.school)
lms <-step(lm(Crime ~ Ale + School + Worship))
lms <-step(lm(Crime ~ Ale + School))
lms <-step(lm(Crime ~ 1))

require(car)
vif(lm(Crime ~ Ale + School + Worship))

layout(matrix(1:2,ncol=2))
plot(Crime ~ Ale)
abline(lm.ale.school)
plot(lm.ale.school, which = 1)

par(mfrow=c(1,2))
tmp <- qqnorm(residuals(lm.ale.school), pch=20,
                main="Normal Q-Q plot, residuals")
qqline(residuals(lm.ale.school))
diff <- (tmp$x - tmp$y)
### label the residuals that deviate too far from the line
text(tmp$x, tmp$y, ifelse((abs(diff) > 3), names(diff), ""), pos=2)
rm(tmp,diff)
### observed vs. fitted
plot(Crime ~ fitted(lm.ale.school), pch=20,
         xlim=c(0,30), ylim=c(0,30),
         xlab="Fitted",ylab="Observed",
         main="Observed vs. Fitted")
abline(0,1); grid(col="black")
par(mfrow=c(1,1))

xtable(Crime)

# make labels and margins smaller
par(mai=c(0.1,0.1,0.2,0.1))
par(mai=c(0.1,0,0.1,0))

par(fig=c(0.1,0.7,0.3,0.9))
pairs( ~ Ale + School + Worship,gap=1/4)
# define area for the boxplot
par(fig=c(0.8,1,0,1), new=TRUE)
boxplot(Temperature)
# define area for the stripchart
par(fig=c(0.1,0.67,0.1,0.25), new=TRUE)
stripchart(Temperature, method="jitter")