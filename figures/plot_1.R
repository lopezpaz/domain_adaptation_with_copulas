library(copula)
library(reshape2)
library(grid)
library(ggplot2)
library(np)

my_color <- rgb(0.02, .4, .7)

ke <- function(x)
{
bw <- npudensbw(as.matrix(x), bwmethod="normal-reference")
list(
pdf = function(x) npudens(bw, edat=x)$dens,
cdf = function(x) npudist(bw, edat=x)$dist,
H   = bw$bw,
x   = x
)
}

NPCopula_fit <- function(sample)
{
sampleR <- qnorm(as.matrix(sample))
kest    <- ke(sampleR)
pdf     <- function(u,v)
{
kest$pdf(cbind(qnorm(u), qnorm(v))) / dnorm(qnorm(u)) / dnorm(qnorm(v))
}

d       <- sampleR
H       <- kest$H
list(pdf = pdf, d = sampleR)
}

x <- read.table("wifi.txt")[,c(4,11)]
x[x==1] <- .9999

# Plot sample
qplot(x[,1],x[,2], colour=I(densCols(x)), size=I(3), xlab="", ylab="")#+
#theme(axis.ticks  = element_blank(),
#      axis.text   = element_blank(),
#      panel.background =  theme_blank(),
#      panel.grid =  element_blank(),
#      plot.margin = unit(c(0,0,-1.5,-1.5),"lines"))
ggsave("wifi_sample.pdf")

cs <- melt(outer(seq(0.01,0.99,0.01),
                 seq(0.01,0.99,0.01),
                 function(x,y) dCopula(cbind(x,y), normalCopula(-0.008483453))))

names(cs) <- c("x", "y", "z")

ggplot(cs, aes(x,y,z=z), xlab="", ylab="")+stat_contour(colour=I(my_color), xlab="", ylab="")+xlab("")+ylab("")#+
#theme(axis.ticks  = element_blank(),
#      axis.text   = element_blank(),
#      panel.background =  theme_blank(),
#      panel.grid =  element_blank(),
#      plot.margin = unit(c(0,0,-1.6,-1.6),"lines"))

ggsave("wifi_gaussian.pdf")

npcop <- NPCopula_fit(x)

kn <- melt(outer(seq(0.01,0.99,0.01),
                 seq(0.01,0.99,0.01),
                 function(x,y) npcop$pdf(x,y)))

names(kn) <- c("x","y","z")

ggplot(kn, aes(x,y,z=z), xlab="", ylab="")+stat_contour(colour=I(my_color), xlab="", ylab="")+xlab("")+ylab("")#+
#theme(axis.ticks  = element_blank(),
#      axis.text   = element_blank(),
#      panel.grid =  element_blank(),
#      panel.background =  theme_blank(),
#      plot.margin = unit(c(0,0,-1.6,-1.6),"lines"))

ggsave("wifi_contour.pdf")

rc <- rCopula(3000, normalCopula(0.8))

qplot(rc[,1],rc[,2], colour=I(densCols(rc)), size=I(3), xlab="", ylab="")+
#theme(axis.ticks  = element_blank(),
#      axis.text   = element_blank(),
#      panel.background =  theme_blank(),
#      panel.grid =  element_blank(),
#      plot.margin = unit(c(0,0,-1.6,-1.6),"lines"))+
      scale_x_continuous(limits=c(min(rc[,1]),max(rc[,1]))) + 
      scale_y_continuous(limits=c(min(rc[,2]),max(rc[,2]))) + 
      geom_rug(col=rgb(.2,0.6,0.8,alpha=.2))
      #geom_rug(col=rgb(.2,.6,.8,alpha=.2))

ggsave("system0.pdf")

rc1 <- cbind(qnorm(rc[,1]), qexp(rc[,2]))

qplot(rc1[,1],rc1[,2], colour=I(densCols(rc1)), size=I(3), xlab="", ylab="")+
#theme(axis.ticks  = element_blank(),
#      axis.text   = element_blank(),
#      panel.background =  theme_blank(),
#      panel.grid =  element_blank(),
#      plot.margin = unit(c(0,0,-1.6,-1.6),"lines"))+
      scale_x_continuous(limits=c(min(rc1[,1]),max(rc1[,1]))) + 
      scale_y_continuous(limits=c(min(rc1[,2]),max(rc1[,2]))) + 
#      #geom_rug(col=rgb(.2,.6,.8,alpha=.2))+
      geom_rug(col=rgb(.2,0.6,0.8,alpha=.2))

ggsave("system1.pdf")

qbimodal <- function(x) {
  for(i in 1:length(x)) {
    if(runif(1) < .5) { 
      x[i] <- (qnorm(x[i], mean=0, sd=0.5))
    } else {
      x[i] <- (qnorm(x[i], mean=5, sd=0.5))
    }
  }
  x
}

rc2 <- cbind(qbimodal(rc[,1]), qbimodal(rc[,2]))

qplot(rc2[,1],rc2[,2], colour=I(densCols(rc2)), size=I(3), xlab="", ylab="")+
#theme(axis.ticks  = element_blank(),
#      axis.text   = element_blank(),
#      panel.background =  theme_blank(),
#      panel.grid =  element_blank(),
#      plot.margin = unit(c(0,0,-1.6,-1.6),"lines"))+
#      scale_x_continuous(limits=c(min(rc2[,1]),max(rc2[,1]))) + 
#      scale_y_continuous(limits=c(min(rc2[,2]),max(rc2[,2]))) + 
      #geom_rug(col=rgb(.2,.6,.8,alpha=.2))
      geom_rug(col=rgb(.2,0.6,0.8,alpha=.2))

ggsave("system2.pdf")
