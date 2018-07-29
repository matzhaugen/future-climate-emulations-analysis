# With this project we attempt to project current temperatures
# into the future using model simulations combined with reanalysis 
# data. The goes is to make a mapping between current and future 
# quantiles. For example, if we want to know the 95th quantile
# temperature in 2080 we first look at the current 95th quantile
# from reanalysis. 
# setwd("~/Extremes/emulation")
# We first raw data and the results are contained in two files,
# 1. everything_1463_1463_310_1010.Rdata which is the cpu 
# intensive stuff and 
# 2. maps_1463_1463_310_1010.Rdata which is the projected temperatures
# which we plot. The first file takes 24 hours with 8 cores/6Gb RAM each
# and the second takes 10 minutes with the same setup.
# We then define constants and produce each figure, which are made to be
# independent of each other.
source('emulationHelper.R')
load('raw_temperatures.RData')
co2 = read.table("RCP85_MIDYR_CONC.DAT", skip=38, header=TRUE)
load('everything_1463_1463_310_1010.Rdata')
# load('backup_1463_1463_310_1010.Rdata')
# IMPORTANT
# Load ONLY if this file exists and maps.R has run successfully.
# OTHERWISE run maps.R to do any plotting
load('maps_1463_1463_310_1010.Rdata')

q = c(.1, .18, 0.25, .35, .42, 0.5, .58, .65, 0.75, .82, 0.9)
q_tail = c(0.01, 0.1, .25, 0.5, .75)
q_low = q_tail
q_high = rev(1-q_tail)
q_all = c(q_low*q[1], q, q[length(q)] + q_high*(1-q[length(q)]))
cex.main=1.5
cex.lab=1.5
q_norm = c(.1, .5, .9)
Sys.setenv(TZ='UTC')
season = list(c(11, 0, 1), c(2,3,4), c(5, 6, 7), c(8,9,10))
location = list(c(30, 360-95), c(42, 360-90), c(48, 360-122), c(34, 360-118),
	c(40, 360-105), c(34, 360-112), c(41, 360-74), c(34, 360-84))
nloc = length(location)
season.names = c('Winter', 'Spring', 'Summer', 'Autumn')
# The first row here is the order in which the cities are stored in qout and map*
# The other rows refer to plotting orders
loc.names = c("Houston", "Chicago", "Seattle", "Los Angeles", "Denver", "Phoenix", "New York City", "Atlanta")
loc.names3 = c("Seattle", "Denver", "Chicago", "New York City", "Los Angeles", "Phoenix", "Houston", "Atlanta")
p.ryan = dim(y.ryan[[1]])[2]
p.lens = dim(y.lens[[1]])[2]
years_i = '1979/2016'
years_f = '2059/2096'


qout.lens = lapply(qout.all.inclusive, function(x) x$lens)
qout.ryan = lapply(qout.all.inclusive, function(x) x$ryan)
map.lens = lapply(maps.all.inclusive, function(x) x$lens)
map.ryan = lapply(maps.all.inclusive, function(x) x$ryan)

map.lens.ryan = lapply(maps.all.inclusive, function(x) x$lens.ryan)
map.ryan.lens = lapply(maps.all.inclusive, function(x) x$ryan.lens)
map.nested.ryan.lens = lapply(maps.nested, function(x) x$ryan.then.lens)


#########################################
########### FIGURE 1 ####################
#########################################
pdf('figures/MarginalComparison.pdf', height=6, width=12)
par(oma=c(0,2,0,2), mfrow=c(2,4))
cex.lab = 1.2
ylab = 'Density'
templab = expression(paste("Temperature [", degree*C, ']'))
for (j in 1:length(loc.names3)) {
	print(j)
	idx = which(loc.names3[j]==loc.names)
	comp.histlike(list(
		as.numeric(y.lens[[idx]][years_i,]), 
		as.numeric(y.ryan[[idx]][years_i,]),
		as.numeric(rea.y[[idx]])), breaks=60, lwd=2, xlab="", main=loc.names3[j], cex.main=2, ylab="")
	if (j==1) legend('topleft', c('LENS', 'Sriver', 'Reanalysis'), col=c(1,2,3), lwd=2)
	if (j %in% c(1,5)) mtext(ylab, side = 2, line = 3,  cex=cex.lab)
	if (j %in% c(5:8)) mtext(templab, side = 1, line = 3,  cex=cex.lab)
}
dev.off()

#########################################
########### FIGURE 2 ####################
#########################################
# Plot difference between LENS and RYAN raw data
templab = expression(paste("Temperature [", degree*C, ']'))
cex.lab=1.2
ylim = c(-5, 36)
j = which(loc.names3[1] == loc.names)
yqhat_real1 = make.quantile.surfaces(qout.ryan[[j]])  
yqhat_real2 = make.quantile.surfaces(qout.lens[[j]])    
mcolors = c(1, 2)
png('figures/sample.png', width=700, height=350)
par(mfrow=c(1,2))
year.str = '1920'
par(mar=c(5,5,4,2))
plot.frac(rep(.indexyday((y.lens[[j]][year.str,])), 40), y.lens[[j]][year.str,], col=alpha(mcolors[1], .3),
	xlab="Day of the year", ylab=templab, cex.lab=cex.lab, ylim=ylim)
points.frac(rep(.indexyday((y.ryan[[j]][year.str,])), 50), y.ryan[[j]][year.str,], col=alpha(mcolors[2], .3))
lines(yqhat_real2[,1,11], col=mcolors[1], lwd=2)
lines(yqhat_real1[,1,11], col=mcolors[2], lwd=2)
title(year.str)
legend('topleft', c('LENS', 'Sriver'), col=mcolors, lwd=2)
year.str = '2099'
par(mar=c(5,5,4,2))
plot.frac(rep(.indexyday((y.lens[[j]][year.str,])), 40), y.lens[[j]][year.str,], col=alpha(mcolors[1], .3),
	xlab="Day of the year", ylab=templab, cex.lab=cex.lab, ylim=ylim)
points.frac(rep(.indexyday((y.ryan[[j]][year.str,])), 50), y.ryan[[j]][year.str,], col=alpha(mcolors[2], .3))
lines(yqhat_real2[,180,11], col=mcolors[1], lwd=2)
lines(yqhat_real1[,180,11], col=mcolors[2], lwd=2)
title(year.str)
dev.off()

#########################################
########### FIGURE 3 ####################
#########################################
# This requires access to a X11 device, and wont work when using
# 'mosh', only 'ssh' if working on a remote.
#### Quantile Maps Figure 
lw <- list(left.padding = list(x = -0.1, units = "inches"))
lw$right.padding <- list(x = -0.2, units = "inches")
lh <- list(bottom.padding = list(x = -.2, units = "inches"))
lh$top.padding <- list(x = -.2, units = "inches")
zlabel = expression(paste("Change [", degree*C, ']'))
cex.lab = 1
zlab = list(label=zlabel, x=.5, rot=90, cex=cex.lab)
ylab = list(label='Percentile', y=.5, rot=-40, cex=cex.lab)
xlab = list(label='Days', y=.5, rot=30, cex=cex.lab)
# plot Quantile Surfaces
logplot = FALSE
if (logplot) {
	nqa = length(q_all)
	lq = log(q_all, .5)-1
	hq = log(1-q_all, .5)-1
	q_axis = c(-lq[1:floor(nqa/2)], hq[(floor(nqa/2)+1): nqa])
	qticks = c(1,3, 6,11,15,19,21)
} else {
	q_axis = q_all
	qticks = c(6, 8, 11, 14, 16)
}

k = 1
p.ryan = list(); length(p.ryan) = nloc
p.lens = list(); length(p.lens) = nloc
cex.tck = .8
mdays = seq(1, 365, 10)
q_all = qout.ryan[[1]]$q_all
for (j in 1:nloc) {
	print(j)
	dy.ryan = get.quantile.change(qout.ryan[[j]])
	dy.lens = get.quantile.change(qout.lens[[j]])
	zlim = range(dy.ryan, dy.lens)
	p.ryan[[j]] = mywireframe(mdays, rev(q_axis), dy.ryan[, length(q_axis):1], 
		scales=list(arrows=FALSE, y=list(arrows=F, at=q_axis[qticks], 
      		labels=rev(q_axis[qticks])), col=1, cex=cex.tck, distance=c(1.3,1.3,1.3), 
			tck=.8),
		main=list(label=loc.names[j], vjust=2.5, hjust=.5),
		zlab=(zlab),
    par.settings = list(axis.line = list(col = 'transparent')),
    lattice.options = list(layout.widths = lw, layout.heights = lh),
    ylab=ylab, 
    xlab=xlab, zlim=zlim)
	
	p.lens[[j]] = mywireframe(mdays, rev(q_axis), dy.lens[, length(q_axis):1], 
		scales=list(arrows=FALSE, y=list(arrows=F, at=q_axis[qticks], 
      		labels=rev(q_axis[qticks])), col=1, cex=cex.tck, distance=c(1.3,1.3,1.3), 
			tck=.8),
		main=list(label=""),
		zlab=(zlab),
    par.settings = list(axis.line = list(col = 'transparent')),
    lattice.options = list(layout.widths = lw, layout.heights = lh),
    ylab=ylab, 
    xlab=xlab, zlim=zlim)
}
loc.ind = unlist(lapply(loc.names3, function(x) {which(x==loc.names)}))
pp = arrangeGrob(grobs=c(p.ryan[loc.ind[1:4]], p.lens[loc.ind[1:4]], p.ryan[loc.ind[5:8]], p.lens[loc.ind[5:8]]), ncol=4, nrow=4)
# pp1 = arrangeGrob(grobs=interlace(p.ryan[1:4], p.lens[1:4]), ncol=4, nrow=2)
# pp2 = arrange 				Grob(grobs=interlace(p.ryan[5:8], p.lens[5:8]), ncol=4, nrow=2)
k = ggsave('figures/multiQuantileSurfAllInclusive.png', pp, width=11, height=11)


#########################################
########### FIGURE 4 ####################
#########################################
j=1
p = plot.model.difference.all.inclusive(qout.ryan[[j]], qout.lens[[j]], mq=c(.01, .5, .99))
pp = arrangeGrob(grobs=p, ncol=3, nrow=1)
ggsave('figures/modelDifferenceAllInclusive.pdf', pp, width=11, height=4)


#########################################
########### FIGURE 5 ####################
#########################################
# Plot double transform first lens and then Ryans back again
# Using the all.inclusive data set (using all 40 simulations for the 
# model building and no resampling)
xlab = list(label='Days', y=.5, rot=30, cex=cex.lab)
zlabel = expression(paste("Change [", degree*C, ']'))
zlab = list(label=zlabel, x=.48, rot=90, cex=cex.lab)
qticks = qticks = c(6, 8, 11, 14, 16)
mdays = seq(1, 365, 10)
scales=list(arrows=FALSE, draw=TRUE, tck=.5, col=1,
          y=list(at=rev(q_all[qticks]), labels=rev(q_all[qticks])))

p = list(); length(p) = 4
for (j in 1:4) {
	idx = which(loc.names3[j] == loc.names)
	p[[j]] = plot.double.transform.all.inclusive(qout.lens[[idx]], qout.ryan[[idx]], 
		main=list(label=loc.names3[j], vjust=2.5), zlim=c(-.5, 6), xlab=xlab, zlab=zlab)
}
pp = arrangeGrob(grobs=p, ncol=4)
ggsave('figures/QuantileDoubleAllInclusive.png', pp, width=13, height=4)


#########################################
########### FIGURE 6 ####################
#########################################
# Plot jacknife variance estimates
qout.lens = lapply(qout.jack, function(x) x$lens)
m=0;null.idx = unlist(lapply(qout.lens[[3]], function(x) {m<<-m+1; if(is.null(x$q)) m}))
res.lens = lapply(qout.lens, model.variance.jackknife, mc.cores=16) # Takes a minute
pp.lens = plot.jackknife(res.lens)
ggsave('figures/JackknifeStdLens.png', pp.lens, width=15, height=10)
	

#########################################
########### FIGURE 7 ####################
#########################################
# Plot marginal plot with full current and future projections
j=1
mseason=1
idx = which(loc.names3[j] == loc.names)
y.lens.f = y.lens[[idx]][years_f,]
y.lens.i = y.lens[[idx]][years_i,]
y.f.hat1 = map.lens[[idx]]$rea.y.f
y.f.hat2 = map.ryan.lens[[idx]]$rea.y.f
png('figures/marginalSingleAllInclusive.png')
marginal.summary.plot.multi(list(
	  y.lens.f[.indexmon(y.lens.f) %in% mseason], 
      y.f.hat1[.indexmon(y.f.hat1) %in% mseason], 
      y.f.hat2[.indexmon(y.f.hat2) %in% mseason], 
      y.lens.i[.indexmon(y.lens.i) %in% mseason]))
cex.lab=1.8
legend('topleft', c('Target', 'Est-LENS', 'Est-SFK15', 'Initial'), lty=1, col=c(1,2, 3, 4), cex=1.3)
mtext(expression(paste("Estimate [", degree*C, ']')), side = 4, line = 3,  cex=cex.lab)
mtext('Density', side = 2, line = 3,  cex=cex.lab)
mtext(expression(paste("Target [", degree*C, ']')), side = 1, line = 3,  cex=cex.lab)
dev.off()

#########################################
########### FIGURE 8,9,10 ####################
#########################################
## Plot seasonal marginal plots
# Winter
pheight = 500
pwidth = 1000
y2lab = expression(paste("Estimate [", degree*C, ']'))
ylab = 'Density'
cex.lab=1.5
i = 1 # winter idx
png('figures/marginalsWinterLensRyanOnLens.png', width=pwidth, height=pheight)
multi.marginal.summary.plot.all.new(location, loc.names3, season[[i]], years_f, 
	list(map.lens, map.ryan.lens), y.lens, 
	include.all.runs=FALSE)
dev.off()

i = 3 #summer
png('figures/marginalsSummerLensRyanOnLens.png', width=pwidth, height=pheight)
multi.marginal.summary.plot.all.new(location, loc.names3, season[[i]], years_f, 
	list(map.lens, map.ryan.lens), y.lens, 
	include.all.runs=FALSE, leg.pos=2)
dev.off()

i = 1 # winter idx
png('figures/marginalsWinterLensRyanOnLensAllRuns.png', width=pwidth, height=pheight)
multi.marginal.summary.plot.all.new(location, loc.names3, season[[i]], years_f, 
	list(map.lens, map.ryan.lens), y.lens, 
	include.all.runs=TRUE)
dev.off()


#########################################
########### FIGURE 11 ####################
#########################################
# Do the same with reanalysis projected into 2099
mseason=1
j=1
idx = which(loc.names3[j] == loc.names)
y.f.hat.lens = maps.all.inclusive.const.year.rea[[idx]]$lens$rea.y.f
y.f.hat.ryan = maps.all.inclusive.const.year.rea[[idx]]$ryan$rea.y.f
y.f.hat.ryan.then.lens = maps.nested[[idx]]$ryan.then.lens.rea$rea.y.f
png('figures/marginalSingleReaNested.png')
par(oma=c(0,1,0,3))
marginal.summary.plot.multi(list(
	  rea.y[[idx]][.indexmon(rea.y[[idx]]) %in% mseason], 
      y.f.hat.lens[.indexmon(y.f.hat.lens) %in% mseason], 
      y.f.hat.ryan[.indexmon(y.f.hat.ryan) %in% mseason],
      y.f.hat.ryan.then.lens[.indexmon(y.f.hat.ryan.then.lens) %in% mseason]))
cex.lab=1.8
legend('topleft', c('Initial', 'Est-LENS', 'Est-SFK15', 'Est-SFK15->Lens'), lty=1, col=c(1,2,3,4), cex=1.3)
mtext(expression(paste("Estimate [", degree*C, ']')), side = 4, line = 3,  cex=cex.lab)
mtext('Density', side = 2, line = 3,  cex=cex.lab)
mtext(expression(paste("Target [", degree*C, ']')), side = 1, line = 3,  cex=cex.lab)
dev.off()


#########################################
########### MISC FIGURE #################
########## NOT IN PAPER #################
#########################################
# Didn't quite make it
i = 1 # winter idx
png('figures/marginalsWinterLensIdentityAllRuns.png', width=pwidth, height=pheight)
multi.marginal.summary.plot.all(location, loc.names3, season[[i]], years_f, 
	list(map.lens, maps.lens.identity), y.lens, 
	include.all.runs=TRUE, leg.str=c('Target', 'Est-LENS', 'Est-Identity'), 1)
dev.off()


i = 1 # winter idx
png('figures/marginalsWinterLensIdentity.png', width=pwidth, height=pheight)
multi.marginal.summary.plot.all(location, loc.names3, season[[i]], years_f, 
	list(map.lens, maps.lens.identity), y.lens, 
	include.all.runs=FALSE, leg.str=c('Target', 'Est-LENS', 'Est-Identity'), 1)
dev.off()



i = 3
png('figures/marginalsSummerLensIdentity.png', width=pwidth, height=pheight)
multi.marginal.summary.plot.all(location, loc.names3, season[[i]], years_f, 
	list(map.lens, maps.lens.identity), y.lens, 
	include.all.runs=FALSE, leg.str=c('Target', 'Est-LENS', 'Est-Identity'), 1)
dev.off()


y.f = lapply(1:length(location), function(j) {
	y.f = y.lens[[j]][years_f,]	
	lapply(season, function(mseason) y.f[.indexmon(y.f) %in% mseason])
})

y.f.hat = lapply(list(map.lens, map.ryan, maps.lens.identity), function(mi) {
	lapply(mi, function(micity) {
		y.all = xts(do.call(cbind, lapply(micity, function(x) 
    	    as.matrix(x$rea.y.f))), order.by=time(micity[[1]]$rea.y.f[,1]))	
		lapply(season, function(mseason) y.all[.indexmon(y.all) %in% mseason])
	})
})   
dfloss = lapply(y.f.hat, function(x) get.percentile.diff(y.f, x)) #Takes a minute

filenames = c('figures/multi/multiBoxLens.png', 
		'figures/multi/multiBoxRyanOnLens.png',
	'figures/multi/multiBoxLensIdentity.png')

for (i in 1:length(dfloss)) {
	print(filenames[i])
	png(filenames[i], width=1000, height=500)
	ggplot(dfloss[[i]]) + 
	geom_boxplot(aes(x=Percentile, y=Loss, fill=Season)) +
	scale_y_continuous(name=expression(paste(tau - hat(tau), " [",degree*C,"]", sep=''))) +  
	coord_flip() + 
	facet_wrap(~City, ncol=4, nrow=2) +
	theme(axis.title.x=element_text(size=25),
		axis.title.y=element_text(size=25),
		strip.text=element_text(size=20, face='bold'),
		legend.title = element_text(size=20),
		axis.text = element_text(size=15))
	dev.off()
}


png('figures/DetailedComparison.png', width=860, height=860)
cex.axis=1.4
cex.main=1.8
idx = 1
# y.ryan = getModelData2(location[[idx]][1], location[[idx]][2])['1920/2099']
# y.lens = getLENSData(location[[idx]][1], location[[idx]][2], 1)['1920/2099']
par(mfrow=c(2,2), oma=c(0,1,0,0), xpd=NA)
plot.quantile.fit(y.lens[[idx]], qout.lens[[idx]], tau=c(.1, .9), cex.main=cex.main,
	plotit=TRUE, normit=F, by.year=T, main="LENS", ylab="Temperature (deg C)", cex.axis=cex.axis)
plot.quantile.fit(y.ryan[[idx]], qout.ryan[[idx]], tau=c(.1, .9), cex.main=cex.main,
	plotit=TRUE, normit=F, by.year=T, main="Sriver", ylab="", cex.axis=cex.axis)
plot.quantile.fit(y.lens[[idx]], qout.lens[[idx]], tau=c(.1, .9), cex.main=cex.main,
	plotit=TRUE, normit=F, by.year=F, main="LENS", ylab="Temperature (deg C)", cex.axis=cex.axis)
plot.quantile.fit(y.ryan[[idx]], qout.ryan[[idx]], tau=c(.1, .9), cex.main=cex.main,
	plotit=TRUE, normit=F, by.year=F, main="Sriver", ylab="", cex.axis=cex.axis)
dev.off()






