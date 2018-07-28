library(ncdf4)
library(quantreg)
library(pbs)
library(splines)
library(parallel)
library(gplots)
library(fields)
library(RColorBrewer)
library(xts)
# library(Matching)
# library(lubridate)
library(reshape2)
library(gridExtra) # also loads grid
library(grid)
library(ggplot2)
library(lattice)
library(Hmisc)
library(scales)
rwb.palette=colorRampPalette(c('darkblue','white','red2'),interpolate='spline')
wr.palette=colorRampPalette(c('white','red2'),interpolate='spline')

export.LENS.RYAN.data <- function(location) {
  y.ryan = list(); length(y.ryan) = length(location)
  y.lens = list(); length(y.lens) = length(location)
  rea.y = list(); length(rea.y) = length(location)
  i = 1
  for (loc in location) {
    print(paste("Exporting location #", i, "/", length(location), sep=""))
    mlat = loc[1]
    mlon = loc[2]
    #Import model data
    y.ryan[[i]] = getModelData2(mlat, mlon)['1920/2099']
    y.lens[[i]] = getLENSData(mlat, mlon, 1)['1920/2099']    
    reanalysis_y = getReanalysisData2(mlat, mlon)
    reanalysis_y = reanalysis_y[.indexhour(reanalysis_y) %in% 12]
    dayofyear = as.numeric(strftime(time(reanalysis_y), format = "%j"))
    rea.y[[i]] = reanalysis_y[dayofyear!=366]
    i = i + 1
  }
  save(y.lens, y.ryan, rea.y, file="raw_temperatures.RData")
}

plot.frac <- function(x, y, frac=.03, ...) {
  x = as.numeric(x)
  y = as.numeric(y)
  n = length(x)
  idx = sample(1:n, n*frac)  
  plot(x[idx], y[idx], ...)
}

points.frac <- function(x, y, frac=.03, ...) {
  x = as.numeric(x)
  y = as.numeric(y)
  n = length(x)
  idx = sample(1:n, n*frac)  
  points(x[idx], y[idx], ...)
}

outer_concat <- function(la, lo) {
  mla = matrix(rep(la, length(lo)), length(la), length(lo))
  mlo = matrix(rep(lo, length(la)), length(lo), length(la))
  latlon = cbind(c(t(mla)), c(mlo))
  latlon
}
# Get volcanic forcing timeseries
# time - the time points on which to interpolate the volcanic time series
# mlat - desired latitude: used to find the latitude
# of volcanic data closest the location.
get.volcanic <- function(mlat) {
  days.per.month = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  cdpm = cumsum(days.per.month)
  ncdata = nc_open('/project/moyer/mahaugen/proj2/hindcast/CCSM4_volcanic_1850-2008_prototype1.nc')
  volc = ncvar_get(ncdata, 'MMRVOLC')
  lev = ncvar_get(ncdata, 'lev')
  lat = ncvar_get(ncdata, 'lat')
  date = ncvar_get(ncdata, 'date')
  endidx = length(time)
  nc_close(ncdata)
  # find the latitude closest to the location
  dlat = abs(lat - mlat)
  volc.lat.idx = which(dlat == min(dlat))
  mean.volc = apply(volc[volc.lat.idx, , ], 2, mean)
  # plot(mean.volc, lty=1)
  volc.date = as.Date(as.character(date), '%Y%m%d')
  # volc.monthly = xts(mean.volc, order.by=volc.date)

  # Convert days since date to numeric YYYYMMDD
  # from.days.to.numeric(time)
  origin.year = 1850
  time = 1:(365*250)
  year = floor(time[-91251] / 365)
  month = rep(0, 250*365)
  d = time[-91251] %% 365

  for (i in 1:12) {
    if (i == 1) {
      temp = ((time[-91251] %% 365) <= cdpm[i])*1
      month = month + temp*i
      day = d
    } else {
      temp = ((time[-91251] %% 365) <= cdpm[i] &
                (time[-91251] %% 365) > cdpm[i-1])*1
      month = month + temp*i
      d[d < cdpm[i] & d >= cdpm[i-1]] = d[d < cdpm[i] & d >= cdpm[i-1]] - cdpm[i-1]
    }
  }

  sim.date = (origin.year + year) * 1e4 + month * 1e2 + d + 1
  volc.daily = rep(NA, 365*250)
  volc.daily[sim.date%in%date] = mean.volc[-1]
  volc.daily.int = na.approx(volc.daily, na.rm=FALSE)
  volc.daily.int[is.na(volc.daily.int)] = 0
  volc.daily.int = (volc.daily.int) / sd(volc.daily.int)
}

printf <- function(...) invisible(print(sprintf(...)))

histlike <- function(breaks, values, col=1, xlim=range(breaks), 
  ylim=c(0, max(values)), add=FALSE, ...) {
  if (!add) {
    plot(NA, ylim=ylim, xlim=xlim, bty='n', ...)
  } else {

  }
  nb = length(breaks)
  nv = length(values)
  segments(x0=breaks[1:(nb-1)], x1=breaks[2:nb], y0=values, y1=values, col=col)
  segments(x0=breaks[2:(nb-1)], x1=breaks[2:(nb-1)], y0=values[1:(nv-1)], y1=values[2:nv], col=col)
  segments(x0=breaks[c(1,nb)], x1=breaks[c(1,nb)], y0=c(0,0), y1=values[c(1,nv)], col=col)
  abline(h=0)
}

sym_log_breaks <- function (n = 5, base = 10)
{
    function(x) {
        rng <- range((log(abs(x)+base, base) - 1)*sign(x))
        print(rng)
        min <- floor(rng[1])
        max <- ceiling(rng[2])
        if (max == min)
            return(base^min)
        by <- floor(min(max, abs(min))/n) + 1
        positive.side = base^(seq(0, max, by = by) + 1) - base
        negative.side = -rev(base^(seq(0, abs(min), by = by) + 1) - base)
        c(round(negative.side), round(positive.side))
    }
}

sym.log.trans <- function(base = exp(1))
{

    trans <- function(x) (log(abs(x)+base, base) - 1)*sign(x)
    inv <- function(x) (base^(abs(x) + 1) - base) * sign(x)
    trans_new(paste0("sym.log-", format(base)), trans, inv, sym_log_breaks(base = base),
        domain = c(-Inf, Inf))
}

multi.marginal.summary.plot.all.new <- function(location, loc.names3, mseason, myears, 
  maps, raw.data, include.all.runs=FALSE, leg.str=c('Target', 'Est-LENS', 'Est-SFK15'), leg.pos=1) {
  cex.lab = 1.8
  y2lab = expression(paste("Estimate [", degree*C, ']'))
  xlab = expression(paste("Target [", degree*C, ']'))
  ylab = 'Density'
  par(oma=c(2,2,2,3))
  par(mfrow=c(2,4))
  for (j in 1:nloc) {   
    print(j)
    idx = which(loc.names3[j] == loc.names)
    print(loc.names3[j])
    y.f = raw.data[[idx]][myears,]
    if (include.all.runs) {
      y.f.hat = lapply(maps, function(mi) {mi[[idx]]$rea.y.f})
    } else {
      y.f.hat = lapply(maps, function(mi) {mi[[idx]]$rea.y.f[,1]})
    }

    y.f.hat = lapply(y.f.hat, function(x) x[.indexmon(x) %in% mseason])
    y.f = y.f[.indexmon(y.f) %in% mseason,]
    marginal.summary.plot.multi(c(list(y.f), y.f.hat))    
    mtext(loc.names3[j], side = 3, line = 3,  cex=cex.lab)
    if (j==leg.pos) legend('topleft', leg.str, lty=1, col=seq(length(leg.str)), cex=1.5)
    if (j %in% c(4,8)) mtext(y2lab, side = 4, line = 4,  cex=cex.lab)
    if (j %in% c(1,5)) mtext(ylab, side = 2, line = 3,  cex=cex.lab)
    if (j %in% c(5:8)) mtext(xlab, side = 1, line = 4,  cex=cex.lab)
  }  
}


multi.marginal.summary.plot.all <- function(location, loc.names3, mseason, myears, 
  maps, model.getter=getLENSData, include.all.runs=FALSE, leg.str=c('Target', 'Est-LENS', 'Est-SFK15'), leg.pos=1) {
  cex.lab = 1.8
  y2lab = expression(paste("Estimate [", degree*C, ']'))
  xlab = expression(paste("Target [", degree*C, ']'))
  ylab = 'Density'
  par(oma=c(2,2,2,3))
  par(mfrow=c(2,4))
  for (j in 1:nloc) {   
    print(j)
    idx = which(loc.names3[j] == loc.names)
    print(loc.names3[j])
    y.f = model.getter(location[[idx]][1], location[[idx]][2], 1)[myears,]
    if (include.all.runs) {
      y.f.hat = lapply(maps, function(mi) {xts(do.call(cbind, lapply(mi[[idx]], function(x) 
        as.matrix(x$rea.y.f))), order.by=time(mi[[idx]][[1]]$rea.y.f[,1]))})   
    } else {
      y.f.hat = lapply(maps, function(x){x[[idx]][[1]]$rea.y.f[,1] })
    }

    y.f.hat = lapply(y.f.hat, function(x) x[.indexmon(x) %in% mseason])
    y.f = y.f[.indexmon(y.f) %in% mseason,]
    marginal.summary.plot.multi(c(list(y.f), y.f.hat))    
    mtext(loc.names3[j], side = 3, line = 3,  cex=cex.lab)
    if (j==leg.pos) legend('topleft', leg.str, lty=1, col=seq(length(leg.str)), cex=1.5)
    if (j %in% c(4,8)) mtext(y2lab, side = 4, line = 4,  cex=cex.lab)
    if (j %in% c(1,5)) mtext(ylab, side = 2, line = 3,  cex=cex.lab)
    if (j %in% c(5:8)) mtext(xlab, side = 1, line = 4,  cex=cex.lab)
  }  
}


multi.marginal.summary.plot.both <- function(location, loc.names3, mseason, myears, 
  map1, map2, model.getter=getLENSData, include.all.runs=FALSE, leg.pos=1, ...) {
  cex.lab = 1.8
  y2lab = expression(paste("Estimate [", degree*C, ']'))
  ylab = 'Density'
  par(oma=c(2,2,0,3))
  layout(matrix(c(1,3,5,7, 2,4,6,8, 9, 11, 13, 15, 10, 12, 14, 16), 
    4, 4, byrow = TRUE))
  for (j in 1:nloc) {   
    print(j)
    idx = which(loc.names3[j] == loc.names)
    print(loc.names3[j])
    y.f = model.getter(location[[idx]][1], location[[idx]][2], ...)[myears,]
    if (include.all.runs) {
      y.f.hat1 = xts(do.call(cbind, lapply(map1[[idx]], function(x) as.matrix(x$rea.y.f))), order.by=time(map1[[1]][[1]]$rea.y.f))  
      y.f.hat2 = xts(do.call(cbind, lapply(map2[[idx]], function(x) as.matrix(x$rea.y.f))), order.by=time(map2[[1]][[1]]$rea.y.f))  
    } else {
      y.f.hat1 = map1[[idx]][[1]]$rea.y.f[,1]
      y.f.hat2 = map2[[idx]][[1]]$rea.y.f[,1]
    }
    
    marginal.summary.plot(y.f[.indexmon(y.f) %in% mseason], 
      y.f.hat1[.indexmon(y.f.hat1) %in% mseason],
      main=loc.names3[j], cex.lab=cex.lab, cex.main=2, ylab="")    
    if (j==leg.pos) legend('topleft', c('Target', 'Estimate'), lty=1, col=c(1,2), cex=1.5)
    if (j %in% c(4,8)) mtext(y2lab, side = 4, line = 4,  cex=cex.lab)
    if (j %in% c(1,5)) mtext(ylab, side = 2, line = 3,  cex=cex.lab)
    marginal.summary.plot(y.f[.indexmon(y.f) %in% mseason], 
      y.f.hat2[.indexmon(y.f.hat2) %in% mseason],
      main="", cex.lab=cex.lab, cex.main=2, ylab="")
    if (j %in% c(4,8)) mtext(y2lab, side = 4, line = 4,  cex=cex.lab)
    if (j %in% c(1,5)) mtext(ylab, side = 2, line = 3,  cex=cex.lab)
    if (j %in% c(5:8)) mtext(expression(paste("Target [", degree*C, ']')), line=4, padj=.1, side=1, cex=cex.lab)
  }  
}

multi.marginal.summary.plot <- function(location, loc.names3, mseason, myears, map, 
  model.getter=getLENSData, include.all.runs=FALSE, ...) {
  cex.lab = 1.7
  y2lab = expression(paste("Estimate [", degree*C, ']'))
  ylab = 'Density'
  par(oma=c(0,2,0,2))
  par(mfrow=c(2,4))
  for (j in 1:nloc) {   
    idx = which(loc.names3[j] == loc.names)
    print(loc.names3[j])
    y.f = model.getter(location[[idx]][1], location[[idx]][2], ...)[myears,]
    if (include.all.runs) {
      y.f.hat = xts(do.call(cbind, lapply(map[[idx]], function(x) as.matrix(x$rea.y.f))), order.by=time(map[[1]][[1]]$rea.y.f))  
    } else {
      y.f.hat = map[[idx]][[1]]$rea.y.f[,1]  
    }
    
    marginal.summary.plot(y.f[.indexmon(y.f) %in% mseason], 
      y.f.hat[.indexmon(y.f.hat) %in% mseason],
      main=loc.names3[j], cex.lab=1.8, cex.main=2.5, ylab="")
    if (j==1) legend('topleft', c('Target', 'Estimate'), lty=1, col=c(1,2), cex=1.3)
    if (j %in% c(4,8)) mtext(y2lab, side = 4, line = 3,  cex=cex.lab)
    if (j %in% c(1,5)) mtext(ylab, side = 2, line = 3,  cex=cex.lab)
  }  
}

# x appears in blue and y appears in red.
comp.histlike <- function(data, breaks=20, ...) {
  
  h.x = lapply(data, function(xi) hist(as.numeric(xi), breaks=breaks, plot=FALSE))   
  dots = list(...)
  dx = min(unlist(lapply(h.x, function(h.xi) {abs(h.xi$breaks[2]-h.xi$breaks[1])})))
  ylim = range(unlist(lapply(h.x, function(h.xi) h.xi$density)))
  u = range(unlist(lapply(h.x, function(h.xi) h.xi$breaks)))
  lwd = 1
  if (any(names(dots) %in% lwd)) lty=dots$lwd
  
  breaks = seq(u[1],u[2]+dx,dx)
  if (max(breaks)<max(ylim)) breaks = c(breaks,max(breaks)+dx)
    if ("plot" %in% names(dots) && dots$plot == FALSE) {
      hist(x,breaks=breaks, plot=FALSE) 
    } else if ("xlim" %in% names(dots)) {
      histlike(breaks=h.x[[1]]$breaks,values=h.x[[1]]$density,ylim=ylim,col="black", ...)
    } else {
      histlike(breaks=h.x[[1]]$breaks,values=h.x[[1]]$density,ylim=ylim,col="black", xlim=u, ...)
    }

  # abline(v=median(x), col="blue")
    k = 1
  lapply(h.x[-1], function(h.xi) {
    k <<- k + 1
    histlike(breaks=h.xi$breaks,values=h.xi$density, col=k, lwd=lwd, add=TRUE)
  })
  h.x
}


seasonal.marginal.summary.plot <- function(y, y.hat, season, season.names) {
  par(mfcol=c(2,2))
  par(mar=c(3,0,0,3))
  for (j in 1:length(season)) {
    par(mar=c(4, 5, 4, 5))
    marginal.summary.plot(y[.indexmon(y)%in%(season[[j]]),], 
      y.hat[.indexmon(y.hat)%in%(season[[j]])], main=season.names[j])
  }
}

marginal.summary.plot.multi <- function(data) {
  cex.points=.3
  cex.axis = 2
  comp.histlike(data, breaks=70,
    xlab="", cex.axis=cex.axis, cex.lab=1.8, ylab='')
  ndata = length(data)
  k=1
  for (i in 1:ndata) {
    sx = sort(as.numeric(data[[1]]))

    lenx = length(sx)
    if (k>1) {
      sy = sort(as.numeric(data[[i]]))
      leny = length(sy)
      if (leny < lenx) 
        sx <- approx(1L:lenx, sx, n = leny)$y
      if (leny > lenx) 
        sy <- approx(1L:leny, sy, n = lenx)$y  
    } 
    if (k == 2) {
      par(new=T)
      plot(sx, sy, axes=F, xlab=NA, ylab=NA, cex=cex.points, col=alpha(i, .5))
    } else if (k > 2) {
      points(sx, sy, cex=cex.points, col=alpha(i, .5))
    }
    k=k+1
  }
  par(xpd=FALSE)
  abline(0,1)
  abline(v=sx[floor(leny*0.01)], lty=2)
  axis(side = 4, cex.axis=cex.axis)
}

marginal.summary.plot <- function(y, y.hat, ...) {
  cex.axis = 2
  comp.histlike(list(y, y.hat), breaks=70,
    xlab="", cex.axis=cex.axis, ...)
  
  x = as.numeric(y)
  y = as.numeric(y.hat)
  sx <- sort(x)
  sy <- sort(y)
  lenx <- length(sx)
  leny <- length(sy)
  if (leny < lenx) 
      sx <- approx(1L:lenx, sx, n = leny)$y
  if (leny > lenx) 
      sy <- approx(1L:leny, sy, n = lenx)$y
  par(new=T)
  plot(sx, sy, axes=F, xlab=NA, ylab=NA, cex=.5, col=alpha(2, .5))
  par(xpd=FALSE)
  abline(0,1)
  abline(v=sx[floor(leny*0.01)], lty=2)
  axis(side = 4, cex.axis=cex.axis)
}

print.memory.usage <- function(min.size=1e6) {
  size = 0
  for (o in objects(envir=.GlobalEnv)) {
    size = size + object.size(get(o))
    if (object.size(get(o)) > min.size) {
        message(o); print(object.size(get(o)), units='auto')}
  }
  print(paste("Total memory usage:", size / 1000 / 1000, "Mb"))
}

get.ncdf.timefirst.ryan <- function(filename, lat_start, lat_count, lon_start, lon_count, t_start, t_count, varname='TREFHT') {
  ncdata <- nc_open(filename) 
  data = ncvar_get(ncdata, 
    varname, 
    start=c(t_start, lon_start, lat_start), 
    count=c(t_count, lon_count, lat_count))
  nc_close(ncdata)
  data
}


get.ncdf.timefirst <- function(filename, lat_start, lat_count, lon_start, lon_count, t_start, t_count, varname='TREFHT') {
  ncdata <- nc_open(filename) 
  data = ncvar_get(ncdata, 
    varname, 
    start=c(t_start, lon_start, lat_start), 
    count=c(t_count, lon_count, lat_count))
  nc_close(ncdata)
  data
}

qqplot <- function(x, y, ...) {
  x = as.numeric(x)
  y = as.numeric(y)
  sx <- sort(x)
  sy <- sort(y)
  lenx <- length(sx)
  leny <- length(sy)
  if (leny < lenx) 
      sx <- approx(1L:lenx, sx, n = leny)$y
  if (leny > lenx) 
      sy <- approx(1L:leny, sy, n = lenx)$y
 
  # plot(sx[1:200], sy[1:200])
  plot(sx, sy, ...)
  par(xpd=FALSE)
  abline(0,1)
}

qqplot2 <- function(y, gpd.fit) {
  n = length(y)
  z = gpd.fit$data
  nz = length(z)
  xz = 1:nz / (nz + 1)
  sig = gpd.fit$mle[1]
  sh = gpd.fit$mle[2]
  x = 1:n/(n+1)
  F_x = 1 - (1 + sh*y/sig)^(-1/sh)
  F_z = 1 - (1 + sh*z/sig)^(-1/sh)
  plot(- log(1 - xz), -log(1 - sort(F_z)), col=alpha(1, .6))
  points(- log(1 - x), -log(1 - sort(F_x)), , col=alpha(3, .6))
  abline(0,1, lwd=2, col=4)
}

# Take two lists and output one ist with entries alternating 
# between the first and last list: x_1, x_2, ... and y_1, y_2, ...
# becomes # x_1, y_1, x_2, y_2, ... etc
interlace <- function(x, y) {
  nx = length(x)
  ny = length(y)
  output = list(); length(output) = nx + ny
  output[seq(1,2*nx, by=2)] = x
  output[seq(2,2*ny, by=2)] = y
  output
}

get.quantile.change <- function(input, year.i=70, year.f=180) {
  yqhat_real = make.quantile.surfaces(input)  
  dy = apply(yqhat_real, c(1, 3), function(y) {
  y_i = y[year.i]
  y_f = y[year.f]
  dy = y_f - y_i
  })
  dy
}

plot.quantile.change <- function(input, year.i=70, year.f=180, ...) {
  yqhat_real = make.quantile.surfaces(input)  
  dy = apply(yqhat_real, c(1, 3), function(y) {
    y_i = y[year.i]
    y_f = y[year.f]
    dy = y_f - y_i
  })
  mdays = seq(1, 365, 20)
  q_all = input$q_all

  mywireframe(mdays, q_all, dy, ...)
}

mywireframe <- function(x, y, z, ...) {
  g <- expand.grid(x = x, y = y)
  g$z <- c(z[x,])
  wireframe(z ~ x * y, g, 
    ...)
}

# Plot the quantile map applied twice, with the last time applied in reverse.
# if the transform is the same you should get the identity transform or a zero surface.
plot.double.transform.all.inclusive <- function(input1, input2, year.i=70, year.f=170, logplot=TRUE, take.abs=FALSE, ...) {  
  ylab = list(label='Percentile', y=.5, rot=-40, cex=cex.lab)
  lw <- list(left.padding = list(x = -0.2, units = "inches"))
  lw$right.padding <- list(x = -0.2, units = "inches")
  lh <- list(bottom.padding = list(x = -.2, units = "inches"))
  lh$top.padding <- list(x = -.2, units = "inches")

  q_all = input1$q_all
  nfolds = length(input1)
  double.dy = matrix(0, 365, length(q_all))
  yqhat_real1 = make.quantile.surfaces(input1)  
  yqhat_real2 = make.quantile.surfaces(input2)  
  dy_i_f = apply(yqhat_real1, c(1, 3), function(y) {
    dy = y[year.f] - y[year.i]
  })
  dy_f_i = apply(yqhat_real2, c(1, 3), function(y) {
      dy = y[year.i] - y[year.f]
  })
  if (take.abs) {
    double.dy = double.dy + abs(dy_i_f + dy_f_i)
  } else {
    double.dy = double.dy + (dy_i_f + dy_f_i)
  }

  trellis.par.set("axis.line",list(col=NA,lty=1,lwd=1))
  p = mywireframe(mdays, q_all, double.dy[,length(q_all):1], 
        par.settings = list(axis.line = list(col = 'transparent')),
        lattice.options = list(layout.widths = lw, layout.heights = lh),
        ylab=ylab, 
         ...)
  p
}

# Plot the quantile map applied twice, with the last time applied in reverse.
# if the transform is the same you should get the identity transform or a zero surface.
plot.double.transform <- function(input1, input2, year.i=70, year.f=170, logplot=TRUE, take.abs=FALSE, ...) {  
  ylab = list(label='Percentile', y=.5, rot=-40, cex=cex.lab)
  lw <- list(left.padding = list(x = -0.2, units = "inches"))
  lw$right.padding <- list(x = -0.2, units = "inches")
  lh <- list(bottom.padding = list(x = -.2, units = "inches"))
  lh$top.padding <- list(x = -.2, units = "inches")

  q_all = input1[[1]]$q_all
  nfolds = length(input1)
  double.dy = matrix(0, 365, length(q_all))
  for (i in 1:(nfolds-1)) {
    for (j in (i+1):nfolds) {
      yqhat_real1 = make.quantile.surfaces(input1[[i]])  
      yqhat_real2 = make.quantile.surfaces(input2[[j]])  
      dy_i_f = apply(yqhat_real1, c(1, 3), function(y) {
        dy = y[year.f] - y[year.i]
      })
      dy_f_i = apply(yqhat_real2, c(1, 3), function(y) {
        dy = y[year.i] - y[year.f]
      })
      if (take.abs) {
        double.dy = double.dy + abs(dy_i_f + dy_f_i)
      } else {
        double.dy = double.dy + (dy_i_f + dy_f_i)
      }
    }    
  }  
  double.dy = double.dy / choose(nfolds,2)

  trellis.par.set("axis.line",list(col=NA,lty=1,lwd=1))
  p = mywireframe(mdays, q_all, double.dy[,length(q_all):1], 
        par.settings = list(axis.line = list(col = 'transparent')),
        lattice.options = list(layout.widths = lw, layout.heights = lh),
        ylab=ylab, 
         ...)
  p
}

plot.jackknife <- function(res) {
  cex.lab=1.5
  xlab = list(label='Days', y=.5, rot=30, cex=cex.lab)
  zlabel = expression(paste("|Error| [", degree*C, ']'))
  zlab = list(label=zlabel, x=.48, rot=90, cex=cex.lab)
  qticks = qticks = c(6, 8, 11, 14, 16)
  mdays = seq(1, 365, 10)
  
  ylab = list(label='Percentile', y=.5, rot=-40, cex=cex.lab)
  lw <- list(left.padding = list(x = -0.2, units = "inches"))
  lw$right.padding <- list(x = -0.2, units = "inches")
  lh <- list(bottom.padding = list(x = -.2, units = "inches"))
  lh$top.padding <- list(x = -.2, units = "inches")
  scales=list(arrows=FALSE, draw=TRUE, tck=.5, col=1,
            y=list(at=rev(q_all[qticks]), labels=(q_all[qticks])),
          z=list(at=c(0, .1, .2, .3, .4), labels=c("0", '0.1  ', '0.2  ', '0.3  ', '0.4  ')))
  mdays = seq(1, 365, 5)
  p = list(); length(p) = 8
  trellis.par.set("axis.line",list(col=NA,lty=1,lwd=1))

  for (j in 1:length(res)) {
    idx = which(loc.names3[j] == loc.names)
    p[[j]] = mywireframe(mdays, q_all, sqrt(res[[idx]]$v_jack[,rev(1:length(q_all))]), 
      par.settings = list(axis.line = list(col = 'transparent')),
      lattice.options = list(layout.widths = lw, layout.heights = lh),
      ylab=ylab, main=list(label=loc.names3[j], vjust=2.5),
      xlab=xlab, zlab=zlab, scales=scales, zlim=c(0, .45))
  }
  pp = arrangeGrob(grobs=p, ncol=4, nrow=2)
  pp
}

model.variance.jackknife <- function(input, year.i=70, year.f=170, mc.cores=1) {
  n = length(input)
  yqhats = mclapply(input, make.quantile.surfaces, mc.cores=mc.cores)
  nq = dim(yqhats[[1]])[3]
  t_j = lapply(yqhats, function(l) {l[,year.f,] - l[,year.i,]})
  l_jack = lapply(t_j[-1], function(t_ji) (n - 1) * (t_j[[1]] - t_ji))
  b_jack = - apply(array(unlist(l_jack), dim=c(365, nq, n - 1)), c(1,2), sum) / n
  v_i = lapply(l_jack, function(l_jack_j) (l_jack_j^2 - b_jack^2))
  v_jack = apply(array(unlist(v_i), dim=c(365, nq, n - 1)), c(1,2), sum) / (n*(n-1))

  res = list(b_jack=b_jack, v_jack=v_jack)
  
}

# Plot the quantile map applied twice, with the last time applied in reverse.
# if the transform is the same you should get the identity transform or a zero surface.
plot.model.difference.all.inclusive <- function(input1, input2, year.i=70, year.f=170, mq=c(.1, .5, .9), ...) {
  zlabel = expression(paste("Change [", degree*C, ']'))
  zlab = list(label=zlabel, x=.48, rot=90, cex=cex.lab)
  ylab = list(label='Years', y=.5, rot=-40, cex=cex.lab)
  xlab = list(label='Days', y=.5, rot=30, cex=cex.lab)

  q_all = input1$q_all
  nfolds = length(input1)
  nfolds = 2
  double.dy = array(0, dim=c(365, 180,length(q_all)))
  yqhat_real1 = make.quantile.surfaces(input1)  
  yqhat_real2 = make.quantile.surfaces(input2)    
  double.dy = double.dy + (yqhat_real2 - yqhat_real1)
  all_equal = unlist(lapply(q_all, function(q) any(abs(q-mq)<0.0001)))
  zlim = range(double.dy[,, all_equal])
  zlim2 = range(yqhat_real2[,, all_equal])
  mdays = seq(1, 365, 5)
  myears = seq(1,180, 10)
  
  p = list(); length(p) = length(mq)
  for (i in 1:length(mq)) {
    z = double.dy[,,abs(q_all-mq[i])<.00001]
    # z = yqhat_real2[,,q_all==mq[i]]
    p[[i]] = mywireframe(mdays, myears, z[,myears], 
        scales=list(arrows=FALSE, col=1),
        zlab=zlab,
        main=list(label=paste("Percentile=", mq[i], sep=''), vjust=2.5),
        par.settings = list(axis.line = list(col = 'transparent')),
        ylab=ylab, 
        xlab=xlab
        ,zlim=zlim
        )
  }
  p
}



# Plot the quantile map applied twice, with the last time applied in reverse.
# if the transform is the same you should get the identity transform or a zero surface.
plot.model.difference <- function(input1, input2, year.i=70, year.f=170, mq=c(.1, .5, .9), ...) {
  zlabel = expression(paste("Change [", degree*C, ']'))
  zlab = list(label=zlabel, x=.48, rot=90, cex=cex.lab)
  ylab = list(label='Years', y=.5, rot=-40, cex=cex.lab)
  xlab = list(label='Days', y=.5, rot=30, cex=cex.lab)

  q_all = input1[[1]]$q_all
  nfolds = length(input1)
  nfolds = 2
  double.dy = array(0, dim=c(365, 180,length(q_all)))
  for (i in 1:(nfolds-1)) {
    for (j in (i+1):nfolds) {
      print(c(i,j))
      yqhat_real1 = make.quantile.surfaces(input1[[i]])  
      yqhat_real2 = make.quantile.surfaces(input2[[j]])    
      double.dy = double.dy + (yqhat_real2 - yqhat_real1)
    }    
  }  
  all_equal = unlist(lapply(q_all, function(q) any(abs(q-mq)<0.0001)))
  double.dy = double.dy / choose(nfolds,2)
  zlim = range(double.dy[,, all_equal])
  zlim2 = range(yqhat_real2[,, all_equal])
  mdays = seq(1, 365, 5)
  myears = seq(1,180, 10)
  
  p = list(); length(p) = length(mq)
  for (i in 1:length(mq)) {
    z = double.dy[,,abs(q_all-mq[i])<.00001]
    # z = yqhat_real2[,,q_all==mq[i]]
    p[[i]] = mywireframe(mdays, myears, z[,myears], 
        scales=list(arrows=FALSE, col=1),
        zlab=zlab,
        main=list(label=paste("Percentile=", mq[i], sep=''), vjust=2.5),
        par.settings = list(axis.line = list(col = 'transparent')),
        ylab=ylab, 
        xlab=xlab
        ,zlim=zlim
        )
  }
  p
}


plot.surface <- function(z) {
  dims = dim(z)
  x = seq(1, dims[1], len=50)
  y = seq(1, dims[2], len=50)
  persp(x=x,y=y,z=z[x,y], phi=45, theta=25)
}

# Returns a matrix whose columns are the quantile fits specified
# by the 'input' parameters.
make.quantile.surfaces <- function(input) {
  norm.x.df = input$norm.x.df
  bulk.x.df = input$bulk.x.df
  tail.x.df = input$tail.x.df
  mlat = input$mlat
  year.range = input$year.range
  nyears = diff(year.range) + 1
  norm.x.one = getPredictors(n_files=1, df.x=norm.x.df[1], df.t=norm.x.df[2], df.xt=norm.x.df[3], year.range=year.range, get.volc=T, mlat=mlat)
  bulk.x.one = getPredictors(n_files=1, df.x=bulk.x.df[1], df.t=bulk.x.df[2], df.xt=bulk.x.df[3], year.range=year.range)
  tail.x.one = getPredictors(n_files=1, df.x=tail.x.df[1], df.t=tail.x.df[2], df.xt=tail.x.df[3], year.range=year.range)
  p.norm = dim(norm.x.one)[2]
  yqhat = norm.x.one %*% input$coef_norm
  yqhat_bulk = bulk.x.one %*% input$coef_bulk
  yqhat_bulk_real = yqhat_bulk * (yqhat[,3] - yqhat[,1]) + yqhat[,2]
  yqhat_low = tail.x.one %*% input$coef_tail[[1]]
  yqhat_low_real = (yqhat_low + yqhat_bulk[,1]) * (yqhat[,3] - yqhat[,1]) + yqhat[,2]
  yqhat_high = tail.x.one %*% input$coef_tail[[2]]
  yqhat_high_real = (yqhat_high + yqhat_bulk[,length(input$q)]) * (yqhat[,3] - yqhat[,1]) + yqhat[,2]
  nq_tot = dim(yqhat_high)[2] + dim(yqhat_low)[2] + dim(yqhat_bulk)[2]
  yqhat_real = cbind(yqhat_low_real, yqhat_bulk_real, yqhat_high_real)
  yqhat_real = array(yqhat_real, dim=c(365, nyears, nq_tot))
  yqhat_real
}

plot.quantile.fit <- function(yi, input, tau=0.5, plotit=TRUE, normit=TRUE, by.year=FALSE, ...) {
  cex.lab=2
  my_yqhat = make.quantile.surfaces(input)
  nq_tot = dim(my_yqhat)[3]
  q_all = input$q_all
  t = time(yi)
  p = dim(yi)[2]
  year = as.numeric(strftime(t, "%Y"))
  year.range = range(year)
  unique.year = unique(year)
  nyears = diff(year.range) + 1
  model.year.i = input$year.range[1]
  year.i = min(year)
  year.f = max(year)
  nruns = length(yi) / 365 / nyears
  doy = as.numeric(strftime(t, "%j"))
  m = my_yqhat[,,q_all==.5]
  r = my_yqhat[,,q_all==q_norm[3]] - my_yqhat[,,q_all==q_norm[1]] 
  
  if (normit) {
    bulk.x.df = input$bulk.x.df
    tail.x.df = input$tail.x.df
    bulk.x.one = getPredictors(n_files=1, df.x=bulk.x.df[1], df.t=bulk.x.df[2], df.xt=bulk.x.df[3], year.range=input$year.range)
    tail.x.one = getPredictors(n_files=1, df.x=tail.x.df[1], df.t=tail.x.df[2], df.xt=tail.x.df[3], year.range=input$year.range)
    bulk = bulk.x.one %*% input$coef_bulk
    diff_hat = array(cbind(tail.x.one %*% input$coef_tail[[1]] + bulk[,1],
                          bulk,
                          tail.x.one %*% input$coef_tail[[2]] + bulk[,length(q)]),
                  dim=c(365, diff(input$year.range) + 1, nq_tot))
  } else {
    r = array(1,dim=c(365,nyears))
    m = array(0,dim=c(365,nyears))
    diff_hat = apply(my_yqhat, 3, function(z) (z - m) / r)
    diff_hat = array(diff_hat, dim=c(365, nyears, nq_tot)) 
  }

  diff = as.numeric(t(apply(cbind(yi, doy, year - model.year.i + 1), 1, function(z) {
        doy_i = z[p+1]
        year_i = z[p+2]
        (z[1:p] - m[doy_i, year_i]) / r[doy_i, year_i]
      })))
  
  ny = 1:length(yi)
  idx = sample(ny, min(10000, length(ny)))
  q_idx = lapply(tau, function(mtau) round(q_all - mtau, 3) == 0)
  if (plotit){
    if (by.year) {
      plot(rep(year, nruns)[idx], diff[idx], col=rgb(0,0,0,.2),
        xlab="Years Elapsed", cex.lab=cex.lab, ...)   
      for (qi in q_idx) {
        lines(unique.year, diff_hat[1, unique.year - year.i + 1, qi], col=2, lwd=2)
        lines(unique.year, diff_hat[182, unique.year - year.i + 1, qi], col=3, lwd=2, lty=1)
      }
      leg.pos = ifelse(tau[1]<.5, 'bottomleft', 'topleft')
      legend(leg.pos, c('Jan1', 'July1'), col=c(2,3), lty=1, cex=1.5)
    
    } else {
      plot(rep(doy, nruns)[idx], diff[idx], col=rgb(0,0,0,.2),
        xlab="Day of the year", cex.lab=cex.lab, ...)
      for (qi in q_idx) {
        lines(1:365, diff_hat[, 1, qi], col=2, lwd=2)
        lines(1:365, diff_hat[, nyears, qi], col=3, lwd=2, lty=1)
      }
      leg.pos = ifelse(tau[1]<.5, 'bottomleft', 'topleft')
      legend(leg.pos, c('1920', '2099'), col=c(2,3), lty=1, cex=1.5)
    
    }
  }
  tauhat = lapply(q_idx, function(qi) {
    sum(apply(cbind(yi, doy, year - year.i + 1), 1, function(z) {
      doy_i = z[p+1]
      year_i = z[p+2]
      z_norm = (z[1:p] - m[doy_i, year_i]) / r[doy_i, year_i]
      sum(z_norm < diff_hat[doy_i, year_i, qi])
    })) / length(t) / p
  })
}

getReducedPredictors <- function(n_files, df.x=15, df.t=4, df.xt=3, 
  year.range=c(1850,2099), get.volc=FALSE, mlat=30, coef=NULL) {
  x = getPredictors(n_files=n_files, df.x=df.x, df.t=df.t, df.xt=df.xt, year.range=year.range, get.volc=get.volc, mlat=mlat)
  reduce_predictors(x, coef)
}

getPredictors <- function(n_files, df.x=14, df.t=6, df.xt=3, 
  year.range=c(1850,2099), get.volc=FALSE, mlat=30) {
  day_of_year = 1:365
  co2 = read.table("RCP85_MIDYR_CONC.DAT", skip=38, header=TRUE)
  b2 = co2[co2$YEARS>=year.range[1] & co2$YEARS<=year.range[2], 2]
  nyears = diff(year.range) + 1
  x = scale(rep(day_of_year, n_files*nyears), center=FALSE)
  x.main.basis = as.matrix(pbs(x, df=df.x))
  t = scale(rep(c(t(matrix(rep(1:nyears, 365), ncol=365))), n_files), center=FALSE)
  if (df.t==-1) {
    b2.trend = rep(c(t(matrix(rep(b2, 365), ncol=365))), n_files)
    b2.trend.squared = rep(c(t(matrix(rep(b2^2, 365), ncol=365))), n_files)
    # t.basis = scale(cbind(b2.trend, b2.trend.squared), center=FALSE)
    t.basis = scale(cbind(t, b2.trend, b2.trend.squared), center=FALSE)
  } else {
    
    t.basis = ns(t, df=df.t)
  }
  if (df.xt != 0) {
    x.int.basis = scale(as.matrix(pbs(x, df=df.xt)), center=FALSE)
    if (df.t==-1) {
      # t = rep(c(t(matrix(rep(1:nyears, 365), ncol=365))), n_files)
      # t.basis.int = cbind(t.basis, ns(t, df=1))
      # t.basis.int = t.basis
      X = model.matrix(~x.main.basis + t.basis + x.int.basis:t.basis)
    } else {
      X = model.matrix(~x.main.basis + t.basis + x.int.basis:t.basis)
    }
  } else {
    X = model.matrix(~x.main.basis + t.basis)
  }
  if (get.volc) {
      volc.time = c(t(matrix(rep(seq(1850, 2099), 365), 250, 365)))
      volc = scale(rep(get.volcanic(mlat)[volc.time %in% seq(year.range[1], year.range[2])], n_files), center=FALSE)
      X = cbind(X, volc)
    }
  rm(x, x.main.basis)
  X
}

fit_rq <- function(x, y, tau=.9, method='lasso', plotit=FALSE,
  lambda=exp(seq(-1,7, len=40))) {
  
  out = list(); length(out) = length(lambda)
  for (i in 1:length(lambda)) {
    out[[i]] = rq(as.numeric(y)~x-1, tau=tau, method=method, lambda = lambda[i])
  }
  coefs = do.call(cbind, lapply(out, function(x) x$coef))

  if (plotit) {
    plot(coefs[2,], ylim=range(coefs), type='l')
    apply(coefs[-c(1, 2),], 1, function(x) lines(x))
  }

  coefs
}

plot.rq.fit <- function(dayofyear, y, fit) {
  year.range = range(as.numeric(strftime(time(y), format = "%Y")))
  shift = 365*1
  shift2 = 365*(diff(year.range))
  yhat_high = fit$fit[1:365+shift,3]
  yhat_mid = fit$fit[1:365+shift,2]
  yhat_low = fit$fit[1:365+shift,1]
  plot(dayofyear, rea.y)
  lines(dayofyear[1:365], yhat_low, col=2, lwd=2, lty=2)
  lines(dayofyear[1:365], yhat_high, col=2, lwd=2, lty=2)
  lines(dayofyear[1:365], yhat_mid, col=2, lwd=2, lty=1)
  lines(dayofyear[1:365], fit$fit[1:365+shift2,1], col=3, lwd=2, lty=2)
  lines(dayofyear[1:365], fit$fit[1:365+shift2,2], col=3, lwd=2, lty=1)
  lines(dayofyear[1:365], fit$fit[1:365+shift2,3], col=3, lwd=2, lty=2)
}


# x appears in blue and y appears in red.
comp.hist <- function(x,y, breaks=20, ...) {
  h.x = hist(x, breaks=breaks, plot=FALSE)
  h.y = hist(y, breaks=breaks, plot=FALSE)
  dots = list(...)
  dx = min(abs(h.x$breaks[2]-h.x$breaks[1]),
           abs(h.y$breaks[2]-h.y$breaks[1]))

  ylim = range(c(h.x$density,h.y$density))
  u = range(c(h.x$breaks,h.y$breaks))
  breaks = seq(u[1],u[2]+dx,dx)
  if (max(breaks)<max(ylim)) breaks = c(breaks,max(breaks)+dx)
    if ("plot" %in% names(dots) && dots$plot == FALSE) {
      h.x = hist(x,breaks=breaks, plot=FALSE) 
    } else if ("xlim" %in% names(dots)) {
      h.x = hist(x,breaks=breaks,density=15,ylim=ylim,col="black",freq=FALSE, ...)
    } else {
      h.x = hist(x,breaks=breaks,density=15,ylim=ylim,col="black",xlim=u,freq=FALSE, ...)
    }

  # abline(v=median(x), col="blue")
  h.y = hist(y,breaks=breaks,density=15,col="red",add=TRUE,freq=FALSE,angle=-45)
  list(h.x=h.x, h.y=h.y)
}


exceedence_loss_1 <- function(y.te, y.te.hat) {
  n = dim(y.te)[1]
  yi.hat = c(y.te.hat)
  tauhat = apply(y.te, 2, function(yi) {
    sum(yi.hat > yi) / n
  })
  tauhat
}


exceedence_loss_i <- function(y.te, y.te.hat, nbins.days=10, nbins.years=5) {
  n = dim(y.te)[1]
  nruns = dim(y.te)[2]
  yi.hat = c(y.te.hat)
  time.day = as.numeric(strftime(time(y.te), format = "%j"))
  time.year = as.numeric(strftime(time(y.te), format = "%Y"))
  tauhat = apply(y.te, 2, function(yi) {
    res = yi - yi.hat
    h = hist2d(time.day[res<0], time.year[res<0], nbins=c(nbins.days, nbins.years), show=FALSE)
    h$count / (n/nbins.days/nbins.years)
  })
  tauhat
}


exceedence_loss <- function(y.te, y.te.hat, nbins.days=10, nbins.years=5) {
  n = length(y.te)
  nruns = dim(y.te)[2]
  yi.hat = c(y.te.hat)
  if (nruns==1 || is.null(dim(y.te))) { 
    res = y.te - c(y.te.hat)
  } else {
    res = c(apply(y.te, 2, function(yi) yi - yi.hat))
  }
  time.day = rep(as.numeric(strftime(time(y.te), format = "%j")), nruns)
  time.year = rep(as.numeric(strftime(time(y.te), format = "%Y")), nruns)
  h = hist2d(time.day[res<0], time.year[res<0], nbins=c(nbins.days, nbins.years), show=FALSE)
  tauhat = h$count / (n/nbins.days/nbins.years)
  tauhat
}

get.map.transform <- function(initial_value, q_initial, q_final, idx) {
  ql = q_initial[idx]
  qh = q_initial[idx+1]
  q_frac = (initial_value - ql) / (qh - ql)
  # my_quantile = q[idx] + q_frac * (q[idx+1] - q[idx])
  final_transform = q_final[idx] + q_frac * (q_final[idx+1] - q_final[idx])
  final_transform
}

model_qqstats <- function(map, y, indexmon=seq(0,11)) {
  k = 0
  idx = .indexmon(map[[1]]$rea.y.f) %in% (indexmon ) 
  output = lapply(map, function(mapout) {
    k <<- k + 1
    apply(mapout$rea.y.f[idx,], 2, function(yi.hat) {
        qqstats(as.numeric(y), as.numeric(yi.hat), standardize=FALSE)  
    })
   })
  qqstats = apply(do.call(rbind, do.call(rbind, output)), 2, unlist)
}

map.data <- function(params, rea.x.df, rea.y, d.year=20, identity.transform=FALSE, 
  const.projection.year=NULL, nested=FALSE) {
  # const.projection.year: if  not NULL this will project every year into the same target year.
  # Otherwise, d.year will be used to project every year a fixed amount.

  q = params$q
  nq = length(params$q)
  q_norm = params$q_norm
  q_low = params$q_tail
  mlat = params$mlat
  q_high = rev(1-params$q_tail)
  year.range = params$year.range
  nyears = diff(year.range) + 1
  norm.x.df = params$norm.x.df
  bulk.x.df = params$bulk.x.df
  tail.x.df = params$tail.x.df
  if (nested) {
    norm.x.one = getReducedPredictors(n_files=1, df.x=norm.x.df[1], df.t=norm.x.df[2], 
      df.xt=norm.x.df[3], year.range=year.range, get.volc=T, coef=params$coef_norm)
    bulk.x.one = getReducedPredictors(n_files=1, df.x=bulk.x.df[1], df.t=bulk.x.df[2], 
      df.xt=bulk.x.df[3], year.range=year.range, coef=params$coef_bulk)
    low.tail.x.one = getReducedPredictors(n_files=1, df.x=tail.x.df[1], df.t=tail.x.df[2], 
      df.xt=tail.x.df[3], year.range=year.range, coef=params$coef_tail[[1]])
    high.tail.x.one = getReducedPredictors(n_files=1, df.x=tail.x.df[1], df.t=tail.x.df[2], 
      df.xt=tail.x.df[3], year.range=year.range, coef=params$coef_tail[[2]])
    p.norm = dim(norm.x.one[[1]])[2]

    yqhat = array(list_matrix_multiply(norm.x.one, params$reduced_coef_norm), dim=c(365, nyears, length(q_norm)))
    yqhat_bulk = array(list_matrix_multiply(bulk.x.one, params$reduced_coef_bulk), dim=c(365, nyears, length(q)))
    yqhat_low = array(list_matrix_multiply(low.tail.x.one, params$reduced_coef_tail[[1]]) + c(yqhat_bulk[,,1]), dim=c(365, nyears, length(q_low)))
    yqhat_high = array(list_matrix_multiply(high.tail.x.one, params$reduced_coef_tail[[2]]) + c(yqhat_bulk[,,length(q)]), dim=c(365, nyears, length(q_high)))
      
  } else {
    norm.x.one = getPredictors(n_files=1, df.x=norm.x.df[1], df.t=norm.x.df[2], df.xt=norm.x.df[3], year.range=year.range, get.volc=T)
    bulk.x.one = getPredictors(n_files=1, df.x=bulk.x.df[1], df.t=bulk.x.df[2], df.xt=bulk.x.df[3], year.range=year.range)
    tail.x.one = getPredictors(n_files=1, df.x=tail.x.df[1], df.t=tail.x.df[2], df.xt=tail.x.df[3], year.range=year.range)
    p.norm = dim(norm.x.one)[2]
  
    yqhat = array(norm.x.one %*% params$coef_norm, dim=c(365, nyears, length(q_norm)))
    yqhat_bulk = array(bulk.x.one %*% params$coef_bulk, dim=c(365, nyears, length(q)))
    yqhat_low = array(tail.x.one %*% params$coef_tail[[1]] + c(yqhat_bulk[,,1]), dim=c(365, nyears, length(q_low)))
    yqhat_high = array(tail.x.one %*% params$coef_tail[[2]] + c(yqhat_bulk[,,length(q)]), dim=c(365, nyears, length(q_high)))
  }
  nq_low = length(q_low)
  nq_high = length(q_high)
  nq_tot = nq + nq_low + nq_high

  t = time(rea.y)
  doy = as.numeric(strftime(t, format = "%j"))
  rea.year = as.numeric(strftime(t, format = "%Y"))
  min.rea.year = min(rea.year)
  nyears = diff(range(rea.year)) + 1
  p = dim(rea.y)[2]
  rea.x = getPredictors(n_files=1, df.x=rea.x.df[1], df.t=rea.x.df[2], 
    df.xt=rea.x.df[3], year.range=range(rea.year), get.volc=TRUE, mlat=mlat)
  rea.coef = apply(rea.y, 2, function(rea.yi) {      
      list(rq(as.numeric(rea.yi)~rea.x-1, tau=q_norm)$coef)
    })
  yqhat_rea = do.call(rbind, lapply(rea.coef, function(coef.i) {
      rea.x %*% coef.i[[1]]
    }))

  m = yqhat_rea[, 2]
  r = yqhat_rea[, 3] - yqhat_rea[, 1]
  rea.y.norm = (rea.y - m) / r
  m_model = yqhat[, , q_norm==0.5]
  r_model = yqhat[, , 3] - yqhat[, , 1]
  
  loop.matrix = as.matrix(cbind(as.numeric(rea.y.norm), doy, r, m, rea.year - year.range[1] + 1)) 
  rea_final = apply(loop.matrix, 1, function(z) {
    initial_value = z[1]
    d = z[2]
    r_d = z[3]
    m_d = z[4]
    y.i = z[5]
    if (!is.null(const.projection.year)) {
      y.f = const.projection.year - year.range[1] + 1  
    } else {
      y.f = y.i + d.year  
    }
    
    q_initial = c(yqhat_low[d, y.i, ], yqhat_bulk[d, y.i, ], yqhat_high[d, y.i, ])
    q_final = c(yqhat_low[d, y.f, ], yqhat_bulk[d, y.f, ], yqhat_high[d, y.f, ])
    nq_tot = length(q_initial)
    idx = min(max(which(initial_value > q_initial), 1), nq_tot)
    idx.at.edge = FALSE
    if (idx==1) {
      idx.at.edge = TRUE
    } else if (idx == nq_tot) {
      idx.at.edge = TRUE
    }
    if (idx.at.edge) {
      # final_transform = q_final[idx]
      final_transform = initial_value + q_final[idx] - q_initial[idx]
    } else {
      final_transform = get.map.transform(initial_value, q_initial, q_final, idx)
    }
    if (identity.transform) {
      final_value = r_d * r_model[d, y.f] / r_model[d, y.i] * initial_value + m_model[d, y.f] - m_model[d, y.i] + m_d
    } else {
      final_value = r_d * r_model[d, y.f] / r_model[d, y.i] * final_transform + m_model[d, y.f] - m_model[d, y.i] + m_d
    }
    # final_value = r_d * final_transform + m_f[d] - m_i[d] + m_d
    
    list(final_value=final_value, final_transform=final_transform, idx=idx)
  })

  rea.y.f = xts(matrix(unlist(lapply(rea_final, function(x) x$final_value)), length(doy), p), order.by=t)
  idx = unlist(lapply(rea_final, function(x) x$idx))
  
  transform = unlist(lapply(rea_final, function(x) x$final_transform))
  n_lowtails = sum(idx==1)
  n_hightails = sum(idx==nq_tot)

  list(rea.y.f=rea.y.f, rea.coef=rea.coef, 
    n_tails=c(n_lowtails, n_hightails), rea.x.df=rea.x.df, d.year=d.year)
}

percentile.diff <- function(y.hat,y, percentiles=c(0.01, 0.1, 0.5, 0.9, 0.99)) {
  x = as.numeric(y)
  y = as.numeric(y.hat)
  sx <- sort(x)
  sy <- sort(y)
  lenx <- length(sx)
  leny <- length(sy)
  if (leny < lenx) 
      sx <- approx(1L:lenx, sx, n = leny)$y
  if (leny > lenx) 
      sy <- approx(1L:leny, sy, n = lenx)$y
  diff = sx - sy
  lenx = length(sx)
  result = unlist(lapply(percentiles, function(per) {diff[lenx*per]}))
  result
}

get.percentile.diff <- function(y.f, y.f.hat) {
  diff = lapply(seq(length(location)), function(i) {
    lapply(seq(length(season)), function(j) {
      apply(y.f.hat[[i]][[j]], 2, percentile.diff, y.f[[i]][[j]])    
    })
  })
  dfloss = melt(diff)
  colnames(dfloss) = c('Percentile', 'Run', 'Loss', 'Season', 'City')
  dfloss[,'Percentile'] = as.factor(c(0.01, 0.1, 0.5, 0.9, 0.99)[dfloss[,'Percentile']])
  dfloss[,'Season'] = as.factor(season.names)[dfloss[,'Season']]
  dfloss[,'City'] = factor(loc.names, levels=loc.names)[dfloss[,'City']]
  dfloss
}

multi.quantile.loss <- function(qout, map, location) {
  # Make boxplots of the loss separated by Winter/Summer, Mid/Late periods,
  # and quantiles
  season.days = list(c(335:365, 1:59), c(60:151), c(152:242), c(243:334))
  mid.range = seq(2019, 2056)
  late.range = seq(2059, 2096)
  season.names = c('Winter', 'Spring', 'Summer', 'Autumn')
  period = list(mid.range, late.range)
  season = list(c(11, 0, 1), c(2,3,4), c(5, 6, 7), c(8,9,10))
  all.days = 1:365; all.years = 1:180;
  tseq= timeBasedSeq('19200101::20991231')
  tseq = tseq[strftime(tseq, format = "%j")!='366'] # remove leep years
  nfolds = length(map[[1]])
  p = dim(map[[1]][[1]]$rea.y.f)[2] * nfolds
  time.month = as.numeric(strftime(tseq, format = "%m")) - 1
  time.year = as.numeric(strftime(tseq, format = "%Y"))
  temploss = array(0, dim=c(length(season), length(map), p, length(q_all)))
  qqstats = array(0, dim=c(length(season), length(map), p, 3))
  for (i in 1:length(season)) {
    print(i)
    for (j in 1:length(location)) {
      idx = time.month %in% season[[i]] & time.year %in% period[[2]]    
      temploss[i,j,,] = do.call(rbind, lapply(seq(nfolds), function(k) {
        y.hat = make.quantile.surfaces(qout[[j]][[k]])
        l = 0
        y.hat.mean = apply(y.hat, 3, function(y.hat.q) mean(c(y.hat.q)[idx]))
        loss.ij = apply(y.hat, 3, function(y.hat.q) {
          l <<- l + 1
          y.hat.sub = c(y.hat.q)[idx]
          y = map[[j]][[k]]$rea.y.f
          y.sub = y[.indexmon(y) %in% season[[i]],] # Test on projected out-of-sample data
          loss.q = exceedence_loss_1(y.sub, y.hat.sub) - q_all[l]
          dq2temp(loss.q, y.hat.mean, q_all[l], q_all)
        })    
        loss.ij
        }))   
      
    }
  }
    # Loss vs. Quantile
  dfloss = melt(temploss, varnames=c('Season', 'City', 'Run', 'Percentile'),
  value.name='Loss')
  dfloss[,'Percentile'] = as.factor(q_all[dfloss[,'Percentile']])
  dfloss[,'Season'] = as.factor(season.names)[dfloss[,'Season']]
  dfloss[,'City'] = factor(loc.names, levels=loc.names)[dfloss[,'City']]
  dfloss = dfloss[dfloss[,'Percentile'] %in% c(0.01, .1, .5, .9, 0.99), ]

  dfloss
}

dq2temp <- function(mloss, f, mq, q_all) {
  q_idx = which(q_all == mq)
  df = diff(f)
  nq = length(q_all)
  df.dq = df / diff(q_all)

  if (q_idx != nq && q_idx != 1) {
    temp =  (df.dq[q_idx] * (mloss > 0) + df.dq[q_idx - 1] * (mloss < 0)) * mloss
  } else if (q_idx == 1) {
    temp = df.dq[q_idx] * mloss
  } else {
    temp = df.dq[q_idx-1] * mloss
  }
  temp
}

plot.rea.fit <- function(params, rea.y, df.x=10, normit=TRUE) {
  rea.x.df = params$rea.x.df
  rea.coef = params$rea.coef
  q_all = params$q_all
  t = time(rea.y)
  rea.x = getPredictors(n_files=1, df.x=df.x, df.t=rea.x.df[2], 
    df.xt=rea.x.df[3], year.range=range(rea.year), get.volc=TRUE, mlat=mlat)
  rea.coef = apply(rea.y, 2, function(rea.yi) {
      list(rq(as.numeric(rea.yi)~rea.x-1, tau=q_norm)$coef)
    })
  yqhat_rea = lapply(rea.coef, function(coef.i) {
      rea.x %*% coef.i[[1]]
    })
  if (normit) {
    y.hat = do.call(rbind, yqhat_rea)
    y = as.numeric(rea.y)
    y.norm = (y - y.hat[,2]) / (y.hat[,3] - y.hat[,1])
    rea.y = xts(matrix(y.norm, length(t), dim(rea.y)[2]), order.by=t)
  }
  doy = as.numeric(strftime(time(rea.y), '%j'))
  p = dim(rea.y)[2]
  n = length(rea.y)
  idx = sample(seq(n), 50000)
  plot(rep(doy, p)[idx], as.numeric(rea.y)[idx])
  lines(yqhat_rea[[2]][1:365, 2], col=2, lwd=2)
  lines(yqhat_rea[[2]][1:365 + (365*37), 2], col=2, lwd=2) 
}

implode <- function(..., sep='') {
     paste(..., collapse=sep)
}

rq.fit.pfn.mult <- function(x, y, tau = 0.5,  Mm.factor = 0.8,
  max.bad.fixup = 3, eps = 1e-6, alpha=1)
# If y is a matrix assume that we have repeated sets of predictors
{
  #rq function for n large --
  py = dim(y)[2]
  n <- length(y)
  y <- as.numeric(y)
  n_one_y = n / py
  if(nrow(x) != n)
    stop("x and y don't match n")
  if(tau < 0 | tau > 1)
    stop("tau outside (0,1)")
  p <- ncol(x)
  m <- py * alpha * round(((p + 1) * n / py)^(2/3))
  not.optimal <- TRUE
  while(not.optimal) {
    if(m < n)
      s <- sample(n, m)
    else {
      b <- rq.fit.fnb(x, y, tau = tau,  eps = eps)$coef
      break
    }
    xx <- x[s,  ]
    yy <- y[s]
    z <- rq.fit.fnb(xx, yy, tau = tau,  eps = eps)
    xxinv <- solve(chol(crossprod(xx)))
    band <- sqrt(((x %*% xxinv)^2) %*% rep(1, p))
    #sqrt(h<-ii)
    r <- y - x %*% z$coef
    M <- Mm.factor * m
    lo.q <- max(1/n, tau - M/(2 * n))
    hi.q <- min(tau + M/(2 * n), (n - 1)/n)
    kappa <- quantile(r/pmax(eps, band), c(lo.q, hi.q))
    sl <- r < band * kappa[1]
    su <- r > band * kappa[2]
    bad.fixup <- 0
    while(not.optimal & (bad.fixup < max.bad.fixup)) {
      xx <- x[!su & !sl,  ]
      yy <- y[!su & !sl]
      if(any(sl)) {
        glob.x <- c(t(x[sl,  , drop = FALSE]) %*% rep(
          1, sum(sl)))
        glob.y <- sum(y[sl])
        xx <- rbind(xx, glob.x)
        yy <- c(yy, glob.y)
      }
      if(any(su)) {
        ghib.x <- c(t(x[su,  , drop = FALSE]) %*% rep(
          1, sum(su)))
        ghib.y <- sum(y[su])
        xx <- rbind(xx, ghib.x)
        yy <- c(yy, ghib.y)
      }
      z <- rq.fit.fnb(xx, yy, tau = tau,  eps = eps)
      b <- z$coef
      r <- y - x %*% b
      su.bad <- (r < 0) & su
      sl.bad <- (r > 0) & sl
      if(any(c(su.bad, sl.bad))) {
        if(sum(su.bad | sl.bad) > 0.10000000000000001 *
          M) {
          warning("Too many fixups:  doubling m")
          m <- 2 * m
          break
        }
        su <- su & !su.bad
        sl <- sl & !sl.bad
        bad.fixup <- bad.fixup + 1
      }
      else not.optimal <- FALSE
    }
  }
  coefficients <- b
  names(coefficients) <- dimnames(x)[[2]]
  residuals <- y - x %*% b
  return(list(coefficients=coefficients, tau=tau,
    residuals=residuals))
}

# This function will perform the quantile mapping. rea.y has to be a 
# time series object, like xts. The other x and y can be arrays and matrices.
# Give q_tail as the fractions of exceedences beyond the upper and lower 
# regular percentile, given in the vector q.
quantile.map <- function(norm.x.df, bulk.x.df, model.y, tail.x.df, q, q_tail, q_norm, mlat=30) {
  nq = length(q)
  n_files = dim(model.y)[2]
  q_low = q_tail
  q_high = rev(1-q_tail)
  q_all = c(q_low*q[1], q, q[length(q)] + q_high*(1-q[length(q)]))
  t = time(model.y)
  doy = as.numeric(strftime(t, format = "%j"))
  year = as.numeric(strftime(t, format = "%Y"))
  min.year = min(year)
  nyears = diff(range(year)) + 1
  year.range=range(year)

  # NORM FIT
  norm.x = getPredictors(n_files=n_files, df.x=norm.x.df[1], 
    df.t=norm.x.df[2], df.xt=norm.x.df[3], 
    year.range=year.range, get.volc=TRUE, mlat=mlat)    
  print('Fitting Norm')
  coef_norm = matrix(0, dim(norm.x)[2], length(q_norm))
  y.vector = as.numeric(model.y)
  try({for (i in 1:length(q_norm)) {
    coef_norm[,i] = rq.fit.pfn(norm.x, y=y.vector, tau=q_norm[i], max.bad.fixup=20)$coefficients 
    # coef_norm[,i] = rq.fit.fnb(norm.x, y=y.vector, tau=q_norm[i])$coefficients 
  }})
  yqhat = norm.x %*% coef_norm
  y_norm = (y.vector - yqhat[,q_norm==0.5]) / (yqhat[,3]  - yqhat[,1])
  rm(norm.x)
  gc()
  # BULK FIT
  bulk.x = getPredictors(n_files=n_files, df.x=bulk.x.df[1], df.t=bulk.x.df[2], df.xt=bulk.x.df[3], year.range=year.range)
  print('Fitting Bulk')
  p.bulk = bulk.x[2]
  coef_bulk = matrix(0, dim(bulk.x)[2], length(q))
  try({for (i in 1:length(q)) {
    if (q[i]==.5) {
      coef_bulk[,i] = 0  
    } else {
      coef_bulk[,i] = rq.fit.pfn(bulk.x, y=y_norm, tau=q[i], max.bad.fixup=20)$coefficients  
    }
  }})
  yhat_low = bulk.x %*% coef_bulk[,1]
  # yhat_low = bulk.x %*% coef_bulk[,q==.25]
  yhat_high = bulk.x %*% coef_bulk[,length(q)] 
  rm(bulk.x)
  gc()
  # For the tails use a simpler model
  tail.x = getPredictors(n_files=n_files, df.x=tail.x.df[1], df.t=tail.x.df[2], df.xt=tail.x.df[3], year.range=year.range)
  lowtail = y_norm - yhat_low
  hightail = y_norm - yhat_high

  # coef_low = rq(lowtail[lowtail<0] ~ tail.x[lowtail<0,] - 1, tau=q_low)$coef
  # coef_high = rq(hightail[hightail>0] ~ tail.x[hightail>0,] - 1, tau=q_high)$coef

  coef_low = matrix(0, dim(tail.x)[2], length(q_low))
  coef_high = matrix(0, dim(tail.x)[2], length(q_high))
  print('Fitting Tail')
  p.tail = tail.x[2]
  print(c(sum(lowtail<0), sum(hightail>0)))
  for (i in 1:length(q_tail)) {
    coef_low[,i] = rq.fit.pfn(tail.x[lowtail<0,], lowtail[lowtail<0], tau=q_low[i])$coef
    coef_high[,i] = rq.fit.pfn(tail.x[hightail>0,], hightail[hightail>0], tau=q_high[i])$coef
  }
  rm(tail.x)
  gc()
  list(coef_norm=coef_norm, coef_bulk=coef_bulk, coef_tail=list(coef_low, coef_high),
    q=q, q_tail=q_tail, q_norm=q_norm, q_all=q_all, norm.x.df=norm.x.df, 
    bulk.x.df=bulk.x.df, tail.x.df=tail.x.df, mlat=mlat, year.range=year.range)
}

get.terms <- function(norm.x) {
  nchars = unlist(lapply(colnames(norm.x), nchar))
  time_idx = grep("t.basis",colnames(norm.x), perl=TRUE)
  intercept_idx = grep("Intercept",colnames(norm.x))
  interaction_idx = grep(":",colnames(norm.x))
  time_idx = time_idx[ !(time_idx %in% interaction_idx)]
  if (any(nchars==0)) {
    time_idx = c(time_idx, which(nchars == 0))
  }
  seasonal_idx = grep("x.main*",colnames(norm.x))
  list(seasonal_idx=seasonal_idx,
    interaction_idx=interaction_idx,
    intercept_idx=intercept_idx,
    time_idx=time_idx)  
}



list_matrix_multiply <- function(m_list, m_matrix) {
  n = length(m_list)
  m = dim(m_matrix)[2]
  if (n != m) stop("Dimensions of list and matrix don't align.")
    j = 0
  do.call(cbind, lapply(m_list, function(mat_i) {
      j <<- j + 1
      mat_i %*% m_matrix[, j]
    }))
}

selective.rq.pfn <- function(x, y, q, which.pred) {
  nq = length(x)
  coefs = matrix(0, dim(x[[1]])[2], nq)
  print(dim(x[[1]]))
  try({for (i in 1:nq) {
    if (length(which.pred) != 0) {
        y.res = y - rowSums(x[[i]][, -which.pred])
        coefs[, which.pred] = rq.fit.pfn(x[[i]][, which.pred], y=y.res, 
          tau=q[i], max.bad.fixup=20)$coefficients 
        coefs[, -which.pred] = 1
    } else {
      stop("Need at least one predictor to regress")
    }
  }})
  coefs
}

selective.rq <- function(x, y, q, which.pred) {
  nq = length(x)
  np = dim(x[[1]])[2]
  coefs = matrix(0, np, nq)

  try({for (i in 1:nq) {
    if (length(which.pred) != 0) {
        if (np-length(which.pred) > 1) {
          y.res = y - rowSums(x[[i]][, -which.pred])
        } else {
          y.res = y - x[[i]][, -which.pred] 
        }
        coefs[which.pred, i] = rq.fit.br(x[[i]][, which.pred], y=y.res, 
          tau=q[i])$coefficients 
        coefs[-which.pred, i] = 1
    } else {
      stop("Need at least one predictor to regress")
    }
  }})
  coefs
}

reduce_predictors <- function(x, coef) {
  nq = dim(coef)[2]
  indeces = get.terms(x)
  reduced.x = lapply(indeces, function(idx) {
      my_idx = c(1:dim(coef)[1]) %in% idx
      if (sum(my_idx) == 1) {
        output = x[,my_idx] %*% t(coef[my_idx,]) 
        
      } else if (any(my_idx)) {
        output = x[,my_idx] %*% coef[my_idx,]
      } else {
        output = NULL  
      }
      output
    })
  null.idx = unlist(lapply(reduced.x, is.null))
  final.output = lapply(1:nq, function(mq) {
    do.call(cbind, lapply(reduced.x[!null.idx], function(mlist) mlist[,mq]))
  })
  final.output
}
# This function will perform the quantile mapping. rea.y has to be a 
# time series object, like xts. The other x and y can be arrays and matrices.
# Give q_tail as the fractions of exceedences beyond the upper and lower 
# regular percentile, given in the vector q. which.pred is indexed from 1 to 4
# as shown in `get.terms` simply remove one of the indeces to remove this
# predictor, e.g. remove index 2 to remove interaction.
quantile.map.nested <- function(input, model.y,  q, q_tail, q_norm, mlat=30,
    which.pred=c(1:4), which.pred.tail=c(2,3)) {

  nq = length(q)
  nq_norm = length(q_norm)
  n_files = dim(model.y)[2]
  q_low = q_tail
  q_high = rev(1-q_tail)
  q_all = c(q_low*q[1], q, q[length(q)] + q_high*(1-q[length(q)]))
  t = time(model.y)
  doy = as.numeric(strftime(t, format = "%j"))
  year = as.numeric(strftime(t, format = "%Y"))
  min.year = min(year)
  nyears = diff(range(year)) + 1
  year.range=range(year)

  # Create predictors 
  coef_norm = input$coef_norm
  coef_bulk = input$coef_bulk
  coef_tail = input$coef_tail
  norm.x.df = input$norm.x.df
  bulk.x.df = input$bulk.x.df
  tail_tail = input$tail_tail
  norm.x = getPredictors(n_files=n_files, df.x=norm.x.df[1], df.t=norm.x.df[2], df.xt=norm.x.df[3], year.range=year.range, get.volc=T, mlat=mlat)
  bulk.x = getPredictors(n_files=n_files, df.x=bulk.x.df[1], df.t=bulk.x.df[2], df.xt=bulk.x.df[3], year.range=year.range)
  tail.x = getPredictors(n_files=n_files, df.x=tail.x.df[1], df.t=tail.x.df[2], df.xt=tail.x.df[3], year.range=year.range)
  reduced_norm = reduce_predictors(norm.x, coef_norm)
  reduced_bulk = reduce_predictors(bulk.x, coef_bulk)
  reduced_tail_low = reduce_predictors(tail.x, coef_tail[[1]])
  reduced_tail_high = reduce_predictors(tail.x, coef_tail[[2]]) 
  
  #NORM FIT
  y.vector = as.numeric(model.y)
  
  reduced_coef_norm = selective.rq(reduced_norm, y.vector, q_norm, which.pred)
  
  yqhat = list_matrix_multiply(reduced_norm, reduced_coef_norm)
  yqhat_low = yqhat[,1]
  yqhat_mid = yqhat[,2]
  yqhat_high = yqhat[,3]
  y_norm = (y.vector - yqhat_mid) / (yqhat_high - yqhat_low)

  # BULK FIT
  reduced_coef_bulk = selective.rq(reduced_bulk, y_norm, q, which.pred)
  reduced_coef_bulk[,q==.5] = 0 # enforce zero here size this is
  # already fitted above in the normalization
  
  yhat_low = reduced_bulk[[1]] %*% reduced_coef_bulk[,1]
  # yhat_low = bulk.x %*% coef_bulk[,q==.25]
  yhat_high = reduced_bulk[[nq]] %*% reduced_coef_bulk[,nq] 
  # debug
  # idx = sample(seq(length(t)), 10000)
  # plot(doy[idx], y_norm[idx], col=alpha(1, .2))
  # lines(doy[idx], yhat_low[idx], col=alpha(2, .5))
  # TAIL FIT
  # For the tails use a simpler model
  lowtail = y_norm - yhat_low
  hightail = y_norm - yhat_high

  # coef_low = rq(lowtail[lowtail<0] ~ tail.x[lowtail<0,] - 1, tau=q_low)$coef
  # coef_high = rq(hightail[hightail>0] ~ tail.x[hightail>0,] - 1, tau=q_high)$coef
  subset_low = lapply(reduced_tail_low, function(x) x[lowtail<0,])
  subset_high = lapply(reduced_tail_high, function(x) x[hightail>0,])
  reduced_coef_low = selective.rq(subset_low, lowtail[lowtail<0], q_low, which.pred.tail)
  reduced_coef_high = selective.rq(subset_high, hightail[hightail>0], q_high, which.pred.tail)
  print("coef high")
  # print(reduced_coef_high)

  list(coef_norm=coef_norm, coef_bulk=coef_bulk, coef_tail=coef_tail,
    q=q, q_tail=q_tail, q_norm=q_norm, q_all=q_all, norm.x.df=norm.x.df, 
    bulk.x.df=bulk.x.df, tail.x.df=tail.x.df, mlat=mlat, year.range=year.range,
    reduced_coef_norm=reduced_coef_norm, reduced_coef_bulk=reduced_coef_bulk,
    reduced_coef_tail=list(reduced_coef_low, reduced_coef_high))
}


# cv - cross-validate anything!
# The idea here is that you give this function some data and 
# specify how you want the data to be trained and tested and out pops
# the loss metric
# Intput:
#   x - data in n x p matrix form
#   y - response in vector form
#   train.func - training function, taking x and y and spitting out some 
#       information to run the test function, or out of sample prediction
#        The output could e.g. be coefficients;
#       The coefficients can be a vector or a matrix of coefficients in each column
#   test.func - test function, taking in x and c, test data and the trained model
#   loss.func - loss function, taking in the test response and the estimate test response.
# Output:
#   loss - a list of length nfolds containing the loss
cv <-function(x, y, lambda, train.func=function(x, y, lambda) {lm(y~x)$coef}, 
  test.func=function(x, c) {x %*% c}, 
  loss.func=function(y, y.hat) {sqrt(mean(y - y.hat)^2)}, 
  temporal.data=TRUE, nfolds=5) {

  n = length(y)
  loss = list(); length(loss) = nfolds
  if (temporal.data) {
    ncols = ceiling(n/nfolds)
    foldid = c(t(matrix(rep(seq(nfolds), ncols), nfolds, ncols)))[1:n]
  } else {
    foldid = sample(seq(nfolds), size = n, replace=TRUE)
  }

  for (i in 1:nfolds) {
  # loss = mclapply(seq(nfolds), function(i, x, y, train.func, test.func, loss.func,
    # lambda) {
    x.tr = x[i != foldid, ]
    x.te = x[i == foldid, ]
    y.tr = y[i != foldid]
    y.te = y[i == foldid]

    # train model
    tr.fit = train.func(x=x.tr, y=y.tr, lambda=lambda)
    # test model
    y.te.hat = test.func(x.te, tr.fit)
    # plot(as.numeric(strftime(time(y.te), format = "%j")), y.te)
    # lines(as.numeric(strftime(time(y.te), format = "%j")), y.te.hat)
    # y.tr.hat = test.func(x.tr, tr.fit)
    # loss function
    loss[[i]] = loss.func(y.te, y.te.hat)
    # loss = loss.func(y.te, te.out)
  # }, x=x, y=y, train.func=train.func, test.func=test.func, 
  # loss.func=loss.func, lambda=lambda) 
  }
  
  loss
}

get.lambda.opt <- function(loss) {
  loss.matrix = do.call(rbind, loss)
  cvm = colMeans(loss.matrix)
  cvsd = apply(loss.matrix, 2, sd)
  lambda.min = max(lambda[cvm==min(cvm)])
  sd.idx = max(which(cvm < cvm[cvm==min(cvm)] + cvsd[cvm==min(cvm)]))
  lambda.1sd = lambda[sd.idx]
  cv.obj = list(lambda=lambda, cvm=cvm, cvsd=cvsd, 
    lambda.min=lambda.min, lambda.1sd=lambda.1sd)
  class(cv.obj) <- "cv.object"
  cv.obj
}

plot.cv.object <- function(obj) {  
  errbar(log(obj$lambda), obj$cvm, obj$cvm + obj$cvsd, obj$cvm - obj$cvsd,
    main="Test Set Loss", xlab=expression(log(lambda)), ylab="Loss")
  abline(v=log(obj$lambda.min))
  abline(v=log(obj$lambda.1sd), col=2)
}

# Pick the data record nearest the specified lat and lon
# points and get the reanalysis data
# Import reanalysis data
getReanalysisData2 <- function(mlat, mlon, npoints=1) {
  filename = '/project/moyer/CMIP5_Raw/ERA-INTERIM/tas_day_ERA-INTERIM_reanalysis_mean_19790101-20101231.nc'
  ncdata <- nc_open(filename) 
  lats = ncvar_get(ncdata, 'lat')
  lons = ncvar_get(ncdata, 'lon')
  time = ncvar_get(ncdata, 'time')
  last.t = max(which(time!=0))
  time = time[1:last.t]

  time.s=as.POSIXct('1979-01-01 12:00',tz='UTC')
  time.e=as.POSIXct('2010-12-31 12:00',tz='UTC')
  tseq=seq(time.s, time.e, by='24 hours')
  times=as.POSIXct(time*60*60, origin='1900-01-01 00:00:0.0', tz='UTC')
  t1=which(times==time.s)
  t2=which(times==time.e)
  
  coords = get_nearest_lonlat(lats, lons, mlat, mlon)

  dt=t2-t1+1
  output = apply(coords, 1, function(x) {
    rea_data = ncvar_get(ncdata, 'tas', start=c(x[1], x[2], t1),
    count=c(1, 1, dt)) - 273.15 
    rea_data
  }) 
  result = xts(output, order.by=tseq)
}

get_nearest_lonlat <- function(lats, lons, mlat, mlon, npoints=1) {
  dlat = abs(lats - mlat)
  dlon = abs(lons - mlon)
  lat_idx = order(dlat)[1:npoints]
  lon_idx = order(dlon)[1:npoints]
  coords = merge(lon_idx, lat_idx)
}

# Pick the data record nearest the specified lat and lon
# points and get the reanalysis data
# Import reanalysis data
getReanalysisData <- function(mlat, mlon, npoints=1) {
  filename = "../reanalysis/t_19790101_19940101.nc"
  ncdata <- nc_open(filename) 
  lats = ncvar_get(ncdata, 'latitude')
  lons = ncvar_get(ncdata, 'longitude')
  time = ncvar_get(ncdata, 'time')
  last.t = max(which(time!=0))
  time = time[1:last.t]

  time.s=as.POSIXct('1979-02-01 00:00',tz='UTC')
  time.e=as.POSIXct('1994-01-01 00:00',tz='UTC')
  tseq=seq(time.s, time.e, by='6 hours')
  times=as.POSIXct(time*60*60, origin='1900-01-01 00:00:0.0', tz='UTC')
  t1=which(times==time.s)
  t2=which(times==time.e)
  
  dlat = abs(lats - mlat)
  dlon = abs(lons - mlon)
  lat_idx = order(dlat)[1:npoints]
  lon_idx = order(dlon)[1:npoints]
  coords = merge(lon_idx, lat_idx)

  dt=t2-t1+1
  output = apply(coords, 1, function(x) {
    rea_data = ncvar_get(ncdata, 't2m', start=c(x[1], x[2], t1),
    count=c(1, 1, dt)) - 273.15 
    rea_data
  }) 
  xts(output, order.by=tseq)
}




# Pick the data record nearest the specified lat and lon
# points
getLENSData <- function(mlat, mlon, npoints=1) {
  KtoC = 273.15
  
  path.to.files = "/project/moyer/mahaugen/lens/timefirst/"
  files = list.files(path.to.files)
  
  # Find which pixel index to import
  ncdata <- nc_open('../lens/timefirst/LENS_001.nc')
  lons = signif(ncvar_get(ncdata, 'lon'), digits=3)
  lats = signif(ncvar_get(ncdata, 'lat'), 3)
  date = ncvar_get(ncdata, 'date')
  nc_close(ncdata)
  dlat = abs(lats - mlat)
  dlon = abs(lons - mlon)
  lat_idx = order(dlat)[1:npoints]
  lon_idx = order(dlon)[1:npoints]
  latlon = outer_concat(lat_idx, lon_idx)
  # Import data  
  one_pixel = lapply(files, function(f) {
        path_f = paste(path.to.files, f, sep="")
        point = rowMeans(apply(latlon, 1, function(lalo) {
                  lat_idx = lalo[1]
                  lon_idx = lalo[2]
                  get.ncdf.timefirst(path_f, lat_idx, 1, lon_idx, 1, 1, -1,
                  varname='TREFHT') - KtoC  
                }))    
      })
  # t = as.Date(as.character(date), format="%Y%m%d")
  t = timeBasedSeq('19200101::21001231')
  t = t[strftime(t, format = "%j")!='366'] # remove leap years

  ts_data = xts(do.call(cbind, one_pixel), order.by=t)
  rm(one_pixel)  
  ts_data 
}

# Pick the data record nearest the specified lat and lon
# points
getModelData2 <- function(mlat, mlon) {
  KtoC = 273.15
  npoints=1
  path.to.files = "../rsriver/timefirst2/"
  files = list.files(path.to.files)
  
  # Find which pixel index to import
  ncdata <- nc_open('../rsriver/timefirst2/trefht_4400.nc')
  lons = signif(ncvar_get(ncdata, 'lon'), digits=3)
  lats = signif(ncvar_get(ncdata, 'lat'), 3)
  date = ncvar_get(ncdata, 'date')
  nc_close(ncdata)
  dlat = abs(lats - mlat)
  dlon = abs(lons - mlon)
  lat_idx = order(dlat)[1:npoints]
  lon_idx = order(dlon)[1:npoints]

  # Import data
  one_pixel = lapply(files, function(f) {
        path_f = paste(path.to.files, f, sep="")
        get.ncdf.timefirst(path_f, lat_idx, 1, lon_idx, 1, 1, -1,
          varname='TREFHT') - KtoC
      })
  # t = as.Date(as.character(date), format="%Y%m%d")
  t = timeBasedSeq('18500101::21000101')
  t = t[strftime(t, format = "%j")!='366'] # remove leap years

  ts_data = xts(do.call(cbind, one_pixel), order.by=t)
  rm(one_pixel)  
  ts_data 
}


# Pick the data record nearest the specified lat and lon
# points
getModelData <- function(mlat, mlon) {
  KtoC = 273.15
  max_time = 365 * 250

  path.to.files = "../rsriver/timefirst/"
  files = list.files(path.to.files)
  
  # Find which pixel index to import
  ncdata <- nc_open('../rsriver/timefirst/trefht_4200.nc')
  lons = signif(ncvar_get(ncdata, 'lon'), digits=3)
  lats = signif(ncvar_get(ncdata, 'lat'), 3)
  nc_close(ncdata)
  dlat = abs(lats - mlat)
  dlon = abs(lons - mlon)
  lat_idx = which(dlat == min(dlat))[1]
  lon_idx = which(dlon == min(dlon))[1]

  # Import data
  one_pixel = lapply(files, function(f) {
        path_f = paste(path.to.files, f, sep="")
        get.ncdf.timefirst(path_f, lat_idx, 1, lon_idx, 1, 1, -1,
          varname='TREFHT') - KtoC
      })
  ts_data = unlist(lapply(one_pixel, function(ts) {
          ts[1:(max_time)]
          }))
  rm(one_pixel)  
  tseq= timeBasedSeq('18500101::20991231')
  tseq = tseq[strftime(tseq, format = "%j")!='366'] # remove leep years
  y.xts = xts(matrix(ts_data, 365*250, 50), order.by=tseq)
  y.xts 
}
