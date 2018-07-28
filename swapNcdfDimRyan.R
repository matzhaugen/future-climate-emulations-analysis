# This code is intended to swap the time dimension with the space dimensions 
# in order to access the time data faster

library(ncdf4)
source('emulationHelper.R')

n_cores=4
# source("helper.R")
varname = 'TREFHT'
path.to.output = "../rsriver/timefirst2/"
path.to.input= "../rsriver/spacefirst/"
files = list.files(path.to.input)
startyear = 1850
r = mclapply(files, function(f) {
	printf("Reading and writing file %s", f)
	path_f = paste(path.to.input, f, sep="")
	ncdata <- nc_open(path_f)
	lons = ncvar_get(ncdata, 'lon')
	lats = ncvar_get(ncdata, 'lat')
	time = ncvar_get(ncdata, 'time')
	lon_idx = which(lons > 180)
	lat_idx = which(lats > 0)
	year = floor(time / 365) + 1850
	# If handling the first file from 1850-2100 only use the 1920-2100 part
	subset = year %in% seq(startyear, 2100) & c(TRUE, diff(time)!=0)
	time_idx = which(subset)
	date = as.Date(paste(as.character(year), as.character(time - (year-1850)*365 + 1), sep='-'), 
		format='%Y-%j')
	newdate = as.integer(format(date[subset], '%Y%m%d'))
	
	data = ncvar_get(ncdata, varname, 
		start=c(lon_idx[1], lat_idx[1], time_idx[1]), 
		count=c(length(lon_idx), length(lat_idx), length(time_idx)))
	latdim = ncdata$dim[[which(names(ncdata$dim) == "lat")]]
	latdim$vals = lats[lat_idx]
	latdim$len = length(lats[lat_idx])
	londim = ncdata$dim[[which(names(ncdata$dim) == "lon")]]
	londim$vals = lons[lon_idx]
	londim$len = length(lons[lon_idx])
	timedim = ncdata$dim[[which(names(ncdata$dim) == "time")]]
	timedim$vals = timedim$vals[subset]
	timedim$len = sum(subset)
	timedim$unlim = FALSE
	datedim = ncdim_def(name='Date', units='YYYYMMDD', calendar='noleap', vals=1:sum(subset), unlim=FALSE)
	nc_close(ncdata)
	

	var3d = ncvar_def(varname, 'K', list(timedim, londim, latdim), 
		longname="Reference height temperature", chunksizes=c(timedim$len, 1, 1))
	datevar = ncvar_def('date', '', list(datedim), 
		longname="current date (YYYYMMDD)", prec='integer')

	path_f = paste(path.to.output, f, sep="")
	nc = nc_create(path_f, list(var3d, datevar))
	# nc = nc_create(path_f, datevar)
	# nc = nc_create(path_f, var3d)
	ncvar_put(nc, datevar, vals=newdate)
	for( i in 1:length(lon_idx)) {
		for (j in 1:length(lat_idx)){
				temp = ncvar_put(nc, var3d, data[i, j,], start=c(1,i,j), count=c(-1,1,1) )
			}
	}
	
	nc_close(nc)
}, mc.cores=n_cores)

# Testing
# ncdata = nc_open(paste(path.to.output, f, sep=""))
# print(ncdata)
# nc_close(ncdata)