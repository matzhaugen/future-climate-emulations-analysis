# This code is intended to swap the time dimension with the space dimensions 
# in order to access the time data faster

library(ncdf4)
source('emulationHelper.R')

n_cores=4
# source("helper.R")
varname = 'TREFHT'
path.to.output = "../lens/timefirst/"
path.to.input= "../lens/spacefirst/"
files = list.files(path.to.input)

r = mclapply(files, function(f) {
	printf("Reading and writing file %s", f)
	path_f = paste(path.to.input, f, sep="")
	ncdata <- nc_open(path_f)
	lons = ncvar_get(ncdata, 'lon')
	lats = ncvar_get(ncdata, 'lat')
	time = ncvar_get(ncdata, 'time')
	date = ncvar_get(ncdata, 'date')
	lon_idx = which(lons > 180)
	lat_idx = which(lats > 0)
	year = floor(date / 1e4)
	# If handling the first file from 1850-2100 only use the 1920-2100 part
	subset = year %in% seq(1920,2100) & c(TRUE, diff(date)!=0)
	time_idx = which(subset)
	newdate = date[subset]
	
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
	datevar = ncdata$var$date
	datedim = ncdata$var$date$dim[[1]]
	datedim$vals = datedim$vals[subset]
	datedim$len = sum(subset)
	datedim$unlim = FALSE
	nc_close(ncdata)
	

	var3d = ncvar_def(varname, 'K', list(timedim, londim, latdim), 
		longname="Reference height temperature", chunksizes=c(timedim$len, 1, 1))
	datevar = ncvar_def('date', 'K', list(datedim), 
		longname="current date (YYYYMMDD)", prec='integer')

	path_f = paste(path.to.output, f, sep="")
	nc = nc_create(path_f, list(var3d, datevar))
	# nc = nc_create(path_f, datevar)
	# nc = nc_create(path_f, var3d)
	ncvar_put(nc, datevar, vals=newdate)
	for( i in 1:length(lon_idx)) {
		print(i)
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