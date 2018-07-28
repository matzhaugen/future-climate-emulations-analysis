source('emulationHelper.R')
load('raw_temperatures.RData')
load('everything_1463_1463_310_1010.Rdata')
mc.cores = 8 
norm.x.df = c(14, 6, 3)
bulk.x.df = c(14, 6, 3)
tail.x.df = c(3, 1, 0)
rea.x.df = c(10, 1, 0)
location = list(c(30, 360-95), c(42, 360-90), c(48, 360-122), c(34, 360-118),
	c(40, 360-105), c(34, 360-112), c(41, 360-74), c(34, 360-84))
years_i = '1979/2016'
years_f = '2059/2096'
d.year = 80
const.projection.year = 2099

# Basic Maps using all simulations
maps.all.inclusive = mclapply(1:length(qout.all.inclusive), function(j) {
	mlat = location[[j]][1]
	mlon = location[[j]][2]

	print('Mapping Ryans model')
	map.ryan = map.data(qout.all.inclusive[[j]]$ryan, 
		rea.x.df, y.ryan[[j]][years_i, ], d.year=d.year)

	print('Mapping Lens model')	
	map.lens = map.data(qout.all.inclusive[[j]]$lens, 
		rea.x.df, y.lens[[j]][years_i, ], d.year=d.year)
	
	print('Mapping RyanOnLens model')
	map.ryan.lens = map.data(qout.all.inclusive[[j]]$ryan, 
			rea.x.df, y.lens[[j]][years_i, ], d.year=d.year)
	
	print('Mapping LensOnRyan model')
	map.lens.ryan = map.data(qout.all.inclusive[[j]]$lens, 
			rea.x.df, y.ryan[[j]][years_i, ], d.year=d.year)
	
	list(ryan=map.ryan, lens=map.lens, 
		ryan.lens=map.ryan.lens, lens.ryan=map.lens.ryan)
}, mc.cores=mc.cores)

# Project data using the identity transform
maps.lens.identity = mclapply(1:length(qout.all.inclusive), function(j) {
	print(j)
	mlat = location[[j]][1]
	mlon = location[[j]][2]
	map.data(qout.all.inclusive[[j]]$lens, 
			rea.x.df, y.lens[[j]][years_i, ], 
			d.year=d.year, 
			identity.transform=TRUE)
}, mc.cores=mc.cores)

# Project every year into 2099 - Reanalysis data
maps.all.inclusive.const.year.rea = mclapply(1:length(qout.all.inclusive), function(j) {
	print(j)
	mlat = location[[j]][1]
	mlon = location[[j]][2]

	map.ryan.rea = map.data(qout.all.inclusive[[j]]$ryan, 
		rea.x.df, rea.y[[j]][years_i], const.projection.year=const.projection.year)
	map.lens.rea = map.data(qout.all.inclusive[[j]]$lens, 
		rea.x.df, rea.y[[j]][years_i], const.projection.year=const.projection.year)
	list(ryan=map.ryan.rea , lens=map.lens.rea)
}, mc.cores=mc.cores)
# Project every year into 2099 - Model data
maps.all.inclusive.const.year = mclapply(1:length(qout.all.inclusive), function(j) {
	mlat = location[[j]][1]
	mlon = location[[j]][2]
	
	print('Mapping Ryans model')	
	map.ryan = map.data(params=qout.all.inclusive[[j]]$ryan, 
		rea.x.df=rea.x.df, rea.y=y.ryan[[j]][years_i, ], const.projection.year=const.projection.year)

	print('Mapping Lens model')	
	y.te = y.lens
	map.lens = map.data(qout.all.inclusive[[j]]$lens, 
		rea.x.df, y.lens[[j]][years_i, ], const.projection.year=const.projection.year)
	
	print('Mapping RyanOnLens model')
	y.te = y.lens
	map.ryan.lens = map.data(qout.all.inclusive[[j]]$ryan, 
			rea.x.df, y.lens[[j]][years_i, ], const.projection.year=const.projection.year)
	
	print('Mapping LensOnRyan model')
	y.te = y.ryan
	map.lens.ryan = map.data(qout.all.inclusive[[j]]$lens, 
			rea.x.df, y.ryan[[j]][years_i, ], const.projection.year=const.projection.year)
	
	list(ryan=map.ryan, lens=map.lens, 
		ryan.lens=map.ryan.lens, lens.ryan=map.lens.ryan)
}, mc.cores=mc.cores)

maps.nested = mclapply(1:length(qout.nested), function(j) {
	mlat = location[[j]][1]
	mlon = location[[j]][2]

	print('Mapping Ryans model')	
	map.ryan.then.lens = map.data(params=qout.nested[[j]]$ryan.then.lens, 
		rea.x.df=rea.x.df, rea.y=y.ryan[[j]][years_i, ], 
		const.projection.year=const.projection.year, nested=TRUE)

	print('Mapping Lens model')	
	map.lens.then.ryan = map.data(qout.nested[[j]]$lens.then.ryan, 
		rea.x.df, y.lens[[j]][years_i, ], const.projection.year=const.projection.year, nested=TRUE)
	
	print('Mapping Ryans-then-lens model with reanalysis data')	
	map.ryan.then.lens.rea = map.data(qout.nested[[j]]$ryan.then.lens, 
		rea.x.df, rea.y[[j]], const.projection.year=const.projection.year, nested=TRUE)

	list(ryan.then.lens=map.ryan.then.lens,
		ryan.then.lens.rea=map.ryan.then.lens.rea,
		lens.then.ryan=map.lens.then.ryan)
}, mc.cores=mc.cores)

save(maps.all.inclusive, 
	maps.lens.identity, 
	maps.nested, 
	maps.all.inclusive.const.year.rea,
	maps.all.inclusive.const.year,
 	file=paste('maps_',
 		implode(norm.x.df), '_',
 		implode(bulk.x.df), '_',
 		implode(tail.x.df), '_',
 		implode(rea.x.df), '.Rdata', sep=''))
