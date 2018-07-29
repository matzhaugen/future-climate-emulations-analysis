source('emulationHelper.R')
load('raw_temperatures.RData')
mc.cores = 8 # Set depending on your CPU power
# If running on bigmem, set the following to 2, to use 16 cores in total.
mc.cores.sub = 3

####################################
# Pre defined constants #
####################################
q = c(.1, .18, 0.25, .35, .42, 0.5, .58, .65, 0.75, .82, 0.9)
q_tail = c(0.01, 0.1, .25, 0.5, .75)
q_low = q_tail
q_high = rev(1-q_tail)
q_all = c(q_low*q[1], q, q[length(q)] + q_high*(1-q[length(q)]))

nq = length(q)
q_norm = c(.1, .5, .9)
Sys.setenv(TZ='UTC')

location = list(c(30, 360-95), c(42, 360-90), c(48, 360-122), c(34, 360-118),
	c(40, 360-105), c(34, 360-112), c(41, 360-74), c(34, 360-84))
years_i = '1979/2016'
years_m = '2019/2056'
years_f = '2059/2096'
p.lens = 40
p.ryan = 50
K = 1 #number of jackknife runs coefficient (n.simulations / K)
nfolds.lens = p.lens / K
nfolds.ryan = p.ryan / K
foldid.ryan = ceil(c(1:p.ryan) / K)
foldid.lens = ceil(c(1:p.lens) / K)
norm.x.df = c(14, 6, 3)
bulk.x.df = norm.x.df
tail.x.df = c(3, 1, 0)
rea.x.df = c(10, 1, 0)

loc.names = c("Houston", "Chicago", "Seattle", "Los Angeles", "Denver", "Phoenix", "New York City", "Atlanta")

qout.all.inclusive = mclapply(1:length(location), function(j) {
	loc = location[[j]]
	mlat = loc[1]
	mlon = loc[2]
	#Import model data
	print('Fitting Ryans model')
	qout.ryan = quantile.map(norm.x.df, bulk.x.df, y.ryan[[j]], tail.x.df, 
				q, q_tail, q_norm, mlat=mlat)
		
	print('Fitting Lens model')
	qout.lens = quantile.map(norm.x.df, bulk.x.df, y.ryan[[j]], tail.x.df, 
				q, q_tail, q_norm, mlat=mlat)
	
	list(ryan=qout.ryan, lens=qout.lens)
}, mc.cores=mc.cores)

j = 0
which.pred = c(3,4)
which.pred.tail = c(2,3)
qout.nested = mclapply(1:length(qout.all.inclusive), function(j) {
	print(j)
	mlat = location[[j]][1]
	mlon = location[[j]][2]
	#Import model data
	print('Fitting Ryan model')
	qout.lens.nested = quantile.map.nested(qout.all.inclusive[[j]]$lens, 
		y.ryan[[j]][,1],  q, q_tail, q_norm, 
		mlat=mlat, which.pred=which.pred, which.pred.tail=which.pred.tail)
		
	print('Fitting Lens model')
	qout.ryan.nested = quantile.map.nested(qout.all.inclusive[[j]]$ryan, 
		y.lens[[j]][,1],  q, q_tail, q_norm, 
		mlat=30, which.pred=which.pred, which.pred.tail=which.pred.tail)

	list(ryan.then.lens=qout.ryan.nested, lens.then.ryan=qout.lens.nested)
}, mc.cores=mc.cores)

qout.jack = mclapply(1:length(location), function(j) {
	loc = location[[j]]
	mlat = loc[1]
	mlon = loc[2]
	#Import model data
	print('Fitting Ryans model')
	# qout.ryan = mclapply(1:2, function(k) {
	qout.ryan = mclapply(1:nfolds.ryan, function(k) {
			print(paste('run #', k))
			quantile.map(norm.x.df=norm.x.df, bulk.x.df=bulk.x.df, 
				model.y=y.ryan[[j]][,foldid.ryan!=k], tail.x.df=tail.x.df, 
				q=q, q_tail=q_tail, q_norm=q_norm, mlat=mlat)
		}, mc.cores=mc.cores.sub)		
	print('Fitting Lens model')
	# qout.lens = mclapply(1:2, function(k) {
	qout.lens = mclapply(1:nfolds.lens, function(k) {
			print(paste('run #', k))
			quantile.map(norm.x.df, bulk.x.df, 
				y.lens[[j]][,foldid.lens!=k], tail.x.df, 
				q, q_tail, q_norm, mlat=mlat)
		}, mc.cores=mc.cores.sub)

	list(ryan=qout.ryan, lens=qout.lens)
}, mc.cores=mc.cores)

save(qout.nested, 
	 qout.jack,
	 qout.all.inclusive, 
 	 file=paste('everything_',
 		implode(norm.x.df), '_',
 		implode(bulk.x.df), '_',
 		implode(tail.x.df), '_',
 		implode(rea.x.df), '.Rdata', sep=''))
