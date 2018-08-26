library(diagram)
library(seraphim)
library(lubridate)
library(OutbreakTools)

e_YFV = extent(-53,-37,-26.5,-12)
models = c("cauchy")

# 1. Extraction of 1,000 trees

burnIn = 0
randomSampling = FALSE
nberOfTreesToSample = 1000
mostRecentSamplingDatum = decimal_date(ymd("2017-04-22"))
coordinateAttributeName = "Locations"
nberOfCores = 1
for (m in 1:length(models))
	{
		localTreesDirectory = paste0("Tree_extractions/YFV1_BR_",models[m])
		allTrees = scan(file=paste0("YFV_skygrid_",models[m],".trees"), what="", sep="\n", quiet=TRUE)
		treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName, nberOfCores)
	}

# 2. Generating graphical spread and plotting the MCC tree

for (m in 1:length(models))
	{
		file1 = scan(paste0("YFV_skygrid_",models[m],".trees"), what="", sep="\n", quiet=TRUE)
		index = which(grepl("tree STATE_",file1))[1]-1
		file2 = c(file1[1:index], file1[seq(index+16,index+16000,16)],"End;")
		write(file2, paste0("YFV_skygrid_",models[m],"_1000.trees"))
	}
# to do : building the MCC tree from these 1,000 selected trees (= those extracted above)

for (m in 1:length(models))
	{
		mostRecentSamplingDatum = decimal_date(ymd("2017-04-22"))
		mcc_tre = read.annotated.nexus(paste0("YFV_skygrid_",models[m],"_MCC.tree"))
		mcc_tab = matrix(nrow=dim(mcc_tre$edge)[1], ncol=11)
		colnames(mcc_tab) = c("node1","node2","length","startLon","startLat","endLon","endLat","startNodeL","endNodeL","startYear","endYear")
		mcc_tab[,c("node1","node2")] = mcc_tre$edge
		mcc_tab[,c("length")] = mcc_tre$edge.length
		for (i in 1:length(mcc_tre$annotations))
			{
				annotations = mcc_tre$annotations[[i]]
				mcc_tab[i,c("endLon","endLat")] = cbind(annotations$Locations2, annotations$Locations1)
				mcc_tab[i,"endNodeL"] = annotations$height
			}
		for (i in 1:length(mcc_tre$annotations))
			{
				index = which(mcc_tab[,"node2"] == mcc_tab[i,"node1"])
				if (length(index) > 0)
					{
						mcc_tab[i,c("startLon","startLat")] = mcc_tab[index,c("endLon","endLat")]
						mcc_tab[i,"startNodeL"] = mcc_tab[index,"endNodeL"]
					}	else		{
						annotations = mcc_tre$root.annotation
						mcc_tab[i,c("startLon","startLat")] = cbind(annotations$Locations2, annotations$Locations1)
						mcc_tab[i,"startNodeL"] = annotations$height
					}
				mcc_tab[i,"startYear"] = mostRecentSamplingDatum - mcc_tab[i,"startNodeL"]
				mcc_tab[i,"endYear"] = mostRecentSamplingDatum - mcc_tab[i,"endNodeL"]
			}
		write.csv(mcc_tab, paste0("YFV_skygrid_",models[m],"_MCC.csv"), row.names=F, quote=F)
	}

nberOfExtractionFiles = 1000
elevation = raster("Environmental_files/YFV_rasters/Elevation_YF_0.008.asc")
emptyRaster = elevation; emptyRaster[is.na(emptyRaster[])] = 0; rast = emptyRaster
prob = 0.95
precision = 0.025
showingPlots = FALSE
nberOfCores = 1
origin = FALSE
for (m in 1:length(models))
	{
		localTreesDirectory = paste0("Tree_extractions/YFV1_BR_",models[m])
		mcc = read.csv(paste0("YFV_skygrid_",models[m],"_MCC.csv"), header=T)
		startDatum = min(mcc[,"startYear"])
		timeLayers = FALSE
		spread = spreadGraphic(localTreesDirectory, nberOfExtractionFiles, rast, prob, startDatum, precision, nberOfCores, origin, timeLayers)
		spread[is.na(elevation[])] = NA
		writeRaster(spread, paste0("YFV_BR_",models[m],"_spread.asc"), overwrite=T)
		timeLayers = TRUE
		spreads = spreadGraphic(localTreesDirectory, nberOfExtractionFiles, rast, prob, startDatum, precision, nberOfCores, origin, timeLayers)
		for (i in 1:length(spreads))
			{
				s = spreads[[i]]; s[is.na(elevation[])] = NA
				writeRaster(s, paste0("YFV_BR_",models[m],"_spreads/YFV_BR_",models[m],"_spread_",i,".asc"), overwrite=T)
			}
	}

e_smaller = extent(-50,-37,-26,-13)
# cols = colorRampPalette(brewer.pal(11,'RdYlGn'))(131)[31:131]
cols = colorRampPalette(brewer.pal(9,'YlGn'))(141)[21:121]
background = raster("Environmental_files/Natural_Earth/Gray_background.tif")
elevation = raster("Environmental_files/YFV_rasters/Elevation_YF_0.008.asc")
emptyRaster = elevation; emptyRaster[!is.na(emptyRaster[])] = 0
emptyRaster = crop(emptyRaster, e_smaller)
borders = getData("GADM",country="BRA",level=1)
borders = crop(borders, e_smaller)
background = crop(background, e_YFV); borders = crop(borders, e_YFV)
background_cols = colorRampPalette(c("grey","white"),bias=1)(max(background[],na.rm=T)-min(background[],na.rm=T))[1:(max(background[],na.rm=T)-min(background[],na.rm=T))]
background = crop(background, e_smaller)
colNA = "grey90"
for (m in 1:length(models))
	{
		mcc = read.csv(paste0("YFV_skygrid_",models[m],"_MCC.csv"), header=T)
		mcc = mcc[order(mcc[,"startYear"],decreasing=F),]
		mcc1 = mcc[1,]; mcc2 = mcc[2:dim(mcc)[1],]
		mcc2 = mcc2[order(mcc2[,"endYear"],decreasing=F),]
		mcc = rbind(mcc1, mcc2)
		spread = raster(paste0("YFV_BR_",models[m],"_spread.asc"))
		spread = crop(spread, e_smaller)
		minYear = min(c(spread[!is.na(spread[])], mcc[,"startYear"]))
		maxYear = max(c(spread[!is.na(spread[])], mcc[,"endYear"]))
		spread[1] = maxYear; spread[2] = minYear
		spreads = list(); maxYears = list(); nberOfFiles = length(list.files(paste0("YFV_BR_",models[m],"_spreads")))
		for (i in 1:nberOfFiles)
			{
				s = raster(paste0("YFV_BR_",models[m],"_spreads/YFV_BR_",models[m],"_spread_",i,".asc"))
				s = crop(s, e_smaller); maxYears[[i]] = max(s[], na.rm=T); s[1] = maxYear; s[2] = minYear; spreads[[i]] = s
			}
		endYearsC = (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
		YFV_logN_cols_1 = cols[endYearsC]
		pdf(paste0("YFV_BR_",models[m],"_spread.pdf"), width=6, height=6); # dev.new(width=6, height=6)
		par(mar=c(0,0,0,0), oma=c(1.2,3.5,1,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
		plot(emptyRaster, col="white", box=F, axes=F, colNA="grey90", legend=F)
		plot(spread, gridded=T, useRaster=T, maxpixels=length(spread[]), col=cols, add=T, legend=F)
		plot(borders, add=T, lwd=0.1, border="gray10")
		for (i in 1:dim(mcc)[1])
			{
				curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
						    arr.width=0, lwd=0.2, lty=1, lcol="gray10", arr.col=NA, arr.pos=FALSE, curve=0.1, dr=NA, endhead=F)
			}
		for (i in dim(mcc)[1]:1)
			{	
				if (i == 1)
					{
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=cols[1], cex=0.8)
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray10", cex=0.8)
					}
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=YFV_logN_cols_1[i], cex=0.8)
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray10", cex=0.8)
			}
		rect(xmin(emptyRaster), ymin(emptyRaster), xmax(emptyRaster), ymax(emptyRaster), xpd=T, lwd=0.2)
		axis(1, c(ceiling(xmin(emptyRaster)), floor(xmax(emptyRaster))), pos=ymin(emptyRaster), mgp=c(0,0.3,0), cex.axis=0.5, lwd=0, lwd.tick=0.2, padj=-0.8, tck=-0.01, col.axis="gray30")
		axis(2, c(ceiling(ymin(emptyRaster)), floor(ymax(emptyRaster))), pos=xmin(emptyRaster), mgp=c(0,0.6,0), cex.axis=0.5, lwd=0, lwd.tick=0.2, padj=1, tck=-0.01, col.axis="gray30")
		plot(spread, legend.only=T, add=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.40,0.81,0.13,0.15),
			 legend.args=list(text="years", cex=0.7, line=0.3, col="gray30"), horizontal=T,
		     axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, tck=-0.4, col.axis="gray30", line=0, mgp=c(0,0.05,0)))
		dev.off()
		for (s in 1:length(spreads))
			{
				spread_maxYear = maxYears[[s]]
				if (s == length(spreads)) spread_maxYear = max(mcc[,"endYear"])
				pdf(file=paste0("YFV_BR_",models[m],"_spreads/YFV_BR_",models[m],"_spread_",round(spread_maxYear,2),".pdf"), width=6, height=6); # dev.new(width=6, height=6)
				par(mar=c(0,0,0,0), oma=c(1.2,3.5,1,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
				plot(emptyRaster, col="white", box=F, axes=F, colNA="grey90", legend=F)
				plot(spreads[[s]], gridded=T, useRaster=T, maxpixels=length(spread[]), col=cols, add=T, legend=F)
				plot(borders, add=T, lwd=0.1, border="gray10")
				for (i in 1:dim(mcc)[1])
					{
						if (mcc[i,"endYear"] < maxYears[[s]]) # print(i)
							{
								curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
						   		    			arr.width=0, lwd=0.2, lty=1, lcol="gray10", arr.col=NA, arr.pos=FALSE, curve=0.1, dr=NA, endhead=F)
						   	}	else		{
								if (s == length(spreads))
									{
										curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
						   		    			arr.width=0, lwd=0.2, lty=1, lcol="gray10", arr.col=NA, arr.pos=FALSE, curve=0.1, dr=NA, endhead=F)
									}
							}
					}
				for (i in dim(mcc)[1]:1)
					{	
						if (i == 1)
							{
								points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=cols[1], cex=0.8)
								points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray10", cex=0.8)
							}
						if (mcc[i,"endYear"] < maxYears[[s]]) # print(i)
							{
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=YFV_logN_cols_1[i], cex=0.8)
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray10", cex=0.8)
							}	else		{
								if (s == length(spreads))
									{
										points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=YFV_logN_cols_1[i], cex=0.8)
										points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray10", cex=0.8)
									}
							}	
					}
				# axis(1, c(ceiling(xmin(emptyRaster)), floor(xmax(emptyRaster))), pos=ymin(emptyRaster), mgp=c(0,0.3,0), cex.axis=0.5, lwd=0, lwd.tick=0.2, padj=-0.8, tck=-0.01, col.axis="gray30")
				# axis(2, c(ceiling(ymin(emptyRaster)), floor(ymax(emptyRaster))), pos=xmin(emptyRaster), mgp=c(0,0.6,0), cex.axis=0.5, lwd=0, lwd.tick=0.2, padj=1, tck=-0.01, col.axis="gray30")
				rect(xmin(emptyRaster), ymin(emptyRaster), xmax(emptyRaster), ymax(emptyRaster), xpd=T, lwd=0.2)
				plot(spread, legend.only=T, add=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.40,0.81,0.13,0.15),
			 		legend.args=list(text="years", cex=0.7, line=0.3, col="gray30"), horizontal=T,
		     		axis.args=list(at=spread_maxYear, labels=c(round(spread_maxYear,2)), cex.axis=0.6, lwd=0, lwd.tick=0.2, tck=-0.4, col.axis="gray30", line=0, mgp=c(0,0.05,0)))
				dev.off()
			}
	}

# 3. Dispersion statistics estimations

for (m in 1:length(models))
	{
		localTreesDirectory = paste0("Tree_extractions/YFV1_BR_",models[m])
		nberOfExtractionFiles = 1000
		timeSlices = 100
		timeInterval = 1
		onlyTipBranches = FALSE
		showingPlots = FALSE
		outputName = paste0("YFV_",models[m])
		nberOfCores = 1
		spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores)
	}

for (m in 1:length(models))
	{
		pdf(paste0("Wavefront_",models[m],".pdf"), width=5, height=3.5); # dev.new(width=5, height=3.5)
		LWD = 0.2; LINE = 0.7; par(mgp=c(10,0.35,0), oma=c(1,0.5,1,3), mar=c(3,4,2,0))
		waveFrontDistances1MedianValue = read.table(paste0("YFV_",models[m],"_median_spatial_wavefront_distance.txt"), header=T)[,"distance"]
		lower_l_1 = read.table(paste0("YFV_",models[m],"_95%HPD_spatial_wavefront_distance.txt"), header=T)[,2]
		upper_l_1 = read.table(paste0("YFV_",models[m],"_95%HPD_spatial_wavefront_distance.txt"), header=T)[,3]
		xLim = c(2016, 2017.304); yLim = c(0, 900)
		text = "Furthest extent of epidemic wavefront (spatial distance from epidemic origin)"
		xLab = "years"; yLab="wavefront distance (km)"
		slicedTimes = read.table(paste0("YFV_",models[m],"_median_spatial_wavefront_distance.txt"), header=T)[,"time"]
		plot(slicedTimes, waveFrontDistances1MedianValue, type="l", axes=F, ann=F, ylim=yLim, xlim=xLim)
		xx_l = c(slicedTimes,rev(slicedTimes)); yy_l = c(lower_l_1,rev(upper_l_1))
		getOption("scipen"); opt = options("scipen"=20)
		polygon(xx_l, yy_l, col=rgb(187/255,187/255,187/255,0.5), border=0)
		lines(slicedTimes, waveFrontDistances1MedianValue, lwd=1)
		axis(side=1, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.020, col.axis="gray30", mgp=c(0,0.0,0))
		axis(side=2, lwd.tick=LWD, cex.axis=0.6, lwd=0, tck=-0.015, col.axis="gray30", mgp=c(0,0.2,0))
		title(xlab=xLab, cex.lab=0.7, mgp=c(1.0,0,0), col.lab="gray30")
		title(ylab=yLab, cex.lab=0.7, mgp=c(1.2,0,0), col.lab="gray30")
		title(main=text, cex.main=0.55, col.main="gray30", line=LINE)
		box(lwd=LWD, col="gray30")
		dev.off()
	}

tab = read.table("YFV_cauchy_median_spatial_wavefront_distance.txt", header=T)
(tab[89,2]-tab[62,2])/((tab[89,1]-tab[62,1])*365) # wavefront velocity between JUl. 2016 and Feb. 2017

stats = read.table("YFV_cauchy_estimated_statistics.txt", header=T); stat = stats[,1]/365
print(paste0(round(median(stat),2)," km/day [",round(quantile(stat,0.025),2),", ",round(quantile(stat,0.975),2),"]"))
