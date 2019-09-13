# 1)

# loading csv file
ParticulateMatterData <- read.csv(file="D:/UPitt/Studies/Sem 2/Spatial DA/Projects/Project-3/ParticulateMatter.csv", 
                                  header=TRUE, sep=",")

# computing the mean of PM25 attribute
valMean <- mean(ParticulateMatterData$PM25)
# creating a new column in Data Frame to store (x-mean(x))
ParticulateMatterData['D'] = ParticulateMatterData$PM25 - valMean
# creating a new column in Data Frame to store square of (x-mean(x))
ParticulateMatterData['D2'] = ParticulateMatterData$D^2
# get the distance matrix by applying dist function.
distance <- as.matrix(dist(cbind(ParticulateMatterData$Lon, ParticulateMatterData$Lat)))
# get inverse distance
distance.inv <- 1/distance
# make diagonals=0 in case they are not equal to 0
diag(distance.inv) <- 0

library(reshape2)
# convert matrix to data frame
distancedf <- melt(distance.inv, varnames = c("row", "col"))
# N = number of records in the data set
N <- nrow(ParticulateMatterData)
# W = sum of all the distances in the distance matrix
W <- sum(distancedf$value)
# denomsum = sum of the newly added D2 column
denomsum <- sum(ParticulateMatterData$D2)
# initializing numsum to 0
numsum <- 0
  
# loop over all the rows in the distance data frame(1764 rows)
for (row in 1:nrow(distancedf)) {
  # rowval = 'row' column value of the specific row
  rowval <- distancedf[row, "row"]
  # columnval = 'col' column value of the specific row
  columnval  <- distancedf[row, "col"]
  # weight = 'value' column value of the specific row
  weight <- distancedf[row, "value"]
  # get the (x-mean(x)) value of the specific row
  prdt1 <- ParticulateMatterData[rowval, "D"]
  # get the (x-mean(x)) value of the specific column
  prdt2 <- ParticulateMatterData[columnval, "D"]
  # numsum = numerator value 
  numsum <- numsum + (weight*prdt1*prdt2) 
}

moranvalue <- (N*numsum)/(W*denomsum)
cat("Moran's I value is: ", moranvalue)

# 2)
library(rgdal)
library(GISTools)  
library(spatstat)

# loading 1st .shp file
OilGasLocnPA <- readOGR(dsn="D:/UPitt/Studies/Sem 2/Spatial DA/Projects/Project-3/OilGasLocationPA/OilGasLocationPA",
                  layer = "OilGasLocationPA") 

# creating a new data frame with only the coordinates
oildf<-data.frame(coordinates(OilGasLocnPA))
# x-minimum coordinate
xleftoil <- apply(OilGasLocnPA@coords,2,min)[1]
# y-minimum coordinate
ybottomoil <- apply(OilGasLocnPA@coords,2,min)[2]
# x-maximum coordinate
xrightoil <- apply(OilGasLocnPA@coords,2,max)[1]
# y-maximum coordinate
ytopoil <- apply(OilGasLocnPA@coords,2,max)[2]

# loading 2nd .shp file
IndMinMiningPA <- readOGR(dsn="D:/UPitt/Studies/Sem 2/Spatial DA/Projects/Project-3/IndustrialMineralMiningPA",
                        layer = "IndustrialMineralMiningOperations2014_10") 

# creating a new data frame with only the coordinates
miningdf<-data.frame(coordinates(IndMinMiningPA))
# x-minimum coordinate
xleftmining <- apply(IndMinMiningPA@coords,2,min)[1]
# y-minimum coordinate
ybottommining <- apply(IndMinMiningPA@coords,2,min)[2]
# x-maximum coordinate
xrightmining <- apply(IndMinMiningPA@coords,2,max)[1]
# y-maximum coordinate
ytopmining <- apply(IndMinMiningPA@coords,2,max)[2]

# to enable things to be drawn outside the plot region
par(xpd=TRUE)

# plot all the events in OilGasLocnPA
plot(oildf, pch=20, col="orange", main="PAOilAndGasLocations",
     xlab="x-coordinate", ylab="y-coordinate")

#add Legend
legend(xrightoil-95000, ytopoil+100000, legend=c("OilLocations"),
       col=c("orange"), pch=c(20, 0), cex=0.8,
       title="Legend", text.font=4)

# add North arrow
north.arrow(xleftoil+10000, ytopoil+60000, len=10000, lab="N", col="red", tcol="red")

# plot all the events in IndMinMiningPA
plot(miningdf, pch=20, col="green", main="PAMineralMiningLocations",
     xlab="x-coordinate", ylab="y-coordinate")

legend(xrightmining-95000, ytopmining+100000, legend=c("MiningLocations"),
              col=c("green"), pch=c(20, 0), cex=0.8,
              title="Legend", text.font=4)


# add North arrow
north.arrow(xleftmining+10000, ytopmining+60000, len=10000, lab="N", col="red", tcol="red")

PAOilppp <- as.ppp(OilGasLocnPA@coords, c(xleftoil,xrightoil,
                                          ybottomoil,ytopoil))
qCA <-  quadratcount(PAOilppp,100,100)

PAMiningppp <- as.ppp(IndMinMiningPA@coords, c(xleftmining,xrightmining,
                                        ybottommining,ytopmining))
qCB <-  quadratcount(PAMiningppp,100,100)


# To compute G-Function of PAOil
GPAOil <- Gest(PAOilppp)

# To compute G-Function of PAMining
GPAMining <- Gest(PAMiningppp)

# To compute F-Function of PAOil
FPAOil <- Fest(PAOilppp)

# To compute F-Function of PAMining
FPAMining <- Fest(PAMiningppp)

# To compute K-Function of PAOil
KPAOil <- Kest(PAOilppp)

# To compute K-Function of PAMining
KPAMining <- Kest(PAMiningppp)

# To compute L-Function of PAOil
LPAOil <- Lest(PAOilppp)

# To compute L-Function of PAMining
LPAMining <- Lest(PAMiningppp)

# To plot the result of G-function of PAOil
plot(GPAOil, main="G-Function of PAOilAndGas")

# To plot the result of G-function of PAMining
plot(GPAMining, main="G-Function of PAMineralMining")

# To plot the result of F-function of PAOil
plot(FPAOil, main="F-Function of PAOilAndGas")

# To plot the result of F-function of PAMining
plot(FPAMining, main="F-Function of PAMineralMining")

# To plot the result of K-function of PAOil
plot(KPAOil, main="K-Function of PAOilAndGas")

# To plot the result of K-function of PAMining
plot(KPAMining, main="K-Function of PAMineralMining")

# To plot the result of L-function of PAOil
plot(LPAOil, main="L-Function of PAOilAndGas")

# To plot the result of L-function of PAMining
plot(LPAMining, main="L-Function of PAMineralMining")


