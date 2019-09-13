
require(rgdal)

# Read .shp file from below path
ogShape <- readOGR(dsn = "D:/UPitt/Studies/Sem 2/Spatial DA/Projects/Project 2/OilGasLocationPA/OilGasLocationPA", layer = "OilGasLocationPA")
# Convert SpatialPointsDataFrame object to DataFrame object
coordsdf<-data.frame(coordinates(ogShape))

(1)

require(spatstat)
# x-minimum coordinate
xleft <- apply(ogShape@coords,2,min)[1]
# y-minimum coordinate
ybottom <- apply(ogShape@coords,2,min)[2]
# x-maximum coordinate
xright <- apply(ogShape@coords,2,max)[1]
# y-maximum coordinate
ytop <- apply(ogShape@coords,2,max)[2]
# Convert SpatialPointsDataFrame object to DataFrame object
ogpp <- as.ppp(ogShape@coords, c(xleft-1,xright+1,ybottom-1,ytop+1))
# number of x-quadrats
qX <- 200
# number of y-quadrats
qY <- 100
# quadrat size along x-axis
xcellsize <- (xright-xleft)/qX
# quadrat size along y-axis
ycellsize <- (ytop-ybottom)/qY
# doing regular quadrat count
qC <-  quadratcount(ogpp,qX,qY)
# total number of events
n <- nrow(ogShape)

require(DataCombine)

# calling random quadrat count function
randquadcount.data <- randomQuadratCount.ppp(ogpp, coordsdf, qX, qY)

# random quadrat count function. To generate random numbers along x and y-axis and create a quadrat. (x,y) generated is the left-top coordinate.
# We get the rest of the points using xcellsize and ycellsize
randomQuadratCount.ppp <- function(X, coordsdf, nx=5, ny=nx)  {
	# total number of quadrats
	numofquadrats <- nx*ny
	# window object of the study region
	W <- as.owin(X)
	# x-minimum coordinate of the study region
	xmin <- W$xrange[1]  
	# x-maximum coordinate of the study region
	xmax <- W$xrange[2]
	# y-minimum coordinate of the study region
	ymin <- W$yrange[1]
	# y-maximum coordinate of the study region
	ymax <- W$yrange[2]
	# quadrat size in X
	xcellsize <- (xmax-xmin)/nx
	# quadrat size in Y
	ycellsize <- (ymax-ymin)/ny
	# changing ymin value to avoid quadrat going outside the study region
	ymin <- ymin+ycellsize
	# New data frame to store random quadrat count result(set of quadrats)
	randquadcount.data <- data.frame(xmin = c(0), xmax = c(0), ymin = c(0), 
									 ymax = c(0), Freq = c(0))
    # looping over the number of quadrats to be generated 									 
	for (i in 1:numofquadrats){
	  # random number generator along x-axis
	  rx <- runif(1, min=xmin, max=xmax)
	  # random number generator along y-axis
	  ry <- runif(1, min=ymin, max=ymax)
	  # stores list of coordinates of quadrat. (rx,ry) corresponds to left top coordinate of quadrat. Generating rest of the co-ordinates using xcellsize & ycellsize
	  quadrat_list <- list(c(rx,ry),c(rx+xcellsize,ry),c(rx,ry-ycellsize),c(rx+xcellsize,ry-ycellsize))
	  # counter
	  quadcount <- 0
	  # getting only those events which lie between the generated quadrat
	  dfquad <- coordsdf[coordsdf$coords.x1>=rx & coordsdf$coords.x1<=(rx+xcellsize) &
	  coordsdf$coords.x2>=(ry-ycellsize) & coordsdf$coords.x2<=ry,]
	  # looping over the events lying within the quadrat
	  for(j in 1:nrow(dfquad)) {
	    # incrementing the counter varible
		quadcount <- quadcount+1
	  }
	  # creating a new row for a dataframe with the quadrat details
	  New <- c(xmin = quadrat_list[[3]][1], xmax = quadrat_list[[4]][1], ymin = quadrat_list[[3]][2], 
			   ymax = quadrat_list[[1]][2], Freq = c(quadcount))
	  # inserting this row into the resultant data frame  
	  randquadcount.data <- InsertRow(randquadcount.data, NewRow = New, RowNum = i)
	}
	# removing the last NA row
	randquadcount.data <- randquadcount.data[-c(nrow(randquadcount.data)), ]
	# return the resultant data frame  
	return(randquadcount.data)
}					
# plot all the events
plot(coordsdf, pch=20, col="green", main="OilGasLocationsPA",
     xlab="x-coordinate", ylab="y-coordinate")

# call plotting function for regular quadrat method
regplotfn()

library(GISTools)  

# plotting function for regular quadrat method
regplotfn <- function() {
  lx<-xleft
  ly<-ybottom
  # looping over the number of quadrats along y-axis. number of lines to be drawn=number of quadrats+1
  for(i in 1:qY+1){
    # draws lines horizontally
    lines(c(xleft,xright),c(ly,ly))
	# increments y-value
    ly<-ly+ycellsize
  }
  # looping over the number of quadrats along x-axis. number of lines to be drawn=number of quadrats+1
  for(i in 1:qX+1){
    # draws lines vertically
    lines(c(lx,lx),c(ybottom,ytop))
	# increments x-value
    lx<-lx+xcellsize
  }
  # add Legend
  legend(120000, 185000, legend=c("OilGasLocations", "Quadrats"),
         col=c("green", "black"), pch=c(20, 0), cex=0.8,
         title="Legend", text.font=4)
  
  # add North Arrow		 
  north.arrow(100000, 150000, len=10000, lab="N", col="red")
}

plot(coordsdf, pch=20, col="orange", main="OilGasLocationsPA",
     xlab="x-coordinate", ylab="y-coordinate")

# call plotting function for regular quadrat method
randplotfn(randquadcount.data)

# plotting function for random quadrat method
randplotfn <- function(randdf) {
  # loop over the number of rows of the resultant dataframe generated after random quadrat count function
  for(index in 1:nrow(randdf)){
    # create a new dataframe of a particular row
    rowdf <- randdf[index,] 
	# draws line base edge of quadrat
    lines(c(rowdf["xmin"],rowdf["xmax"]),c(rowdf["ymin"],rowdf["ymin"]))
	# draws line left edge of quadrat
    lines(c(rowdf["xmin"],rowdf["xmin"]),c(rowdf["ymin"],rowdf["ymax"]))
	# draws line top edge of quadrat
    lines(c(rowdf["xmin"],rowdf["xmax"]),c(rowdf["ymax"],rowdf["ymax"]))
	# draws line right edge of quadrat
    lines(c(rowdf["xmax"],rowdf["xmax"]),c(rowdf["ymax"],rowdf["ymin"]))
	
	# add Legend
    legend(120000, 185000, legend=c("OilGasLocations", "Quadrats"),
           col=c("orange", "black"), pch=c(20, 0), cex=0.8,
           title="Legend", text.font=4)
		   
    # add North arrow
    north.arrow(100000, 150000, len=10000, lab="N", col="red", tcol="red")
  }
}

(2)

# Regular Quadrant Sampling
# convert the resultant data frame after regular quadrat count method to data frame
dfexhsch <- as.data.frame(qC)
# sort 'Frequency of events in quadrat' column in ascending order 
dfexhsch <- dfexhsch[order(dfexhsch$Freq),] 
# mean=(total number of events)/(total number of quadrats)
Mean <- n/(qX*qY) 
# call function to compute statistics table  
df1.data <- computeStats(dfexhsch, Mean)
View(df1.data)

# Random Quadrat Sampling
# create a new data frame
df2.data <- data.frame(numeric(), numeric(),numeric(), numeric(), numeric())
# sort 'Frequency of events in quadrat' column in ascending order 
randquadcount.data <- randquadcount.data[order(randquadcount.data$Freq),] 
# mean=(total number of events)/(total number of quadrats)
Mean <- n/(qX*qY) 
# call function to compute statistics table  
df2.data <- computeStats(randquadcount.data, Mean)
View(df2.data)

# function to compute the statistics table
computeStats <- function(dfexhsch, Mean)  {
  # create a new data frame
  df1.data <- data.frame(numeric(), numeric(),numeric(), numeric(), numeric())
  i <- 1
  # loop until the number of rows of resultant data frame obtained after quadrat count method
  while (i <= nrow(dfexhsch)) {
    # count variable to keep track of the number of quadrats corresponding to each event
    count <- 1
	# To get the number of quadrats having each event. loop till the 2nd-last row and until the n-th row Frequency column=(n+1)-th row Frequency column
    while(i <= nrow(dfexhsch)-1 & dfexhsch$Freq[i] == dfexhsch$Freq[i+1]){
	  # increment count variable
      count <- count+1
	  # increment i
      i <- i+1
    }  
	# compute (number of events-mean)
    Difference <- (dfexhsch$Freq[i]-Mean)
	# create new row 
    New <- c(dfexhsch$Freq[i], count, Difference, (Difference)^2, 
             count*(Difference)^2)
	# insert this row into df1.data data frame 
    df1.data <- InsertRow(df1.data, NewRow = New)
    i <- i+1
  }
  colnames(df1.data) <- c("Number of Events(K)","Number of Quadrants(X)", "K-Mean",
                          "(K-Mean)^2", "X(K-Mean)^2") 
  return(df1.data)
}				

(3)

# call function to compute VMR value in regular quadrat count method
message("VMR Value of Regular Quadrant Sampling Approach is: ", computeVMR(df1.data, n, qX*qY))
# call function to compute VMR value in random quadrat count method
message("VMR Value of Random Quadrant Sampling Approach is: ", computeVMR(df2.data, n, qX*qY))

# compute VMR function
computeVMR <- function(df.data, n, numquads)  {
  # Mean=number of events/number of quadrats
  Mean <- n/numquads
  # Variance=X*(K-Mean)^2
  Variance <- (sum(df.data[, 5]))/(numquads-1)
  # VMR value =Variance/Mean
  return(Variance/Mean)  
}
