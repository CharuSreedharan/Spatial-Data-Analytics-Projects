library(rgdal)

PALocs <- readOGR(dsn="D:/UPitt/Studies/Sem 2/Spatial DA/Mid-Term/PALocs/PALocs",
                  layer = "PALocs") 

coordsdflocs<-data.frame(coordinates(PALocs))
nrow(PALocs)
plot(coordinates(PALocs))

PACoals <- readOGR(dsn="D:/UPitt/Studies/Sem 2/Spatial DA/Mid-Term/PACoals/PACoals",
                   layer = "PACoals") 
coordsdfcoals<-data.frame(coordinates(PACoals))
nrow(PACoals)
plot(coordinates(PACoals))

View(PALocs@coords)

library(spatstat)

xleftlocs <- apply(PALocs@coords,2,min)[1]
ybottomlocs <- apply(PALocs@coords,2,min)[2]
xrightlocs <- apply(PALocs@coords,2,max)[1]
ytoplocs <- apply(PALocs@coords,2,max)[2]

xleftcoals <- apply(PACoals@coords,2,min)[1]
ybottomcoals <- apply(PACoals@coords,2,min)[2]
xrightcoals <- apply(PACoals@coords,2,max)[1]
ytopcoals <- apply(PACoals@coords,2,max)[2]

GPALocsppp <- as.ppp(PALocs@coords, c(xleftlocs,xrightlocs,
                                      ybottomlocs,ytoplocs))
GPACoalsppp <- as.ppp(PACoals@coords, c(xleftcoals,xrightcoals,
                                        ybottomcoals,ytopcoals))
FPALocsppp <- as.ppp(PALocs@coords, c(xleftlocs,xrightlocs,
                                      ybottomlocs,ytoplocs))
FPACoalsppp <- as.ppp(PACoals@coords, c(xleftcoals,xrightcoals,
                                        ybottomcoals,ytopcoals))
#1)
# number of x-quadrats in PALocs
qXlocs <- 30
# number of y-quadrats in PALocs
qYlocs <- 20

# number of x-quadrats in PACoals
qXcoals <- 100
# number of y-quadrats in PACoals
qYcoals <- 40

nlocs <- nrow(PALocs)
ncoals <- nrow(PACoals)

require(DataCombine)

# calling random quadrat count function for PALocs File
randquadcountlocs.data <- randomQuadratCount.ppp(GPALocsppp, coordsdflocs, qXlocs, qYlocs)

# calling random quadrat count function for PACoals File
randquadcountcoals.data <- randomQuadratCount.ppp(GPACoalsppp, coordsdfcoals, qXcoals, qYcoals)

View(randquadcountlocs.data)

View(randquadcountcoals.data)


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

# plot all the events in PALocs
plot(coordsdflocs, pch=20, col="orange", main="PALocations",
     xlab="x-coordinate", ylab="y-coordinate")

# plot all the events in PACoals
plot(coordsdfcoals, pch=20, col="green", main="PACoals",
     xlab="x-coordinate", ylab="y-coordinate")

library(GISTools)  

# call plotting function for PACoals using regular quadrat count method
randplotfnlocs(randquadcountlocs.data)

# call plotting function for PALocs using regular quadrat method
randplotfncoals(randquadcountcoals.data)

# plotting function for PACoals using random quadrat method
randplotfnlocs <- function(randdf) {
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
    
    #add Legend
    legend(xleftlocs-0.2, ytoplocs, legend=c("PALocations", "Quadrats"),
           col=c("orange", "black"), pch=c(20, 0), cex=0.8,
           title="Legend", text.font=4)

    # add North arrow
    north.arrow(xrightlocs, ybottomlocs+0.2, len=0.1, lab="N", col="red", tcol="red")
  }
}

# plotting function for PACoals using random quadrat count method
randplotfncoals <- function(randdf) {
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
    
    #add Legend
    legend(xrightcoals-0.6, ybottomcoals+0.6, legend=c("PACoals", "Quadrats"),
           col=c("orange", "black"), pch=c(20, 0), cex=0.8,
           title="Legend", text.font=4)
    
    # add North arrow
    north.arrow(xrightcoals-1, ybottomcoals+0.2, len=0.1, lab="N", col="red", 
                tcol="red")
  }
}

# Random Quadrat Sampling PALocs
# create a new data frame
df2locs.data <- data.frame(numeric(), numeric(),numeric(), numeric(), numeric())
# sort 'Frequency of events in quadrat' column in ascending order 
randquadcountlocs.data <- randquadcountlocs.data[order(randquadcountlocs.data$Freq),] 
# mean=(total number of events)/(total number of quadrats)
Mean <- nlocs/(qXlocs*qYlocs) 
# call function to compute statistics table for PALocs 
df2locs.data <- computeStats(randquadcountlocs.data, Mean)
View(df2locs.data)

# Random Quadrat Sampling PACoals
# create a new data frame
df2coals.data <- data.frame(numeric(), numeric(),numeric(), numeric(), numeric())
# sort 'Frequency of events in quadrat' column in ascending order 
randquadcountcoals.data <- randquadcountcoals.data[order(randquadcountcoals.data$Freq),] 
# mean=(total number of events)/(total number of quadrats)
Mean <- ncoals/(qXcoals*qYcoals) 
# call function to compute statistics table for PACoals   
df2coals.data <- computeStats(randquadcountcoals.data, Mean)
View(df2coals.data)

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

# call function to compute VMR value for PALocs using random quadrat count method
message("VMR Value of PALocs using Random Quadrat Sampling Approach is: ", computeVMR(df2locs.data, nlocs, qXlocs*qYlocs))

# call function to compute VMR value for PACoals using random quadrat count method
message("VMR Value of PACoals using Random Quadrat Sampling Approach is: ", computeVMR(df2coals.data, ncoals, qXcoals*qYcoals))

# compute VMR function
computeVMR <- function(df.data, n, numquads)  {
  # Mean=number of events/number of quadrats
  Mean <- n/numquads
  # Variance=X*(K-Mean)^2
  Variance <- (sum(df.data[, 5]))/(numquads-1)
  # VMR value =Variance/Mean
  return(Variance/Mean)  
}

#2)
# To compute G-Function of PALocs
GPALocs <- Gest(GPALocsppp)
#nrow(GPALocs)

# To compute G-Function of PACoals
GPACoals <- Gest(GPACoalsppp)
#nrow(GPACoals)

# To compute F-Function of PALocs
FPALocs <- Fest(FPALocsppp)
#nrow(FPALocs)

# To compute F-Function of PACoals
FPACoals <- Fest(FPACoalsppp)

# To plot the result of G-function of PALocs
plot(GPALocs, main="G-Function of PALocs")

# To plot the result of G-function of PACoals
plot(GPACoals, main="G-Function of PACoals")

# To plot the result of F-function of PALocs
plot(FPALocs, main="F-Function of PALocs")

# To plot the result of F-function of PACoals
plot(FPACoals, main="F-Function of PACoals")




