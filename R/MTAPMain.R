#
# Author: Esad Kopru
#
#

library(ggplot2)
library(ompr)
library(sf)
library(sp)
library(spdep)
library(dplyr)
library(ROI.plugin.glpk)
library(ROI.plugin.symphony)
library(ompr.roi)

set.seed(12345)
# Plot integer axis.
integer_breaks <- function(x){
  seq(floor(min(x)), ceiling(max(x)))}


#' Solve a MTAP .
#'
#' This function solves an MTAP instance
#'
#' @param n size of the surface
#' @param neigType queen-rook-hexagon
#' @param upperLimit upper limit for random cost surface
#'
#' @return MTAPResult
#' @export none
solveMTAPInstance <- function(n, neigType, upperLimit){
  # starting from 0 or 1 ?
  numOfTargetPts <- nrow(targetPoints)
  numOfStartPts <- nrow(startPoint)
  minRow <- 1
  minCol <- 1
  maxRow <- n
  maxCol <- n

  if(tolower(neigType)=='hexagon'){
    ########################### Hexagonal ##############################
    hexaCellCount <- function(i){1+sum(0:i*6)}  # number of cells by hexagonal ring

    nOfRings <- n
    hexaCellCount(nOfRings)

    ## Build hexagonal polygon frame
    tri <- cos(pi/6) # sqrt(1-0.5^2)
    bndHexa <- rbind(c(-1,0),c(-0.5,tri),c(0.5,tri),c(1,0),
                     c(0.5,-tri),c(-0.5,-tri),c(-1,0))

    stBnd <- st_polygon(list(bndHexa))
    stBnd <- st_sfc(stBnd, crs=2960)

    ## Export frame to "sp" object
    spBnd <- as(stBnd, Class = "Spatial")
    proj4string(spBnd)

    ## Generate hexagonal tiles within frame
    hexa <- spsample(spBnd,  type="hexagonal", cellsize=1/(nOfRings+1))
    hexa <- HexPoints2SpatialPolygons(hexa)
    df <- coordinates(hexa)
    df <- data.frame(idx=1:nrow(df),x=round(df[,1], 12),y=round(df[,2],12))
    hexaShp <- SpatialPolygonsDataFrame(hexa, df)
    nb <- spdep::poly2nb(hexaShp, queen=F)       # extract first order neighbors links

    allPoints <<- data.frame(idx = 1:nrow(hexaShp@data), x = hexaShp@data$x, y = hexaShp@data$y, grp="Intermediate")

    # Set Start and Target points
    startPoint <- data.frame(x = allPoints[12,2], y = allPoints[12,3])

    targetPoints <- data.frame(x = allPoints[22,2], y = allPoints[22,3])
    targetPoints <- rbind(targetPoints, c(x = allPoints[31,2], y = allPoints[31,3] ))
    targetPoints <- rbind(targetPoints, c(x = allPoints[41,2], y = allPoints[41,3] ))


    # Group Start, Target and Int nodes for plotting

    for(i in 1:numOfStartPts){
      allPoints[which(allPoints$x == startPoint$x[i] & allPoints$y == startPoint$y[i]), "grp"] <<- "Start"
    }

    for(i in 1:numOfTargetPts){
      allPoints[which(allPoints$x == targetPoints$x[i] & allPoints$y == targetPoints$y[i]), "grp"] <<- "Target"

    }

    # Number of Start and Target points
    numOfTargetPts <- nrow(targetPoints)
    numOfStartPts <- nrow(startPoint)


    hexagonalNeighbors <- function(){
      hexaDF <- data.frame()
      for(i in 1:length(nb)){
        a <- nb[[i]]
        for(j in 1:length(a)){
          hexaDF <- rbind(hexaDF, list(i, a[j], hexaShp@data[i, 2], hexaShp@data[i, 3], hexaShp@data[a[j], 2], hexaShp@data[a[j], 3]))
        }
      }

      names(hexaDF)[1] <- "F"
      names(hexaDF)[2] <- "T"
      names(hexaDF)[3] <- "x"
      names(hexaDF)[4] <- "y"
      names(hexaDF)[5] <- "nX"
      names(hexaDF)[6] <- "nY"

      return(hexaDF)

    }

    ########################### Hexagonal ##############################

    hexaDF <- hexagonalNeighbors()


    impedances <<- hexaDF

    # Total number of variables
    totNV <<- nrow(impedances)
    index <<- 1:totNV
    neigType <- 'hexagon'

    # Create random impedance between 1-20
    impValues <<- sample(1:upperLimit, totNV, replace = T)
    impedances <<- cbind(index, impedances, impValues, impValues)

    flowVars <<- cbind(index,"flow", "f", neigType, hexaDF)

    decisionVars <<- cbind(index,"decision", "d", neigType, hexaDF)

    dummyVars <<- cbind(index, "dummy", "u", neigType, hexaDF)


    ########################### CALCULATE MTAP ###########################

    # Create a new MIPModel()
    MTAPModel <- ompr::MIPModel()

    # Flow variables
    MTAPModel <- add_variable(MTAPModel, f[i], i = 1:totNV, type = "integer", lb=0)

    # Decision variables
    MTAPModel <- add_variable(MTAPModel, d[i], i = 1:totNV, type = "binary")

    # Dummy variables
    MTAPModel <- add_variable(MTAPModel, u[i], i = 1:totNV, type = "binary")


    # minimize travel distance
    MTAPModel <- set_objective(MTAPModel, sum_expr(d[i] * impedances[i, 8]  , i = 1:totNV), "min")

    ####################################################
    # Subject TO:
    ####################################################

    # Start point constraints:  all Flow out = # of target points
    for (i in 1:numOfStartPts) {
      # FROM and TO Start point Indexes
      FROMstartPoint <- flowVars[which(flowVars$x == startPoint$x[i] & flowVars$y == startPoint$y[i]), "index" ]
      MTAPModel <- add_constraint(MTAPModel, sum_expr(f[i], i = FROMstartPoint) == numOfTargetPts)

    }

    ####################################################
    # Target point constraints: all Flow in = 1
    ####################################################

    for(i in 1:numOfTargetPts){
      # To Target Point Indexes
      ToTargetPoint <- flowVars[which(flowVars$nX == targetPoints$x[i] & flowVars$nY == targetPoints$y[i]), "index" ]
      MTAPModel <- add_constraint(MTAPModel, sum_expr(f[i], i = ToTargetPoint ) == 1)

    }

    ####################################################
    ################## Intermediate Nodes ##############
    ####################################################
    for(i in 1:numOfStartPts){
      intNodes <- flowVars
      # Remove To Start Point
      intNodes <- intNodes[-which(flowVars$nX==startPoint$x[i] & flowVars$nY==startPoint$y[i]),]
      nrow(intNodes)
    }


    for(i in 1:numOfTargetPts){
      # Remove From Target Points
      intNodes <- intNodes[-which(intNodes$x == targetPoints$x[i] & intNodes$y == targetPoints$y[i]),]
      nrow(intNodes)
    }


    # Intermediate nodes constraints (excluding start and target points), all Flow in == all Flow out
    IntNodesCondition <- ""
    for(s in 1:numOfStartPts){
      IntNodesCondition <- paste(IntNodesCondition, "(r != startPoint$x[",s,"] | c != startPoint$y[",s,"]) &")
    }
    for(t in 1:numOfTargetPts){
      if( t == numOfTargetPts){
        IntNodesCondition <- paste(IntNodesCondition, "(r != targetPoints$x[",t,"] | c != targetPoints$y[",t,"])")
      }
      else{
        IntNodesCondition <- paste(IntNodesCondition, "(r != targetPoints$x[",t,"] | c != targetPoints$y[",t,"]) &")
      }
    }

    for(row in 1:nrow(intNodes)){
      r <- intNodes$x[row]
      c <- intNodes$y[row]
      if(eval(parse(text=IntNodesCondition))){
        fin <- c(intNodes[(which(intNodes$x==r & intNodes$y == c)), "index"])
        fout <- c(intNodes[(which(intNodes$nX==r & intNodes$nY == c)), "index"])
        MTAPModel <- add_constraint(MTAPModel, sum_expr(f[a], a = fin ) == sum_expr(f[b], b = fout ))
      }
    }

    ####################################################
    # Dummy constraints.
    ####################################################
    MTAPModel <- add_constraint(MTAPModel, f[k] == d[k] + u[k], k = 1:totNV)
    MTAPModel <- add_constraint(MTAPModel, d[k] >= u[k], k = 1:totNV)


    # Solve an instance and calculate solution time
    start_time <- Sys.time()

    # Solution with GLPK solver
    result <- solve_model(MTAPModel, with_ROI(solver = "glpk"))

    return(result)

  }
  ## End of hexagon
  else{

    # Create all possible locations
    createPossibleLocations <- function(){
      xyCmb <- gtools::permutations(n, 2, minRow:maxRow, repeats.allowed=TRUE)
      allPoints <<- data.frame(idx = 1:nrow(xyCmb), x = xyCmb[,1], y = xyCmb[,2], grp="Intermediate")

      # Group Start, Target and Int nodes for plotting

      for(i in 1:numOfStartPts){
        allPoints[which(allPoints$x == startPoint$x[i] & allPoints$y == startPoint$y[i]), "grp"] <<- "Start"
      }

      for(i in 1:numOfTargetPts){
        allPoints[which(allPoints$x == targetPoints$x[i] & allPoints$y == targetPoints$y[i]), "grp"] <<- "Target"

      }

      return(allPoints)

    }

    # (Queen or Rook negihbors)
    neighborsWType <- function(nStyle, varLetter, varType, i, j){
      if((i>=minRow) & (i<=maxRow) & (j>=minCol) & (j<=maxCol)){
        neigs <- data.frame()
        locAndneigs <- data.frame()
        for (a in  as.numeric(i-1):as.numeric(i+1)){
          for (b in as.numeric(j-1):as.numeric(j+1)){
            if((b >= minRow) & (b <= maxRow) & (a >= minRow) & (a <= maxRow) & ((a != i) || (b != j))){
              if(tolower(nStyle) == 'rook'){
                if( !(a %in% c(as.numeric(i-1), as.numeric(i+1) ) )  || !(b %in% c(as.numeric(j-1), as.numeric(j+1)))){
                  # print(data.frame(a,b))
                  locAndneigs <- rbind(locAndneigs, list(varLetter, i, j, a, b))
                  neigs <- rbind(neigs, list(a,b))
                }
              }

              else if(tolower(nStyle) =='queen'){
                #print(data.frame(a,b))
                locAndneigs <- rbind(locAndneigs, list(varLetter, i, j, a, b))
                neigs <- rbind(neigs, list(a,b))
              }

            }
          }
        }
        names(locAndneigs)[1] <- varType
        names(locAndneigs)[2] <- "x"
        names(locAndneigs)[3] <- "y"
        names(locAndneigs)[4] <- "nX"
        names(locAndneigs)[5] <- "nY"

        return(locAndneigs)
      }else{
        print("cannot be lower than min and greater than max")
      }
    }

    # Create neighbors
    recordAllneighbors <- function(varType, varLetter, nStyle){
      allCombined <- data.frame()
      for (a in  minRow:maxRow){
        for (b in minCol:maxCol){
          allCombined <- rbind(allCombined, neighborsWType(nStyle, varLetter, varType, a, b))
        }
      }
      return(allCombined)

    }


    # Create DF for flow, decision, dummy and Impedances
    createVariables <- function(neigType){
      impedances <<- recordAllneighbors("impedance", "I", neigType)

      # Total number of variables
      totNV <<- nrow(impedances)
      index <<- 1:totNV

      # Create random impedance between 1-20
      impValues <<- sample(1:upperLimit, totNV, replace = T)
      impedances <<- cbind(index, impedances, impValues, impValues)

      flowVars <<- recordAllneighbors("flow", "f", neigType)
      flowVars <<- cbind(index, flowVars)

      decisionVars <<- recordAllneighbors("decision", "d", neigType)
      decisionVars <<- cbind(index, decisionVars)

      dummyVars <<- recordAllneighbors("dummy", "u", neigType)
      dummyVars <<- cbind(index, dummyVars)

    }

    allPoints <- createPossibleLocations()

    createVariables(neigType)

    ########################### CALCULATE MTAP ###########################

    # Create a new MIPModel()
    MTAPModel <- ompr::MIPModel()

    # Flow variables
    MTAPModel <- add_variable(MTAPModel, f[i], i = 1:totNV, type = "integer", lb=0)

    # Decision variables
    MTAPModel <- add_variable(MTAPModel, d[i], i = 1:totNV, type = "binary")

    # Dummy variables
    MTAPModel <- add_variable(MTAPModel, u[i], i = 1:totNV, type = "binary")


    # minimize travel distance
    MTAPModel <- set_objective(MTAPModel, sum_expr(d[i] * impedances[i, 8]  , i = 1:totNV), "min")

    ####################################################
    # Subject TO:
    ####################################################

    # Start point constraints:  all Flow out = # of target points
    for (i in 1:numOfStartPts) {
      # FROM and TO Start point Indexes
      FROMstartPoint <- flowVars[which(flowVars$x == startPoint$x[i] & flowVars$y == startPoint$y[i]), "index" ]
      # TOstartPoint <- flowVars[which(flowVars$nX == startPoint$x[i] & flowVars$nY == startPoint$y[i]), "index" ]
      MTAPModel <- add_constraint(MTAPModel, sum_expr(f[i], i = FROMstartPoint) == numOfTargetPts)

    }

    ####################################################
    # Target point constraints: all Flow in = 1
    ####################################################

    for(i in 1:numOfTargetPts){
      # To Target Point Indexes
      ToTargetPoint <- flowVars[which(flowVars$nX == targetPoints$x[i] & flowVars$nY == targetPoints$y[i]), "index" ]
      MTAPModel <- add_constraint(MTAPModel, sum_expr(f[i], i = ToTargetPoint ) == 1)

    }

    ####################################################
    ################## Intermediate Nodes ##############
    ####################################################
    for(i in 1:numOfStartPts){
      intNodes <- flowVars
      # Remove To Start Point
      intNodes <- intNodes[-which(flowVars$nX==startPoint$x[i] & flowVars$nY==startPoint$y[i]),]
      nrow(intNodes)
    }


    for(i in 1:numOfTargetPts){
      # Remove From Target Points
      intNodes <- intNodes[-which(intNodes$x == targetPoints$x[i] & intNodes$y == targetPoints$y[i]),]
      nrow(intNodes)
    }


    # Intermediate nodes constraints (excluding start and target points), all Flow in == all Flow out
    IntNodesCondition <- ""
    for(s in 1:numOfStartPts){
      IntNodesCondition <- paste(IntNodesCondition, "(r != startPoint$x[",s,"] | c != startPoint$y[",s,"]) &")
    }
    for(t in 1:numOfTargetPts){
      if( t == numOfTargetPts){
        IntNodesCondition <- paste(IntNodesCondition, "(r != targetPoints$x[",t,"] | c != targetPoints$y[",t,"])")
      }
      else{
        IntNodesCondition <- paste(IntNodesCondition, "(r != targetPoints$x[",t,"] | c != targetPoints$y[",t,"]) &")
      }
    }

    for(r in minRow:maxRow){
      for(c in minCol:maxCol){
        if(eval(parse(text=IntNodesCondition))){
          # print(c(r, c))
          fin <- c(intNodes[(which(intNodes$x==r & intNodes$y == c)), "index"])
          fout <- c(intNodes[(which(intNodes$nX==r & intNodes$nY == c)), "index"])
          MTAPModel <- add_constraint(MTAPModel, sum_expr(f[a], a = fin ) == sum_expr(f[b], b = fout ))
        }
      }
    }

    ####################################################
    # Dummy constraints.
    ####################################################
    MTAPModel <- add_constraint(MTAPModel, f[k] == d[k] + u[k], k = 1:totNV)
    MTAPModel <- add_constraint(MTAPModel, d[k] >= u[k], k = 1:totNV)

    # Solution with GLPK solver
    result <- solve_model(MTAPModel, with_ROI(solver = "glpk"))
    print(MTAPModel)

    return(result)
  }

}


#
#' Plot the result of MTAP .
#'
#' This function plots the solution path
#'
#' @param MTAPResult result of MTAP
#' @param neigType queen-rook-hexagon
#'
#' @return none
#' @export none
#'
plotMTAPResults<- function(MTAPResult, neigType){
  solution <- get_solution(MTAPResult, d[i]) %>%
    filter(value > 0)

  d <- data.frame(impedances[solution$i,])
  d2 <- data.frame(index=d$index, x=d$x-0.5, y= d$y-0.5, nX=d$nX-0.5, nY=d$nY-0.5, impValues = d$impValues)

  ##################
  # Plot the solution
  if(tolower(neigType)=='hexagon'){
    xOffset <- 0
    yOffset <- 0.07
  }else{
    xOffset <- 0.1
    yOffset <- 0.1
  }


  allP2 <- data.frame(idx=allPoints$idx, x=allPoints$x-0.5, y=allPoints$y-0.5, grp=allPoints$grp)

  # Result
  ggplot() +
    geom_point(data=allP2, aes(x, y, color=grp),size = 3) +
    scale_y_continuous(minor_breaks = integer_breaks) +
    scale_x_continuous(minor_breaks = integer_breaks) +
    # use minor_breaks for background grid lines.
    scale_color_manual(values=c("green", "red", "gray"), breaks = c("Start", "Target", "Intermediate")) +
    labs(title = paste0("Optimal Cost : ", sum(d$impValues)), subtitle = paste0("Computation Time : ","compTime", " secs"))+
    geom_segment(data=d2, mapping = aes(x=x, y=y, xend=nX, yend=nY), col='blue', size=0.6)+
    geom_text(data=d2, mapping = aes(x=nX+xOffset, y=nY+yOffset,label=impValues) )+
    theme(panel.background = element_rect(fill = "lightgray"))
}


#
#' Reads a shapefile and creates a topology based on its dataframe .
#'
#'
#' @param MTAPResult result of MTAP
#' @param neigType queen-rook
#'
#' @return none
#' @export none

MTAPonShapefile <- function(shpFilePath, neigType, howManyTargetPts, upperLimit){
  ##
  ## Read Shapefiles
  ##

  inShp <- rgdal::readOGR(dsn=shpFilePath, integer64="warn.loss")
  inShp.centroid <- coordinates(inShp)        # Get province centroids

  # Add lat, long to input shapefile.
  inShp$x <- inShp.centroid[,1]
  inShp$y <- inShp.centroid[,2]


  createNeighbors <- function(inShp, neigType, upperLimit){
    proj4string(inShp)                          # map projection
    inShp.bbox <- bbox(inShp)           # province bounding box for map region

    isQueen = F
    if(tolower(neigType)=='queen'){
      isQueen = T
    }

    inShp.link <- poly2nb(inShp, queen=isQueen)
    # Create neighborhood lists
    neighDF <- data.frame()
    for(i in 1:length(inShp.link)){
      a <- inShp.link[[i]]
      if(a != 0){
        for(j in 1:length(a)){
          neighDF <- rbind(neighDF, list(i, a[j], inShp@data$x[i] ,inShp@data$y[i]
                                         ,inShp@data$x[a[j]] ,inShp@data$y[a[j]]))
        }
      }

    }
    names(neighDF)[1] <- "F"
    names(neighDF)[2] <- "T"
    names(neighDF)[3] <- "x"
    names(neighDF)[4] <- "y"
    names(neighDF)[5] <- "nX"
    names(neighDF)[6] <- "nY"

    return(neighDF)
  }


  # howManyTargetPts <- 4
  # Record all points
  sizeOfinput <- nrow(inShp@data)
  allPoints <- data.frame(idx = 1:sizeOfinput, x = inShp@data$x, y = inShp@data$y, grp="Intermediate")

  # decide on start and target points
  # Set Start and Target points
  # If Texas
  if(length(inShp@polygons) == 254 && howManyTargetPts == 4){
    # 226 - houston
    # 216 - austin
    # 130 - san saba
    # 155 - wichita
    # 239 - dallas

    # Record all points
    sizeOfinput <- nrow(inShp@data)
    allPoints <- data.frame(idx = 1:sizeOfinput, x = inShp@data$x, y = inShp@data$y, grp="Intermediate")

    startPoint <- data.frame()
    startPoint <- data.frame(x = allPoints[226,2], y = allPoints[226,3])
    targetPoints <- data.frame()
    targetPoints <- data.frame(x = allPoints[216,2], y = allPoints[216,3])
    targetPoints <- rbind(targetPoints, c(x = allPoints[130,2], y = allPoints[130,3]))
    targetPoints <- rbind(targetPoints, c(x = allPoints[155,2], y = allPoints[155,3]))
    targetPoints <- rbind(targetPoints, c(x = allPoints[239,2], y = allPoints[239,3]))

  }

  else if(length(inShp@polygons) == 254 && howManyTargetPts == 3){
    # 226 - houston
    # 216 - austin
    # 130 - san saba
    # 155 - wichita
    # 239 - dallas

    # Record all points
    sizeOfinput <- nrow(inShp@data)
    allPoints <- data.frame(idx = 1:sizeOfinput, x = inShp@data$x, y = inShp@data$y, grp="Intermediate")

    startPoint <- data.frame()
    startPoint <- data.frame(x = allPoints[226,2], y = allPoints[226,3])
    targetPoints <- data.frame()
    targetPoints <- data.frame(x = allPoints[216,2], y = allPoints[216,3])
    targetPoints <- rbind(targetPoints, c(x = allPoints[130,2], y = allPoints[130,3]))
    targetPoints <- rbind(targetPoints, c(x = allPoints[155,2], y = allPoints[155,3]))

  }

  else {
    tNum <- 32
    targetPoints <- data.frame()
    startPoint <- data.frame(x = allPoints[5,2], y = allPoints[5,3])
    for(i in 1:howManyTargetPts){
      tps <- floor(sizeOfinput/tNum)
      if(tps == 5){
        tps <- tps + 3
      }
      targetPoints <- rbind(targetPoints, c(x = allPoints[tps,2], y = allPoints[tps,3]))
      tNum <- tNum/2
    }
    names(targetPoints)[1] <- "x"
    names(targetPoints)[2] <- "y"
  }


  # Number of Start and Target points
  numOfTargetPts <- nrow(targetPoints)
  numOfStartPts <- nrow(startPoint)


  # Group Start, Target and Int nodes for plotting
  for(i in 1:numOfStartPts){
    allPoints[which(allPoints$x == startPoint$x[i] & allPoints$y == startPoint$y[i]), "grp"] <- "Start"
  }

  for(i in 1:numOfTargetPts){
    allPoints[which(allPoints$x == targetPoints$x[i] & allPoints$y == targetPoints$y[i]), "grp"] <- "Target"

  }

  shapeDF <- createNeighbors(inShp, neigType, upperLimit)
  nrow(shapeDF)
  impedances <- shapeDF

  # Total number of variables
  totNV <- nrow(impedances)
  index <- 1:totNV

  # Create random impedance between 1-20
  impValues <- sample(1:upperLimit, totNV, replace = T)
  impedances <- cbind(index, impedances, impValues, impValues)

  flowVars <- cbind(index,"flow", "f", neigType, shapeDF)

  decisionVars <- cbind(index,"decision", "d", neigType, shapeDF)

  dummyVars <- cbind(index, "dummy", "u", neigType, shapeDF)


  ########################### CALCULATE MTAP ###########################

  # Create a new MIPModel()
  MTAPModel <- ompr::MIPModel()

  # Flow variables
  MTAPModel <- add_variable(MTAPModel, f[i], i = 1:totNV, type = "integer", lb=0)

  # Decision variables
  MTAPModel <- add_variable(MTAPModel, d[i], i = 1:totNV, type = "binary")

  # Dummy variables
  MTAPModel <- add_variable(MTAPModel, u[i], i = 1:totNV, type = "binary")


  # minimize travel distance
  MTAPModel <- set_objective(MTAPModel, sum_expr(d[i] * impedances[i, 8]  , i = 1:totNV), "min")

  ####################################################
  # Subject TO:
  ####################################################

  # Start point constraints:  all Flow out = # of target points
  for (i in 1:numOfStartPts) {
    # FROM and TO Start point Indexes
    FROMstartPoint <- flowVars[which(flowVars$x == startPoint$x[i] & flowVars$y == startPoint$y[i]), "index" ]
    # TOstartPoint <- flowVars[which(flowVars$nX == startPoint$x[i] & flowVars$nY == startPoint$y[i]), "index" ]
    MTAPModel <- add_constraint(MTAPModel, sum_expr(f[i], i = FROMstartPoint) == numOfTargetPts)

  }

  ####################################################
  # Target point constraints: all Flow in = 1
  ####################################################

  for(i in 1:numOfTargetPts){
    # To Target Point Indexes
    ToTargetPoint <- flowVars[which(flowVars$nX == targetPoints$x[i] & flowVars$nY == targetPoints$y[i]), "index" ]
    MTAPModel <- add_constraint(MTAPModel, sum_expr(f[i], i = ToTargetPoint ) == 1)

  }

  ####################################################
  ################## Intermediate Nodes ##############
  ####################################################
  for(i in 1:numOfStartPts){
    intNodes <- flowVars
    # Remove To Start Point
    intNodes <- intNodes[-which(flowVars$nX==startPoint$x[i] & flowVars$nY==startPoint$y[i]),]
    nrow(intNodes)
  }


  for(i in 1:numOfTargetPts){
    # Remove From Target Points
    intNodes <- intNodes[-which(intNodes$x == targetPoints$x[i] & intNodes$y == targetPoints$y[i]),]
    nrow(intNodes)
  }


  # Intermediate nodes constraints (excluding start and target points), all Flow in == all Flow out
  IntNodesCondition <- ""
  for(s in 1:numOfStartPts){
    IntNodesCondition <- paste(IntNodesCondition, "(r != startPoint$x[",s,"] | c != startPoint$y[",s,"]) &")
  }
  for(t in 1:numOfTargetPts){
    if( t == numOfTargetPts){
      IntNodesCondition <- paste(IntNodesCondition, "(r != targetPoints$x[",t,"] | c != targetPoints$y[",t,"])")
    }
    else{
      IntNodesCondition <- paste(IntNodesCondition, "(r != targetPoints$x[",t,"] | c != targetPoints$y[",t,"]) &")
    }
  }

  for(row in 1:nrow(intNodes)){
    r <- intNodes$x[row]
    c <- intNodes$y[row]
    if(eval(parse(text=IntNodesCondition))){
      fin <- c(intNodes[(which(intNodes$x==r & intNodes$y == c)), "index"])
      fout <- c(intNodes[(which(intNodes$nX==r & intNodes$nY == c)), "index"])
      MTAPModel <- add_constraint(MTAPModel, sum_expr(f[a], a = fin ) == sum_expr(f[b], b = fout ))
    }
  }

  ####################################################
  # Dummy constraints.
  ####################################################
  MTAPModel <- add_constraint(MTAPModel, f[k] == d[k] + u[k], k = 1:totNV)
  MTAPModel <- add_constraint(MTAPModel, d[k] >= u[k], k = 1:totNV)


  # Solve an instance and calculate solution time
  start_time <- Sys.time()

  # Solution with GLPK solver
  result <- solve_model(MTAPModel, with_ROI(solver = "glpk"))
  end_time <- Sys.time()
  compTime <- round(difftime(end_time, start_time, units="secs"), 4)



  solution <- get_solution(result, d[i]) %>%
    filter(value > 0)

  d <- data.frame(impedances[solution$i,])

  d$grp <- allPoints[d$F,]$grp

  d

  ##################
  # Plot the solution
  if(tolower(neigType)=='hexagon'){
    xOffset <- 0
    yOffset <- 0.07
  }else{
    xOffset <- 0.1
    yOffset <- 0.1
  }


  # Alternative plotting with arrows - does not look as attractive
  plot(inShp,col="palegreen3" ,border=grey(0.9), axes=T ) # Second plot areas
  arrows(x0 = d$x, x1 = d$nX, y0 = d$y, y1 = d$nY, length = 0.3, angle = 30 , col="blue")


  # Plot the results
  map <- ggplot() +
    geom_polygon(data = inShp, aes(x = long, y = lat, group = group), colour = "gray", fill = NA) +
    geom_point(data=targetPoints, aes(x, y), color="red",size = 4) +
    geom_point(data=startPoint, aes(x, y), color="green",size = 4) +
    labs( title = paste0("Optimal Cost : ", sum(d$impValues)), subtitle = paste0("Computation Time : ",compTime, " secs"))+
    geom_segment(data=d, mapping = aes(x=x, y=y, xend=nX, yend=nY), arrow=arrow(length=unit(0.21,"cm"), ends="last", type = "closed"), col="blue", size=0.1) +
    theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1))

  map
}

