### Geomorphic Approach ###

#' A function to execute the Geomorphic Approach hydraulic simulation for flows below bankfull.
#'
#' This function executes the hydraulic geometry simulator to evaluate reach-averaged depths and velocities
#' generated at flows less than bankfull conditions.
#' For more information about this model see: McParland et al. (2016) and Gronsdahl et al. (XXXX)

#' @param S channel gradient (m/m)
#' @param wb reach averaged bankfull width (m).
#' @param db reach averaged bankfull depth (m).
#' @param db_max Defaults to NULL. reach averaged maximum bankfull depth (m).  Specifying db_max is preferred to calculate 'b'.
#' @param b User-specified b-value. Defaults to NULL and calculated within model unless specified.
#' @param max_Q maximum discharge (m3/s) to simualted WUA for.  Defaults to 1 m3/s.
#' @param D84 grain size (mm)
#' @param xs_output Defaults to TRUE. An expression specifying whether to produce a .csv and .jpg of the simulated channel cross section.
#' @export
#' @return .csv and .jpeg of channel cross section if specified
#' @return data frame of reach-averaged hydraulics
#' AvgHydraulics()

library(dplyr)
library(zoo)

calcWidth <- function(elev, b, wb, db, dmax) {

  #To calculate width, the channel is split into 3 segments: left bank, channel and right bank
  #using Figure 1 in McParland, D., Eaton, B., and Rosenfeld, J. 2016

  #The widths are calculated using the b value following the 2016 paper
  width.leftbank = b * wb
  width.channel = (0.99-b) * wb
  width.rightbank = 0.01*wb

  #The depths are defined as db for left bank and dmax for right bank
  #The channel depths is the difference between the left and right bank depths

  depth.leftbank = db
  depth.channel = dmax-db
  depth.rightbank = dmax

  #there are two different ways to calculate depth, depending on if the channel is fully submerged

  if(elev >= depth.channel) {

    #channel is fully submerged

    #The wetted width is calculated using trigonometry
    #Each segment of the channel is a triangle
    #We denote the depth a (opposite) in the triangle
    #width is b (adjacent) in the triangle
    #and the channel bottom is c (hypotenuse) in the triangle
    #From the tangent equations in trigonometry, we now that the ratio between the straight lines are constant
    #i.e. a1/b1 = a2/b2
    #in our case, that means the ratio of the depth to width is the same for the geometry and wetted parts
    #so we get segment depth/segment width = wetted depth/wetted width
    #We know the wetted depth, so we reorder to get wetted width = segment width * wetted depth/segment depth
    #The wetted width becomes a fraction of the segment width, the same ratio as the wetted depth to the segment depth

    wetted.depth.leftbank = elev - depth.channel
    wetted.width.leftbank = width.leftbank*wetted.depth.leftbank/depth.leftbank

    wetted.width.channel = width.channel #the entire channel is submerged, so the wetted width is the segment width

    wetted.depth.rightbank = elev
    wetted.width.rightbank = width.rightbank*wetted.depth.rightbank/depth.rightbank

  } else {

    wetted.width.leftbank = 0 #the entire left bank is above water, so wetted width is zero

    wetted.depth.channel = elev
    wetted.width.channel = width.channel*wetted.depth.channel/depth.channel

    wetted.depth.rightbank = elev
    wetted.width.rightbank = width.rightbank*wetted.depth.rightbank/depth.rightbank
  }

  wetted.width = wetted.width.leftbank + wetted.width.channel + wetted.width.rightbank

  return(wetted.width)
}

calcArea <- function(elev, b, wb, db, dmax) {

  #The area calculation is using the same trigonometry as calcWidth to calculated wetted width
  #The area is calculate based on the area of a triangle (width*depth/2)

  width.leftbank = b * wb
  width.channel = (0.99-b) * wb
  width.rightbank = 0.01*wb

  depth.leftbank = db
  depth.channel = dmax-db
  depth.rightbank = dmax

  if(elev >= depth.channel) {

    #channel is fully submerged

    wetted.depth.leftbank = elev - depth.channel
    wetted.width.leftbank = width.leftbank*wetted.depth.leftbank/depth.leftbank
    wetted.area.leftbank = wetted.width.leftbank*wetted.depth.leftbank/2

    #the entire channel is submerged, so the wetted area is the triangle for the width and depth of the channel
    #and the rectangle from the channel depth up to the water elevation
    wetted.area.channel = width.channel*depth.channel/2 + width.channel*(elev-depth.channel)

    wetted.depth.rightbank = elev
    wetted.width.rightbank = width.rightbank*wetted.depth.rightbank/depth.rightbank
    wetted.area.rightbank = wetted.width.rightbank*wetted.depth.rightbank/2

  } else {

    wetted.area.leftbank = 0 #the entire left bank is above water, so zero area

    wetted.depth.channel = elev
    wetted.width.channel = width.channel*wetted.depth.channel/depth.channel
    wetted.area.channel = wetted.width.channel*wetted.depth.channel/2

    wetted.depth.rightbank = elev
    wetted.width.rightbank = width.rightbank*wetted.depth.rightbank/depth.rightbank
    wetted.area.rightbank = wetted.width.rightbank*wetted.depth.rightbank/2

  }

  wetted.area = wetted.area.leftbank + wetted.area.channel + wetted.area.rightbank

  return(wetted.area)
}

calcP <- function(elev, b, wb, db, dmax) {

  #The bottom length is calculated using the Pythagorean theorem: a^2+b^2=c^2
  #For our purpose, that means wetted depth ^ 2 + wetted width ^ 2 = wetted bottom length ^ 2
  #Reorderede get wetted bottom length = sqrt(wetted depth ^ 2 + wetted width ^ 2)

  #wetted width are calculated as in calcWidth

  width.leftbank = b * wb
  width.channel = (0.99-b) * wb
  width.rightbank = 0.01*wb

  depth.leftbank = db
  depth.channel = dmax-db
  depth.rightbank = dmax

  if(elev >= depth.channel) {

    #channel is fully submerged

    wetted.depth.leftbank = elev - depth.channel
    wetted.width.leftbank = width.leftbank*wetted.depth.leftbank/depth.leftbank
    length.bottom.leftbank = sqrt(wetted.width.leftbank^2+wetted.depth.leftbank^2)

    #the entire channel is submerged, so wetted width and depth is the same as channel width and depth
    length.bottom.channel = sqrt(width.channel^2+depth.channel^2)

    wetted.depth.rightbank = elev
    wetted.width.rightbank = width.rightbank*wetted.depth.rightbank/depth.rightbank
    length.bottom.rightbank = sqrt(wetted.width.rightbank^2+wetted.depth.rightbank^2)

  } else {

    length.bottom.leftbank = 0 #the entire left bank is above water

    wetted.depth.channel = elev
    wetted.width.channel = width.channel*wetted.depth.channel/depth.channel
    length.bottom.channel = sqrt(wetted.width.channel^2+wetted.depth.channel^2)

    wetted.depth.rightbank = elev
    wetted.width.rightbank = width.rightbank*wetted.depth.rightbank/depth.rightbank
    length.bottom.rightbank = sqrt(wetted.width.rightbank^2+wetted.depth.rightbank^2)

  }

    p = length.bottom.leftbank + length.bottom.channel + length.bottom.rightbank

  return(p)
}

# Velocity function per Ferguson 2007
calcUi <- function(Ri, D84, S) {
  D.84 <- D84 / 1000 #D84 grain size in m
  g <- 9.81 # gravity
  a1 <- 6.5
  a2 <- 2.5
  Res <- a1 * a2 * (Ri / D.84) /
    (a1^2 + a2^2 * (Ri / D.84) ^ (5/3)) ^ (1/2)
  Ui <- Res * sqrt(g * Ri * S) # Velocity (m/s)
  return(Ui)
}

AvgHydraulics <- function(S, wb, db, db_max = NULL, b_value = NULL, max_Q = 1,
                            D84, xs_output = TRUE) {

  #####################################################
  # input validation
  if(!is.numeric(S)) {stop("AvgHydraulics expects 'S' to be numeric")}
  if(S <= 0 ) {stop("AvgHydraulics expects 'S' to be a positive number")}

  if(!is.numeric(wb)) {stop("AvgHydraulics expects 'wb' to be numeric")}
  if(wb <= 0 ) {stop("AvgHydraulics expects 'wb' to be a positive number")}
  if(wb > 100 ) {warning("Warning: width outside of recommended range")}

  if(!is.numeric(db)) {stop("AvgHydraulics expects 'db' to be numeric")}
  if(db <= 0 ) {stop("AvgHydraulics expects 'db' to be a positive number")}
  if(db > 5 ) {warning("Warning: db outside of recommended range")}

  if(!is.numeric(D84)) {stop("AvgHydraulics expects 'D84' to be numeric")}
  if(D84 <= 1 ) {warning("Warning: Confirm D84 is entered in mm")}
  if(D84 > 400 ) {warning("Warning: D84 may be outside of recommended range")}

  if(is.null(db_max)){
  } else {
    if(db > db_max) {stop("AvgHydraulics expects db < db_max")}
  }

  deltaX = 0.0001
  deltaY = 0.001
  #############################################################
  ##### Simulate Hydraulics #####

  # Ferguson model's shape factor (b): define based on specified inputs
  if(is.null(b_value) == FALSE){
    b <- b_value # option 1: user specified
  } else if (is.null(db_max) == FALSE) {
    b <- 1 - (db / db_max) # option 2: incorporates db_max values
  } else {
    b <- (wb / db) / 100 # option 3: uses mean width
  }
  # b value input validation
  if(b > 0.7 ) {warning("Warning: b_value outside of recommended range")}

  # estimate max depth using b-value
  dmax <- (1 + b) / (1 - b) * db

  # generate xs_corrdinates
  X <- c(0, b * wb, 0.99 * wb, wb)
  Y <- 5 * db- c(0, db, dmax, 0) # depths are relative

  # Interpolate the distribution onto an xs raster
  Xgrid <- wb * seq(0, 1, deltaX)
  Ygrid <- matrix(unlist(approx(X, Y, Xgrid)), ncol = length(Xgrid), byrow = TRUE)[2,]

  # Specify water surface elevations for which to calculate Wi
  Zw <- 5 * db - dmax + seq(0.02 * dmax, dmax, deltaY * dmax)

  ######################################################
  # For loop to calculate the width and discharge for each chosen water level

  # create objects to hold store results
  simulated <- data.frame(Zw=Zw, Q = NA, Ai = NA, Wi = NA, di = NA, Ui = NA, elev = NA)
  results <- list()

  for (j in 1:length(Zw)) {
    #j = 20
    elev = Zw[j] - min(Ygrid)
    Wi = calcWidth(elev, b, wb, db, dmax)
    Ai = calcArea(elev, b, wb, db, dmax)
    Pi = calcP(elev, b, wb, db, dmax)

    di <- Ai / Wi
    Ri <- Ai/Pi
    Ui = calcUi(Ri, D84, S)
    Q = Ui * Ai

    simulated$elev[j] = elev
    simulated$Wi[j] = Wi
    simulated$Ai[j] = Ai
    simulated$Pi[j] = Pi
    simulated$di[j] = di
    simulated$Ri[j] = Ri
    simulated$Ui[j] = Ui
    simulated$Q[j] = Q
  }

  # set up data frame of outputs with varying streamflow intervals
  Q <- c(seq(0.001, 0.1, 0.001), seq(0.11, 1, 0.01), seq(1.1, 10, 0.1),
         seq(11, 100, 1), seq(110, 1000, 10), seq(1100, 10000, 100))

  # add modelled hydraulics to output dataframe
  Ai <- approx(simulated$Q, simulated$Ai, xout = Q)[2]
  Wi <- approx(simulated$Q, simulated$Wi, xout = Q)[2]
  di <- approx(simulated$Q, simulated$di, xout = Q)[2]
  Ui <- approx(simulated$Q, simulated$Ui, xout = Q)[2]

  # filter results
  mod_hyd <- data.frame(Q, Ai = Ai$y, Wi = Wi$y, di = di$y, Ui = Ui$y) %>%
    filter(is.na(Ai) == FALSE) %>% filter(Q <= max_Q)

  #####################################################
  # Prepare graph of cross section
  if(xs_output == TRUE){

    # output coordinates
    if(is.null(db_max) == TRUE){
      plot_y <- c(0, (db * -1), (dmax * - 1), 0)
    } else {
      plot_y <- c(0, (db * -1), (db_max * - 1), 0)
    }

    # set up x-values to plot
    plot_x <- c(0, (b * wb), (0.99 * wb), wb )

    # write channel cross section
    channel_xs <- data.frame(x = plot_x, y = plot_y)
    write.csv(channel_xs, "channel_xs.csv", row.names = FALSE)

    # plot simple figure
    jpeg("channel_xs.jpeg", width = 6, height = 4, units = "in", res = 300)
    par(mar = c(4.5, 4.5, 1, 1))
    plot(plot_x, plot_y, type = "l",
         xlab = "Width (m)",
         ylab = "Depth (m)",
         ylim = c((min(plot_y) * 1.2), 0), cex.lab = 0.8, cex.axis = 0.8,
    )
    abline(a = (db * -1), 0, lty = 2, col = "grey")
    legend("bottomleft", col = c("black", "grey"), bty = "n",
           lty = c(1, 2), cex = 0.86,
           legend = c("Channel cross section", "Average depth"))
    dev.off()
  } else {
  }

  # return modelled hydraulics
  return(mod_hyd)

}

