### Geomorphic Approach ###

#' A function to execute the Geomorphic Approach hydraulic simulation routine
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

AvgHydraulics2 <- function(S, wb, db, db_max = NULL, b_value = NULL, max_Q = 1,
                             D84, xs_output = TRUE) {

  # Frank - I know that loading packages within a function is bad practice, please
  # advise how to address.
  # load libraries
  library(dplyr)
  library(zoo)

  ###########################################

  ##### Define Find_U Function #####
  # this will be later used to simulate velocities at different water laters
  findU <- function(Wb, S, D84, depths) {

    deltaX <- 0.0001
    Xgrid <- Wb * seq(0, 1, deltaX)

    wet.vert <- depths[depths >= 0]
    Wi <- length(wet.vert) * deltaX * Wb
    Ai <- sum(wet.vert * deltaX * Wb)
    di <- Ai / Wi
    Pi <- sum((diff(wet.vert)^2 + (max(Xgrid) * deltaX) ^2) ^ (1/2))
    Ri <- Ai/Pi

    # Ferguson's continuously varying power law. See Ferguson (2007).
    D.84 <- D84 / 1000 #D84 grain size in m
    g <- 9.81 # gravity
    a1 <- 6.5
    a2 <- 2.5
    Res <- a1 * a2 * (Ri / D.84) /
      (a1^2 + a2^2 * (Ri / D.84) ^ (5/3)) ^ (1/2)
    Ui <- Res * sqrt(g * Ri * S) # Velocity (m/s)

    # Formatting the outputs in a dataframe
    df <- data.frame(Ai, Wi, di, Ui)

    return(df)
  }

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

  ### Frank - please change below stop statement
  # stop function execution if error message too high
  try(if(b > 0.7) stop("Error: model will not produce realistic results because b-value unrealistically high"))

  # define grid
  deltaX <- 0.0001
  deltaY <- 0.001  # increment by which to change depths when estimating HG

  # estimate max depth using b-value
  dmax <- (1 + (b / (1 - b))) * db

  # generate xs_corrdinates
  X <- c(0, b * wb, 0.99 * wb, wb)
  Y <- 5 * db- c(0, db, dmax, 0) # depths are relative

  # Interpolate the distribution onto an xs raster
  Xgrid <- wb * seq(0, 1, deltaX)
  Ygrid <- matrix(unlist(approx(X, Y, Xgrid)), ncol = 10001, byrow = TRUE)[2,]

  # Specify water surface elevations for which to calculate Wi
  Zw <- 5 * db - dmax + seq(0.02 * dmax, dmax, deltaY * dmax)

  ######################################################
  # For loop to calculate the width and discharge for each chosen water level

  # create objects to hold store results
  simulated <- data.frame(Q = NA, Ai = NA, Wi = NA, di = NA, Ui = NA)
  results <- list()

  for (j in 1:length(Zw)) {
    #j = 20
    depths <- Zw[j] - Ygrid # Calculate the depths, for each vertical
    results <- findU(wb, S, D84, depths) # calculate hydraulics for each vertical
    results <- c(Q = results[1, 4] * results[1, 1], results) # calculate discarhge
    simulated[j, ] <- results
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


