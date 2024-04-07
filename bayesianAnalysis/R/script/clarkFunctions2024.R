
# TOF

rampBrown2Green <- c('#8c510a','#bf812d','#dfc27d','#80cdc1','#35978f','#01665e') # color ramp
rampRed2Blue    <- c('#d73027','#f46d43','#fdae61','#fee090','#abd9e9','#74add1','#4575b4')
rampYell2Blue   <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#253494','#081d58')
rampRed2Gray    <- c('#b2182b','#d6604d','#f4a582','#fddbc7','#ffffff','#e0e0e0','#bababa',
                     '#878787','#4d4d4d')
rampWhite2Green <- c('#ffffe5','#f7fcb9','#d9f0a3','#addd8e','#78c679',
                     '#41ab5d','#238443','#005a32')


points2grid <- function( x, y, z, FUN = 'mean', nx = 50, ny = 50, 
                         xlim = NULL, ylim = NULL, LOG = F ){
  
  # FUN can be 'mean' or 'sum'
  
  KEEPX <- KEEPY <- T
  
  if( is.null(xlim) ){
    xlim <- quantile(x, c(.01, .99), na.rm=T )
    KEEPX <- F
  }
  if( is.null(ylim) ){
    ylim <- quantile(y, c(.01, .99), na.rm=T )
    KEEPY <- F
  }
  
  xseq <- seq( xlim[1], xlim[2], length = nx )
  yseq <- seq( ylim[1], ylim[2], length = ny )
  
  if( !KEEPX ){
    xnew <- xseq
    nnx <- nx
    d1 <- 10
    
    while( nnx == nx ){
      nnx <- length( unique( round(xnew, d1 ) ) )
      d1  <- d1 - 1
    }
    dx <- min( diff( round( xseq, d1+2 ) ) )
    xseq <- sort( round( seq( xlim[1], xlim[2], by = dx ), d1+2 ) )
    nx   <- length(xseq)
  }
  
  if( !KEEPY ){
    ynew <- yseq
    nny <- ny
    d2 <- 10
    while( nny == ny ){
      nny <- length( unique( round(ynew, d2 ) ) )
      d2  <- d2 - 1
    }
    dy <- min( diff( round( yseq, d2+2 ) ) )
    yseq <- sort( round( seq( ylim[1], ylim[2], by = dy ), d2+2 ) )
    ny   <- length(yseq)
  }
  
  ix <- findInterval( x, xseq, all.inside = T )
  iy <- findInterval( y, yseq, all.inside = T )
  
  if( FUN == 'mean' ){
    z0 <- tapply( z, list(x = ix, y = iy), mean, na.rm=T)
  }else{
    z0 <- tapply( z, list(x = ix, y = iy), sum, na.rm=T)
  }
  
  wx <- sort( unique( which( !c(1:nx) %in% as.numeric( rownames(z0) ) ) ) )
  wy <- sort( unique( which( !c(1:ny) %in% as.numeric( colnames(z0) ) ) ) )
  
  if( length(wx) > 0 ){
    xx <- matrix( NA, length(wx), ncol(z0) )
    rownames(xx) <- wx
    z0 <- rbind( z0, xx )
    z0 <- z0[order( as.numeric( rownames(z0)) ),]
  }
  if( length(wy) > 0 ){
    xx <- matrix( NA, nrow(z0), length(wy) )
    colnames(xx) <- wy
    z0 <- cbind( z0, xx )
    z0 <- z0[, order( as.numeric( colnames(z0)) )]
  }
  
  rownames(z0) <- xseq[ as.numeric(rownames(z0)) ]
  colnames(z0) <- yseq[ as.numeric(colnames(z0)) ]
  
  kn <- 4
  
  wna <- which( is.na( z0 ), arr.ind = T )
  if( length(wna) > 0 ){
    wf <- which( is.finite( z0 ), arr.ind = T )
    tmp <- RANN::nn2( wf, wna, k = kn, radius = 2 )
    mm <- tmp[[1]]
    dd <- tmp[[2]] 
    dd[ dd > 3 ] <- NA
    wm <- which( rowSums( (dd*0 + 1), na.rm = T ) > 1 )
    
    fix <- wna[wm,]
    dd  <- dd[wm,]
    mm  <- mm[wm,]
    
    ss <- matrix( z0[ wf[mm,] ], ncol = kn )
    ss[ is.na(dd) ] <- 0
    
    z0[ fix ] <- rowSums( ss/dd^2, na.rm = T )/rowSums( dd^2, na.rm = T )
  }
  
  xseq <- xseq[-1] - .5*diff(xseq)[1] # midpoints
  yseq <- yseq[-1] - .5*diff(yseq)[1]
  
  
  if( LOG )z0 <- log10( z0 )
  
  z0[ !is.finite(z0) ] <- NA
  
  xseq <- as.numeric(rownames(z0))
  yseq <- as.numeric(colnames(z0))
  
  list( x = xseq, y = yseq, z = z0 )
}



grid2map <- function( xy, z, nx = 20, ny = 20, 
                      colRamp = rampWhite2Green,
                      xlim = NULL, ylim = NULL, zlim = NULL,
                      xaxt = 's', yaxt = 's', units = '', 
                      MAP = T, LEGEND = T, 
                      xleg = 'bottomright', 
                      yleg = c( ylim[1] - .2, ylim[1] + 1 ) ){
  
  # xy comes from plot file or expand.grid( x, y )
  # z is x by y matrix
  
  if( is.null( xlim ) ){
    xlim = range( xy[,1] ) + 1.2*c(-1, 1)
    ylim = range( xy[,2] ) + 1.2*c(-1, 1)
  }
  if( is.null( zlim ) )zlim <- range( z, na.rm = T ) + 1.2*c(-1, 1)
  
  wx <- which( xy[,1] >= xlim[1] & xy[,1] <= xlim[2] &
                 xy[,2] >= ylim[1] & xy[,2] <= ylim[2] )
  xy <- xy[wx,]
  z  <- z[wx]
  
  tgrid <- points2grid( x = xy[,1], y = xy[,2], z = z, nx = nx, ny = ny )
  z <- tgrid$z
  z[ !is.finite(z) ] <- NA
  x <- tgrid$x
  y <- tgrid$y
  
  surf2map( x, y, z, colRamp = colRamp,
            xlim = xlim, ylim = ylim, zlim = zlim,
            xaxt = xaxt, yaxt = yaxt, 
            MAP = MAP, LEGEND = LEGEND, 
            units = units, xleg = xleg, yleg = yleg )
}

surf2map <- function( x, y, z, colRamp = c('#a6611a','#dfc27d','#80cdc1','#018571'),
                      xlim = NULL, ylim = NULL, zlim = NULL,
                      xaxt = 's', yaxt = 's', 
                      MAP = T, LEGEND = T,
                      units = '', xleg = 'bottomright', 
                      yleg = c( ylim[1] - .2, ylim[1] + 1 ) ){
  
  if( is.null( xlim ) ){
    xlim = range( x ) + .1*c(-1, 1)
    ylim = range( y ) + .1*c(-1, 1)
  }
  wx <- which( x >= xlim[1] & x <= xlim[2] )
  x  <- x[ wx ]
  wy <- which( y >= ylim[1] & y <= ylim[2] )
  y  <- y[ wy ]
  z  <- z[wx, wy]
  
  if( is.null( zlim ) )zlim <- range( z, na.rm = T ) + 1.2*c(-1, 1)
  
  z[ z > zlim[2] ] <- zlim[2]
  
  image( x, y, z, xlim = xlim, ylim = ylim, col = colRamp, bty = 'n', 
         zlim = zlim, xlab = '', ylab = '', asp = 1, xaxt = xaxt, yaxt = xaxt,
         add = F )
  if( MAP ){
    maps::map('world', xlim = xlim, ylim = ylim, add = T, lwd = 2.5, col = 'white' )
    maps::map('world', xlim = xlim, ylim = ylim, add = T, lwd = 2, col = 'grey' )
  }
  
  zlim <- round( zlim, -2 )
  
  if( LEGEND )
    colorScaleLegend( xleg = xleg, yleg = yleg, 
                      xlim = xlim, ylim = ylim, 
                      zlim = zlim, ncol = 15, units = units, 
                      colRamp = colRamp ) # in mastifFunctions.R
}



poisReg <- function(n = 100, E = 1, nsim = 1000, slope = .5, PLOT = T, 
                    col = 'black'){
  
  # simulate data for intercept and slope with a
  # Poisson likelihood, log link, and effort E
  
  beta <- matrix( c(.1, slope), ncol = 1)
  x    <- cbind(1, rnorm(n) )
  mu   <- exp(x%*%beta)
  y    <- rpois( n, mu*E )
  
  # estimate beta, predict counts, note offset E for log link
  fit  <- glm( y ~ x[,2], "poisson", offset = rep(log(E), n) )
  pred <- predict.glm(fit, type = 'response', se.fit = T)      # contains error in lambda
  
  # combine parameter error with Poisson observations
  lambda <- rnorm(n*nsim, pred$fit, pred$se)
  lambda[ lambda < .001 ] <- .001
  py <- rpois( n*nsim, lambda )
  py <- matrix(py, n, nsim)
  pred$se <- apply( py, 1, sd )                                # lambda and sampling error
  
  if(PLOT){
    rmspe <- sqrt( mean( fit$residuals^2 ) )
    z <- jitter(y)
    plot(z, pred$fit, xlab = 'counts', ylab = 'predicted', 
         pch = 14, cex=.3, col = col)
    abline(0, 1, lty=2)
    abline( h = mean(y), lty = 2 )
    segments( z, pred$fit - pred$se.fit, 
              z, pred$fit + pred$se.fit, col = col )
    title( paste('rmspe =', round(rmspe, 2)) )
  }
  out <- cbind(x[,2], y, pred$fit, pred$se)
  colnames(out) <- c('x2', 'y', 'yhat','yse')
  invisible( list( fit = fit, out = out) )
}

monthAxis <- function( at = c(1, 4, 7, 10 ), cex.axis = 1.2, tick = T ){
  
  mnames <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec' )
  
  axis( 1, at = at, labels = mnames[at], cex.axis = cex.axis, tick = tick ) 
  axis( 1, at = c(1,12), labels = F )
}

.tnormMVNmatrix <- function(avec, muvec, smat, 
                            lo=matrix(-1000,nrow(muvec),ncol(muvec)), 
                            hi=matrix(1000,nrow(muvec),ncol(muvec)),
                            whichSample = c(1:nrow(smat)) ){
  
  #lo, hi must be same dimensions as muvec,avec
  
  if( !is.matrix( avec ) )avec <- matrix( avec, 1 )
  if( !is.matrix( muvec ) )muvec <- matrix( muvec, 1 )
  if( !is.matrix( lo ) )lo <- matrix( lo, 1 )
  if( !is.matrix( hi ) )hi <- matrix( hi, 1 )
  
  lo[lo < -1000] <- -1000
  hi[hi > 1000]  <- 1000
  
  if(max(whichSample) > length(muvec))
    stop('whichSample outside length(muvec)')
  
  r <- avec
  a <- trMVNmatrixRcpp(avec, muvec, smat, lo, hi, whichSample, 
                       idxALL = c(0:(nrow(smat)-1)) )  
  r[,whichSample] <- a[,whichSample]
  r
}

.getPlotLayout <- function( np, WIDE = TRUE ){
  
  # np - no. plots
  
  if( np == 1 )return( list( mfrow = c( 1, 1 ), left = 1, bottom = c( 1, 2 ) ) )
  if( np == 2 ){
    if( WIDE )return( list( mfrow = c( 1, 2 ), left = 1, bottom = c( 1, 2 ) ) )
    return( list( mfrow = c( 2, 1 ), left = c( 1, 2 ), bottom = 2 ) )
  }
  
  if( np == 3 ){
    if( WIDE )return( list( mfrow = c( 1, 3 ), left = 1, bottom = c( 1:3 ) ) )
    return( list( mfrow = c( 3, 1 ), left = 1:3, bottom = 3 ) )
  }
  if( np <= 4 )return( list( mfrow = c( 2, 2 ), left = c( 1, 3 ), bottom = c( 3:4 ) ) )
  if( np <= 6 ){
    if( WIDE )return( list( mfrow = c( 2, 3 ), left = c( 1, 4 ), bottom = c( 4:6 ) ) )
    return( list( mfrow = c( 3, 2 ), left = c( 1, 3, 5 ), bottom = 5:6 ) )
  }
  if( np <= 9 )return( list( mfrow = c( 3, 3 ), left = c( 1, 4, 7 ), bottom = c( 7:9 ) ) )
  if( np <= 12 ){
    if( WIDE )return( list( mfrow = c( 3, 4 ), left = c( 1, 5, 9 ), bottom = c( 9:12 ) ) )
    return( list( mfrow = c( 4, 3 ), left = c( 1, 4, 7, 10 ), bottom = 10:12 ) )
  }
  if( np <= 16 )return( list( mfrow = c( 4, 4 ), left = c( 1, 5, 9, 13 ), 
                              bottom = c( 13:16 ) ) )
  if( np <= 20 ){
    if( WIDE )return( list( mfrow = c( 4, 5 ), left = c( 1, 6, 11, 15 ), 
                            bottom = c( 15:20 ) ) )
    return( list( mfrow = c( 5, 4 ), left = c( 1, 5, 9, 13 ), bottom = 17:20 ) )
  }
  if( np <= 25 )return( list( mfrow = c( 5, 5 ), left = c( 1, 6, 11, 15, 20 ), 
                              bottom = c( 20:25 ) ) )
  if( np <= 25 ){
    if( WIDE )return( list( mfrow = c( 5, 6 ), left = c( 1, 6, 11, 15, 20, 25 ), 
                            bottom = c( 25:30 ) ) )
    return( list( mfrow = c( 6, 5 ), left = c( 1, 6, 11, 16, 21, 26 ), bottom = 26:30 ) )
  }
  if( np <= 36 ){
    return( list( mfrow = c( 6, 6 ), left = c( 1, 7, 13, 19, 25, 31 ), bottom = c( 31:36 ) ) )
  }
  return( list( mfrow = c( 7, 6 ), left = c( 1, 7, 13, 19, 25, 31, 37 ), bottom = c( 37:42 ) ) )
}


wideFormat2waterBalance <- function( prec1, pet1, years, ylim = range( c(prec1, pet1) ),
                                     xlab = '', ylab = c( 'P, PET, (mm)'), 
                                     yaxt = rep( 's', nrow(prec1) ) ){
  sites <- rownames(prec1)
  nsite <- length(sites)
  
  if( length( yaxt ) == 1 )yaxt <- rep( yaxt, nsite )
  
  ym <- columnSplit( colnames(prec1), '_' )
  yi <- ym[,1]
  mi <- ym[,2]
  
  ww <- which( yi %in% years )
  prec1 <- prec1[,ww]
  pet1  <- pet1[, ww]
  yi   <- yi[ww]
  mi   <- mi[ww]
  
  siteIndex <- rep( sites, ncol(prec1) )
  monIndex  <- rep(mi, each = nsite)
  siteMonth <- list( site = siteIndex, month = monIndex )
  pr <- round( tapply( as.vector ( prec1 ), siteMonth, mean ), 1 )[sites,]
  pe  <- round( tapply( as.vector ( pet1 ), siteMonth, mean ), 1 )[sites,]
  
  for( i in 1:nsite ){
    plot( NA, xlim = c(0, 13), ylim = ylim, bty = 'n', xaxt = 'n', las = 1,
          ylab = ylab, yaxt = yaxt[i] )
    monthAxis()
    
    shadeThreshold( 1:12, pr[i,], tmin = pe[i,], tmax = NULL, 
                    border = 'darkblue', col = '#9ecae1', LINES = T )
    shadeThreshold( 1:12, pe[i,], tmin = pr[i,], tmax = NULL, 
                    border = '#bd0026', col = '#fed976', LINES = T )
    title( sites[i], cex.main = 1.2, cex.main = 1.5)
  }
}

shadeInterval <- function( xvalues, loHi, col = 'grey', PLOT = TRUE, add = TRUE, 
                            xlab = ' ', ylab = ' ', xlim = NULL, ylim = NULL, 
                            LOG = FALSE, trans = .5 ){
  tmp <- NULL
  
  #draw shaded interval
  
  loHi <- as.matrix( loHi )
  tmp  <- smooth.na( xvalues, loHi )
  
  xvalues <- tmp[, 1]
  loHi    <- tmp[, -1]
  
  xbound <- c( xvalues, rev( xvalues ) )
  ybound <- c( loHi[, 1], rev( loHi[, 2] ) )
  if( is.null( ylim ) )ylim <- range( as.numeric( loHi ) )
  if( is.null( xlim ) )xlim <- range( xvalues )
  
  if( !add ){
    if( !LOG )plot( NULL, xlim = xlim, ylim = ylim, 
                    xlab = xlab, ylab = ylab, bty = 'n' )
    if( LOG )suppressWarnings( plot( NULL, xlim = xlim, ylim = ylim, 
                                     xlab = xlab, ylab = ylab, log = 'y', bty = 'n' ) )
  }
  
  
  if( PLOT )polygon( xbound, ybound, border = NA, col = .getColor( col, trans ) )
  
  invisible( cbind( xbound, ybound ) )
}


smooth.na <- function(x,y){   
  
  #remove missing values
  #x is the index
  #y is a matrix with rows indexed by x
  
  if(!is.matrix(y))y <- matrix(y,ncol=1)
  
  wy <- which(!is.finite(y),arr.ind =T)
  if(length(wy) == 0)return(cbind(x,y))
  wy <- unique(wy[,1])
  ynew <- y[-wy,]
  xnew <- x[-wy]
  
  return(cbind(xnew,ynew))
}

.getColor <- function( col, trans ){
  
  # trans - transparency fraction [ 0, 1]
  
  tmp <- col2rgb( col )
  rgb( tmp[ 1, ], tmp[ 2, ], tmp[ 3, ], maxColorValue = 255, 
       alpha = 255*trans, names = paste( col, trans, sep = '_' ) )
}


shadeThreshold <- function( x, y, tmin = NULL, tmax = NULL, ylim = range( c(y, tmin, tmax) ),
                            border = 'brown', col = '#fed976', xaxt = 's', yaxt = 's', 
                            add = T, LINES = F ){
  
  # LINES draws vertical lines at data points
  
  if( length(y) == 1 ){
    if( !is.null( tmin ) )y <- rep( y, length(tmin) )
    if( !is.null( tmax ) )y <- rep( y, length(tmax) )
  }
  if( !is.null( tmin ) )
    if( length(tmin) == 1) tmin <- rep( tmin, length(y) )
  if( !is.null( tmax ) )
    if( length(tmax) == 1) tmax <- rep( tmax, length(y) )
  
  if( !add ){
    plot( x, y, type = 'l', xaxt = 'n', xlab = '', 
          xaxt = xaxt, yaxt = yaxt,
          ylab = '', col = 'white', lwd = 6, ylim = ylim,
          bty = 'n' )
  }
  
  if( !is.null( tmin ) ){   # shade above tmin
    pts <- intersectLines( x, y1 = y, y2 = tmin )
    
    if( length(pts$poly1) > 0 ){
      for( j in 1:length( pts$poly1 ) )
        polygon(pts$poly1[[j]][,1], pts$poly1[[j]][,2], border = border, col = col)
    }
  }
  if( !is.null( tmax ) ){  # shade below tmax
    
    pts <- intersectLines( x, y1 = y, y2 = tmax )
    
    if( length(pts$poly2) > 0 ){
      for( j in 1:length( pts$poly2 ) )
        polygon(pts$poly2[[j]][,1], pts$poly2[[j]][,2], border = border, col = col)
    }
  }
  if( LINES ){
    
    if( is.null(tmax) ){
      ww <- which( y > tmin )
      if( length(ww) > 0 )segments( x[ww], tmin[ww], x[ww], y[ww], col = border )
    }
    if( is.null(tmin) ){
      ww <- which( y < tmax )
      if( length(ww) > 0 )segments( x[ww], y[ww], x[ww], tmax[ww], col = border )
    }
  }
}


points2map <- function( x, inputList, verbose = F, traitPath = '../../' ){
  
  require( maps )
  require(mapdata)
  require(maptools)
  require(ggplot2)
  require(rgeos)
  
  if( length(x) > 1 )stop('x longer than 1')
  
  # ipath = '../inventories/'
  
  # on log10 scale:
  logVars <- c("fecGmMu", "fecGmSd", "disp10gmMu", "disp10gmSd",
               "fecMu","fecSd", "disp10seedMu", "disp10seedSd", "ISP","BA", "SIZE" )
  climVars <- c("temp","prec","def","pet")
  soilVars <- c("ph","cec30","slope")
  
  species <- genus <- xlim <- ylim <- zlim <- years <- xlegend <- ylegend <- NULL
  continent <- 'northAmerica'
  sigx  <- sigy   <- 1
  CV    <- RASTER <- F
  ANOMALY <- PLOT  <- SMOOTH <- LEGEND <- T
  CIRCLES <- F
  cex     <- .8
  scaleSymbol <- 2
  EOF <- 0
  outputDir  <- NULL
  outputType <- 'png'
  colorRamp <- colGreens
  
  for(k in 1:length(inputList) )assign( names(inputList)[k], inputList[[k]] )
  
  if( !is.null(years) )years  <- years[years %in% 2008:2022]
  
  kvars <- c("fecMu","fecSd", "fecGmMu","fecGmSd","disp10seedMu",  # gm to kg
             "disp10seedSd", "disp10gmMu", "disp10gmSd") 
  ivars <- 'ISP'        # per tree BA variables
  svars <- 'SIZE'       # community-weighted mean
  avars <- c( kvars, "SIZE", "BA", "ISP", "mastifData", climVars, soilVars)  
  tcols <- c("plot","lon","lat","species", "year","diam","shade","wtPerHa","ecoReg",
             "tempSite","deficitSite","cec30","slope","aspect", "betaYr")
  if( !x %in% avars ){
    stop( paste( 'x must be one of these:', paste0( avars, collapse = ', ' ) ) )
  }
  
  if( is.null(years) ){
    ANOMALY <- F
    if( EOF > 0 )stop( 'EOF only possible if there are multiple years' )
    nyr <- length(years)
    if( EOF > nyr )EOF <- nyr
  }
  if( ANOMALY | EOF )colorRamp <- colBlu2Red
  
  # g/ha to kg/ha
  if( x %in% c("fecGmMu", "fecGmSd", "disp10gmMu", "disp10gmSd") )units <- 'log10 kg/ha'
  if( x %in% c("fecMu","fecSd", "disp10seedMu", "disp10seedSd" ) )units <- 'log10 per ha'
  if( x %in% c("ba", "BA") )units <- 'log10 m2/ha'
  if( x == 'ISP' ) units <- 'log10 g/m2'  # per m2 basal area
  if( x == 'SIZE' )units <- 'log10 g'
  if( CV )         units <- 'log10 CV'
  if( x == 'mastifData' )units = 'tree-years'
  
  if( !is.null(outputDir) ){
    cat('\nPlots will be written to outputDir')
    cat('\nSend plots to console by removing outputDir from inputList')
    if( !dir.exists(outputDir) )dir.create( outputDir )
  }
  
  sfile <- paste( sumPlotPath, continent, '/summary/', x, '.Rdata', sep = '')
  plotFile <- .replaceString( sfile, x, 'plotData' )
  
  # fitted data
  if( x == 'mastifData' ){
    
    sfile    <- paste( sumPlotPath, continent, '/summary/mastifData.rdata', sep = '')
    plotFile <- paste( dataPlotPath, 'mastPlotsAll.csv', sep = '')
    
    EOF   <- 0
    tmp   <- getSummaryColumns( x, sfile, continent, xlim, ylim, plotFile, 
                                species, genus, years, EOF, CV, ANOMALY, verbose ) 
    summary <- tmp$summary
    taxa    <- tmp$taxa
    xlim    <- tmp$xlim
    ylim    <- tmp$ylim
  }
  
  # summary files are plots by species, years averaged or summed or individual files
  
  fe <- file.exists( sfile )
  
  if( !fe & !x %in% c('mastifData', climVars, soilVars) ){
    
    if( length( list.files( ipath ) ) == 0 )stop( 'files in ipath not available' )
    
  }else if( x %in% c(climVars, soilVars) ){
    
    load( plotFile, verbose = T )
    
    if( x %in% climVars ){
      
      clfile = paste( '../inventories/', continent, '/', x, '.rdata', sep = '' )
      load( clfile, verbose = T )
      assign( 'tmp', get( x ) )
      
      tmp <- tmp[rownames(plotData), ]
      
      clim <- wide2Anom( vname = x, wideForm = tmp, normYr = 2000:2020 )
      summary <- cbind( plotData[, c('lon','lat')], clim$norm )
      if( !is.null(years) ){
        
        if( ANOMALY ){
          summary   <- clim$anom
          summary <- summary[, colnames(summary) %in% as.character(years)]
          summary <- cbind( plotData[, c('lon','lat')], summary )
        }else{
          summary[,3] <- clim$annual
        }
      }
      rm( clim )
    }
    
    if( x %in% soilVars ){
      summary <- plotData[, c('lon','lat', x)] 
    }
    
  }else{
    
    if( x != 'mastifData' ){
      
      tmp <- getSummaryColumns( x, sfile, continent, xlim, ylim, plotFile, 
                                species = species, genus = genus, 
                                years = years, EOF = EOF, CV, ANOMALY, verbose = verbose )
      summary <- tmp$summary
      varEOF  <- tmp$varEOF
      taxa    <- tmp$taxa
      xlim    <- tmp$xlim
      ylim    <- tmp$ylim
    }
    
    #  z <- log10( summary[,3] + .01 )
    #  qz <- seq( -2, 3, by = 1 )
    #  colm <- findInterval( z, qz, all.inside = T )
    #  symbols( summary[,1], summary[,2], 
    #           rectangles = cbind( rep(pointsWide, nrow(summary)), rep(pointsHigh, nrow(summary))), 
    #           inches = F, fg = colf[colm], bg = colf[colm], add = F )
    
  }
  
  if( !PLOT )return( summary )
  
  
  xseq <- floor( min(summary[,1]) ):ceiling( max(summary[,1]) )
  yseq <- floor( min(summary[,2]) ):ceiling( max(summary[,2]) )
  
  #dd <- summary
  
  if( !is.null(years) ){
    
    tmp <- smoothGrid( x = summary[,1], y = summary[,2], z = summary[,3], 
                       sigx, sigy, xlim = xlim, ylim = ylim,
                       zmax = NULL, SMOOTH = T)
    z    <- tmp$z; 
    xseq <- tmp$xseq; yseq <- tmp$yseq
    dd   <- tmp$dd
    nyr  <- length(years)
  }else{
    if( x != 'mastifData' ){
      
      tmp <- smoothGrid( x = summary[,1], y = summary[,2], z = summary[,3], 
                         sigx, sigy, xlim = xlim, ylim = ylim,
                         zmax = NULL, SMOOTH = T)
      z  <- tmp$z; xlim = tmp$xlim; ylim = tmp$ylim; xseq = tmp$xseq; yseq = tmp$yseq
      dd <- tmp$dd
    }
  }
  
  #dd <- dd[ dd[,3] != 0, ]
  
  pointsWide <- diff(xlim)/length(xseq)
  pointsHigh <- diff(ylim)/length(yseq)
  
  
  if( x == 'mastifData' & !CIRCLES ){
    
    pointsWide <- pointsHigh <- 1
    qz <- c(0, 10^(0:6) )
    summary[ ,-c(1:2)] <- findInterval( summary[,-c(1:2)],  qz , all.inside = T ) - 1      
    qz <- qz[-1]
    qz <- qz[ qz <= 10^max( summary[ ,-c(1:2)] ) ]
    colf <- colorRampPalette( colReds[-1] )( length(qz) )
  }
  
  tt <- getLims( continent )
  mdat <- tt$world
  sdat <- tt$state
  
  if( x != 'mastifData' ){
    ylist <- NULL
    if( !is.null(years) ){
      ylist <- summary
      z     <- NULL
    }else{
      z <- z[z > 0]
    }
    #   yz   <- ylist
    LOGZ <- F
    if( x %in% logVars )LOGZ = T
    
    
    tmp  <- getZscale( z, zlim = zlim, ylist = ylist, cramp = colorRamp, LOGZ = LOGZ )
    qz   <- tmp$qz
    colf <- tmp$cols
    ncol <- length(qz)
  }
  
  nmap <- 1
  if( !is.null(years) )nmap <- nyr
  if( EOF > 0 )nmap <- EOF
  
  col <- 1     # for mastifData
  if( x != 'mastifData' & is.null(years) ){
    
    d3 <- dd[,3]
    if( x %in% logVars)d3 <- log10(d3)
    
    col <- findInterval( d3, qz, all.inside = T)  # log10 scale
  }
  if( x == 'mastifData' ){
    if( nmap == 1 ){
      maxobs <- max(summary[,'z'], na.rm = T)
    }else{
      maxobs <- max( summary[,as.character(years)], na.rm = T )
    }
  }
  
  cnames <- taxa
  if( !is.null(years) )cnames <- colnames(summary)[-c(1,2)] 
  
  for(m in 1:nmap){
    
    tname <- cnames[m]
    if( EOF > 0)tname <- paste( colnames(varEOF)[m], round( varEOF['varFrac',m], 3 ) ,sep = ', ' )
    
    if( !is.null( outputDir ) ){
      oname <- .replaceString( tname, ' ', '' )
      oname <- .replaceString( oname, ',', '' )
      mfile <- paste( outputDir, '/', x, '_', oname, '_',  continent, '.', outputType, sep = '' )
      if( outputType == 'jpeg')jpeg( mfile, width = 350, height = 350)
      if( outputType == 'png')png( mfile, width = 350, height = 350)
      if( outputType == 'pdf')pdf( mfile, width = 350, height = 350)
    }
    
    colm <- col
    if( x == 'mastifData' ){
      dd   <- points2degree( summary[,'lon'], summary[,'lat'], summary[,3] )
      colm <- findInterval( dd[,3], qz, all.inside = T)
    }
    
    if( !is.null(years) ){
      
      dd <- summary[, c('lon', 'lat', cnames[m]) ]
      
      if( x == 'mastifData' ){
        dd <- points2degree( summary[,'lon'], summary[,'lat'], summary[,cnames[m]] )
      }
      colm <- findInterval( dd[,3], qz, all.inside = T)
    }else{
      if( x != 'mastifData' )colm <- findInterval( log10(dd[,3]), qz, all.inside = T)
    }
    
    if( x == 'mastifData' ){
      
      map( 'world', xlim = xlim, ylim = ylim, fg = '#7bccc4', col = '#ccebc5', 
           border = '#7bccc4', fill = T )
      map( sdat, col = '#7bccc4', add=T )
      
      if( CIRCLES ){
        zm <- scaleSymbol*2*(dd[,3]/maxobs/pi)^.5
        zm[ zm < scaleSymbol*.3 ] <- scaleSymbol*.3
        points( dd[,1], dd[,2], pch = 16, cex = scaleSymbol*.2 + zm, col = 'white' )
        points( dd[,1], dd[,2], pch = 16, cex = zm, col = .getColor( '#67000d', .3) )
      }else{
        symbols( dd[,1], dd[,2], 
                 rectangles = cbind( 1.8*rep(pointsWide, nrow(dd)), 1.8*rep(pointsHigh, nrow(dd))), 
                 inches = F, fg = 'white', bg = 'white', add = T )
        symbols( dd[,1], dd[,2], 
                 rectangles = cbind( rep(pointsWide, nrow(dd)), rep(pointsHigh, nrow(dd))), 
                 inches = F, fg = colf[colm], bg = colf[colm], add = T )
      }
      
    }else{
      
      map( 'world', xlim = xlim, ylim = ylim, bg = 'white', col = '#f7fcf0', 
           border = '#a8ddb5', fill = T )
      map( sdat, col = '#a8ddb5', add=T )
      symbols( dd[,1], dd[,2], 
               rectangles = cbind( rep(pointsWide, nrow(dd)), rep(pointsHigh, nrow(dd))), 
               inches = F, fg = colf[colm], bg = colf[colm], add = T )
      map( 'world', fill = F, col = '#a8ddb5', add=T )
    }
    if( !is.null(years) )title( tname )
    
    if( nmap > 1 & m < nmap & !is.null( outputDir ) )dev.off()
  }
  
  if( length(species) == 0 )species <- NULL
  
  f1 <- ' '
  if( !is.null(genus) & length(genus) == 1 )f1 <- upperFirstLetter( genus )
  if( !is.null(species) ){
    tpath <- paste( traitPath, 'traitsByGroup/', sep='') 
    gen <- getTraits( specNames = species, traitName =  'genus', path2Traits = tpath ) 
    spe <- getTraits( specNames = species, traitName =  'specEpith', path2Traits = tpath )
    f1  <- paste( gen, spe )
  }
  
  if( nchar(f1)[1] > 0 & is.null(years) ){
    title( main=f1, font.main = 3, font.sub = 4 )
  }
  
  if( LEGEND ){
    
    rr <- 1
    ytick <- round(qz, rr)
    while( length(which(duplicated(ytick))) > 0 ){
      rr <- rr + 1
      ytick <- round(qz, rr)
    }
    
    if(is.null(xlegend)){
      dx <- .02*diff(par('xaxp')[1:2])
      xlegend <- par('xaxp')[1] + dx*c(1.5, 3)
    }
    
    if(is.null(ylegend)){
      dy <- .5*diff(par('yaxp')[1:2])
      ylegend <- par('yaxp')[1] + dy*c(.35, .85)
    }
    
    dx <- diff(xlegend)
    dy <- diff(ylegend)
    mx <- diff(xlim)
    my <- diff(ylim)
    
    xtext <- xlegend[1] - .15*mx
    x0    <- min( xtext, xlegend[1] - 8 )
    
    rect( x0, ylegend[1] - 1, xlegend[2] + 1, ylegend[2] + 1, 
          col = 'white', border = NA )
    
    par( xpd = T )
    colorLegend(xlegend, ylegend, ytick=ytick,
                scale = qz, cols=colf, labside='right', cex = cex)
    xl <- xtext
    if( x == 'mastifData' )xl <- xl - log10( max(qz) ) - 1
    text( xl, .01*dy + mean(ylegend), units, pos=3, 
          srt=90, cex = cex)
    par( xpd = F )
  }
  
  if( !is.null( outputDir ) ){
    dev.off()
  }
  
  out <- dd
  if( !is.null(years) )out <- summary
  invisible( out )
}


intersectLines <- function (x, y1, y2){
  
  # intersection points and polygons for y1 > y2 and y2 > y1
  
  DATE <- F
  
  if( inherits(x[1], 'Date') ){
    # assumes YYYY-MM-DD
    DATE <- T
    date <- x
    xx <- columnSplit( as.character( x ), '-' )
    x  <- as.numeric(xx[,1]) + as.POSIXlt(x)$yday/365
    #   origin <- as.POSIXlt( as.Date(paste( min( xx[,1] ), '-01-01', sep = '' )) 
  }
  
  tiny <- 1e-12
  y1 <- y1 + rnorm( length(y1), 0, tiny )
  
  n <- length(x)
  above <- y1 > y2
  intersectPts <- which(diff(above) != 0) 
  
  y1.diff <- y1[intersectPts+1] - y1[intersectPts]
  y2.diff <- y2[intersectPts+1] - y2[intersectPts]
  x.diff  <- x[intersectPts+1] - x[intersectPts]
  
  slope1 <- y1.diff/x.diff
  slope2 <- y2.diff/x.diff
  intercept1 <- y1[intersectPts] - slope1*x[intersectPts]
  intercept2 <- y2[intersectPts] - slope2*x[intersectPts]
  xPts <- ifelse(slope1 == slope2, NA, 
                 (intercept2-intercept1)/(slope1-slope2))
  yPts <- ifelse(slope1 == slope2, NA,
                 slope1*xPts+intercept1)
  
  jointPts <- which(y1 == y2)
  xPts <- c(xPts, x[jointPts])
  yPts <- c(yPts, y1[jointPts])
  
  ipoints <- cbind( xPts,  yPts)
  
  pt1  <- rbind( c(x[1], y1[1]), ipoints, c(x[n], y1[n]) )
  pt2  <- rbind( c(x[1], y2[1]), ipoints, c(x[n], y2[n]) )
  
  pt2 <- pt2[ !duplicated(pt1[,1]), ]
  pt1 <- pt1[ !duplicated(pt1[,1]), ]
  
  p1 <- p2 <- numeric(0)
  
  for( k in 2:nrow(pt1) ){
    
    w <- which( x >= pt1[k-1,1] & x <= pt1[k,1] )
    if( length(w) == 0 )next
    if( above[w[1]] ){
      pk <- rbind( pt1[k-1,], cbind(x[w], y1[w]), pt1[k,] )
      pb <- rbind( pt2[k-1,], cbind(x[w], y2[w]), pt2[k,] ) #bottom
      pk <- rbind( pk, pb[ nrow(pb):1, ] )
      pk <- pk[ !duplicated(pk[,1:2]), ]
      pk <- rbind( pk, pk[1,] )
      p1 <- append( p1, list( pk ) )
    }else{
      pk <- rbind( pt2[k-1,], cbind(x[w], y2[w]), pt2[k,] ) #top
      pb <- rbind( pt1[k-1,], cbind(x[w], y1[w]), pt1[k,] ) #bottom
      pk <- rbind( pk, pb[ nrow(pb):1, ] )
      pk <- pk[ !duplicated(pk[,1:2]), ]
      pk <- rbind( pk, pk[1,] )
      p2 <- append( p2, list( pk ) )
    }
  }
  if( length(p1) > 0 ){
    for( k in 1:length(p1)){  # all above y2 
      p1[[k]][ abs(p1[[k]][,2]) < 1e-10, 2] <- 0
      pk <- p1[[k]]
      pn   <- RANN::nn2( x, pk[,1], k = 1 )[[1]]
      ymin <- y2[ pn ]
      whi  <- which( pk[, 2] < ymin )
      if( length( whi ) == 0 )next
      pk[ whi, 2] <- ymin[ whi ]
      p1[[k]] <- pk
    }
  }
  if( length(p2) > 0 ){
    for( k in 1:length(p2)){  # all below y2 
      p2[[k]][ abs(p2[[k]][,2]) < 1e-10, 2] <- 0
      pk   <- p2[[k]]
      pn   <- RANN::nn2( x, pk[,1], k = 1 )[[1]]
      ymax <- y2[ pn ]
      whi  <- which( pk[, 2] > ymax )
      if( length( whi ) == 0 )next
      pk[ whi, 2] <- ymax[ whi ]
      p2[[k]] <- pk
    }
  }
  if( length( p1 ) > 1 ){
    p1 <- p1[ which( !sapply( p1, var )[3,] == 0 ) ]
  }else{
    if( length(p1) == 1)p1 <- p1[ which( !sapply( p1, var )[4] == 0 ) ]
  }
  if( length( p1 ) > 1 ){
    p2 <- p2[ which( !sapply( p2, var )[3,] == 0 ) ]
  }else{
    if( length(p2) == 1)p2 <- p2[ which( !sapply( p2, var )[4] == 0 ) ]
  }
  
  list(ipoints = ipoints, poly1 = p1, poly2 = p2) 
}

extractPlot <- function( inputs, plot, taxa, plotData ){
  
  xytree <- NULL
  inputs$treeData <- inputs$treeData[ inputs$treeData$plot == plot, ]
  specNames <- inputs$specNames <- sort( unique( inputs$treeData$species ) )
  
  inputs$priorTable <- inputs$priorTable[ drop=F, specNames, ]
  
  if( !is.null( inputs$xytree ) ){
    
    g1 <- stri_locate_first_regex(specNames, "[A-Z]")
    gen <- substr( specNames, 1, g1 - 1 )
    unkn <- unique( paste( gen, 'UNKN', sep = '' ) )
    
    inputs$xytree   <- inputs$xytree[ inputs$xytree$plot == plot, ]
    inputs$xytrap   <- inputs$xytrap[ inputs$xytrap$plot  == plot, ]
    inputs$seedData <- inputs$seedData[ inputs$seedData$plot == plot, ]
    
    dcols <- colnames( inputs$seedData )
    seeds <- c( specNames, unkn )
    scols <- colnames(inputs$seedData)
    scols <- scols[ scols %in% seeds ]
    
    inputs$seedData <- inputs$seedData[,c('plot','trap','year','area','active', scols ) ]
    inputs$seedNames <- scols
  }
  taxa <- taxa[specNames, ]
  plotData <- plotData[ drop = F, match( plot, plotData$plot), ]
  list( inputs = inputs, plotData = plotData, taxa = taxa )
}

extractGenus <- function( inputs, taxa, plotData, genus = NULL ){
  
  xytree <- NULL
  
  taxa$genus <- tolower( taxa$genus )
  if( is.null( genus ) )genus <- sort( unique( taxa$genus ) )
  genus      <- tolower( genus )
  specNames  <- rownames( taxa )[taxa$genus %in% genus ]
  tplots     <- sort( unique( inputs$treeData$plot ) )
  plotData   <- plotData[ plotData$plot %in% tplots, ]
  
  inputs$priorTable <- inputs$priorTable[ drop=F, specNames, ]
  inputs$treeData   <- inputs$treeData[ inputs$treeData$species %in% specNames, ]
  if( !is.null( inputs$xytree ) ){
    inputs$xytree   <- inputs$xytree[ inputs$xytree$plot %in% tplots, ]
    inputs$xytrap   <- inputs$xytrap[ inputs$xytrap$plot %in% tplots, ]
    inputs$seedData <- inputs$seedData[ inputs$seedData$plot %in% tplots, ]
    
    dcols <- colnames( inputs$seedData )
    scols <- character(0)
    seeds <- specNames
    for( k in 1:length(genus) ){
      seeds <- c( seeds, paste( genus, 'UNKN', sep = '' ) )
    }
    for( k in 1:length(seeds)){
      scols <- c( scols, dcols[ startsWith( dcols, seeds[k] ) ] )
    }
    inputs$seedData <- inputs$seedData[,c('plot','trap','year','area','active', scols ) ]
    inputs$seedNames <- scols
  }
  inputs$specNames <- sort( unique( inputs$treeData$species ) )
  taxa <- taxa[specNames, ]
  list( inputs = inputs, plotData = plotData, taxa = taxa )
}



row2Mat <- function(vec){
  
  if(is.matrix(vec))return(vec)
  vn  <- names(vec)
  vec <- matrix(vec,1)
  colnames(vec) <- vn
  vec
}

myrmultinom <- function(size, p, ASVECTOR=F){  
  
  # n multinomial r.v. for a n by ncol(p) matrix of probs
  # each row of p is a probability vector
  # size is one integer or a length-n vector of integers
  # if ASVECTOR = T all size == 1, returns a vector of columns, otherwise a matrix
  
  p <- row2Mat(p)
  
  n     <- nrow(p)
  J     <- ncol(p)
  
  if(length(size) == 1)size <- rep(size,n)
  
  jord  <- sample(J,J)    #randomize order
  
  rs <- rowSums(p)
  ws <- which(rs != 1)
  if(length(ws) > 0){
    p[ws,] <- p[ws,]/rs[ws]
  }
  
  p <- row2Mat(p[,jord])
  
  
  sizej <- size
  sumj  <- rep(0,n)
  dpj   <- rep(1,n)
  pj    <- p
  wj    <- c(1:n)
  
  if(ASVECTOR){        #  only if all size == 1
    
    yy <- size*0
    
    for(j in 1:(J-1)){
      a     <- round(pj[wj,1],10)
      tmp  <- rbinom(length(wj),sizej[wj],a)
      yy[wj[tmp == 1]] <- j
      sumj[wj]  <- sumj[wj] + tmp
      sizej <- size - sumj                       # no. remaining to choose
      dpj   <- dpj - p[,j]                       # Pr for remainder
      pj    <- matrix(p[,c((j+1):J)]/dpj,nrow(p))
      wj    <- which(sumj < size,arr.ind=T) 
    }
    
    yy[yy == 0] <- J
    
    return(yy)
  }
  
  yy  <- matrix(0,n,J)
  
  for(j in 1:(J-1)){
    a     <- round(pj[wj,1],10)
    yy[wj,j] <- rbinom(length(wj),sizej[wj],a)
    sumj  <- sumj + yy[,j]
    sizej <- size - sumj                       # no. remaining to choose
    dpj   <- dpj - p[,j]                       # Pr for remainder
    pj    <- matrix(p[,c((j+1):J)]/dpj,nrow(p))
    wj    <- which(sumj < size,arr.ind=T) 
  }
  
  if(n == 1)yy[,J] <- size - sum(yy)
  if(n > 1) yy[,J] <- size - rowSums(yy)
  
  yy[,jord] <- yy
  yy
  
}

obsGrid <- function( treeData, plotData, colRamp = c('#fdbb84','#fc8d59','#7f0000') ){
  
  # location of each observation
  mm <- match( inputs$treeData$plot, plotData$plot )
  ll <- plotData[mm, c('lon','lat') ]
  
  # aggregate to 1 degree grid
  loc  <- round( ll, 1 )
  ltab <- table( loc[,1], loc[,2] ) # no. observation at each grid point
  wtab <- which( ltab > 0, arr.ind = T )
  
  nobs <- ltab[ wtab ]
  xy   <- cbind( rownames( ltab )[wtab[,1]], colnames( ltab )[wtab[,2]] )
  effort <- data.frame( lon = as.numeric( xy[,1] ),
                        lat = as.numeric( xy[,2] ),
                        nobs )
  
  # assign colorRamp, low to high
  zlevs <- 1:10
  icol  <- findInterval( log10(1 + effort$nobs), zlevs ) + 1
  zlevs <- zlevs[ zlevs <= max(icol) ]
  cvals <- colorRampPalette( colRamp )( max( icol ) )
  effort$cols  <- cvals[ icol ]
  effort
}

upperFirstLetter <- function( xx ){
  
  # FIRSTONLY - only first letter of first word when a string has multiple words
  
  f <- toupper( substring( xx, 1, 1 ) )
  l <- substring( xx, 2 )
  paste( f, l, sep='' )
}

cor2cov <- function(sigvec, cormat){ 
  
  #correlation matrix and variance vector to covariance
  
  n <- nrow(cormat)
  d <- length(sigvec)
  if(d == 1)sigvec <- rep( sigvec, n )
  
  s <- matrix(sigvec,n,n)
  cormat*sqrt(s*t(s))
}

cov2cor <- function( covmat ){  
  
  # covariance matrix to correlation matrix
  # if covInv provided, return inverse correlation matrix
  
  d    <- nrow(covmat)
  di   <- diag(covmat)
  s    <- matrix(di,d,d)
  covmat/sqrt(s*t(s))
}

cov2dist <- function(sigma){ #distance induced by covariance
  
  n <- nrow(sigma)
  matrix(diag(sigma),n,n) + matrix(diag(sigma),n,n,byrow=T) - 2*sigma
}


predVsObsPlot <- function( obs, pred, col = 'black', 
                           ybins = NULL, ADD = F ){
  
  if( is.null( ybins ) )ybins <- unique( quantile( obs, seq(0, 1, by = .02 ) ) )
  yi    <- findInterval( obs, ybins, all.inside = T )
  q     <- tapply( pred, yi, quantile, c(.5, .156, .841, .025, .975 ) )
  yb    <- matrix( unlist(q), ncol = 5, byrow = T )
  mids  <- ybins[-length(ybins)] + diff(ybins)/2
  if( !ADD ){
  plot( mids, yb[,1], ylim = range( yb ), bty = 'n',
        xlab = 'Observed', ylab = 'Predicted' )
  }else{
    points( mids, yb[,1] )
  }
  for( j in 1:length(mids)){
    rect( mids[j] - .1, yb[j,2], mids[j] + .1, yb[j,3],
          col = .getColor( col, .5 ), border = col )
    lines( mids[ c(j,j) ], yb[j,4:5], col = col )
    lines( mids[j] + .12*c(-1,1), yb[j,c(1,1)], col = col )
  }
  abline(0, 1, lty = 2, lwd = 2, col = 'grey' )
}

colorScaleLegend <- function( xleg, yleg = NULL, zlim, 
                              zlimLabs = zlim,
                              xlim = NULL, ylim = NULL, ncol = 10, units = '',
                              colRamp ){
  
  xchar <- max( strwidth( zlimLabs ) )
  
  xusr <- par('usr')[1:2]
  yusr <- par('usr')[3:4]
  
  if( is.null( xlim ) )xlim <- xusr
  if( is.null( ylim ) )ylim <- yusr
  
  dx    <- diff( xlim )
  dy    <- diff( ylim )
  xwide <- dx/25
  yhigh <- dy/5
  
  
  if( is.character( xleg ) ){
    pos <- c('bottomleft', 'bottomright','topleft','topright')
    if( !xleg %in% pos )
      stop( paste( 'legend position must be:', paste0( pos, collapse = ', ' ) ) )
    
    pos <- xleg
    
    xleg  <- xlim[1] + c(-2*xwide, -xwide ) # bottomleft
    yleg  <- ylim[1] + c(0, yhigh )
    tpos  <- 4                              # right of scale
    xtext <- xleg[2] + dx/150
    
    if( pos == 'topleft' )yleg <- ylim[2] + c(-yhigh, 0 )
    
    if( xusr[1] < xlim[1] ){
      xleg  <- xleg - xwide
      xtext <- xtext - xwide
    }
    
    if( pos == 'bottomright' ){
      xleg <- xlim[2] + c(0, xwide )
      yleg <- ylim[1] + c(0, yhigh )
      tpos <- 2                     # left of scale
      xtext <- xlim[2] - dx/150
      if( xusr[2] > xlim[2] ){
        xleg <- xleg + xwide
        xtext <- xtext + xwide
      }
    }
    if( pos == 'topright' ){
      xleg <- xlim[2] + c(0, xwide )
      yleg <- ylim[2] + c( -yhigh, 0 )
      tpos <- 2                     # left of scale
      xtext <- xleg[1] - dx/150
      if( xusr[2] > xlim[2] ){
        xleg <- xleg + xwide
        xtext <- xtext + xwide
      }
    }
  }else{
    
    if( xleg[1] < ( xlim[1] + diff(xlim)/2 ) ){ # left side
      tpos  <- 4                                # right of scale
      xtext <- xleg[2] + dx/150
    }else{
      tpos <- 2                     # left of scale
      xtext <- xlim[2] - dx/150
    }
  }
  
  xbox <- xleg
  ybox <- yleg + dy/20*c(-1, 1)
  if( nchar( units ) > 0 )ybox[2] <- ybox[2] + dy/20
  if( tpos == 2 )xbox[1] <- xbox[1] - xchar - dx/10
  if( tpos == 4 )xbox[2] <- xbox[2] + xchar + dx/10
  
  lcol <- colorRampPalette(colRamp)( ncol )
  zseq <- seq( yleg[1], yleg[2], length = ncol + 1 )
  
  par( xpd = T )
  rect( xbox[1], ybox[1], xbox[2], ybox[2], col = 'white', border = NA )
  for( i in 1:ncol ){
    yi <- zseq[i:(i+1)]
    rect(xleg[1], yi[1], xleg[2], yi[2], col = lcol[i], border = NA )
  }
  
  textColor <- colorRampPalette(colRamp)( 2 )
  text( rep( xtext, 2 ), yleg, zlimLabs, pos = tpos, col = textColor, 
        cex = .9)
  
  if( nchar( units ) > 0 )text( mean( xleg ), yleg[2] + diff(yleg)/5, units,
                                cex = .9 )
  
  par( xpd = F )
}



imageMap <- function( lonLat, z, colRamp = NULL, 
                      xlim = range(lonLat[,1]), ylim = range(lonLat[,2]),
                      zlim = range( z, na.rm = T), ncol = 20, xaxt = 's', yaxt = 's', 
                      xlab = '', ylab = '', add = F ){
  
  # lonLat has rownames Elon_Nlat, e.g., 'E30.5_N-22.067797'
  # z is a climate vector with names(z) == rownames(lonLat)
  
  if( is.null( colRamp ) )
    colRamp <- c('#8c510a','#bf812d','#dfc27d','#80cdc1','#35978f','#01665e')
  
  cfun  <- colorRampPalette( colRamp )
  cfn   <- cfun(ncol)
  
  ww <- which( lonLat[,1] >= xlim[1] & lonLat[,1] <= xlim[2] &
                 lonLat[,2] >= ylim[1] & lonLat[,2] <= ylim[2] )
  lonLat <- lonLat[ww,]
  z <- z[ww]
  
  x <- sort( unique( lonLat[,1] ) )
  y <- sort( unique( lonLat[,2] ) )
  names(x) <- paste( 'E', x, sep = '' )
  names(y) <- paste( 'N', y, sep = '' )
  zmat <- matrix( NA, length(x), length(y), dimnames = list( names(x), names(y) ) )
  zmat[ columnSplit( names(z), '_' ) ] <- z
  
  image( x, y, zmat, xlim = xlim, ylim = ylim, zlim = zlim, xlab = xlab, ylab = ylab,
         asp = 1, col = colRamp, xaxt = xaxt, yaxt = yaxt )
}

polyMap <- function( polyList, type, lwd = 1, xlim = NULL, ylim = NULL,
                     border = terrain.colors(length(unique(type)), alpha = .5 ), 
                     col = NULL, xaxt = 's', yaxt = 's', 
                     xlab = '', ylab = '', add = F ){
  
  # polyList is a list of polygons with columns x, y
  # type is used for the color legend
  # border and col have length equal to the number of unique values in type
  
  subs <- unique( type )
  if( length( border ) == 1 )border <- rep( border, length( subs ) )
  
  names( border ) <- subs
  cc <- fill <- NA
  if( !is.null( col ) ){
    fill <- col
    names( fill ) <- subs
  }
  
  if( !add ){
    if( is.null( xlim ) ){
      sapply( polyList, range )
      rr <- sapply( polyList, function(x) apply( x, 2, range ) )
      xlim <- range( rr[1:2,] )
      ylim <- range( rr[3:4,] )
    }
    plot( NA, xlim = xlim, ylim = ylim, asp = 1, bty = 'n',
          xlab = xlab, ylab = ylab, xaxt = xaxt, yaxt = yaxt )
  }
  
  for( k in 1:length( polyList ) ){
    kbord <- border[ type[k] ]
    if( !is.null( col ) )cc <- fill[ type[k] ]
    
    polygon( polyList[[k]][,1], polyList[[k]][,2], border = kbord, col = cc,
             lwd = lwd )
  }
}



rampColor <- function( values, 
                       ramp = c('#d53e4f','#fc8d59','#e6f598','#99d594','#3288bd'),
                       nlevs = 10 ){
  # assign colors to the magnitude of values
  
  s  <- seq( min( values, na.rm = T ), max( values, na.rm = T ), length = nlevs )
  si <- findInterval( values, s, all.inside = T )
  cfun <- colorRampPalette( ramp )
  cfun( nlevs )[ si ]
}

colorSegment <- function(xx, yy, colors, nc = 20, lwd=1, lty=1,
                          clear=T){
  
  # add line to plot, colored by magnitude
  # colors - colors to interpolate, see colorRampPalette
  #          interpolate for range(y)
  # x, y   - vectors to plot
  # clear  - clear background of line
  
  yr    <- range(yy)
  sc    <- seq(yr[1],yr[2],length=nc)
  yc    <- c(1:nc)
  y1    <- 1:(nc-1)
  y2    <- y1 + 1
  sc    <- round( yr[1] + yc/nc*diff(yr) ,0)
  
  cf   <- colorRampPalette(colors)(nc)
  ycol <- findInterval(yy,sc,all.inside=T)
  
  if(clear) lines( new$density, new$mids, lwd=lwd*3, col='white')
  segments(xx[y1], yy[y1], xx[y2], yy[y2], col=cf[ycol],lwd=lwd,lty=lty)
}

vec2Mat <- function( x, i, j){
  I <- sort(unique(i))
  J <- sort(unique(j))
  mat <- matrix( NA, length(I), length(J), dimnames = list( I, J) )
  mm <- cbind( match(i, I), match(j, J) )
  mat[ mm ] <- x
  mat
}

mat2Vec <- function( mat, include.na = F ){  # reverse: matrix back to vector
  
  if(is.null(colnames(mat)))colnames(mat) <- c(1:ncol(mat))
  if(is.null(rownames(mat)))rownames(mat) <- c(1:nrow(mat))
  
  i   <- rownames(mat)
  j   <- colnames(mat)
  ii  <- rep( i, length(j) )
  jj  <- rep( j, each = length(i) )
  id  <- paste(ii, jj, sep='_')
  vec <- mat[ cbind(ii, jj) ]
  names( vec ) <- id
  w   <- which(is.finite(vec))
  df  <- data.frame( i = ii, j = jj, x = vec)
  if( !include.na )df <- df[ is.finite(df[,3]), ]
  df
}

.dMVN <- function(xx,mu,smat=NULL,sinv=NULL,log=F){ #MVN density for mean 0
  
  if(!is.matrix(xx))xx <- matrix(xx,1)
  if(!is.matrix(mu))mu <- matrix(mu,1)
  
  xx <- xx - mu
  if(ncol(xx) != nrow(smat))xx <- t(xx)
  
  if(!is.null(sinv)){
    distval <- diag( xx%*%sinv%*%t(xx) )
    ev      <- eigen(sinv, only.values = TRUE)$values
    logd    <- -sum(log(ev))
  }
  
  if(is.null(sinv)){
    testv <- try(chol(smat),T)
    if(inherits(testv,'try-error')){
      tiny  <- min(abs(xx))/100 + 1e-5
      smat  <- smat + diag(diag(smat + tiny))
      testv <- try(chol(smat),T)
    }
    covm    <- chol2inv(testv)
    distval <- rowSums((xx %*% covm) * xx)
    ev      <- eigen(smat, only.values = TRUE)$values 
    logd    <- sum(log( ev ))
  }
  
  z <- -(ncol(xx) * log(2 * pi) + logd + distval)/2
  if(!log)z <- exp(z)
  z
}

biVarMoments <- function(x1, x2, wt=1,PLOT = F, POINTS=F, color=1, pr = .95, ADD=F, lty=1,lwd=1,
                         xylab=c('x','y')){  
  
  require(cluster)
  
  #x1, x2 - vectors for variables, wt is weight
  
  if(length(pr) > 1){
    if(length(lty) == 1)lty = rep(lty,length(pr))
    if(length(lwd) == 1)lwd = rep(lwd,length(pr))
    if(length(color) == 1)color = rep(color,length(pr))
  }
  
  if(length(wt) == 1)wt <- rep(wt,length(x1))
  
  ww <- which(is.finite(x1) & is.finite(x2))
  
  x1 <- x1[ww]
  x2 <- x2[ww]
  wt <- wt[ww]
  
  w1 <- x1*wt
  w2 <- x2*wt
  m1 <- sum(w1)/sum(wt)
  m2 <- sum(w2)/sum(wt)
  
  v1  <- sum(wt*x1^2)/sum(wt) - m1^2
  v2  <- sum(wt*x2^2)/sum(wt) - m2^2
  c   <- sum(wt*(x1 - m1)*(x2 - m2))/sum(wt)
  
  for(k in 1:length(pr)){
    tmp <- list(loc = c(m1,m2), cov = matrix(c(v1,c,c,v2),2,2), d2 = qchisq(pr[k],1) )
    tmp <- predict.ellipsoid(tmp)
    if(PLOT & !ADD & k == 1)plot(tmp[,1],tmp[,2],type='l',col=color[k],lty=lty[k],
                                 lwd=lwd[k],xlab=xylab[1], ylab=xylab[2])
    if(POINTS)points(x1,x2,cex=.3,col='grey')
    if(ADD | k > 1)lines(tmp[,1],tmp[,2],type='l',col=color[k],lwd=lwd[k],lty=lty[k])
  }
  invisible( list(loc = c(m1,m2), 
                  cov = matrix(c(v1,c,c,v2),2,2) ,
                  d2 = qchisq(pr[1],1), ellipse = tmp) )
}

metRatio <- function(beta, x, y, likelihood,
                     priorB=beta*0, priorVB=diag(1000,length(beta))){
  B <- T
  if(likelihood %in% c('poisson', 'dpois') )B <- F
  
#  link  <- match.fun(link)
#  like  <- match.fun(likelihood)
  xb <- x %*% beta
  
  if(B){
    theta <- exp( xb )/( 1 + exp( xb ) )
    p1    <- dbinom(y, 1, theta, log = T)
  }else{
    theta <- exp( xb )
    p1 <- dpois(y, theta, log = T)
  }
  sum( p1 ) + .dMVN(beta, priorB, priorVB, log=T)
}

updateBetaGLM <- function(bg, cmat, x, y, likelihood='binom',
                          lo = NULL, hi=NULL, ...){
  
  # metropolis sampler for GLM reg parameters
  # lo, hi - truncation points for MVN proposals
  
  require(TruncatedNormal)
  
  ac <- 0
  if(is.null(lo)){             # proposal
    bs <- t( .rMVN(1,bg,cmat) )      
  }else{
    bs <- matrix( rtmvnorm(1, bg, cmat, lb = lo, ub = hi), ncol = 1)
  }
  
  pnow <- metRatio(bg, x, y, likelihood )
  pnew <- metRatio(bs, x, y, likelihood )
  a    <- exp(pnew - pnow)
  z    <- runif(1,0,1)
  if(a > z){
    bg <- bs
    ac <- 1
  }
  list(beta = bg, accept = ac)
}


columnSplit <- function( vec, sep = '_', ASFACTOR = F, ASNUMERIC = FALSE, 
                         LASTONLY = FALSE ){
  
  vec <- as.character( vec )
  nc  <- length( strsplit( vec[ 1], sep, fixed = TRUE )[[ 1]] )
  
  mat <- matrix( unlist( strsplit( vec, sep, fixed = TRUE ) ), ncol = nc, byrow = TRUE )
  if( LASTONLY & ncol( mat ) > 2 ){
    rnn <- mat[, 1]
    for( k in 2:( ncol( mat )-1 ) ){
      rnn <- columnPaste( rnn, mat[, k] )
    }
    mat <- cbind( rnn, mat[, ncol( mat )] )
  }
  if( ASNUMERIC ){
    mat <- matrix( as.numeric( mat ), ncol = nc )
  }
  if( ASFACTOR ){
    mat <- data.frame( mat )
  }
  if( LASTONLY )mat <- mat[, 2]
  mat
}

 
varPrior <- function(mu, wt){ 
  # mean and wt of prior IG distribution, returns prior parameter values for IG
  c(wt, mu*(wt - 1))
}



initialStatesSS <- function(y){
  
  require(stats)
  
  if(!is.matrix(y))y <- matrix(y,ncol=1)
  r <- ncol(y)
  
  n    <- length(y)
  time <- c(1:n)
  wm   <- which(is.na(y))
  notMiss <- c(1:n)
  
  x <- y
  
  if(length(wm) > 0){
    notMiss <- notMiss[-wm]
    x[wm]   <- predict(smooth.spline(time[-wm],y[-wm]),time)$y[wm]
  }
  
  list(x = as.vector(x), notMiss = notMiss, miss = wm)
}


updateSSRW <- function( states, yy, zb=rep(0, length(y)), missing, tg, sg ){        
  #state-space random walk 
  #update continuous states, random walk
  #missing times, obs y, obs error tg, process error sg
  #zb = z%*%beta
  
  nn <- length(states)
  
  for(t in 1:nn){
    
    VI <- 0
    v  <- 0
    
    if(!t %in% missing){          #observations
      v  <- yy[t]/tg
      VI <- 1/tg
    }
    
    if(t < nn){              #t+1 term excluded for last 
      v  <- v + (states[t+1] - zb[t])/sg 
      VI <- VI + 1/sg
    }
    
    if(t > 1){                #t-1 term excluded for 1st 
      v  <- v + (states[t-1] + zb[t-1])/sg
      VI <- VI + 1/sg
    }
    
    V     <- 1/VI
    states[t] <- rnorm(1,V*v,sqrt(V))
  }
  states
}


updateVariance <- function(yy,mu,s1=1,s2=1){
  
  #y and mu are n by 1
  
  u1 <- s1 + length(yy)/2
  u2 <- s2 + .5*crossprod(yy - mu)
  
  1/rgamma(1,u1,u2) 
}

updateSSbeta <- function(states,z,sg,priorB,priorIVB,addStates=F){  
  #covariates z
  #if addStates = T:  states[t+1] <- states[t] + zbeta[t]
  
  require(mvtnorm)
  
  nt <- length(states)
  xx <- states[-1]
  zz <- z
  if(nrow(zz) == nt)zz <- zz[-nt,]     # last value already clipped?
  
  V  <- solve( crossprod(z[-nt,])/sg + priorIVB )
  
  if(addStates)xx <- xx - states[-nt]
  
  v  <- 1/sg*crossprod(zz,xx)
  t( rmvnorm(1,V%*%v,V) )
}




bayesReg <- function(formula, data, ng = 3000, burnin = 100, TOBIT=NULL){
  
  fc <- as.character(formula)
  
  if(missing(data)){
    data <- environment(formula)
    yx   <- match.call()
    m    <- match(c("formula", "data"), names(yx), 0L)
    yx   <- yx[c(1L, m)]
    yx[[1L]] <- quote(stats::model.frame)
    yx   <- eval(yx, parent.frame())
    y    <- yx[,1]
  } else {
    yy <- unlist( strsplit( fc, '~' ) )
    yy <- yy[ nchar(yy) > 0]
    y  <- data[,yy[1]]
  }
  
  yzero <- which (y == 0)
  nzero <- length(yzero)
  
  if( is.null(TOBIT) )TOBIT <- F
  
  if(TOBIT)message('fitted as Tobit model')
  
  tmp <- model.frame(formula, data, na.action=NULL)
  x   <- model.matrix(formula, data=tmp)
  
  colnames(x)[1] <- 'intercept'
  
  xnames    <- colnames(x)
  snames    <- colnames(y)
  Q         <- ncol(x)
  n         <- nrow(x)
  predXcols <- 2:Q
  
  ymiss <- which(is.na(y))
  
  yy <- y
  xx <- x
  if(length(ymiss) > 0){
    yy <- yy[-ymiss]
    xx <- xx[-ymiss,]
  }
  
  XX  <- crossprod(xx)
  IXX <- solve(XX)
  bg  <- IXX%*%crossprod(xx,yy)
  py  <- x%*%bg
  w   <- y
  w[ymiss] <- mean(y,na.rm=T)
  if(TOBIT){
    p0 <- py[y == 0]
    p0[p0 > 0] <- -p0[p0 > 0]
    py[y == 0] <- p0
    w[w == 0] <- py[y == 0]
  }
  wi  <- c(which(y == 0),ymiss)   
  
  priorB   <- matrix(0,Q,1)         #prior mean regression parameters
  priorIVB <- solve(diag(1000,Q))   #inverse prior covariance matrix
  s1       <- .1                    #variance prior values
  s2       <- .1
  
  bchains <- matrix(NA,ng,Q)
  schains <- rep(0,ng)              #store parameters
  colnames(bchains) <- xnames
  ychains <- matrix(0,ng,n)
  mu <- x%*%bg
  
  for(g in 1:ng){
    
    sg <- updateSigma(w, mu, s1, s2)
    bg <- updateBeta(xx, w, sg, priorIVB, priorB )
    mu <- x%*%bg
    if(TOBIT)w[yzero] <- .tnorm(nzero,-Inf,0,mu[yzero],sqrt(sg))
    w[ymiss] <- py[ymiss]
    
    py <- rnorm(n,mu,sqrt(sg))
    py[py <= 0] <- 0
    
    bchains[g,] <- bg   # store estimates
    schains[g]  <- sg
    ychains[g,] <- py
  }
  
  beta <- signif( t( apply(bchains,2,quantile,c(.5,.025,.975)) ), 4)
  bhat <- beta[drop=F,,1]
  bsd  <- signif( apply(bchains,2,sd ), 4)
  beta <- data.frame(beta[,1], bsd, beta[,2:3])
  betaVar <- var(bchains)
  
  shat <- mean(schains)
  sse  <- sd(schains)
  
  zero <- which(beta[,3] < 0 & beta[,4] > 0)
  notZero <- rep('*',Q)
  notZero[ zero ] <- ' '
  beta <- data.frame( cbind(beta,notZero) )
  
  colnames(beta) <- c('median','std error','0.025','0.975','not zero')
  
  py <- signif( t( apply( ychains, 2, quantile,c(.5,.025,.975) ) ), 3)
  py <- cbind(y,py)
  
  rmspe <- signif( sqrt( mean((y - py[,2])^2, na.rm=T) ), 4)
  sig   <- signif(sqrt(shat),4)
  
  ssr   <- crossprod(y - x%*%bhat)/(n - Q)
  rsq   <- signif(1 - ssr/var(y), 3)
  
  cat("\nCoefficients:\n")
  print( beta )
  
  cat("\n * indicates that 95% predictive interval does not include zero\n")
  
  out <- paste("\nResidual standard error ", sig, ", with ",n - Q, 
               " degrees of freedom, \n root mean sq prediction error ",
               rmspe, ".", sep='')
  cat(out,"\n")
  
  
  list(beta = beta, betaVar = betaVar, predictY = py, sigma = median(schains), 
       rmspe = rmspe)
}

#################################################################
simX <- function(n,loX,hiX){                #generate design matrix
  
  k <- length(loX)
  x <- matrix(1,n,k)
  for(j in 1:k)x[,j] <- runif(n,loX[j],hiX[j])
  x
}
####################################################
simY <- function(x,b,LIKE,r = 1,size=rep(1,nrow(x)),sigma = 0,Effort = 1){     #simulate response
  
  u <- x%*%b
  
  if(LIKE == 'norm')return( rnorm(length(u),u,sqrt(sigma)) )
  if(LIKE == 'pois'){
    u <- exp(u)*Effort
    return( rpois(nrow(x),u) )
  } 
  if(LIKE == 'binom'){
    u <- invlogit(u)
    return( rbinom(nrow(x),1,u) )
  }
  if(LIKE == 'multinom'){
    zs <- rowSums(exp(u))
    z1 <- 1/(1 + zs)
    zm <- exp(u)/ (1 + zs)
    u  <- cbind(zm,z1)
    return( myrmultinom(size,u) )
  }
  if(LIKE == 'mvnorm'){
    u <- myrmvnorm(n,u,sigma)
    return( u )
  }
  if(LIKE == 'mvnorm-multinom'){
    u <- myrmvnorm(n,u,sigma)
    zs <- rowSums(exp(u))
    z1 <- 1/(1 + zs)
    zm <- exp(u)/ (1 + zs)
    u  <- cbind(zm,z1)
    return( myrmultinom(size,u) )
  }
  numeric(0)
}

#######################################################
deviance <- function(y,x,b,s=0,LIKE='norm'){
  
  if(LIKE == 'norm')  dv <- dnorm(y,x%*%b,sqrt(s),log=T) 
  if(LIKE == 'pois')  dv <- dpois(y,exp(x%*%b),log=T)
  if(LIKE == 'binom') dv <- dbinom(y,1,invlogit(x%*%b),log=T)
  if(LIKE == 'mvnorm')dv <- dmvnormZeroMean(y - x%*%b,s)
  if(LIKE == 'multinom')dv <- multinomLike(y,x,b)$like
  -2*dv
}

gibbsLoop <- function(LIKE, ng, x, y, b, priorB=b*0, priorVB=diag(1000,length(b)),
                      sigma = 0, z = numeric(0), burnin=10,
                      loB=NULL, hiB=NULL ){
  
  if(ng < 100)ng <- 100
  if(burnin > ng/2) burnin <- round(ng/2)
  
  tiny <- 1e-13
  
  k <- ncol(x)
  r <- br <- 1                       #responses
  if(is.matrix(y)){
    r <- ncol(y)
    if(LIKE == 'multinom')br <- r - 1
  }
  n  <- nrow(x)
  kx <- length(b)
  if(!is.matrix(b))b <- as.matrix(b)
  
  priorIVB <- solve(priorVB)
  
  bgibbs <- matrix(NA,ng,kx)
  if(is.null(rownames(b))) rownames(b) <- paste('q',c(1:k),sep='-')
  if(is.null(colnames(b))) colnames(b) <- paste('s',c(1:br),sep='-')
  colnames(bgibbs) <- as.vector(outer(rownames(b),colnames(b),paste,sep='_') )
  
  sgibbs <- matrix(NA,ng,max(1,length(sigma)))
  
  if(sigma[1] == 0)sgibbs <- numeric(0) #variances
  
  sg <- sigma
  
  if(is.matrix(sigma)){                  #Wishart prior
    sg        <- prior.W
    colnames(sgibbs) <- outer(rownames(sigma),rownames(sigma),paste,sep='_') 
    colnames(bgibbs) <- as.vector(outer(colnames(x),colnames(b),paste,sep='_'))
  }
  
  pBVar <- solve(crossprod(x))
  
  pred <- pred2 <- rep(0,nrow(x)*r)
  
  bg <- b
  yg <- y
  y1 <- yg
  
  if(LIKE == 'pois') bg <- pBVar%*%crossprod(x,log(y + .1))
  if(LIKE == 'binom'){
    yl <- y
    yl[yl == 0] <- .1
    yl[yl == 1] <- .9
    bg <- pBVar%*%crossprod(x,logit(yl))
  }
  
  if(LIKE == "mvnorm-multinom"){
    y1   <- prob2Logit(yg)
  }
  if(LIKE == 'multinom'){
    size <- rowSums(y)
    pBVar   <- diag(.0001,k*(r-1))
    priorB  <- b*0
    priorVB <- diag(100,k*(r-1))
    priorIVB <- solve(priorVB)
    sigma <- sg <- sgibbs <- numeric(0)
  }
  
  dev <- 0
  np  <- 0
  
  pbar <- txtProgressBar(min=1,max=ng,style=1)
  
  for(g in 1:ng){
    
    bg <- bUpdateGibbs(LIKE,x,y1,bg,priorB,priorIVB,loB,hiB,sigma=sg,pBVar)
    
    if( sigma[1] > 0 & length(sigma) == 1 ){
      sg <- updateSigma( y1, x%*%bg, 1, 1)
    }
    
    if(length(grep('mvnorm',LIKE)) > 0){
      sinv <- wishsamp(x,y1,bg)
      sg   <- solve(sinv)
    }
    
    if(LIKE == "mvnorm-multinom"){
      y1 <- ysampMvnormMultinom(x,y1,z,bg,sg)$y
      yg <- logit2Prob(y1)
    }
    
    dev <- dev + sum( deviance( y1, x, bg, sg, LIKE) )
    
    if(g > burnin){
      py    <- as.vector(simY(x,bg,LIKE,r,size,sg))
      pred  <- pred + py
      pred2 <- pred2 + py^2
      np    <- np + 1
    }
    
    bgibbs[g,] <- bg
    if(length(sgibbs) > 0)sgibbs[g,] <- sg
    
    if(g %in% c(100,200,500,1000)){
      pBVar <- .5*cov(bgibbs[(g-50):g,]) + diag(tiny,ncol(bgibbs))
    }
    
    setTxtProgressBar(pbar,g)
  }
  
  ymean <- pred/np
  yse   <- sqrt(pred2/np - ymean^2)
  
  if(r > 1){
    ymean <- matrix(ymean,n,r)
    yse   <- matrix(yse,n,r)
  }
  bmean <- colMeans(bgibbs)
  bmean <- matrix(bmean,nrow(bg),ncol(bg))
  smean <- numeric(0)
  if(length(sg) > 1) smean <- matrix(colMeans(sgibbs),nrow(sg),ncol(sg))
  if(length(sg) == 1)smean <- mean(sgibbs)
  
  meanDev <- dev/ng
  
  pd  <- meanDev - sum(deviance(y1,x,bmean,smean,LIKE))
  dic <- 2*pd + meanDev
  
  colnames(ymean) <- colnames(yse) <- colnames(y)
  
  list(bgibbs = bgibbs,sgibbs = sgibbs, ymean = ymean, yse = yse, dic = dic)
}


logit2Prob <- function(y){     #multivar logit to fractions
  
  if(!is.matrix(y)){
    z <- invlogit(y)
    return(cbind(z,1 - z))
  }
  
  zs   <- rowSums(exp(y))
  z1   <- 1/(1 + zs)
  zm   <- exp(y)/ (1 + zs)
  cbind(zm,z1)
}

prob2Logit <- function(z){     #fractions to multivar logit
  
  r <- ncol(z)
  if(r == 2){
    ss <- z[,1]/rowSums(z)
    return(log(ss/(1 - ss) ))
  }
  
  log(z[,-r]/(1 - rowSums(z[,-r])))
}

bUpdateMNom <- function(x,y,b,priorB,priorIVB,
                        pBVar=diag(.1,nrow(b)*(ncol(b)-1)),loB=NULL,hiB=NULL){
  
  bvec <- as.vector(b)
  priorVB <- solve(priorIVB)
  
  tmp  <- bmultiProp(ncol(b)+1,nrow(b),b,pBVar,loB,hiB)
  cc   <- tmp$cc
  cvec <- tmp$cvec
  
  pnow <- multinomLike(y,x,b)$like
  pnew <- multinomLike(y,x,cc)$like
  
  pnow <- sum(pnow) + .dMVN(bvec,as.vector(priorB),priorVB,log=T)
  pnew <- sum(pnew) + .dMVN(cvec,as.vector(priorB),priorVB,log=T)
  
  a <- exp(pnew - pnow)
  z <- runif(1,0,1)
  if(z < a)b <- cc
  b
}

#######################################################
multinomLike <- function(y,x,b){  #log likelihood multinomial logit
  
  tiny <- 1e-20
  huge <- 1 - tiny
  
  z2 <- logit2Prob(x%*%b)
  
  z2[z2 < tiny] <- tiny
  z2[z2 > huge] <- huge
  list(like = y*log(z2), theta = z2)
}

bUpdateGibbs <- function(LIKE,x,y,b,priorB,priorIVB,
                         loB=NULL,hiB=NULL,sigma = 0,pBVar=0){
  
  if(LIKE == 'norm')    return( bUpdateNorm(x,y,b,priorB,priorIVB,loB,hiB,sigma) )
  if(LIKE == 'multinom')return( bUpdateMNom(x,y,b,priorB,priorIVB,pBVar,loB,hiB) )
  if(LIKE %in% c('mvnorm','mvnorm-multinom'))return( bUpdateMVNorm(x,y,b,sigma) )
  
  b <- matrix(b,length(b),1)
  if( is.null(loB)) c <- t(myrmvnorm(1,t(b),pBVar))   #proposal
  if(!is.null(loB)) c <- tnorm.mvt(b,b,pBVar,loB,hiB,times=1)
  
  znow <- x%*%b
  znew <- x%*%c
  
  if(LIKE == 'pois'){
    pnow <- dpois(y,exp(znow),log=T)
    pnew <- dpois(y,exp(znew),log=T)
  }
  if(LIKE == 'binom'){
    pnow <- dbinom(y,1,invlogit(znow),log=T)
    pnew <- dbinom(y,1,invlogit(znew),log=T)
  }
  
  priorVB <- solve(priorIVB)
  
  pnow <- sum(pnow) + mydmvnorm(t(b),priorB,priorVB,log=T)
  pnew <- sum(pnew) + mydmvnorm(t(c),priorB,priorVB,log=T)
  
  a <- exp(pnew - pnow)
  z <- runif(1,0,1)
  if(z < a)b <- c
  b
}


bmultiProp <- function( r, k, b = matrix( 0, k, r-1),
                       pBVar = diag(.1,k*(r-1) ), loB=NULL, hiB=NULL ){  
  
  bvec <- as.vector(b)
  cvec <- .rMVN(1,t(bvec),pBVar)
  cc   <- matrix(cvec,nrow(b),ncol(b))
  
  if(!is.null(loB)){                         #if lob and hib available, use tnorm.mvt
    cvec <- tnorm.mvt(bvec,bvec,pBVar,loB,hiB,times=4)
    cc   <- matrix(cvec,nrow(b),ncol(b))
  }
  
  list(cc = cc, cvec = cvec)
}

updateWishart <- function(xx, yy, df, beta=NULL, IXX=NULL){
  #df  <- n - Q + S - 1
  S     <- ncol(yy)
  index <- 0
  if(is.null(IXX)){
    XX  <- crossprod(xx)
    IXX <- solve(XX)
  }
  
  D  <- diag(1,nrow(xx)) - xx%*%IXX%*%t(xx)
  SS  <-  t(yy)%*%D%*%yy
  testv <- try(chol(SS),T)
  
  if( inherits(testv,'try-error') ){
    tiny <- 1e-8
    SS[SS < tiny] <- tiny
    message('warning: updateWishartNoPrior')
    
    if( is.null( beta ) ){
      XC <- yy - xx
    }else{
      XC <- yy - xx%*%beta
    }
    
    SS <- crossprod( XC ) +  diag(diag(SS)*.001)#*nrow(SS)
    SS <- SS + diag(diag(SS)*.1)
    testv <- try(chol(SS),T)
    index <- 1
  }
  SI <- chol2inv(testv)
  z  <- matrix(rnorm(df*S),df,S)%*%chol(SI)
  
  sinv  <- crossprod(z)
  sigma <- solve(sinv)
  list( sigma = sigma, sinv = sinv, indicator = index )
}

getVv <- function(x, y, sigma2, priorIVB, priorB){
  
  si <- 1/sigma2
  V  <- solve( si*crossprod(x) + priorIVB )
  v  <- si*crossprod(x,y) + priorIVB%*%priorB
  list( V = V, mu = V%*%v )
}


updateZ <- function( x, y, z, alpha, beta, Sigma, tau ){
  
  # update latent states for 
  #     z[t+1,] = alpha z[t,] + beta x[t,] + epsilon
  #     y[t,] = z[t,] + obs error
  
  SI <- solve( Sigma )
  TI <- solve( tau )
  
  for( j in 1:n ){
    
    VI <- TI
    v  <- TI%*%y[j,]
    
    if( j > 1 ){
      v  <- v + SI%*%( alpha%*%z[j-1,] + t(beta)%*%x[j-1,] )
      VI <- VI + SI
    }
    if( j < n ){
      v  <- v + t( alpha )%*%SI%*%( z[j+1,] - t(beta)%*%x[j,] )
      VI <- VI + t( alpha )%*%SI%*%alpha 
    }
    V <- solve( VI )
    z[j,] <- t(V%*%v) + .rMVN(1, 0, V )
  }
  z
}


updateBetaMVN <- function( x, y, sg){
  
  XX  <- crossprod( x )
  IXX <- chol2inv(chol( XX ) )
  WX  <- crossprod( x, y)
  WIX <- IXX%*%WX
  matrix( .rMVN(1,as.vector(WIX),
                kronecker(sg,IXX)),nrow(IXX),ncol(WIX) )
}

updateBeta <- function( x, y, sigma, 
                        priorIVB = diag(.001, ncol(x) ), priorB = rep(0, ncol(x) ) ){  # random vector of coefficients
  
  Vv <- getVv(x, y, sigma, priorIVB, priorB)
  V  <- Vv$V
  mu <- Vv$mu
  t( .rMVN(1, mu, V) )                     
}

updateSigma <- function(y, mu, s1, s2){ # random value for residual variance
  n  <- length(y)
  u1 <- s1 + n/2
  u2 <- s2 + 0.5*crossprod(y - mu)
  1/rgamma(1,u1,u2)                          
}

imputeX <- function(missx, xx, yy, bg, sg, priorx ){
  
  ix <- sort( unique(missx[,2]) )
  IV <- 1/1000
  
  for(i in ix){
    
    wmi <- missx[missx[,2] == i,1]
    
    V <- 1/( bg[i]^2/sg + IV )
    v <- bg[i]*(yy[wmi+1] - xx[wmi,-i,drop=F]%*%bg[drop=F,-i])/sg + 
      priorx[i]*IV
    xx[wmi,i] <- v*V
  }
  xx
}

getMissX <- function(x, y, first=NULL, last=NULL){
  
  prows <- rrows <- numeric(0)
  missX <- missY <- numeric(0)
  missx <- which(is.na(x),arr.ind=T)
  missy <- which(is.na(y))
  
  if( !is.null(first) ){        # first, last bound groups
    ngroup <- length(first)
    for(j in 1:ngroup){           #rows to predict and respond
      prows <- c(prows,c( first[j]:(last[j]-1)) )
      rrows <- c(rrows,c( (first[j]+1):last[j] ) )
      jr <- first[j]:last[j]
      wna <- which( is.na(x[jr,]), arr.ind=T ) 
      if( length( wna ) > 0 )missX <- append(missX, list(wna) )
      wna   <- which( is.na( y[jr] ), arr.ind = T )
      missY <- append(missY, list( wna ) )
    }
    missXbyGroup <- which(is.na(x[prows,]),arr.ind=T)
    missYbyGroup <- which(is.na(y[rrows]))
  }
  
  list(prows = prows, rrows = rrows, missx = missx, 
       missy = missy,
       missX = missX, 
       missY = missY,
       missXbyGroup = missXbyGroup, 
    missYbyGroup = missYbyGroup)
}


.rMVN <- function (nn, mu, sigma){
  
  # nn - no. samples from one mu vector or nrow(mu) for matrix
  
  if(!is.matrix(mu)) mu <- matrix(mu,1)
  if(length(mu) == 1)mu <- matrix(mu,1,nrow(sigma))
  if(ncol(mu) == 1)  mu <- t(mu)
  
  m <- ncol(sigma)
  
  if(ncol(mu) != m)stop('dimension mismatch mu, sigma')
  
  if(nn > 1 & nrow(mu) == 1)mu <- matrix(mu,nn,m,byrow=T)
  
  if(nn != nrow(mu))stop('sample size does not match mu')
  
  testv <- try(svd(sigma),T)
  
  if( inherits(testv,'try-error') ){
    ev <- eigen(sigma, symmetric = T)
    testv <- t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
  } else {
    testv <- t(testv$v %*% (t(testv$u) * sqrt(testv$d)))
  }
  
  p <- matrix(rnorm(nn * m), nn) %*% testv
  p + mu
}


logit <- function(x) log(x/(1 - x))
invlogit <- function(x) exp(x)/(1 + exp(x))

.tnorm <- function(n,lo,hi,mu,sig, tiny = .000001){   
  
  #normal truncated lo and hi
  
  tiny <- 10e-6
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig)
  q2 <- pnorm(hi,mu,sig) 
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  
  z[z == Inf]  <- lo[z == Inf] + tiny
  z[z == -Inf] <- hi[z == -Inf] - tiny
  z
}

chainPlot <- function( chain, trueValues = NULL, burnin = 1 ){
  
  if( !is.matrix( chain ) )chain <- matrix( chain, ncol = 1 )
  nc     <- ncol( chain )
  nr     <- nrow( chain )
  cnames <- colnames( chain )
  if( is.null( cnames ) )cnames <- paste( 'p', 1:nc ) 
  
  sx <- c( burnin, nr )
  
  stats <- chain2tab( chain[ sx[1]:sx[2],] )
  
  mf <- .getPlotLayout( nc )
  par( mfrow = mf$mfrow, bty = 'n', mar = c( 3, 3, 2, 1),
       omi = c(.5, .1, .1,.1) )
  
  for( k in 1:ncol(chain) ){
    plot( chain[,k], type = 'l', xlab = '', ylab = '' )
    title( cnames[k] )
    lines( sx, stats[k, c(1,1) ], lwd = 2, col = 'grey' )
    lines( sx, stats[k, c(3,3) ], lwd = 3, col = 'white' )
    lines( sx, stats[k, c(3,3) ], lty = 2 )
    lines( sx, stats[k, c(4,4) ], lwd = 3, col = 'white' )
    lines( sx, stats[k, c(4,4) ], lty = 2 )
    if( !is.null( trueValues ) )
      lines( sx, trueValues[ c(k,k) ], col = 'brown', lwd = 2 )
  }
  mtext( 'Iteration', 1, outer = T, line = -1 )
}


chain2density <- function( chain, trueValues = NULL, burnin = 1 ){

  if( !is.matrix( chain ) )chain <- matrix( chain, ncol = 1 )
  nc     <- ncol( chain )
  nr     <- nrow( chain )
  cnames <- colnames( chain )
  if( is.null( cnames ) )cnames <- paste( 'p', 1:nc ) 
  
  sx <- c( burnin, nr )
  xd <- sx[1]:sx[2]
  
  stats <- chain2tab( chain[ xd,] )
  
  mf <- .getPlotLayout( nc )
  par( mfrow = mf$mfrow, bty = 'n', omi = c(.1, .5, .1,.1) )
  
  for( k in 1:ncol(chain) ){
    xy <- density( chain[ xd, k ] )
    plot( xy$x, xy$y, type = 'l', xlab = cnames[k], ylab = '', lwd = 2 )
    abline( v = stats[k,1], col = 'grey', lwd = 2 )
    abline( v = stats[k,3:4], lty = 2 )
    if( !is.null(trueValues) )abline( v = trueValues[k], col = 'brown', lwd = 2 )
  }
  mtext( 'Density', 2, outer = T, line = -1, las = 0 )
}

chain2tab <- function( chain, sigfigs = 4 ){
  
  if( !is.matrix( chain ) )chain <- matrix( chain, ncol = 1 )
  
  mu <- colMeans( chain )    
  
  SE <- apply( chain, 2, sd, na.rm=T )
  CI <- apply( chain, 2, quantile, c( .025, .975 ), na.rm=T )
  splus <- rep( ' ', length = length( SE ) )
  splus[ CI[ 1, ] > 0 | CI[ 2, ] < 0] <- '*'
  
  tab <- cbind( mu, SE, t( CI ) )
  tab <- signif( tab, sigfigs )
  colnames( tab ) <- c( 'estimate', 'SE', 'CI_025', 'CI_975' )
  tab <- as.data.frame( tab )
  tab$sig95 <- splus
  attr( tab, 'note' ) <- '* indicates that zero is outside the 95% CI'
  
  tab
}
