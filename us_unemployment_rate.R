library(raster)
library(maptools)
library(rgeos)
library(rgdal)
library(akima)
library(reshape)
library(prevR)

options(scipen=999)

source("transition-map-helpers.R")


# Albers projection, proj4string
county_proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"



#
# Load data.
#

# Annual unemployment
unemp <- read.csv("data/unemployment-combined.tsv", sep="\t", colClasses=c(state_fips="character", county_fips="character"), stringsAsFactors=FALSE)
unemp <- unemp[ unemp$state_fips != "02" & unemp$state_fips != "15" & unemp$state_fips != "72",  ]

# Lower US boundaries
usa <- readShapePoly("data/nation/nationp010g/nationp010g.shp", proj4string=CRS("+proj=longlat"))
usa_lower <- usa[2,]
usa_lower_albers <- spTransform(usa_lower, CRS(county_proj))

# US states
states <- readShapePoly("data/states/cb_2013_us_state_20m/cb_2013_us_state_20m.shp", proj4string=CRS("+proj=longlat"))
states_albers <- spTransform(states, CRS(county_proj))

# Lakes
lakes_world <- readShapePoly("data/lakes/ne_50m_lakes/ne_50m_lakes.shp", proj4string=CRS("+proj=longlat"))
lakes_world_albers <- spTransform(lakes_world, CRS(county_proj))
lakes_albers <- intersect(lakes_world_albers, usa_lower_albers)

# Counties
counties <- readShapePoly("data/counties/cb_2013_us_county_20m/cb_2013_us_county_20m.shp", proj4string=CRS("+proj=longlat"))
counties <- counties[counties$STATEFP != "02" & counties$STATEFP != "15" & counties$STATEFP != "72", ]
counties_albers <- spTransform(counties, CRS(county_proj))

# Approximate center of each county
centroids <- SpatialPointsDataFrame(gCentroid(counties_albers, byid=TRUE), counties@data, match.ID=FALSE)



#
# For filtering later.
# 

xo <- seq(usa_lower_albers@bbox["x", "min"], usa_lower_albers@bbox["x", "max"], length = 600)
yo <- seq(usa_lower_albers@bbox["y", "min"], usa_lower_albers@bbox["y", "max"], length = 400)
metric_mat <- matrix(NA, nrow=length(xo), ncol=length(yo), dimnames = list(xo, yo))
metric_df <- melt(metric_mat)
metric_df$in_simple <- point.in.SpatialPolygons(metric_df[,1], metric_df[,2], usa_lower_albers)
usa_lower_mat <- matrix(metric_df$in_simple, nrow=length(xo), ncol=length(yo), dimnames = list(xo, yo)) 


# Full FIPS code
unemp$full_fips <- paste(unemp$state_fips, unemp$county_fips, sep="")


# Map color scale
breaks <- c(0, 2, 4, 6, 8, 10, 12, 14, 100)
pal <- colorRampPalette(c("#b9d7ef", "#102d44"))
col <- pal(length(breaks)-1)

years <- unique(unemp$year)

#
# Use Animation package
#

library(animation)


# Animation options
ani.options(outdir = paste(getwd(), "/images", sep=""), interval=0.4, ani.width=700, ani.height=420)

# Animated GIF
saveGIF({
  for (yi in 1:length(years)) {
    
    unemp_curr <- unemp[unemp$year == years[yi],]
    centroid_join <- match(unemp_curr$full_fips, as.vector(centroids@data$GEOID))
    unemp_curr[,c("x", "y")] <- centroids@coords[centroid_join,]
    
    # Remove missing counties.
    unemp_curr <- na.omit(unemp_curr)
    unemp_curr$unemp_rate <- as.numeric(unemp_curr$unemp_rate)
    
    # Interpolate and filter.
    unemp_curr_smooth <- interp(unemp_curr$x, unemp_curr$y, unemp_curr$unemp_rate, xo = xo, yo = yo, extrap = TRUE, linear=FALSE)
    unemp_curr_smooth$z[which(usa_lower_mat == FALSE)] <- NA
    
    
    # Not on first year in series.
    if (yi != 1) {
      
      frames_per_year <- 2
      frames_so_far <- 0
      
      one_tick <- (unemp_curr_smooth[["z"]] - unemp_prev_smooth[["z"]]) / frames_per_year
      ticker <- unemp_prev_smooth
      
      for (j in 1:frames_per_year) {
        
        # Increment one tick
        ticker[["z"]] <- ticker[["z"]] + one_tick
        
        par(mar=c(1,1,1,1))
        image(ticker, col=col, breaks=breaks, asp=1, axes=FALSE, useRaster=TRUE)
        plot(states_albers, add=TRUE, border="#ffffff", lwd=0.4)
        plot(lakes_albers, add=TRUE, border=NA, col="#ffffff")
        
        par(list(new=TRUE, plt=c(.015, .285, .015, .260)))
        plot(0, 0, type="n", axes=FALSE, xlim=c(0,10), ylim=c(0,10))
        
        # Year label
        text(5, 6, labels = years[yi], family="sans", cex=3, pos=3)
        
        # Legend limits
        legend_xlim <- c(0, 10)
        legend_ylim <- c(0, 10)
        legend_height <- 1
        legend_cat_width <- (legend_xlim[2] - legend_xlim[1]) / length(col)
        
        # Color scale
        text(5, 3, "Unemployment Rate (%)", pos=3, cex=0.9)
        for (k in 1:length(col)) {
          rect((k-1)*legend_cat_width, 2, k*legend_cat_width, legend_height+2, border="black", col=col[k], lwd=0.6)
        }
        text(0:length(col)*legend_cat_width, 1.2, labels=breaks, family="sans", cex=0.75)
        
        
        frames_so_far <- frames_so_far + 1
      }

    }
  unemp_prev_smooth <- unemp_curr_smooth
  }
}, movie.name="US_unemployment_rate.gif")