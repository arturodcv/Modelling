library(corrplot)
library(gplots)
library(ggplot2)
library(MASS)
library(ggmap)
library(maps)
library(dplyr)
library(GGally)
library(geosphere)

##################################
##################DATA PREPARATION
##################################
load(file = "data.Rda")
load(file = "data_non_aftershocks.Rda")

summary(data)




################################
###################DATA CLUSTERS
################################

num_data = data[,-c(1,2)]
mag = data$mag
mag_data = num_data[,c(4,11,12)]
depth_data = num_data[,c(3,10,5)]
dist_data = num_data[,c(11,6)]

################################
#############################MAP
################################

register_google(key = "AIzaSyC4TejqvRFWKSGaTIJFsFrGTOpbqeb1AdQ")
japan.map = get_map(location = c(lon = 138.077031, lat = 35.801589 ),
                    maptype = "satellite", zoom = 5)

data_from_2011 = data[5000:11780,]
ggmap(japan.map, extent = "device")+
  geom_point( data = data_from_2011,aes(x=longitude,y = latitude),
              colour = "yellow", alpha = 0.1)+
  scale_alpha(range = c(0.1,0.4))

################################
####################CORRELATIONS
################################

pairs(mag_data)
pairs(depth_data)
pairs(dist_data)

M = cor(na.omit(mag_data, na.action = "omit", fill = NULL))
M = cor(na.omit(depth_data, na.action = "omit", fill = NULL))
M = cor(na.omit(dist_data, na.action = "omit", fill = NULL))

corrplot.mixed(M)

################################
#######################MAGNITUDE
################################

ts_mag = ts(as.data.frame(mag))

plot(data.frame(data$time_series,data$mag)[1:80,],type = "b", 
     main = "Magnitude among time in events", ylab = "Magnitude", xlab = "Time (hours)")
hist_mag = hist(mag,breaks = 47 , plot=T, xlim = c(4,9), 
                xlab = "Magnitude", main = "Histrogram of magnitude")
hist_fit = hist(data$mag[data$mag <= 7.5],breaks = 31 , plot=T)

#Logarithmic fit
hist_fit$counts = log(1+hist_fit$counts)
yy <- hist_fit$counts[1:50]
xx <- hist_fit$mids[1:50]
fit <- lm(yy ~ xx)
fitround <- round(coef(fit),3)
eq <- paste("y=",fitround[1],"+",fitround[2],"t")
eq = paste("y = 19.01 - 2.56")
print(fit$coefficients)
hist_mag$counts = log(1+hist_mag$counts)

#Plot log histogram
plot(hist_mag, xlim=c(4.4,9.1),main="Histogram of magnitude", xlab = "Magnitude")
lines(hist_fit$mids, fit$coefficients[1]+hist_fit$mids*fit$coefficients[2])
mtext(eq,3,line=-5, cex = 1.5)

#Autocorrelation function for magnitude
corr <- acf(ts_mag, type = "correlation", plot = T,
            main = "Autocorrelation values for magnitude")
plot(corr, ylim=c(0.001,1), log="y", main = "Autocorrelation values for magnitude")

#Plot periodogram for magnitude
a = spec.pgram(ts_mag, kernel("daniell",m=5),
               plot = T, main = "Series: Magnitude\nSmoothed Periodogram")

#Autocorrelation function
corr <- acf(data$time_btw_events, type = "correlation")

#Plot acf
plot(corr, ylim=c(0.001,1), log="y", 
     main = "Autocorrelation values for time between events")

spec.pgram(data$time_btw_events, kernel("daniell",m=5),
           plot = T, main = "Series: Time between events\nSmoothed Periodogram")


########################################
#############################AFTERSHOCKS
########################################

#Haversine distance given latitude and longitude of two points
haversine = function(lat1,long1,lat2,long2){
  lat1 = lat1*pi/180
  lat2 = lat2*pi/180
  long1 = long1*pi/180
  long2 = long2*pi/180
  
  R = 6371
  d = 2*R*asin( sqrt(   sin(0.5*(lat2-lat1))^2+cos(lat1)*cos(lat2)*sin(0.5*( long2-long1))^2 ) )
  d
}

#Remove aftershocks iteratively
data_no_aftershocks = function(data){
  data = data[,c(3,4,6,15,16)]
  data$time_btw_events = data$time_btw_events/24
  a = c()
  for(k in 1:(length(data$latitude)-1)){
    if (data$mag[k] >  data$mag[k+1]){
      if (data$mag[k] < 5){
        if ((data$time_series[k+1]-data$time_series[k]) <= 83 || haversine(data$latitude[k],data$longitude[k],data$latitude[k+1],data$longitude[k+1]) < 35){
          a = c(a,k+1)}}
      else if(data$mag[k] > 5 && data$mag[k] < 5.5){
        if ((data$time_series[k+1]-data$time_series[k]) <= 155 || haversine(data$latitude[k],data$longitude[k],data$latitude[k+1],data$longitude[k+1])< 40){
          a = c(a,k+1)
        } 
      }
      else if(data$mag[k] > 5.5 && data$mag[k] < 6){
        if ((data$time_series[k+1]-data$time_series[k]) <= 290 || haversine(data$latitude[k],data$longitude[k],data$latitude[k+1],data$longitude[k+1])< 47){
          a = c(a,k+1)
        } 
      }
      else if(data$mag[k] > 6 && data$mag[k] < 6.5){
        if ((data$time_series[k+1]-data$time_series[k]) <= 510 || haversine(data$latitude[k],data$longitude[k],data$latitude[k+1],data$longitude[k+1]) < 54){
          a = c(a,k+1)
        } 
      }
      else if(data$mag[k] > 6.5 && data$mag[k] < 7){
        if ((data$time_series[k+1]-data$time_series[k]) <= 790 || haversine(data$latitude[k],data$longitude[k],data$latitude[k+1],data$longitude[k+1]) < 61){
          a = c(a,k+1)
        } 
      }
      else if(data$mag[k] > 7 && data$mag[k] < 7.5){
        if ((data$time_series[k+1]-data$time_series[k])<= 915 || haversine(data$latitude[k],data$longitude[k],data$latitude[k+1],data$longitude[k+1]) < 70){
          a = c(a,k+1)
        } 
      }
      else if(data$mag[k] > 7.5 && data$mag[k] < 8){
        if ((data$time_series[k+1]-data$time_series[k]) <= 960 || haversine(data$latitude[k],data$longitude[k],data$latitude[k+1],data$longitude[k+1])< 81){
          a = c(a,k+1)
        } 
      }
      else if(data$mag[k] > 8){
        if ((data$time_series[k+1]-data$time_series[k]) <= 985 || haversine(data$latitude[k],data$longitude[k],data$latitude[k+1],data$longitude[k+1]) < 94){
          a = c(a,k+1)
        } 
      }
    }}
  a
} 
iter_aftershocks = function(data){
  b  = data_no_aftershocks(data)
  while (!is.null(b)) {
    data = data[-b,]
    rownames(data) <- 1:nrow(data)
    b = data_no_aftershocks(data)
  }
  data
}

data_non_aftershocks = iter_aftershocks(data)


#Autocorrelation function for magnitude
corr <- acf(data_non_aftershocks$mag, type = "correlation",
            main = "Autocorrelation values for magnitude with no aftershocks")
plot(corr, ylim=c(0.001,1), log="y", main = "Autocorrelation values for magnitude")

#Plot periodogram for magnitude
a = spec.pgram(data_non_aftershocks$mag, kernel("daniell",m=3),
               plot = T, main = "Series: Magnitude for non aftershocks earthquakes\nSmoothed Periodogram")

#Autocorrelation function for time between events
corr <- acf(data_non_aftershocks$time_btw_events, type = "correlation")
plot(corr, ylim=c(0.001,1), log="y", 
     main = "Autocorrelation values for time between events with no aftershocks")

#Plot periodogram for time between events
spec.pgram(data_non_aftershocks$time_btw_events, kernel("daniell",m=3),
           plot = T, main = "Series: Time between events for non aftershocks\nSmoothed Periodogram")









