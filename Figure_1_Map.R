## This code creates a map of the locations sampled.
## This code requires a .csv file of the latitude and longitude of the sampling site (Samplingsite.csv file).
## Due to poaching concerns, please contact clare.im.adams@gmail.com if you would like the location file.

library(maps)
library(mapdata)
library(mapproj)
library(maptools)
library(plotly)

map('nzHires', xlim=c(165,179), ylim=c(-50,-35))
nz<-map('nzHires')
visual <- ggplotly(nz, height = 1.2 * 600, width = 600, tooltip=c("text"), 
                   hoverinfo='hide', dynamicTicks = FALSE) 

samps <- read.csv("SamplingSite.csv", header=TRUE, stringsAsFactors=T)

#plot all city names of NZ onto
#map.cities(country="New Zealand", label=TRUE, cex=1,xlim=c(165,179), ylim=c(-50,-35), pch=20)

map('nzHires', xlim=c(164,179), ylim=c(-50,-35))
points(samps$Long, samps$Lat, pch=19, col="red",cex=1)
pointLabel(samps$Long, samps$Lat, samps$Sampling_site)



### 
##alternatively, ggmap

library(ggmap)

map.nz <- map_data(map = "nz")
p1 <- ggplot(map.nz, aes(x = long, y = lat, group=group))
p1 <- p1 + geom_polygon()
p1 <- p1
p1


p2=p1+geom_point(data = samps, aes(x = Longitude, y = Latitude,shape = Sampling_site, colour = Sampling_site), size = 7,inherit.aes = FALSE)+theme_bw()

p2

pal9<-c("#920000", "#009292")

p2+ geom_text(data = samps, aes(x = Longitude, y = Latitude, label = Sampling_site), hjust = -0.2, colour = "black", size = 3,inherit.aes = FALSE)+scale_color_manual(values=pal9)+xlab("Longitude")+ylab("Latitude")+theme(axis.title=element_text(size=15),legend.title=element_text(size=15))


#now with scalebar
p3<-p2+ geom_text(data = samps, aes(x = Longitude, y = Latitude, label = Sampling_site), hjust = -0.2, colour = "black", size = 3,inherit.aes = FALSE)+scale_color_manual(values=pal9)+xlab("Longitude")+ylab("Latitude")+theme(axis.title=element_text(size=15),legend.title=element_text(size=15))+ggsn::scalebar(map.nz, dist = 100, transform=TRUE, dist_unit = "km", model="WGS84", st.size=3)

# and this outputs the file into a PDF
pdf(file = "MAP_01.pdf")
plot(p3)
dev.off()



