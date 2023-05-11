
micro <- read.csv("/Users/ianbrennan/Google.Drive/ANU/AHE/T428_Microhylidae/Microhylidae_Range_Elevation.csv")
    micro$logRange <- log(micro$Range.size)



ggplot(micro, aes(x=Range.size, y=Elevation_min)) +
  geom_point(aes(color=EPBC), size=3) + coord_trans(x="log") +
  theme_bw() + theme(legend.position = "below")

library(plotly)

plot_ly(micro,
        x = ~Range.size,
        y = ~Elevation_min,
        text = ~paste("Species:", Species),
        type = 'scatter',
        mode = 'markers',
        color = ~EPBC,
        marker = list(size = 10,
                      line = list(width = 1,
                                  color = "black"))) %>%
  layout(xaxis = list(type = "log"))


mtree <- read.tree("/Users/ianbrennan/Google.Drive/ANU/AHE/T428_Microhylidae/Australian_Microhylids_Elevation.tre")
plot(mtree); axisPhylo()


elev <- as.matrix(micro[,c( "logRange", "Elevation_min")])
rownames(elev) <- micro$Species
phylomorphospace(mtree, elev, label="horizontal")
