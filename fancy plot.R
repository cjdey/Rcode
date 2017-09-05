z_theme <- function() {
  library(RColorBrewer)
  # Generate the colors for the chart procedurally with RColorBrewer
  palette <- brewer.pal("Greys", n=9)
  color.background = "white"
  color.grid.major = palette[3]
  color.axis.text = palette[7]
  color.axis.title = palette[7]
  color.title = palette[8]
  # Begin construction of chart
  theme_bw(base_size=9) +
    # Set the entire chart region to a light gray color
    theme(panel.background=element_rect(fill=color.background, color=color.background)) +
    theme(plot.background=element_rect(fill=color.background, color=color.background)) +
    theme(panel.border=element_rect(color=color.background)) +
    # Format the grid
    theme(panel.grid.major=element_line(color=color.grid.major,size=.25)) +
    theme(panel.grid.minor=element_blank()) +
    theme(axis.ticks=element_blank()) +
    # Format the legend, but hide by default
    theme(legend.position="none") +
    theme(legend.background = element_rect(fill=color.background)) +
    theme(legend.text = element_text(size=7,color=color.axis.title)) +
    # Set title and axis labels, and format these and tick marks
    theme(plot.title=element_text(color=color.title, size=20, vjust=1.25)) +
    theme(axis.text.x=element_text(size=9,color=color.axis.text)) +
    theme(axis.text.y=element_text(size=9,color=color.axis.text, vjust = 0)) +
    theme(axis.title.x=element_text(size=12,color=color.axis.title)) +
    theme(axis.title.y=element_text(size=12,color=color.axis.title, vjust=1.25))
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

###Plot code

require(ggjoy)
require(viridis)

A<-ggplot(AVD, (aes(x=Research_effort+1, y=group2, fill = ..x.., color="white")))+
  geom_joy_gradient(scale =1.5, rel_min_height=0.01)+
  xlab("Research effort")+ 
  z_theme()+  
  scale_x_log10(
    breaks = c(1,10,100, 1000, 10000),
    labels = c(1, 10, 100, 1000, 10000)
    
  )+
  theme(axis.title.y = element_blank(), legend.position = "none",
        axis.title.x=element_text(hjust=0.38))+
  scale_y_discrete(expand = c(0.01, 0))+
  scale_fill_viridis(option="inferno", direction = -1)+
  scale_color_manual(values=rep("white", 7))

B<-ggplot(dfa)+
  geom_col(aes(g, y=val*100), fill ="grey")+
  z_theme()+
  scale_y_continuous(limits =c(0,105), breaks=c(0,25,50,75,100))+
  theme(legend.position="none")+
  coord_flip()+xlab("")+ylab("% IUCN evaluated")+
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_blank())

C<-ggplot(df)+
  geom_col(aes(y=(traitperc*100), x=group2), fill ="grey")+
  z_theme()+
  theme(legend.position="none")+
  scale_y_continuous(limits =c(0,105), breaks=c(0,25,50,75,100))+
  coord_flip()+xlab("")+ylab("% known traits")+
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_blank())



###MULTIPANEL PLOT
plots<-list(A,C,B)
layout <- matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE)
multiplot(plotlist = plots, layout = layout)
