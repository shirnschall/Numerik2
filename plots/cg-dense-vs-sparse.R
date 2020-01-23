library(ggplot2)
require(scales)
require(dplyr)

data=read.csv("/Users/shirnschall/Desktop/Numerik2/plots/cg-dense-vs-sparse",header = TRUE ,sep = "\t")


#vergleichsfunktionen
n <- seq(from=0.1,to=1250,by=0.1)
f <- function(a){
  a*a
}
g <- function(a){
  a
}
t<-c(f(n),g(n))
type<-c(rep("x*x",times=length(n)), 
        rep("x",times=length(n)))
density<-c(rep("n",times=length(n)), 
           rep("1",times=length(n)))
n<-c(n,n)
d = data.frame(n,t,type,density)





p <- ggplot(data,aes(x=n,y=t,color=type,group=type))+
  geom_point(aes(shape = type)) + 
  #geom_path(aes(group = type))+
  geom_smooth()+ # argument se=F schaltet konvidenzintervall aus
  theme_bw() +
  labs(color = "Art der Matrix",group="Art der Matrix",linetype="Art der Matrix",shape="Art der Matrix")+
  theme(
    legend.position = c(.97, .03),
    legend.justification = c("right", "bottom"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )+
scale_y_log10()+
  ylab("Zeit [\u03bcs]") +
  xlab("Matrix (n\u00d7n)")+
  #vergleichsfunktionen
  geom_line(data = d, aes(x=n, y=t, group=density, colour=density),
            show_guide = FALSE)


p


#ggsave("cg-dense-vs-sparse.png", units="in", width=5, height=4, dpi=300)