library(ggplot2)
library(ggpmisc)
require(scales)
require(dplyr)


data=read.csv("/Users/shirnschall/Desktop/Numerik2/plots/precond-vs-cg-complexity",header = TRUE ,sep = "\t")
#data=read.csv("C:\\Users\\shirnschall\\Documents\\GitHub\\Numerik2\\plots\\cg-sparsity-sparse",header = TRUE ,sep = "\t")

#data=data[!(data$density==7),]
#data=data[!(data$density==5),]
data=data[(data$n>20),]
#data=data[!(data$density=='n'),]


#vergleichsfunktionen
lm_eqn = function(df){
  m = lm(t~poly(n,2,raw=TRUE),df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


n <- seq(from=60,to=2000,by=0.1)
f <- function(a){
  a*2
}
g <- function(a){
  a*2
}
t<-c(f(n),g(n))
type<-c(rep("x",times=length(n)), 
        rep("x",times=length(n)))
density<-c(rep("n",times=length(n)), 
           rep("1",times=length(n)))
n<-c(n,n)
d = data.frame(n,t,type,density)

#data = rbind(d,data)



scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
}


p <- ggplot(data,aes(x=n,y=t))+
  geom_point(aes(shape = factor(type),color = factor(type))) + 
  #geom_path(aes(group = factor(density)))+
  geom_smooth(aes(color=factor(type)),size=0.5,method="lm", se=FALSE, formula = y~poly(x,1,raw=TRUE))+ # argument se=F schaltet konvidenzintervall aus
  
  theme_bw() +
  #umlaut a = \u00e4
  labs(linetype="Vergleichsgerade", color = "Verfahren",shape="Verfahren")+
  theme(
    legend.position = c(.03, .97),
    legend.justification = c("left", "top"),
    legend.box.just = "top",
    legend.margin = margin(6, 6, 6, 6),
    legend.box = "horizontal"
  )+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  ylab("Zeit [\u03bcs] pro Iterationsschritt") +
  xlab("Matrixgr\u00f6\u00dfe (n)")+
  scale_color_discrete(labels = c("CG", "CG mit Vorkonditionierung"))+
  scale_shape_discrete(labels = c("CG", "CG mit Vorkonditionierung"))+
  #vergleichsfunktionen
  geom_line(data = d, aes(x=n, y=t,linetype=type))+
  scale_linetype_discrete(labels = c("\u039f(n)","\u039f(n\u00b2)"))





p


ggsave("precond-vs-cg-complexity.png", units="in", width=7, height=5, dpi=500)