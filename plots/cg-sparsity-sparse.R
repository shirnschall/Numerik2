library(ggplot2)
require(scales)
require(dplyr)


data=read.csv("/Users/shirnschall/Desktop/Numerik2/plots/cg-sparsity-sparse",header = TRUE ,sep = "\t")
#data=read.csv("C:\\Users\\shirnschall\\Documents\\GitHub\\Numerik2\\plots\\cg-sparsity-sparse",header = TRUE ,sep = "\t")

p <- ggplot(data,aes(x=n,y=t,color=factor(density),group=factor(density)))+
  geom_point(aes(shape = factor(density))) + 
  #geom_path(aes(group = factor(density)))+
  geom_smooth()+ # argument se=F schaltet konvidenzintervall aus
  theme_bw() +
  #umlaut a = \u00e4
  labs(color = "Eintr\u00e4ge ungleich null",group="Eintr\u00e4ge ungleich null",linetype="Eintr\u00e4ge ungleich null",shape="Eintr\u00e4ge ungleich null")+
  theme(
    legend.position = c(.97, .03),
    legend.justification = c("right", "bottom"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )+
scale_y_log10()+
  ylab("Zeit [\u03bcs]") +
  xlab("Matrix (n\u00d7n)")

p


ggsave("cg-sparsity-sparse.png", units="in", width=5, height=4, dpi=300)