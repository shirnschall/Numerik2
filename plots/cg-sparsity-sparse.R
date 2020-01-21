library(ggplot2)
require(scales)
require(dplyr)

data=read.csv("C:\\Users\\shirnschall\\Documents\\GitHub\\Numerik2\\plots\\cg-sparsity-sparse",header = TRUE ,sep = "\t")

p <- ggplot(data,aes(x=n,y=t,color=factor(density),group=factor(density)))+
  geom_point(aes(shape = factor(density))) + 
  #geom_path(aes(group = factor(density)))+
  geom_smooth()+ # argument se=F schaltet konvidenzintervall aus
  theme_bw() +
  labs(color = "Einträge ungleich null",group="Einträge ungleich null",linetype="Einträge ungleich null",shape="Einträge ungleich null")+
  theme(
    legend.position = c(.03, .97),
    legend.justification = c("left", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )+
scale_y_log10()+
  ylab("Zeit [us]") +
  xlab("Matrix (nxn)")

p


ggsave("cg-sparsity-sparse.png", units="in", width=5, height=4, dpi=300)