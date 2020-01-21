library(ggplot2)
require(scales)
require(dplyr)

data=read.csv("/Users/shirnschall/Desktop/Numerik2/plots/cg-dense-vs-sparse",header = TRUE ,sep = "\t")

p <- ggplot(data,aes(x=n,y=t,color=type,group=type))+
  geom_point(aes(shape = type)) + 
  geom_path(aes(group = type))+
  theme_bw() +
  labs(color = "Typ",group="Typ",linetype="Typ",shape="Typ")+
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


ggsave("cg-dense-vs-sparse.png", units="in", width=5, height=4, dpi=300)