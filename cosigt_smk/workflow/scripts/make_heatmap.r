library(data.table)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

df<-data.frame(fread(args[1]))
selected_paths<-df$path.name
df<-df[,c(2:ncol(df))]
df_out<-matrix(0, nrow=nrow(df)*ncol(df), ncol=3)

n<-0
for (i in c(1:nrow(df))) {

    print(i)
    for (j in c(1:ncol(df))) {

        n<-n+1
        df_out[n,1]<- selected_paths[i]
        df_out[n,2]<- j
        df_out[n,3]<- df[i,j]

    }

}
df<-data.frame(df_out)
colnames(df)<-c('x', 'y', 'z')

df$y<-as.numeric(df$y)
df$z<-as.numeric(df$z)


p_all <- ggplot(df, aes(x = y, y = x, fill = z)) +
  geom_raster() +
  scale_fill_gradient2(
      oob = scales::squish
  ) +
  theme_void() +
  theme(axis.text.x = element_text(), axis.title.x=element_text())+
  theme(axis.text.y = element_text(size=5))+
  theme(legend.position = "none") +
  labs(x="node #")

ggsave("all.pdf", width=15, height=12)

selected_paths<-c("chr1_hg19_103998686_104406594_0", "chr1_hg38_103456064_103863972_0", "chr1_chm13_103304997_103901127_0", "HG01175#1#JAHAMA010000003.1_33204707_33808232_1", "HG02109#1#JAHEPG010000026.1_16673035_17080816_0", "HG002#2#JAHKSD010000021.1_102757493_103165350_0", "HG02257#2#JAGYVH010000034.1_98190388_98857961_0")
df<-subset(df, x%in%selected_paths)
#df$x<-factor(df$x,levels=selected_paths)
df$y<-as.numeric(df$y)
df$z<-as.numeric(df$z)

p <- ggplot(df, aes(x = y, y = x, fill = z)) +
  geom_raster() +
  scale_fill_gradient2(
      oob = scales::squish
  ) +
  theme_void() +
  theme(axis.text.x = element_text(), axis.title.x=element_text())+
  theme(axis.text.y = element_text())+
  theme(legend.position = "none") +
  labs(x="node #")

ggsave("zoom1.pdf", width=15, height=2)

#zoom on selected nodes
#example_nodes
nodes<-c(9800:9850)
df_sub<-subset(df, y%in%nodes)

p_zoom <- ggplot(df_sub, aes(x = y, y = x, fill = z)) +
  geom_tile(color = "white",
            lwd = .5,
            linetype = 1) +
  geom_text(aes(label = z), color = "white", size = 4) +
  scale_fill_gradient2(
      oob = scales::squish
  ) +
  theme_void() +
  theme(axis.text.x = element_text(), axis.title.x=element_text())+
  theme(axis.text.y = element_text(size=5))+
  theme(legend.position = "none") +
  labs(x="node #")

ggsave("zoom2.pdf", width=15, height=2)
