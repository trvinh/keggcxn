library(ggplot2)

benchmark <- as.data.frame(read.table("/Users/trvinh/work/R_projects/keggcxn/data/benchmark.txt", sep='\t',header=TRUE,check.names=FALSE,comment.char="",stringsAsFactors=FALSE))

p <- ggplot(benchmark,aes(x=Time)) +  
  geom_point(aes(y=Nodes, colour="nodes")) +
  geom_smooth(aes(y=Nodes, colour="nodes"), method = "lm", se = FALSE, size=.5) +
  geom_point(aes(y=Edges, colour="edges")) +
  geom_smooth(aes(y=Edges, colour="edges"), method = "lm", se = FALSE, size=.5) +
  geom_point(aes(y=Degree*100, colour="degree")) +
  geom_smooth(aes(y=Degree*100, colour="degree"), method = "lm", se = FALSE, size=.5) +
  scale_color_manual(name = "", values = c(nodes = "#ffcc5c", edges = "#ff6f69", degree = "#96ceb4")) +
  theme_minimal()
p <- p + xlab("Time (seconds)") +
  ylab("Count") +
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15))
p