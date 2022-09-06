## part6. compare all methods of prognostic value of cell-subtype

## library
library(ggplot2)
library(grid)
library(ggthemes)

## load data
data_df <- read.table(file = "./compare.txt", sep = "\t", header = T, fileEncoding = "utf-8")
colnames(data_df) <- c("Cell_type", "Methods", "Accuracy", "Flase rate", "Contradiction rate")
head(data_df)

## plot
for (i in 1:(ncol(data_df) - 2)) {
    data_df_tem <- data_df
    colnames(data_df_tem)[i + 2] <- "y"

    p <- ggplot(data_df_tem, aes(x = Methods, y = y, fill = Cell_type, label = y)) +
        geom_bar(stat = "identity", position = "dodge") +
        theme_wsj(color = "white") +
        scale_fill_wsj("colors6", "") +
        facet_grid(Cell_type ~ .) +
        coord_flip() +
        guides(fill = guide_legend(title = NULL)) +
        ggtitle(paste0("The ", colnames(data_df)[i + 2], " of all methods")) +
        geom_text(aes(y = y + 0.05), position = position_dodge(0.9), hjust = 1.1, colour = "black", size = 3.5) +
        theme(
            axis.title = element_blank(),
            legend.position = "none",
            panel.grid.major.x = element_line(linetype = "dashed", colour = "grey60"),
            panel.grid.major.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(),
            axis.ticks.length = unit(0.25, "cm"),
            axis.line.x = element_blank(),
            axis.line.y = element_line(),
            axis.text.y = element_text(size = 10),
            axis.text.x = element_text(size = 10)
        )

    pdf(paste0("./figure/Compare all methods via ", colnames(data_df)[i + 2], ".pdf"), width = 12, height = 8)
    print(p)
    dev.off()
}
