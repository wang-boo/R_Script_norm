library(ggplot2)
library(openxlsx)



setwd("D:/WD_1711125/combined/txt")
dd1 <- read.xlsx("KEGG.xlsx")



bar_plot <- function(df1){
    nrow_df <- nrow(df1)
    # limit the number of iterms
    iterm_count <- function(dataframe, MAX = 10){
        if(nrow_df == 0){
            itermcount = NULL
        }else{
            itermcount <- 1:min(MAX, nrow_df)
            return(dataframe[itermcount,])
        }
    }
    cols <- brewer.pal(8, "Paired")[2:4]


    dd1 <- iterm_count(df1)
    dd1$Description <- factor(dd1$Description, levels = rev(x = dd1$Description))


    fg <- ggplot(dd1, aes(x = Description, y = -log10(p.adjust))) + 
        coord_flip() +
        geom_bar(
            fill = cols[1], 
            stat = "identity",
            width = 0.7,
            size = 0.25,
            alpha = 0.8
            ) +
        ggtitle("Enrich KEGG") + 
        xlab("")+
        geom_text(
            aes(label = paste("[", Count, " genes]", sep=" ")),
            nudge_y = 2,
            size = 3
            )
        ggsave(
            dpi = 300,
            units = "cm",
            scale = 0.7,
            paste("KEGG","tiff",sep="."),
            plot = fg,
            width = 35,
            height = 10 + nrow_df * 0.5
            )
        ggsave(
            dpi = 300,
            units = "cm",
            scale = 0.7,
            paste("KEGG","pdf",sep="."),
            plot = fg,
            width = 35,
            height = 10 + nrow_df * 0.5
        )
}
