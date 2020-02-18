
PGCExpr <- GeneExprPan("ENSG00000096088.15", cancerid)

PGCSumrise <- PGCExpr %>% mutate(logPGC = log10(ENSG00000096088.15 + 1)) %>% 
    group_by(CacnerType, Group) %>% summarise(mean = mean(logPGC), sd = sd(logPGC),
                                              median = median(logPGC)) %>% 
    arrange(CacnerType)
write.csv(PGCSumrise, "PGCExpr.csv", row.names = F, quote = F)
ggplot(PGCExpr, aes(CacnerType, log10(ENSG00000096088.15 + 1), color = Group)) + 
    geom_boxplot() + 
    labs(x = NULL, y = "log10(PGC + 1)") + theme_bw() + 
    theme(axis.text.x = element_text(angle = 60, hjust = 1),
          legend.position = "top")
ggsave("PGCPan.pdf", width = 8, height = 4)