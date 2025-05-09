setwd("/home/sangjun.lim/2025bi/binfo1-work")

library(dplyr)
library(ggplot2)

# 데이터 불러오기 및 처리 함수
process_file <- function(filename, sample_name) {
  df <- read.table(filename, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  colnames(df) <- c("chr", "read_start", "read_end", "count", "exon_chr", "exon_start", "exon_end", "transcript_id", "start_codon", "strand")
  
  df$relative_position <- ifelse(
    df$strand == "+",
    df$read_start - df$start_codon,
    df$start_codon - df$read_start
  )
  
  df_filtered <- subset(df, relative_position >= -50 & relative_position <= 50)
  
  df_counts <- df_filtered %>%
    group_by(relative_position) %>%
    summarise(total_count = sum(count)) %>%
    mutate(sample = sample_name)
  
  return(df_counts)
}

# 파일 처리
df_siluc <- process_file("fivepcounts-filtered-RPF-siLuc.txt", "siLuc")

# plot
p <- ggplot(df_siluc, aes(x=relative_position, y=total_count / 1000)) +
  annotate("rect", xmin=-10, xmax=10, ymin=0, ymax=Inf, alpha=0.1, fill="orange") +
  geom_segment(x=0, xend=0, y=0, yend=100, data=NULL, color="red", linewidth=0.5) +
  coord_cartesian(ylim = c(0, 100)) +
  geom_hline(yintercept=0, color="black", linewidth=0.5) +
  geom_bar(stat="identity", width=0.6, fill="black") +
  theme_bw() +
  labs(
    title = "Ribosome footprint density near start codons",
    x = "Relative position to start codon of 5'-end of reads",
    y = "Raw read count (x1000)"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size=12),
    axis.title = element_text(size=14),
    plot.title = element_text(hjust=0.5, face="bold", size=14)
  )

ggsave("mission2.png", plot = p, width=8, height=4, dpi=300)
