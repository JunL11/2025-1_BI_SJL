##CLIP-seq data processing
#awk '$3 == "exon" && $0 ~ /gene_type "miRNA"/' gencode.gtf > mirna_only.gtf

#featureCounts -a mirna_only.gtf \
#-o clip_miRNA_counts.txt \
#-T 4 \
#-g gene_name \
#-t exon \
#CLIP-35L33G.bam


library(dplyr)

process_mirna_file <- function(file_path) {
  df <- read.table(file_path,
                   sep = "\t",
                   header = TRUE,
                   skip = 9,
                   quote = "",
                   fill = TRUE,
                   stringsAsFactors = FALSE)
  
#QC
  df_filtered <- df %>%
    filter(IsManualFlag == 0, gIsPosAndSignif == 1)
  
#필요한 열만 선택
  mirna_expr <- df_filtered[, c("GeneName", "gProcessedSignal")]
  colnames(mirna_expr) <- c("miRNA", "Signal")
  
#같은 miRNA에서 probe 평균
  mirna_avg <- mirna_expr %>%
    group_by(miRNA) %>%
    summarise(Signal = mean(Signal, na.rm = TRUE))
  
#별표 제거 (e.g., mmu-let-7d*)
  mirna_avg$base_miRNA <- gsub("\\*", "", mirna_avg$miRNA)
  
#base 이름 기준으로 signal 합산
  mirna_final <- mirna_avg %>%
    group_by(base_miRNA) %>%
    summarise(Signal = sum(Signal, na.rm = TRUE))
  
  return(mirna_final)
}

# 파일별 처리
gfp1 <- process_mirna_file("miRNA_siGFP1.txt")
gfp2 <- process_mirna_file("miRNA_siGFP2.txt")
lin1 <- process_mirna_file("miRNA_siLin28a1.txt")
lin2 <- process_mirna_file("miRNA_siLin28a2.txt")

# 병합
library(purrr)
library(tidyr)

all_data <- list(gfp1, gfp2, lin1, lin2) %>%
  purrr::reduce(full_join, by = "base_miRNA")

# 컬럼 이름 지정
colnames(all_data) <- c("miRNA", "siGFP1", "siGFP2", "siLin28a1", "siLin28a2")

# 결측값 처리
all_data <- na.omit(all_data)

##로그화
log_data <- log_data %>%
  mutate(
    log2_siGFP1     = log2(siGFP1 + 1),
    log2_siGFP2     = log2(siGFP2 + 1),
    log2_siLin28a1  = log2(siLin28a1 + 1),
    log2_siLin28a2  = log2(siLin28a2 + 1)
  )

#정규화
log2_mat <- as.matrix(log_data[, c("log2_siGFP1", "log2_siGFP2", "log2_siLin28a1", "log2_siLin28a2")])

# 정규화 (quantile normalization)
library(limma)
norm_mat <- normalizeBetweenArrays(log2_mat, method = "quantile")

# 정규화된 값을 다시 log_data에 넣기
log_data$log2_siGFP1     <- norm_mat[, "log2_siGFP1"]
log_data$log2_siGFP2     <- norm_mat[, "log2_siGFP2"]
log_data$log2_siLin28a1  <- norm_mat[, "log2_siLin28a1"]
log_data$log2_siLin28a2  <- norm_mat[, "log2_siLin28a2"]


# 2. 그룹별 평균 계산 (log2 평균)

log_data <- log_data %>%
  mutate(
    mean_log2_GFP = rowMeans(dplyr::select(., log2_siGFP1, log2_siGFP2), na.rm = TRUE),
    mean_log2_Lin28a = rowMeans(dplyr::select(., log2_siLin28a1, log2_siLin28a2), na.rm = TRUE)
  )

# 3. log2 Fold Change 계산
log_data <- log_data %>%
  mutate(
    log2FC = mean_log2_Lin28a - mean_log2_GFP
  )

##################norm_data 이름 형식 clip_data처럼 바꾸기
log_data$miRNA_simple <- gsub("miR-", "Mir", log_data$miRNA)
log_data$miRNA_simple <- gsub("^mmu-", "", log_data$miRNA_simple)
log_data$miRNA_simple <- ifelse(
  grepl("let-", log_data$miRNA_simple),
  paste0("Mir", log_data$miRNA_simple),
  log_data$miRNA_simple
)
log_data$miRNA_simple <- gsub("let-", "let", log_data$miRNA_simple)

##################Variation 큰 애들 cutoff
log_data$sd_GFP <- apply(log_data[, c("log2_siGFP1", "log2_siGFP2")], 1, sd, na.rm = TRUE)
log_data$sd_Lin28a <- apply(log_data[, c("log2_siLin28a1", "log2_siLin28a2")], 1, sd, na.rm = TRUE)

# 특정 기준 이상이면 플래그 후 제거
log_data$high_variation <- (log_data$sd_GFP > 1 | log_data$sd_Lin28a > 1)
log_data_filtered <- log_data[log_data$high_variation == FALSE, ]


# 1. base miRNA 이름 생성 (-3p/-5p 이름 통합)
log_data_filtered$base_miRNA <- gsub("-[35]p$", "", log_data_filtered$miRNA_simple)

# 2. base_miRNA 기준으로 샘플별 합산 → log2 변환
collapsed_log_data <- log_data_filtered %>%
  group_by(base_miRNA) %>%
  summarise(
    siGFP1     = sum(siGFP1, na.rm = TRUE),
    siGFP2     = sum(siGFP2, na.rm = TRUE),
    siLin28a1  = sum(siLin28a1, na.rm = TRUE),
    siLin28a2  = sum(siLin28a2, na.rm = TRUE)
  ) %>%
  mutate(
    log2_siGFP1    = log2(siGFP1 + 1),
    log2_siGFP2    = log2(siGFP2 + 1),
    log2_Lin28a1   = log2(siLin28a1 + 1),
    log2_Lin28a2   = log2(siLin28a2 + 1)
  )

collapsed_log_data <- collapsed_log_data %>%
  mutate(
    mean_log2_GFP     = (log2_siGFP1 + log2_siGFP2) / 2,
    mean_log2_Lin28a  = (log2_Lin28a1 + log2_Lin28a2) / 2,
    log2FC            = mean_log2_Lin28a - mean_log2_GFP
  )

colnames(collapsed_log_data)[colnames(collapsed_log_data) == "base_miRNA"] <- "miRNA"

#######################figure 위해 a =a-1, f = f-1, c = c-2로 변환

collapsed_log_data <- collapsed_log_data %>%  mutate(miRNA = case_when(
    miRNA == "Mirlet7a" ~ "Mirlet7a-1",
    miRNA == "Mirlet7f" ~ "Mirlet7f-1",
    miRNA == "Mirlet7c" ~ "Mirlet7c-2",
    TRUE ~ miRNA
  ))

##############clip_count 합치기
collapsed_log_data <- collapsed_log_data %>%
  left_join(clip_data %>% dplyr::select(miRNA, clip_count), by = "miRNA")


###############clip_count NA나 0으로 나오는거 삭제
collapsed_log_data <- collapsed_log_data %>%
  filter(!is.na(clip_count) & clip_count != 0)

###################x축 위해 clip count/mean_siGFP (log2)
collapsed_log_data <- collapsed_log_data %>%
  mutate(clip_norm = log2(clip_count / mean_log2_GFP))

#################그래프 그리기
library(ggplot2)
library(dplyr)

# 관심 있는 miRNA (let-family 및 mir98)
highlight_miRNAs <- collapsed_log_data %>%
  filter(grepl("let", miRNA, ignore.case = TRUE) | grepl("mir98", miRNA, ignore.case = TRUE))

# 파란색 강조할 miRNA 목록
blue_miRNAs <- c("Mir703", "Mir682", "Mir677", "Mir708", "Mir344")

highlight_blue <- collapsed_log_data %>%
  filter(miRNA %in% blue_miRNAs) %>%
  mutate(color = "blue")

# 나머지 (배경 점들)
background_miRNAs <- collapsed_log_data %>%
  filter(!(miRNA %in% highlight_miRNAs$miRNA))

# 그래프 그리기
ggplot() +
  geom_point(data = background_miRNAs,
             aes(x = clip_norm, y = log2FC),
             color = "black", size = 2) +
  
  geom_point(data = highlight_miRNAs,
             aes(x = clip_norm, y = log2FC),
             color = "red", size = 2) +
  geom_text(data = highlight_miRNAs,
            aes(x = clip_norm, y = log2FC, label = miRNA),
            vjust = -0.5, hjust = 0.5, color = "red", size = 3) +
  
  geom_point(data = highlight_blue,
             aes(x = clip_norm, y = log2FC),
             color = "blue", size = 2) +
  geom_text(data = highlight_blue,
            aes(x = clip_norm, y = log2FC, label = miRNA),
            vjust = -0.5, hjust = 0.5, color = "blue", size = 3) +
  
  scale_x_continuous(limits = c(-5, 10)) +  # 👈 x축 범위 지정
  scale_y_continuous(limits = c(-1, 3)) +  # 👈 x축 범위 지정
  
  labs(
    title = "LIN28A binding and miRNA regulation",
    x = "LIN28A CLIP tag count relative to mature miRNA amount (log2)",
    y = "Mature microRNA change upon Lin28a knockdown (log2)"
  ) +
  theme_minimal()
