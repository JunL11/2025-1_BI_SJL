# 1. pandas → read_csv에 해당: 파일 불러오기
pileup <- read.table("CLIP-let7d-gene.pileup",
                     sep = "\t", header = FALSE,
                     col.names = c("chrom", "pos", "_ref", "count", "basereads", "quals"),
                     stringsAsFactors = FALSE)

# 2. 정규표현식 제거: basereads에서 < > $ * ^ 문자 제거하여 matches 열 생성
pileup$matches <- gsub("[<>\\$\\*\\^]", "", pileup$basereads)

# 3. 확인: 필요한 컬럼만 보기
pileup_subset <- pileup[, c("chrom", "pos", "matches")]
print(head(pileup_subset))  # 또는 View(pileup_subset) in RStudio

# A/T/G/C만 남기고 모두 제거 (정규표현식 사용)
pileup_subset$matches <- gsub("[^ATGC]", "", pileup_subset$matches)

# 필요한 패키지
library(dplyr)
library(stringr)

# 1. 파일 불러오기
pileup <- read.table("CLIP-let7d-gene.pileup", sep = "\t", header = FALSE,
                     col.names = c("chrom", "pos", "_ref", "count", "basereads", "quals"),
                     stringsAsFactors = FALSE)

# 2. basereads에서 분석에 불필요한 문자 제거 → matches 열 생성
pileup <- pileup %>%
  mutate(matches = gsub("[^ATGC]", "", basereads))  # A/T/G/C만 남기기

# 3. Shannon entropy 계산 함수 정의
shannon_entropy <- function(seq) {
  bases <- str_split(seq, "")[[1]]
  total <- length(bases)
  if (total == 0) return(0)
  probs <- table(bases) / total
  -sum(probs * log2(probs))
}

# 4. 각 position별 Shannon entropy 계산
pileup <- pileup %>%
  rowwise() %>%
  mutate(entropy = shannon_entropy(matches)) %>%
  ungroup()

# 5. bedGraph 형식으로 변환: chrom, start, end, entropy
pileup <- pileup %>%
  mutate(start = pos, end = pos + 1)

bedgraph <- dplyr::select(pileup, chrom, start, end, entropy)

# 6. 파일로 저장
write.table(bedgraph, "shannon_entropy.bedgraph_let7d", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
