df3 <- read.table(
  "miRNA_siGFP1.txt",
  sep = "\t",
  header = TRUE,
  skip = 9,
  quote = "",
  fill = TRUE,
  stringsAsFactors = FALSE
)

# GeneName에 "let-7a"가 포함된 행 전체 필터링
let7c_rows <- df3[grepl("let-7c", df$GeneName), ]

# 결과 출력
print(let7c_rows)
