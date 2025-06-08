#資料.output(base data) .model  genotype matrix(target data)
#test_data<-fread("test_data.txt")
#raw_model <- fread("single_model_2D.txt", fill = TRUE, header = FALSE)
#output_2D <- fread("single_output_2D.txt",
#                   header = FALSE,
#                   skip = 3,
#                   fill = FALSE,
#                   col.names = c("snp1", "snp2", "F_value", "p_value")) %>%
#  mutate(pair = paste(snp1, snp2, sep = "_"))
#load_all(".")
#setwd("D:/桌面/mrs_rewrite/lassosumMRS")
rm(list = ls())
library(data.table)
library(dplyr)
library(tibble)
library(lassosumMRS)
library(dplyr)

#等效devtools::load_all()
library(devtools)
load_all(".")




f_test_data <- system.file("extdata", "test_data.txt", package = "lassosumMRS")
f_model     <- system.file("extdata", "model_2D.txt", package = "lassosumMRS")
f_output    <- system.file("extdata", "output_2D.txt", package = "lassosumMRS")

test_data <- data.table::fread(f_test_data)
raw_model <- data.table::fread(f_model, fill = TRUE, header = FALSE)
output_2D <- data.table::fread(f_output,
                               header = FALSE,
                               skip = 3,
                               fill = FALSE,
                               col.names = c("snp1", "snp2", "F_value", "p_value")) %>%
  dplyr::mutate(pair = paste(snp1, snp2, sep = "_"))

#.model轉格式
model<- parse_model_to_mdr(raw_model)

#LASSO soft-thresholding lambda預設是預設為 log scale 的 20 個值
beta_obj <- mrs_indeplasso(output_2D)

#snp_data 需要是row是樣本 column是snp id
#表型共變數需要在外部先移除這些欄位，只留下純粹的 SNP 資料
#beta、pair_names是mrs_indeplasso輸出 需要先跑過mrs_indeplasso
#目前的 compute_one() 假設每個 SNP 有三種基因型：0、1、2 僅支援 3x3 的 HLO matrix
risk_mat_test <- mrs.default(
  snp_data = test_data %>% select(-Class),
  mdr_model = model,
  beta = beta_obj$beta,
  pair_names = beta_obj$pair_names
)





#將 target_data 切成2半（validation + test）
#一半從從risk score matrix中選出最好的lambda，
#另一半在test set評估預測力。
#每個人隨機被分配到 group 1（validation）或 group 2（test）
#是獨立抽樣，有可能不平均，也可能有人都被分配到某一組


#risk_mat由mrs.default()得到n個樣本風險分數矩陣 row=個體，col= lambda
#test_data 同樣有 n 個樣本，必須含 Class 欄（0/1）
#lambda要跟ncol(risk_mat)一致，來自 mrs_indeplasso() 的 lambda

res <- mrs_splitvalidate(
  risk_mat = risk_mat_test,
  test_data = test_data,
  lambda = beta_obj$lambda,
  seed = 42
)

print(res$best_lambda)
#直接算每個 lambda 的分數與phenotype的pearson相關
print(res$test_correlation)
print(res$results_table)
print(res$validation_correlation)


#畫圖輸入為
#1.手動組合 mrs_indeplasso() ➜ mrs.default() ➜ mrs_splitvalidate() 的結果
#直接用mrs.pipeline()調用上面三個函數

#lambda tuning 過程
plot_mrs_validate(res, metric = "correlation")
#ROC曲線
plot_mrs_roc(res)
#風險分數的分布比較
plot_mrs_score_dist(res)



#使用 Wilcoxon rank-sum test檢查測試資料的mrs是否在兩個 phenotype 群體之間有顯著差異
test <- mrs_group_test(res)
paste0("Group p-value = ", signif(test$p_value, 3))


#mrs.pipeline整合以下
#風險加權 shrinkage（mrs_indeplasso()）+分數計算（mrs.default()）+lambda選擇與模型評估（mrs_splitvalidate()）
#輸入:MB-MDR 結果（F 值 + HLO）、測試資料（genotype + phenotype）、參數 λ（可選）、隨機種子

#回傳:風險分數建模結果（最佳 λ、分數、準確率、AUC…）、results_table供後續視覺化
res2 <- mrs.pipeline(
  base_output = output_2D,
  test_data = test_data,
  raw_model = raw_model,
  lambda = exp(seq(log(0.001), log(0.1), length.out = 20)),
  seed = 42
)
#比對結果，應該相同
all.equal(res, res2)



#file.edit()
#file.remove("R/mrs.default.R")
#rm(list = c("mrs_indeplasso"))
#devtools::install()
#devtools::check()
#devtools::document()
#devtools::load_all()
#>getwd() 必須是R package目錄
#>1./R要放.R 函數檔案2.DESCRIPTION 檔案：標明 package 名稱,作者,依賴 3.NAMESPACE


#usethis::create_package(".")
#usethis::use_description()
#usethis::use_namespace()


usethis::use_package("dplyr")
usethis::use_package("magrittr")
usethis::use_package("ggplot2")
usethis::use_package("pROC")
