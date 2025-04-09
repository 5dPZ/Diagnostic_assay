###SepsetER Source Code
###Copyright to Sepset Biosciences Inc., subsidy of Asep Medical Holdings Inc.
###Confidential
###The software is copyrighted and the genes referred to are patented to prevent unauthorized commercial exploitation


#A dynamic libpath/setwd is used, so the code can be run anywhere
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(digits=14)
options(scipen=999)

# ###Train data with Caret for Classification ML
# ###input1: expr_data (raw counts), column1 has gene annotation
# ###input2: clin_data_train, column1 has sample id that matches all or training portion of column names of expr_data
# ###input3: x, column name of clin_data with group info
# ###input4: sig, list of signature genes that matches annotation in expr_data column1
# ###input5: hkg, list of housekeeping genes that matches annotation in expr_data column1
# ###optional6: train_name, default="Classification", name of the training set
# ###optional7: signature_name - default - "", name of the signature
# ###optional7: clin_data_validate - default - 0.  If numeric between 0 to 1, % Of input data reserved as validation.
# ###if clin_data_validate=0, then the entire training set is used for validation ***beware of overfit***
# ###if clin_data_validate is a data frame, then it is processed as a seperate dataset the same way as clin_data_train
# ###output1: df, validation report
# ###output2: df, training report
# ###output3: figure ROC curve of best algorithm
# ###output4: algorithm object
# ###optional output5-7: univariate predictor figures

if (!require("tidyverse")) install.packages("tidyverse")
if (!require("readxl")) install.packages("readxl")
if (!require("scales")) install.packages("scales")
if (!require("ggpubr")) install.packages("ggpubr")
if (!require("recipes")) install.packages("recipes")
if (!require("caret")) install.packages("caret")
if (!require("MLeval")) install.packages("MLeval")
if (!require("doParallel")) install.packages("doParallel")
if (!require("AppliedPredictiveModeling")) install.packages("AppliedPredictiveModeling")
if (!require("ComplexHeatmap")) install.packages("ComplexHeatmap")


message("Fetching required expression and clinical data")

###User_defined functions
###User_defined functions
###User_defined functions

REMOVE_EMPTYCOL <- function (df) {
  # Convert POSIXct column to character
  df <- df %>%
    mutate(across(where(is.POSIXct), as.character))
  
  output<-df%>%
    dplyr::select_if(function(x) !(all(is.na(x)) | all(x=="")))
  return(output)
}

# ##remove a row if all are NA or ""
# ##input 1: df, a data frame
# ##output: df, a data frame with empty row removed
# clin_data<-read_csv("data/clin_data_N589_n267.csv")
# clin_data<-REMOVE_EMPTYROW(clin_data)

REMOVE_EMPTYROW <- function (df) {
  # Convert POSIXct column to character
  df <- df %>%
    mutate(across(where(is.POSIXct), as.character))
  
  output<- df[!apply(is.na(df) | df == "", 1, all),]
  return(output)
}

TRANSPOSE_DF<-function (df, col_names="input_colnames") {
  message("All columns (except first column) must be same type")
  df[[1]]<-make.unique(df[[1]])
  if (any(is.na(df[[1]]))) {
    message("Removed NAs from first column")
    df<-filter(df, !is.na(df[[1]]))
  }
  output<-df%>%
    column_to_rownames(var = colnames(df)[1])%>%
    t()%>%
    as.data.frame()%>%
    rownames_to_column(var = col_names)%>%
    as_tibble(.name_repair = 'unique')
  return(output)
}

EXPRESSION_FETCH<-function (expr_data, clin_data, genes, groups) {
  
  if (missing(groups)) {
    clin_data_use<-clin_data
  } else {
    clin_data_use<-dplyr::select(clin_data, 1, any_of(groups))
  }
  
  ###filter by clin data
  expr_data_use<-expr_data %>%
    filter(expr_data[[1]] %in% genes)%>%
    dplyr::select(1, any_of(clin_data_use[[1]]))
  
  for (i in setdiff(genes, expr_data_use[[1]])) {
    warning(paste0("Gene expression data is missing for ", i))}
  
  output<-TRANSPOSE_DF(expr_data_use, col_names = colnames(clin_data_use)[1])%>%
    ###merge with clin_data
    left_join(clin_data_use, by=colnames(clin_data_use)[1])
  return(output)
}

GEOMETRIC_MEAN <- function(x) {
  if (is.numeric(x)) {
    return(exp(mean(log(x))))
  } else if (is.data.frame(x)) {
    output<-x %>% 
      mutate(geometric_mean = apply(dplyr::select_if(., is.numeric), 1, function(row) {
        exp(mean(log(row[row != 0])))
      }))
    return(output)
  } else {
    stop("Input must be a numeric vector or a data frame.")
  }
}

###Start of SepsetER
###Start of SepsetER
###Start of SepsetER

###Fetch expression
expr_data_use<-EXPRESSION_FETCH(expr_data = expr_data,
                                clin_data = clin_data_train,
                                genes = c(sig, hkg))

sig_use<-sort(unique(sig[sig%in% colnames(expr_data_use)]))

filename<-paste0(train_name, 
                 "_N", nrow(clin_data_train), 
                 "_hkg", paste0(str_sub(hkg, 1, 1), collapse = ""),
                 "_sig_", signature_name, 
                 "_n", length(sig_use))

###add 1 to expression data to prevent 0 in calculating geometric mean
expr_data_use[,c(sig_use, hkg)]<-expr_data_use[,c(sig_use, hkg)]+1

###calculate geometric mean for hkg
hkg_use<-dplyr::select(expr_data_use, any_of(hkg))%>%
  GEOMETRIC_MEAN()

hkg_summary<-hkg_use%>%
  colMeans()%>%
  tibble()%>%
  t()
colnames(hkg_summary)<-colnames(hkg_use)
rownames(hkg_summary)<-NULL

expr_data_use<-expr_data_use%>%
  bind_cols(hkg_use[,"geometric_mean"])

###calculate ratios
ratios<-character()

for (i in sig_use) {
  current_ratio<-paste0(i, "_hkg", paste0(str_sub(hkg, 1, 1), collapse = ""))
  expr_data_use[[current_ratio]]<-expr_data_use[[i]]/expr_data_use[["geometric_mean"]]
  ratios<-c(ratios, current_ratio)
}

expr_data_use<-expr_data_use%>%
  dplyr::select(1, any_of(ratios), !!as.symbol(group))%>%
  column_to_rownames(var = colnames(expr_data_use)[1])

ratios<-colnames(expr_data_use)[-ncol(expr_data_use)]

if ("Positive" %in% unique(expr_data_use[[group]])) {
  expr_data_use[[group]]<-factor(expr_data_use[[group]], levels = c("Positive", "Negative"))
} else {
  expr_data_use[[group]]<-as.factor(expr_data_use[[group]])
}

summary(expr_data_use)

###end of data preparation
###end of data preparation
###end of data preparation

###data pre-processing
nzv <- nearZeroVar(expr_data_use, saveMetrics = T)

message(paste0("Removing ", nrow(filter(nzv, nzv==T)), " predictors with non-zero variance"))
###filter out nzv predictor
expr_data_use<-expr_data_use%>%
  dplyr::select(-rownames(filter(nzv, nzv==T)))

###Only perform correlation check when there are more than 1 predictors
if (sum(str_detect(colnames(expr_data_use), "_hkg"))>1) {
  ###filter out correlated predictors
  descrCor<-cor(expr_data_use[, 1:(ncol(expr_data_use)-1)])
  summary(descrCor[upper.tri(descrCor)])
  highlyCor<- findCorrelation(descrCor,
                              cutoff = 0.99)
  
  message(paste0("Removing ", length(highlyCor), " predictors with high correlation"))
  expr_data_use<-dplyr::select(expr_data_use, -any_of(highlyCor))
}

###filter out linear dependencies
comboInfo<- findLinearCombos(expr_data_use[, 1:(ncol(expr_data_use)-1)])

message(paste0("Removing ", length(comboInfo$remove), " predictors with linear dependencies"))
expr_data_use<-dplyr::select(expr_data_use, -any_of(comboInfo$remove))

summary(expr_data_use)

if (is.numeric(clin_data_validate)) {
  message("Spliting data")
  ###split data into training and validation
  set.seed(666)
  trainIndex  <- createDataPartition(expr_data_use[[group]], #stratified with outcome
                                     p= 1-clin_data_validate, 
                                     list=F,
                                     times = 1)
  training <- expr_data_use[trainIndex,]
  validation<- expr_data_use[-trainIndex,]
  
  if (nrow(validation)==0) {
    validation<- training
  }} else {  ###use input validation id to fetch gene expression for validation
    training <- expr_data_use
    
    ###filter by clin data
    validation<-EXPRESSION_FETCH(expr_data = expr_data,
                                 clin_data = clin_data_validate,
                                 genes = c(sig, hkg))
    
    sig_use<-sort(unique(sig[sig%in% colnames(validation)]))
    
    ###add 1 to expression data to prevent 0 in calculating geometric mean
    validation[,c(sig_use, hkg)]<-validation[,c(sig_use, hkg)]+1
    
    ###calculate geometric mean for hkg
    hkg_use<-dplyr::select(validation, any_of(hkg))%>%
      GEOMETRIC_MEAN()
    
    validation<-validation%>%
      bind_cols(hkg_use[,"geometric_mean"])
    
    ###calculate ratios
    ratios<-character()
    
    for (i in sig_use) {
      current_ratio<-paste0(i, "_hkg", paste0(str_sub(hkg, 1, 1), collapse = ""))
      validation[[current_ratio]]<-validation[[i]]/validation[["geometric_mean"]]
      ratios<-c(ratios, current_ratio)
    }
    
    validation<-validation%>%
      dplyr::select(any_of(ratios), !!as.symbol(group))
    
    ratios<-colnames(validation)[-ncol(validation)]
    
    if ("Positive" %in% unique(validation[[group]])) {
      validation[[group]]<-factor(validation[[group]], levels = c("Positive", "Negative"))
    } else {
      validation[[group]]<-as.factor(validation[[group]])
    }
    
    summary(validation)
  }

message(paste0("Training dataset contains ", nrow(training), " samples"))
message(paste0("Validation dataset contains ", nrow(validation), " samples"))

###check outcomes in training and validation
table(training[[group]])
table(validation[[group]])

###transform the predictors
set.seed(666)
preProcValues <- preProcess(as.data.frame(training), method = "knnImpute")

set.seed(666)
training<-predict(preProcValues, training)
set.seed(666)
validation<-predict(preProcValues, validation)

try(rm(list=ls()[str_detect(ls(), "(F|f)it")]))  

###Start parallel computing with all available cores
cl <- makePSOCKcluster(detectCores())
clusterEvalQ(cl, .libPaths())
registerDoParallel(cl)

set.seed(666)
fitControl <- trainControl(
  ## 10-fold CVs
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10, 
  savePredictions = "all",
  summaryFunction = twoClassSummary,
  classProbs = T,
  sampling = "smote", ##SMOTE sampling
  allowParallel = T)

###Use AUROC as metric for model selection
metric="ROC"

set.seed(666)
try(xgbFit <- train(as.formula(paste0(colnames(training)[ncol(training)], " ~ .")), 
                    data=training, 
                    method="xgbTree", 
                    trControl=fitControl,
                    metric=metric,
                    maximize=TRUE))


###Stop parallel computing
stopCluster(cl)

###generate a summary for algorithms
algorithms<-tibble(name=c(
  "XGB"
))

algorithms<-mutate(algorithms, fit=paste0(tolower(name), "Fit"),
                   working=NA,
                   method=NA,
                   label=NA,
                   library=NA,
                   type=NA)

for (n in seq_len(nrow(algorithms))) {
  algorithms$working[n]<-exists(algorithms$fit[n])
  if (algorithms$working[n]) {
    current_fit<-get(algorithms$fit[n])
    algorithms$method[n]<-current_fit$method
    algorithms$label[n]<-current_fit$modelInfo$label
    algorithms$library[n]<-paste0(sort(current_fit$modelInfo$library), collapse = ", ")
    algorithms$type[n]<-paste0(sort(current_fit$modelInfo$type), collapse=", ")
  }}

algorithms<-filter(algorithms, working)
###compare results of 10 algorithms

algorithm.list<-list()

n<- 1
for (i in algorithms$fit) {
  algorithm.list[[n]] <- get(i)
  n <- n + 1
}

names(algorithm.list)<-algorithms$name

set.seed(666)
results <- resamples(algorithm.list)

message("Generating summary figures for algorithms")
# summarize the distributions
summary(results)

pdf(file = paste0(filename, "_bwplot.pdf"), width = 10, height = 10)
# boxplots of results
try(print(bwplot(results)))
dev.off()

pdf(file = paste0(filename, "_dotplot.pdf"), width = 10, height = 10)
# dot plots of results
try(print(dotplot(results)))
dev.off()

###compare results of 10 algorithms

algorithm.list<-list()

n<- 1
for (i in algorithms$fit) {
  algorithm.list[[n]] <- get(i)
  n <- n + 1
}

names(algorithm.list)<-algorithms$name

set.seed(666)
###generate training report
g<-evalm(algorithm.list, 
         gnames = names(algorithm.list), 
         rlinethick = 3/length(algorithm.list),
         cols=hue_pal()(12)[seq_along(algorithm.list)],
         positive = levels(training[[group]])[1])

g$roc+
  theme(legend.position = c(0.88, 0.28))

ggsave(paste0(filename, "_AUCplot.pdf"), device = "pdf", width = 10, height = 10)

training_report<-tibble()
for (n in seq_len(nrow(algorithms))) {
  current.results<-g$stdres[[n]]%>%
    separate(CI, into = c("_lower", "_upper"), sep = "-")
  
  current.results[]<-sapply(current.results, as.numeric)
  
  current.results<-current.results%>%
    rownames_to_column(var = "stat")%>%
    pivot_longer(cols=2:4, values_to = "value")
  
  current.results[[2]][current.results[[2]]=="Score"]<-NA
  
  current.results<-current.results%>%
    unite(col = "stat", sep = "", remove = T, 1:2, na.rm = T)%>%
    TRANSPOSE_DF()%>%
    dplyr::select(-1)%>%
    dplyr::mutate(name=names(g$probs)[n],
                  .before=1)
  
  training_report<-bind_rows(training_report, current.results)
}

training_report<-algorithms%>%
  left_join(training_report, by="name")%>%
  mutate(accuracy_balanced=(SENS+SPEC)/2,
         accuracy_overall=(TP+TN)/(TP+TN+FP+FN))%>%
  arrange(desc(`AUC-ROC`),
          desc(`AUC-PRG`),
          desc(`AUC-PR`))

training_report<-REMOVE_EMPTYCOL(training_report)

write_csv(training_report, 
          paste0(filename, "_training_report_",
                 str_replace_all(as.character(Sys.Date()), "\\-", ""),
                 ".csv"))

###20211206 Generate a ROC figure
###dplyr::select prediction result for best algorithm
###generate training report
set.seed(666)
g<-evalm(algorithm.list[training_report$name[1]], 
         gnames = training_report$label[1],
         positive = levels(training[, ncol(training)])[1])

g$roc+
  theme(legend.position = c(0.8, 0.08))

ggsave(paste0(filename, "_AUCplot_top.pdf"), device = "pdf", width = 7, height = 7)


###Generating a validation report
###Generating a validation report
###Generating a validation report

###remove RF algorithm if clin_data_validate=0, i.e use training for validation, due to overfit 
if (is.numeric(clin_data_validate)) {
  if (clin_data_validate==0) {
    algorithm.list["RF"]<-NULL
    message("Removed RF algorithm to prevent overfit")
  }
  clin_data_pred<-filter(clin_data_train, clin_data_train[[1]]%in%rownames(validation))
} else {    
  clin_data_pred<-clin_data_validate
}

validation_report<-tibble()
for (n in seq_len(nrow(algorithms))) {
  current.fit<-get(algorithms$fit[n])
  set.seed(666)
  current.matrix <- predict(current.fit, newdata = validation, type = "prob")
  
  ###Remove any rownames, drop any empty rows, then fill any NaN as 1
  rownames(current.matrix) <- NULL
  current.matrix <- current.matrix %>%
    REMOVE_EMPTYROW( ) %>%
    mutate_all(~ifelse(is.nan(.), 1, .))
  
  current.table <- evalm(data.frame(current.matrix, 
                                    validation[[group]][as.numeric(rownames(current.matrix))]),
                         positive = levels(validation[[group]])[1])
  
  current.results<-current.table$stdres[[1]]%>%
    separate(CI, into = c("_lower", "_upper"), sep = "-")
  
  current.results[]<-sapply(current.results, as.numeric)
  
  current.results<-current.results%>%
    rownames_to_column(var = "stat")%>%
    pivot_longer(cols=2:4, values_to = "value")
  
  current.results[,2][current.results[,2]=="Score"]<-NA
  
  current.results<-current.results%>%
    unite(col = "stat", sep = "", remove = T, 1:2, na.rm = T)%>%
    TRANSPOSE_DF()%>%
    dplyr::select(-1)%>%
    dplyr::mutate(name=algorithms$name[n],
                  .before=1)
  
  validation_report<-bind_rows(validation_report, current.results)
  
  set.seed(666)
  clin_data_pred<-clin_data_pred%>%
    bind_cols(predict(current.fit, newdata = validation, type= "prob"))%>%
    mutate(prediction=predict(current.fit, newdata = validation))
  colnames(clin_data_pred)[(length(colnames(clin_data_pred))-2):length(colnames(clin_data_pred))]<-c(paste0(names(algorithm.list)[n], "_", levels(validation[[group]])[1]),
                                                                                                     paste0(names(algorithm.list)[n], "_", levels(validation[[group]])[2]),
                                                                                                     names(algorithm.list)[n])
}

validation_report<-algorithms%>%
  left_join(validation_report, by="name")%>%
  mutate(accuracy_balanced=(SENS+SPEC)/2,
         accuracy_overall=(TP+TN)/(TP+TN+FP+FN))%>%
  arrange(desc(`AUC-ROC`),
          desc(`AUC-PRG`),
          desc(`AUC-PR`))

validation_report<-REMOVE_EMPTYCOL(validation_report)

write_csv(validation_report, 
          paste0(filename, "_validation_report_",
                 str_replace_all(as.character(Sys.Date()), "\\-", ""),
                 ".csv"))

write_csv(as_tibble(clin_data_pred, .name_repair = 'unique'), 
          paste0(filename, "_clin_data_pred_",
                 str_replace_all(as.character(Sys.Date()), "\\-", ""),
                 ".csv"))

clin_data_pred<-bind_cols(clin_data_pred, validation[,-ncol(validation)])%>%
  arrange(!!as.symbol(group))
rownames(clin_data_pred)<-NULL

###heatmap for prediction vs. reference 
df<-data.frame(dplyr::select(clin_data_pred, any_of(group), any_of(names(algorithm.list))))
ha <- HeatmapAnnotation(
  df=df
)

validation<-TRANSPOSE_DF(dplyr::select(clin_data_pred, 1, any_of(ratios)))%>%
  column_to_rownames(colnames(.)[1])
# Combine the heatmap and the annotation
pdf(paste0(filename, "_heatmap.pdf"),
    width=20, height=20)
print(Heatmap(validation,
              top_annotation = ha,
              column_split = clin_data_pred[[group]],
              show_column_dend = F,
              show_row_dend = F,
              show_column_names = F,
              show_row_names = T,
))
dev.off()

training_report<-training_report%>%
  dplyr::select(name, `AUC-ROC`, `AUC-PRG`, `AUC-PR`)%>%
  dplyr::rename(model=name, ROC_training=`AUC-ROC`, PRG_training=`AUC-PRG`, PR_training=`AUC-PR`)

validation_report<-validation_report%>%
  dplyr::select(name, `AUC-ROC`, `AUC-PRG`, `AUC-PR`)%>%
  dplyr::rename(model=name, ROC_validation=`AUC-ROC`, PRG_validation=`AUC-PRG`, PR_validation=`AUC-PR`)

model_summary<-training_report%>%
  left_join(validation_report, by="model")%>%
  mutate(ROC_mean=round((ROC_training*ROC_validation)^0.5, 4),
         PRG_mean=round((PRG_training*PRG_validation)^0.5, 4),
         PR_mean=round((PR_training*PR_validation)^0.5, 4))%>%
  arrange(desc(ROC_mean),
          desc(PRG_mean),
          desc(PR_mean))

write_csv(model_summary,paste0(filename, "_model_summary_",
                               str_replace_all(as.character(Sys.Date()), "\\-", ""),
                               ".csv"))

model_summary<-model_summary%>%
  pivot_longer(cols = 2:10, names_to = "metric", values_to = "performance")%>%
  separate(metric, into = c("metric", "dataset"), sep = "_")%>%
  filter(metric == "ROC")

model_summary$dataset <- factor(model_summary$dataset, levels = c("training", "validation", "mean"))

ggplot(model_summary, aes(x=dataset, y=performance, group=model, color=model))+
  geom_line()

ggsave(paste0(filename, "_model_summary.pdf"), device = "pdf", width = 6, height = 10)

save(list=c("training", "hkg_summary", "preProcValues", algorithms$fit), 
     file=paste0(filename, ".Rda"))

