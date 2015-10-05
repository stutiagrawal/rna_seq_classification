library(e1071)
library(preprocessCore)
library(qvalue)
library(VGAM)

data = read.table("/Users/stuti/Data/gec22/out_900.txt", header=T, sep="\t")
data = data.frame(data)

#Perform the normalization
upper_quantile = apply(data, 2, function(x) quantile(x, 0.75))
upper_quantile <- data.frame(upper_quantile)
for (i in 1:nrow(upper_quantile)){
    data[,i] <- (data[,i]/upper_quantile[i,1]) * 1000
}

#Perform log transformation
data_log = log2(data + 0.0000000001)
data_trans = t(data_log)
temp = lapply(data.frame(data_trans), var)
w = which(temp==0)
data_trans = data_log[-w, ]

#remove low-expression genes
threshold = -20
t <- apply(data_trans, 1, mean)
selected_genes <-  which(t > threshold)
data_trans <- data_trans[selected_genes,]

#Perform student's t-test
n = 482
p_val <- apply(data_trans, 1, function(x) t.test(x[1:n], x[n+1:length(x)])$p.value)
p_val <- data.frame(p_val)
p_val$gene <- rownames(p_val)
p <- p_val[order(p_val$p_val),]
2015-1961
#Select discordant samples
discordant <- read.table("/Users/stuti/Data/gec22/star_cufflinks_pipeline/discordant_lusc", header=F)
discordant <- gsub("-",".", discordant$V1)
v <- numeric()
for (i in 1:length(discordant)){
    v <- c(v, which(grepl(discordant[i], colnames(data_trans))))
}

lung_discordant <- data.frame(t(data_trans[,v]))
lung_other <- data.frame(t(data_trans[,-v]))

#Apply labels
label = rep("LUAD", length(rownames(lung_other)))
lusc <- which(grepl("LUSC",rownames(lung_other)))
label[lusc] <- "LUSC"
lung_other$condition <- as.factor(label)
test_set_labels = rep("LUSC", length(rownames(lung_discordant)))
lung_discordant$condition <- as.factor(test_set_labels)

#Create intial train and test set
get_train_set <- function(dataset, parts=3){
    index <- 1:nrow(dataset)
    testindex <- sample(index, trunc(length(index)/parts))
    testset <- dataset[testindex,]
    trainset <- dataset[-testindex,]
    return(list(trainset, testset))
}

#Perform SVM based classification
classify_svm <- function(train, test){
    label="condition"
    
    test_label_index <- which(colnames(test) == label)
    svm.model <- best.tune(svm, condition~., data=train, 
                     ranges = list(gamma = c(0, 0.2, 0.4, 0.6, 0.8, 1), cost = 2^(0:4)),
                     tunecontrol = tune.control(sampling = "fix"), kernel="linear",
                     type="C-classification", cross=5)
    
    svm.pred <- predict(svm.model, test[,-test_label_index])
    t <- table(pred = svm.pred, true = test[,test_label_index])
    classAgreement(t)$diag   
}

result <- matrix(nrow=4, ncol=2)
i = 0
#get accuracy for different numbers of genes
for(n in c(50, 100, 500, 1000)){
    i <- i + 1
    w <- p$gene[1:n]
    w <- gsub("-", ".", w)
    w <- c(w, "condition")
    lung_other_sel <- lung_other[,w]
    t <- get_train_set(lung_other_sel)
    train <- data.frame(t[1])
    test <- data.frame(t[2])
    x <- classify_svm(train, test)
    result[i,1] <- n
    result[i,2] <- x
}

