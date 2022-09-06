
## the interface of Seurat object and SingleCellNet
Interface_Seurat_SCN <- function(data, ClassList) {
    data <- subset(data, (Class == ClassList[1] | Class == ClassList[2]))

    exp <- data@assays$RNA@counts
    anno <- data@meta.data

    ### training ScingleCellNet classifier
    delete_result <- remove_celltype_less_n(exp, anno, 5)
    exp <- delete_result[[1]]
    anno <- delete_result[[2]]
    rm(delete_result)

    anno<- cbind("cellid" = rownames(anno), anno)
    anno <- droplevels(anno)
    
    return_list=list()
    return_list[["exp_matrix"]]=exp
    return_list[["anno_matrix"]]=anno
    return(return_list)
}

## trainging a singlecellnet classifier
training_SCN_classifier <- function(exp_matrix,anno_matrix, dataset_id,
                                    nTopGenes, nTopGenePairs, ncells,
                                    nrand = 0) {

    ## training set
    set.seed(619)
    stList <- splitCommon(sampTab = anno_matrix, ncells = ncells, dLevel = "Cell_subtype") ## ncell can adjust
    stTrain <- stList[[1]]
    expTrain <- exp_matrix[, rownames(stTrain)]
    expTrain[1:5, 1:5]

    ## test-set
    stTestList <- splitCommon(sampTab = stList[[2]], ncells = ncells, dLevel = "Cell_subtype") # normalize validation data so that the assessment is as fair as possible
    stTest <- stTestList[[1]]
    expTest <- exp_matrix[, rownames(stTest)]

    ## cell subtype classifier model
    class_info <- scn_train(
        stTrain = stTrain, expTrain = expTrain,
        nTopGenes = nTopGenes, nRand = 70,
        nTrees = 1000, nTopGenePairs = nTopGenePairs,
        dLevel = "Cell_subtype", colName_samp = "cellid"
    )

    ## save the classifier
    saveRDS(class_info, file = paste0(
        "./result/", dataset_id,
        "_class_info(ntopgenes=", nTopGenes,
        ",ntoppaies=", nTopGenePairs,
        ",ncells=", ncells, ").rds"
    ))

    ## assess classifier
    classRes_val_all <- scn_predict(cnProc = class_info[["cnProc"]], expDat = expTest, nrand = nrand)

    tm_heldoutassessment <- assess_comm(
        ct_scores = classRes_val_all,
        stTrain = stTrain, stQuery = stTest,
        dLevelSID = "cellid",
        classTrain = "Cell_subtype", classQuery = "Cell_subtype", nRand = 0
    )

    pdf(file = paste0(
        "./figure/", dataset_id,
        " Assess cell-subtype classifier(ntopgenes=", nTopGenes,
        ",ntoppairs=", nTopGenePairs,
        ",ncells=", ncells, ").pdf"
    ), width = 15, height = 7.5)
    p <- plot_PRs(tm_heldoutassessment)
    print(p)
    dev.off()

    numPairs <- length(class_info$cnProc$xpairs)
    gpTab <- compareGenePairs(
        query_exp = expTest, training_exp = expTrain, training_st = stTrain,
        classCol = "Cell_subtype", sampleCol = "cellid",
        RF_classifier = class_info$cnProc$classifier, numPairs = numPairs, trainingOnly = TRUE
    )
    saveRDS(gpTab, file = paste0("./result/", dataset_id, "(ncells=", ncells, ")", "top_genepairs.rds"))

    return(NULL)
}

## debug some functions in SCN
if (T) {
    assess_comm <- function(ct_scores, # matrix of classification scores, rows = classifiers, columns = samples, colnames=sampleids || where query cells is in the training
                            stTrain, # sample table where cells in query are in the training
                            stQuery,
                            resolution = 0.005, # increment at which to evalutate classification
                            nRand = 0,
                            dLevelSID = "sample_name",
                            classTrain = "cell_ontology_class",
                            classQuery = "description2", # query data
                            AUCmethod = "trapezoid") {
        shared_cell_type <- intersect(unique(stTrain[, classTrain]), unique(stQuery[, classQuery]))
        stVal_com <- stQuery[which(stQuery[, classQuery] %in% shared_cell_type), ]

        if (nRand > 0) {
            tmp <- as.data.frame(matrix("rand", nrow = nRand, ncol = (ncol(stVal_com))))
            colnames(tmp) <- colnames(stVal_com)
            tmp[, dLevelSID] <- colnames(ct_scores)[(ncol(ct_scores) - nRand + 1):ncol(ct_scores)]
            rownames(tmp) <- tmp[, dLevelSID]
            stVal_com <- rbind(stVal_com, tmp)
        }

        cells_sub <- as.character(stVal_com[, dLevelSID])

        # subsetting the ct_scores where the cells' true identity is within the range of the classifiers
        ct_score_com <- ct_scores[, cells_sub]
        report <- list()
        ct_scores_t <- t(ct_score_com)
        true_label <- as.character(stVal_com[, classQuery])
        # multiLogLoss
        names(true_label) <- rownames(ct_scores_t)
        if (is.matrix(true_label) == FALSE) {
            y_true <- model.matrix(~ 0 + ., data.frame(as.character(true_label)))
        }
        eps <- 1e-15
        y_pred <- pmax(pmin(ct_scores_t, 1 - eps), eps)
        multiLogLoss <- (-1 / nrow(ct_scores_t)) * sum(t(y_true) %*% log(y_pred)) # want columns to be the cell types for y_pred
        # cohen's kappa, accuracy
        pred_label <- c()
        pred_label <- colnames(ct_scores_t)[max.col(ct_scores_t, ties.method = "random")]

        cm <- as.matrix(table(Actual = true_label, Predicted = pred_label))

        # in case of misclassfication where there are classifiers that are not used
        if (length(setdiff(unique(true_label), unique(pred_label))) != 0) {
            misCol <- setdiff(unique(true_label), unique(pred_label))
            for (i in 1:length(misCol)) {
                cm <- cbind(cm, rep(0, nrow(cm)))
            }
            colnames(cm)[(ncol(cm) - length(misCol) + 1):ncol(cm)] <- misCol
        }

        if (length(setdiff(unique(pred_label), unique(true_label))) != 0) {
            misRow <- setdiff(unique(pred_label), unique(true_label))
            for (i in 1:length(misRow)) {
                cm <- rbind(cm, rep(0, ncol(cm)))
            }
            rownames(cm)[(nrow(cm) - length(misRow) + 1):nrow(cm)] <- misRow
        }

        cm <- cm[, colnames(cm)[match(rownames(cm), colnames(cm))]]


        # sort table names accordigly

        n <- sum(cm) # number of instances
        nc <- nrow(cm) # number of classes
        diag <- diag(cm) # number of correctly classified instances per class
        rowsums <- apply(cm, 1, sum) # number of instances per class
        colsums <- apply(cm, 2, sum) # number of predictions per class
        p <- rowsums / n # distribution of instances over the actual classes
        q <- colsums / n # distribution of instances over the predicted classes
        expAccuracy <- sum(p * q)
        accuracy <- sum(diag) / n

        # PR
        confusionMatrix <- cn_classAssess(ct_score_com, stVal_com, classLevels = classQuery, dLevelSID = dLevelSID, resolution = resolution)
        PR_ROC <- cal_class_PRs(confusionMatrix)
        nonNA_PR <- PR_ROC[which(!is.nan(PR_ROC$recall)), ]
        nonNA_PR[which((nonNA_PR$TP == 0 & nonNA_PR$FP == 0)), "precision"] <- 1

        w <- c()
        areas <- c()
        for (i in 1:length(unique(nonNA_PR$ctype))) {
            tmp <- nonNA_PR[which(nonNA_PR$ctype %in% unique(nonNA_PR$ctype)[i]), ]
            area <- DescTools::AUC(tmp$recall, tmp$precision, method = AUCmethod)
            areas <- c(areas, area[1])
            w <- c(w, sum(stVal_com[, classQuery] %in% unique(nonNA_PR$ctype)[i]) / nrow(stVal_com))
        }

        report[["accuracy"]] <- accuracy
        report[["kappa"]] <- (accuracy - expAccuracy) / (1 - expAccuracy)
        report[["AUPRC"]] <- areas
        report[["AUPRC_w"]] <- mean(areas)
        report[["AUPRC_wc"]] <- weighted.mean(areas, w)
        report[["multiLogLoss"]] <- multiLogLoss
        report[["cm"]] <- cm
        report[["confusionMatrix"]] <- confusionMatrix
        report[["nonNA_PR"]] <- nonNA_PR
        report[["PR_ROC"]] <- PR_ROC

        return(report)
    }

    cal_class_PRs <- function(assessed) {
        ctts <- names(assessed)
        df <- data.frame()
        for (ctt in ctts) {
            tmp <- assessed[[ctt]]
            tmp <- cbind(tmp, ctype = ctt)
            df <- rbind(df, tmp)
        }

        prsAll <- transform(df,
            TP = as.numeric(as.character(TP)),
            TN = as.numeric(as.character(TN)),
            FN = as.numeric(as.character(FN)),
            FP = as.numeric(as.character(FP))
        )

        precfunc <- function(df) {
            ans <- vector()
            for (i in 1:nrow(df)) {
                ans <- append(ans, df[i, "TP"] / (df[i, "TP"] + df[i, "FP"]))
            }
            ans
        }

        sensfunc <- function(df) {
            ans <- vector()
            for (i in 1:nrow(df)) {
                ans <- append(ans, df[i, "TP"] / (df[i, "TP"] + df[i, "FN"]))
            }
            ans
        }

        precs <- precfunc(prsAll)
        sens <- sensfunc(prsAll)
        prsAll2 <- cbind(prsAll, data.frame(recall = sens, precision = precs))
    }

    ## plot PRs and calculate AUPR
    plot_PRs <- function(assessed, collapse = F) {
        if (collapse) {
            ggplot(data = assessed$nonNA_PR, aes(x = as.numeric(as.vector(recall)), y = as.numeric(as.vector(precision)), colour = ctype)) +
                geom_point(size = 0.5, alpha = 0.5) +
                geom_path(size = 0.5, alpha = 0.75) +
                theme_bw() +
                xlab("Recall") +
                ylab("Precision") +
                theme(axis.text = element_text(size = 5)) +
                ggtitle("Classification performance_PR Curve")
        } else {
            library(DescTools)
            areas <- c()
            for (i in 1:length(unique(assessed$nonNA_PR$ctype))) {
                tmp <- assessed$nonNA_PR[which(assessed$nonNA_PR$ctype %in% unique(assessed$nonNA_PR$ctype)[i]), ]
                area <- DescTools::AUC(tmp$recall, tmp$precision, method = "trapezoid")
                areas <- c(areas, area)
            }
            plist <- list()
            for (i in 1:table(table(assessed$nonNA_PR$ctype) == 0)[1]) {
                a <- assessed$nonNA_PR[(1 + 201 * (i - 1)):(201 + 201 * (i - 1)), ]
                celltype <- as.character(a$ctype[1])[1]
                tem <- ggplot(data = a, aes(x = as.numeric(as.vector(recall)), y = as.numeric(as.vector(precision)))) +
                    geom_point(size = 0.5, alpha = 0.5) +
                    geom_path(size = 0.5, alpha = 0.75) +
                    theme_bw() +
                    labs(
                        x = "Recall", y = "Precision",
                        title = paste0(celltype),
                        subtitle = paste0("AUPR=", round(areas[i], digits = 4)), vjust = 1.5, size = 4.0
                    ) +
                    theme(axis.text = element_text(size = 5))
                plist[[i]] <- tem
            }
            p <- plist[[1]]
            for (i in 2:length(plist)) {
                p <- p + plist[[i]]
            }
            return(p)
        }
    }

    hm_gpa_sel <- function(expDat,
                           genes,
                           grps, ## vector of cellnames -> grp label
                           maxPerGrp = 100,
                           cRow = FALSE,
                           cCol = FALSE,
                           limits = c(0, 10),
                           toScale = FALSE,
                           fontsize_row = 4,
                           reOrderCells = FALSE) {
        allgenes <- rownames(expDat)
        missingGenes <- setdiff(genes, allgenes)
        if (length(missingGenes) > 0) {
            cat("Missing genes: ", paste0(missingGenes, collapse = ","), "\n")
            genes <- intersect(genes, allgenes)
        }

        value <- expDat[genes, ]
        if (toScale) {
            if (class(value)[1] != "matrix") {
                value <- t(scale(Matrix::t(value)))
            } else {
                value <- t(scale(t(value)))
            }
        }

        value[value < limits[1]] <- limits[1]
        value[value > limits[2]] <- limits[2]

        groupNames <- unique(grps)
        if (reOrderCells) {
            grps <- grps[order(grps)]
            groupNames <- sort(unique(grps))
        }

        cells <- names(grps)

        cells2 <- vector()
        for (groupName in groupNames) {
            xi <- which(grps == groupName)
            if (length(xi) > maxPerGrp) {
                tmpCells <- sample(cells[xi], maxPerGrp)
            } else {
                tmpCells <- cells[xi]
            }
            cells2 <- append(cells2, tmpCells)
        }
        value <- value[, cells2]

        xcol <- colorRampPalette(rev(brewer.pal(n = 12, name = "Paired")))(length(groupNames))
        names(xcol) <- groupNames
        anno_colors <- list(group = xcol)

        xx <- data.frame(group = as.factor(grps))
        rownames(xx) <- cells

        val_col <- colorRampPalette(rev(brewer.pal(n = 12, name = "Spectral")))(25)
        # val_col <- colorRampPalette(brewer.pal(n = 12,name = "Spectral"))(100)

        p <- pheatmap(value,
            cluster_rows = cRow, cluster_cols = cCol, color = val_col,
            show_colnames = FALSE, annotation_names_row = FALSE, show_rownames = FALSE,
            ##        annotation_col = annTab,
            annotation_col = xx,
            annotation_names_col = FALSE, annotation_colors = anno_colors, fontsize_row = fontsize_row
        )
        return(p)
    }
}

## applied lasso-cox model to bulid signature
LassoCox_signature <- function(stable_pairs_list, celltype_list, iteration_times, all_exp, all_clinical, training_size, test_size, validation_sets) {
    dir.create(paste0("./model/"))
    for (celltype in celltype_list) {
        t1 <- proc.time()
        setwd(paste0("./model/"))
        dir.create(path = paste0("./", celltype, "/"))
        setwd(paste0("./", celltype, "/"))

        ## each ncells train n times
        set.seed(619)
        for (iteration in c(1:iteration_times)) {
            cat("Interation:", iteration, "\n")
            dir.create(path = paste0("./iteration", iteration, "/"))
            setwd(paste0("./iteration", iteration, "/"))
            ## traing sets
            traing_index <- sample(1:nrow(all_clinical), training_size, replace = F)
            training_exp <- all_exp[, traing_index]
            training_clinical <- all_clinical[traing_index, ]

            ## testing set
            testing_index <- c(1:nrow(all_clinical))[-traing_index]
            testing_exp <- all_exp[, testing_index]
            testing_clinical <- all_clinical[testing_index, ]

            rm(traing_index, testing_index)

            celltype_specific_genes_pairs <- get_stable_specific_genepairs(stable_pairs_list[[celltype]])

            training_exp_genepairs_rank <- calculate_genepair_rank(training_exp, celltype_specific_genes_pairs)

            ## traning model
            training_exp_genepairs_rank <- remove_all_1or0_features(training_exp_genepairs_rank)

            ## lasso-cox model with L1 regulation
            train_model <- tryCatch(
                {
                    training_Lasso_Cox_model(training_exp_genepairs_rank, training_clinical)
                },
                error = function(x) {
                    cat("Interation:", iteration, " caused errors.", "\n")
                    return("ERROR_label")
                }
            )
            rm(training_exp_genepairs_rank)
            if (train_model == "ERROR_label" & iteration == iteration_times) {
                break
            }
            if (train_model == "ERROR_label" & iteration != iteration_times) {
                setwd("../")
                next
            }

            if (length(table(as.matrix(train_model$`gene-pair weights`))) < 6) {
                cat("Interation:", iteration, " model features are too less.", "\n")
                setwd("../")
                next
            }

            ## search the best group decision boundary
            cutoff <- search_group_boundary(training_exp, training_clinical, train_model, way = "cutoff")
            if (cutoff == "ERROR_label") {
                cat("Interation:", iteration, " surv_cutpoint caused errors.", "\n")
                setwd("../")
                next
            }

            ## save model
            train_model[["cutoff"]] <- cutoff
            saveRDS(train_model,
                file = paste0(
                    "Traning set(size=",
                    nrow(training_clinical), ") prognostic-signature of ",
                    celltype, "(ncells=", ncells, ")_model.rds"
                )
            )

            ## traning set  K-M plot
            group_label <- grouping(train_model, training_exp, training_clinical, cutoff)
            if (length(table(group_label)) < 2) {
                setwd("../")
                next
            }
            cat("The training set:GSE39582 cutoff: ", cutoff, "\n")
            cat("The training set:GSE39582 labels: ", table(group_label), "\n")
            KMplot_list <- KM_plot(training_clinical, group_label)
            pdf(
                file = paste0("Traning set(size=", nrow(training_clinical), ") prognostic-signature of ", celltype, "(ncells=", ncells, ").pdf"),
                width = 8, height = 6
            )
            p <- KMplot_list$p
            print(p)
            dev.off()

            ## test
            group_label <- grouping(train_model, testing_exp, testing_clinical, cutoff)
            if (length(table(group_label)) < 2) {
                setwd("../")
                next
            }
            KMplot_list <- KM_plot(testing_clinical, group_label)
            cat("The testing set labels: ", table(group_label), "\n")
            if (KMplot_list$logrank_p > 0.05) {
                cat(paste0("Testing set of iteration ", iteration, " is not significant."), "\n")
                setwd("../")
                next
            }

            p <- KMplot_list$p
            pdf(
                file = paste0("Testing set(size=", nrow(testing_clinical), ") prognostic-signature of ", celltype, "(ncells=", ncells, ").pdf"),
                width = 8, height = 6
            )
            print(p)
            dev.off()

            ## validate
            for (i in validation_sets) {
                group_label <- grouping(train_model, exp_matrix_list[[i]], clinic_list[[i]], cutoff)
                if (length(table(group_label)) < 2) {
                    next
                }
                cat("The validate set(", i, ") labels: ", table(group_label), "\n")
                KMplot_list <- KM_plot(clinic_list[[i]], group_label)
                if (KMplot_list$logrank_p > 0.05) {
                    cat(paste0("Validation set ", i, " of iteration ", iteration, " is not significant."), "\n")
                    break
                }
                pdf(
                    file = paste0("Validation set(", i, " size=", nrow(clinic_list[[i]]), ") prognostic-signature of ", celltype, "(ncells=", ncells, ").pdf"),
                    width = 8, height = 6
                )
                p <- KMplot_list$p
                print(p)
                dev.off()
            }
            setwd("../")
        }
        t2 <- proc.time()
        t <- t2 - t1
        cat(paste0("Cell-subtype: ", celltype, "\t", "ncells= ", ncells, "\t", "using time:", t[3][[1]], " seconds.", "\n"))
        setwd("../../")
    }
    
    return(NULL)
}

## get the gene-pairs from stable cellsubtype-specific gene pairs
get_stable_specific_genepairs <- function(stable_pairs_ele) {
    gene_pairs <- stable_pairs_ele

    gene_1 <- c()
    gene_2 <- c()
    for (i in 1:length(gene_pairs)) {
        tem_1 <- strsplit(gene_pairs[i], split = "_")[[1]][1]
        tem_2 <- strsplit(gene_pairs[i], split = "_")[[1]][2]
        gene_1 <- c(gene_1, tem_1)
        gene_2 <- c(gene_2, tem_2)
    }

    allPairs <- data.frame(Gene1 = gene_1, Gene2 = gene_2, stringsAsFactors = F)
    return(allPairs)
}

## calculate the gene-pairs rank in bulk exp data
calculate_genepair_rank <- function(bulkdata, celltype_specific_genes_pairs) {
    ## select the gene pairs invovle in bulk exp
    celltype_specific_genes_pairs <- celltype_specific_genes_pairs[celltype_specific_genes_pairs[, 1] %in% rownames(bulkdata), ]
    celltype_specific_genes_pairs <- celltype_specific_genes_pairs[celltype_specific_genes_pairs[, 2] %in% rownames(bulkdata), ]

    gene_pair_rank <- matrix(data = NA, ncol = ncol(bulkdata), nrow = nrow(celltype_specific_genes_pairs))
    colnames(gene_pair_rank) <- colnames(bulkdata)
    rownames(gene_pair_rank) <- paste0(celltype_specific_genes_pairs[, 1], "__", celltype_specific_genes_pairs[, 2])

    for (i in 1:ncol(bulkdata)) {
        patient_exp <- bulkdata[, i]
        names(patient_exp) <- rownames(bulkdata)

        ## calculate the gene pair rank
        gene1 <- celltype_specific_genes_pairs[, 1]
        gene2 <- celltype_specific_genes_pairs[, 2]

        patient_rank <- patient_exp[gene1] - patient_exp[gene2]
        patient_rank <- ifelse(patient_rank > 0, 1, 0)

        gene_pair_rank[, i] <- patient_rank
    }
    return(gene_pair_rank)
}

## removing some features which all equal 1 or 0
remove_all_1or0_features <- function(exp_genepairs_rank) {
    exp_genepairs_rank <- as.data.frame(exp_genepairs_rank)

    index <- c()

    for (i in 1:nrow(exp_genepairs_rank)) {
        tem <- as.numeric(exp_genepairs_rank[i, ])
        if (length(table(tem)) > 1) {
            index <- c(index, i)
        }
    }
    rm(i)
    exp_genepairs_rank <- exp_genepairs_rank[index, ]
    return(exp_genepairs_rank)
}

## select prognostic-related gene pairs and build lasso-cox model with L1 regulation
training_Lasso_Cox_model <- function(exp, clinical) {

    ## survival time equal to zero should be remove
    zero_time <- clinical[, 1] %in% 0

    ## x matrix(features)
    exp <- exp[, !zero_time]

    ## y matrix(survival time and status)
    y <- as.matrix(Surv(time = as.double(clinical[, 1]), event = as.double(clinical[, 2])))
    y <- y[!zero_time, ]

    ## features
    x <- as.matrix(t(exp))

    ## 生成模型的特征
    ## 挑选正则化程度最高的lamda
    ## C-index 衡量error, alpha=0, L1 regulation, cox family
    cindex_lamda_set <- cv.glmnet(x, y, family = "cox", type.measure = "C", alpha = 1, nfolds = 5, maxit = 100000)

    pdf(file = "lasso-cox model cindex lamda_set.pdf", width = 8, height = 6)
    print(plot(cindex_lamda_set))
    dev.off()

    feature_wights <- coef(cindex_lamda_set, s = "lambda.1se")

    model <- list()
    model[["lasso-cox model"]] <- cindex_lamda_set
    model[["gene-pair weights"]] <- feature_wights
    model[["cindex"]] <- cindex_lamda_set$cvm[cindex_lamda_set$index["1se", ]]

    cat("Lasso-Cox model training success!", "\n")

    coxph_model <- coxph(Surv(time = as.double(clinical[, 1]), event = as.double(clinical[, 2])) ~ x)

    pdf(file = "Check the coxph and glmnet coef.pdf", width = 8, height = 6)
    plot(coef(cindex_lamda_set, s = 0), coef(coxph_model))
    abline(0, 1)
    dev.off()
    return(model)
}

## search the best group decision boundary
search_group_boundary <- function(exp, clincal, train_model, way) {
    if (way == "cutoff") {
        pairs <- train_model$`gene-pair weights`@Dimnames[[1]]
        coef <- train_model$`gene-pair weights`
        score <- riskScore(exp, pairs, coef)

        data <- cbind(time = clincal[, 1], event = clincal[, 2], score = score)
        data <- as.data.frame(data)

        cutoff <- tryCatch(
            {
                surv_cutpoint(data, time = "time", event = "event", variables = "score")
            },
            error = function(x) {
                return("ERROR_label")
            }
        )
        if (cutoff == "ERROR_label") {
            return(cutoff)
        } else {
            p <- plot(cutoff, "score", palette = "npg")

            pdf(file = "Determine the optimal cutpoint of training set.pdf", width = 8, height = 6)
            print(p)
            dev.off()

            return(cutoff$cutpoint$cutpoint)
        }
    }
    if (way == "mean") {
        pairs <- train_model$`gene-pair weights`@Dimnames[[1]]
        coef <- train_model$`gene-pair weights`
        score <- riskScore(exp, pairs, coef)

        return(mean(score))
    }
}

## calculate risk score
riskScore <- function(exp, pairs, coef) {
    gene1 <- c()
    gene2 <- c()

    score_all_patients <- c()

    for (i in 1:length(pairs)) {
        gene1 <- c(gene1, strsplit(pairs[i], split = "__")[[1]][1])
        gene2 <- c(gene2, strsplit(pairs[i], split = "__")[[1]][2])
    }
    for (i in 1:ncol(exp)) {
        tem_exp <- exp[, i]
        names(tem_exp) <- rownames(exp)
        score <- ifelse(tem_exp[gene1] - tem_exp[gene2] > 0, 1, 0)
        score <- score %*% coef
        score <- as.numeric(score)
        score <- exp(score)

        score_all_patients <- c(score_all_patients, score)
    }
    return(score_all_patients)
}

## grouping
grouping <- function(train_model, exp, clincal, cutoff) {
    pairs <- train_model$`gene-pair weights`@Dimnames[[1]]
    coef <- train_model$`gene-pair weights`
    score <- riskScore(exp, pairs, coef)
    label <- ifelse(score >= cutoff, "high risk", "low risk")
    return(label)
}

## K-M plot
KM_plot <- function(clinical, label) {
    clinical[, 1] <- as.numeric(clinical[, 1])
    clinical[, 2] <- as.numeric(clinical[, 2])

    data <- cbind(clinical, score = label)

    fit <- surv_fit(Surv(clinical[, 1], clinical[, 2]) ~ label, data = data)
    logrank_pvalue <- survdiff(Surv(clinical[, 1], clinical[, 2]) ~ label, data = data)
    logrank_pvalue <- 1 - pchisq(logrank_pvalue$chisq, length(logrank_pvalue$n) - 1)

    fit2 <- coxph(Surv(clinical[, 1], clinical[, 2]) ~ label, data = data)
    fit2_sum <- summary(fit2)
    c_index <- fit2_sum$concordance[1]
    c_index <- round(digits = 3, c_index)

    p <- ggsurvplot(fit,
        data = data,
        linetype = c("solid", "solid"),
        surv.median.line = "hv", surv.scale = "percent",
        pval = T, risk.table = T,
        conf.int = T, conf.int.alpha = 0.1, conf.int.style = "ribbon",
        risk.table.y.text = T,
        palette = c("#CC3300", "#3300CC"),
        xlab = "Disease-Free survival time(month)"
    )

    p$plot <- p$plot + labs(caption = paste0("C-idnex: ", c_index))

    return_list <- list()
    return_list[["p"]] <- p
    return_list[["logrank_p"]] <- logrank_pvalue
    return_list[["c-index"]] <- c_index

    return(return_list)
}