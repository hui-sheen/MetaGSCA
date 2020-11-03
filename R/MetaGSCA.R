
#' @export
#############################
## 1.Function: GSAR_boot ####
## group: patient subgroup or condition (values: 1 or 2)
## genelist: a pre-specified gene list or signature or defined by pathway
## object: gene-expression matrix wehre rows represent genes and columns represent patient/sample. The first column should be gene name.

GSAR_boot <- function(projectname, gematrix, group, genelist, nperm = 500, cor.method = "pearson",
                      check.sd = TRUE, min.sd = 0.001, max.skip =50, R = 200, level = 0.05){

    call <- match.call()

    # in case some genes in the provided gene list cannot be found in the dataset
    genes.not.found <- paste(genelist[!(genelist %in% as.character(gematrix[,1]))], collapse = ",")
    if(sum(!(genelist %in% as.character(gematrix[,1]))) == 0) genes.not.found <- "-"
    #current.time <- Sys.time()
    #print(paste0(current.time, " Genes not found: ", genes.not.found))

    mat<-gematrix[as.character(gematrix[,1]) %in% genelist, ]
    object <- as.matrix(mat[,2:ncol(mat)])  # be careful: starting with which column? Might be 4. Also the first column should be gene name

    if (!(is.matrix(object)))
        stop("'object' must be a matrix where rows are features and columns are samples")
    if (is.null(group))
        stop("'group' must be a vector indicating group association. Possible values are 1 and 2")

    nv <- ncol(object)
    group <- as.numeric(group)
    if (length(group) != nv)
        stop("length of 'group' must equal the number of columns in 'object'")
    if (sum(group %in% c(1, 2)) < nv)
        stop("all members in 'group' must have values 1 or 2")
    if ((sum(group == 1) < 3) || (sum(group == 2) < 3))
        stop("there are less than 3 samples in at least one group")
    # if (!(cor.method %in% c("pearson", "spearman", "kendall")))
    #   stop("'cor.method' must be a character string indicating which \n correlation coefficient to be calculated. \n One of 'pearson' (default), 'spearman' or 'kendall'")
    #

    ### 1, Original Data
    object1 <- object[, c(which(group == 1))]; objt1 <- aperm(object1, c(2, 1))
    object2 <- object[, c(which(group == 2))]; objt2 <- aperm(object2, c(2, 1))

    if (check.sd == TRUE) {
        sd1 <- apply(objt1, 2, "sd")
        sd2 <- apply(objt2, 2, "sd")

        delcol <- (sd1 < min.sd) | (sd2 < min.sd)
        genes.removed <- paste(as.character(mat[,1])[delcol], collapse = ",") ## genes removed cause sd less than min.sd
        if(sum(delcol) == 0) genes.removed <- "-"
        #current.time <- Sys.time()
        #print(paste0(current.time, " Genes removed(sd<", min.sd, "): ", genes.removed))

        if(sum(!delcol) < 2)
            stop(paste0("There must be at least 2 genes with standard deviation of gene expression data bigger than ", min.sd, " in group 1 and group 2"))

        objt1 <- objt1[, !delcol]
        objt2 <- objt2[, !delcol]
        object <- object[!delcol,]
        object1 <- object1[!delcol,]
        object2 <- object2[!delcol,]

        sd1 <- apply(objt1, 2, "sd")
        sd2 <- apply(objt2, 2, "sd")

        if (sum(sd1 < min.sd, na.rm=TRUE) > 0)
            stop(paste("There are ", sum(sd1 < min.sd), " Original Data: features with standard deviation smaller than ", min.sd, " in group 1", sep = ""))
        if (sum(sd2 < min.sd, na.rm=TRUE) > 0)
            stop(paste("There are ", sum(sd2 < min.sd), " Original Data: features with standard deviation smaller than ", min.sd, " in group 2", sep = ""))
    }

    cormat1 <- abs(cor(objt1, method = cor.method))
    cormat2 <- abs(cor(objt2, method = cor.method))
    e1 <- eigen(cormat1)        ### eigen: Computes eigenvalues and eigenvectors of numeric (double, integer, logical) or complex matrices.
    e2 <- eigen(cormat2)
    p1 <- abs(e1$vectors[, 1])
    p2 <- abs(e2$vectors[, 1])
    D_obs <- sum(abs((p1 * norm(matrix(p1))) - (p2 * norm(matrix(p2)))))       ### Compute the Norm of a Matrix

    domain <- seq_len(ncol(object))
    D_perm <- array(0, c(1, nperm))
    skip.counter <- 0

    for (itr in seq_len(nperm)) {
        randperm <- sample(domain, replace = FALSE)      ### Random Samples and Permutations
        objt <- aperm(object[, randperm], c(2, 1))
        nv1 <- sum(group == 1)
        group1 <- objt[seq_len(nv1), ]
        group2 <- objt[(nv1 + 1):nv, ]
        if (check.sd == TRUE) {
            sd1 <- apply(group1, 2, "sd")
            sd2 <- apply(group2, 2, "sd")
            while (((sum(sd1 < min.sd) > 0) || (sum(sd2 < min.sd) > 0)) && (skip.counter <= max.skip)) {
                if (skip.counter == max.skip){
                    stop("number of skipped permutations exceeded 'max.skip'")
                }

                skip.counter <- skip.counter + 1
                randperm <- sample(domain, replace = FALSE)
                objt <- aperm(object[, randperm], c(2, 1))
                group1 <- objt[seq_len(nv1), ]
                group2 <- objt[(nv1 + 1):nv, ]
                sd1 <- apply(group1, 2, "sd")
                sd2 <- apply(group2, 2, "sd")
            }
        }
        cormat1 <- abs(cor(group1, method = cor.method))
        cormat2 <- abs(cor(group2, method = cor.method))
        e1 <- eigen(cormat1)
        e2 <- eigen(cormat2)
        p1 <- abs(e1$vectors[, 1])
        p2 <- abs(e2$vectors[, 1])
        D_perm[itr] <- sum(abs((p1 * norm(matrix(p1))) - (p2 * norm(matrix(p2)))))
    }

    pvalue <- (sum(D_perm >= D_obs) + 1)/(length(D_perm) + 1)


    ### 2, Bootstrap Sample

    n1 <- nrow(objt1); n2 <- nrow(objt2)
    D_obs_boot <- c()
    pvalue_boot <- c()


    for (i in seq_len(R)){

        #print(i)

        output1 <- sample.int(n1, n1, replace=TRUE)
        output2 <- sample.int(n2, n2, replace=TRUE)
        object_boot1 <- object1[, output1]; objt_boot1 <- aperm(object_boot1, c(2, 1))
        object_boot2 <- object2[, output2]; objt_boot2 <- aperm(object_boot2, c(2, 1))
        object_boot <- cbind(object_boot1, object_boot2)

        if (check.sd == TRUE) {
            sd1 <- apply(objt_boot1, 2, "sd")
            sd2 <- apply(objt_boot2, 2, "sd")
            while (((sum(sd1 < min.sd) > 0) || (sum(sd2 < min.sd) > 0)) && (skip.counter <= max.skip)) {
                if (skip.counter == max.skip)
                    stop("number of skipped bootstraps exceeded 'max.skip'")

                output1 <- sample.int(n1, n1, replace=TRUE)
                output2 <- sample.int(n2, n2, replace=TRUE)
                object_boot1 <- object1[, output1]; objt_boot1 <- aperm(object_boot1, c(2, 1))
                object_boot2 <- object2[, output2]; objt_boot2 <- aperm(object_boot2, c(2, 1))
                object_boot <- cbind(object_boot1, object_boot2)
                sd1 <- apply(objt_boot1, 2, "sd")
                sd2 <- apply(objt_boot2, 2, "sd")
            }
        }

        cormat1 <- abs(cor(objt_boot1, method = cor.method))
        cormat2 <- abs(cor(objt_boot2, method = cor.method))
        e1 <- eigen(cormat1)        ### eigen: Computes eigenvalues and eigenvectors of numeric (double, integer, logical) or complex matrices.
        e2 <- eigen(cormat2)
        p1 <- abs(e1$vectors[, 1])
        p2 <- abs(e2$vectors[, 1])
        D_obs_boot[i] <- sum(abs((p1 * norm(matrix(p1))) - (p2 * norm(matrix(p2)))))       ### Compute the Norm of a Matrix

        domain <- seq_len(ncol(object))
        D_perm_boot <- array(0, c(1, nperm))
        skip.counter <- 0

        for (itr in seq_len(nperm)) {
            randperm <- sample(domain, replace = FALSE)      ### Random Samples and Permutations
            objt <- aperm(object[, randperm], c(2, 1))
            nv1 <- sum(group == 1)
            group1 <- objt[seq_len(nv1), ]
            group2 <- objt[(nv1 + 1):nv, ]
            if (check.sd == TRUE) {
                sd1 <- apply(group1, 2, "sd")
                sd2 <- apply(group2, 2, "sd")
                while (((sum(sd1 < min.sd) > 0) || (sum(sd2 < min.sd) > 0)) && (skip.counter <= max.skip)) {
                    if (skip.counter == max.skip)
                        stop("number of skipped permutations exceeded 'max.skip'")

                    skip.counter <- skip.counter + 1
                    randperm <- sample(domain, replace = FALSE)
                    objt <- aperm(object[, randperm], c(2, 1))
                    group1 <- objt[seq_len(nv1), ]
                    group2 <- objt[(nv1 + 1):nv, ]
                    sd1 <- apply(group1, 2, "sd")
                    sd2 <- apply(group2, 2, "sd")
                }
            }
            cormat1 <- abs(cor(group1, method = cor.method))
            cormat2 <- abs(cor(group2, method = cor.method))
            e1 <- eigen(cormat1)
            e2 <- eigen(cormat2)
            p1 <- abs(e1$vectors[, 1])
            p2 <- abs(e2$vectors[, 1])
            D_perm_boot[itr] <- sum(abs((p1 * norm(matrix(p1))) - (p2 * norm(matrix(p2)))))
        }

        pvalue_boot[i] <- (sum(D_perm_boot >= D_obs_boot[i]) + 1)/(length(D_perm) + 1)


    }


    ## Original Part
    D_obs <- round(D_obs, 6)
    pvalue <- round(pvalue, 6)

    ## Bootstrap Part
    D_obs_boot <- round(D_obs_boot, 6)
    pvalue_boot <- round(pvalue_boot, 6)

    ## 1, Test statistics/ Mean
    ts.mean <- round(mean(D_obs_boot), 6)
    ts.bias <- round((mean(D_obs_boot)-D_obs), 6)
    ts.std <- round(sd(D_obs_boot), 6)
    ts.norm.ci <- paste((1-level)*100, "% CI: [", round(mean(D_obs_boot)-qnorm(1-level/2)*sd(D_obs_boot), 6), ", ", round(mean(D_obs_boot)+qnorm(1-level/2)*sd(D_obs_boot), 6), "]", sep = "")

    ## 2, Test statistics/ Median
    ts.percentile <- round(quantile(D_obs_boot, probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)), 6)
    ts.pt.ci <- paste((1-level)*100, "% CI: [", round(as.numeric(quantile(D_obs_boot,  probs = level/2)), 6), ", ", round(as.numeric(quantile(D_obs_boot,  probs = (1-level/2))), 6), "]",  sep = "")

    ## 3, P value/ Mean
    p.mean <- round(mean(pvalue_boot), 6)
    p.bias <- round((mean(pvalue_boot)-pvalue), 6)
    p.std <- round(sd(pvalue_boot), 6)
    p.norm.ci <- paste((1-level)*100, "% CI: [", round(mean(pvalue_boot)-qnorm(1-level/2)*sd(pvalue_boot), 6), ", ", round(mean(pvalue_boot)+qnorm(1-level/2)*sd(pvalue_boot), 6), "]", sep = "")

    ## 4, P value/ Meidan
    p.percentile <- round(quantile(pvalue_boot, probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)), 6)
    p.pt.ci <- paste((1-level)*100, "% CI: [", round(as.numeric(quantile(pvalue_boot,  probs = level/2)), 6), ", ", round(as.numeric(quantile(pvalue_boot,  probs = (1-level/2))), 6), "]",  sep = "")

    ## Combined matrix
    Coexpression <- data.frame(Original.TS = D_obs,
                               Original.Pvalue = pvalue,
                               Boot.TS.Mean = ts.mean,
                               Boot.TS.Bias = ts.bias,
                               Boot.TS.Std = ts.std,
                               Boot.TS.Norm.CI = ts.norm.ci,
                               Boot.TS.PT.CI = ts.pt.ci,
                               Boot.Pvalue.Mean = p.mean,
                               Boot.Pvalue.Bias = p.bias,
                               Boot.Pvalue.Std = p.std,
                               Boot.Pvalue.Norm.CI = p.norm.ci,
                               Boot.Pvalue.PT.CI = p.pt.ci,
                               Genes.not.found = genes.not.found,
                               Genes.removed = genes.removed)
    colnames(Coexpression)[which(colnames(Coexpression) == "Genes.removed")] <- paste0("Genes.removed", "(sd<", min.sd, ")")

    ## this is only for the report
    #write.csv(Coexpression, paste("Coexpression_", projectname, ".csv"))

    output<-list(Original.TS = D_obs,
                 Original.Pvalue = pvalue,
                 Boot.TS = D_obs_boot,
                 Boot.Pvalue = pvalue_boot,
                 Boot.TS.Mean = ts.mean,
                 Boot.TS.Bias = ts.bias,
                 Boot.TS.Std = ts.std,
                 Boot.TS.Norm.CI = ts.norm.ci,
                 Boot.TS.Percentile = ts.percentile,
                 Boot.TS.PT.CI = ts.pt.ci,
                 Boot.Pvalue.Mean = p.mean,
                 Boot.Pvalue.Bias = p.bias,
                 Boot.Pvalue.Std = p.std,
                 Boot.Pvalue.Norm.CI = p.norm.ci,
                 Boot.Pvalue.Percentile = p.percentile,
                 Boot.Pvalue.PT.CI = p.pt.ci,
                 call = call,
                 nperm = nperm,
                 R = R,
                 Coexpression = Coexpression)

    return(output)
}








#' @export
############################
### 2.Function: MetaGSCA ###
### single gene set
MetaGSCA <- function(list.geneset,  ## a pre-specified gene list from a gene set or pathway
                     list.dataset,  ## a list of datasets, first column is gene name
                     list.group,  ## a list of samples/patients subgroup or condition (e.g. (1,1,1,2,2,2))
                     name.geneset = "Geneset",  ## the name of the gene set, used for output file name
                     name.dataset,  ## a list of dataset names corresponding to list.dataset, used for forest plot
                     nperm = 500,
                     nboot = 200,
                     method.Inverse = FALSE,
                     method.GLMM = TRUE,
                     fixed.effect = FALSE,  ## A logical indicating whether a fixed effect meta-analysis should be conducted.
                     random.effect = TRUE  ## A logical indicating whether a random effects meta-analysis should be conducted.
){
    if(!xor(method.Inverse, method.GLMM))
        stop("Only one of Inverse method and GLMM method can be TRUE")
    if(!xor(fixed.effect, random.effect))
        stop("Only one of random effect and fixed effect can be TRUE")
    if (length(list.group) != length(list.dataset))
        stop("The length of the 'list.dataset' and the length of the 'list.group' must be equal")
    if (length(name.dataset) != length(list.dataset))
        stop("The length of the 'list.dataset' and the length of the 'name.dataset' must be equal")


    ### GSAR_boot Analysis
    list <- list()
    res <- NULL
    coexp <- NULL
    for (j in seq_len(length(list.dataset))){
        current.time <- Sys.time()
        cat(paste0(current.time, " Processing dataset: ", name.dataset[j], "..."))
        cat('\n')

        L <- GSAR_boot(projectname=name.dataset[j], gematrix=list.dataset[[j]], group=list.group[[j]], genelist=list.geneset, nperm=nperm, R=nboot)
        list[[j]] <- L
        coexp <- L$Coexpression
        rownames(coexp) <- name.dataset[j]
        res<-rbind(res, coexp)
    }
    write.csv(res, paste0(name.geneset, "_Individual Dataset GSAR Results with Bootstrapping.csv"))

    ### Meta Analysis
    Event <- c()
    for (j in seq_len(length(list))){
        event <- ceiling(list[[j]]$Original.Pvalue * (nperm + 1) - 1)
        Event <- c(Event, event)
        #print(event)
    }
    N <- rep(nperm, length(list))


    ### APPROACH 1: INVERSE VARIANCE ###
    if(method.Inverse == TRUE){
        ## STEP 1, Meta Analysis for original p-value
        meta.inv <- metaprop(
            event = Event,
            n = N,
            studlab = name.dataset,
            comb.fixed = fixed.effect,
            comb.random = random.effect,
            prediction = TRUE,
            method = "inverse",
            method.tau = "ML")


        ## STEP 2, Meta Analysis for bootstrap p-value
        meta.p.inv.fixed <- c()
        meta.p.inv.random <- c()
        for (i in seq_len(nboot)) {
            event.boot <- c()
            for(j in seq_len(length(list))){
                event.boot <- c(event.boot, ceiling(list[[j]]$Boot.Pvalue[i] * (nperm + 1) - 1))
            }

            N.boot <- rep(nperm, length(list))

            m.boot <- tryCatch(metaprop(
                event = event.boot,
                n = N.boot,
                studlab = name.dataset,
                comb.fixed = fixed.effect,
                comb.random = random.effect,
                prediction = TRUE,
                method = "inverse",
                method.tau = "ML"),
                error = function(e) return(NA))

            if(is.na(m.boot)){
                meta.p.inv.fixed[i] <- NA
                meta.p.inv.random[i] <- NA
            }else{
                meta.p.inv.fixed[i] <- exp(m.boot$TE.fixed)
                meta.p.inv.random[i] <- exp(m.boot$TE.random)
            }
        }

        pdf(file = paste0(name.geneset, "_Inverse_Meta Analysis for bootstrap p-value_Forest plot.pdf"), width=10, height=3+length(list)*0.25);
        forest(
            meta.inv,
            leftlabs = c("Datasets", "No. of sig", "NPerm"),
            rightlabs = c("P statistic", "95%-CI", "Weight"),
            xlim = c(0, max(meta.inv$upper)+0.1),
            digits = 4,
            col.diamond = "red",
            col.diamond.lines = "black",
            print.tau2 = FALSE,
            print.I2 = FALSE,
            prediction = FALSE,
            col.predict = "black",
            smlab = "")

        if (fixed.effect == TRUE){
            grid.text(label = paste("Fixed effect model bootstrapping result (nperm=", nperm, ", nboot=", nboot, ")", sep=""),
                      x = unit(0.12, "npc"), y = unit(0.15, "npc"),
                      just = "left", gp = gpar(fontsize = 12, fontface = "bold"))

            grid.text(label = paste(format(round(median(meta.p.inv.fixed, na.rm = TRUE),4), nsmall=4), " [",
                                    format(round(quantile(meta.p.inv.fixed, probs=0.025, na.rm = TRUE),4), nsmall=4), ", ",
                                    format(round(quantile(meta.p.inv.fixed, probs=0.975, na.rm = TRUE),4), nsmall=4), "]", sep=""),
                      x = unit(0.70, "npc"), y = unit(0.15, "npc"),
                      just = "left", gp = gpar(fontsize = 12, fontface = "bold"))
        }

        if (random.effect == TRUE){
            grid.text(label = paste("Random effect model bootstrapping result (nperm=", nperm, ", nboot=", nboot, ")", sep=""),
                      x = unit(0.12, "npc"), y = unit(0.15, "npc"),
                      just = "left", gp = gpar(fontsize = 12, fontface = "bold"))

            grid.text(label = paste(format(round(median(meta.p.inv.random, na.rm = TRUE),4), nsmall=4), " [",
                                    format(round(quantile(meta.p.inv.random, probs=0.025, na.rm = TRUE),4), nsmall=4), ", ",
                                    format(round(quantile(meta.p.inv.random, probs=0.975, na.rm = TRUE),4), nsmall=4), "]", sep=""),
                      x = unit(0.70, "npc"), y = unit(0.15, "npc"),
                      just = "left", gp = gpar(fontsize = 12, fontface = "bold"))
        }
        dev.off()
    } ## method.Inverse end


    ### APPROACH 2: GLMM ###
    if(method.GLMM == TRUE){
        ## sTEP 1, Meta Analysis for original p-value
        if (sum(Event) / length(list) == Event[1]) {
            print(length(list))
            meta.glmm <- Event[1] / nperm

            cat('\n')
            cat("---------------------- GLMM Warning--------------------")
            cat('\n')
            cat("If all input datasets have same event, the GLMM approach cannot fit ML model.")
            cat('\n')
            cat("Error in rma.glmm(xi = event[!exclude], ni = n[!exclude], method = method.tau")
            cat('\n')
            cat(paste("P-value (poroportion) = ", meta.glmm))
            cat('\n')

        } else {
            meta.glmm <- metaprop(
                event = Event,
                n = N,
                studlab =  name.dataset,
                comb.fixed = fixed.effect,
                comb.random = random.effect,
                prediction = TRUE,
                method = "GLMM",
                method.tau = "ML")
        }

        ## STEP 2, Meta Analysis for bootstrap p-value
        meta.p.glmm.fixed <- c()
        meta.p.glmm.random <- c()
        for (i in seq_len(nboot)) {
            event.boot <- c()
            for(j in seq_len(length(list))){
                event.boot <- c(event.boot, ceiling(list[[j]]$Boot.Pvalue[i] * (nperm + 1) - 1))
            }

            N.boot <- rep(nperm, length(list))

            if (sum(event.boot) / length(list) == event.boot[1]) {
                meta.p.glmm.fixed[i] <- event.boot[1] / nperm
                meta.p.glmm.random[i] <- event.boot[1] / nperm
            } else {

                m.boot <- tryCatch(metaprop(
                    event = event.boot,
                    n = N.boot,
                    studlab = name.dataset,
                    comb.fixed = fixed.effect,
                    comb.random = random.effect,
                    prediction = TRUE,
                    method = "GLMM",
                    method.tau = "ML"),
                    error = function(e) return(NA))

                if(is.na(m.boot)){
                    meta.p.glmm.fixed[i] <- NA
                    meta.p.glmm.random[i] <- NA
                }else{
                    meta.p.glmm.fixed[i] <- exp(m.boot$TE.fixed)
                    meta.p.glmm.random[i] <- exp(m.boot$TE.random)
                }
            }
        }

        if (is.numeric(meta.glmm)==TRUE){
            print(paste("glmm.fixed bootstrap result (permutation time=", nperm, ", bootstrap time=", nboot, ")",
                        round(median(meta.p.glmm.fixed),4), " [", round(quantile(meta.p.glmm.fixed, probs=0.025),4),
                        ", ", round(quantile(meta.p.glmm.fixed, probs=0.975),4), "]",
                        round(median(meta.p.glmm.random),4), " [", round(quantile(meta.p.glmm.random, probs=0.025),4),
                        ", ", round(quantile(meta.p.glmm.random, probs=0.975),4), "]", sep=""))
        } else if (is.numeric(meta.glmm)==FALSE){
            pdf(file = paste0(name.geneset, "_GLMM_Meta Analysis for bootstrap p-value_Forest plot.pdf"), width=10, height=3+length(list)*0.25);
            forest(
                meta.glmm,
                leftlabs = c("Datasets", "No. of sig", "NPerm"),
                rightlabs = c("P statistic", "95%-CI"),
                xlim = c(0, max(meta.glmm$upper)+0.1),
                digits = 4,
                col.diamond = "red",
                col.diamond.lines = "black",
                print.tau2 = FALSE,
                print.I2 = FALSE,
                prediction = FALSE,
                col.predict = "black",
                smlab = "")

            if (fixed.effect == TRUE){
                grid.text(label = paste("Fixed effect model bootstrapping result (nperm=", nperm, ", nboot=", nboot, ")", sep=""),
                          x = unit(0.12, "npc"), y = unit(0.15, "npc"),
                          just = "left", gp = gpar(fontsize = 12, fontface = "bold"))

                grid.text(label = paste(format(round(median(meta.p.glmm.fixed, na.rm = TRUE),4), nsmall=4), " [",
                                        format(round(quantile(meta.p.glmm.fixed, probs=0.025, na.rm = TRUE),4), nsmall=4), ", ",
                                        format(round(quantile(meta.p.glmm.fixed, probs=0.975, na.rm = TRUE),4), nsmall=4), "]", sep=""),
                          x = unit(0.70, "npc"), y = unit(0.15, "npc"),
                          just = "left", gp = gpar(fontsize = 12, fontface = "bold"))
            }

            if (random.effect == TRUE){
                grid.text(label = paste("Random effect model bootstrapping result (nperm=", nperm, ", nboot=", nboot, ")", sep=""),
                          x = unit(0.12, "npc"), y = unit(0.15, "npc"),
                          just = "left", gp = gpar(fontsize = 12, fontface = "bold"))

                grid.text(label = paste(format(round(median(meta.p.glmm.random, na.rm = TRUE),4), nsmall=4), " [",
                                        format(round(quantile(meta.p.glmm.random, probs=0.025, na.rm = TRUE),4), nsmall=4), ", ",
                                        format(round(quantile(meta.p.glmm.random, probs=0.975, na.rm = TRUE),4), nsmall=4), "]", sep=""),
                          x = unit(0.70, "npc"), y = unit(0.15, "npc"),
                          just = "left", gp = gpar(fontsize = 12, fontface = "bold"))
            }
            dev.off()
        }

    } ## method.GLMM end

    ## bootstrap result
    bootstrap <- list()
    bootstrap$Geneset.Name <- name.geneset
    bootstrap$Number.of.Genes <- length(list.geneset)

    if(method.Inverse == TRUE){
        if(fixed.effect == TRUE){
            inv.fixed <- c(meta.inv$TE.fixed, meta.inv$lower.fixed, meta.inv$upper.fixed)
            bootstrap$meta.p.fixed = meta:::backtransf(inv.fixed, sm="PLOGIT")[1]
            bootstrap$bootstrap.p.fixed = median(meta.p.inv.fixed, na.rm = TRUE)
        }
        if(random.effect == TRUE){
            inv.random <- c(meta.inv$TE.random, meta.inv$lower.random, meta.inv$upper.random)
            bootstrap$meta.p.random = meta:::backtransf(inv.random, sm="PLOGIT")[1]
            bootstrap$bootstrap.p.random = median(meta.p.inv.random, na.rm = TRUE)
        }
    }

    if(method.GLMM == TRUE & is.numeric(meta.glmm) == FALSE){
        if(fixed.effect == TRUE){
            glmm.fixed <- c(meta.glmm$TE.fixed, meta.glmm$lower.fixed, meta.glmm$upper.fixed)
            bootstrap$meta.p.fixed = meta:::backtransf(glmm.fixed, sm="PLOGIT")[1]
            bootstrap$bootstrap.p.fixed = median(meta.p.glmm.fixed, na.rm = TRUE)
        }
        if(random.effect == TRUE){
            glmm.random <- c(meta.glmm$TE.random, meta.glmm$lower.random, meta.glmm$upper.random)
            bootstrap$meta.p.random = meta:::backtransf(glmm.random, sm="PLOGIT")[1]
            bootstrap$bootstrap.p.random = median(meta.p.glmm.random, na.rm = TRUE)
        }
    }

    if(method.GLMM == TRUE & is.numeric(meta.glmm) == TRUE){
        bootstrap$meta.glmm <- meta.glmm
    }

    bootstrap[(length(bootstrap)+1):(length(list.dataset)+length(bootstrap))] <- res$Original.Pvalue
    names(bootstrap)[(length(bootstrap)-length(list.dataset)+1):length(bootstrap)] <- name.dataset

    return(bootstrap)
}






#' @export
##########################################
## 3.Function: MetaGSCA_Multi_Geneset ####
## Multiple gene set Meta-Analysis
MetaGSCA_Multi_Geneset <- function(list.geneset,  ## list of geneset
                                   list.dataset,  ## a list of datasets, first column is gene name
                                   list.group,  ## a list of samples/patients subgroup or condition (e.g. (1,1,1,2,2,2))
                                   name.geneset = "Geneset",  ## the name of the gene set, used for output file name
                                   name.dataset, ## a list of dataset names corresponding to list.dataset, used for forest plot
                                   nperm = 500,
                                   nboot = 200,
                                   method.Inverse = FALSE,
                                   method.GLMM = TRUE,
                                   fixed.effect = FALSE,
                                   random.effect = TRUE){

    res <- NULL
    for (i in seq_len(length(list.geneset))) {
        cat(name.geneset[i])
        cat('\n')

        bootstrap <- tryCatch(MetaGSCA(list.geneset = list.geneset[[i]],
                                       list.dataset = list.dataset,
                                       list.group = list.group,
                                       name.geneset = name.geneset[[i]],
                                       name.dataset = name.dataset,
                                       nperm = nperm,
                                       nboot = nboot,
                                       method.Inverse = method.Inverse,
                                       method.GLMM = method.GLMM,
                                       fixed.effect = fixed.effect,
                                       random.effect = random.effect),
                              error = function(e) return(NULL))
        if(!is.null(bootstrap)) res <- rbind(res, as.data.frame(bootstrap))
    }

    if(!is.null(res))
        write.csv(res, "_Meta Analysis bootstrap result.csv", row.names = FALSE)
}


