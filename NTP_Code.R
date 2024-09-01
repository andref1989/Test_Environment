#### Function for identifying the most variable genes across a set of patients/cell lines.
### mat = expression matrix
### margin = whether individual genes are across the rows or columns
### val_cutoff An option to use an absolute value cutoff instead of variance
### Whether to replace NA values with zero (or some other value)
get_most_variable <- function(mat, margin=c('row','column'),val_cutoff=NULL,quantile=NULL,verbose=F,na_replace=0){
    if(any(is.na(mat)) == TRUE){
        mat[is.na(mat)] <- na_replace
        if(verbose){print("There are NA values in the matrix. Working anyway")}}

    else{ if(verbose){print("Working")}}


    if(is.null(quantile) & !is.null(val_cutoff)){
        cutoff= val_cutoff
    } else if(is.null(quantile) & is.null(val_cutoff)) { cutoff = 0} else if(!is.null(quantile) & is.null(val_cutoff)){ cutoff = quantile
                                                                   } else { print("Can't have an absolute val_cutoff *AND* quantile-based cutoff; setting to 0"); cutoff =0}

    if(verbose){print(cutoff)}
    if( margin =="row"){
        vec <- apply(mat, 1, var)
        vec <- sort(vec,decreasing=TRUE)

    }
    else if( margin== "column"){
        vec <- apply(mat,2,var)
    }

    if(!is.null(val_cutoff)){ variable <- which(vec >= cutoff)} else { variable <- which(vec > quantile(vec,cutoff,na.rm=TRUE))}

    return(variable)}




#### Function for assessing which group a given patient or cell line belongs to using either a reference expression dataset and labels or a named list of genes
### reference_df = reference expression matrix for identifying group-representative genes
### alt_df = the expression matrix being used for predicting group membership. (Ideally you'd do some normalization of the data"
### tsne_df = group membership dataframe for use with the reference expression dataset
### is_ranked = whether the values in either the reference or test dataset is using rank values instead of numeric.
### quantile cutoff = the cutoff to use to identify the top 1-quantile_cutoff variable genes
### method used to identify group-representative genes if none are provided
### random_sample = how many random genesets to use to assess statistical significance.
### variable length = whether to allow group_representative genesets to be of different lengths (or not)
### zscore cutoff = value for determining which genes are included in the genesets (if none are provided)
### templates = named list of group representative genes.
nearest_template_prediction_v3 <- function(reference_df=NULL,alt_df, tsne_df=NULL, cluster_column="Cluster",is_ranked=FALSE,quantile_cutoff=0.5, method=c("simple","randomforest"),template_length=NULL,verbose=TRUE,random_sample=100,variable_length=TRUE,zscore_cutoff=1.4,templates=NULL,sample_col="Sample") {
    require(lsa)
    require(dplyr)



    if(!is.null(reference_df)){tsne_df <- dplyr::filter(tsne_df, get(cluster_column) !="")

        predicted_template_df <- data.frame(stringsAsFactors=F)
        template_df <- unique(tsne_df[unique(c(cluster_column))])
        template_df <- apply(template_df,2, as.character)


    } else{ print("Not using a reference")
        predicted_template_df <- data.frame(stringsAsFactors=F)}

##    gene_set <- unique(unlist(lapply(unique(tsne_df[,TF_column]), function(x) unlist(strsplit(x,"-")))))
    gene_set <- rownames(alt_df)


    if(is.null(templates) & !is.null(reference_df)){
        templates <- get_template_genes_v2(reference_df, tsne_df,cluster_column,quantile_cutoff=quantile_cutoff,method=method, template_length=template_length, is_ranked=is_ranked,value_returned="Templates",variable_length=variable_length, zscore_cutoff=zscore_cutoff)
    } else if(!is.null(templates)) { templates <- lapply(templates, function(x) intersect(x, rownames(alt_df)))}
    else{ print("Need to provide either template genes and/or a reference dataset/groups")
    stop()}

    if(!is.null(reference_df)){
    if(!all(unique(unlist(templates)) %in% union(rownames(reference_df), rownames(alt_df)))){ print("There are template genes missing from the query dataframe. Dropping missing genes")
    templates <- lapply(templates, function(x) intersect(x, intersect(rownames(reference_df),rownames(alt_df))))} else{ templates <- templates}
    print("Finished templates")


###########

    samples <- lapply(sort(unique(template_df[,cluster_column])), function(x) dplyr::filter(tsne_df, get(cluster_column)==x)[,sample_col])
    names(samples) <- sort(unique(template_df[,cluster_column]))
    samples <- lapply(samples, function(x) intersect(x, colnames(reference_df)))



    sample_template_list <- lapply(names(samples), function(x) rowMeans(reference_df[templates[[x]],samples[[x]]],na.rm=TRUE))
    names(sample_template_list) <- names(samples)

    } else if(is.null(reference_df)){

                                                                                                                                  templates2 <- lapply(templates, function(x) unique(c(intersect(x, rownames(alt_df)), setdiff(x, x=unlist(templates)))))
                                                                                                                                  names(templates2) <- names(templates)

                                                                                                                                  sample_template_list <- list()
                                                                                                                                  for(i in names(templates)){
                                                                                                                                      vec <- templates2[[i]]
                                                                                                                                      index <- which(vec %in% templates[[i]])
                                                                                                                                      out_vec <- rep(0,length(vec)); names(out_vec) <- vec
                                                                                                                                      out_vec[index] <- 1
                                                                                                                                      sample_template_list[[i]] <- out_vec}
                                                                                                                                  templates <- templates2}



##    str(names(templates))
##    str(templates)
##    str(sample_template_list)


#######################
#######################

    for(i in colnames(alt_df)){

        if(verbose){
        print(i)}
        distance_list <- unlist(lapply(sort(unique(names(templates))), function(x) 1-cosine(alt_df[templates[[x]],i],sample_template_list[[x]])))
###        str(distance_list)

        names(distance_list) <- sort(unique(names(templates)))


        random_list <- list()
        for(k in names(templates)){
            random_list[[k]] <- unlist(lapply(1:random_sample, function(x) 1-cosine(alt_df[sample(gene_set, length(sample_template_list[[k]])),i],sample_template_list[[k]])))
            random_list[[k]] <- random_list[[k]][!is.na(random_list[[k]])]



        }



        pval_list <- unlist(lapply(names(distance_list), function(x) 1-pnorm(distance_list[[x]],mean(random_list[[x]]),sd(random_list[[x]]),lower.tail=F)))
##        str(pval_list)
        names(pval_list) <- names(distance_list)

        candidate_template_df <- data.frame(pval_list,distance_list,names(distance_list),i,stringsAsFactors=F)
        colnames(candidate_template_df) <- c("P_val","Distance","Template","Cell_Line")
        candidate_template_df$Q_val <- p.adjust(candidate_template_df$P_val,method="BH")
        candidate_template_df$Signif <- candidate_template_df$Q_val <= 0.1
        sub_df <- dplyr::filter(candidate_template_df, Signif==TRUE)

        if(nrow(sub_df) >=1){

            candidate_template_df$Prediction <- sub_df[which(sub_df$Distance == min(sub_df$Distance)),"Template"]} else if (nrow(sub_df) == 0){

                candidate_template_df$Prediction <- candidate_template_df[which(candidate_template_df$Distance == min(candidate_template_df$Distance,na.rm=TRUE)),"Template"]}



        predicted_template_df <- rbind(predicted_template_df,candidate_template_df)

        if(verbose==TRUE){
            print(paste0("Finished ", grep(paste0("^",i,"$"), colnames(alt_df)), " of ",ncol(alt_df)))}

    }

    output <- list(templates,predicted_template_df,sample_template_list)
    names(output) <- c("Templates","Predictions","Reference_Templates")
    print("Finished")
        return(output)
}


### Function for converting the output of prior function into a more readable dataframe.
### NTP_list = a list or single dataframe corresponding to the "Predictions" output from prior function.
### metadata_df = a metadata file describing the cell lines or samples (optional)
### signif_cutoff = Statistical significance threshold
### is_patient = either T/F or skip. Impacts whether any additional metadata columns are considered. Default is to skip
parse_NTP <- function(NTP_list, metadata_df=NULL,signif_cutoff=0.1,is_patient="skip"){
    require(dplyr)

    if(is.data.frame(NTP_list)){

        NTP_list <- list(NTP_list)

        names(NTP_list) <- "NTP" } else{
            NTP_list <- NTP_list
    }

    for(i in names(NTP_list)){
        sub <- NTP_list[[i]]
        ##str(sub)
        sub$Signif <- sub$Q_val <= signif_cutoff
        sub_df <- dplyr::filter(sub, Signif==TRUE)
        for(k in unique(sub$Cell_Line)){
##            str(k)
            index <- which(sub$Cell_Line==k)
            int <- dplyr::filter(sub_df, Cell_Line==k)
##            str(int)
            if(nrow(int) >=1){
            ##print(int[which(int$Distance == min(int$Distance)),"Template"] == unique(sub[index,"Prediction"]))
            sub[index,"Prediction"] <- int[which(int$Distance == min(int$Distance)),"Template"]} else if (nrow(int) == 0){
#                print("Changing")
                sub[index,"Prediction"] <- sub[index,"Prediction"]




                                                                                               }}

        if(is_patient==FALSE){
        sub$Subtype <- metadata_df[sub$Cell_Line,"TCGA"]
        sub$Detailed_Cancer <- metadata_df[sub$Cell_Line,"lineage_subtype"]
        sub$Detailed_Cancer_Ext <- metadata_df[sub$Cell_Line,"disease_subtype"]
        sub$Primary <- metadata_df[sub$Cell_Line,"primary_or_metastasis"]
        sub <- sub %>% group_by(Cell_Line, Prediction) %>% mutate(Accurate= grepl(unique(Subtype), Prediction)) %>% data.frame
        sub$Method <- i
        ##str(sub)
        sub <- sub %>% group_by(Cell_Line) %>% mutate(Num_Predictions=sum(as.numeric(Signif)),Signif=ifelse(any(Signif)==TRUE,TRUE,FALSE)) %>% data.frame
        sub$Q_val_cutoff <- signif_cutoff
        } else if(is_patient==TRUE){
            sub$Subtype <- metadata_df[sub$Cell_Line,"Subtype"]
            sub$Detailed_Cancer <- "Unknown"
            sub$Detailed_Cancer_Ext <- "Unknown"
            sub$Primary <- "Unknown"
            sub <- sub %>% group_by(Cell_Line, Prediction) %>% mutate(Accurate= grepl(unique(Subtype), Prediction)) %>% data.frame
        sub$Method <- i
        ##str(sub)
        sub <- sub %>% group_by(Cell_Line) %>% mutate(Num_Predictions=sum(as.numeric(Signif)),Signif=ifelse(any(Signif)==TRUE,TRUE,FALSE)) %>% data.frame
        sub$Q_val_cutoff <- signif_cutoff

                                                                         } else if(is_patient =="skip"){
                                                                             sub <- sub %>% group_by(Cell_Line) %>% mutate(Num_Predictions=sum(as.numeric(Signif)),Signif=ifelse(any(Signif)==TRUE,TRUE,FALSE)) %>% data.frame
                                                                             sub$Q_val_cutoff <- signif_cutoff
                                                                             }

        NTP_list[[i]] <- sub
    }
    if(is_patient =="skip"){ print("Skipping accuracy assessment")
            output <- do.call("rbind",lapply(NTP_list, function(x) unique(x[c(4,6:ncol(x))])))
            final <- list(NTP_list,output)
            return(final)
            stop()} else{
##    str(NTP_list)
##str(NTP_list)

#print(colnames(output))

    output <- output%>% group_by(Subtype, Method, Primary) %>% mutate(Accuracy_fraction=sum(as.numeric(Accurate))/length(Accurate)) %>% data.frame
output <- output%>% group_by(Subtype, Method) %>% mutate(Accuracy_fraction_overall=sum(as.numeric(Accurate))/length(Accurate)) %>% data.frame
output <- output%>% group_by(Subtype, Method,Signif) %>% mutate(Accuracy_fraction_signif=sum(as.numeric(Accurate))/length(Accurate)) %>% data.frame
output$Primary <- gsub("^$","Unknown",output$Primary)
##output <- as.data.frame(output, stringsAsFactors=F)
##str(output)
##str(NTP_list)
final <- list(NTP_list,output)
##str(final)
    return(list(NTP_list,output))
        }}


### Function to normalize expression across a given dataset.
zscore_df <- function(df, margin=c("row","column"),central_tendency=c("mean","median","median_deviation")){
    if(length(margin) >1){ margin <- "row"}
    if(length(central_tendency>1)){central_tendency <- "mean"}
    if(margin == "row"){
        if(central_tendency=="mean"){df <- t(apply(df,1,function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)))} else if(central_tendency=="median"){ df <- t(apply(df,1,function(x) (x-median(x,na.rm=TRUE))/sd(x,na.rm=TRUE)))} else if (central_tendency=="median_deviation"){ df <- t(apply(df,1,function(x) (x-median(x,na.rm=TRUE))))}

    } else if(margin=="column") { if(central_tendency=="mean"){df <- apply(df,2,function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))} else if (central_tendency=="median"){ df <- apply(df,2,function(x) (x-median(x,na.rm=TRUE))/sd(x,na.rm=TRUE))} else if (central_tendency=="median_deviation"){ df <- apply(df,2,function(x) (x-median(x,na.rm=TRUE)))}
                                                                                                                                                                                                                                                         }
    df[is.na(df)] <- 0
    df <- as.data.frame(df)

    return(df)}



### Function to get group representative genesets using an expression matrix and a metadata/group membership dataframe.

### expresion_df = reference expression matrix for identifying group-representative genes
### tsne_df = group membership dataframe for use with the reference expression dataset
### cluster_column = column where the group membership data is stored
### is_ranked = whether the values in either the reference or test dataset is using rank values instead of numeric.
### quantile cutoff = the cutoff to use to identify the top 1-quantile_cutoff variable genes
### method used to identify group-representative genes if none are provided
### value returned = whether to return the actual genesets or an importance matrix (only relevant for the random forest approach
### variable length = whether to allow group_representative genesets to be of different lengths (or not)
### zscore cutoff = value for determining which genes are included in the genesets (if none are provided)
### sig_limit = maximum number of genes to be considered for a given geneset.

get_template_genes_v2 <- function(expression_df, tsne_df, cluster_column="Cluster",is_ranked=FALSE,quantile_cutoff=0.5, method=c("simple","randomforest"),template_length=50,value_returned=c("Templates","Matrix"),variable_length=TRUE,zscore_cutoff=1.8,abs_val=FALSE,sig_limit=1000){

    if(length(method)>1){ method <- "simple"} else{ method <- method}
    if(length(value_returned)>1){ method <- "Templates"} else { method <- method}
    require(dplyr)


    sample_column <- grep("sample", colnames(tsne_df),ignore.case=TRUE, value=TRUE)
    tsne_df$Sample <- tsne_df[,sample_column]
    tsne_df <- tsne_df[which(tsne_df$Sample %in% colnames(expression_df)),]


    sub_df <- expression_df[,intersect(tsne_df$Sample, colnames(expression_df))]
    sub_df <- sub_df[names(get_most_variable(sub_df, "row",quantile=quantile_cutoff)),]

    if(method == "simple"){
        df <- zscore_df(do.call("cbind", lapply(unique(tsne_df[,cluster_column]), function(x) rowMeans(sub_df[,dplyr::filter(tsne_df, get(cluster_column)==x)$Sample]))),"row")
##        str(df)

        df_abs <- do.call("cbind", lapply(unique(tsne_df[,cluster_column]), function(x) rowMeans(sub_df[,dplyr::filter(tsne_df, get(cluster_column)==x)$Sample])))
        colnames(df) <- unique(tsne_df[,cluster_column])
        colnames(df_abs) <- unique(tsne_df[,cluster_column])
        if(is_ranked == TRUE){
            print("Working 1")
            if(variable_length==FALSE){
                templates_rel <- lapply(colnames(df), function(x) names(sort(df[,x], decreasing=FALSE))[1:sig_limit])
                templates_abs <- lapply(colnames(df_abs), function(x) names(sort(df_abs[,x], decreasing=FALSE))[1:sig_limit])
                names(templates_abs) <- colnames(df)
                names(templates_rel) <- colnames(df)
                templates <- lapply(names(templates_abs), function(x) intersect(templates_rel[[x]], templates_abs[[x]])[1:template_length])} else if(variable_length==TRUE){

                    if(abs_val ==TRUE){ df <- abs(df)*-1


                                    } else{ df <- df}
                    templates_rel <- lapply(colnames(df), function(x) names(sort(df[which(df[,x]<=zscore_cutoff),x], decreasing=FALSE))[1:sig_limit])
                templates_abs <- lapply(colnames(df_abs), function(x) names(sort(df_abs[,x], decreasing=FALSE))[1:sig_limit])
                names(templates_abs) <- colnames(df)
                names(templates_rel) <- colnames(df)
                templates <- lapply(names(templates_abs), function(x) intersect(templates_rel[[x]], templates_abs[[x]]))}

            names(templates) <- names(templates_abs)



        } else if (is_ranked==FALSE){
            print("Working 2")
            if(variable_length==FALSE){

                templates_rel <- lapply(colnames(df), function(x) names(sort(df[,x], decreasing=TRUE))[1:sig_limit])
                templates_abs <- lapply(colnames(df_abs), function(x) names(sort(df_abs[,x], decreasing=TRUE))[1:sig_limit])
                names(templates_abs) <- colnames(df)
                names(templates_rel) <- colnames(df)
                templates <- lapply(names(templates_abs), function(x) intersect(templates_rel[[x]], templates_abs[[x]])[1:template_length])
            } else if(variable_length==TRUE){
                                if(abs_val ==TRUE){ df <- abs(df)} else{ df <- df}
                                          templates_rel <- lapply(colnames(df), function(x) names(sort(df[which(df[,x]>=zscore_cutoff),x], decreasing=TRUE))[1:sig_limit])
                templates_abs <- lapply(colnames(df_abs), function(x) names(sort(df_abs[,x], decreasing=TRUE))[1:sig_limit])
                names(templates_abs) <- colnames(df)
                names(templates_rel) <- colnames(df)
                templates <- lapply(names(templates_abs), function(x) intersect(templates_rel[[x]], templates_abs[[x]]))




                                      }
        names(templates) <- colnames(df)

                                  }}

    else if(method == "randomforest"){
                                          require(randomForest)
        rownames(sub_df) <- gsub("-","_",rownames(sub_df))
        ##str(sub_df)
##        sub_df <- sub_df[names(get_most_variable(sub_df, "row",quantile=quantile_cutoff)),]
    transposed_df <- as.data.frame(t(sub_df))
    transposed_df[,cluster_column] <- as.factor(tsne_df[colnames(expression_df),cluster_column])
    uniq_classes <- unique(tsne_df[,cluster_column])
##    print(index)
##   print(cluster_column)


##    str(transposed_df)
    print("Starting Random Forest")
    rf <- randomForest(as.formula(paste0(cluster_column,"~.")),data=transposed_df,importance=TRUE)
        importance_df <- as.data.frame(importance(rf))[setdiff(colnames(importance(rf)),uniq_classes)]
##        str(importance_df)

    importance_df <- importance_df[order(-importance_df[,2]),]
    final_importance_matrix <- data.frame(importance(rf)[rownames(importance_df),uniq_classes],stringsAsFactors=F)


    colnames(final_importance_matrix) <- uniq_classes
##        str(final_importance_matrix)
##        str(colnames(final_importance_matrix))

    templates_abs <- lapply(colnames(final_importance_matrix), function(x) rownames(sort_df(final_importance_matrix,x, sort_order=TRUE))[1:sig_limit])

    templates_rel <- lapply(colnames(final_importance_matrix), function(x) rownames(sort_df(zscore_df(final_importance_matrix),x,sort_order=TRUE))[1:sig_limit])
#    str(templates_rel)
        ##str(final_importance_matrix)
    names(templates_abs) <- colnames(final_importance_matrix)
    names(templates_rel) <- colnames(final_importance_matrix)

  ##      str(templates_abs)
  ##      str(templates_rel)
     if(variable_length==FALSE){
                templates <- lapply(names(templates_abs), function(x) intersect(templates_rel[[x]], templates_abs[[x]])[1:template_length])
            } else if(variable_length==TRUE){
                zscore_importance_matrix <- zscore_df(final_importance_matrix,"row",central_tendency="median")
                templates_rel <- lapply(colnames(zscore_importance_matrix), function(x) names(sort(zscore_importance_matrix[which(zscore_importance_matrix[,x]>=zscore_cutoff),x], decreasing=TRUE))[1:sig_limit])
                names(templates_rel) <- colnames(zscore_importance_matrix)
##                str(templates_rel)
##                str(templates_abs)
                templates <- lapply(names(templates_abs), function(x) intersect(templates_rel[[x]], templates_abs[[x]]))      }


#    templates <- lapply(colnames(final_importance_matrix), function(x) rownames(sort_df(final_importance_matrix,x, decreasing=TRUE))[1:template_length])
    templates <- lapply(templates, function(x) gsub("_","-",x))
        names(templates) <- names(templates_abs)


}
#    str(templates)
        templates <- lapply(templates, function(x) x[!is.na(x)])
    if(value_returned=="Templates"){ return(templates)} else if(value_returned=="Matrix"){ return(final_importance_matrix)}

}
