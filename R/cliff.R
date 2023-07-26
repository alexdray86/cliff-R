#' CLIFF deconvolution of drug sensitivity data
#'
#' Takes the output of CLIMB deconvolution together with drug sensitivity data as AUC values,
#' and optionnaly somatic mutation data that will be attributed to cancer cell-type containing
#' the pattern defined by 'cancer_pattern' parameter. 
#'
#' @param climb_output CLIMB output object, containing deconvoluted expression and proportions
#' @param drug_data data.frame containing a column named 'auc' with AUC values of drug sensitivity
#' @param mutation_data matrix optional parameter containing a 1/0 matrix with mutation per bulk sample
#' @param min.mutation integer minimum number of patient having a mutation to include it in the model
#' @param em_steps integer maximum number of steps in the EM algorithm
#' @param mode string 'overall' will uses the same average cell-type expression across patients. 'high-res' uses the deconvoluted expression at high-resolution produced by CLIMB
#' @param regularization string 'none' does not apply normalization. L2 applied ridge regularization

sigmoid_ <- function(x, a){
    return( 1 / (1 + exp(-x + a)))
}
rmse <- function(tr, pr){ return( sqrt( sum((tr - pr)^2) / length(tr) ) ) }

cliff <- function(climb_output, drug_data, mutation_data=NULL, min.mutation=0, em_steps = 30, 
                  mode='highres', regularization='none'){
    message('Select overlapping samples')
    climb_expr = climb_output$expr.highres
    climb_expr_overall = climb_output$expr.overall
    climb_prop = climb_output$props
    bulk_mat = climb_output$bulk_mat
    
    mean_auc = mean(drug_data$auc)
    drug_data$auc = drug_data$auc + (0.5 - mean(mean(drug_data$auc)))
    N = dim(bulk_mat)[2]
    if(is.null(mutation_data)){
        message('No mutation data provided')
        mutation_data = matrix(0, ncol=2, nrow=N)
        colnames(mutation_data) = c('a', 'b')
        rownames(mutation_data) = 1:N
    }
    sel.sample = intersect(dimnames(climb_expr)[[1]], intersect(rownames(mutation_data), drug_data$sample))
    order.sample.climb = match(sel.sample, dimnames(climb_expr)[[1]])
    ## Subset Sample to get samples available for a given drug
    drug_data = drug_data[match(sel.sample, drug_data$sample),]
    mutation_data = mutation_data[sel.sample,]
    bulk_mat = bulk_mat[,sel.sample]
    climb_expr = climb_expr[order.sample.climb,,]
    climb_prop = climb_prop[order.sample.climb,]
    # Select mutation that spans at least 5 patients
    sel.mutation = colSums(mutation_data) >= min.mutation
    mutation_data = mutation_data[,sel.mutation]
    mutation_data = as.matrix(mutation_data)
    # Check that we have the same sample names in good order
    stopifnot(all(dimnames(climb_expr)[[1]] == rownames(mutation_data)))
    stopifnot(all(dimnames(climb_expr)[[1]] == drug_data$sample))
    
    ## PREPARE CLIFF INPUT FROM CLIMB OUTPUT
    message('Prepare CLIFF input from CLIMB output, mutation data, and drug sensitivity data')
    K = num(dim(climb_expr)[3]) ; N = num(dim(climb_expr)[1])
    sample_names = dimnames(climb_expr)[[1]]
    tabs_ = list()
    for(n in 1:N){
        if(mode == 'overall'){
            climb_expr[n,,] = climb_expr_overall
        } 
        df_expr = data.frame(t(log2(climb_expr[n,,]+1))/rowSums(t(log2(climb_expr[n,,]+1)))) * climb_prop[n,]
        df_expr[is.na(df_expr)] <- 0
        #df_expr = data.frame(mean.expr.sc.es) - If one wants to use average cell-type expression
        rownames(df_expr) = paste0(rownames(df_expr),'_', sample_names[n])
        df_expr$props = climb_prop[n,]
        mat.mut = matrix(0, ncol=length(mutation_data[n,]), nrow=length(colnames(climb_prop)))
        mat.mut[grepl(cancer_pattern,colnames(climb_prop)),] = mutation_data[n,] # attribute mutation to all cancer cell types
        colnames(mat.mut) = paste0('mut.',colnames(mutation_data))
        df_expr_0 = cbind(df_expr,mat.mut) # Add mutations as covariates
        df_expr_1 = cbind(df_expr,mat.mut) # Add mutations as covariates
        rownames(df_expr_0) = paste0(rownames(df_expr_0), '.0')
        rownames(df_expr_1) = paste0(rownames(df_expr_1), '.1')
        tabs_[[n]] = rbind(df_expr_0, df_expr_1)
    }
    cliff_input = do.call(rbind, tabs_)
    prop_col = num(cliff_input[,grepl('props',colnames(cliff_input))])
    cliff_input = cliff_input[,!grepl('props',colnames(cliff_input))]
    # Expend Drug sensitivity data to match CLIFF input
    drug_data = drug_data[rep(seq_len(nrow(drug_data)), each = 2*K), ]
    drug_data$y = drug_data$auc
    drug_data$y_bin = drug_data$auc
    for(i in 0:(2*N-1)){
        start_ = i*K+1
        end_ = (i+1)*K
        if ( i %% 2 == 1){
            drug_data$y_bin[start_:end_] <- 0
        } else {
            drug_data$y_bin[start_:end_] <- 1
        }
    }

    
    ### CLIFF model with Expectation-Maximization Algorithm ###
    message('Launch CLIFF algorithm and iterate the EM algorithm')
    drug_data$pi_hat = mean_auc ; drug_data$y_hat = 0
    drug_data$q_0 = 0 ; drug_data$q_1 = 0
    drug_data$w_0 = 0 ; drug_data$w_1 = 0
    prev_rmse = 1e5 ; lambda = 0
    drug_data$y = drug_data$auc
    it.increasing.rmse = 0
    min_rmse = 1000
    # EM algorithm with 50 iterations
    for(e in 1:em_steps){
        # E-step - compute values for q_0 and q_1 
        for(i in 0:(2*N-1)){
            start_ = i*K+1
            end_ = (i+1)*K  
            div_factor_0 = sum(prop_col[start_:end_]*(1-drug_data$pi_hat[start_:end_]))
            div_factor_1 = sum(prop_col[start_:end_]*drug_data$pi_hat[start_:end_])
            drug_data$q_0[start_:end_] <- ((prop_col[start_:end_]*(1-drug_data$pi_hat[start_:end_]))/div_factor_0)
            drug_data$q_1[start_:end_] <- ((prop_col[start_:end_]*drug_data$pi_hat[start_:end_])/div_factor_1)
            drug_data$w_0[start_:end_] <- (1-drug_data$y[start_:end_])*drug_data$q_0[start_:end_]
            drug_data$w_1[start_:end_] <- drug_data$y[start_:end_]*drug_data$q_1[start_:end_]
        }
        # M-step - define weigths and optimize beta coefficients
        y = drug_data$y_bin
        drug_data$weights = drug_data$w_0
        drug_data$weights[drug_data$y_bin == 1] <- drug_data$w_1[drug_data$y_bin == 1]
        w_ = drug_data$weights
        # We launch glmnet with binomial mode (logistic regression) with weights
        set.seed(1)
        if(regularization == 'L2'){
            fit = glmnet(as.matrix(cliff_input), y, weights = w_, family = 'binomial', 
                      intercept = F, scale=F, alpha=0) # lambda=0: no regularization
        } else { 
            fit = glmnet(as.matrix(cliff_input), y, weights = w_, family = 'binomial', 
                         intercept = F, scale=F, lambda=0) # alpha=0: L2-regularization
        }
        cliff_coefs = coef(fit)[-1,dim(coef(fit))[2]]
        # Compute predicted values for pi - cell-type drug sensitivity
        drug_data$pi_hat = sigmoid_(as.matrix(cliff_input) %*% cliff_coefs, log((1-mean_auc)/mean_auc))
        for(i in 0:(2*N-1)){
            start_ = i*K+1
            end_ = (i+1)*K  
            drug_data$y_hat[start_:end_] <- sum(prop_col[start_:end_]*drug_data$pi_hat[start_:end_])
        }
        # Record error and metrics
        drug_data.sub = drug_data[drug_data$y_bin == 1,]
        drug_data.sub$cellType = rep(colnames(climb_prop), N)
        pred_mean_pi = aggregate(drug_data.sub$pi_hat, list(drug_data.sub$cellType), mean)
        # Recover values for pi
        pi_hat_nk = data.frame(drug_data.sub$pi_hat)
        pi_hat_nk$sample = gsub('.*_','',rownames(pi_hat_nk))
        pi_hat_nk$cellType = gsub('_.*','',rownames(pi_hat_nk))
        colnames(pi_hat_nk) = c('pi_hat_nk', 'sample', 'celltype')
        PI_hat_nk = dcast(pi_hat_nk, sample ~ celltype, value.var="pi_hat_nk") 
        rownames(PI_hat_nk) = PI_hat_nk[,1] ; PI_hat_nk = PI_hat_nk[,-1]
        rownames(PI_hat_nk) = gsub('\\.0','',rownames(PI_hat_nk))
        PI_hat_nk = PI_hat_nk[rownames(mutation_data),]
        rmse_ = rmse(rowSums(as.matrix(PI_hat_nk) * as.matrix(climb_prop)), drug_data.sub$auc)
        if(rmse_ < min_rmse){
            min_rmse = rmse_
            min_coefs = cliff_coefs
            min_PI = PI_hat_nk
            it.increasing.rmse = 0
        } else {
            it.increasing.rmse = it.increasing.rmse + 1
            if( it.increasing.rmse >= 2 ){
                message(paste0('early stopping of EM algorithm at step ', e))
                break
            }
        }
    }
    return(list(min_PI, mutation_data, climb_prop))
}

