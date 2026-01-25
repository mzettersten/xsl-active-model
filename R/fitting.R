order_dir = "orders/"
model_dir = "models/"
data_dir = "data/"


#load(here::here("data/combined_data.RData"))

# imported from fitting_functions.R

#' Run specified model on provided condition(s)
#'
#'
#' @return Returns model association matrix after training, conditional probability of selecting correct referent given each word, and sum of squared error (SSE)
#' @export
run_model <- function(conds, model_name, parameters, SSE_only=F, print_perf=F) {
  source(here::here(paste0(model_dir,model_name,".R"))) # sourcing not allowed in packages
  if(!is.null(conds$train)) {
    mod = list(perf = model(parameters, ord=conds$train)$perf)
    SSE = sum( (mod$perf - conds$HumanItemAcc)^2 )
  } else {
    mod = list()
    SSE = 0
    unweighted_SSE = 0
    totSs = 0
    for(i in 1:length(conds)) {
      #print(conds[[i]]$Condition)
      #print(parameters)
      mp = model(parameters, conds[[i]]$train)
      if(!is.null(conds[[i]]$test)) {
        mperf = mafc_test(mp$matrix, conds[[i]]$test)
        mod[[names(conds)[i]]] = mperf
      } else {
        mperf = mp$perf
        mod[[names(conds)[i]]] = mperf
      }
      SSE = SSE + conds[[i]]$Nsubj * sum( (mperf - conds[[i]]$HumanItemAcc)^2 )
      unweighted_SSE = unweighted_SSE + sum( (mperf - conds[[i]]$HumanItemAcc)^2 )
      totSs = totSs + conds[[i]]$Nsubj
    }
    SSE = SSE / totSs
  }
  if(print_perf) {
    print(mod)
    print(paste0("SSE: ",SSE))
    #print(paste0("unweighted SSE: ", unweighted_SSE)) # baseline: 31.84
  }
  mod$SSE = SSE

  if(SSE_only) return(SSE)
  return(mod)
}

fit_model <- function(model_name, conds, lower, upper) {
  fit = DEoptim::DEoptim(run_model, lower=lower, upper=upper,
                         DEoptim::DEoptim.control(reltol=.001, NP=100, itermax=100),
                model_name=model_name, conds=conds, SSE_only=T)
  return(fit)
}

stochastic_dummy <- function(n, parameters, ord) {
  return(model(parameters, ord)$perf)
}

stochastic_matrix_dummy <- function(n, parameters, ord) {
  return(model(parameters, ord)$matrix)
}

stochastic_traj_dummy <- function(n, parameters, ord) {
  return(model(parameters, ord)$traj)
}


run_stochastic_model <- function(conds, model_name, parameters, SSE_only=F, print_perf=F, get_resp_matrix=F, Nsim=500) {
  source(here::here(paste0(model_dir,"stochastic/",model_name,".R")))
  # fitting a single condition
  if(!is.null(conds$train)) {
    if(get_resp_matrix) {
      mp = lapply(1:Nsim, stochastic_matrix_dummy, parameters=parameters,
                  ord=conds$train)
      mod = Reduce('+', mp)
    } else {
      mp = sapply(1:Nsim, stochastic_dummy, parameters=parameters, ord=conds$train)
      mod = get_perf(mp)
      SSE = sum( (mod - conds$HumanItemAcc)^2 )
    }
  } else {
    mod = list()
    SSE = 0
    totSs = 0
    for(i in 1:length(conds)) {
      # just want response matrix, not AFC performance / SSE
      if(get_resp_matrix) {
        mp = lapply(1:Nsim, stochastic_matrix_dummy, parameters=parameters,
                    ord=conds[[i]]$train)
        mod[[names(conds)[i]]] = Reduce('+', mp)
      } else {
        #print(conds[[i]]$Condition)
        # 4AFC conditions use different testing
        if(!is.null(conds[[i]]$test)) {
          mp = lapply(1:Nsim, stochastic_matrix_dummy, parameters=parameters,
                      ord=conds[[i]]$train)
          mperf = Reduce('+', mp)
          mperf = mafc_test(mperf, conds[[i]]$test)
        } else {
          mp = sapply(1:Nsim, stochastic_dummy, parameters=parameters, ord=conds[[i]]$train)
          mperf = get_perf(mp)
        }
        mod[[names(conds)[i]]] = mperf
        SSE = SSE + conds[[i]]$Nsubj * sum( (mperf - conds[[i]]$HumanItemAcc)^2 )
        totSs = totSs + conds[[i]]$Nsubj
        SSE = SSE / totSs
      }

    }
    if(print_perf) {
      print(mod)
      print(paste0("SSE: ",SSE))
    }
    mod$SSE = SSE

    if(SSE_only) return(SSE)
    return(mod)
  }
}



fit_stochastic_model <- function(model_name, conds, lower, upper) {
  fit = DEoptim::DEoptim(run_stochastic_model, lower=lower, upper=upper,
                         DEoptim::DEoptim.control(reltol=.001, NP=100, itermax=30),
                model_name=model_name, conds=conds, SSE_only=T)
  return(fit)
}


# for group fits (all conditions per model)
get_model_dataframe <- function(fits, conds, cvdat=c()) {
  mdf = tidyr::tibble()
  for(model_name in names(fits)) {
    if(length(cvdat)!=0) {
      mdat = cvdat
    } else {
      if(is.element(model_name, stochastic_models)) {
        pars = fits[[model_name]]$optim$bestmem
        mdat = run_stochastic_model(conds, model_name, pars)
      } else {
        pars = fits[[model_name]]$optim$bestmem
        mdat = run_model(conds, model_name, pars)
      }
    }
    for(c in names(mdat)) {
      if(c!="SSE") {
        Nitems = length(mdat[[c]])
        HumanPerf = conds[[c]]$HumanItemAcc
        if(length(HumanPerf)==0) HumanPerf = rep(NA, Nitems)
        tmp = tidyr::tibble(Model=rep(model_name, Nitems),
                     condnum=rep(c, Nitems),
                     Condition=rep(conds[[c]]$Condition, Nitems),
                     ModelPerf=as.vector(mdat[[c]]),
                     HumanPerf=HumanPerf,
                     Nsubj=rep(conds[[c]]$Nsubj, Nitems))
        mdf = rbind(mdf, tmp)
      }
    }
  }
  return(mdf)
}


get_model_dataframe_cond_fits <- function(fits, conds) {
  mdf = tidyr::tibble()
  SSE = 0
  for(model_name in names(fits)) {
    print(model_name)
    for(c in names(fits[[model_name]])) {
      print(c)
      pars = fits[[model_name]][[c]]$optim$bestmem
      SSE = fits[[model_name]][[c]]$optim$bestval

      if(is.element(model_name, stochastic_models)) {
        mdat = run_stochastic_model(conds, model_name, pars)
      } else {
        mdat = run_model(conds[[c]], model_name, pars)
      }
      Nitems = length(mdat$perf)
      HumanPerf = conds[[c]]$HumanItemAcc
      if(length(HumanPerf)==0) HumanPerf = rep(NA, Nitems)
      tmp = tidyr::tibble(Model=rep(model_name, Nitems),
                   condnum=rep(c, Nitems),
                   Condition=rep(conds[[c]]$Condition, Nitems),
                   ModelPerf=as.vector(mdat$perf),
                   HumanPerf=HumanPerf,
                   Nsubj=rep(conds[[c]]$Nsubj, Nitems))
      mdf = rbind(mdf, tmp)
    }
  }
  return(mdf)
}



# use for fitting given model to each condition
fit_by_cond <- function(mname, conds, lower, upper) {
  mod_fits = list()
  for(cname in names(conds)) {
    print(cname)
    mod_fits[[cname]] = fit_model(mname, conds[[cname]], lower, upper)
  }
  return(mod_fits)
}


fit_stochastic_by_cond <- function(mname, conds, lower, upper) {
  mod_fits = list()
  for(cname in names(conds)) {
    print(cname)
    mod_fits[[cname]] = fit_stochastic_model(mname, conds[[cname]], lower, upper)
  }
  return(mod_fits)
}



# given a vector of indices and a list, remove all indexed items
get_train_test_split <- function(test_inds, conds) {
  dat = list(train = conds, test = list())
  for(i in test_inds) {
    dat$train[[i]] = NULL
    dat$test[[names(conds)[i]]] = conds[[i]]
  }
  return(dat)
}


# run fit_model / fit_stochastic_model for each of 5 subsets of combined_data conditions
# save the parameters...and whole dataframe?
cross_validated_group_fits <- function(model_name, combined_data, lower, upper) {
  #model_name = "kachergis"
  dat = list()
  test = list()
  testdf = tidyr::ibble()
  set.seed(123)
  folds <- caret::createFolds(names(combined_data), k = 5, list = T)
  cv_group_fits = list()
  for(i in 1:length(folds)) {
    conds = get_train_test_split(folds[[i]], combined_data)
    if(is.element(model_name, stochastic_models)) {
      opt = fit_stochastic_model(model_name, conds$train, lower, upper)
      test[[i]] = run_stochastic_model(conds$test, model_name, opt$optim$bestmem)
    } else {
      opt = fit_model(model_name, conds$train, lower, upper)
      test[[i]] = run_model(conds$test, model_name, opt$optim$bestmem)
    }
    dat[["pars"]] = rbind(dat[["pars"]], opt$optim$bestmem)
    dat[["train_acc"]] = c(dat[["train_acc"]], opt$optim$bestval)

    # add the data frame of test data? much more convenient..
    tmp = list()
    tmp[[model_name]] = opt
    testdf = rbind(testdf, get_model_dataframe(tmp, conds$test, cvdat=test[[i]]))
    print(paste("fold",i,"train SSE:",round(opt$optim$bestval,3),"test SSE:",round(test[[i]]$SSE,3)))
  }
  # add the folds ??
  dat[["test"]] = test
  dat[["testdf"]] = testdf # get_cv_test_df(test)
  return(dat)
}



#' Evaluate m-alternative forced choice test
#'
#' This function evaluates a given set of test trials using the provided model memory matrix
#' (word x referent). Each test trial is assumed to present one word and a set of referents
#' of size less than the width of the model memory matrix.
#'
#' @return A vector with the probability of choosing the correct object, given each word.
#' @export
mafc_test <- function(mperf, test) {
  perf = rep(0, length(test$trials))
  for(i in 1:length(test$trials)) {
    w = test$trials[[i]]$word
    denom = sum(mperf[w, test$trials[[i]]$objs])
    perf[i] = mperf[w,w] / denom
  }
  return(perf)
}




