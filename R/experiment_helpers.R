get_experiment_cooc_matrix <- function(dat, expt_name, iterations=0) {
  this_exp = dat |> filter(experiment_name==expt_name)
  xd_s <- subset(this_exp, subject==unique(this_exp$subject)[1])
  voc = sort(na.omit(unique(unlist(xd_s |> select(starts_with("label"))))))
  refs = sort(na.omit(unique(unlist(xd_s |> select(starts_with("image"))))))
  coocs <- matrix(0, nrow=length(voc), ncol=length(refs))
  row.names(coocs) = voc
  colnames(coocs) = refs
  if(iterations>0) {
    for(i in 1:iterations) {
      words <- na.omit(unlist(xd_s[i,] |> select(starts_with("label"))))
      refs <- unlist(xd_s[i,] |> select(starts_with("image")))
      refs <- refs[which(!is.na(refs))]
      coocs[words,refs] = coocs[words,refs] + 1
    }
  }
  return(coocs)
}


run_subject_trials <- function(xd_s) {
  coocs <- get_experiment_cooc_matrix(xd_s, xd_s[1,]$experiment_name)
  # iterate over trials, handle training and selection/test trials
  for(i in 1:nrow(xd_s)) {
    # should ignore
    refs = setdiff(unlist(xd_s[i,] |> select(starts_with("image"))), NA)
    labs = setdiff(unlist(xd_s[i,] |> select(starts_with("label"))), NA)
    coocs[labs,refs] = coocs[labs,refs] + 1
  }
  return(coocs)
}


# specific to crossact
get_subject_trial_ordering <- function(xd_s, keep_selection_trials = FALSE) {
  if(keep_selection_trials) {
    xd_train <- xd_s |> filter(trial_type!="test") # can result in a 1x1 learning trial - break any models?
  } else {
    xd_train <- xd_s |> filter(trial_type=="learning") # remove "selection" and "test" trials
  }
  voc = sort(na.omit(unique(unlist(xd_s |> select(starts_with("label"))))))
  refs = sort(na.omit(unique(unlist(xd_s |> select(starts_with("image"))))))
  w_inds = which(startsWith(names(xd_train), "label")) # columns with label info
  o_inds = which(startsWith(names(xd_train), "image")) # columns with object info
  words = list() # words per trial
  objects = list() # objs per trial
  for(i in 1:nrow(xd_train)) {
    if(xd_train[i,]$trial_type=="learning") {
      words[[i]] = as.vector(unlist(xd_train[i,w_inds][which(!is.na(xd_train[i,w_inds]))]))
      objects[[i]] = as.vector(unlist(xd_train[i,o_inds][which(!is.na(xd_train[i,o_inds]))]))
    } else if(xd_train[i,]$trial_type=="selection") {
      words[[i]] = unlist(xd_train[i,]$choice_label)
      objects[[i]] = unlist(xd_train[i,]$choice_image)
    }
  }
  return(list(words = words, objs = objects,
              voc=voc, refs=refs)) # for indexing
}



# fit an entire group of subjects with a single best-fitting set of parameters
fit_model_to_group <- function(model_name, exp_name, xd, lower, upper) {
  source(paste0("R/", model_name, ".R"))
  xd_exp <- xd |> filter(experiment_name == exp_name)
  ords <- list()
  for(s in unique(xd_exp$subject)) {
    ords[[s]] <- get_subject_trial_ordering(xd_exp |> filter(subject == s))
  }

  testdat_group <- xd |> filter(experiment_name == exp_name,
                                trial_type == "test", trial_index < 37)

  fit <- DEoptim::DEoptim(eval_model_group, lower = lower, upper = upper,
                          DEoptim::DEoptim.control(reltol = .001, NP = 200, itermax = 200),
                          ords = ords, testdat = testdat_group, logLik_only = TRUE,
                          start_matrix = crossact_mats[[exp_name]])
  return(fit$optim$bestmem)
}


eval_model_group <- function(parms, ords, testdat, start_matrix, logLik_only = TRUE) {
  total_logLik = 0

  for(s in names(ords)) {
    testdat_s <- testdat |> filter(subject == s)
    res <- model(parms, ord = ords[[s]], start_matrix = start_matrix)
    mat <- res$matrix[testdat_s$target_label, testdat_s$target_image]
    mod_perf <- diag(mat) / rowSums(mat)

    is_right <- testdat_s$is_right
    subject_logLik <- -sum(is_right * log(mod_perf) + (1 - is_right) * log(1 - mod_perf))
    total_logLik <- total_logLik + subject_logLik
  }

  return(total_logLik)
}
