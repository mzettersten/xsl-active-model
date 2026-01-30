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
    if(xd_s[i,]$trial_type=="selection") {
      refs = na.omit(c(xd_s[i,]$choice_image, xd_s[i,]$random_image))
      labs = na.omit(c(xd_s[i,]$choice_label, xd_s[i,]$random_label))
    } else {
      refs = setdiff(unlist(xd_s[i,] |> select(starts_with("image"))), NA)
      labs = setdiff(unlist(xd_s[i,] |> select(starts_with("label"))), NA)
    }
    coocs[labs,refs] = coocs[labs,refs] + 1
  }
  return(coocs)
}


# specific to crossact
get_subject_trial_ordering <- function(xd_s, keep_selection_trials = TRUE) {
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
      words[[i]] = na.omit(c(xd_train[i,]$choice_label, xd_train[i,]$random_label))
      objects[[i]] = na.omit(c(xd_train[i,]$choice_image, xd_train[i,]$random_image))
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


# take a single subject's selection trials and their association matrix, and evaluate given selection strategy
eval_selection_strategy <- function(mat, selections, strategy, temperature = 1.0) {
  # mat: current association matrix (words x objects)
  # selections: data frame with choice_image and choice_label columns
  # strategy: which selection strategy to evaluate
  # temperature: softmax temperature (higher = more random)

  nlogLik = 0
  probOfSelections = rep(NA, nrow(selections))

  for(i in 1:nrow(selections)) {
    # Get available options (all objects in matrix)
    available_objects = colnames(mat)
    n_options = length(available_objects)

    # Calculate selection probabilities based on strategy
    if(strategy == "entropy") {
      # Choose objects with highest entropy (most uncertain associations)
      selection_scores = apply(mat, 2, shannon.entropy)

      # Handle NAs and zeros
      selection_scores[is.na(selection_scores)] = 0

      if(sum(selection_scores) == 0) {
        # If no entropy anywhere, choose uniformly
        predictedSelections = rep(1/n_options, n_options)
      } else {
        # Softmax over entropy scores
        predictedSelections = exp(selection_scores / temperature)
        predictedSelections = predictedSelections / sum(predictedSelections)
      }

    } else if(strategy == "random") {
      # Baseline: uniform selection
      predictedSelections = rep(1/n_options, n_options)

    } else if(strategy == "margin") {
      # Choose objects where top 2 associations are most similar (least confident)
      predictedSelections = rep(NA, n_options)

      for(o in 1:n_options) {
        word_assocs = mat[, o]

        if(sum(word_assocs > 0) < 2) {
          # Can't compute margin with fewer than 2 non-zero associations
          predictedSelections[o] = 1  # High uncertainty
        } else {
          top2 = sort(word_assocs, decreasing = TRUE)[1:2]
          margin = top2[1] - top2[2]
          # Small margin = uncertain = high selection probability
          predictedSelections[o] = 1 / (margin + 0.01)  # Add small constant to avoid division by zero
        }
      }

      # Normalize
      predictedSelections = predictedSelections / sum(predictedSelections)

    } else if(strategy == "confirmatory") {
      # Choose objects with strongest single association (most certain)
      max_strength = apply(mat, 2, max)

      if(sum(max_strength) == 0) {
        predictedSelections = rep(1/n_options, n_options)
      } else {
        predictedSelections = exp(max_strength / temperature)
        predictedSelections = predictedSelections / sum(predictedSelections)
      }

    } else if(strategy == "novelty") {
      # Choose objects with weakest associations (least familiar)
      max_strength = apply(mat, 2, max)

      # Novelty is inverse of familiarity
      novelty_scores = 1 / (max_strength + 0.01)  # Add constant to avoid division by zero

      predictedSelections = exp(novelty_scores / temperature)
      predictedSelections = predictedSelections / sum(predictedSelections)

    } else {
      stop("Unknown strategy: ", strategy)
    }

    # Get actual choice
    actual_choice = selections[i,]$choice_image

    # Find probability of actual choice
    choice_idx = which(available_objects == actual_choice)

    if(length(choice_idx) == 0) {
      warning("Choice ", actual_choice, " not found in available objects for subject on trial ", i)
      actual_selection_prob = 1e-10  # Small probability as penalty
    } else {
      actual_selection_prob = predictedSelections[choice_idx]

      # Ensure probability is not too small (for numerical stability)
      actual_selection_prob = max(actual_selection_prob, 1e-10)
    }

    # Accumulate negative log-likelihood
    # This is CATEGORICAL choice, not binary
    nlogLik = nlogLik - log(actual_selection_prob)
    probOfSelections[i] = actual_selection_prob

    # Update matrix with the selection that was made
    # (for next trial's predictions)
    mat[selections[i,]$choice_label, selections[i,]$choice_image] =
      mat[selections[i,]$choice_label, selections[i,]$choice_image] + 1
  }

  return(list(
    nlogLik = nlogLik,
    probOfSelections = probOfSelections,
    mat = mat,
    mean_prob = mean(probOfSelections)  # Average predicted probability
  ))
}

# iterates through participants, gets their model-based association matrix,
eval_selection_trials_group <- function(ords, xd, strategy, use_model = FALSE,
                                        model_params = NULL, temperature = 1.0) {
  total_logLik = 0
  n_selections = 0

  for(s in names(ords)) {
    # Get this subject's selection trials
    selections_s = xd |> filter(subject == s, trial_type == "selection")

    if(nrow(selections_s) == 0) {
      next  # Skip subjects with no selection trials
    }

    if(use_model && !is.null(model_params)) {
      # Use model-based association matrix
      start_mat = crossact_mats[[xd |> filter(subject == s) |>
                                   pull(experiment_name) |> first()]]
      res = model(model_params, ord = ords[[s]], start_matrix = start_mat)
      mat = res$matrix
    } else {
      # Use co-occurrence matrix from training
      mat = create_cooc_matrix(ords[[s]])
    }

    # Evaluate selection strategy
    sel = eval_selection_strategy(mat, selections_s, strategy = strategy,
                                  temperature = temperature)

    total_logLik = total_logLik + sel$nlogLik
    n_selections = n_selections + nrow(selections_s)
  }

  # Return both total and per-selection average
  return(list(
    total_nlogLik = total_logLik,
    mean_nlogLik = total_logLik / n_selections,
    n_selections = n_selections
  ))
}
