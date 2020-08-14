
new_ensemble <- function(learners = NULL) {
  sl3::make_learner(sl3::Lrnr_sl,
                    learners = learners,
                    metalearner = sl3::make_learner("Lrnr_nnls", convex = TRUE),
                    keep_extra = FALSE)
}

new_sl3 <- function(data, Y, X, id = NULL) {
  sw(sl3::sl3_Task$new(
    data = data,
    covariates = X,
    outcome = Y,
    outcome_type = "binomial",
    id = id,
    drop_missing_outcome = FALSE
  ))
}

run_ensemble <- function(ensemble, task) {
  if (!is.null(ensemble)) {
    ensemble <- sl3::delayed_learner_train(ensemble, task)
    if (check_future_job()) {
      sched <- delayed::Scheduler$new(ensemble, delayed::FutureJob, nworkers = get_workers())
    } else {
      sched <- delayed::Scheduler$new(ensemble, delayed::SequentialJob)
    }
    sched$compute()
  }
}

predict_sl3 <- function(object, task) {
    return(object$predict(task))
}
