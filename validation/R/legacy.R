## Load the original SIGA scripts (../legacy) into an environment. library()
## and require() are stubbed so that the scripts do not attach packages; the
## gatbxr operators and survival functions are made visible explicitly.
load_legacy <- function(legacy_dir = file.path(validation_root(), "..", "legacy")) {
  env <- new.env(parent = globalenv())
  env$library <- function(...) invisible(TRUE)
  env$require <- function(...) invisible(TRUE)
  for (f in c("crtrp", "ranking", "reclin", "recombin", "reins", "mutbga", "sus"))
    assign(f, getExportedValue("gatbxr", f), envir = env)
  for (f in c("coxph", "Surv", "survfit")) assign(f, getExportedValue("survival", f), envir = env)
  env$rmvnorm <- mvtnorm::rmvnorm
  for (f in c("decodsubset.R", "objectiveFunctions.R", "VariableScreen.R", "SIGA.tte.R",
              "Best.R", "simul.Cox.R"))
    sys.source(file.path(legacy_dir, f), envir = env)
  env
}

validation_root <- function() {
  root <- Sys.getenv("SIGA_VALIDATION_DIR", "")
  if (nzchar(root)) return(normalizePath(root, winslash = "/"))
  normalizePath(".", winslash = "/")
}
