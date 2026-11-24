
suppressPackageStartupMessages(library(jsonlite))

if (!file.exists("figures/data/manifest.json"))
  stop("Missing figures/data/manifest.json. Restore this file from the code archive.")

.MF      <- fromJSON("figures/data/manifest.json", simplifyVector = FALSE)
.MF_META <- .MF[["_meta"]]
.mf_used <- new.env(parent = emptyenv())

mf_entry <- function(key) {
  e <- .MF[[key]]
  if (is.null(e))
    stop("manifest key does not exist: ", key,
         "\nRestore the matching manifest from the code archive.")
  assign(key, TRUE, envir = .mf_used)
  e
}

## Does the key exist? Use before a fallback chain; mf() deliberately stops on an
## unknown key so that a typo cannot silently return NA and leave a blank in a panel.
mf_has <- function(key) !is.null(.MF[[key]])

mf_status <- function(key) {
  e <- mf_entry(key)
  if (identical(e$status, "OWED")) "OWED" else if (identical(e$status, "NA")) "NA" else "ok"
}

## The value. Numeric scalars come back as numeric; lists come back as vectors.
## "OWED" (value not yet supplied) and "NA" (not applicable) entries return NA, so
## a placeholder cannot be plotted as a number.
mf <- function(key, quiet = FALSE) {
  e <- mf_entry(key)
  st <- if (identical(e$status, "OWED")) "OWED" else if (identical(e$status, "NA")) "NA" else "ok"
  if (st != "ok") {
    if (!quiet) message(sprintf("  [%s] %s -- %s", st, key, e$reason))
    return(NA)
  }
  v <- e$value
  if (is.list(v)) v <- unlist(v)
  if (!quiet) {
    src <- if (!is.null(e$source_file)) sprintf("%s :: %s", e$source_file, e$source_key)
           else sprintf("computed <- %s", paste(unlist(e$inputs), collapse = ", "))
    message(sprintf("  [mf] %-46s %s   [%s]", key,
                    if (is.numeric(v) && length(v) == 1) sprintf("%.6g", v)
                    else sprintf("<%s[%d]>", class(v)[1], length(v)), src))
  }
  v
}

mf_num <- function(key, quiet = TRUE) {
  v <- mf(key, quiet = quiet)
  if (length(v) != 1) return(NA_real_)
  as.numeric(v)
}

## What to draw where a number would have gone.
mf_label <- function(key, fmt = "%.3g", na_text = "n/a", owed_text = "[VALUE OWED]") {
  switch(mf_status(key),
         OWED = owed_text,
         "NA" = na_text,
         sprintf(fmt, mf_num(key)))
}

mf_is_owed <- function(key) mf_status(key) == "OWED"
mf_is_na   <- function(key) mf_status(key) == "NA"

mf_keys <- function() sort(ls(.mf_used))

## Write the list of manifest keys used, beside the figure. The final check reads this
## rather than grepping the script: keys built in a loop or by paste() are invisible
## to a regex but are recorded here because the file is produced by running the plot code.
mf_write_keys <- function(path_noext) {
  k <- mf_keys()
  writeLines(k, paste0(path_noext, ".keys.txt"))
  owed <- k[vapply(k, mf_is_owed, logical(1))]
  nas  <- k[vapply(k, mf_is_na,   logical(1))]
  message(sprintf("  [keys] %d manifest keys used -> %s.keys.txt", length(k),
                  basename(path_noext)))
  if (length(owed)) {
    message("  [keys] VALUE OWED, rendered as a visible placeholder:")
    for (x in owed) message("           ", x)
  }
  if (length(nas))
    message(sprintf("  [keys] %d keys are NOT APPLICABLE and are drawn as such, never as zero",
                    length(nas)))
  invisible(k)
}
