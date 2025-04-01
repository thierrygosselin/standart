#' @title clean_fq_filename
#' @description Removes the last part of the fq filename
#' @rdname clean_fq_filename
#' @export
#' @keywords internal
clean_fq_filename <- function(x) {
  x %<>%
    stringi::stri_replace_all_fixed(
      str = .,
      pattern = unique(fq_file_type(.)),
      replacement = "",
      vectorize_all = FALSE
    )
}#End clean_fq_filename


#' @title stats_standart
#' @description Generate useful stats
#' @rdname stats_standart
#' @export
#' @keywords internal

stats_standart <- function(data, x, group.by = NULL, digits = NULL) {

  if (!is.null(group.by)) {
    data <- dplyr::group_by(.data = data, .data[[group.by]])
  }

  s <- dplyr::summarise(
    .data = data,
    N = n(),
    SUM = sum(.data[[x]], na.rm = TRUE),
    MEAN = mean(.data[[x]], na.rm = TRUE),
    SE = sqrt(stats::var(.data[[x]]) / length(.data[[x]])),
    SD = stats::sd(.data[[x]], na.rm = TRUE),
    MEDIAN = stats::median(.data[[x]], na.rm = TRUE),
    Q25 = stats::quantile(.data[[x]], 0.25, na.rm = TRUE),
    Q75 = stats::quantile(.data[[x]], 0.75, na.rm = TRUE),
    IQR = stats::IQR(.data[[x]], na.rm = TRUE),
    # IQR = abs(diff(stats::quantile(.data[[x]], probs = c(0.25, 0.75), na.rm = TRUE))),
    MIN = min(.data[[x]], na.rm = TRUE),
    MAX = max(.data[[x]], na.rm = TRUE),
    OUTLIERS_LOW = Q25 - (1.5 * IQR),
    OUTLIERS_HIGH = Q75 + (1.5 * IQR),
    OUTLIERS_LOW = ifelse(OUTLIERS_LOW < MIN, MIN, OUTLIERS_LOW), # don'T use dplyr::if_else here... you don't want to preserve types
    OUTLIERS_LOW_N = length(.data[[x]][.data[[x]] < OUTLIERS_LOW]),
    OUTLIERS_HIGH = ifelse(OUTLIERS_HIGH > MAX, MAX, OUTLIERS_HIGH),
    OUTLIERS_HIGH_N = length(.data[[x]][.data[[x]] > OUTLIERS_HIGH]),
    .groups = "keep"
  )

  if (!is.null(digits)) {
    s %<>% dplyr::mutate(.data = ., dplyr::across(.cols = where(is.numeric), .fns = round, digits = digits))
  }

  return(s)
}#End stats_standart


#' @title list_sample_file
#' @description List sample file in folder
#' @rdname list_sample_file
#' @export
#' @keywords internal
list_sample_file <- function(f, full.path = FALSE, recursive = FALSE, paired.end = FALSE) {
  sample_file <- function(x, f) {
    sample.file <- list.files(
      path = f,
      pattern = x,
      full.names = full.path,
      recursive = recursive
    )

    # fq files with .rem.
    not.wanted <- list.files(
      path = f,
      pattern = ".rem.",
      full.names = full.path,
      recursive = recursive
    )

    sample.file <- purrr::keep(.x = sample.file, .p = !sample.file %in% not.wanted)

    if (!paired.end) {
      # reverse paired-end files
      not.wanted <- list.files(
        path = f,
        pattern = "\\.2\\.",
        full.names = full.path,
        recursive = recursive
      )

      sample.file <- purrr::keep(.x = sample.file, .p = !sample.file %in% not.wanted)
    }

    if (length(sample.file) > 0) {
      return(sample.file)
    } else {
      return(NULL)
    }
  }
  sample.list <- purrr::map(
    .x = c("fq.gz", "fq", "fasta", "fastq", "gzfasta", "gzfastq", "fastq.gz", "FASTQ.gz", "FASTQ.GZ"),
    .f = sample_file, f = f) %>%
    purrr::flatten_chr(.) %>%
    unique
  return(sample.list)
}#End list_sample_file



#' @title fq_file_type
#' @description Detect fq file type
#' @rdname fq_file_type
#' @export
#' @keywords internal
fq_file_type <- function(x) {
  fq.file.type <-  suppressWarnings(
    stringi::stri_match_all_regex(
      str = x,
      omit_no_match = TRUE,
      pattern = c( ".fq", ".fq.gz", ".fasta", ".gzfasta", ".gzfastq", ".fastq", ".fastq.gz", ".FASTQ.gz", ".FASTQ.GZ")
    ) %>%
      purrr::flatten_chr(.)
  ) %>%
    unique

  # if (identical(x = c(".fastq", ".fastq.gz", ".FASTQ.gz", ".FASTQ.GZ"), y = fq.file.type)) fq.file.type <- ".fastq.gz"
  # if (identical(x = c(".FASTQ.GZ"), y = fq.file.type)) fq.file.type <- ".fastq.gz"
  # if (identical(x = c(".FASTQ.gz"), y = fq.file.type)) fq.file.type <- ".fastq.gz"
  if (identical(x = c(".fastq", ".fastq.gz"), y = fq.file.type)) fq.file.type <- ".fastq.gz"
  if (identical(x = c(".fq", ".fq.gz"), y = fq.file.type)) fq.file.type <- ".fq.gz"
  if (identical(x = c(".fasta", ".gzfasta"), y = fq.file.type)) fq.file.type <- ".gzfasta"
  return(fq.file.type)
}#End fq_file_type


# standart_future --------------------------------------------------------------
#' @name standart_future
#' @title radiator parallel function
#' @description Updating radiator to use future
# @inheritParams future::plan
# @inheritParams future::availableCores
#' @inheritParams future.apply::future_apply
#' @rdname standart_future
#' @keywords internal
standart_future <- function(
    .x,
    .f,
    flat.future = c("int", "chr", "dfr", "dfc", "walk", "drop"),
    split.vec = FALSE,
    split.with = NULL,
    split.chunks = 4L,
    parallel.core = parallel::detectCores() - 1,
    forking = FALSE,
    ...
) {
  os <- Sys.info()[['sysname']]
  if (os == "Windows") forking <- FALSE

  opt.change <- getOption("width")
  options(width = 70)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(if (parallel.core > 1L && !forking) future::plan(strategy = "sequential"), add = TRUE)

  # argument for flattening the results
  flat.future <- match.arg(
    arg = flat.future,
    choices = c("int", "chr", "dfr", "dfc", "walk", "drop"),
    several.ok = FALSE
  )

  # splitting into chunks-------------------------------------------------------
  if (split.vec && is.null(split.with)) {
    # d: data, data length, data size
    # sv: split vector
    d <- .x
    df <- FALSE
    if (any(class(d) %in% c("tbl_df","tbl","data.frame"))) {
      d <- nrow(d)
      df <- TRUE
    }
    if (length(d) > 1L) d <- length(d)
    stopifnot(is.integer(d))
    sv <- as.integer(floor((split.chunks * (seq_len(d) - 1) / d) + 1))
    # sv <- as.integer(floor((parallel.core * cpu.rounds * (seq_len(d) - 1) / d) + 1))
    stopifnot(length(sv) == d)

    # split
    if (df) {
      .x$SPLIT_VEC <- sv
      .x %<>% dplyr::ungroup(.) %>% dplyr::group_split(.tbl = ., "SPLIT_VEC", .keep = FALSE)
    } else {
      .x %<>% split(x = ., f = sv)
    }
  }
  if (!is.null(split.with)) {
    # check
    if (length(split.with) != 1 || !is.character(split.with)) {
      rlang::abort(message = "Contact author: problem with parallel computation")
    }
    .data <- NULL
    stopifnot(rlang::has_name(.x, split.with))
    if (split.vec) {
      sv <- dplyr::distinct(.x, .data[[split.with]])
      d <- nrow(sv)
      sv$SPLIT_VEC <- as.integer(floor((split.chunks * (seq_len(d) - 1) / d) + 1))
      .x %<>%
        dplyr::left_join(sv, by = split.with) %>%
        dplyr::ungroup(.) %>%
        dplyr::group_split(.tbl = ., "SPLIT_VEC", .keep = FALSE)
    } else {
      .x %<>% dplyr::ungroup(.) %>% dplyr::group_split(.tbl = ., .data[[split.with]], .keep = TRUE)
    }
  }


  if (parallel.core == 1L) {
    future::plan(strategy = "sequential")
  } else {
    # parallel.core <- parallel_core_opt(parallel.core = parallel.core)
    lx <- length(.x)
    if (lx < parallel.core) {
      future::plan(strategy = "multisession", workers = lx)
    } else {
      if (!forking) future::plan(strategy = "multisession", workers = parallel.core)
    }
  }

  # Run the function in parallel and account for dots-dots-dots argument

  if (forking) {
    if (length(list(...)) == 0) {
      rad_map <- switch(
        flat.future,
        int = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, mc.cores = parallel.core) %>%
            purrr::flatten_int(.)
        },
        chr = {.x %<>%
            parallel::mclapply(X = ., FUN = .f, mc.cores = parallel.core) %>%
            purrr::flatten_chr(.)
        },
        dfr = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, mc.cores = parallel.core) %>%
            purrr::flatten_dfr(.)
        },
        dfc = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, mc.cores = parallel.core) %>%
            purrr::flatten_dfc(.)
        },
        walk = {furrr::future_walk},
        drop = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, mc.cores = parallel.core)
        }
      )
    } else {
      rad_map <- switch(
        flat.future,
        int = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, ..., mc.cores = parallel.core) %>%
            purrr::flatten_int(.)
        },
        chr = {.x %<>%
            parallel::mclapply(X = ., FUN = .f, ..., mc.cores = parallel.core) %>%
            purrr::flatten_chr(.)
        },
        dfr = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, ..., mc.cores = parallel.core) %>%
            purrr::flatten_dfr(.)
        },
        dfc = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, ..., mc.cores = parallel.core) %>%
            purrr::flatten_dfc(.)
        },
        walk = {furrr::future_walk},
        drop = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, ..., mc.cores = parallel.core)
        }
      )
    }
  } else {
    rad_map <- switch(flat.future,
                      int = {furrr::future_map_int},
                      chr = {furrr::future_map_chr},
                      dfr = {furrr::future_map_dfr},
                      dfc = {furrr::future_map_dfc},
                      walk = {furrr::future_walk},
                      drop = {furrr::future_map}
    )
    p <- NULL
    p <- progressr::progressor(along = .x)
    opts <- furrr::furrr_options(globals = FALSE, seed = TRUE)
    if (length(list(...)) == 0) {
      .x %<>% rad_map(.x = ., .f = .f, .options = opts)
    } else {
      .x %<>% rad_map(.x = ., .f = .f, ..., .options = opts)
    }
  }
  return(.x)
}#End standart_future
