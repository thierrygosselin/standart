#' @name reads_length_distribution
#' @title Generate the read length distribution


#' @description This function reads the fastq file of an individual, lane/chip or
#' entire folder of fastq files.
#' The read length distribution is generated and saved in the working directory.
#' This helps to decide the threshold to cut the reads to a specific length.

#' @description
#' This function reads FASTQ files corresponding to a single individual, a sequencing lane/chip,
#' or an entire folder containing multiple FASTQ files. It generates a distribution of read lengths,
#' which is saved to the working directory.
#'
#' This distribution can be used to determine an appropriate threshold for keeping or
#' trimming reads of a specific length.

#' @param fq.files (character, path) A character vector of paths to FASTQ files or a folder containing FASTQ files.
#' Default: \code{fq.files = "my-sample.fq.gz"}.

#' @param with.future (logical) Whether to use parallel processing via the future package.
#' Default: \code{with.future = FALSE}.

#' @param parallel.core (integer) Number of cores to use when with.future is TRUE.
#' Default: \code{parallel.core = parallel::detectCores() - 1}.

#' @param chunk.size (integer) Number of reads per chunk to stream from each FASTQ file.
#' Default: \code{chunk.size = 1e6}.

#' @param show.progress (logical) Whether to display a progress bar.
#' Default: \code{show.progress = TRUE}.



#' @details
#'
#' coming soon, just try it in the meantime...
#'

#' @rdname reads_length_distribution
#' @export

#' @return A tibble with columns: FILE and READ_LENGTH.


#' @examples
#' \dontrun{
#' # Example 1: Process a folder of FASTQ files with progress bar
#'
#' read_lengths <- read_length_many(
#'   fq.files = fq.files,
#'   with.future = TRUE,
#'   parallel.core = 4,
#'   chunk.size = 500000,
#'   show.progress = TRUE
#' )
#'
#' # View summary
#' dplyr::glimpse(read_lengths)
#'
#' # Optional: Plot read length distribution
#' library(ggplot2)
#' ggplot(read_lengths, aes(x = READ_LENGTH)) +
#'   geom_histogram(binwidth = 1) +
#'   facet_wrap(~ FILE, scales = "free_y") +
#'   theme_minimal()
#' }


reads_length_distribution <- function(
    fq.files,
    parallel.core = parallel::detectCores() - 1,
    with.future = FALSE,
    chunk.size = 1e6,
    show.progress = TRUE
) {
  timing <- proc.time()

  # Required package -----------------------------------------------------------
  # vroom turns out to be slower to do this kind of stuff...

  # if (!"vroom" %in% utils::installed.packages()[,"Package"]) {
  #   rlang::abort('Please install vroom for this option:\n
  #                install.packages("vroom")')
  # }
  if (with.future) {
    if (!"listenv" %in% utils::installed.packages()[,"Package"]) {
      rlang::abort('Please install listenv for this option:\n
                 install.packages("listenv")')
    }
    if (!"future" %in% utils::installed.packages()[,"Package"]) {
      rlang::abort('Please install future for this option:\n
                 install.packages("future")')
    }
    if (!"future.apply" %in% utils::installed.packages()[,"Package"]) {
      rlang::abort('Please install future.apply for this option:\n
                 install.packages("future.apply")')
    }
  }

  # ShortRead dependency
  if (!"ShortRead" %in% utils::installed.packages()[,"Package"]) {
    rlang::abort('Please install ShortRead:\n
                 BiocManager::install("ShortRead")')
  }

  # check if argument is a file or a directory...
  if (fs::is_dir(fq.files)) fq.files <- list.files(path = fq.files, full.names = TRUE)

  read_length_many <- function(
    fq.files,
    with.future = FALSE,
    parallel.core = parallel::detectCores() - 1,
    chunk.size = 1e6,
    show.progress = TRUE
    ) {
    # with future
    if (with.future) future::plan(strategy = "multisession", workers = parallel.core)

    # Progress handler
    if (show.progress) {
      progressr::handlers(global = TRUE)
      p <- progressr::progressor(along = fq.files)
    } else {
      p <- function(...) NULL
    }

    # Process a single FASTQ file
    process_file <- function(one.fq.file) {
      f <- ShortRead::FastqStreamer(one.fq.file, n = chunk.size, verbose = FALSE)
      result <- list()
      counter <- 0L

      while (length(fq <- ShortRead::yield(f))) {
        counter <- counter + 1L
        result[[counter]] <- as.integer(fq@sread@ranges@width)
      }
      close(f)

      sample.clean.name <- standart::clean_fq_filename(x = one.fq.file)
      p(sprintf("Finished %s", sample.clean.name))
      message("Finished counting:", sample.clean.name)

      # tibble
      fq <- tibble::tibble(
        FILE = sample.clean.name,
        READ_LENGTH = as.integer(unlist(purrr::list_flatten(result)))
      )
      p(sprintf("Number of reads:  %s", nrow(fq)))
      message("Number of reads: ", nrow(fq))

      # Visualisation ----
      cum_length <- function(threshold, x) x <- length(x$READ_LENGTH[x$READ_LENGTH >= threshold])
      read.breaks <- read.seq <- seq(from = (max(min(fq$READ_LENGTH), 50)), to = (min(max(fq$READ_LENGTH), 200)), by = 10L)

      names(read.seq) <- read.seq
      reads.info <- tibble::tibble(
        READS_LENGTH = read.seq,
        N = purrr::map_int(.x = read.seq, .f = cum_length, x = fq)
      )
      reads.plot <- ggplot2::ggplot(
        data = reads.info,
        ggplot2::aes(x = READS_LENGTH, y = N)) +
        ggplot2::geom_line() +
        ggplot2::geom_point(size = 2, shape = 21, fill = "white") +
        ggplot2::scale_x_continuous(name = "Read length maximum size", breaks = read.breaks) +
        ggplot2::scale_y_continuous(name = "Number of reads") +
        ggplot2::theme_bw() +
        ggplot2::theme(
          axis.title.x = ggplot2::element_text(size = 10, face = "bold"),
          axis.title.y = ggplot2::element_text(size = 10, face = "bold"),
          axis.text.x = ggplot2::element_text(size = 8)
        )
      filename.plot <- stringi::stri_join(sample.clean.name, "_reads_length_dist.png")
      ggplot2::ggsave(
        plot = reads.plot,
        filename = filename.plot,
        width = 25,
        height = 15,
        dpi = 300,
        units = "cm"
      )
      return(fq)
    }# END process_file

    # Process all files (in parallel or not)
    if (with.future) {
      res.list <- progressr::with_progress({
        future.apply::future_lapply(X = fq.files, FUN = process_file)
      })
    } else {
      res.list <- progressr::with_progress({
        lapply(X = fq.files, FUN = process_file)
      })
    }

    # Combine all results into one tibble
    res.list <- dplyr::bind_rows(res.list)
    return(res.list)

  }#End read_length_many

  fq <- read_length_many(fq.files = fq.files, with.future = with.future, parallel.core = parallel.core)

  timing <- proc.time() - timing
  message("\nComputation time: ", round(timing[[3]]), " sec")
  return(fq)
}# End reads_length_distribution
