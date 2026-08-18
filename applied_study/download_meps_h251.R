# Download and prepare the public MEPS 2023 Full-Year Consolidated file (HC-251).
#
# This script makes the MEPS arm of the applied study reproducible from a clean
# checkout. The file is released by the Agency for Healthcare Research and
# Quality and is freely redistributable; it is fetched here rather than stored
# in the repository only because it is roughly 60 MB as CSV.
#
# CONDITIONS OF USE. Two restrictions attach to the MEPS public-use files and
# survive redistribution: the data may be used for statistical reporting and
# analysis only, and they may not be linked with individually identifiable
# records from any source other than MEPS itself and the National Health
# Interview Survey from which the MEPS sample is drawn.
#
# Usage:  Rscript applied_study/download_meps_h251.R

if (!requireNamespace("haven", quietly = TRUE)) {
  stop("Package 'haven' is required to read the Stata release of HC-251.")
}

OUTDIR <- "applied_study"
URL <- "https://meps.ahrq.gov/mepsweb/data_files/pufs/h251/h251dta.zip"
CSV_PATH <- file.path(OUTDIR, "h251.csv")

if (file.exists(CSV_PATH)) {
  cat("Already present:", CSV_PATH, "\n")
  cat("Delete it first if you want to re-download.\n")
  quit(save = "no")
}

work_dir <- tempfile("meps_h251_")
dir.create(work_dir)
on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
zip_path <- file.path(work_dir, "h251dta.zip")

cat("Downloading", URL, "\n")
status <- utils::download.file(URL, zip_path, mode = "wb", quiet = FALSE)
if (status != 0L || !file.exists(zip_path)) {
  stop("Download failed. The AHRQ release page is https://meps.ahrq.gov/",
       "mepsweb/data_stats/download_data_files_detail.jsp?cboPufNumber=HC-251")
}

utils::unzip(zip_path, exdir = work_dir)
dta_files <- list.files(work_dir, pattern = "\\.dta$", full.names = TRUE,
                        ignore.case = TRUE)
if (length(dta_files) != 1L) {
  stop("Expected exactly one .dta file in the archive; found ", length(dta_files))
}

cat("Reading", basename(dta_files), "\n")
meps <- haven::read_dta(dta_files[1L])
meps <- as.data.frame(lapply(meps, function(column) {
  if (inherits(column, "haven_labelled")) haven::zap_labels(column) else column
}), stringsAsFactors = FALSE)

# The analysis addresses variables by their upper-case public-use names.
names(meps) <- toupper(names(meps))

required <- c("AGE23X", "SEX", "RACETHX", "HISPANX", "FAMINC23", "POVCAT23",
              "REGION23", "MARRY23X", "DNTINS23_M23", "DVTEXP23",
              "PERWT23F", "VARSTR", "VARPSU")
missing_columns <- setdiff(required, names(meps))
if (length(missing_columns) > 0L) {
  stop("The downloaded file is missing: ", paste(missing_columns, collapse = ", "))
}

write.csv(meps, CSV_PATH, row.names = FALSE)
cat(sprintf("Wrote %s: %d records, %d variables.\n",
            CSV_PATH, nrow(meps), ncol(meps)))
cat("Next: Rscript applied_study/run_study.R\n")
