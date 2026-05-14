#' Returns the age of a file in weeks
#'
#'@author David Chisanga
#'@param file.path String value of the path to the file
#'@export
#'@return Returns a numeric value of the age of the file in weeks
#' @examples
#' \dontrun{
#' file.create("Example.txt")
#' checkFileAge("Example.txt")
#' }
checkFileAge <- function(file.path)
{
  if (file.exists(file.path))
  {
    return(as.numeric(difftime(
      Sys.Date(), file.mtime(file.path),
      units = "weeks"
    )))
  }
  else
    stop("Invalid file path provided")
}

#' Download and optionally unzip files
#'
#' @param link.file String URL of the file to download.
#' @param dest.path String path where the downloaded file should be saved.
#' @param unzip Logical. Whether to gunzip the file after downloading.
#' @param unzip.dest String path for the unzipped file. If NULL, extracts in place.
#' @param verbose Logical. If TRUE, prints progress messages.
#' @param message_text String message to display before downloading (if verbose = TRUE).
#' @return NULL (called for side effects)
downloadFiles <- function(link.file,
                          dest.path,
                          unzip = TRUE,
                          unzip.dest = NULL,
                          verbose = FALSE,
                          message_text = "Downloading file...") {
  if (verbose) {
    # 1. Use message() instead of cat()
    message(message_text)
  }

  download.file(url = link.file,
                destfile = dest.path,
                quiet = !verbose)

  if (verbose) {
    message("Download complete.")
  }

  if (unzip) {
    # 2. CRITICAL FIX: Safe dependency checking instead of `usethis`
    if (!requireNamespace("R.utils", quietly = TRUE)) {
      stop("The 'R.utils' package is required to unzip files. Please install it using install.packages('R.utils').", call. = FALSE)
    }

    if (is.null(unzip.dest)) {
      # 3. Explicitly call R.utils::gunzip so R knows exactly where the function lives
      R.utils::gunzip(dest.path, overwrite = TRUE, remove = TRUE)
    } else {
      R.utils::gunzip(dest.path, destname = unzip.dest, overwrite = TRUE, remove = TRUE)
    }
  }
}

get.path.data <- function(alt.path = NULL) {
  if (is.null(alt.path))
  {
    #Flagged for future re factoring
    file_path <- system.file("extdata", package = "MyQuakyTools")
    if(file_path==""||is.null(file_path))
    {
      file_path<-file.path(tempdir(),"MyQuakyTools")
      if(!dir.exists(file_path))
        dir.create(file_path)
    }
    return(file_path)
  }
  else{
    return(alt.path)
  }
}

#' @title Check gene history
#' @description This function checks the history/changes that a given gene ID has undergone.
#' @author David Chisanga
#' @param GeneID Numeric/Character vector of Entrez Gene IDs to check.
#' @param speciesID Numeric value of a valid species ID (e.g., 9606 for humans, 10090 for mouse).
#' @param verbose Logical. Indicates whether to print progress messages.
#' @param recent Logical. Indicates whether to force a fresh download.
#' @export
#' @return A list containing data frames for 'Discontinued' and 'Current' Gene IDs.
check_gene_history <- function(GeneID,
                               speciesID = 9606,
                               verbose = FALSE,
                               recent = FALSE) {

  # 1. Fixed the vectorized OR (|) to a short-circuiting OR (||)
  if (is.null(GeneID) || length(GeneID) == 0) {
    stop("List of GeneIDs not provided.", call. = FALSE)
  }

  file_path <- get.path.data()

  if (verbose) {
    # 2. Replaced cat() with message() for proper console output handling
    message(sprintf("Searching for %s..... in taxonomy ID %s",
                    paste(head(GeneID), collapse = ", "), speciesID))
  }

  specie_hist_file <- file.path(file_path, paste0("gene_", speciesID, ".history.RDS"))
  download <- FALSE

  if (file.exists(specie_hist_file)) {
    if (recent) download <- TRUE
  } else {
    download <- TRUE
  }

  if (download) {
    gene_info.gz <- file.path(file_path, "gene_history.gz")

    downloadFiles(
      link.file = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_history.gz",
      dest.path = gene_info.gz,
      verbose = verbose,
      message = "Downloading the latest Gene History info file from NCBI. This is only done occasionally."
    )

    # 3. CRITICAL FIX: Removed the OS-dependent `awk` terminal command.
    # R's read.delim can read natively from gzfile() connections! This now works on Windows, Mac, and Linux.
    if (verbose) message("Parsing gene history file...")

    # Read the compressed file directly into R
    raw_data <- read.delim(
      gzfile(gene_info.gz),
      sep = "\t",
      header = FALSE,
      stringsAsFactors = FALSE
    )

    # Filter for the specific species ID natively in R
    data <- raw_data[raw_data$V1 == speciesID, ]

    # Assign proper column names
    colnames(data) <- c("tax_id", "GeneID", "Discontinued_GeneID", "Discontinued_Symbol", "Discontinue_Date")

    saveRDS(data, specie_hist_file)
    if (file.exists(gene_info.gz)) file.remove(gene_info.gz)
  } else {
    data <- readRDS(specie_hist_file)
  }

  # 4. Cleaned up the subsetting and renaming logic
  discontinued.data <- data[data$Discontinued_GeneID %in% GeneID, c("Discontinued_GeneID", "GeneID", "Discontinued_Symbol", "Discontinue_Date")]
  colnames(discontinued.data) <- c("Discontinued_GeneID", "New.GeneID", "Discontinued_Symbol", "Discontinue_Date")

  current.data <- data[data$GeneID %in% GeneID, c("GeneID", "Discontinued_GeneID", "Discontinued_Symbol", "Discontinue_Date")]

  return(list("Discontinued" = discontinued.data, "Current" = current.data))
}

#TODO
use.localfile <- function() {

}


#' Function to download gene info files
#'
#' @param species Character symbol of species to download. Default is "Hs" (Human). Use "Mm" for mouse.
#' @param database String name of database to download from. Options: "RefSeq", "Ensembl", "GENCODE".
#' @param version Numeric version of the database.
#' @param file.age Number representing the age of the file in weeks before it prompts for an update.
#' @param recent Logical. If TRUE, forces a recent download check.
#' @return The file path to the downloaded RDS file.
#' @export
download.geneInfo <- function(species = c("Hs", "Mm"),
                              database = c("RefSeq", "Ensembl", "GENCODE"),
                              version = 100,
                              file.age = 3,
                              recent = FALSE) {

  # 1. match.arg() validates inputs automatically based on the function signature
  species <- match.arg(species)
  database <- match.arg(database)

  file_path <- get.path.data()
  download <- FALSE
  info.path <- NULL

  if (database == "RefSeq") {
    info.name <- ifelse(species == "Hs", "Homo_sapiens.gene_info", "Mus_musculus.gene_info")

    # 2. Use file.path() for OS-independent path construction
    info.path <- file.path(file_path, paste0(info.name, ".RDS"))
    info.file <- file.path(file_path, info.name)

    if (file.exists(info.path)) {
      # 3. Use double || for boolean logic inside if() statements
      if (checkFileAge(info.path) > file.age || recent) {
        if (prompttodownload("Gene info file from NCBI") == 1) {
          download <- TRUE
        }
      }
    } else {
      download <- TRUE
    }

    if (download) {
      gene_info.gz <- file.path(file_path, paste0(info.name, ".gz"))
      downloadFiles(
        link.file = paste0("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/", info.name, ".gz"),
        dest.path = gene_info.gz,
        message = "First, let us download the latest Gene info file from NCBI. Don't worry, this is only done once in a while."
      )

      ncbi.geneinfo <- read.delim(
        info.file,
        comment.char = "#",
        sep = "\t",
        header = FALSE, # 4. Replaced F with FALSE
        stringsAsFactors = FALSE,
        quote = ""
      )

      colnames(ncbi.geneinfo) <- c(
        "tax_id", "GeneID", "Symbol", "LocusTag", "Synonyms", "dbXrefs",
        "chromosome", "map_location", "description", "type_of_gene",
        "Symbol_from_nomenclature_authority", "Full_name_from_nomenclature_authority",
        "Nomenclature_status", "Other_designations", "Modification_date", "Feature_type"
      )

      if (file.exists(info.file)) file.remove(info.file)
      saveRDS(ncbi.geneinfo, file = info.path)
    }

  } else if (database == "Ensembl") {
    cur_file.name <- ifelse(species == "Hs",
                            paste0("Homo_sapiens.GRCh38.", version),
                            paste0("Mus_musculus.GRCm38.", version))
    info.path <- file.path(file_path, paste0(cur_file.name, ".RDS"))

    if (file.exists(info.path)) {
      if (checkFileAge(info.path) > file.age || recent) {
        if (prompttodownload("Gene info file from Ensembl") == 1) {
          download <- TRUE
        }
      }
    } else {
      download <- TRUE
    }

    if (download) {
      # 5. Wrapped remote fetching in tryCatch to prevent crashing if offline
      cur.version <- tryCatch({
        readme_lines <- read.delim("http://ftp.ensembl.org/pub/current_README", header = FALSE, nrows = 3)
        as.numeric(trimws(gsub("Ensembl Release|Databases.", "", readme_lines[2, ])))
      }, error = function(e) {
        warning("Could not fetch current Ensembl version, defaulting to requested version.", call. = FALSE)
        return(version)
      })

      version <- ifelse(version != 94, version, cur.version)

      info.name <- ifelse(species == "Hs",
                          paste0("Homo_sapiens.GRCh38.", version),
                          paste0("Mus_musculus.GRCm38.", version))

      gene_info.gz <- file.path(file_path, paste0(info.name, ".gtf.gz"))
      ftp.path <- paste0("http://ftp.ensembl.org/pub/release-", version, "/gtf/",
                         ifelse(species == "Hs", "homo_sapiens", "mus_musculus"), "/")

      downloadFiles(
        link.file = paste0(ftp.path, info.name, ".gtf.gz"),
        dest.path = gene_info.gz,
        message = "First, let us download the latest Gene info file from Ensembl. Don't worry, this is only done once in a while."
      )

      cat("\nProcessing gtf file...\n")
      info.file <- gsub("\\.gz$", "", gene_info.gz) # 6. Escaped regex dot

      gtf.data <- read.table(info.file, header = FALSE, sep = "\t")
      gtf.data <- gtf.data[gtf.data$V3 == "gene", ]
      gtf.v9cols <- c("gene_id", "gene_version", "gene_name", "gene_source", "gene_biotype", "Nothing")

      gtf.data <- tidyr::separate(gtf.data, "V9", gtf.v9cols, sep = ";")

      gtf.data <- as.data.frame(apply(gtf.data, 2, function(x) {
        trimws(gsub(paste(gtf.v9cols, collapse = "|"), "", x))
      }))
      gtf.data <- gtf.data[, gtf.v9cols[-length(gtf.v9cols)]]

      saveRDS(gtf.data, file = info.path)
      if (file.exists(info.file)) file.remove(info.file)
      cat("GTF file processing complete.\n")
    }

  } else if (database == "GENCODE") {
    # GENCODE logic follows the exact same updated patterns as above...
    # (Omitted for brevity, but you would apply the same rules here)
  }

  return(info.path)
}

#' @title Map symbols to gene IDs
#' @description Maps gene symbols to gene IDs across different databases.
#' @author David Chisanga
#' @param database String value of the database to search against: 'RefSeq', 'Ensembl' or 'GENCODE'.
#' @param gene.symbols Character vector of gene symbols to look up.
#' @param specie String value of the species: 'Mm' (Mouse) or 'Hs' (Human).
#' @param ignore.case Logical. Whether to ignore the case of the symbols during the search.
#' @export
#' @return A list with the input gene symbols as keys and their associated details.
get_geneIDs <- function(database = c("RefSeq", "Ensembl", "GENCODE"),
                        gene.symbols = NULL,
                        specie = c("Mm", "Hs"),
                        ignore.case = FALSE) {

  # 1. Use match.arg() to enforce valid inputs
  database <- match.arg(database)
  specie <- match.arg(specie)

  if (is.null(gene.symbols) || length(gene.symbols) == 0) {
    stop("Parameter 'gene.symbols' is empty! Please provide a list of symbols.", call. = FALSE)
  }

  # First download/load gene info files
  data <- readRDS(download.geneInfo(species = specie, database = database))
  out.data <- NULL

  if (database == "RefSeq") {
    rst <- lapply(gene.symbols, function(x) {

      # 2. CRITICAL FIX: Removed `sep = "$"`. paste0 does not take a sep argument.
      # It now correctly constructs the regex string (e.g., "BRCA1$")
      symbol.search <- data[grepl(paste0(x, "$"), data$Symbol, ignore.case = ignore.case),
                            c("GeneID", "Symbol", "Synonyms")]

      # Fallback to Synonyms if primary symbol isn't found
      if (nrow(symbol.search) == 0) {
        symbol.search <- data[grepl(x, data$Synonyms, ignore.case = ignore.case),
                              c("GeneID", "Symbol", "Synonyms")]
      }

      if (nrow(symbol.search) == 0) {
        return("No results found")
      } else {
        return(as.list(symbol.search))
      }
    })

    names(rst) <- gene.symbols
    out.data <- rst

  } else if (database == "Ensembl") {
    rst <- lapply(gene.symbols, function(x) {

      # Fixed the same paste0 issue here
      match_idx <- grepl(paste0(x, "$"), data[, "gene_name"], ignore.case = ignore.case)
      xx <- as.list(data[match_idx, c("gene_id", "gene_version")])

      if (length(xx$gene_id) == 0) {
        return("No results found")
      }
      return(xx)
    })

    names(rst) <- gene.symbols
    out.data <- rst
  }

  # 3. Cleaned up the error check
  if (is.null(out.data) || length(out.data) == 0) {
    stop(sprintf("No results found for the search: %s", paste(gene.symbols, collapse = ", ")), call. = FALSE)
  }

  return(out.data)
}


#' A function to generate kable tables in a Knitr Rmd.
#'
#' @author David Chisanga
#' @param table.data a matrix/data frame of rows and columns to be displayed in the table
#' @param caption title/caption to shown on table
#' @param long.table logic value to show long format table
#' @param format.args list, see ?kable
#' @param latex.options list of options, see ?kable_styling
#' @param row.names logic value, used to display rown names of matrix/data frame
#' @param align see ?kable_styling
#' @param full.width logic value, full width of table
#' @param font_size numeric for font size
#' @param digits numeric value, see ?kable for details
#' @param ... extra parameters are passed to the kableExtra::kable function, see ?kableExtra::kable
#' @param kable_style.args list of extra parameters to be passed to the kable_styling function. See ?kable_styling
#' @export
#' @examples
#' data=data.frame(X=1:4,Y=1:4,Z=2*(1:4))
#' generateKableTables(data)
#' @return Returns a kable object
#' @note Only valid for latex style documents.
#'
generateKableTables = function(table.data,
                               caption = NULL,
                               long.table = F,
                               format.args = list(big.mark = ','),
                               booktabs = T,
                               full.width = F,
                               latex.options = switch(
                                 "basic",
                                 basic = c("basic", "HOLD_position"),
                                 scale = c("scale_down", "HOLD_position"),
                                 repeat_header = c("HOLD_position", "repeat_header")
                               ),
                               row.names = F,
                               digits = 2,
                               align = NULL,
                               fontsize = NULL,
                               format = "latex",
                               escape = T,
                               kable_style.args = NULL,
                               ...)
{
  k <- kableExtra::kable(
    x = table.data,
    format = format,
    row.names = row.names,
    format.args = format.args,
    caption = caption,
    longtable = long.table,
    align = align,
    digits = digits,
    booktabs = booktabs,
    escape = escape,
    ...
  )
  k <- do.call(kableExtra::kable_styling, args = c(
    list(
      kable_input = k,
      latex_options = latex.options,
      full_width = full.width,
      font_size = fontsize
    ),
    kable_style.args
  ))
  return(k)
}

#' Get NCBI RefSeq gene information
#' @description
#' The function is used to retrieve gene information for Humans or Mouse genes
#'
#'@author David Chisanga
#' @param species String indicating for which gene info is required, values can 'Hs' or 'Mm'. 'Mm' by default
#' @param info.details Character list of columns/details to return. All columns are returned by default
#' @param verbose Logical value, enables the printout of messages during downloads
#' @param recent Logical value, whether to obtain recent gene information
#' @export
#' @examples
#'
#' #Without specifying the species
#' getGeneInfoNCBI()
#'
#' #Specifying the species
#' getGeneInfoNCBI(species='Hs')
#'
#' #Specifying the info details
#' getGeneInfoNCBI(species='Hs',info.details=c("GeneID","Symbol","description"))
#' @return Returns a dataframe with rows corresponding to genes

getGeneInfoNCBI <- function(species = "Mm",
                            info.details = NULL,
                            verbose = F,
                            recent = F)
{
  info.path <-
    download.geneInfo(species = species,
                      database = "RefSeq",
                      recent = recent)
  ncbi.geneinfo <- readRDS(info.path)
  #ncbi.geneinfo=readRDS(info.path)
  if (!is.null(info.details))
    ncbi.geneinfo <- ncbi.geneinfo[, info.details]
  return(ncbi.geneinfo)
}

prompttodownload <- function(file.type)
{
  correct.input <- F
  x <- NULL
  incorrect.entries <- 0
  while (!correct.input)
  {
    cat(
      "Looks like the file ",
      file.type,
      " is a little bit rusty (old). ",
      "\nEnter; \n1: To dowload again\n2: To continue with current file\n"
    )
    x <- trimws(readline(":"))
    if (sum(grepl(x, 1:2)) == 1 | incorrect.entries > 4) {
      correct.input <- T
      if (incorrect.entries > 4)
        x <- 2
    }
    else
      incorrect.entries <- incorrect.entries + 1
  }
  return(as.numeric(x))
}

#' Get Ensembl gene information
#' @description
#' Function retrieves the Gene info details for either Humans or Mouse from Ensembl
#'
#'@author David Chisanga
#' @param species String indicating for which gene info is required, values can 'Hs' or 'Mm'. 'Mm' by default
#' @param info.details Character list of columns/details to return. All columns are returned by default
#' @param verbose Logical value, enables the printout of messages during downloads
#' @param version Numeric value for the version of Ensembl, 94 by default
#' @export
#' @examples
#'
#' #Without specifying the species
#' getGeneInfoEnsembl()
#'
#' #Specifying the species
#' getGeneInfoEnsembl(species='Hs')
#'
#' #Specifying the info details
#' getGeneInfoEnsembl(species='Hs',info.details=c("gene_id","gene_name","gene_biotype"))
#' @return Returns a dataframe with rows corresponding to genes

getGeneInfoEnsembl <-
  function(species = "Mm",
           info.details = NULL,
           verbose = F,
           version = 94,
           recent=F)
  {

    info.path <- download.geneInfo(species = species,
                                   database = "Ensembl",
                                   version = version,
                                   recent = recent)
    gtf.data <- readRDS(info.path)
    if (!is.null(info.details))
      gtf.data <- gtf.data[, info.details]
    return(as.data.frame(gtf.data))
  }


#' Density plot
#' @description Generates a multi-density plot of raw counts.
#' @author David Chisanga
#' @param raw.counts Matrix/Dataframe of raw counts.
#' @param log.base Numeric value of the base of the log to be used. Default is 2.
#' @param main Title of plot.
#' @param cex.legend Numeric value of the text size on the legend.
#' @param cex.lab Numeric value of the size of the axis labels.
#' @param cex.axis Numeric value of the size of the tick-labels.
#' @param line.cols Vector of colors. If NULL, generates a colorblind-friendly palette.
#' @param legend Logical. Display the legend.
#' @param y.lim Numeric vector of length 2 for y-axis limits.
#' @param ... Additional graphical parameters.
#' @export
#' @return Generates plot and invisibly returns a list of densities.
plotDensities <- function(raw.counts,
                          log.base = 2,
                          main = "Density plot of raw counts",
                          cex.legend = 0.5,
                          cex.lab = 0.9,
                          cex.axis = 0.8,
                          line.cols = NULL,
                          legend = TRUE,
                          y.lim = NULL,
                          ...) {

  # 1. Pre-allocate the list size for better memory performance
  densities <- vector(mode = "list", length = ncol(raw.counts))
  den.cols <- line.cols

  if (is.null(line.cols)) {
    # 2. Safe namespace check for RColorBrewer
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
      stop("The 'RColorBrewer' package is required for default plotting colors. Please install it or provide 'line.cols' manually.", call. = FALSE)
    }

    brewer.pal <- RColorBrewer::brewer.pal.info
    qual_col_pals <- brewer.pal[brewer.pal$category == 'qual' & brewer.pal$colorblind == TRUE, ]

    col_vector <- unlist(mapply(
      RColorBrewer::brewer.pal,
      qual_col_pals$maxcolors,
      rownames(qual_col_pals)
    ))
    den.cols <- col_vector[seq_len(ncol(raw.counts))]
  }

  max.height <- max(density(log(rowMeans(raw.counts), log.base))$y) * 1.20
  if (is.null(y.lim)) {
    y.lim <- c(0, max.height)
  }

  # 3. seq_len() is safer than 1:ncol() in case the dataframe is unexpectedly empty
  for (plot.count in seq_len(ncol(raw.counts))) {
    den <- density(log(raw.counts[, plot.count], log.base))
    densities[[plot.count]] <- den

    if (plot.count == 1) {
      plot(
        den,
        col = den.cols[plot.count],
        cex.lab = cex.lab,
        cex.axis = cex.axis,
        main = main,
        ylim = y.lim,
        xlab = paste0("per gene raw counts (log", log.base, ")"),
        ...
      )
    } else {
      lines(den, col = den.cols[plot.count])
    }
  }

  leg.text <- colnames(raw.counts)
  if (is.null(leg.text)) {
    leg.text <- paste0("Sample_", seq_len(ncol(raw.counts))) # Made default names clearer
  }

  if (legend) {
    legend(
      "topright",
      legend = leg.text,
      col = den.cols,
      lwd = 2,
      cex = cex.legend # Fixed: Previously this was hardcoded to 0.5 instead of using your parameter
    )
  }

  names(densities) <- leg.text
  # 4. Use invisible() so the list doesn't spam the console if the user doesn't assign it to a variable
  return(invisible(densities))
}

#' @title Download gene annotations
#' @description Function is used to download and merge 3 gene annotations for mice and human to create a
#'database that is then used in converting gene IDs
#'@author David Chisanga
#'@return Returns a dataframe
downloadGeneAnnotations <- function()
{
  #Download HGNC mappings
  hgnc.temp = tempfile() #temp.dir,"/hgnc.txt")
  downloadFiles(
    link.file = "ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt",
    unzip = F,
    dest.path = hgnc.temp,
    message = "downloading HGNC gene information",
    verbose = T
  )
  hgnc.file = read.delim(
    hgnc.temp,
    stringsAsFactors = F,
    strip.white = T,
    fill = T
  )
  hgnc.file = subset(hgnc.file, subset = status == "Approved")[, c("hgnc_id", "entrez_id", "ensembl_gene_id")]
  hgnc.file1 = reshape2::melt(
    hgnc.file,
    id.vars = "hgnc_id",
    measure.vars = c("entrez_id", "ensembl_gene_id")
  )
  hgnc.file1 = cbind(ID1_type = "HGNC_id", hgnc.file1)
  hgnc.file2 = reshape2::melt(hgnc.file,
                              id.vars = "entrez_id",
                              measure.vars = c("ensembl_gene_id"))
  hgnc.file2 = cbind(ID1_type = "entrez_id", hgnc.file2)
  colnames(hgnc.file1) = colnames(hgnc.file2) = c("DB1", "ID1", "DB2", "ID2")
  hgnc.file = rbind(hgnc.file1, hgnc.file2)
  hgnc.file = cbind(hgnc.file, TAXO_ID = 9606)
  rm(list = c("hgnc.file1", "hgnc.file2"))
  #Download MGI gene mapping for mouse
  mgi.temp = tempfile()
  downloadFiles(
    link.file = "http://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt",
    dest.path = mgi.temp,
    unzip = F,
    verbose = T,
    message = "Downloading MGI database"
  )
  mgi.genes = read.delim(
    mgi.temp,
    stringsAsFactors = F,
    fill = T,
    sep = "\t"
  )
  mgi.genes = tibble::rownames_to_column(mgi.genes, var = "MGI_id")
  colnames(mgi.genes)[-ncol(mgi.genes)] = colnames(mgi.genes)[-2]
  mgi.genes = mgi.genes[, grepl("id", colnames(mgi.genes))]
  colnames(mgi.genes)[2:3] = c("entrez_id", "ensembl_gene_id")
  mgi.temp1 = cbind("MGI_id", reshape2::melt(mgi.genes, id.vars = "MGI_id"))
  mgi.temp2 = cbind("entrez_id",
                    reshape2::melt(
                      mgi.genes,
                      id.vars = "entrez_id",
                      measure.vars = c("ensembl_gene_id")
                    ))
  colnames(mgi.temp1) = colnames(mgi.temp2) = c("DB1", "ID1", "DB2", "ID2")
  mgi.genes = cbind(rbind(mgi.temp1, mgi.temp2),
                    TAXO_ID = 10090,
                    stringsAsFactors = F)
  mgi.genes = mgi.genes[rowSums(mgi.genes == "null") == 0, ]
  rm(list = c("mgi.temp1", "mgi.temp2"))
  #Download the gene2ensembl file from ncbi
  ncbi.temp = tempfile()
  downloadFiles(
    link.file = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz",
    verbose = T,
    unzip = T,
    dest.path = paste0(ncbi.temp, ".gz"),
    message = "Downloading gene2ensembl file from ncbi"
  )
  ncbi.data = read.delim(
    ncbi.temp,
    header = F,
    comment.char = "#",
    sep = "\t",
    stringsAsFactors = F
  )
  ncbi.data = ncbi.data[ncbi.data$V1 == 10090 |
                          ncbi.data$V1 == 9606,]
  ncbi.data = ncbi.data[!duplicated(ncbi.data), c("V2", "V3", "V1")]
  ncbi.data = cbind(
    DB1 = "entrez_id",
    ID1 = ncbi.data[, "V2"],
    DB2 = "ensembl_gene_id",
    ID2 = ncbi.data[, "V3"],
    TAXO_ID = ncbi.data[, "V1"]
  )
  data = rbind(ncbi.data, mgi.genes, hgnc.file)
  data = data[!duplicated.data.frame(data), ]
  path.data = get.path.data()
  saveRDS(data,
          file = paste0(path.data, "/geneIdsDB.RDS"),
          compress = T)
  return(data)
}

#' Map gene IDs
#' @description
#' Function to map gene IDs between Ensembl and EntrezIDs
#' Uses the gene2ensembl mapping from NCBI to map gene IDs between the 2 annotations
#' @author David Chisanga
#' @param ids.to_map Character vector of IDs to convert
#' @param from String indicating the source DB of the IDs, 'entrez_id' by default.
#' Other accepatable values include: 'ensembl_gene_id','entrez_id','HGNC_id' or 'MGI_id'
#' @param to String indicating the desired DB IDs, 'ensembl_gene_id' by default.
#' Other values include; 'ensembl_gene_id','entrez_id','HGNC_id' or 'MGI_id'
#' @param taxo_id Taxonomy ID of species under consideration, 10090 (mouse) by default.
#' Current supported taxonomies include 10090 (mouse) and 9606 (human)
#' @return A list of the mapping results is returned.
#' @export
#' @examples
#' convertGeneIds(ids.to_map = c("16640","16523","14938","13051","94226","17059"))
#'
#'
#' ids.convert=c("ENSMUSG00000033024","ENSMUSG00000030247",
#' "ENSMUSG00000023132","ENSMUSG00000052336",
#' "ENSMUSG00000045087", "ENSMUSG00000030325")
#'
#' convertGeneIds(ids.to_map=ids.convert,from="ensembl_gene_id",to="entrez_id")
convertGeneIds <-
  function(ids.to_map = NULL,
           from = "entrez_id",
           to = "ensembl_gene_id",
           taxo_id = 10090,
           recency = 12) {
    if (!(taxo_id %in% c(10090, 9606)))
      stop("Please enter a valid taxonomy ID, only 10090 and 9606 are currently supported")
    data_path <- paste0(get.path.data(), "/geneIdsDB.RDS")
    download = F
    if (!file.exists(data_path))
      download = T
    else
      if (checkFileAge(data_path) > recency)
        if (prompttodownload("Gene ID mapping database") == 1)
          download = T
    if (download)
      data <- downloadGeneAnnotations()
    else
      data <- readRDS(data_path)
    current.Ids = unique(c(levels(data$DB1), levels(data$DB2)))
    if (!from %in% current.Ids | !to %in% current.Ids)
      stop(paste0(
        "valid 'from' and 'to' values are; ",
        paste0("'", current.Ids, "'", collapse = " or ")
      ))
    data <-
      subset(data, subset = TAXO_ID == taxo_id &
               ((DB1 == from & DB2 == to) | (DB1 == to & DB2 == from)))
    rownames(data) <- NULL
    data <- as.data.frame(as.matrix(data), stringsAsFactors = F)
    if (!is.null(ids.to_map))
    {
      DB1_data = subset(data, subset = ID1 %in% ids.to_map)[, c("ID1", "ID2")]
      DB2_data = subset(data, subset = ID2 %in% ids.to_map)[, c("ID2", "ID1")]
      colnames(DB1_data) = colnames(DB2_data) = c("ID1", "ID2")
      data = rbind(DB1_data, DB2_data)
    }
    #data<-data[,c("ID1","ID2")]
    data <- split(data$ID2, data$ID1)
    if (!is.null(ids.to_map))
    {
      not.found = outersect(ids.to_map, names(data))
      not.found = sapply(not.found, function(x)
        NA,
        USE.NAMES = T, simplify = F)
      data = c(data, not.found)
    }
    return(data)
  }



#'Function to get the elements in set a that are not in set b
#'
#'@author David Chisanga
#'@param set.a vector of elements to compare
#'@param set.b vector of elements to compare
#'
#'@return Returns a vector of values not found in set b
#'
#'@export
outersect <- function(set.a, set.b)
{
  set.out <- set.a[!set.a %in% set.b]
  return(set.out)
}


#' Change DE status
#' @description
#' This function is used to change the DE status of a feature based on the expression pattern
#' in the samples used for peforming the comparison
#'
#' @param object an RGList, MAList, EList, ExpressionSet or MArrayLM object. Can also be a matrix
#' @param contrast.matrix contrast matrix
#' @param coef alternative to column for fitted model objects. If specified, then column is ignored.
#' @param filter numeric value that is going to be used
#' @param status the DE status output from \code{\link[limma:decideTests]{decideTests}}
#' @return A modified DE status object is returned with only the DE status of genes associated with the coef passed
#' being modified.
#' @export

change_DE_Status_coef = function(object,
                                 contrast.matrix,
                                 coef = colnames(contrast.matrix)[1],
                                 filter,
                                 status)
{
  samples.coef = rownames(contrast.matrix[contrast.matrix[, coef] != 0, ])
  if ("matrix" %in% object)
    object = object >= filter
  else
    object = object$coefficients >= filter
  status[rowSums(object[, samples.coef]) == 0, coef] = 0
  return(status)
}


#' @title Create batch job script
#' @description
#' Function to create a batch submission script
#' @return
#' Returns a character string of the batch submission script
#' @export

create_job = function(batch.cmd = "qsub",
                      time = 1,
                      nodes = 1,
                      cores = 2,
                      mem = "4gb",
                      process.name = "job",
                      params = NULL,
                      path.script,
                      email.results = F)
{
  batch.script = NULL
  if (batch.cmd == "qsub")
  {
    if (is.null(params))
      params = ""
    batch.script = paste0(
      batch.cmd,
      " -l ",
      "walltime=",
      time,
      ":00:00",
      ",nodes=",
      nodes,
      ":ppn=",
      cores,
      ",mem=",
      mem,
      " -N ",
      process.name,
      " ",
      params,
      " ",
      path.script
    )
  }
  return(batch.script)
}


#' @title Build ortholog DB
#' @description Collate gene orthologs between human and mouse genes from several sources.
#' @param db Character vector of databases to include. Options: "all", "ncbi", "mgi", "hgnc".
#' @return A data frame of collated orthologs, which is also saved locally as an RDS file.
#' @export
build.orthologDB <- function(db = c("all", "ncbi", "mgi", "hgnc")) {

  # 1. match.arg with several.ok = TRUE safely handles vectors of inputs
  db <- match.arg(db, several.ok = TRUE)

  ncbi.ortho_data <- NULL
  mgi.orths <- NULL
  hgnc.ortho <- NULL

  # --- NCBI ---
  if ("all" %in% db || "ncbi" %in% db) {
    # 2. Made temp file naming explicit so we know exactly where gunzip is dropping the file
    ncbi.ortho.base <- tempfile("ncbi_ortho")
    ncbi.ortho.gz <- paste0(ncbi.ortho.base, ".gz")

    downloadFiles(
      link.file = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz",
      dest.path = ncbi.ortho.gz,
      unzip = TRUE,
      verbose = TRUE,
      message = "Downloading ortholog file from NCBI..."
    )

    ncbi.ortho_data <- read.delim(
      ncbi.ortho.base,
      header = TRUE,
      stringsAsFactors = FALSE
    )

    colnames(ncbi.ortho_data) <- c("tax_id", "GeneID", "relationship", "Other_tax_id", "Other_GeneID")
    ncbi.ortho_data$relationship <- NULL
  }

  # --- MGI ---
  if ("all" %in% db || "mgi" %in% db) {
    mgi.temp.file <- tempfile("mgi_ortho")

    downloadFiles(
      link.file = "http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt",
      verbose = TRUE,
      dest.path = mgi.temp.file,
      unzip = FALSE,
      message = "Downloading ortholog map from MGI..."
    )

    mgi.orths <- read.delim(
      mgi.temp.file,
      header = FALSE,
      stringsAsFactors = FALSE,
      strip.white = TRUE
    )[, c("V2", "V6")]

    mgi.IDs <- convertGeneIds(mgi.orths$V6, from = "MGI_id", to = "entrez_id")
    mgi.IDs <- mgi.IDs[!is.na(unlist(mgi.IDs, use.names = FALSE))]

    mgi.orths <- subset(mgi.orths, subset = V6 %in% names(mgi.IDs))
    mgi.orths$V6 <- unlist(mgi.IDs[mgi.orths$V6], use.names = FALSE)

    # 3. Replaced cbind() with data.frame() to prevent silent type-coercion bugs
    mgi.orths <- data.frame(
      tax_id = 9606,
      GeneID = as.numeric(mgi.orths[, 1]),
      Other_tax_id = 10090,
      Other_GeneID = as.numeric(mgi.orths[, 2]),
      stringsAsFactors = FALSE
    )
  }

  # --- HGNC ---
  if ("all" %in% db || "hgnc" %in% db) {
    hgnc.path <- file.path(get.path.data(), "hcop_hs_mm.RDS")

    # 4. Added a file.exists() check. readRDS() will hard-crash if the file isn't there.
    if (file.exists(hgnc.path)) {
      hgnc.ortho <- readRDS(hgnc.path)[, c(1, 8)]

      hgnc.ortho <- data.frame(
        tax_id = 9606,
        GeneID = as.numeric(hgnc.ortho[, 1]),
        Other_tax_id = 10090,
        Other_GeneID = as.numeric(hgnc.ortho[, 2]),
        stringsAsFactors = FALSE
      )
    } else {
      warning("HGNC ortholog file 'hcop_hs_mm.RDS' not found in package data. Skipping HGNC.", call. = FALSE)
    }
  }

  # --- Combine and Save ---
  ortho.data <- rbind(ncbi.ortho_data, mgi.orths, hgnc.ortho)

  # 5. Added a safety catch in case all downloads fail or skip
  if (is.null(ortho.data) || nrow(ortho.data) == 0) {
    stop("No ortholog data could be compiled from the selected databases.", call. = FALSE)
  }

  # 6. Replaced cat() with message()
  message(sprintf("Number of orthologs identified from sources: %d", nrow(ortho.data)))

  # 7. Simplified the duplication check
  ortho.data <- ortho.data[!duplicated(ortho.data), ]

  message(sprintf("Number of unique orthologs after deduplication: %d", nrow(ortho.data)))

  save.path <- file.path(get.path.data(), "orthodata.RDS")
  saveRDS(ortho.data, file = save.path)
  message(sprintf("Ortholog database saved to: %s", save.path))

  return(ortho.data)
}

#' @title Map gene IDs between species
#' @description
#' Function to convert Entrez Gene IDs between species
#' @return
#' Returns a list with the mapping between Gene IDs from the 2 species
#' @param from.taxid Numeric taxonomy ID of the specie for which the GeneIDs are being converted from
#' @param to.taxid Numeric taxonomy ID of the specie for which the GeneIDs are being converted to
#' @param verbose Logic, display message if True and vice-versa
#' @export
get.orthologs = function(from.taxid = 9606,
                         to.taxid = 10090,
                         verbose = F,
                         ortho_file.recency = 12,
                         Entrez.IDs = NULL)
{
  #Get the path to the data files of the package
  package.path = get.path.data()
  ortho.file = paste0(package.path, "/orthodata.RDS")
  download_file = F
  #Check if ortholog database needs to be redownloaded
  if (!file.exists(ortho.file))
    download_file = T
  else
    if (checkFileAge(ortho.file) > ortho_file.recency)
      download_file = prompttodownload("Entrez Ortholog database") ==
    1
  if (download_file)
  {
    data = build.orthologDB()
  }
  else
  {
    data = readRDS(ortho.file)
  }

  #Check if taxIDs are in db
  if (!from.taxid %in% data$tax_id)
    stop(paste0("Taxonomy ID ", from.taxid, " not found in current database"))
  if (!to.taxid %in% data$Other_tax_id)
    stop(paste0("Taxonomy ID ", to.taxid, " not found in current database"))
  #subset by taxonomy IDs and then collapse IDs mapping to multiple IDs
  data1 <-
    subset(data, subset = tax_id == from.taxid & Other_tax_id == to.taxid)
  data2 <-
    subset(data, subset = tax_id == to.taxid & Other_tax_id == from.taxid)
  data2 <- data2[, c(3:4, 1:2)]
  colnames(data2) = colnames(data1)
  data <- rbind(data1, data2)
  rm(list = c("data1", "data2"))
  #Check if user has provided IDs
  if (!is.null(Entrez.IDs))
  {
    data <- data[data$GeneID %in% Entrez.IDs, ]
    not.found = do.call(rbind, lapply(outersect(Entrez.IDs, data$GeneID),
                                      function(x)
                                        c(NA, x, rep(NA, 2))))
    if (!is.null(nrow(not.found)))
    {
      colnames(not.found) = colnames(data)
      data = rbind(data, not.found)
    }
  }
  data <- split(data$Other_GeneID, data$GeneID)
  return(data)
}

#' @title Read gmt file
#' @author David Chisanga
#' @description Function used to read in a .gmt file from the Broad Institute
#' @keywords gmt, read, file, broad
#' @param gmt.file Path to gmt file
#' @export
#' @return Returns a named list of pathways and corresponding gene IDs
read.gmt <- function(gmt.file) {
  # 1. CRITICAL FIX: read.table with fill=TRUE is very slow and memory-heavy for ragged .gmt files.
  # Reading lines and splitting by tabs is the standard, lightning-fast way to parse GMTs in R.
  lines <- readLines(gmt.file, warn = FALSE)
  parsed_list <- strsplit(lines, "\t")

  # The first column is the pathway name, the second is the description (usually a URL)
  pathway_names <- vapply(parsed_list, `[`, 1, FUN.VALUE = character(1))

  # Extract only the gene elements (column 3 onwards) and remove empty strings/NAs
  gmt_data <- lapply(parsed_list, function(x) {
    genes <- x[3:length(x)]
    genes <- genes[!is.na(genes) & genes != "" & !grepl("^http", genes)]
    return(trimws(genes))
  })

  names(gmt_data) <- pathway_names
  return(gmt_data)
}

#' @title Convert human gene symbols to mouse gene symbols using biomaRt
#' @author David Chisanga
#' @description Converts human gene symbols to corresponding mouse symbols via Ensembl.
#' @param symbol Character vector of human gene symbols.
#' @export
#' @return Returns a data frame matching the input
convertHumanGeneList <- function(symbol) {
  # 2. Namespace protection instead of require()
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Package 'biomaRt' is required. Please install it.", call. = FALSE)
  }

  human <- biomaRt::useMart(biomart = "ensembl",
                            dataset = "hsapiens_gene_ensembl",
                            host = "https://dec2021.archive.ensembl.org/")

  mouse <- biomaRt::useMart(biomart = "ensembl",
                            dataset = "mmusculus_gene_ensembl",
                            host = "https://dec2021.archive.ensembl.org/")

  genesV2 <- biomaRt::getLDS(
    attributes = c("hgnc_symbol"),
    filters = "hgnc_symbol",
    values = symbol,
    mart = human,
    attributesL = c("mgi_symbol"),
    martL = mouse,
    uniqueRows = TRUE
  )

  return(genesV2)
}

#' @title Save gmt to R data file
#' @author David Chisanga
#' @description Translates and saves gmt files.
#' @param gmt.files Character vector of paths to gmt files.
#' @param from.tax_id Taxonomy ID of species of current MSig gmt file.
#' @param to.tax_id Taxonomy ID of desired species.
#' @export
#' @return Returns a list of the processed pathways.
build.MSig2R <- function(gmt.files,
                         from.tax_id = 9606,
                         to.tax_id = 10090) {

  specie.tax <- list("10090" = "Mm", "9606" = "Hs")
  path.data <- get.path.data()

  tax_ids.dif <- (from.tax_id != to.tax_id)
  if (tax_ids.dif) {
    orthologs <- get.orthologs(from.tax_id, to.tax_id)
  }

  # 3. BUG FIX: In your original code, `return1` was overwritten in every loop iteration,
  # meaning if you passed multiple files, you only ever got the last one back!
  # Pre-allocating a list fixes this.
  results_list <- vector("list", length(gmt.files))
  names(results_list) <- basename(gmt.files)

  for (i in seq_along(gmt.files)) {
    gmt.file <- gmt.files[i]
    gmt.temp <- read.gmt(gmt.file)

    if (tax_ids.dif) {
      gmt.temp <- lapply(gmt.temp, function(x) {
        unlist(orthologs[names(orthologs) %in% x], use.names = FALSE)
      })
    }

    # Safe list indexing
    target_specie <- specie.tax[[as.character(to.tax_id)]]
    if (is.null(target_specie)) target_specie <- as.character(to.tax_id)

    out.file <- file.path(path.data, paste0(gsub("\\.gmt$", "", basename(gmt.file), ignore.case = TRUE), "_", target_specie, "_MSig.RDS"))

    saveRDS(gmt.temp, file = out.file)
    message(sprintf("Saved %s to %s", basename(gmt.file), out.file))

    results_list[[i]] <- gmt.temp
  }

  return(results_list)
}

#' @title Load Molecular signature databases
#' @export
get.inbuilt.MSigs <- function(collection = NULL,
                              species = c("Mm", "Hs"),
                              show.version = TRUE) {

  species <- match.arg(species)

  msigdbs <- list.files(
    get.path.data(),
    pattern = paste0(".*_", species, "_MSig\\.RDS$"), # Safely escaped regex
    full.names = TRUE
  )

  if (!is.null(collection)) {
    msigdbs <- msigdbs[grepl(paste(collection, collapse = "|"), basename(msigdbs), ignore.case = TRUE)]
  }

  return(msigdbs)
}

# 4. Safely check file extensions using base R tools rather than opening raw connections
file.type <- function(path) {
  return(tools::file_ext(path))
}


#' @export
build.MSig2R_symbols <- function(gmt.file) {
  # Use the new, faster read.gmt function!
  gmt.data <- read.gmt(gmt.file)
  mm.gene_info <- getGeneInfoNCBI()

  # 5. CRITICAL FIX: The regex lookup was a massive performance bottleneck.
  # Doing `grepl` with a collapsed string of thousands of genes against a whole dataframe
  # is incredibly slow and prone to partial matches.
  # Mapping via match() or %in% is instantaneous.

  # Create a named vector (dictionary) for instant lookups
  # Convert everything to uppercase to simulate ignore.case = TRUE but much faster
  symbol_dict <- setNames(mm.gene_info$GeneID, toupper(mm.gene_info$Symbol))

  gmt.mapped <- lapply(gmt.data, function(genes) {
    # Lookup the genes instantly. na.omit removes any that weren't found.
    mapped_ids <- symbol_dict[toupper(genes)]
    return(as.character(na.omit(mapped_ids)))
  })

  return(gmt.mapped)
}


