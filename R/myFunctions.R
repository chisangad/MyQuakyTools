#' Returns the age of a file in weeks
#'
#'@author David Chisanga
#'@param file.path String value of the path to the file
#'@export
#'@return Returns a numeric value of the age of the file in weeks
#'@examples
#'file.create("Example.txt")
#'checkFileAge("Example.txt")
#'
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

downloadFiles <-
  function(link.file,
           dest.path,
           unzip = T,
           unzip.dest = NULL,
           verbose = F,
           message)
  {
    if (verbose)
      cat("\n", message, "\n")
    download.file(url = link.file,
                  destfile = dest.path,
                  quiet = !verbose)
    if (verbose)
      cat("\n", "Download complete", "\n")
    if (unzip)
    {
      if (is.null(unzip.dest))
        R.utils::gunzip(dest.path, overwrite = T, remove = T)
      else
        R.utils::gunzip(dest.path,
                        unzip.dest,
                        overwrite = T,
                        remove = T)
    }
  }

get.path.data <- function(alt.path = NULL) {
  if (is.null(alt.path))
  {
    file_path <- system.file("data", package = "MyQuakyTools")
    return(file_path)
  }
  else{
    return(alt.path)
  }
}

#' @title A check gene history function
#' @description This function checks the history/changes that a given gene ID has undergone
#' @author David Chisanga
#' @param GeneID list of gene Ids to check for that are likely to have been discontinued/replaced, has to be EntrezIDs
#' @param speciesID Numeric value of any valid species ID such as 9606 for humans and 10090 for mouse. 9606 by default
#' @param verbose Logical value. Indicates whether to print to the standard screen output messages
#' @keywords Genes, History
#' @export
#' @examples
#' check_gene_history(c(243187,319599,319785),10090)
#'
#' @return Returns a vector list of Discontinued and Current GeneIds from the list entered

check_gene_history <- function(GeneID,
                               speciesID = 9606,
                               verbose = F,
                               recent = F)
{
  if (is.null(GeneID) | length(GeneID) == 0)
  {
    stop("List of GeneIDs not provided")
  }
  file_path <- get.path.data()
  if (verbose)
    cat("\n",
        "Searching for ",
        paste(head(GeneID), sep = ","),
        "..... in ",
        speciesID,
        "\n")
  specie_hist_file <-
    paste0(file_path, "/gene_", speciesID, ".history.RDS")
  download = F
  if (file.exists(specie_hist_file))
  {
    if (recent)
      download = T
  } else
    download = T
  if (download)
  {
    gene_info_file <- paste0(file_path, "/gene_history")
    gene_info.gz <- paste0(file_path, "/gene_history.gz")
    downloadFiles(
      link.file = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_history.gz",
      dest.path = gene_info.gz,
      verbose = verbose,
      message = "First, let us download the latest Gene History info file\
      from NCBI. Don't worry, this is only done once in a while"
    )
    if (.Platform$OS.type != "Windows")
    {
      file.hist <- paste0(file_path, "/gene_", speciesID, ".history")
      awk_command <-
        paste0(
          'cat ',
          gene_info_file,
          ' |awk \'{split($0,b);if (b[1]==',
          speciesID,
          ') print $0}\'>',
          file.hist
        )
      system(command = awk_command)
      data <-
        read.delim(
          file.hist,
          sep = "\t",
          header = F,
          stringsAsFactors = F,
          col.names = c(
            "tax_id",
            "GeneID",
            "Discontinued_GeneID",
            "Discontinued_Symbol",
            "Discontinue_Date"
          )
        )
      saveRDS(data, specie_hist_file)
      file.remove(file.hist)
      file.remove(gene_info_file)
    }
    else{
      stop("This function is currently on available on non-windows platforms")
    }

  }
  data = readRDS(specie_hist_file)

  #search_ids = paste(paste0("^", GeneID, "$"), collapse = "|")
  discontinued.data = data[data$Discontinued_GeneID %in% GeneID,]
  discontinued.data = discontinued.data[, c("Discontinued_GeneID",
                                            "GeneID",
                                            "Discontinued_Symbol",
                                            "Discontinue_Date")]
  colnames(discontinued.data)[-1] = c("New.GeneID", "Discontinued_Symbol", "Discontinue_Date")
  current.data = data[data$GeneID %in% GeneID,]
  current.data = current.data[, c("GeneID",
                                  "Discontinued_GeneID",
                                  "Discontinued_Symbol",
                                  "Discontinue_Date")]
  return(list(Discontinued = discontinued.data, "Current" = current.data))
}

use.localfile <- function() {

}
#'Function to download gene info files
#'
#'@param species Character symbol of species to download, default is Hs. Use Mm for mouse and
#'@param database String name of database to download from. RefSeq by default. Use Ensembl
#'@param version Numeric of the version  of the database to download from
#'@param file.age Number representing the age of the file before it can be updated
download.geneInfo <- function(species = "Hs",
                              database = "RefSeq",
                              version = 100,
                              file.age = 3,
                              recent = F)
{
  file_path = get.path.data()
  download = F
  info.path = NULL
  if (database == "RefSeq")
  {
    info.name = ifelse(species == "Hs",
                       "Homo_sapiens.gene_info",
                       "Mus_musculus.gene_info")
    info.path = paste0(file_path, "/", info.name, ".RDS")
    info.file = gsub(".RDS", "", info.path)
    if (file.exists(info.path))
    {
      if (checkFileAge(info.path) > file.age | recent)
      {
        if (prompttodownload("Gene info file from NCBI") == 1)
          download = T
      }
    }
    else
      download = T
    if (download)
    {
      gene_info.gz = paste0(file_path, "/", info.name, ".gz")
      downloadFiles(
        link.file = paste0(
          "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/",
          info.name,
          ".gz"
        ),
        dest.path = gene_info.gz,
        message = "First, let us download the latest Gene info file from NCBI.,\
          Don't worry, this is only done once in a while"
      )
      ncbi.geneinfo = read.delim(
        info.file,
        comment.char = "#",
        sep = "\t",
        header = F,
        stringsAsFactors = F,
        quote = ""
      )
      colnames(ncbi.geneinfo) = c(
        "tax_id",
        "GeneID",
        "Symbol",
        "LocusTag",
        "Synonyms",
        "dbXrefs",
        "chromosome",
        "map_location",
        "description",
        "type_of_gene",
        "Symbol_from_nomenclature_authority",
        "Full_name_from_nomenclature_authority",
        "Nomenclature_status",
        "Other_designations",
        "Modification_date",
        "Feature_type"
      )
      file.remove(info.file)
      saveRDS(ncbi.geneinfo, file = info.path)
    }
  }
  else if (database == "Ensembl")
  {
    cur_file.name = ifelse(
      species == "Hs",
      paste0("Homo_sapiens.GRCh38.", version),
      paste0("Mus_musculus.GRCm38.", version)
    )
    info.path = paste0(file_path, "/", cur_file.name, ".RDS")
    if (file.exists(info.path))
    {
      if (checkFileAge(info.path) > file.age|recent)
      {
        if (prompttodownload("Gene info file from Ensembl") == 1)
          download = T
      }
    }
    else
      download = T
    if (download)
    {

      #Check current version of Ensembl
      cur.version <-
        as.numeric(trimws(gsub(
          "Ensembl Release|Databases.",
          "",
          read.delim(
            "http://ftp.ensembl.org/pub/current_README",
            header = F,
            nrows = 3
          )[2, ]
        )))
      version<-ifelse(version!=94,version,cur.version)
      info.name<-ifelse(
        species == "Hs",
        paste0("Homo_sapiens.GRCh38.", version),
        paste0("Mus_musculus.GRCm38.", version)
      )
      gene_info.gz = paste0(file_path, "/", info.name, ".gtf.gz")
      ftp.path = paste0(
        "http://ftp.ensembl.org/pub/release-",
        version,
        "/gtf/",
        ifelse(species == "Hs", "homo_sapiens", "mus_musculus"),
        "/")
      downloadFiles(
        link.file = paste0(ftp.path, info.name, ".gtf.gz"),
        dest.path = gene_info.gz,
        message = "First,let us download the latest Gene info file from Ensembl.\
          Don't worry, this is only done once in a while"
      )
      cat("\n", "Processing gtf file", "\n")
      info.file <- gsub(".gz", "", gene_info.gz)
      gtf.data = read.table(info.file, header = F, sep = "\t")
      gtf.data = gtf.data[gtf.data$V3 == "gene",]
      gtf.v9cols = c(
        "gene_id",
        "gene_version",
        "gene_name",
        "gene_source",
        "gene_biotype",
        "Nothing"
      )
      gtf.data = tidyr::separate(gtf.data, "V9", gtf.v9cols, sep = ";")
      gtf.data = apply(gtf.data, 2, function(x) {
        return(trimws(gsub(
          paste(gtf.v9cols, collapse = "|"), "", x
        )))
      })
      gtf.data = gtf.data[, gtf.v9cols[-length(gtf.v9cols)]]
      saveRDS(gtf.data, file = info.path)
      file.remove(paste0(info.file, ".gtf"))
      cat("\n", "gtf file processing complete", "\n")
    }
  }
  else if (database == "GENCODE")
  {
    info.name = paste0("gencode.v", version, ".", "annotation")
    info.path = paste0(file_path, "/", info.name, "_", species, ".RDS")
    if (file.exists(info.path))
    {
      if (checkFileAge(info.path) > file.age)
      {
        if (prompttodownload("Gene info file from GENCODE") == 1)
          download = T
      }
    }
    else
      download = T
    if (download)
    {
      release.version = paste0(ifelse(species == "Hs", "", "M"), version)
      info.name = ifelse(
        species == "Hs",
        paste0("GENCODE_Humans.", release.version),
        paste0("GENCODE_Mouse.", version)
      )
      gene_info.gz = paste0(file_path, "/", info.name, ".gtf.gz")
      ftp.path = paste0(
        "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_",
        ifelse(species == "Hs", "human", "mouse"),
        "/release_",
        release.version,
        "/gencode.v",
        release.version,
        ".annotation.gtf.gz"
      )
      info.file = paste0(file_path, "", info.name, ".gtf")
      downloadFiles(
        link.file = ftp.path,
        dest.path = gene_info.gz,
        message = "First,let us download the latest Gene info file from GENCODE\
          Don't worry, this is only done once in a while"
      )
      cat("\n", "Processing gtf file", "\n")
      gtf.v9cols = c("gene_id", "gene_type", "gene_name")
      gtf.data = read.table(info.file, header = F, sep = "\t")
      gtf.data = gtf.data[gtf.data$V3 == "gene",]
      gtf.data = tidyr::separate(gtf.data, "V9", gtf.v9cols, sep = ";")
      gtf.data = apply(gtf.data, 2, function(x) {
        return(trimws(gsub(
          paste(gtf.v9cols, collapse = "|"), "", x
        )))
      })
      gtf.data = gtf.data[, c(gtf.v9cols, "V1", "V4", "V5", "V7")]
      colnames(gtf.data) = c(gtf.v9cols, "Chr", "Start", "End", "Strand")
      gtf.data = tidyr::separate(as.data.frame(gtf.data),
                                 "gene_id",
                                 c("gene_id", "gene.version"))
      saveRDS(gtf.data, file = info.path)
      file.remove(info.file)
    }
  }
  else
  {
    stop("No valid database provided to download GeneInfo file from")
  }
  return(info.path)
}

#' @title Map symbols to gene IDs
#' @description
#' The function maps gene symbols to gene IDs
#' @author David Chisanga
#' @param database a string value of the database to search against, accepted values include; 'RefSeq', 'Ensembl' and 'GENCODE','RefSeq' by default
#' @param gene.symbols list of gene symbols to look up, NULL by default
#' @param specie a string value of the species to look up, accepted values are 'Mm' for Mouse and 'Hs' for Homo sapiens, 'Mm' by default
#' @export
#' @param ignore.case logic value indicates whether or not to ignore the case of the symbols during the search, FALSE by default
#' @keywords Genes, Gene Symbol, EntrezID, EnsemblID
#' @examples
#' rst=get.geneIDs(gene.symbols = c("mcad","Kras"),specie = "Hs",ignore.case = T)
#' @return Returns a list with the input gene symbols as keys and their associated details including, GeneID, Symbol and Synonyms.

get_geneIDs <-
  function(database = "RefSeq",
           gene.symbols = NULL,
           specie = "Mm",
           ignore.case = F)
  {
    if (is.null(gene.symbols))
      stop("parameter 'symbols' is empty! Make sure to provide a list of symbols")
    #First download gene info files
    data = readRDS(download.geneInfo(species = specie, database = database))
    out.data = NULL
    if (database == "RefSeq")
    {
      rst = lapply(gene.symbols, function(x) {
        symbol.search = data[grepl(paste0(x, sep = "$"), data$Symbol,
                                   ignore.case = ignore.case),
                             c("GeneID", "Symbol", "Synonyms")]
        if (nrow(symbol.search) == 0)
          symbol.search = data[grepl(x, data$Synonyms,
                                     ignore.case = ignore.case),
                               c("GeneID", "Symbol", "Synonyms")]
        if (nrow(symbol.search) == 0)
          symbol.search = "No results found"
        else
          symbol.search = as.list(symbol.search)
        return(symbol.search)
      })
      names(rst) = gene.symbols
      #check that the geneID is not matched to another gene symbol
      out.data = rst
    }
    else if (database == "Ensembl")
    {
      rst <- lapply(gene.symbols, function(x)
      {
        xx = as.list(data[grepl(paste0(x, "$"), data[, "gene_name"], ignore.case = ignore.case),
                          c("gene_id", "gene_version")])
        if (length(xx) == 0)
          xx = "No results found"
        return(xx)
      })
      names(rst) = gene.symbols
      out.data <- rst
    }
    if (is.null(out.data) | length(out.data) == 0)
      stop(paste0(
        "No results found for the search ",
        paste(gene.symbols, collapse = ",")
      ))
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
#' @description
#' Function to generate a multi-density plot
#'
#' @author David Chisanga
#' @param raw.counts Matrix/Dataframe of raw counts
#' @param log.base Numeric value of the base of the log to be used. 2 by default
#' @param main Title of plot
#' @param cex.main Numeric value of the size of the title
#' @param cex.lab Numeric value of the size of the axis labels
#' @param cex.axis Numeric value of the size of the tick-labels
#' @param legend Logic value to display the legend
#' @param cex.legend Numeric value of the text size on the legend
#'
#' @export
#'
#' @examples
#' raw.counts<-matrix(rnbinom(100,size=1,mu=10),25,4)
#' plotDensities(raw.counts)
#' plotDensities(raw.counts,10)
#' @return Generates plot of density values for each column (sample) and returns a list of densities

plotDensities <-
  function(raw.counts,
           log.base = 2,
           main = "Density plot of raw counts",
           cex.legend = 0.5,
           cex.lab = 0.9,
           cex.axis = 0.8,
           line.cols = NULL,
           legend = T,
           y.lim = NULL,
           ...)
  {
    den.cols <- line.cols
    densities <- vector(mode = "list")
    if (is.null(line.cols))
    {
      brewer.pal <- RColorBrewer::brewer.pal.info
      qual_col_pals <-
        brewer.pal[brewer.pal$category == 'qual' &
                     brewer.pal$colorblind == T, ]
      col_vector <- unlist(mapply(
        RColorBrewer::brewer.pal,
        qual_col_pals$maxcolors,
        rownames(qual_col_pals)
      ))
      den.cols <- col_vector[1:ncol(raw.counts)]
    }
    max.height <-
      max(density(log(rowMeans(raw.counts), log.base))$y) * 1.20
    if (is.null(y.lim))
      y.lim <- c(0, max.height)
    for (plot.count in 1:ncol(raw.counts))
    {
      den <- density(log(raw.counts[, plot.count], log.base))
      densities[[plot.count]] <- den
      if (plot.count == 1)
      {
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
      }
      else{
        lines(den, col = den.cols[plot.count])
      }
    }
    leg.text <- colnames(raw.counts)
    if (is.null(leg.text))
      leg.text <- paste0("V", 1:ncol(raw.counts))
    if (legend)
    {
      legend(
        "topright",
        legend = leg.text,
        col = den.cols,
        lwd = 2,
        cex = 0.5
      )
    }
    names(densities) <- leg.text
    return(densities)
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


#'@title Build ortholog DB
#'@description This function is used to collate gene orthologs between human and mouse genes
#' from several sources.
#' @param db Character vector of databases; currently supported options; "all" or a mix of ncbi,mgi,hgnc
#' @return The database is written and saved to a R data file and can be read back to R for further inspection or analysis
#'@export
#'
#'
build.orthologDB = function(db = "all")
{
  #Download the gene orthologs file from NCBI
  ncbi.ortho_data = NULL
  if ("ncbi" %in% db | "all" %in% db)
  {
    ncbi.ortho = tempfile()
    downloadFiles(
      "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz",
      dest.path = paste0(ncbi.ortho, ".gz"),
      unzip = T,
      verbose = T,
      message = "Downloading ortholog file from NCBI"
    )
    ncbi.ortho_data = read.delim(ncbi.ortho,
                                 header = T,
                                 stringsAsFactors = F)
    colnames(ncbi.ortho_data) = c("tax_id",
                                  "GeneID",
                                  "relationship",
                                  "Other_tax_id",
                                  "Other_GeneID")
    ncbi.ortho_data$relationship = NULL
  }

  mgi.orths = NULL
  if ("mgi" %in% db | "all" %in% db)
  {
    # Download Human-Mouse orthologs from MGI
    temp.file = tempfile()
    downloadFiles(
      "http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt",
      verbose = T,
      dest.path = temp.file,
      unzip = F,
      message = "Downloading ortholog map from MGI"
    )
    mgi.orths = read.delim(
      temp.file,
      header = F,
      stringsAsFactors = F,
      strip.white = T
    )[, c("V2", "V6")]
    mgi.IDs = convertGeneIds(mgi.orths$V6, from = "MGI_id", to = "entrez_id")
    mgi.IDs = mgi.IDs[!is.na(unlist(mgi.IDs, use.names = F))]
    mgi.orths = subset(mgi.orths, subset = V6 %in% names(mgi.IDs))
    mgi.orths$V6 = unlist(mgi.IDs[mgi.orths$V6], use.names = F)
    mgi.orths = cbind(
      tax_id = 9606,
      GeneID = mgi.orths[, 1],
      Other_tax_id = 10090,
      Other_GeneID = mgi.orths[, 2]
    )
  }
  hgnc.ortho = NULL
  if ("hgnc" %in% db | "all" %in% db)
  {
    # Get ortho file from HGNC Comparison of Orthology Predictions (HCOP)genenames.org/tools/hcop/#!/
    hgnc.ortho = readRDS(paste0(get.path.data(), "/hcop_hs_mm.RDS"))[, c(1, 8)]
    hgnc.ortho = cbind(
      tax_id = 9606,
      GeneID = hgnc.ortho[, 1],
      Other_tax_id = 10090,
      Other_GeneID = hgnc.ortho[, 2]
    )
  }
  ortho.data = rbind(ncbi.ortho_data, mgi.orths, hgnc.ortho)
  cat("Number of orthologs identified ", nrow(ortho.data), " \n")
  ortho.data = ortho.data[!duplicated.data.frame(ortho.data), ]
  cat("Number of orthologs identified ", nrow(ortho.data), " \n")
  saveRDS(ortho.data, file = paste0(get.path.data(), "/orthodata.RDS"))
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
#' @description Function used to read in gmt file from the Broad Institute
#' @keywords gmt,read,file,broad
#' @param gmt.file Path to gmt file
#' @export
#' @return Returns a list of pathways and corresponding gene IDs
read.gmt = function(gmt.file)
{
  num_cols = max(count.fields(gmt.file, sep = "\t"))
  gmt_data = read.table(
    gmt.file,
    fill = T,
    header = F,
    sep = "\t",
    row.names = 1,
    col.names = 1:num_cols,
    strip.white = T,
    stringsAsFactors = F
  )
  gmt_data = apply(gmt_data, 1, function(x)
  {
    x = trimws(as.vector(x))
    x = x[!is.na(x) & !grepl("http", x)]
  })
  return(gmt_data)
}


#' @title Function to convert human gene symbols to mouse gene symbols using the biomaRt package
#' @author David Chisanga
#' @description The function is used to convert human gene symbols to their corresponding mouse symbols based on information in Ensembl database
#' @param species Character vector of human gene symbols
#' @export
#' @return Returns a data frame matching the input
convertHumanGeneList = function(symbol) {
  require("biomaRt")
  human = useMart(biomart = "ensembl",
                  dataset = "hsapiens_gene_ensembl",
                  host = "https://dec2021.archive.ensembl.org/")
  mouse = useMart(biomart = "ensembl",
                  dataset = "mmusculus_gene_ensembl",
                  host = "https://dec2021.archive.ensembl.org/")
  genesV2 = getLDS(
    attributes = c("hgnc_symbol"),
    filters = "hgnc_symbol",
    values = symbol ,
    mart = human,
    attributesL = c("mgi_symbol"),
    martL = mouse,
    uniqueRows = T
  )
  return(genesV2)
}


#' @title Save gmt 2 R data file
#' @author David Chisanga
#' @description Function used to save gmt file from the Broad Institute
#' @keywords gmt,save,file,broad
#' @param gmt.file Path to gmt file
#' @param from.tax_id Taxonomy ID of specie of current MSig gmt file
#' @param to.tax_id Taxonomy ID of desires
#' @export
#' @return Returns a list of pathways and corresponding gene IDs
build.MSig2R = function(gmt.files,
                        from.tax_id = 9606,
                        to.tax_id = 10090)
{
  specie.tax = list("10090" = "Mm", "9606" = "Hs")
  path.data = get.path.data()
  return1 = NULL
  tax_ids.dif = from.tax_id != to.tax_id
  if (tax_ids.dif)
    orthologs = get.orthologs(from.tax_id, to.tax_id)
  for (gmt.file in gmt.files)
  {
    gmt.temp = read.gmt(gmt.file)
    if (from.tax_id != to.tax_id)
    {
      gmt.temp = lapply(gmt.temp, function(x)
        unlist(orthologs[names(orthologs) %in% x], use.names = F))
    }
    out.file = paste0(path.data,
                      "/",
                      gsub("gmt", specie.tax[[as.character(to.tax_id)]],
                           basename(gmt.file)),
                      "_MSig.RDS")
    saveRDS(gmt.temp, file = out.file)
    return1 = gmt.temp
  }
  cat(basename(gmt.file), " saved in ", out.file, "\n")
  return(return1)
}

#' @title Function to load Molecular signature databases
#' @author David Chisanga
#' @description The function is used to load molecular signature databases provided by the Broad Institute.
#' When run for the first time, the function first
#' @param collection character or list of collections to be loaded, by default is NULL meaning all collections are loaded.
#' Options include, "c1","c2","c3","c4","c5","c6","c7" or h or a combined list of these
#' @param species Character string of the species, currently supports human (hs) and mouse (mm)
#' @return Returns a list of paths to the inbuilt collections
#' @export
get.inbuilt.MSigs = function(collection = NULL,
                             species = "Mm",
                             show.version = T)
{
  if (!grepl("mm|hs", species, ignore.case = T))
    stop(paste0("The species ", species, " is currently not supported"))
  msigdbs = list.files(
    get.path.data(),
    pattern = paste0("*", species, "_MSig.RDS$"),
    full.names = T
  )

  if (!is.null(collection))
  {
    msigdbs = msigdbs[grepl(paste0(collection, collapse = "|"), msigdbs, ignore.case = T)]
  }
  return(msigdbs)
}

file.type <- function(path) {
  f = file(path)
  ext = summary(f)$class
  close.connection(f)
  ext
}

#'@export
build.MSig2R_symbols = function(gmt.file)
{
  num.cols = max(count.fields(gmt.file, sep = "\t"))
  gmt.data = read.table(
    gmt.file,
    fill = T,
    header = F,
    sep = "\t",
    row.names = 1,
    col.names = 1:num.cols,
    strip.white = T,
    stringsAsFactors = F
  )[, -1]
  mm.gene_info = getGeneInfoNCBI()
  gmt.data = sapply(rownames(gmt.data), function(x) {
    x = unlist(gmt.data[x, ], use.names = F)
    x = x[x != ""]
    x = mm.gene_info$GeneID[grepl(paste0(x, "$", collapse = "|"),
                                  mm.gene_info$Symbol,
                                  ignore.case = T)]
    return(x)
  }, USE.NAMES = T, simplify = F)
  return(gmt.data)
}



