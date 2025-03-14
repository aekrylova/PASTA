#' The polyAsiteAssay class
#'
#' The polyAsiteAssay object is an extended \code{\link[Signac]{ChromatinAssay}}
#' and \code{\link[Seurat]{Assay}} for the storage and analysis of single-cell
#' polyAsite data.
#'
#' @slot ranges A \code{\link[GenomicRanges]{GRanges}} object describing the
#' genomic location of features in the object
#' @slot motifs A \code{\link{Motif}} object
#' @slot fragments A list of \code{\link{Fragment}} objects.
#' @slot seqinfo A \code{\link[GenomeInfoDb]{Seqinfo}} object containing basic
#' information about the genome sequence used.
#' @slot annotation A  \code{\link[GenomicRanges]{GRanges}} object containing
#' genomic annotations
#'
#' @import Signac
#' @import Seurat
#' @importFrom methods as new
#' @name polyAsiteAssay-class
#' @rdname polyAsiteAssay-class
#' @exportClass polyAsiteAssay
#' @concept assay
#'
polyAsiteAssay <- setClass(
  Class = "polyAsiteAssay",
  contains = "ChromatinAssay"
)


#' Create polyA object
#'
#' Create a \code{\link{polyAsiteAssay}} object from a count matrix.
#' The expected format of the input matrix is features x
#' cells. A set of genomic ranges must be supplied along with the matrix, with
#' the length of the ranges equal to the number of rows in the matrix. If a set
#' of genomic ranges are not supplied, they will be extracted from the
#' row names of the matrix.
#'
#' @param counts Unnormalized data (raw counts)
#' @param min.cells Include features detected in at least this many cells.
#' Will subset the counts matrix as well.
#' To reintroduce excluded features, create a new object with a lower cutoff.
#' @param max.cells Include features detected in less than this many cells.
#' Will subset the counts matrix as well.
#' To reintroduce excluded features, create a new object with a higher cutoff.
#' This can be useful for chromatin assays where certain artefactual loci
#' accumulate reads in all cells. A percentage cutoff can also be set using
#' 'q' followed by the percentage of cells, for example 'q90' will discard
#' features detected in 90 percent of cells.
#' If NULL (default), do not apply any maximum value.
#' @param min.features Include cells where at least this many features are
#' detected.
#' @param ranges A set of \code{\link[GenomicRanges]{GRanges}} corresponding to
#' the rows of the input matrix
#' @param motifs A Motif object (not required)
#' @param fragments Path to a tabix-indexed fragments file for the data
#' contained in the input matrix. If multiple fragment files are required,
#' you can add additional \code{\link{Fragment}} object to the assay after it is
#' created using the \code{\link{CreateFragmentObject}} and
#' \code{\link{Fragments}} functions. Alternatively, a list of
#' \code{\link{Fragment}} objects can be provided.
#' @param genome A \code{\link[GenomeInfoDb]{Seqinfo}} object containing basic
#' information about the genome used. Alternatively, the name of a UCSC genome
#' can be provided and the sequence information will be downloaded from UCSC.
#' @param annotation A set of \code{\link[GenomicRanges]{GRanges}} containing
#' annotations for the genome used
#' @param sep Separators to use for strings encoding genomic coordinates.
#' First element is used to separate the chromosome from the coordinates,
#' second element is used to separate the start from end coordinate. Only
#' used if \code{ranges} is NULL.
#' @param validate.fragments Check that cells in the assay are present in the
#' fragment file.
#' @param verbose Display messages
#' @param ... Additional arguments passed to \code{\link{CreateFragmentObject}}
#'
#' @importFrom Seurat CreateAssayObject
#' @importFrom Matrix rowSums colSums
#' @importFrom GenomicRanges isDisjoint strand start end
#' @concept assay
#'
#' @export

CreatePolyAAssay <- function(
    counts,
    min.cells = 0,
    min.features = 0,
    max.cells = NULL,
    ranges = NULL,
    motifs = NULL,
    fragments = NULL,
    genome = NULL,
    annotation = NULL,
    sep = c(":",",",":"),
    validate.fragments = TRUE,
    verbose = TRUE,
    ...
) {
  if (!is.null(x = ranges)) {
    if (length(x = ranges) != nrow(x = counts)) {
      stop("Length of supplied genomic ranges does not match number
           of rows in matrix")
    }
    if (sum(!(as.character(strand(ranges)) %in% c("-", "+") ))) {
      stop("At least one feature does not have a strand, please only use features
      that have strand information")
    }
  } else {
    ranges <- FeaturesToGRanges(regions = rownames(x = counts), sep = sep)
    ranges$rownames = paste0(as.character(seqnames(ranges)),"-",
                                    as.character(start(ranges)),"-",
                                    as.character(end(ranges)))
    if (sum(!(as.character(strand(ranges)) %in% c("-", "+") ))) {
      stop("At least one feature does not have a strand, please make sure input is formatted correctly
           or provide strand infromation through the 'ranges' arguement")
    }
  }
  if (!isDisjoint(x = ranges)) {
    warning("Overlapping ranges supplied. Ranges should be non-overlapping.")
  }
  if ( length( which(duplicated(rownames(counts)) == TRUE)) > 0) {
    stop("Features must be unique.")
  }
  chrom.assay <- CreateChromatinAssay(counts = counts,
                                      ranges = ranges,
                                      motifs = motifs,
                                      fragments = fragments,
                                      genome = genome,
                                      annotation = annotation,
                                      sep = sep,
                                      min.cells = min.cells,
                                      min.features = min.features,
                                      max.cells = max.cells,
                                      validate.fragments = validate.fragments,
                                      verbose = verbose)

  # requires input to be a ChromatinAssay
  pA.assay <- as(chrom.assay, Class = "polyAsiteAssay")
  pA.assay <- AddMetaData( object = pA.assay , metadata = as.character(strand(pA.assay@ranges)), col.name = "strand" )
  nCells_feature = rowSums( GetAssayData(pA.assay) > 0 )
  pA.assay <- AddMetaData( object = pA.assay , metadata = nCells_feature, col.name = "nCells_feature" )
  return(pA.assay)
}

#' Merge polyA site assays
#'
#' @export
#' @concept objects
#' @method merge polyAsiteAssay
#'
merge.polyAsiteAssay <- function(x = NULL,
                                 y = NULL,
                                 add.cell.ids = NULL,
                                 cells= NULL, ...) {
  chromatin.x <- as(object = x, Class = 'ChromatinAssay')
  if (is.list(y)) {
    chromatin.y <- list()
    for (i in 1:length(y)) {
      chromatin.y[[i]] <- as(object = y[[i]], Class = 'ChromatinAssay')
    }
  } else {
    chromatin.y <- as(object = y, Class = 'ChromatinAssay')
  }
  chromatin.m <- merge(x = chromatin.x, y = chromatin.y,
                       add.cells.ids = add.cell.ids, ...)
  polyA.m <- as(object = chromatin.m, Class = 'polyAsiteAssay')

  #add back in meta features
  meta.x <- data.frame(strand = chromatin.x@meta.features$strand)
  meta.x$peak.tmp <- rownames(chromatin.x$counts)
  if (is.list(chromatin.y)) {
    meta.y <- lapply(chromatin.y, function(chromatin) {
      df <- data.frame(strand = chromatin@meta.features$strand)
      df$peak.tmp <- rownames(chromatin$counts)
      return(df)
    })
    meta.y <- do.call(rbind, meta.y)
    meta.y <- unique(meta.y)
  } else {
    meta.y <- data.frame(strand = chromatin.y@meta.features$strand)
    meta.y$peak.tmp <- rownames(chromatin.y$counts)
  }

  meta.merge <- merge(meta.x, meta.y, by="peak.tmp", all=TRUE)
  if (any(!is.na(meta.merge$strand.x) & !is.na(meta.merge$strand.y) &
           meta.merge$strand.x != meta.merge$strand.y)) {
     warn(message = "Mismatch in strand values for the same feature when merging,
          converting strand to * for that feature")
  }


  meta.merge$strand <- ifelse(is.na(meta.merge$strand.x), as.character(meta.merge$strand.y),
                              as.character(meta.merge$strand.x))
  meta.merge$strand[meta.merge$strand.x != meta.merge$strand.y] <- "*"
  meta.merge <- meta.merge[,c("peak.tmp", "strand")]

  #check for duplicates
  duplicates <- duplicated(meta.merge$peak.tmp) | duplicated(meta.merge$peak.tmp, fromLast = TRUE)
  meta.merge$strand[duplicates] <- "*"
  meta.merge <- unique(meta.merge)

  rownames(meta.merge) <- meta.merge$peak.tmp
  meta.merge <- meta.merge[rownames(chromatin.m),]
  meta.merge$peak.tmp <- NULL
  BiocGenerics::strand(polyA.m@ranges) <- meta.merge$strand
  polyA.m <- AddMetaData(polyA.m, meta.merge)
  return(polyA.m)
}


#' Subset a polyA site assay
#'
#' @param x A polyAsiteAssay
#' @param features Which features to retain
#' @param cells Which cells to retain
#'
#' @export
#' @concept objects
#' @method subset polyAsiteAssay
#'
subset.polyAsiteAssay <- function(x,
                                  features = NULL,
                                  cells= NULL, ...) {
  chromatin <- as(object = x, Class = 'ChromatinAssay')
  chromatin <- subset(x = chromatin, cells = cells, features = features, ...)
  chromatin <- as(object = chromatin, Class = 'polyAsiteAssay')

  # Do center.scale.data slot subsetting
  #if (dim(x@center.scale.data)[1] >0 ) {
  #  center.scale <- GetAssayData(x, slot = "center.scale.data" )
  #  center.scale <- center.scale[features, cells]
  #  chromatin <- SetAssayData(chromatin, slot = "center.scale.data", new.data = center.scale)
  #}
  return(chromatin)
}

## S4 methods

setMethod(
  f = "show",
  signature = "polyAsiteAssay",
  definition = function(object) {
    cat(
      "polyAsiteAssay data with",
      nrow(x = object),
      "features for",
      ncol(x = object),
      "cells\n"
    )
    cat(
      "Variable features:",
      length(x = VariableFeatures(object = object)),
      "\n"
    )
    cat(
      "Genome:",
      unique(x = genome(x = object)),
      "\n"
    )
    cat(
      "Annotation present:",
      ifelse(
        test = is.null(x = Annotation(object = object)), yes = FALSE, no = TRUE
      ),
      "\n"
    )
    cat(
      "Motifs present:",
      ifelse(
        test = is.null(x = Motifs(object = object)),
        yes = FALSE,
        no = TRUE
      ),
      "\n"
    )
    cat(
      "Fragment files:",
      length(x = Fragments(object = object)),
      "\n"
    )
  }
)







