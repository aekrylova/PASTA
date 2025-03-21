
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Find polyA sites that are differentially used across cells, utilizing
#' polyA residuals.
#' @param object An object
#' @param assay name of polyAsite assay to test. Default is polyA.
#' @param ident.1 Identity class to find polyA sites for
#' @param ident.2 A second identity class for comparison.
#' @param features Features to test. Default is all features. If a subset of
#' features is provided, pct.1 and pct.2 will be calculated using all features
#' with the same gene annotation where residuals have been calculated.
#' @param covariates Vector of covariates to include in linear model.
#' @param gene.names Column name providing gene annotation of each polyA site.
#' Default is "Gene_Symbol"
#'
#' @importFrom stats lm relevel
#'
#' @rdname FindDifferentialPolyA
#' @concept differential_polyA
#' @export
#'
FindDifferentialPolyA <- function(
    object,
    assay = "polyA",
    ident.1,
    ident.2,
    features = NULL,
    covariates = NULL,
    gene.names = "Gene.Symbol") {

  if( !inherits(object[[assay]], "polyAsiteAssay")){
    stop(paste0(assay," assay is not a polyAsiteAssay"))
  }

  if (dim(LayerData(object, layer="scale.data", assay = assay))[1] == 0)  {
    stop ("No features found in scale.data layer for specified assay. Run CalcPolyAResiduals prior to FindPolyASites")
  }

  if (!(gene.names %in% colnames(object[[assay]][[]]))) {
    stop("Gene.names column not found in meta.features, please make sure
         you are specific gene.names correctly")
  }

  if (!(ident.1 %in% unique(Idents(object)))) {
    stop("ident.1 not found in object")
  }
  df <- data.frame(ident = Idents(object))

  if (!is.null(covariates)) {
    for (i in 1:length(covariates)) {
      if (is.na(match(covariates[[i]], colnames(object[[]])))) {
        stop("Covariates are not found in meta data, please make sure you are specifying correctly.")
      }
      df[,i+1] <- object[[]][,match(covariates[[i]], colnames(object[[]]))]
    }
    colnames(df) <- c("ident", covariates)

  }

  features <- features %||% rownames(x = object[[assay]]@scale.data)
  r.matrix <- object[[assay]]@scale.data
  df$ident <- relevel(df$ident, ref = ident.2)

  sub <- subset(df, ident %in% c(ident.1, ident.2))
  r.matrix.sub <- r.matrix[features,rownames(sub)]

  all.models <- lapply(
    X = 1:nrow(x = r.matrix.sub),
    FUN = function(i) {
      sub$residuals <- as.numeric(r.matrix.sub[i, ])
      model <- lm(residuals ~  .,  data=sub)

      to_return <- data.frame(summary(model)$coefficients)
      to_return$coefficients <- rownames(to_return)
      colnames(to_return) <- c("Estimate", "std_error", "t", "p.value", "coefficient")
      to_return$peak <- rownames(r.matrix.sub)[i]
      return(to_return)
    }
  )
  results <- do.call(rbind, all.models)
  main.effects <- results[grep("ident", results$coefficient),]

  gene.idx = match(gene.names, colnames(object[[assay]][[]]))
  main.effects$symbol <- object[[assay]][[]][main.effects$peak,gene.idx]

  #get all features within genes to calculate percentage.usage
  #also require that they have residuals calculated
  features.all.genes <- rownames(object[[assay]][[]][object[[assay]][[]][[gene.names]] %in% unique(main.effects$symbol), ])
  features.all.genes <- intersect(features.all.genes, rownames(LayerData(object[[assay]], layer="scale.data")))
  #calculate percentage usage in group 1 and group 2
  percent.1 <- percentage.usage(object,
                                assay = assay,
                                cells = WhichCells(object, idents = ident.1),
                                features = features.all.genes,
                                gene.names = gene.names)
  main.effects$percent.1 <- percent.1[main.effects$peak]
  percent.2 <- percentage.usage(object,
                                assay = assay,
                                cells = WhichCells(object, idents = ident.2),
                                features = features.all.genes,
                                gene.names = gene.names)
  main.effects$percent.2 <- percent.2[main.effects$peak]
  main.effects$p_val_adj <- main.effects$p.value * nrow(object[[assay]]@scale.data)
  main.effects$p_val_adj[main.effects$p_val_adj  > 1] <- 1

  rownames(main.effects) <- main.effects$peak
  main.effects.return <- main.effects[,c("Estimate", "p.value", "p_val_adj", "percent.1", "percent.2", "symbol")]

  #order by p-value
  main.effects.return <- main.effects.return[ with(main.effects.return, order(p_val_adj, -Estimate)),]
  return(main.effects.return)
}


percentage.usage <- function( object,
                              assay = "polyA",
                              cells,
                              features,
                              gene.names = "Gene_Symbol") {
  df <- data.frame(peak=features)
  meta <- object[[assay]]@meta.features
  df$symbol <- meta[features,gene.names]
  df$counts1 <- rowSums(object[[assay]]@counts[df$peak, cells])
  sum1 <- aggregate(df$counts1, by=list(gene=df$symbol), FUN=sum)
  colnames(sum1) <- c("symbol", "sum")
  df <- merge(df, sum1, by="symbol")
  df$frac <- df$counts1/df$sum
  rownames(df) <- df$peak
  df <- df[features,]
  df$frac[df$sum==0] <- 0
  v <- as.vector(df$frac)
  names(v) <- features
  return(v)
}



