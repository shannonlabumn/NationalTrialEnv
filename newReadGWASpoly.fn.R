read.GWASpoly.object <- function (ploidy, phenos, genos, format, n.traits, delim = ",") 
{
  if (format == "ACTG") {
    format <- "ACGT"
  }
  if (!is.element(format, c("AB", "numeric", "ACGT"))) {
    stop("Invalid genotype format.")
  }
  bases <- c("A", "C", "G", "T")
  get.ref <- function(x, format) {
    if (format == "numeric") {
      ref.alt <- c(0, 1)
    }
    if (format == "AB") {
      ref.alt <- c("A", "B")
    }
    if (format == "ACGT") {
      y <- paste(na.omit(x), collapse = "")
      ans <- apply(array(bases), 1, function(z, y) {
        length(grep(z, y, fixed = T))
      }, y)
      if (sum(ans) > 2) {
        stop("Error in genotype matrix: More than 2 alleles")
      }
      if (sum(ans) == 2) {
        ref.alt <- bases[which(ans == 1)]
      }
      if (sum(ans) == 1) {
        ref.alt <- c(bases[which(ans == 1)], NA)
      }
    }
    return(ref.alt)
  }
  if(is_tibble(genos)){
    genos <- as.data.frame(genos)
  }
  if(is.data.frame(genos)){
    geno <- genos
  }
  else{
    geno <- read.table(file = geno.file, header = T, as.is = T, 
                        check.names = F, sep = delim)
  }
  map <- data.frame(Marker = geno[, 1], 
                    Chrom = factor(geno[, 2], ordered = T), 
                    Position = geno[, 3], stringsAsFactors = F)
  markers <- as.matrix(geno[, -(1:3)])
  rownames(markers) <- geno[, 1]
  tmp <- apply(markers, 1, get.ref, format)
  map$Ref <- tmp[1, ]
  map$Alt <- tmp[2, ]
  if (is.element(format, c("AB", "ACGT"))) {
    M <- apply(cbind(map$Ref, markers), 1, function(x) {
      y <- gregexpr(pattern = x[1], text = x[-1], fixed = T)
      ans <- as.integer(lapply(y, function(z) {
        ifelse(z[1] < 0, ploidy, ploidy - length(z))
      }))
      return(ans)
    })
  }
  else {
    M <- t(markers)
  }
  gid.geno <- colnames(geno)[-(1:3)]
  rownames(M) <- gid.geno
  stopifnot(na.omit(M <= ploidy & M >= 0))
  MAF <- apply(M, 2, function(x) {
    AF <- mean(x, na.rm = T)/ploidy
    MAF <- ifelse(AF > 0.5, 1 - AF, AF)
  })
  polymorphic <- which(MAF > 0)
  M <- M[, polymorphic]
  map <- map[polymorphic, ]
  map <- map[order(map$Chrom, map$Position), ]
  M <- M[, map$Marker]
  m <- nrow(map)
  cat(paste("Number of polymorphic markers:", m, "\n"))
  impute.mode <- function(x) {
    ix <- which(is.na(x))
    if (length(ix) > 0) {
      x[ix] <- as.integer(names(which.max(table(x))))
    }
    return(x)
  }
  impute.mean <- function(x) {
    ix <- which(is.na(x))
    if (length(ix) > 0) {
      x[ix] <- mean(x, na.rm = T)
    }
    return(x)
  }
  missing <- which(is.na(M))
  if (length(missing) > 0) {
    if (any(as.integer(M) != as.numeric(M), na.rm = T)) {
      cat("Missing marker data imputed with population mean \n")
      M <- apply(M, 2, impute.mean)
    }
    else {
      cat("Missing marker data imputed with population mode \n")
      M <- apply(M, 2, impute.mode)
    }
  }
  if(is_tibble(phenos)){
    phenos <- as.data.frame(phenos)
  }
  if(is.data.frame(phenos)){
  pheno <- phenos
  }
  else{
    pheno <- read.table(file = pheno.file, header = T, as.is = T, 
                       check.names = F, sep = delim)
  }
  gid.pheno <- unique(pheno[, 1])
  gid <- intersect(gid.pheno, gid.geno)
  pheno <- pheno[is.element(pheno[, 1], gid), ]
  M <- M[gid, ]
  N <- length(gid)
  cat(paste("N =", N, "individuals with phenotypic and genotypic information \n"))
  n.fixed <- ncol(pheno) - n.traits - 1
  if (n.fixed > 0) {
    fixed <- data.frame(pheno[, (n.traits + 2):ncol(pheno)], 
                        stringsAsFactors = F)
    fixed.names <- colnames(pheno)[(n.traits + 2):ncol(pheno)]
    colnames(fixed) <- fixed.names
    pheno <- data.frame(pheno[, 1:(1 + n.traits)], stringsAsFactors = F)
    cat(paste("Detected following fixed effects:\n", paste(fixed.names, 
                                                           collapse = "\n"), "\n", sep = ""))
  }
  else {
    fixed <- data.frame(NULL)
  }
  traits <- colnames(pheno)[-1]
  cat(paste("Detected following traits:\n", paste(traits, 
                                                  collapse = "\n"), "\n", sep = ""))
  return(new("GWASpoly.data", map = map, pheno = pheno, fixed = fixed, 
             geno = M, ploidy = ploidy))
}
