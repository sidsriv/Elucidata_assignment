pkgname <- "cmapR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('cmapR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("align_matrices")
### * align_matrices

flush(stderr()); flush(stdout())

### Name: align_matrices
### Title: Align the rows and columns of two (or more) matrices
### Aliases: align_matrices

### ** Examples

# construct some example matrices
m1 <- matrix(rnorm(20), nrow=4)
rownames(m1) <- letters[1:4]
colnames(m1) <- LETTERS[1:5]
m2 <- matrix(rnorm(20), nrow=5)
rownames(m2) <- letters[1:5]
colnames(m2) <- LETTERS[1:4]
m1
m2

# align them, padding with NA and returning a 3D array
align_matrices(m1, m2)

# align them, not padding and retuning a list
align_matrices(m1, m2, na.pad=F, as.3D=F)




cleanEx()
nameEx("annotate.gct")
### * annotate.gct

flush(stderr()); flush(stdout())

### Name: annotate.gct
### Title: Add annotations to a GCT object
### Aliases: annotate.gct

### ** Examples

## Not run: 
##D  g <- parse.gctx('/path/to/gct/file')
##D  g <- annotate.gct(g, '/path/to/annot')
## End(Not run)




cleanEx()
nameEx("append.dim")
### * append.dim

flush(stderr()); flush(stdout())

### Name: append.dim
### Title: Append matrix dimensions to filename
### Aliases: append.dim
### Keywords: internal

### ** Examples

(filename <- cmapR:::append.dim("my.gctx.filename", matrix(nrow=10, ncol=15)))
  



cleanEx()
nameEx("check_colnames")
### * check_colnames

flush(stderr()); flush(stdout())

### Name: check_colnames
### Title: Check whether 'test_names' are columns in the 'data.frame' df
### Aliases: check_colnames

### ** Examples

check_colnames(c("pert_id", "pert_iname"), cdesc_char)            # TRUE
check_colnames(c("pert_id", "foobar"), cdesc_char, throw_error=FALSE) # FALSE, suppress error



cleanEx()
nameEx("check_dups")
### * check_dups

flush(stderr()); flush(stdout())

### Name: check_dups
### Title: Check for duplicates in a vector
### Aliases: check_dups

### ** Examples

check_dups(c("a", "b", "c", "a", "d"))



cleanEx()
nameEx("distil")
### * distil

flush(stderr()); flush(stdout())

### Name: distil
### Title: Collapse the rows or columns of a matrix using weighted
###   averaging
### Aliases: distil

### ** Examples

m <- matrix(rnorm(30), ncol=3)
distil(m)




cleanEx()
nameEx("extract.gct")
### * extract.gct

flush(stderr()); flush(stdout())

### Name: extract.gct
### Title: Exract elements from a GCT matrix
### Aliases: extract.gct

### ** Examples

# get the values for all targeted genes from a 
# dataset of knockdown experiments 
res <- extract.gct(kd_gct, row_field="pr_gene_symbol", col_field="pert_mfc_desc")
str(res)
stats::quantile(res$vals)




cleanEx()
nameEx("fix.datatypes")
### * fix.datatypes

flush(stderr()); flush(stdout())

### Name: fix.datatypes
### Title: Adjust the data types for columns of a meta data frame
### Aliases: fix.datatypes
### Keywords: internal

### ** Examples

# meta data table with all character types
str(cdesc_char)
fixed <- cmapR:::fix.datatypes(cdesc_char)
# note how some column classes have changed
str(fixed)




cleanEx()
nameEx("is.wholenumber")
### * is.wholenumber

flush(stderr()); flush(stdout())

### Name: is.wholenumber
### Title: Check if x is a whole number
### Aliases: is.wholenumber

### ** Examples

is.wholenumber(1)
is.wholenumber(0.5)



cleanEx()
nameEx("melt.gct")
### * melt.gct

flush(stderr()); flush(stdout())

### Name: melt.gct
### Title: Transform a GCT object in to a long form 'data.table' (aka
###   'melt')
### Aliases: melt.gct

### ** Examples

# simple melt, keeping both row and column meta
head(melt.gct(ds))

# update row/colum suffixes to indicate rows are genes, columns experiments
head(melt.gct(ds, suffixes = c("_gene", "_experiment")))

# ignore row/column meta
head(melt.gct(ds, keep_rdesc = FALSE, keep_cdesc = FALSE))




cleanEx()
nameEx("merge.gct")
### * merge.gct

flush(stderr()); flush(stdout())

### Name: merge.gct
### Title: Merge two GCT objects together
### Aliases: merge.gct

### ** Examples

# take the first 10 and last 10 rows of an object
# and merge them back together
(a <- subset.gct(ds, rid=1:10))
(b <- subset.gct(ds, rid=969:978))
(merged <- merge.gct(a, b, dimension="row"))




cleanEx()
nameEx("merge_with_precedence")
### * merge_with_precedence

flush(stderr()); flush(stdout())

### Name: merge_with_precedence
### Title: Merge two 'data.frame's, but where there are common fields those
###   in 'x' are retained and those in 'y' are dropped.
### Aliases: merge_with_precedence
### Keywords: internal

### ** Examples

(x <- data.table(foo=letters[1:10], bar=1:10))
(y <- data.table(foo=letters[1:10], bar=11:20, baz=LETTERS[1:10]))
# the 'bar' column from y will be dropped on merge
cmapR:::merge_with_precedence(x, y, by="foo")




cleanEx()
nameEx("na_pad_matrix")
### * na_pad_matrix

flush(stderr()); flush(stdout())

### Name: na_pad_matrix
### Title: Pad a matrix with additional rows/columns of NA values
### Aliases: na_pad_matrix

### ** Examples

m <- matrix(rnorm(10), nrow=2)
rownames(m) <- c("A", "B")
colnames(m) <- letters[1:5]
na_pad_matrix(m, row_universe=LETTERS, col_universe=letters)




cleanEx()
nameEx("parse.gctx")
### * parse.gctx

flush(stderr()); flush(stdout())

### Name: parse.gctx
### Title: Parse a GCTX file into the workspace as a GCT object
### Aliases: parse.gctx

### ** Examples

gct_file <- system.file("extdata", "modzs_n272x978.gctx", package="cmapR")
(ds <- parse.gctx(gct_file))

# matrix only
(ds <- parse.gctx(gct_file, matrix_only=TRUE))

# only the first 10 rows and columns
(ds <- parse.gctx(gct_file, rid=1:10, cid=1:10))




cleanEx()
nameEx("parse.gmt")
### * parse.gmt

flush(stderr()); flush(stdout())

### Name: parse.gmt
### Title: Read a GMT file and return a list
### Aliases: parse.gmt

### ** Examples

gmt_path <- system.file("extdata", "query_up.gmt", package="cmapR")
gmt <- parse.gmt(gmt_path)
str(gmt)




cleanEx()
nameEx("parse.gmx")
### * parse.gmx

flush(stderr()); flush(stdout())

### Name: parse.gmx
### Title: Read a GMX file and return a list
### Aliases: parse.gmx

### ** Examples

gmx_path <- system.file("extdata", "lm_probes.gmx", package="cmapR")
gmx <- parse.gmx(gmx_path)
str(gmx)




cleanEx()
nameEx("parse.grp")
### * parse.grp

flush(stderr()); flush(stdout())

### Name: parse.grp
### Title: Read a GRP file and return a vector of its contents
### Aliases: parse.grp

### ** Examples

grp_path <- system.file("extdata", "lm_epsilon_n978.grp", package="cmapR")
values <- parse.grp(grp_path)
str(values)



cleanEx()
nameEx("process_ids")
### * process_ids

flush(stderr()); flush(stdout())

### Name: process_ids
### Title: Return a subset of requested GCTX row/colum ids out of the
###   universe of all ids
### Aliases: process_ids
### Keywords: internal

### ** Examples

gct_file <- system.file("extdata", "modzs_n272x978.gctx", package="cmapR")
ids <- read.gctx.ids(gct_file)
processed_ids <- cmapR:::process_ids(ids[1:10], ids)
str(processed_ids)




cleanEx()
nameEx("rank.gct")
### * rank.gct

flush(stderr()); flush(stdout())

### Name: rank.gct
### Title: Convert a GCT object's matrix to ranks
### Aliases: rank.gct

### ** Examples

(ranked <- rank.gct(ds, dim="column"))
# scatter rank vs. score for a few columns
plot(ds@mat[, 1:3], ranked@mat[, 1:3],
  xlab="score", ylab="rank")




cleanEx()
nameEx("read.gctx.ids")
### * read.gctx.ids

flush(stderr()); flush(stdout())

### Name: read.gctx.ids
### Title: Read GCTX row or column ids
### Aliases: read.gctx.ids

### ** Examples

gct_file <- system.file("extdata", "modzs_n272x978.gctx", package="cmapR")
# row ids
rid <- read.gctx.ids(gct_file)
head(rid)
# column ids
cid <- read.gctx.ids(gct_file, dimension="column")
head(cid)




cleanEx()
nameEx("read.gctx.meta")
### * read.gctx.meta

flush(stderr()); flush(stdout())

### Name: read.gctx.meta
### Title: Parse row or column metadata from GCTX files
### Aliases: read.gctx.meta

### ** Examples

gct_file <- system.file("extdata", "modzs_n272x978.gctx", package="cmapR") 
# row meta
row_meta <- read.gctx.meta(gct_file)
str(row_meta)
# column meta
col_meta <- read.gctx.meta(gct_file, dimension="column")
str(col_meta)
# now for only the first 10 ids
col_meta_first10 <- read.gctx.meta(gct_file, dimension="column", ids=col_meta$id[1:10])
str(col_meta_first10)




cleanEx()
nameEx("robust.zscore")
### * robust.zscore

flush(stderr()); flush(stdout())

### Name: robust.zscore
### Title: Compoute robust z-scores
### Aliases: robust.zscore

### ** Examples

(x <- rnorm(25))
(robust.zscore(x))

# with min_mad
(robust.zscore(x, min_mad=1e-4))




cleanEx()
nameEx("subset.gct")
### * subset.gct

flush(stderr()); flush(stdout())

### Name: subset.gct
### Title: Subset a gct object using the provided row and column ids
### Aliases: subset.gct

### ** Examples

# first 10 rows and columns by index
(a <- subset.gct(ds, rid=1:10, cid=1:10))

# first 10 rows and columns using character ids
(b <- subset.gct(ds, rid=ds@rid[1:10], cid=ds@cid[1:10]))

identical(a, b) # TRUE




cleanEx()
nameEx("threshold")
### * threshold

flush(stderr()); flush(stdout())

### Name: threshold
### Title: Threshold a numeric vector
### Aliases: threshold

### ** Examples

x <- rnorm(20)
threshold(x, -0.1, -0.1)




cleanEx()
nameEx("transpose.gct")
### * transpose.gct

flush(stderr()); flush(stdout())

### Name: transpose.gct
### Title: Transpose a GCT object
### Aliases: transpose.gct

### ** Examples

transpose.gct(ds)




cleanEx()
nameEx("update.gctx")
### * update.gctx

flush(stderr()); flush(stdout())

### Name: update.gctx
### Title: Update the matrix of an existing GCTX file
### Aliases: update.gctx

### ** Examples

## Not run: 
##D m <- matrix(rnorm(20), nrow=10)
##D # update by integer indices
##D update.gctx(m, ofile="my.gctx", rid=1:10, cid=1:2)
##D # update by character ids
##D row_ids <- letters[1:10]
##D col_ids <- LETTERS[1:2]
##D update.gctx(m, ofile="my.gctx", rid=row_ids, cid=col_ids)
## End(Not run)



cleanEx()
nameEx("write.gct")
### * write.gct

flush(stderr()); flush(stdout())

### Name: write.gct
### Title: Write a GCT object to disk in GCT format
### Aliases: write.gct

### ** Examples

## Not run: 
##D write.gct(ds, "dataset", precision=2)
## End(Not run)



cleanEx()
nameEx("write.gctx")
### * write.gctx

flush(stderr()); flush(stdout())

### Name: write.gctx
### Title: Write a GCT object to disk in GCTX format
### Aliases: write.gctx

### ** Examples

## Not run: 
##D # assume ds is a GCT object
##D write.gctx(ds, "my/desired/outpath/and/filename")
## End(Not run)



cleanEx()
nameEx("write.gctx.meta")
### * write.gctx.meta

flush(stderr()); flush(stdout())

### Name: write.gctx.meta
### Title: Write a 'data.frame' of meta data to GCTX file
### Aliases: write.gctx.meta
### Keywords: internal

### ** Examples

## Not run: 
##D # assume ds is a GCT object
##D cmapR:::write.gctx.meta("/my/file/path", cdesc_char, dimension="col")
## End(Not run)



cleanEx()
nameEx("write.gmt")
### * write.gmt

flush(stderr()); flush(stdout())

### Name: write.gmt
### Title: Write a nested list to a GMT file
### Aliases: write.gmt

### ** Examples

## Not run: 
##D write.gmt(gene_set, "gene_set.gmt")
## End(Not run)




cleanEx()
nameEx("write.grp")
### * write.grp

flush(stderr()); flush(stdout())

### Name: write.grp
### Title: Write a vector to a GRP file
### Aliases: write.grp

### ** Examples

## Not run: 
##D write.grp(letters, "letter.grp")
## End(Not run)




cleanEx()
nameEx("write.tbl")
### * write.tbl

flush(stderr()); flush(stdout())

### Name: write.tbl
### Title: Write a 'data.frame' to a tab-delimited text file
### Aliases: write.tbl

### ** Examples

## Not run: 
##D write.tbl(cdesc_char, "col_meta.txt")
## End(Not run)




### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
