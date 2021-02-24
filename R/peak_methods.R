#' @name peak-methods
#'
#' @title Peak into \code{\link[phyloseq]{phyloseq}} objects
#'
#' @description
#' These functions work on \code{otu_table}, \code{tax_table}, \code{sample_data}
#' or on \code{data.frame} and \code{matrix}.
#'
#' \code{peak_abundance} returns, user specified number of rows and columns
#' for \code{otu_table}.
#' 
#' #' \code{peak_taxonomy} returns, user specified number of rows and columns
#' for \code{tax_table}.
#' 
#' \code{peak_sample} returns, user specified number of rows and columns
#' for \code{sample_data}.  
#' 
#' \code{peak_base} returns, user specified number of rows and columns
#' for \code{data.frame} and \code{matrix}.
#'
#' @param x a
#'   \code{\link[phyloseq]{phyloseq}} or \code{data.frame} or
#'   \code{matrix} object 
#'   
#' @param nrows number of rows, to be specified as numeric e.g. 1, or sequence 
#'              of numeric specified as 1:5. to return first to fifth row.
#'
#' @param ncols number of cols, to be specified as numeric e.g. 1, or sequence 
#'              of numeric specified as 1:5 to return first to fifth col.
#'             
#' @return Print user specified rows and columns
#' 
#' @examples
#' data("zackular2014")
#' 
#' peak_abundance(zackular2014, nrows=1:3, ncols = 1:3)
#' 
#' peak_taxonomy(zackular2014, nrows=1:3, ncols = 1:3)
#' 
#' peak_sample(zackular2014, nrows=1:3, ncols = 1:3)
#' 
#' dat.frm <- meta(zackular2014)
#' # specify specific columns
#' peak_base(dat.frm, nrows=1:3, ncols = c(1, 3, 4))
#' 
#' matrix_ab <- abundances(zackular2014)
#' peak_base(matrix_ab, nrows=1:3, ncols = 1:3)
#' 
NULL

#' @rdname peak-methods
setGeneric("peak_abundance", signature = c("x"),
           function(x, 
                    nrows=1:5,
                    ncols=1:5)
             standardGeneric("peak_abundance"))  


#' @rdname peak-methods
setGeneric("peak_taxonomy", signature = c("x"),
           function(x, 
                    nrows=1:5,
                    ncols=1:5)
             standardGeneric("peak_taxonomy"))

#' @rdname peak-methods
setGeneric("peak_sample", signature = c("x"),
           function(x, 
                    nrows=1:5,
                    ncols=1:5)
             standardGeneric("peak_sample"))

#' @rdname peak-methods
setGeneric("peak_base", signature = c("x"),
           function(x, 
                    nrows=1:5,
                    ncols=1:5)
             standardGeneric("peak_base"))


#' @rdname peak-methods
#' @aliases peak_abundance
#'
#' @export
setMethod("peak_abundance", signature = c(x = "phyloseq"),
          function(x, 
                   nrows=1:5,
                   ncols=1:5){
            
            if(!.check_pseq(x)){
              stop("input must be a phyloseq object for peak_abundance.",
                   call. = FALSE)
              #abundances(x)[nrows,ncols]
            }
            
            abundances(x)[nrows,ncols]
            
          }
)


#' @rdname peak-methods
#' @aliases peak_taxonomy
#'
#' @export
setMethod("peak_taxonomy", signature = c(x = "phyloseq"),
          function(x, 
                   nrows=1:5,
                   ncols=1:5){
            
            if(!.check_pseq(x)){
              stop("input must be a phyloseq object for peak_taxonomy",
                   call. = FALSE)
              #abundances(x)[nrows,ncols]
            }
            
            tax_table(x)[nrows,ncols]
            
          }
)


#' @rdname peak-methods
#' @aliases peak_sample
#'
#' @export
setMethod("peak_sample", signature = c(x = "phyloseq"),
          function(x, 
                   nrows=1:5,
                   ncols=1:5){
            
            if(!.check_pseq(x)){
              stop("input must be a phyloseq object for peak_sample",
                   call. = FALSE)
              #abundances(x)[nrows,ncols]
            }
            
            meta(x)[nrows,ncols]
            
          }
)


#' @rdname peak-methods
#' @aliases peak_base
#'
#' @export
setMethod("peak_base", signature = c(x = "ANY"),
          function(x, 
                   nrows=1:5,
                   ncols=1:5){
            
            if(!is.data.frame(x) & !is.matrix(x) ){
              
              stop("input for peak_base must be data.frame or matrix",
                   call. = FALSE)
             
            }
            
            x[nrows,ncols]
            
          }
)


#' @param x object to test
#' @noRd
.check_pseq <- function (x) {
  length(x) == 1 && is(x) == "phyloseq"
}
