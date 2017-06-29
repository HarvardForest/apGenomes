### https://cran.r-project.org/web/packages/xtable/vignettes/xtableGallery.pdf
italic <- function(x){paste0('{\\emph{',x,'}}')}

### Example.
## mtcars$cyl <- factor(mtcars$cyl, levels = c("four","six","eight"),
##                      labels = c("four",italic("six"),"eight"))
## tbl <- ftable(mtcars$cyl, mtcars$vs, mtcars$am, mtcars$gear,
##               row.vars = c(2, 4),
##               dnn = c("Cylinders", "V/S", "Transmission", "Gears"))
## xftbl <- xtableFtable(tbl, method = "row.compact")
## print.xtableFtable(xftbl,
##                    sanitize.rownames.function = large,
##                    sanitize.colnames.function = bold,
##                    rotate.colnames = TRUE,
##                    rotate.rownames = TRUE)
