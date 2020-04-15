#Functions for teffects class

#' @export
print.teffect <- function(object){
  a = object
  cat("\nCall:\n", paste(deparse(a$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  coefs = data.frame(x1 = a$TE,x2 = a$sd)
  names(coefs)<- c("Mean","Std.err")
  rownames(coefs) <- a$te
  print(coefs)

}
