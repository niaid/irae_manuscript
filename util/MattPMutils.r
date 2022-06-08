# some utility funs 
# adjust p values from correlation matrix 
# given some correlation matrix object returned by hmisc::rcorr(d)
# d containd columns; all pairwise correlations computed 
p.adjust.cormat = function(hmisc.cor, method = 'fdr'){ 
  stopifnot(isTRUE(isSymmetric(hmisc.cor$P)))
  p.adj =  p.adjust(hmisc.cor$P[lower.tri(hmisc.cor$P)], method = method)
  p.adj.mx <- matrix(rep(0,ncol(hmisc.cor$P)*ncol(hmisc.cor$P)), nrow = ncol(hmisc.cor$P))
  p.adj.mx[lower.tri(p.adj.mx)] <- p.adj
  p.adj.mx[upper.tri(p.adj.mx)] <- p.adj
  diag(p.adj.mx) = 1
  colnames(p.adj.mx) = rownames(p.adj.mx) = colnames(hmisc.cor$P)
  return(p.adj.mx)
}



# Visualizations

# add alpha transparency to color 
# taken from rethinking package
# https://github.com/rmcelreath/rethinking/
col.alpha <- function( acol , alpha=0.2 ) {
  acol <- grDevices::col2rgb(acol)
  acol <- grDevices::rgb(acol[1]/255,acol[2]/255,acol[3]/255,alpha)
  acol
}

# blue red palette white in middle for low mid high diverging at mid point
blue.red = c("#053061", "#1E61A5", "#3C8ABE", "#7CB7D6", "#BAD9E9", "#E5EEF3", 
             "#F9EAE1", "#F9C7AD", "#EB9273", "#CF5246", "#AB1529", "#67001F")