F_PrintVec <- function(vec, linebreak=T){
   p = length(vec)
   invisible(sapply(1:(p-1), function(j){cat(vec[j], ' & ')}))
   if(linebreak){cat(vec[p], '\\\\ \n')}else{cat(vec[p], '\n')}
}

F_PrintTab <- function(tab){
   n = nrow(tab); p = ncol(tab)
   cat(' & ')
   F_PrintVec(colnames(tab))
   cat('\\hline \n')
   if(n > 1){
      invisible(sapply(1:(n-1), function(i){cat(rownames(tab)[i], ' & '); F_PrintVec(tab[i, ])}))
   }
   cat(rownames(tab)[n], ' & '); F_PrintVec(tab[n, ], linebreak=F)
}
