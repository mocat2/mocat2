# DESC sums feature values per sample
# PRODUCES xls
# PRODUCES csv

sum <- t(as.matrix(rowSums(DATA)))                     # Transpose
rownames(sum) <- TYPE                                  # Add rownames
RESULTS[['sum']] <- rbind(RESULTS[['sum']], sum)       # Add to results list
