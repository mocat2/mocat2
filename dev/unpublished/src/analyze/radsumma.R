
# DESC summerar rader
# PRODUCES csv

res <- colSums(DATA)

file.result <- t(as.matrix(res))                       # Transpose
rownames(file.result) <- TYPE                                  # Add rownames
results[['FINradSUMMA']] <- rbind(results[['FINradSUMMA']], file.result)   # Add to results list

