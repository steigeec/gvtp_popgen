#!/usr/bin/env Rscript
# Writes N single-column phenotype files, each a random permutation of one
# column of a GEMMA phenotype file (row order = individual order, matching
# the .fam/bfile), for use as a permutation null in gemma -lmm.
#
# Usage: Rscript make_permutations.R <env.txt> <col_idx> <n_perm> <out_dir> <seed_base>

args <- commandArgs(trailingOnly = TRUE)
env_path  <- args[1]
col_idx   <- as.integer(args[2])
n_perm    <- as.integer(args[3])
out_dir   <- args[4]
seed_base <- as.integer(args[5])

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
env <- read.table(env_path, header = FALSE)
y <- env[[col_idx]]

for (i in seq_len(n_perm)) {
  set.seed(seed_base + i)
  y_perm <- sample(y)
  write.table(y_perm, file.path(out_dir, paste0("perm_", i, ".txt")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}
cat("Wrote", n_perm, "permutations to", out_dir, "\n")
