


variables <- c("var1", "var2", "var3", "var4", "var5", "var6", "var7", "var8", "var9", "var10")

# Generate all combinations (e.g., for size 2)
combinations_size_2 <- combn(variables, 2)

# Display the combinations
print(combinations_size_2)


# Generate all possible combinations for different sizes
all_combinations <- lapply(2:length(variables), function(size) combn(variables, size))

# Display combinations of all sizes
all_combinations




all_combinations_flat <- unlist(all_combinations, recursive = FALSE)

# Find the maximum combination size to handle varying lengths
max_size <- max(sapply(all_combinations_flat, length))

# Convert each combination into a data frame row and fill missing values with NA
combinations_df <- do.call(rbind, lapply(all_combinations_flat, function(x) {
  c(x, rep(NA, max_size - length(x)))
}))

# Convert to data.frame and set column names
combinations_df <- as.data.frame(combinations_df, stringsAsFactors = FALSE)
colnames(combinations_df) <- paste0("Var", 1:ncol(combinations_df))

# Display the final data frame
print(combinations_df)

###############


# Generate all possible combinations for different sizes
all_combinations <- lapply(2:length(variables), function(size) combn(variables, size, simplify = FALSE))

# Flatten the list of lists into a single list
all_combinations_flat <- unlist(all_combinations, recursive = FALSE)

# Find the maximum combination size to handle varying lengths
max_size <- max(sapply(all_combinations_flat, length))

# Convert each combination into a data frame row and fill missing values with NA
combinations_df <- do.call(rbind, lapply(all_combinations_flat, function(x) {
  c(x, rep(NA, max_size - length(x)))
}))

# Convert to data.frame and set column names
combinations_df <- as.data.frame(combinations_df, stringsAsFactors = FALSE)
colnames(combinations_df) <- paste0("Var", 1:ncol(combinations_df))

# Display the final data frame
print(combinations_df)




