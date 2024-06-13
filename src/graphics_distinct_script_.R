geno_umap_2d$ID <- 1:nrow(geno_umap_2d)

ggplot(geno_umap_2d, aes(x = X1, y = X2)) +
  geom_text(aes(label = ID), vjust = -0.5, hjust = 0.5, 
            color = "blue") +
  ggtitle("UMAP 2D plot for REFPOP genotype data") +
  xlab("first component") +
  ylab("second component") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

geno_pca_mat_$ID <- 1:nrow(geno_pca_mat_)

ggplot(geno_pca_mat_, aes(x = PC1, y = PC2)) +
  geom_text(aes(label = ID), vjust = -0.5, hjust = 0.5, 
            color = "blue") +
  ggtitle("PCA 2D plot for REFPOP genotype data") +
  xlab(paste0(
    names(geno_pca_exp_var_)[1], ": ",
    signif(100 * as.numeric(geno_pca_exp_var_)[1], 2), "%"
  )) +
  ylab( paste0(
    names(geno_pca_exp_var_)[2], ": ",
    signif(100 * as.numeric(geno_pca_exp_var_)[2], 2), "%"
  )) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))