setwd("path_to_working_directory")
library("seqinr") 
library("NMF")
library("reshape2")
library("RColorBrewer")
library("ggtree")
library("ape")
library("phangorn")
library("plyr")
library("tidyverse")
library("ggtreeExtra")
library("ggnewscale")
library("fastbaps")
library("rPinecone")
library("pairsnp")
library("ggh4x")
library("adegenet")
library("TreeTools")
library("PanVizGenerator")
library("ggstance")
devtools::install_github("gtonkinhill/pairsnp-r")
library("pairsnp-r")

## SNPs analysis

## Use of alignment done with reference M avium subsp avium H87 (formerly hominissuis) Closest genome from ncbi
## We import the alignment from the map2ref workflow into a sparse matrix
sparse.data <- import_fasta_sparse("clean_core.aln")  ## The clean_core.aln file is the output from snp-sites after gubbins

## Calculate SNp matrix
snp.matrix <- snp_dist(sparse.data) 

check_tree <- read.tree("RAxML_bipartitions.SNP") ## RAxML tree from snp-sites
ggtree(check_tree) + geom_tiplab()

p16 <- read.tree("RAxML_bipartitions.SNP")
p16 <- midpoint(p16)


# Import annotation files --> csv file with sample names same as fasta headers, patient clusters, SKA results, etc.
annotation <- read_csv(file = "info.csv",
                       col_names = T,
                       na = "")
order_tree <- p16$tip.label
sort(order_tree)
sort(annotation$ID)
annotation

## check
dd <- annotation
dd4 <- dd %>% arrange(desc(factor(ID, levels = order_tree)))

# orig_tiplabs <- c(p16$tip.label) ## This is in case the names are different
# new_tiplab <- rev(c(dd4$New_name))

p16$tip.label <- new_tiplab

ggtree(p16) + geom_tiplab()


# Export matrix to csv
write.table(snp.matrix,file="snp_difference_matrix.tsv",sep = "\t")

## Import tree
tree_wgs_i <- ggtree(p16) + 
  geom_text2(aes(subset=!is.na(as.numeric(label)) & as.numeric(label) < 75, 
                 label=label, 
                 hjust=1.5,
                 vjust=-.5), size=2) + 
  geom_treescale(width = .05)

tree_wgs_i + geom_tiplab(align = TRUE)

annotation
snp.matrix
# Merge and assign patients to samples
snp_dif_wgs_i <- melt(data = snp.matrix)
snp_dif_wgs_i$Var1 <- as.character(snp_dif_wgs_i$Var1)
snp_dif_wgs_i$Var2 <- as.character(snp_dif_wgs_i$Var2)
snp_dif_wgs_i <- filter(snp_dif_wgs_i, !Var1 == Var2)

snp_dif_wgs_i_edit <- bind_cols(snp_dif_wgs_i, select(
  left_join(snp_dif_wgs_i, annotation, by = c("Var1" = "ID")), 
  "Patient"))
snp_dif_wgs_i_edit <- rename(snp_dif_wgs_i_edit, Var1Patient = Patient)

snp_dif_wgs_i_edit <- bind_cols(snp_dif_wgs_i_edit, select(
  left_join(snp_dif_wgs_i_edit, annotation, by = c("Var2" = "ID")), 
  "Patient"))
snp_dif_wgs_i_edit <- rename(snp_dif_wgs_i_edit, Var2patient = Patient)

snp_dif_wgs_i_edit$patient <- if_else(
  condition = snp_dif_wgs_i_edit$Var1Patient == snp_dif_wgs_i_edit$Var2patient,
  true = "within",
  false = "between")

#### You can plot the differences in SNVs Plots 

ggplot(na.omit(snp_dif_wgs_i_edit),aes(value,fill=patient,group=patient)) + 
  geom_histogram(binwidth = 2) + 
  theme_bw() + scale_fill_brewer(type = "qual", palette = 2) 


######## Inhouse Seq Cluster Analysis ###############

assigned_cluster_wgs_i <- hclust(as.dist(snp.matrix)) ## using parsnp we have to transform matrix to a dist() object

order_cluster_labels <- assigned_cluster_wgs_i$labels
dd3 <- dd %>% arrange(desc(factor(ID, levels = order_cluster_labels)))

orig_cluster_labs <- c(assigned_cluster_wgs_i$labels)
new_cluster_labs <- rev(c(dd3$New_name))

assigned_cluster_wgs_i$labels <- new_cluster_labs

# Plot clustering tree to see relationships
plot(assigned_cluster_wgs_i, hang = -1)
# Draw tree as a dendrogram
acd_wgs_i <- as.dendrogram(assigned_cluster_wgs_i)
# Find the largest SNP difference within same patient samples
tail(sort(pull(filter(snp_dif_wgs_i_edit, patient == "within"), value)))
# Plot the dendrogram with SNP cutoff in the y axis limits to manually inspect the relationships. the y limit c(x, y) y = the upper threshold difference within the same patient, 
# in our example ty = 4
plot(acd_wgs_i, ylim = c(0, 4))


# Assign clusters by cutting tree and snp difference you have chosen to define as transmission "h = " the same number as before
clusters_wgs_i <- cutree(assigned_cluster_wgs_i, h = 4)
# Convert cluster assignments to a data frame
wgs_clusters_wgs_i <- data.frame(clusters_wgs_i)
# Add sample names from rownames
wgs_clusters_wgs_i$sample <- rownames(wgs_clusters_wgs_i)
# Remove isolates with no predicted transmission
wgs_clusters_wgs_i_trim <- wgs_clusters_wgs_i[wgs_clusters_wgs_i$clusters %in% 
                                                names(table(wgs_clusters_wgs_i$clusters))[
                                                  table(wgs_clusters_wgs_i$clusters) > 1],]
# Convert clusters to a factor
wgs_clusters_wgs_i_trim$clusters <- as.factor(wgs_clusters_wgs_i_trim$clusters)
## Annotate cluster numbers
wgs_clusters_wgs_i_trim$clusters
# Revalue the factor from numbers to Roman numerals
wgs_clusters_wgs_i_trim$clusters <- revalue(wgs_clusters_wgs_i_trim$clusters, c(
  "1"="I",
  "2"="II",
  "4"= "III"))

# Add the groups to the annotation metadata table
annotation_w_cluster <- full_join(annotation,wgs_clusters_wgs_i_trim, by = c("New_name" = "sample"))
annotation_w_cluster <- select(annotation_w_cluster, -clusters_wgs_i)
annotation_w_cluster <- rename(annotation_w_cluster, clusters_wgs_i = clusters)

# Create a list of the clusters
annotation_w_cluster$ID <- annotation_w_cluster$New_name ### To change the name
cluster_list <- split(annotation_w_cluster$ID, annotation_w_cluster$clusters_wgs_i)
# Use ggtree to add the group metadata to the tree
tree_groups_i <- groupOTU(tree_wgs_i, cluster_list) + geom_tiplab(align = T,linesize=.1,size=2)
# Replot the tree colouring the branches and tips with the transmission groups
tree_groups_i <- tree_groups_i + aes(color=group) + 
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("black", brewer.pal(7,"Dark2")))


## Annotate with patient number
p <- groupOTU(tree_wgs_i, cluster_list) + geom_tiplab(align = T,linesize=.1,size=4, offset = .01) 

### rpinecone 
library("rPinecone")
library("adegenet")

ggtree(read.tree("pyjar_SNP.joint.tre"))

pine_tree <- read.tree("pyjar_SNP.joint.tre")

ggtree(pine_tree) + geom_tiplab()
### Order the tree
order_pinetree <- pine_tree$tip.label
dd6 <- dd %>% arrange(desc(factor(ID, levels = order_pinetree)))

orig_pinetiplabs <- c(order_pinetree)
new_pinetiplab <- rev(c(dd6$New_name))

pine_tree$tip.label <- new_pinetiplab
##############################


pinerooted_tree <- midpoint(pine_tree)

pinecone.results <- pinecone(pinerooted_tree, 4, 23, quiet = TRUE)

pinecone.results$table

pine_table <- pinecone.results$table

plot.pinecone <- data.frame(pine_table)

plot.pinecone$Major.Sub.lineage <- NULL

plot.pinecone$Sub.lineage <- str_remove(plot.pinecone$Sub.lineage, "singleton_.{1,}")

plot.pinecone$Sub.lineage[plot.pinecone$Sub.lineage == ""] <- NA # We remove Singletons to clarify, so there is not too much letters
order.patients <- get_taxa_name(p)
plot.pinecone2 <- plot.pinecone %>% arrange(desc(factor(Taxa, levels = order.patients))) ## It works but is not smart


## This works to mantain aesthetics:
col.4.tree <- c("black", brewer.pal(8,"Dark2"), brewer.pal(8, "Dark2"))

names(col.4.tree) <- c(0, "I", "II", "III", 1:4) ### Very important!!

snp.matrix
annotation

dd2 <- dd %>% arrange(desc(factor(New_name, levels = order.patients))) 

dd2$ID <- dd2$New_name
dd3 <- cbind(dd2[2,], dd2[1,])

## To create the annotated tree with all information
p12 <- p + aes(colour = group) +  guides(color = guide_legend(override.aes = aes(label = ""))) + xlim_tree(xlim = .30) + 
  
  geom_facet(data = dd2, panel = "Patient", geom = geom_label2, mapping = aes(x = 0), label = dd2$Patient, size = 4,
             fill = "white", show.legend = FALSE) + 
  
  geom_facet(data = plot.pinecone2, panel = "rPinecone", geom = geom_label, mapping = aes(x = 0, color = as.factor(plot.pinecone2$Sub.lineage)), 
             label = plot.pinecone2$Sub.lineage, size = 4,  fontface = "bold", fill = "white", show.legend = FALSE) +
  
  geom_facet(data = dd2, panel = "ST", geom = geom_label, mapping = aes(x = 0, color = as.factor(dd2$ST)), 
             label = dd2$ST, size = 4,  fontface = "bold", fill = "white", show.legend = FALSE) +
  
  theme(legend.position = "right") + labs(colour = "Cluster")  + 
  
  scale_color_manual(values = col.4.tree, breaks = c("I", "II", "III")) +
  
  force_panelsizes(cols = c(1.9, .15, .15, .15), respect = F)

p12

pdf(file = "results.pdf", 
    width = 17, 
    height = 12)
p12
dev.off()

## PNG
png(file = "results.png",
    width = 17,
    height = 12)
p12
dev.off()

