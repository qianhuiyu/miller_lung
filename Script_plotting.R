# read all fetal data and generate the dotplot for all fetal data
library(Seurat)
source("~/Work/commonScript/Script_functions.R")
fetal <- readRDS("../Res_fetal_1500.rds")
# rearrange the order of clusters for clearer presentation
current.cluster.ids <- 26:0
new.cluster.ids <- c(7,6,23,21,26,1,27,25,5,10,9,20,8,19,18,17,4,3,24,22,2,16:11)
cell.idx <- rep(NA, nrow(fetal@meta.data))
for(new.idx in new.cluster.ids){
        old.idx=current.cluster.ids[which(new.cluster.ids==new.idx)]
                cell.idx[which(fetal@meta.data$res.0.8==old.idx)] <- new.idx

}
fetal@meta.data[,"cl_idx_2"] <- cell.idx
lineage <- rep("RBC", nrow(fetal@meta.data))
lineage[which(fetal@meta.data$cl_idx_2%in%seq(10))] <- "Epithelial"
lineage[which(fetal@meta.data$cl_idx_2%in%c(11:21))] <- "Mesenchymal"
lineage[which(fetal@meta.data$cl_idx_2==22)] <- "Endothelial"
lineage[which(fetal@meta.data$cl_idx_2%in%c(23:26))] <- "Immune"
fetal@meta.data[,"lineage"] <- lineage
saveRDS(fetal, file="Res_fetal_1500.rds")
 
cl.ave.expr <- getAveExpr(meta.data=fetal@meta.data, feature.to.calc="cl_idx_2", expr=fetal@data, genes=NULL, core.num=20, colname.prefix="Cluster")
saveRDS(cl.ave.expr, file="Res_fetal_all_cl_idx_2_cluster_average_expr.rds")

pan.epi.markers <- rev(c("EPCAM", "KRT8", "KRT18", "FXYD3", "PERP", "CDH1"))
pan.meso.markers <- c("COL6A2","MFAP4","DCN","FHL1","COL1A2","COL3A1")
pan.imm.markers <- c("CD37","CORO1A","LCP1","PTPRC","CD53","LAPTM5")
endothelial.markers <- c("FLT1", "CDH5", "CLDN5", "EGFL7", "ESAM","IFI27")
neuron.mix.markers <- c("CHGA","ASCL1","NNAT","STMN2","GRP","MPZ")
rbc.markers <- c("HBB", "HBA1", "HBA2", "HBG1", "HBG2")
g1 <- c(neuron.mix.markers, pan.epi.markers, pan.meso.markers, endothelial.markers, pan.imm.markers, rbc.markers)

scale.cl.expr <- t(scale(t(cl.ave.expr[g1, ])))
colorPal <- grDevices::colorRampPalette(c("#d9d9d9", "#252525"))
color.break.num <- 30
cellColor <- as.vector(apply(scale.cl.expr, 2, function(values){
	adjustcolor(colorPal(color.break.num), alpha=.8)[as.numeric(cut(values, breaks=color.break.num, right=F, include.lowest=T))]
}))
cl.num <- length(unique(fetal@meta.data$"cl_idx_2"))
expr.prop <- sapply(sort(unique(fetal@meta.data$"cl_idx_2")), function(i){
	cell.idx <- which(fetal@meta.data$"cl_idx_2"==i)
	mat <- as.matrix(fetal@data[g1, cell.idx])
	apply(mat, 1, function(vec){sum(vec>0)/length(vec)})
})
point.size <- as.vector(expr.prop)*6
pdf("Dotplot_fetal_all_lineage_marker.pdf", height=20, width=20)
par(mar=c(11,3,2,2))
plot(rep(seq(length(g1)), cl.num), rep(seq(cl.num), each=length(g1)), pch=16, col=cellColor, cex=point.size, bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
mtext(g1, side=1, at=seq(length(g1)), las=2, cex=3)
mtext(sort(unique(fetal@meta.data$"cl_idx_2")), side=2, at=seq(cl.num), las=1, cex=3)
dev.off()

p <- c(0, 0.25, 0.5, 0.75, 1)*6
pdf("Plot_dot_size_legend.pdf")
plot(seq(length(p)), seq(length(p)), pch=1, cex=p, bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
dev.off()

#plot fetal tSNE showing sample information
g.cols <- c("#a4d371", "#3c7ab6", "#4daa99", "#a64a97", "#c86879", "#1f7539", "#322a84", "#852655")
names(g.cols) <- paste0("S", seq(8))
sample.names <- c("11.5 week distal", "11.5 week airway", "15 week distal", "15 week airway", "15 week trachea", "18 week distal", "18 week airway", "21 week trachea")
tsne.coor <- fetal@dr$tsne@cell.embeddings
group.res <- as.character(fetal@meta.data$orig.ident)
png("Plot_tsne_fetal_all_no_sample_identity.png", height=20, width=40, unit="cm", res=500)
par(mfrow=c(1,2))
plotFeature2(tsne.coor, values=group.res, xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", gCols=g.cols, cex=1)
plotFeature2(tsne.coor, values=group.res, xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", gCols=g.cols, type="n")
legend("topleft", legend=sample.names, text.col=g.cols, bty="n")
dev.off()

# plot fetal tSNE showing cluster information
unique(fetal@meta.data$lineage)
epi.cl <- unique(fetal@meta.data$cl_idx_2[which(fetal@meta.data$lineage=="Epithelial")])
epi.cols <- colorRampPalette(c("#fee0d2", "#fb6a4a","#a50f15"))(length(epi.cl))
names(epi.cols) <- paste0("Cluster", epi.cl)

mes.cl <- unique(fetal@meta.data$cl_idx_2[which(fetal@meta.data$lineage=="Mesenchymal")])
mes.cols <- colorRampPalette(c("#e5f5e0", "#74c476","#006d2c"))(length(mes.cl))
names(mes.cols) <- paste0("Cluster", mes.cl)

imm.cl <- unique(fetal@meta.data$cl_idx_2[which(fetal@meta.data$lineage=="Immune")])
imm.cols <- colorRampPalette(c("#deebf7", "#6baed6","#08519c"))(length(imm.cl))
names(imm.cols) <- paste0("Cluster", imm.cl)

endo.cols <- "#6a51a3"
names(endo.cols) <- paste0("Cluster", unique(fetal@meta.data$cl_idx_2[which(fetal@meta.data$lineage=="Endothelial")]))
rbc.cols <- "#bdbdbd"
names(rbc.cols) <- paste0("Cluster", unique(fetal@meta.data$cl_idx_2[which(fetal@meta.data$lineage=="RBC")]))
cell.type.cols <- c(epi.cols, mes.cols, imm.cols, endo.cols, rbc.cols)
group.res <- paste0("Cluster", as.character(fetal@meta.data$cl_idx_2))
png("Plot_tSNE_fetal_all_cluster_sequential_cols.png", height=20, width=20, unit="cm", res=500)
plotFeature2(tsne.coor, values=group.res, xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", cex=1,gCols=cell.type.cols)
for(i in seq(27)){
        idx <- which(fetal@meta.data[,"cl_idx_2"]==i)
        coor <- apply(tsne.coor[idx,],2,median)
        text(coor[1], coor[2], labels=i, cex=1.5)
}
dev.off()

gene.cols <- c("#ffffcc", "#a1dab4", "#41b6c4", "#2c7fb8", "#253494")
#gene.cols <- c("#feebe2", "#fbb4b9", "#f768a1", "#c51b8a", "#7a0177")
g2 <- intersect(c("EPCAM", "NKX2-1", "FOXI1", "ASCL3", "DCLK1", "TMPR5"), rownames(fetal@data))
for(x in g2){
	file.name <- paste("Plot_fetal_all_cell_", x, "-2.png")
	png(file.name, height=20, width=20, unit="cm", res=500)
	plotFeature2(tsne.coor, values=fetal@data[x,], xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", cex=1,nCols=gene.cols)
	dev.off()
}

epi <- readRDS("../Res_epi2_noCorrection_10PC_2reso.rds")
# color scheme
cl.cols <- c("#197636", "#7BCAA4", "#A4D371", "#EED951", "#DCCB7C", "#FDD884", "#4DA7B0", "#CA6778","#9E0142","#5E4FA2", "#A84798", "#3C7AB6")
names(cl.cols) <- paste0("Cluster", c(2,1,11,7,10,12,3,6,5,4,8,9))
cell.types <- c("Multiciliated cell", "Multiciliated precursor", "Intermediate ciliated", "Basal cell", "Goblet-like secretory", "Club-like secretory", "Neuroendocrine", "Bud tip adjacent", "Bud tip progenitor", "Hub cell", "Submucosal gland",  "Submucosal gland basal")
names(cell.types) <- paste0("Cluster", c(2,1,11,7,10,12,3,6,5,4,8,9))

# tSNE color coded by clusters
plot.coor=epi@dr$tsne@cell.embeddings
group.res=paste0("Cluster", epi@meta.data$merge.cl.idx)
group.cols <- cl.cols
png("Plot_tSNE_epi_noCorrection_merged_cluster_info.png", height=20, width=20, unit="cm", res=500)
plotFeature2(coor=plot.coor, values=group.res, main="", cex=1, cex.main=3, xaxt="n", yaxt="n", bty="n", xlab="", ylab="",gCols=cl.cols)
#for(x in unique(group.res)){
#        cell.idx <- which(group.res==x)
#        vec <- apply(plot.coor[cell.idx, ], 2, median)
#        text(vec[1], vec[2], label=x, cex=2)
#}
dev.off()

# tSNE color coded by samples
g.cols <- c("#a4d371", "#3c7ab6", "#4daa99", "#a64a97", "#c86879", "#1f7539", "#322a84", "#852655")
names(g.cols) <- paste0("S", seq(8))
sample.names <- c("11.5 week distal", "11.5 week airway", "15 week distal", "15 week airway", "15 week trachea", "18 week distal", "18 week airway", "21 week trachea")
plot.coor=epi@dr$tsne@cell.embeddings
group.res=epi@meta.data$orig.ident
png("Plot_tSNE_epi_noCorrection_sample_info.png", height=20, width=20, unit="cm", res=500)
plotFeature2(coor=plot.coor, values=group.res, main="", cex=1, cex.main=3, xaxt="n", yaxt="n", bty="n", xlab="", ylab="",gCols=g.cols)
dev.off()

# heatmap for expression patterns for selected top markers across clusters
cm <- readRDS("../Res_epi_merged_cluster_top50_removeCC_rbc_gene.rds")
cl.ave.expr <- readRDS("../Res_all_genes_merged_cluster_average_expression.rds")
colors <- colorRampPalette(c("#f7f7f7", "#d9d9d9", "#252525"), space="rgb")(50)
file.name <- "Plot_epi_denovo_cluster_markers_by_cluster.pdf"
pdf(file.name, height=8)
heatmap.2(cl.ave.expr[cm,],trace="none",density.info="none",scale="row",dendrogram="none", col=colors,Rowv=FALSE, Colv=FALSE, labRow=FALSE, labCol=names(cl.cols)[colnames(cl.ave.expr)], ColSideColors=cl.cols[colnames(cl.ave.expr)], srtCol=30, margins = c(8, 5))
dev.off()

known.markers <- c("FOXJ1", "CDK1", "CHGA", "SFTPB", "SCGB3A2", "SFTPC", "AGER", "TP63", "SCGB1A1", "MUC5B", "LTF", "FHL2") 
cluster.order <- sub("Cluster", "", colnames(cl.ave.expr))
e1 <- lapply(known.markers, function(x){
	lapply(cluster.order, function(i){
		epi@data[x,which(epi@meta.data$merge.cl.idx==i)]
	})	
})
pdf("Plot_boxplot_known_and_denovo_markers_epi_clusters.pdf", height=50)
par(mfrow=c(length(known.markers),1))
for(i in seq(length(e1))){
	boxplot(e1[[i]], main=known.markers[i], outline=F, bty="n", col=cl.cols[colnames(cl.ave.expr)], border=cl.cols[colnames(cl.ave.expr)], lwd=1, names=cluster.order, frame=F, xaxt="n", yaxt="n")
}
dev.off()

# feature plots 
plot.coor=epi@dr$tsne@cell.embeddings
gene.cols <- c("#ffffcc", "#a1dab4", "#41b6c4", "#2c7fb8", "#253494")
g2 <- intersect(c("SFTPC", "AGER", "TP63", "IL33", "KRT5", "KRT15", "MUC5B", "MUC5AC", "SPDEF", "SCGB1A1", "SCGB3A2", "LTF", "FOXJ1", "MKI67", "CDK1",  "DCLK1", "NKX2-1", "EPCAM", "FOXI1"), rownames(fetal@data))
for(x in g2){
	file.name <- paste("Plot_fetal_epi_", x, ".png")
	png(file.name, height=20, width=20, unit="cm", res=500)
	plotFeature2(plot.coor, values=epi@data[x,], xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", cex=1,nCols=gene.cols)
	dev.off()
}

# heatmap for basal cell and bud tip progenitor cluster 
cl.ave.quantile <- t(apply(cl.ave.expr, 1, function(vec){
	(vec-min(vec))/(max(vec)-min(vec))
}))
markers <- read.table("../for_manuscript/Table_fetal_cluster_markers_for_fetal_vitro_analysis.txt", sep="\t", stringsAsFactors=F, head=T)
mat <- markers[which(markers$cluster==7),]
basal.cm <- setdiff(c(mat$gene_name[order(mat$avg_logFC, decreasing=T)[1:21]], "EGFR", "F3"), "IGF2")

mat <- markers[which(markers$cluster==5),]
budtip.cm <- setdiff(mat$gene_name[order(mat$avg_logFC, decreasing=T)[1:22]], c("CPM", "COL9A3"))

colors <- colorRampPalette(gene.cols, space="rgb")(50)
file.name <- "Plot_heatmap_basal_cluster_top20_markers.pdf"
file.name <- "Plot_heatmap_budTipProgenitor_cluster_top20_markers.pdf"
pdf(file.name, height=8)
heatmap.2(cl.ave.quantile[basal.cm,],trace="none",density.info="none",scale="none",dendrogram="none", col=colors,Rowv=FALSE, Colv=FALSE, labCol=names(cl.cols)[colnames(cl.ave.expr)], ColSideColors=cl.cols[colnames(cl.ave.expr)], srtCol=30, margins = c(8, 5))
dev.off()

# SPRING color coded by clusters
fetal.spring.coor <- read.table("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/SPRING_cm/epi_noCrrection_10PC_cm_k50/coordinates.txt",sep=",",row.names=1)
plot.coor=cbind(fetal.spring.coor[,2], -fetal.spring.coor[,1])
rownames(plot.coor) <- rownames(epi@meta.data)
saveRDS(plot.coor, file="Res_epi_SPRING_coordinate.rds")
hvg.expr <- read.csv("../SPRING_cm/Table_data.csv", head=F, row.names=1)
disMat <- 1-cor(hvg.expr)
# get the kNN network
k <- 15
idx1 <- apply(disMat, 2, function(vec){
	order(vec)[2:(k+1)]
})
knn.idx <- cbind(as.vector(idx1), rep(seq(ncol(disMat)), each=k))
group.res=paste0("Cluster", epi@meta.data$merge.cl.idx)
group.cols <- cl.cols
png("Plot_SPRING_epi_noCorrection_merged_cluster_info_with_edges.png", height=20, width=20, unit="cm", res=500) 
plotFeature2(coor=plot.coor, values=group.res, knn.pairs=knn.idx, main="", cex=1, cex.main=3, xaxt="n", yaxt="n", bty="n", xlab="", ylab="",gCols=cl.cols)
dev.off()

# SPRING color coded by samples
g.cols <- c("#a4d371", "#3c7ab6", "#4daa99", "#a64a97", "#c86879", "#1f7539", "#322a84", "#852655")
names(g.cols) <- paste0("S", seq(8))
sample.names <- c("11.5 week distal", "11.5 week airway", "15 week distal", "15 week airway", "15 week trachea", "18 week distal", "18 week airway", "21 week trachea")
group.res=epi@meta.data$orig.ident
png("Plot_SPRING_epi_noCorrection_sample_info_with_edges.png", height=20, width=20, unit="cm", res=500)
plotFeature2(coor=plot.coor, values=group.res, main="", cex=1, cex.main=3, xaxt="n", yaxt="n", bty="n", xlab="", ylab="",gCols=g.cols, knn.pairs=knn.idx)
dev.off()

# feature plots on SPRING coordinate
g2 <- c("TP63", "SFTPC", "PDPN", "SCGB3A2", "SFTPB", "FOXJ1")
g2 <- "ATOH1"
gene.cols <- c("#ffffcc", "#a1dab4", "#41b6c4", "#2c7fb8", "#253494")
for(x in g2){
	file.name <- paste0("Plot_SPRING_fetal_epi_", x, ".png")
	png(file.name, height=20, width=20, unit="cm", res=500)
	plotFeature2(plot.coor, values=epi@data[x,], knn.pairs=knn.idx, xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", cex=1,nCols=gene.cols)
	dev.off()
}


# give cell type name for each cluster
cell.type.vec <- rep(NA, nrow(epi@meta.data))
for(x in names(cell.types)){
	ct.name <- cell.types[x]
	cl.idx <- sub("Cluster", "", x)
	cell.idx <- which(epi@meta.data$merge.cl.idx==cl.idx)
	cell.type.vec[cell.idx] <- ct.name
}
epi@meta.data$Cell_type <- cell.type.vec
saveRDS(epi, file="Res_epi2_noCorrection_10PC_2reso.rds")

# in vitro analysis
d0 <- readRDS("../Res_d0.rds")
d3 <- readRDS("../Res_d3.rds")
d21 <- readRDS("../Res_d21.rds")
combined.vitro <- merge(x=d0, y=list("d3"=d3, "d21"=d21))
combined.vitro <- FindVariableFeatures(object = combined.vitro, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
combined.vitro <- ScaleData(object = combined.vitro, vars.to.regress = c("nCount_RNA", "percent.mito"), verbose = T)
combined.vitro <- RunPCA(object = combined.vitro, features = VariableFeatures(combined.vitro), verbose = T, npcs = 20, ndims.print = 2, nfeatures.print = 2)
usefulPCs <- 1:20
combined.vitro <- FindNeighbors(object = combined.vitro, dims = usefulPCs, force.recalc = T, k.param = 15)
combined.vitro <- FindClusters(object = combined.vitro, resolution = 0.8)
combined.vitro <- RunTSNE(object = combined.vitro, dims = usefulPCs)
saveRDS(combined.vitro, file="Res_inVitro_combined.rds")

cm <- findAllMarkers(seu.obj=combined.vitro, seu.version=3, selected.column="RNA_snn_res.0.8")
#plot fetal tSNE showing sample information
tsne.coor <- Embeddings(combined.vitro, reduction="tsne")
group.res <- paste0("Cluster", combined.vitro@meta.data$RNA_snn_res.0.8)
cluster.cols <- c("#C0392B", "#E74C3C", "#9B59B6", "#8E44AD", "#2980B9", "#3498DB", "#1ABC9C", "#16A085", "#27AE60", "#2ECC71", "#F1C40F", "#F39C12", "#E67E22", "#D35400")
names(cluster.cols) <- sort(unique(group.res))
png("Plot_tsne_invitro_combined_cluster_identity.png", height=20, width=20, unit="cm", res=500)
plotFeature2(tsne.coor, values=group.res, xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", cex=1,gCols=cluster.cols)
for(i in unique(group.res)){
        idx <- which(group.res==i)
        coor <- apply(tsne.coor[idx,],2,median)
        text(coor[1], coor[2], labels=sub("Cluster", "", i), cex=1.5)
}
dev.off()

group.res <- combined.vitro@meta.data$Sample
sample.cols <- c("#BB8FCE", "#85C1E9", "#82E0AA")
names(sample.cols) <- c("D0", "D3", "D21")
png("Plot_tsne_invitro_samples.png", height=20, width=20, unit="cm", res=500)
plotFeature2(tsne.coor, values=group.res, xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", cex=1,gCols=sample.cols)
legend("topleft", legend=names(sample.cols), bty="n", text.col=sample.cols)
dev.off()

gene.cols <- c("#ffffcc", "#a1dab4", "#41b6c4", "#2c7fb8", "#253494")
g2 <- c("EPCAM", "NKX2-1", "SFTPC", "SOX9", "NPC2", "TESC", "CA2","ETV5", "ID2", "AGER", "TP63", "IL33", "KRT5", "KRT15", "MUC5B", "MUC5AC", "SPDEF", "SCGB1A1", "SCGB3A2", "LTF", "FOXJ1", "MKI67", "CDK1", "EGFR", "F3")
g2 <- grep("MT-", rownames(combined.vitro), value=T)
column.num <- 5
row.num <- ceiling(length(g2)/column.num)
#png("Plot_tsne_invitro_known_cell_type_marker_expression.png", height=500*row.num, width=500*column.num)
png("Plot_tsne_invitro_mt_gene_expression.png", height=500*row.num, width=500*column.num)
par(mfrow=c(row.num, column.num))
for(x in g2){
	plotFeature2(tsne.coor, values=combined.vitro@assays$RNA@data[x,], xaxt="n", yaxt="n", bty="n",main=x, xlab="", ylab="", cex=1, cex.main=3, nCols=gene.cols)
}
dev.off()

top.cm <- cm %>% group_by(cluster) %>% top_n(n=5, wt=avg_logFC)
selected.markers <- unique(top.cm$gene_name)
column.num <- 5
row.num <- ceiling(length(selected.markers)/column.num)
png("Plot_tsne_invitro_denovo_cluster_marker_expression.png", height=500*row.num, width=500*column.num)
par(mfrow=c(row.num, column.num))
for(x in selected.markers){
	plotFeature2(tsne.coor, values=combined.vitro@assays$RNA@data[x,], xaxt="n", yaxt="n", bty="n",main=x, xlab="", ylab="", cex=1, cex.main=3, nCols=gene.cols)
}
dev.off()

cl6 <- cm$gene_name[which(cm$cluster==6)]
cl7 <- cm$gene_name[which(cm$cluster==7)]
cl6.specific <- setdiff(cl6, cl7)
cl7.specific <- setdiff(cl7, cl6)
g2 <- cl7.specific[1:20]
column.num <- 5
row.num <- ceiling(length(g2)/column.num)
png("Plot_tsne_invitro_cluster7_marker_expression.png", height=500*row.num, width=500*column.num)
par(mfrow=c(row.num, column.num))
for(x in g2){
	plotFeature2(tsne.coor, values=combined.vitro@assays$RNA@data[x,], xaxt="n", yaxt="n", bty="n",main=x, xlab="", ylab="", cex=1, cex.main=3, nCols=gene.cols)
}
dev.off()

png("Plot_tsne_invitro_basic_statistics.png", height=500*2, width=500*2)
par(mfrow=c(2,2))
plotFeature2(tsne.coor, values=combined.vitro@meta.data$percent.mito, xaxt="n", yaxt="n", bty="n",main="MT (%)", xlab="", ylab="", cex=1, cex.main=3, nCols=gene.cols)
plotFeature2(tsne.coor, values=combined.vitro@meta.data$nCount_RNA, xaxt="n", yaxt="n", bty="n",main="Count number", xlab="", ylab="", cex=1, cex.main=3, nCols=gene.cols)
plotFeature2(tsne.coor, values=combined.vitro@meta.data$nFeature_RNA, xaxt="n", yaxt="n", bty="n",main="Gene number", xlab="", ylab="", cex=1, cex.main=3, nCols=gene.cols)
dev.off()

v1 = combined.vitro@meta.data$percent.mito[which(combined.vitro@meta.data$Sample=="D21")]
v2 = combined.vitro@meta.data$percent.mito[which(combined.vitro@meta.data$Sample=="D3")]
v3 = combined.vitro@meta.data$percent.mito[which(combined.vitro@meta.data$Sample=="D0")]
pdf("Plot_density_mt_percentage.pdf")
plot(density(v3, from=0, to=0.1, bw=0.003), lwd=2, col=sample.cols["D0"], ylim=c(0, 50), bty="n", main="")
lines(density(v2, from=0, to=0.1, bw=0.003), lwd=2, col=sample.cols["D3"])
lines(density(v1, from=0, to=0.1, bw=0.003), lwd=2, col=sample.cols["D21"])
legend("topright", legend=names(sample.cols), bty="n", text.col=sample.cols)
dev.off()

vitroCell2vivoCluster.top.cm <- list()
for(x in c("d0", "d3", "d21")){
	rds.file <- paste0("../in_vitro/Res_cell2Cluster_PCC_", x, ".rds")
	res <- readRDS(rds.file)
	vitroCell2vivoCluster.top.cm[[x]] <- res
}
saveRDS(vitroCell2vivoCluster.top.cm, file="Res_vitroCell2vivoCluster_pcc.rds")

# plot the heatmap showing distribution of correlation to different in vivo clusters for each day in vitro cells
selected.organoid <- "d21"
sample.cor.mat <- t(vitroCell2vivoCluster.top.cm[[selected.organoid]]$cell2cluster.cor.mat)
bestVivoCluster2VItroCell.vec <- vitroCell2vivoCluster.top.cm[[selected.organoid]]$best.cor.cl

# color d21 cells according to best correlated cell types
d21 <- readRDS("../Res_d21.rds")
group.res=bestVivoCluster2VItroCell.vec
plot.coor=Embeddings(d21, reduction="tsne")
file.name <- "Plot_tSNE_d21_bestCorrelatedFetalClusters.png"
png(file.name, height=20, width=20, res=500, unit="cm")
plotFeature2(coor=plot.coor, values=group.res,xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", cex=1, cex.main=3, gCols=cl.cols)
dev.off()

# color d3 cells according to best correlated cell types
sample.cor.mat <- t(vitroCell2vivoCluster.top.cm[["d3"]]$cell2cluster.cor.mat)
bestVivoCluster2VItroCell.vec <- vitroCell2vivoCluster.top.cm[["d3"]]$best.cor.cl
group.res=bestVivoCluster2VItroCell.vec
plot.coor=Embeddings(d3, reduction="tsne")
file.name <- "Plot_tSNE_d3_bestCorrelatedFetalClusters.png"
png(file.name, height=20, width=20, res=500, unit="cm")
plotFeature2(coor=plot.coor, values=group.res,xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", cex=1, cex.main=3, gCols=cl.cols)
dev.off()

# get the PCC range for the best correlated clusters
best.cl.pcc.range <- sapply(unique(bestVivoCluster2VItroCell.vec), function(i){
	cell.idx <- which(bestVivoCluster2VItroCell.vec==i)
	vec <- round(sample.cor.mat[cell.idx, i],1)
	return(c(min(vec), max(vec)))
})

selected.cell.idx <- sample(seq(nrow(sample.cor.mat)), 2000) 
selected.cor.mat <- sample.cor.mat[selected.cell.idx,]
bestVivoCluster2VItroCell.top.cm <- bestVivoCluster2VItroCell.vec[selected.cell.idx]

values=bestVivoCluster2VItroCell.top.cm
cl.cols <- c("#197636", "#7BCAA4", "#A4D371", "#EED951", "#DCCB7C", "#FDD884", "#4DA7B0", "#CA6778","#9E0142","#5E4FA2", "#A84798", "#3C7AB6")
names(cl.cols) <- paste0("Cluster", c(2,1,11,7,10,12,3,6,5,4,8,9))
vitro.cell.idx <- unlist(lapply(paste0("Cluster",seq(cl.cols)), function(x){
	which(values==x)
}))
colors <- colorRampPalette(c("#2166ac","#67a9cf","#d1e5f0","#f7f7f7","#fddbc7","#ef8a62","#b2182b"), space="rgb")(50)
cell.col <- cl.cols[values[vitro.cell.idx]]

library(gplots)
file.name <- paste0("Plot_heatmap_", selected.organoid, "_cell_vivoCell2VitroCluster_topCM.png")
png(file.name, height=20, width=20, res=500, unit="cm")
par(mar=c(6,6,5,3))
heatmap.2(selected.cor.mat[vitro.cell.idx, ],trace="none",density.info="none",scale="none",dendrogram="none", col=colors,Rowv=FALSE, Colv=FALSE,labRow=NA, RowSideColors=cell.col, ColSideColors=cl.cols[colnames(selected.cor.mat)])
dev.off()

d3 <- readRDS("../Res_d3.rds")
cm <- findAllMarkers(seu.obj=d3, seu.version=3, selected.column="RNA_snn_res.0.8", do.par=F, cluster.to.test=setdiff(unique(d3@meta.data$RNA_snn_res.0.8), 2))
mat <- FindMarkers(seu.obj, ident.1=2, only.pos=T, min.pct=0.25, logfc.threshold=0.1)
dat <- data.frame(mat, "cluster"=rep(2, nrow(mat)), "gene_name"=rownames(mat), stringsAsFactors=F)
rownames(dat) <- paste(rownames(dat), 2, sep="_")
cm <- rbind(cm, dat)
saveRDS(cm, file="Res_d3_cluster_markers.rds")
top.cm <- cm %>% group_by(cluster) %>% top_n(n=5, wt=avg_logFC)
selected.markers <- unique(top.cm$gene_name)
column.num <- 5
row.num <- ceiling(length(selected.markers)/column.num)
tsne.coor=Embeddings(d3, reduction="tsne")
png("Plot_tsne_d3_denovo_cluster_marker_expression.png", height=500*row.num, width=500*column.num)
par(mfrow=c(row.num, column.num))
for(x in selected.markers){
	plotFeature2(tsne.coor, values=d3@assays$RNA@data[x,], xaxt="n", yaxt="n", bty="n",main=x, xlab="", ylab="", cex=1, cex.main=3, nCols=gene.cols)
}
dev.off()

top.cm <- cm %>% group_by(cluster) %>% top_n(n=50, wt=avg_logFC)
selected.genes <- unique(top.cm$gene_name)
selected.genes <- unique(cm$gene_name)
ref <- d3@assays$RNA@data["TP63",]
top.cm.cor.basal <- sapply(selected.genes, function(x){
	que <- d3@assays$RNA@data[x,]
	cor(que, ref)
})

cor.by.cl <- lapply(unique(cm$cluster), function(i){
	g0 <- intersect(cm$gene_name[which(cm$cluster==i)], selected.genes)
	top.cm.cor.basal[g0]
})
pdf("Plot_boxplot_cm_cor_to_basal_distribution_by_cluster.pdf")
boxplot(cor.by.cl, names=unique(cm$cluster))
dev.off()

#plot tSNE showing sample information
tsne.coor <- Embeddings(d3, reduction="tsne")
group.res <- paste0("Cluster", d3@meta.data$RNA_snn_res.0.8)
cluster.cols <- c("#C0392B", "#9B59B6", "#2980B9", "#1ABC9C", "#27AE60", "#F1C40F", "#E67E22")
names(cluster.cols) <- sort(unique(group.res))
png("Plot_tsne_d3_cluster_identity_no_idx.png", height=20, width=20, unit="cm", res=500)
plotFeature2(tsne.coor, values=group.res, xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", cex=1,gCols=cluster.cols)
#for(i in unique(group.res)){
#        idx <- which(group.res==i)
#        coor <- apply(tsne.coor[idx,],2,median)
#        text(coor[1], coor[2], labels=sub("Cluster", "", i), cex=1.5)
#}
dev.off()

# feature plots
markers <- read.table("../for_manuscript/Table_fetal_cluster_markers_for_fetal_vitro_analysis.txt", sep="\t", stringsAsFactors=F, head=T)
basal.cm <- intersect(markers$gene_name[which(markers$cluster==7)], rownames(combined.vitro))
score <- apply(apply(combined.vitro@assays$RNA@data[basal.cm, ], 1, scale), 1, sum)
names(score) <- colnames(combined.vitro)
colorPal <- grDevices::colorRampPalette(c("dark green", "yellow", "red"))
cellColor <- adjustcolor(colorPal(30), alpha=.8)[as.numeric(cut(score, breaks=30, right=F,include.lowest=T))]
names(cellColor) <- names(score)
tsne.d0 <- Embeddings(d0, reduction="tsne")
tsne.d3 <- Embeddings(d3, reduction="tsne")
tsne.d21 <- Embeddings(d21, reduction="tsne")
png("Plot_tSNE_cl7_marker_sum_scaled_expr_each_organoid-2.png", height=20, width=60, unit="cm", res=500)
par(mfrow=c(1,3))
plot(tsne.d0, col=cellColor[rownames(tsne.d0)], pch=16,xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", cex=2)
plot(tsne.d3, col=cellColor[rownames(tsne.d3)], pch=16,xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", cex=2)
plot(tsne.d21, col=cellColor[rownames(tsne.d21)], pch=16,xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", cex=2)
dev.off()

# check fetal cluster 9 expression
d3.cm <- readRDS("Res_d3_cluster_markers.rds")
g0 <- intersect(d3.cm$gene_name, markers$gene_name[which(markers$cluster==9)])
png("Plot_feature_plot_d3_cluster_marker_fetal_cluster_marker_overlap.png", height=7*500, width=8*500)
par(mfrow=c(7, 8))
for(x in g0){
	plotFeature2(tsne.d3, values=d3@assays$RNA@data[x,], xaxt="n", yaxt="n", bty="n",main=x, xlab="", ylab="", cex=1,nCols=gene.cols, cex.main=5)
}
dev.off()

# GO enrichment for d3 cluster markers
dir.create("d3_GO/")
link <- read.table("~/Work/Annotation/Table_10X_gene_symbols_and_ensemblID_hg19.tsv", stringsAsFactors=F, sep="\t")
for(i in unique(d3.cm$cluster)){
	g0 <- d3.cm$gene_name[which(d3.cm$cluster==i)]
	g0.id <- unique(link[which(link[,2]%in%g0),1])
	file.name <- paste0("d3_GO/List_d3_cluster",i,"_markers.txt")
	write.table(g0.id, file=file.name, row.names=F, col.names=F, quote=F)
}
g0 <- rownames(d3)
g0.id <- unique(link[which(link[,2]%in%g0),1])
file.name <- "d3_GO/List_d3_expressed_genes.txt"
write.table(g0.id, file=file.name, row.names=F, col.names=F, quote=F)

# read in the GO enrichment analysis results
go.res <- c()
for(i in sort(unique(d3.cm$cluster))){
	file.name <- paste0("d3_GO/GO_d3_cl",i,".txt")
	dat <- read.table(file.name, sep="\t", stringsAsFactors=F, fill=T)
	idx <- which(dat[,12]<0.05 & dat[,12] != "")
	if(length(idx)>0){
		for(x in idx){
			term <- dat[x,2]
			#genes <- strsplit(dat[x,6], ", ")
			genes <- dat[x,6]
			fold.enrichment <- dat[x,10]
			pval <- dat[x,5]
			prop <- dat[x,4]
			bonf <- dat[x,11]
			bh <- dat[x,12]
			tmp <- c(paste0("Cluster",i), term, genes, fold.enrichment, pval, bonf, bh)
			go.res <- rbind(go.res, tmp)
		}
	}
}
go.res.2 <- data.frame("Cluster"=go.res[,1], "Term"=go.res[,2], "Gene"=go.res[,3], "Fold_enrichment"=as.numeric(go.res[,4]), "P"=as.numeric(go.res[,5]), "Bonferroni"=as.numeric(go.res[,6]), "BH"=as.numeric(go.res[,7]), stringsAsFactors=F)

gene.symbol <- sapply(go.res.2$Gene, function(vec){
	id <- strsplit(vec, ", ")[[1]]
	symbol <- unique(link[link[,1]%in%id, 2])
	symbol <- paste(symbol, collapse=", ")
	return(symbol)
})
go.res.2[["Gene_symbol"]] <- gene.symbol
saveRDS(go.res.2, file="Res_d3_cm_with_enriched_GO_terms.rds")
go.enriched.cm <- unique(unlist(strsplit(go.res.2$Gene_symbol, ", ")))
expr.mat <- as.matrix(d3@assays$RNA@data[go.enriched.cm, ])
tp63.pattern <- d3@assays$RNA@data["TP63",]
cm.tp63.cor <- cor(t(expr.mat), tp63.pattern)
cor.vec <- sapply(go.res.2$Gene_symbol, function(vec){
	g0 <- strsplit(vec, ", ")[[1]]
	pcc <- cm.tp63.cor[g0,]
	pcc.median <- median(pcc) 
	return(pcc.median)
})
go.res.2[["PCC_median_to_TP63_day3"]] <- cor.vec
saveRDS(go.res.2, file="Res_d3_cm_with_enriched_GO_terms.rds")

d0.total.mean <- sapply(go.enriched.cm, function(x){
	if(sum(rownames(d0)%in%x)>0){
		return(mean(d0@assays$RNA@data[x,]))
	}else{
		return(0)
	}
})
d3.total.mean <- apply(d3@assays$RNA@data[go.enriched.cm,], 1, mean) 
d3.diff <- d3.total.mean-d0.total.mean
logfc.vec <- sapply(go.res.2$Gene_symbol, function(vec){
	g0 <- strsplit(vec, ", ")[[1]]
	diff <- d3.total.mean[g0]
	diff.median <- median(diff) 
	return(diff.median)
})
go.res.2[["day3_logfc"]] <- logfc.vec
saveRDS(go.res.2, file="Res_d3_cm_with_enriched_GO_terms.rds")

go.res.used <- go.res.2[grep("GO:", go.res.2$Term),]

names(cluster.cols) <- sapply(names(cluster.cols), function(x){
	paste(strsplit(x, "_")[[1]], collapse="")
})
pdf("Plot_d3_cm_GO_enrichment_fold_enrichment.pdf")
barplot(rev(log(go.res.used$Fold_enrichment)), col=rev(cluster.cols[go.res.used$Cluster]), border=NA, names=rev(terms), horiz=TRUE, las=2)
dev.off()

pdf("Plot_d3_cm_enriched_GO_pcc2TP63_day3.pdf")
barplot(rev(go.res.used$PCC_median_to_TP63_day3), col=rev(cluster.cols[go.res.used$Cluster]), border=NA, names=rev(terms), horiz=TRUE, las=2)
dev.off()

pdf("Plot_d3_cm_enriched_GO_logfc_day3vs0.pdf")
barplot(rev(go.res.used$day3_logfc), col=rev(cluster.cols[go.res.used$Cluster]), border=NA, names=rev(terms), horiz=TRUE, las=2)
dev.off()

pos.go <- go.res.used[go.res.used$PCC_median_to_TP63_day3>0,c("Term", "PCC_median_to_TP63_day3")]
neg.go <- go.res.used[go.res.used$PCC_median_to_TP63_day3<0,c("Term", "PCC_median_to_TP63_day3")]

mat <- go.res.used[,c("PCC_median_to_TP63_day3", "day3_logfc")]
terms <- sapply(go.res.used$Term, function(x){
	strsplit(x, "~")[[1]][2]
})
group.res <- paste0("Cluster", d3@meta.data$RNA_snn_res.0.8)
cluster.cols <- c("#C0392B", "#9B59B6", "#2980B9", "#1ABC9C", "#27AE60", "#F1C40F", "#E67E22")
names(cluster.cols) <- sort(unique(group.res))
cols <- cluster.cols[go.res.used$Cluster]
size <- go.res.used$Fold_enrichment^0.25 
for(x in unique(go.res.used$Cluster)){
	idx <- which(go.res.used$Cluster==x)
	file.name <- paste0("Plot_scatter_cm_enriched_GO_pcc_vs_fc_", x, ".pdf")
	pdf(file.name)
	plot(mat[idx,], pch=16, bty="n", cex=size[idx], col=cols[idx])
	abline(v=0)
	text(mat[idx,], label=terms[idx])
	dev.off()
}

group.res <- paste("Cluster", d3@meta.data$RNA_snn_res.0.8, sep="_")
cluster.cols <- c("#C0392B", "#9B59B6", "#2980B9", "#1ABC9C", "#27AE60", "#F1C40F", "#E67E22")
names(cluster.cols) <- sort(unique(group.res))

cl.ave.expr <- getAveExpr(meta.data=d3@meta.data, feature.to.calc="RNA_snn_res.0.8", expr=d3@assays$RNA@data, genes=NULL, core.num=20, colname.prefix="Cluster")
saveRDS(cl.ave.expr, file="Res_vitro_d3_cluster_ave_expr.rds")
cl.ave.expr <- cl.ave.expr[,c(2,6,5,4,1,7,3)]
# heatmap for expression patterns for selected top markers across clusters
top.d3.cm.mat <- d3.cm %>% group_by(cluster) %>% top_n(n=50, wt=avg_logFC)
top.d3.cm <- unique(top.d3.cm.mat$gene_name)
colors <- colorRampPalette(c("#f7f7f7", "#d9d9d9", "#252525"), space="rgb")(50)
file.name <- "Plot_vitro_d3_denovo_cluster_markers_by_cluster.pdf"
pdf(file.name, height=8)
heatmap.2(cl.ave.expr[top.d3.cm,],trace="none",density.info="none",scale="row",dendrogram="none", col=colors,Rowv=FALSE, Colv=FALSE, labRow="", ColSideColors=cluster.cols[colnames(cl.ave.expr)], srtCol=30, margins = c(8, 5))
dev.off()

# boxplot for expression patterns for selected top markers across clusters
top.d3.cm.mat <- d3.cm %>% group_by(cluster) %>% top_n(n=5, wt=avg_logFC)
top.d3.cm <- c(unique(top.d3.cm.mat$gene_name), "FHL2", "TP63")
names(cluster.cols) <- sub("Cluster", "Cluster_", names(cluster.cols))
cluster.order <- sub("Cluster_", "", colnames(cl.ave.expr))
e1 <- lapply(top.d3.cm, function(x){
	lapply(cluster.order, function(i){
		d3@assays$RNA@data[x,which(d3@meta.data$RNA_snn_res.0.8==i)]
	})	
})
pdf("Plot_boxplot_denovo_markers_d3-tp63.pdf", height=100)
par(mfrow=c(length(top.d3.cm),1))
for(i in seq(length(e1))){
	boxplot(e1[[i]], main=top.d3.cm[i], outline=F, bty="n", col=cluster.cols[colnames(cl.ave.expr)], border=cluster.cols[colnames(cl.ave.expr)], lwd=1, names=cluster.order, frame=F, xaxt="n", yaxt="n")
}
dev.off()

marker.genes <- c("EPCAM", "SFTPC", "SOX9", "TP63", "SCGB3A2", "AGER", "PDPN", "FOXJ1", "FHL2", "CHGA", "SFTPB", "LTF", "MUC5B", "SCGB1A1", "KRT5", "KRT15")
for(x in marker.genes){
	score=combined.vitro@assays$RNA@data[x,]
	colorPal <- grDevices::colorRampPalette(gene.cols)
	cellColor <- adjustcolor(colorPal(30), alpha=.8)[as.numeric(cut(score, breaks=30, right=F,include.lowest=T))]
	cellColor[which(score==0)] <- "#bdbdbd30"
	names(cellColor) <- names(score)
	file.name <- paste0("Plot_tSNE_vitro_", x, "_expr_each_organoid-2.png")
	png(file.name, height=20, width=60, unit="cm", res=500)
	par(mfrow=c(1,3))
	plot(tsne.d0, col=cellColor[rownames(tsne.d0)], pch=16,xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", cex=2)
	plot(tsne.d3, col=cellColor[rownames(tsne.d3)], pch=16,xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", cex=2)
	plot(tsne.d21, col=cellColor[rownames(tsne.d21)], pch=16,xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", cex=2)
	dev.off()
}

de.res.list <- readRDS("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/basal_include_cl9/Res_basal_lineage_pairwise_DE_comparison_and_pathway_enrichment.rds")
pathway.res <- sapply(seq(length(de.res.list)), function(i){
	res <- de.res.list[[i]]$pathway
	or.vec <- res[,2]
	pval.vec <- res[,4]
	or.vec[which(pval.vec>0.05)] <- 0
	return(or.vec)
})
colnames(pathway.res) <- names(de.res.list)
pathway.res <- pathway.res[order(pathway.res[,"Cl7_Cl5"]),]
rownames(pathway.res) <- sub("HALLMARK_", "", rownames(pathway.res))
rownames(pathway.res) <- gsub("_", " ", rownames(pathway.res))
cols <- colorRampPalette(c("#d9d9d9", "#252525"))(50)
col.pos <- seq(from=0, to=1, length=ncol(pathway.res))
row.pos <- seq(from=0, to=1, length=nrow(pathway.res))
pdf("Plot_pairwise_pathway_enrichment.pdf")
par(mar=c(6,10,2,2))
image(t(pathway.res), col=cols, xaxt="n", yaxt="n", bty="n")
mtext(text=colnames(pathway.res), side=1, at=col.pos, line=0.5)
mtext(text=rownames(pathway.res), side=2, las=2, at=row.pos, cex=0.4, line=0.5)
dev.off()


