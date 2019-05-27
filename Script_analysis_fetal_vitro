library(Seurat)
library(dplyr)
source("~/Work/commonScript/Script_functions.R")

# color scheme
#cl.cols <- c("#197636", "#7BCAA4", "#A4D371", "#EED951", "#DCCB7C", "#FDD884", "#4DA7B0", "#CA6778","#9E0142","#5E4FA2", "#A84798", "#3C7AB6","#FDB164")
cl.cols <- c("#197636", "#7BCAA4", "#A4D371", "#EED951", "#DCCB7C", "#FDD884", "#4DA7B0", "#CA6778","#9E0142","#5E4FA2", "#A84798", "#3C7AB6")
names(cl.cols) <- paste0("Cluster", c(2,1,11,7,10,12,3,6,5,4,8,9))
cell.types <- c("Multiciliated", "Multiciliated precursor", "Intermediate", "Basal cell", "Goblet-like secretory", "Club-like secretory", "Neuroendocrine", "Bud tip adjacent", "Bud tip progenitor", "Hub progenitor", "Submucosal gland",  "Undefined")
names(cell.types) <- paste0("Cluster", c(2,1,11,7,10,12,3,6,5,4,8,9))

# load all fetal data
########## Week 11.5 ##########
distal.11.5 <- prepareSeuratObject(rawData = readRDS("/home/qianhui_yu/Work/Lung/Data/Dat_HT234_11.5w_distal_count_matrix.rds"),
                               namePrefix = "S1",
                               age = 11.5, 
                               tissue = "Distal", 
                               sex = NA,
                               sampleName = "W11.5_distal",
                               condition = NA,
                               isolation = NA,
							   min.gene.cutoff=1500,
							   gene.num.low=1500
							   )

airway.11.5 <- prepareSeuratObject(rawData = readRDS("/home/qianhui_yu/Work/Lung/Data/Dat_HT234_11.5w_airway_count_matrix.rds"),
                               namePrefix = "S2",
                               age = 11.5, 
                               tissue = "Airway", 
                               sex = NA,
                               sampleName = "W11.5_airway",
                               condition = NA,
                               isolation = NA,
							   min.gene.cutoff=1500,
							   gene.num.low=1500
							   )

########## Week 15 #############
distal.15 <- prepareSeuratObject(rawData = Read10X(data.dir = "/mnt/SingleCellGenomics/scg_projects/Spence_Michigan_collab/Lung/Run_2321_Miller_HT187lung_HT145DSA-DSI_CzerwinskiHT188AdultDuo/Sample_HT-187-Distal-lung/outs/HT-187-Distal-lung-filtered_gene_bc_matrices/hg19"),
                               namePrefix = "S3",
                               age = 15, 
                               tissue = "Distal", 
                               sex = NA,
                               sampleName = "W15_distal",
                               condition = NA,
                               isolation = NA,
							   min.gene.cutoff=1500,
							   gene.num.low=1500
							   )
                              
airway.15 <- prepareSeuratObject(rawData = Read10X(data.dir = "/mnt/SingleCellGenomics/scg_projects/Spence_Michigan_collab/Lung/Run_2321_Miller_HT187lung_HT145DSA-DSI_CzerwinskiHT188AdultDuo/Sample_HT-187-Small-airway/outs/HT-187-Small-airway-filtered_gene_bc_matrices/hg19"),
                               namePrefix = "S4",
                               age = 15, 
                               tissue = "Airway", 
                               sex = NA,
                               sampleName = "W15_airway",
                               condition = NA,
                               isolation = NA,
							   min.gene.cutoff=1500,
							   gene.num.low=1500
                               )							  
							 						 
trachea.15 <- prepareSeuratObject(rawData = Read10X(data.dir = "/mnt/SingleCellGenomics/scg_projects/Spence_Michigan_collab/Lung/Run_2321_Miller_HT187lung_HT145DSA-DSI_CzerwinskiHT188AdultDuo/Sample_HT-187-Tracheal-epi/outs/HT-187-Tracheal-epi-filtered_gene_bc_matrices/hg19"),
                               namePrefix = "S5",
                               age = 15, 
                               tissue = "Trachea", 
                               sex = NA,
                               sampleName = "W15_trachea",
                               condition = NA,
                               isolation = NA,
							   min.gene.cutoff=1500,
							   gene.num.low=1500
                               )	

########## Week 18 #############							   
distal.18 <- prepareSeuratObject(rawData = Read10X(data.dir = "/mnt/SingleCellGenomics/scg_projects/Spence_Michigan_collab/Lung/Run_2288_Miller_D125_fetalLung_LungBudCulture/Sample_HT-182-d125-lung-Distal/outs/HT-182-d125-lung-Distal-filtered_gene_bc_matrices/hg19"),
                               namePrefix = "S6",
                               age = 18, 
                               tissue = "Distal", 
                               sex = NA,
                               sampleName = "W18_distal",
                               condition = NA,
                               isolation = NA,
							   min.gene.cutoff=1500,
							   gene.num.low=1500
                               )

airway.18 <- prepareSeuratObject(rawData = Read10X(data.dir = "/mnt/SingleCellGenomics/scg_projects/Spence_Michigan_collab/Lung/Run_2288_Miller_D125_fetalLung_LungBudCulture/Sample_HT-182-d125-lung-Prox/outs/HT-182-d125-lung-Prox-filtered_gene_bc_matrices/hg19"),
                               namePrefix = "S7",
                               age = 18, 
                               tissue = "Airway", 
                               sex = NA,
                               sampleName = "W18_airway",
                               condition = NA,
                               isolation = NA,
							   min.gene.cutoff=1500,
							   gene.num.low=1500
                               )

########## week 21 #############
trachea.21 <- prepareSeuratObject(rawData = Read10X(data.dir = "/mnt/SingleCellGenomics/scg_projects/Spence_Michigan_collab/Lung/Run_2330_Miller_HT189-Tracheal-Epithelium/Sample_HT-189-Tracheal-Epi/outs/HT189-TrachealEpi-filtered_gene_bc_matrices/hg19/"),
                               namePrefix = "S8",
                               age = 21, 
                               tissue = "Trachea", 
                               sex = NA,
                               sampleName = "W21_trachea",
                               condition = NA,
                               isolation = NA,
							   min.gene.cutoff=1500,
							   gene.num.low=1500
                               )

vivo.samples <- c("distal.11.5", "airway.11.5", "distal.15", "airway.15", "trachea.15", "distal.18", "airway.18", "trachea.21")

fetal <- MergeSeurat(object1 = get(vivo.samples[1]), object2 = get(vivo.samples[2]))
for(i in 3:length(vivo.samples)){
	fetal <- MergeSeurat(object1 = fetal, object2 = get(vivo.samples[i]))
}
saveRDS(fetal, file="Dat_inVivo_combined_data.rds")

fetal <- NormalizeData(object = fetal)
# Read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines("/home/qianhui_yu/Work/Annotation/cellCycle/regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
fetal <- CellCycleScoring(object = fetal, s.genes = s.genes, g2m.genes = g2m.genes, set.ident=FALSE)
fetal <- ScaleData(object=fetal, vars.to.regress=c("nUMI", "percent.mito", "S.Score", "G2M.Score"), do.par=T)
fetal <- FindVariableGenes(object=fetal, do.plot=FALSE)
saveRDS(fetal, file="Dat_inVivo_combined_data.rds")

# identify variable genes in each sample separately
samples <- as.character(unique(fetal@meta.data$orig.ident))
seu.by.sample <- list()
for(idx in samples){
	cat(paste(idx, "start\n"))
	work.dir <- paste0(idx, "/")
	dir.create(file.path(work.dir))
	
	seu.obj <- SubsetData(fetal, ident.use = idx, do.clean=T)
	seu.obj <- ScaleData(object=seu.obj, vars.to.regress=c("nUMI", "percent.mito", "S.Score", "G2M.Score"), do.par=T)
	seu.obj <- FindVariableGenes(object=seu.obj, do.plot=FALSE)

	
	seu.obj <- RunPCA(object=seu.obj, pc.genes=seu.obj@var.genes, do.print=FALSE)
	seu.obj <- FindClusters(object=seu.obj, reduction.type="pca", dims.use=1:20, resolution=0.8, print.output=0, save.SNN=TRUE)
	seu.obj <- RunTSNE(object=seu.obj, dims.use=1:20, do.fast=TRUE)

	saveRDS(seu.obj, file=paste0(work.dir,"Res_",idx,".rds"))
	
	seu.by.sample[[idx]] <- seu.obj
}

hvg.gene.idx <- matrix(F, nrow=nrow(fetal@data), ncol=length(samples))
rownames(hvg.gene.idx) <- rownames(fetal@data)
colnames(hvg.gene.idx) <- samples
for(j in seq(length(samples))){
	hvg.gene.idx[which(rownames(hvg.gene.idx)%in%seu.by.sample[[j]]@var.genes),j] <- T
}
# selected highly variable genes should be hvg in at least 2 samples
selected.hvg <- rownames(hvg.gene.idx)[apply(hvg.gene.idx, 1, sum)>1]
save(hvg.gene.idx, selected.hvg, file="Res_hvg_gene_by_sample.rdata")

fetal <- RunPCA(object=fetal, pc.genes=selected.hvg, do.print=FALSE)
fetal <- FindClusters(object=fetal, reduction.type="pca", dims.use=1:20, resolution=0.8, print.output=0, save.SNN=TRUE)
fetal <- RunTSNE(object=fetal, dims.use=1:20, do.fast=TRUE)

p2 <- TSNEPlot(object=fetal, pt.size=1, do.return=T, do.label=T,label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
p3 <- TSNEPlot(object=fetal, pt.size=1, group.by="orig.ident", do.return=T, label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
png("Plot_tSNE_res0.8.png",width=1000, height=500)
plot_grid(p2, p3, ncol=2, nrow=1)
dev.off()
saveRDS(fetal, file="Res_fetal_1500.rds")

# classify different lineages based on canonical lineage markers or novel markers well correlated with canonical markers
pan.epi.markers <- rev(c("EPCAM", "KRT8", "KRT18", "FXYD3", "PERP", "CDH1")) 
pan.meso.markers <- c("COL6A2","MFAP4","DCN","FHL1","COL1A2","COL3A1")
pan.imm.markers <- c("CD37","CORO1A","LCP1","PTPRC","CD53","LAPTM5")
endothelial.markers <- c("FLT1", "CDH5", "CLDN5", "EGFL7", "ESAM","IFI27")
neuron.mix.markers <- c("CHGA","ASCL1","NNAT","STMN2","GRP","MPZ")

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

pdf("DotPlot_fetal_lineages.pdf", width=10, height=7)
DotPlot(fetal, genes.plot=c(pan.imm.markers, endothelial.markers, pan.meso.markers, pan.epi.markers,neuron.mix.markers), x.lab.rot=TRUE, group.by="cl_idx_2", plot.legend=TRUE, cols.use=c("#d9d9d9", "#252525"))
dev.off()

# plot fetal tSNE showing sample information
g.cols <- c("#A4D371", "#3C7AB6", "#4DAA99", "#A64A97", "#C86879", "#1F7539", "#322A84", "#852655")
names(g.cols) <- paste0("S", seq(8))
sample.names <- c("11.5 week distal", "11.5 week airway", "15 week distal", "15 week airway", "15 week trachea", "18 week distal", "18 week airway", "21 week trachea")
tsne.coor <- fetal@dr$tsne@cell.embeddings
group.res <- as.character(fetal@meta.data$orig.ident)
pdf("Plot_tSNE_fetal_all_no_sample_identity.pdf", height=15, width=30)
par(mfrow=c(1,2))
plotFeature2(tsne.coor, values=group.res, xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", gCols=g.cols, cex=1.8)
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
pdf("Plot_tSNE_fetal_all_cluster_sequential_cols.pdf", height=15, width=15)
plotFeature2(tsne.coor, values=group.res, xaxt="n", yaxt="n", bty="n",main="", xlab="", ylab="", cex=1.8,gCols=cell.type.cols)
for(i in seq(27)){
	idx <- which(fetal@meta.data[,"cl_idx_2"]==i)
	coor <- apply(tsne.coor[idx,],2,median)
	text(coor[1], coor[2], labels=i, cex=2)
}
dev.off()

# extract epithelial cells, use highly variable genes identified in each cluster and perform subclustering
epi <- SubsetData(fetal, cells.use=fetal@cell.names[which(fetal@meta.data$lineage=="Epithelial")])
epi.by.sample <- list()
wd <- "/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/"
for(idx in samples){
	cat(paste(idx, "start\n"))
	work.dir <- paste0("Epi_", idx, "/")
	dir.create(file.path(work.dir))
	setwd(work.dir)
	epi.by.sample[[idx]] <- preprocessSubset(epi, field.to.select="orig.ident", subset.idx=idx, clean=TRUE, regressed.vars=c("percent.mito", "nUMI", "S.Score", "G2M.Score"),reso=0.8, pc.num=20, hvg.to.exclude=cc.genes, tsne.plot.name=NULL, cluster.only=TRUE, plot.known.markers=TRUE, marker.list.file="~/Work/Annotation/cellTypeMarker/Lung/Table_major_cell_type_markers_from_literature_search.txt", marker.plot.name=NULL, save.rds=TRUE, rds.file.name=NULL, do.return=TRUE)
	setwd(wd)
}

hvg.gene.idx <- matrix(F, nrow=nrow(fetal@data), ncol=length(samples))
rownames(hvg.gene.idx) <- rownames(fetal@data)
colnames(hvg.gene.idx) <- samples
for(j in seq(length(samples))){
	hvg <- setdiff(epi.by.sample[[j]]@var.genes, cc.genes)
	hvg.gene.idx[which(rownames(hvg.gene.idx)%in%hvg),j] <- T
}
# selected highly variable genes should be hvg in at least 2 samples
selected.hvg <- rownames(hvg.gene.idx)[apply(hvg.gene.idx, 1, sum)>1]
save(hvg.gene.idx, selected.hvg, file="Res_epi_hvg_gene_by_sample_ccGeneRemoved.rdata")

epi <- RunPCA(object=epi, pc.genes=selected.hvg, do.print=FALSE)
epi <- FindClusters(object=epi, reduction.type="pca", dims.use=1:20, resolution=0.8, print.output=0, save.SNN=TRUE, force.recalc=TRUE)
epi <- RunTSNE(object=epi, dims.use=1:20, do.fast=TRUE)

p2 <- TSNEPlot(object=epi, pt.size=1, do.return=T, do.label=T,label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
p3 <- TSNEPlot(object=epi, pt.size=1, group.by="orig.ident", do.return=T, label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
png("Plot_tSNE_epi_res0.8.png",width=1000, height=500)
plot_grid(p2, p3, ncol=2, nrow=1)
dev.off()

markers <- read.table("~/Work/Annotation/cellTypeMarker/Lung/Table_major_cell_type_markers_from_literature_search.txt",sep="\t",stringsAsFactors=F)
markers <- markers[markers[,2]%in%rownames(epi@data),]
g1 <- markers[,2]
marker.plot.name <- "Plot_marker_expression_on_epi_tSNE.png"
png(marker.plot.name, width=1600, height=4800)
FeaturePlot(object=epi, features.plot=g1, cols.use=c("gray", "blue"), no.legend=F)
dev.off()
saveRDS(epi, file="Res_epi_1.rds")

# Further filter EpCAM- clusters and use highly variable genes identified in each cluster and perform subclustering
epi <- SubsetData(epi, cells.use=epi@cell.names[which(!epi@meta.data$res.0.8%in%c(18,19,21))])
epi.by.sample <- list()
wd <- "/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/"
for(idx in samples){
	cat(paste(idx, "start\n"))
	work.dir <- paste0("Epi_", idx, "_2/")
	dir.create(file.path(work.dir))
	setwd(work.dir)
	epi.by.sample[[idx]] <- preprocessSubset(epi, field.to.select="orig.ident", subset.idx=idx, clean=TRUE, regressed.vars=c("percent.mito", "nUMI", "S.Score", "G2M.Score"),reso=0.8, pc.num=20, hvg.to.exclude=cc.genes, tsne.plot.name=NULL, cluster.only=TRUE, plot.known.markers=TRUE, marker.list.file="~/Work/Annotation/cellTypeMarker/Lung/Table_major_cell_type_markers_from_literature_search.txt", marker.plot.name=NULL, save.rds=TRUE, rds.file.name=NULL, do.return=TRUE)
	setwd(wd)
}

confound.genes <- c(cc.genes, c("HBB", "HBA1", "HBA2", "HBG1", "HBG2"))
hvg.gene.idx <- matrix(F, nrow=nrow(fetal@data), ncol=length(samples))
rownames(hvg.gene.idx) <- rownames(fetal@data)
colnames(hvg.gene.idx) <- samples
for(j in seq(length(samples))){
	hvg <- setdiff(epi.by.sample[[j]]@var.genes, confound.genes)
	hvg.gene.idx[which(rownames(hvg.gene.idx)%in%hvg),j] <- T
}
# selected highly variable genes should be hvg in at least 2 samples
selected.hvg <- rownames(hvg.gene.idx)[apply(hvg.gene.idx, 1, sum)>1]
save(hvg.gene.idx, selected.hvg, file="Res_epi2_hvg_gene_by_sample_ccGeneRemoved.rdata")

epi <- RunPCA(object=epi, pc.genes=selected.hvg, do.print=FALSE)
pc.num <- 10
reso <- 2
epi <- FindClusters(object=epi, reduction.type="pca", dims.use=1:pc.num, resolution=reso, print.output=0, save.SNN=TRUE, force.recalc=TRUE)
epi <- RunTSNE(object=epi, dims.use=1:pc.num, do.fast=TRUE)

p2 <- TSNEPlot(object=epi, pt.size=1, do.return=T, do.label=T,label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
p3 <- TSNEPlot(object=epi, pt.size=1, group.by="Sample", do.return=T, label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
plot.file <- paste0("Plot_tSNE_epi2_res",reso,"_",pc.num,"PC.png")
png(plot.file,width=1000, height=500)
plot_grid(p2, p3, ncol=2, nrow=1)
dev.off()
saveRDS(epi, file="Res_epi2_noCorrection_10PC_2reso.rds")

markers <- read.table("~/Work/Annotation/cellTypeMarker/Lung/Table_manuscript_used_cell_type_markers.txt",sep="\t",stringsAsFactors=F)
markers <- markers[markers[,2]%in%rownames(epi@data),]
g1 <- markers[,2]
marker.plot.name <- paste0("Plot_marker_expression_on_epi2_tSNE_res",reso,"_",pc.num,"PC.png")
png(marker.plot.name, width=1600, height=1600)
FeaturePlot(object=epi, features.plot=g1, cols.use=c("gray", "blue"), no.legend=F)
dev.off()

# for data without integration, 10PC
## find cluster markers
### merge similar clusters based on hierarchical clustering
epi.origin <- readRDS("Res_epi2_noCorrection_10PC_2reso.rds")
cm <- findAllMarkers(epi.origin, selected.column="res.2", core.num=20)
saveRDS(cm, file="Res_epi2_noCorrection_10PC_2reso_cluster_markers.rds")
top.cm <- cm %>% group_by(cluster) %>% top_n(n=5, wt=avg_logFC)
selected.markers <- unique(top.cm$gene_name)
row.num <- ceiling(length(selected.markers)/4)
plot.name <- "Plot_tSNE_epi2_noCorrection_top5ClusterMarkers.png"
png(plot.name, height=row.num*400, width=4*400)
par(mfrow=c(row.num,4))
for(x in selected.markers){
	plotFeature2(coor=epi.origin@dr$tsne@cell.embeddings, values=epi.origin@data[x,], main=x)
}
dev.off()
top.cm <- cm %>% group_by(cluster) %>% top_n(n=50, wt=avg_logFC)
selected.markers <- unique(top.cm$gene_name)
cm.expr <- getAveExpr(meta.data=epi.origin@meta.data, feature.to.calc="res.2", expr=epi.origin@data, genes=selected.markers, core.num=20)
hc <- hclust(as.dist(1-cor(cm.expr)))
pdf("Plot_hcTree_epi_noCorrection_10PC_2reso.pdf")
plot(hc)
dev.off()
# merge clusters based on hierarchical clustering results
old.cl.idx <- hc$order-1
new.cl.idx <- c(1, rep(2,3), 3, rep(4,3), rep(5,3), rep(6,3), rep(7,3), 8, 9, 10, 11, rep(12,2))
merge.cl.idx <- rep(NA, nrow(epi.origin@meta.data))
for(i in unique(new.cl.idx)){
	idx <- old.cl.idx[which(new.cl.idx==i)]
	merge.cl.idx[which(epi.origin@meta.data$res.2%in%idx)] <- i
}
epi.origin@meta.data$merge.cl.idx <- merge.cl.idx
epi.origin@ident <- as.factor(epi.origin@meta.data$merge.cl.idx)
names(epi.origin@ident) <- epi.origin@cell.names
cm <- findAllMarkers(epi.origin, selected.column="merge.cl.idx", core.num=15, cluster.to.test=sort(unique(epi.origin@meta.data$merge.cl.idx)))
saveRDS(cm, file="Res_epi2_noCorrection_10PC_2reso_merged_cluster_markers.rds")
top.cm <- cm %>% group_by(cluster) %>% top_n(n=50, wt=avg_logFC)
selected.markers <- unique(top.cm$gene_name)
selected.markers <- setdiff(selected.markers, confound.genes)
saveRDS(selected.markers, file="Res_epi_merged_cluster_top50_removeCC_rbc_gene.rds")
saveRDS(epi.origin@meta.data, file="Res_epi2_meta_data.rds")
raw.count <- as.matrix(epi.origin@raw.data[,rownames(epi.origin@meta.data)])
saveRDS(raw.count, file="Res_epi2_raw_data.rds")
normalized.count <- as.matrix(epi.origin@data[,rownames(epi.origin@meta.data)])
saveRDS(normalized.count, file="Res_epi2_normalized_data.rds")
saveRDS(epi.origin@scale.data, file="Res_epi2_scaled_data.rds")

## prepare SPRING input
### use scaled expression of top 50 cluster markers, removing cell cycle related and red blood cell marker genes
spring.dir <- "SPRING_cm/"
dir.create(file.path(spring.dir))
hvg.expr <- as.matrix(epi.origin@scale.data[selected.markers, ])
write.table(hvg.expr, file=paste0(spring.dir, "Table_data.csv"), sep=",", col.names=F,quote=F)
hvg.info <- t(epi.origin@meta.data)
write.table(hvg.info, file=paste0(spring.dir, "Table_meta_data.csv"), sep=",", col.names=F,quote=F)
#saveRDS(epi.origin, file="Res_epi2_noCorrection_10PC_2reso.rds")

## generate SPRING plot 
## quantify the significance of between cluster linkage
# input data is scaled expression levels of cluster markers, distance is calcualted based on correlation distance, k=50
dis.mat <- 1-cor(hvg.expr)
plotSPRING(spring.coor.dir="SPRING_cm/epi_noCrrection_10PC_cm_k50/", dis.mat=NULL, genes=markers[,2], marker.plot.name=NULL, expr.mat=epi.origin@scale.data, hvg.info=hvg.info, plot.sample.info=TRUE, sample.plot.name=NULL, plot.cluster.info=TRUE, cluster.plot.name=NULL)
res.origin <- quantifyAndPlotClusterLinkage(data.dir="SPRING_cm/", cell.coor.dir="SPRING_cm/epi_noCrrection_10PC_cm_k50/", distance.type="cor", k=50, cluster.info.idx="merge.cl.idx", to.plot=TRUE, return=TRUE)


# extract cells from 18 week airway sample
## this batch of command is run in Seurat 3.0
home <- "/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data"
work.dir <- "W18_airway"
dir.create(file.path(work.dir))
setwd(work.dir)
cells <- rownames(epi.origin@meta.data)[which(epi.origin@meta.data$Sample=="W18_airway")]
raw.data=epi.origin@raw.data[,cells]
meta.data=epi.origin@meta.data[cells,]
airway.18w <- CreateSeuratObject(counts=raw.data, meta.data=meta.data)
airway.18w <- NormalizeData(object=airway.18w, verbose=FALSE)
airway.18w <- ScaleData(object=airway.18w, verbose=FALSE, vars.to.regress=c("percent.mito", "nUMI"), do.par=T)
airway.18w <- FindVariableFeatures(object=airway.18w, selection.method="vst", verbose=FALSE)
pc.num=10
reso=0.8
airway.18w <- RunPCA(object=airway.18w, npcs=pc.num, verbose=FALSE)
airway.18w <- RunTSNE(object=airway.18w, dims=1:pc.num)
airway.18w <- FindNeighbors(object=airway.18w, dims=1:pc.num)
airway.18w <- FindClusters(object=airway.18w, resolution=reso)
p1 <- DimPlot(object=airway.18w, reduction="tsne", group.by="RNA_snn_res.0.8")
png(paste0("Plot_tSNE_airway18w_cluster_",pc.num,"PC_reso",reso,".png"))
plot(p1)
dev.off()

markers <- read.table("~/Work/Annotation/cellTypeMarker/Lung/Table_manuscript_used_cell_type_markers.txt",sep="\t",stringsAsFactors=F)
markers <- markers[markers[,2]%in%rownames(airway.18w@assays$RNA@data),]
g1 <- markers[,2]
marker.plot.name <- paste0("Plot_tSNE_marker_expression_on_airway18w_",pc.num,"PC_reso",reso,".png")
row.num <- ceiling(length(g1)/4)
png(marker.plot.name, height=row.num*400, width=4*400)
par(mfrow=c(row.num,4))
for(x in g1){
	plotFeature2(coor=Embeddings(object=airway.18w, reduction="tsne"), values=airway.18w@assays$RNA@data[x,], main=x)	
}
dev.off()

cm <- findAllMarkers(airway.18w, selected.column="RNA_snn_res.0.8", core.num=10)
saveRDS(cm, file="Res_airway18w_10PC_reso0.8_cluster_markers.rds")
top.cm <- cm %>% group_by(cluster) %>% top_n(n=50, wt=avg_logFC)
selected.markers <- unique(top.cm$gene_name)
selected.markers <- setdiff(selected.markers, confound.genes)
# prepare SPRING input
spring.dir <- "SPRING_airway18w_cm/"
dir.create(file.path(spring.dir))
hvg.expr <- as.matrix(airway.18w@assays$RNA@scale.data[selected.markers, ])
write.table(hvg.expr, file=paste0(spring.dir, "Table_data.csv"), sep=",", col.names=F,quote=F)
hvg.info <- t(airway.18w@meta.data)
write.table(hvg.info, file=paste0(spring.dir, "Table_meta_data.csv"), sep=",", col.names=F,quote=F)
# generate SPRING plots
## use cluster marker expression levels
#### correlation distance based
disMat=1-cor(hvg.expr)
plotSPRING(spring.coor.dir="SPRING_airway18w_cm/airway18w_10PC_cm_k20/", dis.mat=disMat, genes=markers[,2], marker.plot.name=NULL, expr.mat=airway.18w@assays$RNA@data, hvg.info=hvg.info, plot.sample.info=FALSE, sample.info.column="Sample", sample.plot.name=NULL, plot.cluster.info=TRUE, cluster.info.column="RNA_snn_res.0.8", cluster.plot.name=NULL)
#### Euclidean distance based
disMat=as.matrix(dist(t(hvg.expr)))
plotSPRING(spring.coor.dir="SPRING_airway18w_cm/airway18w_10PC_cm_k20_eucl/", dis.mat=disMat, genes=markers[,2], marker.plot.name=NULL, expr.mat=airway.18w@assays$RNA@data, hvg.info=hvg.info, plot.sample.info=FALSE, sample.info.column="Sample", sample.plot.name=NULL, plot.cluster.info=TRUE, cluster.info.column="RNA_snn_res.0.8", cluster.plot.name=NULL)
### use PC values
spring.dir <- "SPRING_pc/"
dir.create(file.path(spring.dir))
hvg.expr <- t(Embeddings(object=airway.18w, reduction="pca"))
write.table(hvg.expr, file=paste0(spring.dir, "Table_data.csv"), sep=",", col.names=F,quote=F)
hvg.info <- t(airway.18w@meta.data)
write.table(hvg.info, file=paste0(spring.dir, "Table_meta_data.csv"), sep=",", col.names=F,quote=F)
#### Euclidean distance based, k=20
disMat=as.matrix(dist(t(hvg.expr)))
plotSPRING(spring.coor.dir="SPRING_pc/airway18w_10_pc_k20_eucl/", dis.mat=disMat, genes=markers[,2], marker.plot.name=NULL, expr.mat=airway.18w@assays$RNA@data, hvg.info=hvg.info, plot.sample.info=FALSE, sample.info.column="Sample", sample.plot.name=NULL, plot.cluster.info=TRUE, cluster.info.column="RNA_snn_res.0.8", cluster.plot.name=NULL)
#################################################################################################

#####################################################################
# Use results based on non-corrected values for following analysis ##
#####################################################################
epi <- readRDS("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/Res_epi2_noCorrection_10PC_2reso.rds")
dir.create("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/for_manuscript")
setwd("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/for_manuscript")
# tSNE color coded by clusters
plot.coor=epi@dr$tsne@cell.embeddings
group.res=paste0("Cluster", epi@meta.data$merge.cl.idx)
group.cols <- cl.cols
png("Plot_tSNE_noCorrection_merged_cluster_info-2.png", width=1000, height=1000)
plotFeature2(coor=plot.coor, values=group.res, main="Cluster", cex=1.3, cex.main=3, xaxt="n", yaxt="n", bty="n", xlab="", ylab="",gCols=cl.cols)
for(x in unique(group.res)){
	cell.idx <- which(group.res==x)
	vec <- apply(plot.coor[cell.idx, ], 2, median)
	text(vec[1], vec[2], label=cell.types[x], cex=2)
}
dev.off()

# SPRING color coded by clusters
fetal.spring.coor <- read.table("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/SPRING_cm/epi_noCrrection_10PC_cm_k50/coordinates.txt",sep=",",row.names=1)
plot.coor=fetal.spring.coor
group.res=paste0("Cluster", epi@meta.data$merge.cl.idx)
group.cols <- cl.cols
png("Plot_SPRING_noCorrection_merged_cluster_info-2.png", width=1000, height=1000)
plotFeature2(coor=plot.coor, values=group.res, main="Cluster", cex=1.3, cex.main=3, xaxt="n", yaxt="n", bty="n", xlab="", ylab="",gCols=cl.cols)
for(x in unique(group.res)){
	cell.idx <- which(group.res==x)
	vec <- apply(plot.coor[cell.idx, ], 2, median)
	text(vec[1], vec[2], label=sub("Cluster","",x), cex=3)
}
dev.off()

# tSNE color coded by expression of basal cell markers
g1 <- c("TP63", "KRT5", "KRT15")
plot.coor=epi@dr$tsne@cell.embeddings
png("Plot_tSNE_noCorrection_basal_cell_marker_expression.png", width=3000, height=1000)
par(mfrow=c(1,3))
for(x in g1){
	plotFeature2(coor=plot.coor, values=epi@data[x,], main=x, cex=2, cex.main=3, xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
}
dev.off()

expr.mat=as.matrix(epi@data)
ave.expr <- getAveExpr(meta.data=epi@meta.data, feature.to.calc="merge.cl.idx", expr=expr.mat)
saveRDS(ave.expr, file="Res_all_genes_merged_cluster_average_expression.rds")

# identify DE genes between basal cell cluster and bud tip progenitors
library(qusage)
hallmark <- read.gmt("~/Work/Annotation/MSigDB/h.all.v6.1.symbols.gmt")
fetal <- readRDS("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/Res_epi2_normalized_data.rds")
fetal.meta <- readRDS("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/Res_epi2_meta_data.rds")
# cluster 5: bud tip progenitors; cluster 6: bud tip adjacent; cluster 4: hub progenitor; cluster 7: basal cell
cl5.idx <- which(fetal.meta$merge.cl.idx==5)
cl6.idx <- which(fetal.meta$merge.cl.idx==6)
cl4.idx <- which(fetal.meta$merge.cl.idx==4)
cl7.idx <- which(fetal.meta$merge.cl.idx==7)
compare.mat <- matrix(c(6,5,4,6,7,4,7,5), nrow=2)

# DE between bud tip progenitor and basal
de.res.list <- foreach(j=seq(ncol(compare.mat)), .combine="list", .multicombine=T)%dopar%{

	cl.1 <- compare.mat[1,j]
	cl.2 <- compare.mat[2,j]
	cat(paste("Start comparison between cluster", cl.1, "and cluster", cl.2, "\n"))
	cl.idx.1 <- which(fetal.meta$merge.cl.idx==cl.1)
	cl.idx.2 <- which(fetal.meta$merge.cl.idx==cl.2)
	de.res <- t(apply(fetal, 1, function(vec){
		pval <- wilcox.test(vec[cl.idx.1], vec[cl.idx.2])$p.value
		logfc <- mean(vec[cl.idx.1])-mean(vec[cl.idx.2])
		p.cl5 <- sum(vec[cl.idx.1]>0)/length(cl.idx.2)
		p.cl7 <- sum(vec[cl.idx.1]>0)/length(cl.idx.2)
		return(c(pval, logfc, p.cl5, p.cl7))
	}))
	colnames(de.res) <- c("P_value", "logFC", "Prop_1", "Prop_2")
	padj <- p.adjust(de.res[,"P_value"], method="bonferroni")
	de.res <- cbind(de.res, padj)
	saveRDS(de.res, file=paste0("Res_DE_Cl",cl.1,"_Cl",cl.2,".rds"))
	
	## specific interest are placed to genes upregulated in cluster 1
	cl1.up <- rownames(de.res)[which(de.res[,"padj"]<0.05 & de.res[,"logFC"]> 0.1 & de.res[,"Prop_1"]>0.25)]
	cl1_up_genes <- rownames(de.res)%in%cl1.up
	out <- data.frame(de.res, "cl1_up_genes"=cl1_up_genes, stringsAsFactors=F)
	write.table(out, file=paste0("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/for_manuscript/Table_DE_genes_between_Cl",cl.1,"_and_Cl",cl.2,".txt"), quote=F, sep="\t")
	b <- length(cl1.up)
	d <- nrow(fetal)
	hallmark.cl1 <- t(sapply(seq(length(hallmark)), function(j){
		pathway.genes <- intersect(rownames(fetal), hallmark[[j]])
		c <- length(pathway.genes)
		overlap <- intersect(cl1.up, pathway.genes)
		a <- length(overlap)
		pval <- fisher.test(matrix(c(a,b,c,d),c(2,2)), alternative="g")$p.value
		or <- fisher.test(matrix(c(a,b,c,d),c(2,2)), alternative="g")$estimate
		return(c(pval, or, a))
	}))
	rownames(hallmark.cl1) <- names(hallmark)
	hallmark.cl1.adj <- p.adjust(hallmark.cl1[,1], method="bonferroni")
	hallmark.res <- cbind(hallmark.cl1, hallmark.cl1.adj)
	colnames(hallmark.res) <- c("Nonimal_P_value", "Odds_ratio", "Hit_count", "Adjusted_P_value")
	out <- list("DE"=de.res, "pathway"=hallmark.res)
	write.table(hallmark.res, file=paste0("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/for_manuscript/Table_hallmark_geneset_enrichment_for_DE_Cl",cl.1,"_and_Cl",cl.2,".txt"), sep="\t", quote=F)
	return(out)						   

}
names(de.res.list) <- paste0("Cl", compare.mat[1,], "_Cl", compare.mat[2,])
saveRDS(de.res.list, file="Res_basal_lineage_pairwise_DE_comparison_and_pathway_enrichment.rds")
# visualize the pathway enrichment results using heatmap
pathway.res <- sapply(seq(length(de.res.list)), function(i){
	res <- de.res.list[[i]]$pathway
	or.vec <- res[,2]
	pval.vec <- res[,4]
	or.vec[which(pval.vec>0.05)] <- 0
	return(or.vec)

})
colnames(pathway.res) <- names(de.res.list)
cols <- colorRampPalette(c("#d9d9d9", "#252525"))(50)
col.pos <- seq(from=0, to=1, length=ncol(pathway.res))
row.pos <- seq(from=0, to=1, length=nrow(pathway.res))
pdf("Plot_pairwise_pathway_enrichment.pdf")
par(mar=c(6,10,2,2))
image(t(pathway.res), col=cols, xaxt="n", yaxt="n", bty="n")
mtext(text=colnames(pathway.res), side=1, at=col.pos, line=0.5)
mtext(text=sub("HALLMARK_", "", rownames(pathway.res)), side=2, las=2, at=row.pos, cex=0.4, line=0.5)
dev.off()
################################################################################################
#################################### In vitro part #############################################

## no integration
library(Seurat)
library(dplyr)
source("~/Work/commonScript/Script_functions.R")

# load all fetal data
########## Week 11.5 ##########
d0 <- prepareSeuratObject(rawData = Read10X(data.dir = "/mnt/SingleCellGenomics/scg_projects/Spence_Michigan_collab/Lung/Run_2288_Miller_D125_fetalLung_LungBudCulture/Sample_HT-145-3F/outs/HT145-3F-filtered_gene_bc_matrices/hg19/"),
			   namePrefix = "S9",
			   age = 0, 
			   tissue = "Organoid", 
			   sex = NA,
			   sampleName = "D0",
			   condition = NA,
			   isolation = NA,
			   min.gene.cutoff=1000,
			   gene.num.low=1000,
			   seu.version=3
)

d3 <- prepareSeuratObject(rawData=Read10X(data.dir = "/mnt/SingleCellGenomics/scg_projects/Spence_Michigan_collab/Lung/Run_2288_Miller_D125_fetalLung_LungBudCulture/Sample_HT-145-3daysDSA/outs/HT-145-3daysDSA-filtered_gene_bc_matrices/hg19/"),
			   namePrefix = "S10",
			   age = 3, 
			   tissue = "Organoid", 
			   sex = NA,
			   sampleName = "D3",
			   condition = NA,
			   isolation = NA,
			   min.gene.cutoff=1500,
			   gene.num.low=1500,
			   seu.version=3
)

d21 <- prepareSeuratObject(rawData=Read10X(data.dir = "/mnt/SingleCellGenomics/scg_projects/Spence_Michigan_collab/Lung/Run_2321_Miller_HT187lung_HT145DSA-DSI_CzerwinskiHT188AdultDuo/Sample_HT-145-3dDSA-18dFANY/outs/HT145-3dDSA-18dFANY-filtered_gene_bc_matrices/hg19/"),
			   namePrefix = "S11",
			   age = 0, 
			   tissue = "Organoid", 
			   sex = NA,
			   sampleName = "D21",
			   condition = NA,
			   isolation = NA,
			   min.gene.cutoff=1500,
			   gene.num.low=1500,
			   seu.version=3
)

vivo <- prepareSeuratObject(rawData = readRDS("Res_epi2_raw_data.rds"),
			   namePrefix = NULL,
			   age = NA, 			   
			   tissue = "Fetal_combined", 
			   sex = NA,
			   sampleName = "Fetal",
			   condition = NA,
			   isolation = NA,
			   min.gene.cutoff=0,
			   gene.num.low=0,
			   min.cell.cutoff=0,
			   seu.version=3
)

## analyze separately
organoid.list <- list("d0"=d0, "d3"=d3, "d21"=d21)
for(i in seq(length(organoid.list))){
	print(paste(i, "start"))
	seu.obj <- organoid.list[[i]]
	seu.obj <- NormalizeData(object=seu.obj, verbose=FALSE)
	seu.obj <- FindVariableFeatures(object = seu.obj, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
	seu.obj <- ScaleData(object = seu.obj, vars.to.regress = c("nCount_RNA", "percent.mito"), verbose = T, nthreads=20, block.size = 1280000, min.cells.to.block = 100, do.par = T, features = rownames(seu.obj))
	seu.obj <- RunPCA(object = seu.obj, features = VariableFeatures(seu.obj), verbose = T, npcs = 50, ndims.print = 2, nfeatures.print = 2)
	usefulPCs <- 1:10
	seu.obj <- FindNeighbors(object = seu.obj, dims = usefulPCs, force.recalc = T, k.param = 15)
	seu.obj <- FindClusters(object = seu.obj, resolution = 0.8)
	seu.obj <- RunTSNE(object = seu.obj, dims = usefulPCs)
	organoid.list[[i]] <- seu.obj
	file.name <- paste0("Res_", names(organoid.list)[i], ".rds")
	saveRDS(seu.obj, file=file.name)
}

dir.create("in_vitro/")
setwd("in_vitro/")
# plot sample information on tSNE coordiante for each sample separately
for(x in names(organoid.list)){
	seu.obj <- organoid.list[[x]]
	plot.coor=Embeddings(seu.obj, reduction="tsne")
	cluster.vec=paste0("Cluster", seu.obj@meta.data$RNA_snn_res.0.8)
	plot.file <- paste0("Plot_tSNE_",x,"_cluster_info.png")
	png(plot.file, height=1000, width=1000)
	plotFeature2(coor=plot.coor, values=cluster.vec, xaxt="n", bty="n", yaxt="n", main="Cluster", xlab="", ylab="", cex.main=3, cex=2)
	for(x in unique(cluster.vec)){
		cell.idx <- which(cluster.vec==x)
		vec <- apply(plot.coor[cell.idx, ], 2, median)
		text(vec[1], vec[2], label=sub("Cluster","",x), cex=2.5)
	}
	dev.off()
} 

# calculate the transcriptome similarity between in vitro cell and in vivo cluster 
vivo.cl.expr <- readRDS("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/Res_all_genes_merged_cluster_average_expression.rds")
top.vivo.cm <- readRDS("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/Res_epi_merged_cluster_top50_removeCC_rbc_gene.rds")
vitro.detected.genes <- intersect(intersect(rownames(organoid.list[[1]]), rownames(organoid.list[[2]])), rownames(organoid.list[[3]]))
top.vivo.cm <- intersect(top.vivo.cm, vitro.detected.genes)
saveRDS(top.vivo.cm, file="/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/in_vitro/Res_selected_vitro_expressed_fetal_cluster_markers_for_PCC.rds")

cm <- readRDS("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/Res_epi2_noCorrection_10PC_2reso_merged_cluster_markers.rds")
selected.markers <- readRDS("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/Res_epi_merged_cluster_top50_removeCC_rbc_gene.rds")
top.vivo.cm <- readRDS("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/in_vitro/Res_selected_vitro_expressed_fetal_cluster_markers_for_PCC.rds")
fetal_SPRING_selected <- cm$gene_name%in%selected.markers
vitro_vivo_comparison_selected <- cm$gene_name%in%top.vivo.cm
out <- data.frame(cm, "fetal_SPRING_selected"=fetal_SPRING_selected, "vitro_vivo_comparison_selected"=vitro_vivo_comparison_selected, stringsAsFactors=F)
saveRDS(out, file="/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/in_vitro/Res_selected_fetal_cluster_markers_for_fetal_vitro_analysis.rds")
write.table(out, file="/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/for_manuscript/Table_fetal_cluster_markers_for_fetal_vitro_analysis.txt", sep="\t", quote=F)

library(doParallel)
registerDoParallel(20)
vitroCell2vivoCluster.top.cm <- list()
for(x in names(organoid.list)){
	print(paste(x, "start"))
	seu.obj <- organoid.list[[x]]
	cor.res <- foreach(k=seq(ncol(seu.obj@assays$RNA@data)), .multicombine=T, .combine='cbind')%dopar%{
		if(k%%500==0){cat(paste0(k,"start\n"))}
		cor.vec <- sapply(seq(ncol(vivo.cl.expr)), function(cl){
			cor(seu.obj@assays$RNA@data[top.vivo.cm,k], vivo.cl.expr[top.vivo.cm,cl])
		})
		return(cor.vec)
	}
	colnames(cor.res) <- colnames(seu.obj)
	rownames(cor.res) <- colnames(vivo.cl.expr)
	best.cl <- apply(cor.res, 2, function(vec){order(vec, decreasing=T)[1]})
	best.cl=paste0("Cluster", best.cl)
	names(best.cl) <- colnames(cor.res)
	output=list("cell2cluster.cor.mat"=cor.res, "best.cor.cl"=best.cl)
	rds.file <- paste0("Res_cell2Cluster_PCC_", x, ".rds")
	vitroCell2vivoCluster.top.cm[[x]] <- output
	saveRDS(output, file=rds.file)
}
stopImplicitCluster()

# plot the best correlated in vivo cluster distribution on the in vitro tSNE coordinate
for(x in names(organoid.list)){
	print(paste(x, "start"))
	seu.obj <- organoid.list[[x]]
	plot.coor <- Embeddings(seu.obj, reduction="tsne")
	group.res=vitroCell2vivoCluster.top.cm[[x]]$best.cor.cl
	group.cols <- cl.cols[intersect(group.res, names(cl.cols))]
	plot.file <- paste0("Plot_tSNE_vitroCell_vs_vivoCluster_topCM_",x,"-2.png")
	png(plot.file, height=1000, width=1000)
	plotFeature2(coor=plot.coor, values=group.res, xaxt="n", bty="n", yaxt="n", main="", xlab="", ylab="", cex.main=3, cex=2, gCols=group.cols)
	legend("topleft", bty="n", text.col=group.cols, legend=names(group.cols))
	dev.off()
}

# plot the correlation distribution in vivo cluster 
for(x in c("d0", "d3", "d21")){
	print(paste(x, "start"))
	cor.distribution <- vitroCell2vivoCluster.top.cm[[x]]$cell2cluster.cor.mat
	cor.distribution <- t(cor.distribution)
	seu.obj <- organoid.list[[x]]
	#g1 <- intersect(rownames(seu.obj), c("TP63", "KRT5", "KRT15"))
	#expr.mat <- as.matrix(seu.obj@assays$RNA@data[g1,])
	#basal.cell.idx <- colnames(expr.mat)[which(apply(expr.mat, 2, sum)>0)]
	plot.name <- paste0("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/in_vitro/Plot_boxplot_corDistr_vitroCell_vs_vivoCluster_topCM_",x,"-2.png")
	png(plot.name, height=1000, width=2000)
	par(mar=c(10,6,6,6))
	#boxplot(cor.distribution[basal.cell.idx,], main=x, cex.axis=2.5, cex.main=4, frame=F)
	boxplot(cor.distribution, main=x, cex.axis=2.5, cex.main=4, frame=F)
	dev.off()
}


x <- "d3"
cor.distribution <- vitroCell2vivoCluster.top.cm[[x]]$cell2cluster.cor.mat
cor.distribution <- t(cor.distribution)
seu.obj <- organoid.list[[x]]
plot.name <- paste0("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/in_vitro/Plot_boxplot_corDistr_vitroCell_vs_vivoCluster_topCM_",x,"-2.png")
png(plot.name, height=500, width=3000)
par(mar=c(10,6,6,6), mfrow=c(1,3))
for(i in c(0, 5, 6)){
	boxplot(cor.distribution[which(seu.obj@meta.data$RNA_snn_res.0.8==i),], main=paste("Cluster",i,"@",x), cex.axis=2.5, cex.main=4, frame=F, ylim=c(0.3, 0.8))
}
dev.off()


x <- "d21"
cor.distribution <- vitroCell2vivoCluster.top.cm[[x]]$cell2cluster.cor.mat
cor.distribution <- t(cor.distribution)
seu.obj <- organoid.list[[x]]
plot.name <- paste0("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/in_vitro/Plot_boxplot_corDistr_vitroCell_vs_vivoCluster_topCM_",x,"-2.png")
png(plot.name, height=500, width=3000)
par(mar=c(10,6,6,6), mfrow=c(1,3))
for(i in c(0, 9, 3)){
	boxplot(cor.distribution[which(seu.obj@meta.data$RNA_snn_res.0.8==i),], main=paste("Cluster",i,"@",x), cex.axis=2.5, cex.main=4, frame=F, ylim=c(0.3, 0.8))
}
dev.off()

# plot correlation to bud tip progenitors, basal cells, undefined cluster
x <- "d0"
cor.distribution <- vitroCell2vivoCluster.top.cm[[x]]$cell2cluster.cor.mat
cor.distribution <- t(cor.distribution)
seu.obj <- organoid.list[[x]]
#basal.idx <- which(apply(seu.obj@assays$RNA@counts[c("TP63","KRT15"),],2,sum)>2)
basal.idx <- which(seu.obj@assays$RNA@counts["TP63",]>2)
cell.order <- order(cor.distribution[,"Cluster5"], decreasing=T)
y.max <- max(cor.distribution[,c("Cluster5", "Cluster7", "Cluster9")])
y.min <- min(cor.distribution[,c("Cluster5", "Cluster7", "Cluster9")])
plot.name <- paste0("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/in_vitro/Plot_scatter_corDistr_vitroCell_vs_vivoCluster_topCM_",x,".png")
png(plot.name, height=1000, width=3000)
par(mar=c(6,6,5,3))
plot(seq(length(cell.order)), cor.distribution[cell.order, "Cluster5"], pch=16, col=ifelse(cell.order%in%basal.idx, cl.cols["Cluster5"], paste0(cl.cols["Cluster5"],"20")), ylab="PCC", xlab="Order", cex.lab=2, cex.axis=1.5, ylim=c(y.min, y.max), main=x, cex.main=2.5)
points(seq(length(cell.order)), cor.distribution[cell.order, "Cluster7"], pch=16, col=ifelse(cell.order%in%basal.idx, cl.cols["Cluster7"], paste0(cl.cols["Cluster7"],"20")))
points(seq(length(cell.order)), cor.distribution[cell.order, "Cluster9"], pch=16, col=ifelse(cell.order%in%basal.idx, cl.cols["Cluster9"], paste0(cl.cols["Cluster9"],"20")))
legend("topright", legend=c("PCC to Cl5 & TP63+","PCC to Cl5 & TP63-","PCC to Cl7 & TP63+","PCC to Cl7 & TP63-","PCC to Cl9 & TP63+","PCC to Cl9 & TP63-"), bty="n", text.col=c(cl.cols["Cluster5"], paste0(cl.cols["Cluster5"],"70"), cl.cols["Cluster7"], paste0(cl.cols["Cluster7"],"70"),cl.cols["Cluster9"], paste0(cl.cols["Cluster9"],"70")), cex=2)
dev.off()

plot.coor=Embeddings(seu.obj, reduction="tsne")
plot.name <- paste0("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/in_vitro/Plot_selected_basal_cell_",x,".png")
png(plot.name, height=1000, width=1000)
plot(plot.coor, pch=19, col="gray")
points(plot.coor[basal.idx,], pch=19, col="blue")
dev.off()

vitroCell2vivoCluster.top.cm <- list()
for(x in names(organoid.list)){
	rds.file <- paste0("Res_cell2Cluster_PCC_", x, ".rds")
	res <- readRDS(rds.file)
	vitroCell2vivoCluster.top.cm[[x]] <- res
}

# plot the cluster 9 and basal cell markers for those day 3 cells with highest correlation with cluster 9
vivo.cm.mat <- readRDS("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/Res_epi2_noCorrection_10PC_2reso_merged_cluster_markers.rds")
selected.day <- "d0"
cl9 <- intersect(vivo.cm.mat$gene_name[which(vivo.cm.mat$cluster==9)], top.vivo.cm)
x <- paste(cl9, "9", sep="_")
cl9.mat <- vivo.cm.mat[x,]
cl9.markers <- intersect(cl9.mat$gene_name[order(cl9.mat$avg_logFC, decreasing=T)][1:10], rownames(organoid.list[[selected.day]]))

cl7 <- intersect(vivo.cm.mat$gene_name[which(vivo.cm.mat$cluster==7)], top.vivo.cm)
x <- paste(cl7, "7", sep="_")
cl7.mat <- vivo.cm.mat[x,]
cl7.markers <- intersect(cl7.mat$gene_name[order(cl7.mat$avg_logFC, decreasing=T)][1:10], rownames(organoid.list[[selected.day]]))

selected.cl <- 5
cm0 <- intersect(vivo.cm.mat$gene_name[which(vivo.cm.mat$cluster==selected.cl)], top.vivo.cm)
x <- paste(cm0, as.character(selected.cl), sep="_")
cm.mat <- vivo.cm.mat[x,]
cl5.markers <- intersect(cm.mat$gene_name[order(cm.mat$avg_logFC, decreasing=T)][1:10], rownames(organoid.list[[selected.day]]))

g1 <- c(cl9.markers, cl7.markers, cl5.markers, c("TP63", "KRT15"))
idx <- c(rep("Cluster9", length(cl9.markers)), rep("Cluster7", length(cl7.markers)), rep("Cluster5", length(cl5.markers)),rep("Basal cell", 2))

setwd("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/in_vitro")
plot.file <- paste0("Plot_tSNE_", selected.day, "_vivo_cluster9_cluster7_cluster5_marker_expression.png")
plot.coor=Embeddings(organoid.list[[selected.day]], reduction="tsne")
column.num <- 8
row.num <- ceiling(length(g1)/column.num)
png(plot.file, height=1000*row.num, width=1000*column.num)
par(mfrow=c(row.num, column.num))
for(i in seq(length(g1))){
	g.to.plot <- g1[i]
	g.idx <- idx[i]
	plotFeature2(coor=plot.coor, values=organoid.list[[selected.day]]@assays$RNA@data[g.to.plot,], xaxt="n", bty="n", yaxt="n", main=paste(g.to.plot, g.idx, sep="@"), xlab="", ylab="", cex.main=5, cex=2)
}
dev.off()

vitro.combined <- merge(x=organoid.list[["d0"]], y=list(organoid.list[["d3"]], organoid.list[["d21"]]))
saveRDS(vitro.combined, file="Res_vitro_combined.rds")

vitro.expressed.genes <- rownames(vitro.combined)
cl9.markers <- intersect(cl9.mat$gene_name[order(cl9.mat$avg_logFC, decreasing=T)][1:20], vitro.expressed.genes)
cl7.markers <- intersect(cl7.mat$gene_name[order(cl7.mat$avg_logFC, decreasing=T)][1:20], vitro.expressed.genes)
cl5.markers <- intersect(cm.mat$gene_name[order(cm.mat$avg_logFC, decreasing=T)][1:20], vitro.expressed.genes)
#g1 <- c(cl9.markers, cl7.markers, cl5.markers, c("TP63", "KRT15"))
#idx <- c(rep("Cluster9", length(cl9.markers)), rep("Cluster7", length(cl7.markers)), rep("Cluster5", length(cl5.markers)),rep("Basal cell", 2))

marker.list <- list("Cl5"=cl5.markers, "Cl7"=cl7.markers, "Cl9"=cl9.markers, "Basal"=c("TP63", "KRT15"))
cell.type.markers <- read.table("~/Work/Annotation/cellTypeMarker/Lung/Table_manuscript_used_cell_type_markers.txt",sep="\t",stringsAsFactors=F)

expr.mat <- as.matrix(vitro.combined@assays$RNA@data)
age.vec <- rep(c("Day0", "Day3", "Day21"), c(ncol(organoid.list[["d0"]]), ncol(organoid.list[["d3"]]), ncol(organoid.list[["d21"]])))

d0.coor=Embeddings(organoid.list[["d0"]], reduction="tsne")
d3.coor=Embeddings(organoid.list[["d3"]], reduction="tsne")
d21.coor=Embeddings(organoid.list[["d21"]], reduction="tsne")
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))
column.num <- 9
row.num <- ceiling(length(g1)*3/column.num)

#for(cl.idx in names(marker.list)){
	#g1 <- marker.list[[cl.idx]]
	#column.num <- 9
	#plot.file <- paste0("Plot_tSNE_vivo_",cl.idx,"_marker_expression.png")
	g1 <- cell.type.markers[,2]
	column.num <- 6
	row.num <- ceiling(length(g1)*3/column.num)
	plot.file <- "Plot_tSNE_canonical_cell_type_marker_expression.png"
	png(plot.file, height=1000*row.num, width=1000*column.num)
	par(mfrow=c(row.num, column.num))
	for(g.to.plot in g1){
		values=expr.mat[g.to.plot,]
		all.cell.cols <- adjustcolor(colorPal(30), alpha=.8)[as.numeric(cut(values, breaks=30, right=F, include.lowest=T))]
		plot(d0.coor, col=all.cell.cols[which(age.vec=="Day0")], xaxt="n", bty="n", yaxt="n", main=paste(g.to.plot, cl.idx, "Day 0", sep="@"), xlab="", ylab="", cex.main=5, cex=3, pch=16)
		plot(d3.coor, col=all.cell.cols[which(age.vec=="Day3")], xaxt="n", bty="n", yaxt="n", main=paste(g.to.plot, cl.idx, "Day 3", sep="@"), xlab="", ylab="", cex.main=5, cex=3, pch=16)
		plot(d21.coor, col=all.cell.cols[which(age.vec=="Day21")], xaxt="n", bty="n", yaxt="n", main=paste(g.to.plot, cl.idx, "Day 21", sep="@"), xlab="", ylab="", cex.main=5, cex=3, pch=16)
	}
	dev.off()
#}

fetal <- readRDS("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/Res_epi2_normalized_data.rds")
fetal.spring.coor <- read.table("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/SPRING_cm/epi_noCrrection_10PC_cm_k50/coordinates.txt",sep=",",row.names=1)
rownames(fetal.spring.coor) <- colnames(fetal)
for(cl.idx in names(marker.list)){
	g1 <- marker.list[[cl.idx]]
	plot.file <- paste0("Plot_tSNE_fetal_vivo_", cl.idx, "_marker_expression.png")
	plot.coor=fetal.spring.coor
	expr.mat <- fetal
	column.num <- 10
	row.num <- ceiling(length(g1)/column.num)
	png(plot.file, height=1000*row.num, width=1000*column.num)
	par(mfrow=c(row.num, column.num))
	for(g.to.plot in g1){
		plotFeature2(coor=plot.coor, values=expr.mat[g.to.plot,], xaxt="n", bty="n", yaxt="n", main=paste(g.to.plot, cl.idx, sep="@"), xlab="", ylab="", cex.main=5, cex=2)
	}
	dev.off()
}

# plot gene expression patterns in fetal and organoid data
# specific gene sets
expressed.gene <- intersect(rownames(vitro.combined), rownames(fetal))
tgf.genes <- setdiff(grep("TGF", expressed.gene, value=T),c("PTGFRN","CTGF","PTGFR"))
bmp.genes <- grep("BMP", expressed.gene, value=T)
acvr.genes <- grep("ACV",expressed.gene, value=T)
smad.genes <- grep("SMAD",expressed.gene, value=T)
signal.genes <- intersect(c(tgf.genes, bmp.genes, acvr.genes, smad.genes, "NODAL", "NOG","CHRD"), expressed.gene)

column.num <- 4
row.num <- length(signal.genes)
png("Plot_feature_plot_SMAD_signaling_genes.png", height=500*row.num, width=500*column.num)
par(mfrow=c(row.num, column.num))
for(x in signal.genes){
	# fetal data
	plotFeature2(coor=fetal.spring.coor, values=fetal[x,], xaxt="n", bty="n", yaxt="n", main=paste(x, "Fetal", sep="@"), xlab="", ylab="", cex.main=5, cex=2)

	# organoid data
	values=expr.mat[x,]
	all.cell.cols <- adjustcolor(colorPal(30), alpha=.8)[as.numeric(cut(values, breaks=30, right=F, include.lowest=T))]
	all.cell.cols[which(values==0)] <- "#bdbdbd30"
	plot(d0.coor, col=all.cell.cols[which(age.vec=="Day0")], xaxt="n", bty="n", yaxt="n", main=paste(x, "Day 0", sep="@"), xlab="", ylab="", cex.main=5, cex=3, pch=16)
	plot(d3.coor, col=all.cell.cols[which(age.vec=="Day3")], xaxt="n", bty="n", yaxt="n", main=paste(x, "Day 3", sep="@"), xlab="", ylab="", cex.main=5, cex=3, pch=16)
	plot(d21.coor, col=all.cell.cols[which(age.vec=="Day21")], xaxt="n", bty="n", yaxt="n", main=paste(x, "Day 21", sep="@"), xlab="", ylab="", cex.main=5, cex=3, pch=16)
}
dev.off()

# barplots for average expression levels of TP63, KRT5 and KRT15 in fetal data
pdf("Plot_barplot_basal_cell_marker_expr_across_in_vivo_clusters.pdf")
par(mfrow=c(3,1))
barplot(vivo.cl.expr["TP63",], main="TP63", bty="n")
barplot(vivo.cl.expr["KRT15",], main="KRT15", bty="n")
barplot(vivo.cl.expr["KRT5",], main="KRT5", bty="n")
dev.off()

# plot the sum expression patterns of bud tip progenitor, basal cell cluster in fetal and organoid data
vivo.cm.mat <- readRDS("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/Res_epi2_noCorrection_10PC_2reso_merged_cluster_markers.rds")
vitro.expressed.genes <- rownames(vitro.combined)
cl7.markers <- intersect(vivo.cm.mat$gene_name[which(vivo.cm.mat$cluster==7)], vitro.expressed.genes)
cl5.markers <- intersect(vivo.cm.mat$gene_name[which(vivo.cm.mat$cluster==5)], vitro.expressed.genes)

marker.list <- list("Cl5"=cl5.markers, "Cl7"=cl7.markers)

expr.mat <- as.matrix(vitro.combined@assays$RNA@data)
age.vec <- rep(c("Day0", "Day3", "Day21"), c(ncol(organoid.list[["d0"]]), ncol(organoid.list[["d3"]]), ncol(organoid.list[["d21"]])))

fetal <- readRDS("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/Res_epi2_normalized_data.rds")
fetal.spring.coor <- read.table("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/SPRING_cm/epi_noCrrection_10PC_cm_k50/coordinates.txt",sep=",",row.names=1)
rownames(fetal.spring.coor) <- colnames(fetal)
d0.coor=Embeddings(organoid.list[["d0"]], reduction="tsne")
d3.coor=Embeddings(organoid.list[["d3"]], reduction="tsne")
d21.coor=Embeddings(organoid.list[["d21"]], reduction="tsne")

colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))
column.num <- 4
row.num <- length(marker.list)
png("Plot_fetal_and_organoid_Cl5_Cl7_marker_scale_sum_expression.png", height=1000*row.num, width=1000*column.num)
par(mfrow=c(row.num, column.num))
for(cl.idx in names(marker.list)){
	g1 <- marker.list[[cl.idx]]
	fetal.values=apply(t(scale(t(fetal[g1,]))), 2, sum)
	plotFeature2(coor=fetal.spring.coor, values=fetal.values, xaxt="n", bty="n", yaxt="n", main=paste(cl.idx, "Fetal", sep="@"), xlab="", ylab="", cex.main=5, cex=2)
	
	organoid.values=apply(t(scale(t(expr.mat[g1,]))), 2, sum)
	organoid.cell.cols <- adjustcolor(colorPal(30), alpha=.8)[as.numeric(cut(organoid.values, breaks=30, right=F, include.lowest=T))]
	plot(d0.coor, col=organoid.cell.cols[which(age.vec=="Day0")], xaxt="n", bty="n", yaxt="n", main=paste(cl.idx, "Day 0", sep="@"), xlab="", ylab="", cex.main=5, cex=3, pch=16)
	plot(d3.coor, col=organoid.cell.cols[which(age.vec=="Day3")], xaxt="n", bty="n", yaxt="n", main=paste(cl.idx, "Day 3", sep="@"), xlab="", ylab="", cex.main=5, cex=3, pch=16)
	plot(d21.coor, col=organoid.cell.cols[which(age.vec=="Day21")], xaxt="n", bty="n", yaxt="n", main=paste(cl.idx, "Day 21", sep="@"), xlab="", ylab="", cex.main=5, cex=3, pch=16)
}
dev.off()

# plot the difference of average scaled expression levels of Cl7 and Cl5
png("Plot_fetal_and_organoid_Cl7_Cl5_marker_average_scale_sum_expression_difference.png", height=1000, width=1000*4)
par(mfrow=c(1, 4))
cm.5 <- marker.list[["Cl5"]]
fetal.cl5=apply(t(scale(t(fetal[cm.5,]))), 2, mean)
cm.7 <- marker.list[["Cl7"]]
fetal.cl7=apply(t(scale(t(fetal[cm.7,]))), 2, mean)
fetal.values=fetal.cl7-fetal.cl5
plotFeature2(coor=fetal.spring.coor, values=fetal.values, xaxt="n", bty="n", yaxt="n", main="Fetal", xlab="", ylab="", cex.main=5, cex=2)

organoid.cl5=apply(t(scale(t(expr.mat[cm.5,]))), 2, mean)
organoid.cl7=apply(t(scale(t(expr.mat[cm.7,]))), 2, mean)
organoid.values=organoid.cl7-organoid.cl5
organoid.cell.cols <- adjustcolor(colorPal(30), alpha=.8)[as.numeric(cut(organoid.values, breaks=30, right=F, include.lowest=T))]
plot(d0.coor, col=organoid.cell.cols[which(age.vec=="Day0")], xaxt="n", bty="n", yaxt="n", main="Day 0", xlab="", ylab="", cex.main=5, cex=3, pch=16)
plot(d3.coor, col=organoid.cell.cols[which(age.vec=="Day3")], xaxt="n", bty="n", yaxt="n", main="Day 3", xlab="", ylab="", cex.main=5, cex=3, pch=16)
plot(d21.coor, col=organoid.cell.cols[which(age.vec=="Day21")], xaxt="n", bty="n", yaxt="n", main="Day 21", xlab="", ylab="", cex.main=5, cex=3, pch=16)
dev.off()

# plot the correlation difference 
vitro.pcc <- cbind(vitroCell2vivoCluster.top.cm[["d0"]]$cell2cluster.cor.mat, vitroCell2vivoCluster.top.cm[["d3"]]$cell2cluster.cor.mat, vitroCell2vivoCluster.top.cm[["d21"]]$cell2cluster.cor.mat)
png("Plot_organoid_to_Cl7_Cl5_PCC_difference.png", height=1000, width=1000*3)
par(mfrow=c(1, 3))
organoid.cl5=vitro.pcc["Cluster5",]
organoid.cl7=vitro.pcc["Cluster7",]
organoid.values=organoid.cl7-organoid.cl5
organoid.cell.cols <- adjustcolor(colorPal(30), alpha=.8)[as.numeric(cut(organoid.values, breaks=30, right=F, include.lowest=T))]
plot(d0.coor, col=organoid.cell.cols[which(age.vec=="Day0")], xaxt="n", bty="n", yaxt="n", main="Day 0", xlab="", ylab="", cex.main=5, cex=3, pch=16)
plot(d3.coor, col=organoid.cell.cols[which(age.vec=="Day3")], xaxt="n", bty="n", yaxt="n", main="Day 3", xlab="", ylab="", cex.main=5, cex=3, pch=16)
plot(d21.coor, col=organoid.cell.cols[which(age.vec=="Day21")], xaxt="n", bty="n", yaxt="n", main="Day 21", xlab="", ylab="", cex.main=5, cex=3, pch=16)
dev.off()

# plot the heatmap showing distribution of correlation to different in vivo clusters for each day in vitro cells
# for day 3
sample.cor.mat <- t(vitroCell2vivoCluster.top.cm[["d3"]]$cell2cluster.cor.mat)
bestVivoCluster2VItroCell.vec <- vitroCell2vivoCluster.top.cm[["d3"]]$best.cor.cl
tp63.cell.idx <- which(organoid.list[["d3"]]@assays$RNA@data["TP63",]>0)
cl5.cl6.cell.idx <- which(organoid.list[["d3"]]@meta.data$RNA_snn_res.0.8%in%c(5,6))
selected.cell.idx <- cl5.cl6.cell.idx
selected.cor.mat <- sample.cor.mat[selected.cell.idx,]
bestVivoCluster2VItroCell.top.cm <- bestVivoCluster2VItroCell.vec[selected.cell.idx]

# for day 21
sample.cor.mat <- t(vitroCell2vivoCluster.top.cm[["d21"]]$cell2cluster.cor.mat)
bestVivoCluster2VItroCell.vec <- vitroCell2vivoCluster.top.cm[["d21"]]$best.cor.cl
cl3.cell.idx <- which(organoid.list[["d21"]]@meta.data$RNA_snn_res.0.8==3)
selected.cell.idx <- cl3.cell.idx
selected.cor.mat <- sample.cor.mat[selected.cell.idx,]
bestVivoCluster2VItroCell.top.cm <- bestVivoCluster2VItroCell.vec[selected.cell.idx]

normed.mat <- t(apply(selected.cor.mat, 1, function(vec){
	(vec-min(vec))/(max(vec)-min(vec))
}))
#hc.col <- hclust(as.dist(1-cor(normed.mat)), method="ward.D2")
hc.col <- hclust(as.dist(1-cor(vivo.cl.expr[top.vivo.cm,])), method="ward.D2")
#hc.col$order <- c(1,2,7,8,9,10,11,13)
hc.row <- hclust(as.dist(1-cor(t(normed.mat), method="spearman")), method="ward.D2")
vitro.cell.idx <- hc.row$order
values=bestVivoCluster2VItroCell.top.cm
vitro.cell.idx.2 <- unlist(lapply(paste0("Cluster",1:13), function(x){
	intersect(vitro.cell.idx, which(values==x))
}))
colors <- colorRampPalette(c("#2166ac","#67a9cf","#d1e5f0","#f7f7f7","#fddbc7","#ef8a62","#b2182b"), space="rgb")(50)

#cl.cols <- setNames(scales::hue_pal()(13), paste0("Cluster", seq(13)))
#cl.cols <- c("#197636", "#7BCAA4", "#A4D371", "#EED951", "#DCCB7C", "#FDD884", "#CA6778","#9E0142","#5E4FA2", "#A84798", "#3C7AB6","#4DA7B0")
cl.cols <- c("#197636", "#7BCAA4", "#A4D371", "#EED951", "#DCCB7C", "#FDD884", "#4DA7B0", "#CA6778","#9E0142","#5E4FA2", "#A84798", "#3C7AB6")
names(cl.cols) <- paste0("Cluster", seq(12))
cell.col <- cl.cols[values[vitro.cell.idx.2]]
#cell.types <- c("Ciliated", "Proliferating ciliated", "Undefined", "Differentiating basal cell", "Secretory progenitor", "Mixture", "Neuroendocrine", "Bud tip adjacent", "Bud tip progenitor", "Hub progenitor", "Undefined proliferating cell", "Submucosal gland",  "Basal cell")
#names(cell.types) <- paste0("Cluster", seq(13))

library(gplots)
#file.name <- "Plot_heatmap_d3_cl5_cl6_cell_vivoCell2VitroCluster_topCM.png"
file.name <- "Plot_heatmap_d21_cl3_cell_vivoCell2VitroCluster_topCM.png"
png(file.name, height=500, width=500)
par(mar=c(6,6,5,3))
#heatmap.2(normed.mat[vitro.cell.idx.2, hc.col$order],trace="none",density.info="none",scale="row",dendrogram="none", col=colors,Rowv=FALSE, Colv=FALSE,labRow=NA,labCol=cell.types[hc.col$order],ColSideColors=cl.cols[hc.col$order], RowSideColors=cell.col)
heatmap.2(normed.mat[vitro.cell.idx.2, hc.col$order],trace="none",density.info="none",scale="row",dendrogram="none", col=colors,Rowv=FALSE, Colv=FALSE,labRow=NA)
dev.off()

# feature plots
markers <- c("FOXJ1","TMEM190","CDK1","MKI67","ASCL1","SFTPC", "SCGB3A2", "SFTPB", "TP63", "KRT15")
markers <- c("MUC5AC", "MUC5B", "SCGB1A1", "CHGA", "SPDEF", "CFTR", "EPCAM", "SOX9", "KRT5", "EGFR", "F3", "PDPN")

colorPal <- grDevices::colorRampPalette(c("navy", "darkorange1"))
png("Plot_tSNE_marker_gene_expression_bo-additional.png", height=4500, width=6000)
par(mfrow=c(3,4))
for(x in markers){
	plotFeature2(tsne.coor, values=d21@data[x, ], xaxt="n", yaxt="n", bty="n", main=x, xlab="", ylab="",cex.main=5,nCols=c("navy", "darkorange1"),cex=4)
}
#for(i in 1:30) rect(20+0.5*(i-1),30,20+0.5*i,33, border=NA, col=colorPal(30)[i])
dev.off()

# calculate correlation between the day 3 in vitro cell and fetal cells
top.vivo.cm <- readRDS("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/Res_epi_merged_cluster_top50_removeCC_rbc_gene.rds")
vitro.detected.genes <- intersect(intersect(rownames(organoid.list[[1]]), rownames(organoid.list[[2]])), rownames(organoid.list[[3]]))
top.vivo.cm <- intersect(top.vivo.cm, vitro.detected.genes)


library(doParallel)
registerDoParallel(20)
vitroCell2vivoCell.top.cm <- list()
for(x in names(organoid.list)){
	print(paste(x, "start"))
	seu.obj <- organoid.list[[x]]
	vitro.expr <- as.matrix(seu.obj@assays$RNA@data[top.vivo.cm,])
	cor.res <- cor(vitro.expr, fetal[top.vivo.cm,])
	rownames(cor.res) <- colnames(seu.obj)
	colnames(cor.res) <- colnames(fetal)
	
	rds.file <- paste0("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/in_vitro/Res_cell2Cell_PCC_", x, ".rds")
	vitroCell2vivoCell.top.cm[[x]] <- cor.res
	saveRDS(cor.res, file=rds.file)
}
stopImplicitCluster()
# project the in vitro cells to in vivo SPRING coordinate
best.cell.idx <- apply(cor.res, 1, function(vec){
	order(vec, decreasing=T)[1:10]
})
projected.coor <- t(apply(best.cell.idx, 2, function(vec){
	apply(fetal.spring.coor[vec,], 2, mean)
}))

png("Plot_SPRING_projection_d3.png", height=1000, width=1000)
plot(fetal.spring.coor, pch=16, col="gray")
points(projected.coor, pch=16, col="blue")
dev.off()

########################### characterize basal cell heterogeneity ##################################
# perform subclustering for cluster 7 and 9, the potential basal cell cluster
dir.create("basal_include_cl9")
setwd("/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/basal_include_cl9")
cell.idx <- rownames(epi@meta.data)[which(epi@meta.data$merge.cl.idx%in%c(7,9))]
count <- as.matrix(epi@raw.data[, cell.idx])
metadata <- epi@meta.data[cell.idx,]
basal <- CreateSeuratObject(counts = count, meta.data = metadata)
basal@assays$RNA@data <- epi@data[, cell.idx]
basal@assays$RNA@scale.data <- epi@scale.data[, cell.idx]
basal <- FindVariableFeatures(object = basal, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
pc.num=5
reso=0.6
cc.genes <- readLines("/home/qianhui_yu/Work/Annotation/cellCycle/regev_lab_cell_cycle_genes.txt")
confound.genes <- c(cc.genes, c("HBB", "HBA1", "HBA2", "HBG1", "HBG2"))
selected.hvg <- setdiff(VariableFeatures(basal), confound.genes)
basal <- RunPCA(object=basal, npcs=pc.num, verbose=FALSE, features=selected.hvg)
basal <- RunTSNE(object=basal, dims=1:pc.num)
basal <- FindNeighbors(object=basal, dims=1:pc.num)
basal <- FindClusters(object=basal, resolution=reso)
saveRDS(basal, file="Res_epi2_cluster7_9_basal_cell.rds")

png("Plot_tSNE_basal_cell_5PC_reso0.6.png", height=1000, width=3000)
par(mfrow=c(1,3))
plotFeature2(Embeddings(basal, reduction="tsne"), values=paste0("Cluster", basal@meta.data$RNA_snn_res.0.6), bty="n", xaxt="n", yaxt="n", xlab="", ylab="", main="Cluster", add.legend=T, add.label=T, cex=2, label.cex=2, legend.cex=2, cex.main=3)
plotFeature2(Embeddings(basal, reduction="tsne"), values=paste0("Cluster", basal@meta.data$merge.cl.idx), bty="n", xaxt="n", yaxt="n", xlab="", ylab="", main="Original cluster", add.legend=T, add.label=T, cex=2, label.cex=2, legend.cex=2, cex.main=3,gCols=cl.cols)
plotFeature2(Embeddings(basal, reduction="tsne"), values=basal@meta.data$Sample, bty="n", xaxt="n", yaxt="n", xlab="", ylab="", main="Sample", add.legend=T, add.label=F, cex=2, label.cex=2, legend.cex=2, cex.main=3)
dev.off()

basal.cm <- findAllMarkers(seu.obj=basal, seu.version=3, selected.column="RNA_snn_res.0.6")
saveRDS(basal.cm, file="Res_basal_cell_no_integration_subclustering_markers_5PC_0.6_removeConfoundGenes.rds")
top.cm <- basal.cm %>% group_by(cluster) %>% top_n(n=50, wt=avg_logFC)
selected.markers <- setdiff(unique(top.cm$gene_name), confound.genes)
mat <- read.table("~/Work/Annotation/cellTypeMarker/Lung/Table_manuscript_used_cell_type_markers.txt",sep="\t",stringsAsFactors=F)
selected.genes <- unique(mat[,2])
selected.genes <- selected.markers

# calculate the transcriptome similarity between individual basal cells
basal.expr <- as.matrix(basal@assays$RNA@data)
basal.c2c.cor <- cor(basal.expr, basal.expr)
basal.c2c.cor.hvg <- cor(basal.expr[VariableFeatures(basal),], basal.expr[VariableFeatures(basal),])
basal.c2c.cor.cm <- cor(basal.expr[selected.markers,], basal.expr[selected.markers,])
#basal.ave.expr <- getAveExpr(meta.data=basal@meta.data, feature.to.calc="RNA_snn_res.0.6", expr=basal.expr)


cl.pair.idx <- c()
for(cl.1 in sort(unique(basal@meta.data$RNA_snn_res.0.6))){
	for(cl.2 in sort(unique(basal@meta.data$RNA_snn_res.0.6))){
		#if(as.numeric(as.character(cl.2))>=as.numeric(as.character(cl.1))){
			cl.pair.idx <- rbind(cl.pair.idx, c(as.numeric(as.character(cl.1)), as.numeric(as.character(cl.2))))
		#}
	}
}

selected.cor.mat <- basal.c2c.cor
selected.cor.mat <- basal.c2c.cor.cm
selected.cor.mat <- basal.c2c.cor.hvg
cl.pair.cor <- t(sapply(seq(nrow(cl.pair.idx)), function(i){
	cl.1 <- cl.pair.idx[i,1]
	cl.2 <- cl.pair.idx[i,2]
	cl.1.cell <- which(basal@meta.data$RNA_snn_res.0.6==cl.1)
	cl.2.cell <- which(basal@meta.data$RNA_snn_res.0.6==cl.2)
	vec <- as.vector(selected.cor.mat[cl.1.cell, cl.2.cell])
	return(c(mean(vec), sd(vec)))
}))
cols <- colorRampPalette(c("#d9d9d9", "#252525"))(nrow(cl.pair.cor))
cor.mean.mat <- matrix(cl.pair.cor[,1], nrow=length(unique(basal@meta.data$RNA_snn_res.0.6)))
colnames(cor.mean.mat) <- paste0("Cluster", sort(unique(basal@meta.data$RNA_snn_res.0.6)))
rownames(cor.mean.mat) <- paste0("Cluster", sort(unique(basal@meta.data$RNA_snn_res.0.6)))
hc <- hclust(as.dist(1-cor.mean.mat))

ordered.cor.mean.mat <- cor.mean.mat[hc$order,hc$order]
values=round(as.vector(ordered.cor.mean.mat),2)
cl.num = length(unique(basal@meta.data$RNA_snn_res.0.6))
pos = seq(from=0, to=1, length=cl.num)
pdf("Plot_basal_cell_cluster_transcriptome_PCC_cm.pdf")
image(t(ordered.cor.mean.mat), col=cols, xaxt="n", yaxt="n", bty="n")
text(rep(pos, cl.num), rep(pos, each=cl.num), labels=values)
mtext(text=colnames(ordered.cor.mean.mat), side=1, at=pos)
mtext(text=rownames(ordered.cor.mean.mat), side=2, at=pos)
dev.off()

plot.coor=Embeddings(basal, reduction="tsne")
column.num <- 8
row.num <- ceiling(length(selected.genes)/column.num)
png(paste0("Plot_tSNE_basal_cell_subclustering_markers_",pc.num,"PC_reso",reso,".png"), height=row.num*500, width=column.num*500)
par(mfrow=c(row.num, column.num))
for(x in selected.genes){
	plotFeature2(coor=plot.coor, values=basal@assays$RNA@data[x,], main=x, cex=2, cex.main=3, xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
}
dev.off()

# plot top 5~10 subcluster markers
top.cm <- basal.cm %>% group_by(cluster) %>% top_n(n=8, wt=avg_logFC)
selected.markers <- setdiff(unique(top.cm$gene_name), confound.genes)
selected.top.cm <- intersect(selected.markers, intersect(rownames(d3), rownames(d21)))
plot.coor=Embeddings(basal, reduction="tsne")
column.num <- 10
row.num <- ceiling(length(selected.top.cm)/column.num)
png("Plot_tSNE_basal_cell_subclustering_markers.png", height=row.num*1000, width=column.num*1000)
par(mfrow=c(row.num, column.num))
for(x in selected.top.cm){
	plotFeature2(coor=plot.coor, values=basal@assays$RNA@data[x,], main=x, cex=5, cex.main=7, xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
}
dev.off()


basal.expr.mat <- as.matrix(basal@assays$RNA@data)
basal.cl.expr <- getAveExpr(meta.data=basal@meta.data, feature.to.calc="RNA_snn_res.0.8", expr=basal.expr.mat)
colnames(basal.cl.expr) <- paste("Basal", colnames(basal.cl.expr), sep="_")
saveRDS(basal.cl.expr, file="Res_basal_cell_no_integration_cluster_average_expr.rds")

# pull out TP63+ cells in day 3 and day 21, and compared them to different subtypes
g1 <- intersect(selected.markers, intersect(rownames(d3), rownames(d21)))
basal.subtype.ref <- basal.cl.expr[g1,]
d3.basal.cells <- colnames(d3)[which(d3@assays$RNA@counts["TP63",]>2)]
d21.basal.cells <- intersect(names(which(vitroCell2vivoCluster.top.cm[["d21"]]$best.cor.cl=="Cluster7")), rownames(d21@meta.data)[which(d21@meta.data$RNA_snn_res.0.8==3)])
png("Plot_tSNE_selected_organoid_basal_cells.png", height=1000, width=2000)
par(mfrow=c(1,2))
plot(Embeddings(d3, reduction="tsne"), pch=16, col="gray", main="Day 3")
points(Embeddings(d3, reduction="tsne")[d3.basal.cells,], pch=16, col="blue")
plot(Embeddings(d21, reduction="tsne"), pch=16, col="gray", main="Day 21")
points(Embeddings(d21, reduction="tsne")[d21.basal.cells,], pch=16, col="blue")
dev.off()

selected.for.vitro.compare <- basal.cm$gene_name%in%rownames(basal.subtype.ref)
out <- data.frame(basal.cm, "selected_for_vitro_compare"=selected.for.vitro.compare, stringsAsFactors=F)
write.table(out, file="/r1/people/qianhui_yu/Work/Lung/include_new_fetal_data/for_manuscript/Table_basal_cell_type_cluster_markers.txt", quote=F, sep="\t")

vitro.basal.cell.expr <- as.matrix(vitro.combined@assays$RNA@data[,c(d3.basal.cells, d21.basal.cells)])
basal.cor.res <- cor(vitro.basal.cell.expr[g1,], basal.subtype.ref)
basal.best.cl <- colnames(basal.cor.res)[apply(basal.cor.res, 1, which.max)]
basal.age.vec <- rep(c("Day3", "Day21"), c(length(d3.basal.cells), length(d21.basal.cells)))
basal.subtype.best <- cbind(basal.best.cl, basal.age.vec)
rownames(basal.subtype.best) <- rownames(basal.cor.res)
table(basal.subtype.best[which(basal.subtype.best[,2]=="Day21"),1])

basal.cm.list <- lapply(unique(basal.cm$cluster), function(x){
	intersect(basal.cm$gene_name[which(basal.cm$cluster==x)], g1)
})
names(basal.cm.list) <- paste0("Cluster", unique(basal.cm$cluster))
basal.cm.expr <- sapply(seq(length(basal.cm.list)), function(i){
	apply(t(scale(t(basal.expr.mat[basal.cm.list[[i]], ]))), 2, sum)
})
png("Plot_tSNE_fetal_basal_cm_expr.png", height=2000, width=3000)
par(mfrow=c(2,3))
for(j in seq(ncol(basal.cm.expr))){
	plotFeature2(Embeddings(basal, reduction="tsne"), values=basal.cm.expr[,j], main=names(basal.cm.list)[j], cex=2, bty="n", cex.main=3, xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
}
dev.off()

# use heatmap to visualize the fetal basal cell heterogeneity and in vitro basal cell development
vitro.basal.1 <- d3.basal.cells
e1 <- apply(vitro.basal.cell.expr[g1, vitro.basal.1], 1, mean)
vitro.basal.2 <- rownames(basal.subtype.best)[which(basal.subtype.best[,"basal.age.vec"]=="Day21" & basal.subtype.best[,"basal.best.cl"]=="Basal_Cluster5")]
e2 <- apply(vitro.basal.cell.expr[g1, vitro.basal.2], 1, mean)
vitro.basal.3 <- rownames(basal.subtype.best)[which(basal.subtype.best[,"basal.age.vec"]=="Day21" & basal.subtype.best[,"basal.best.cl"]=="Basal_Cluster2")]
e3 <- apply(vitro.basal.cell.expr[g1, vitro.basal.3], 1, mean)
basal.sub.expr <- cbind(basal.subtype.ref, e1, e2, e3)

library(gplots)
colors <- colorRampPalette(c("#2166ac","#67a9cf","#d1e5f0","#f7f7f7","#fddbc7","#ef8a62","#b2182b"), space="rgb")(50)
file.name <- "Plot_heatmap_basal_cell_heterogeneity.png"
#input <- basal.sub.expr
#mat0 <- t(scale(t(input)))
#sorted.g1 <- unlist(lapply(seq(length(basal.cm.list)), function(i){
#	genes <- basal.cm.list[[i]]
#	sorted.genes <- genes[order(mat0[genes,i], decreasing=T)]
#	return(sorted.genes)
}))
output <- mat0[sorted.g1,]
png(file.name, height=500, width=500)
par(mar=c(6,6,5,3))
heatmap.2(output, trace="none",density.info="none",scale="none", dendrogram="none", Rowv=FALSE, Colv=FALSE, col=colors,labRow=NA)
dev.off()

values=paste0("Cluster", basal@meta.data$RNA_snn_res.0.6)
basal.subtype.cols <- setNames(scales::hue_pal()(length(unique(values))), unique(values))
# Pie Chart with Percentages
slices.1 <- length(vitro.basal.1)
cl.1 <- c("Cluster5")
pct.1 <- round(slices.1/sum(slices.1)*100)
lbls.1 <- paste(cl.1, pct.1) # add percents to labels 
lbls.1 <- paste(lbls.1,"%",sep="") # ad % to labels 

slices.2 <- c(length(vitro.basal.2), length(vitro.basal.3))
cl.2 <- c("Cluster5", "Cluster2")
pct.2 <- round(slices.2/sum(slices.2)*100)
lbls.2 <- paste(cl.2, pct.2) # add percents to labels 
lbls.2 <- paste(lbls.2,"%",sep="") # ad % to labels 
pdf("Plot_piechart_vitro_basal_cell_heterogeneity.pdf", height=5, width=10)
par(mfrow=c(1,2))
pie(slices.1,labels = lbls.1, col=basal.subtype.cols[cl.1], main="Day3")
pie(slices.2,labels = lbls.2, col=basal.subtype.cols[cl.2], main="Day21")
dev.off()

# show the sample distribution of each basal cell cluster
sample.num <- sapply(sort(unique(basal@meta.data$RNA_snn_res.0.6)), function(i){
	sapply(unique(basal@meta.data$Sample), function(x){
		sum(basal@meta.data$Sample==x & basal@meta.data$RNA_snn_res.0.6==i)
	})
})
colnames(sample.num) <- paste0("Cluster", sort(unique(basal@meta.data$RNA_snn_res.0.6)))
sample.freq <- t(t(sample.num)/apply(sample.num, 2, sum))
sample.vec <- basal@meta.data$Sample
sample.col <- setNames(scales::hue_pal()(length(unique(sample.vec))), unique(sample.vec))
pdf("Plot_barplot_fetal_basal_cluster_sample_distribution.pdf")
barplot(sample.freq, col=sample.col, border=NA)
dev.off()



# clean up contaminating hub progenitors
hub.genes <- intersect(vivo.cm.mat$gene_name[which(vivo.cm.mat$cluster==4)], top.vivo.cm)
basal.genes <- intersect(vivo.cm.mat$gene_name[which(vivo.cm.mat$cluster==7)], top.vivo.cm)
# plot the hub gene enrichment scores in hub progenitor cluster and basal cell cluster 2
## define enrichment score as sum of scaled expression levels of hub cluster markers
vec0 <- apply(t(scale(t(fetal[hub.genes, which(fetal.meta$merge.cl.idx==4)]))), 2, sum)
vec1 <- apply(t(scale(t(fetal[hub.genes, which(fetal.meta$merge.cl.idx==7)]))), 2, sum)

pdf("Plot_hist_hub_gene_sum_scaled_expression_in_hub_basal_cluster.pdf", height=7, width=15)
par(mfrow=c(1,2))
hist(vec0, main="Hub progenitor cluster", xlab="Sum of scaled expression levels of hub progenitor markers")
hist(vec1, main="Basal cell cluster", xlab="Sum of scaled expression levels of hub progenitor markers")
dev.off()
cutoff <- quantile(vec0, 0.95)
names(vec1) <- colnames(fetal)[which(fetal.meta$merge.cl.idx==7)]
cell.idx <- names(which(vec1>cutoff))

dev.off()
## define enrichment score as ratio of mean of scaled expression levels of hub cluster markers versus that of basal cell markers
expr.basal.cl4 <- apply(t(scale(t(fetal[basal.genes, which(fetal.meta$merge.cl.idx==4)]))), 2, sum)
expr.basal.cl7 <- apply(t(scale(t(fetal[basal.genes, which(fetal.meta$merge.cl.idx==7)]))), 2, sum)
expr.hub.cl4 <- apply(t(scale(t(fetal[hub.genes, which(fetal.meta$merge.cl.idx==4)]))), 2, sum)
expr.hub.cl7 <- apply(t(scale(t(fetal[hub.genes, which(fetal.meta$merge.cl.idx==7)]))), 2, sum)
r.cl4 <- expr.hub.cl4/expr.basal.cl4
r.cl7 <- expr.hub.cl7/expr.basal.cl7
names(r.cl7) <- colnames(fetal)[which(fetal.meta$merge.cl.idx==7)]
cutoff <- quantile(r.cl4, 0.80)
cell.idx <- names(which(r.cl7>cutoff))

# highlight the contaminating hub cells in basal cell tSNE coordiante
plot.coor=Embeddings(basal, reduction="tsne")
png("Plot_tSNE_contaminated_hub_progenitor_cells_ratio.png", height=1000, width=1000)
plot(plot.coor, pch=16, col="gray")
points(plot.coor[cell.idx,], pch=16, col="blue")
dev.off()

# these hub signature expression cells may actually represent intermediate cell types from hub to basal cells, given that one could observe a expression gradient of hub cell features and basal cell features in cluster 2, and these potential contaminating hub cells does not form a distince clusters when performing subclustering within cluster 2
