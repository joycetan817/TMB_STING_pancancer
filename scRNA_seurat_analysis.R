rm(list = ls())


sr_QC = function (sc_obj, resdir, proFlag, patient, tissue) {
	sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^MT-")
	sc_obj@meta.data$percent.mt[is.na(sc_obj@meta.data$percent.mt)] <- 0
	featplot <- VlnPlot(sc_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	#if (class(featplot) == "try-error") {tmp_mt <- sc_obj@meta.data[["percent.mt"]]
		#tmp_mt[is.na(tmp_mt)] = 0
		#sc_obj@meta.data[["percent.mt"]] <- tmp_mt
		#featplot <- VlnPlot(sc_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
		#}
	ggsave(plot=featplot, file = paste(resdir, patient, tissue, "_feature_plot_bf_filter.png", sep = ""),  width = 9, height = 6)
	plot1 <- FeatureScatter(sc_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
	plot2 <- FeatureScatter(sc_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
	featscat <- CombinePlots(plots = list(plot1, plot2))
	ggsave(plot = featscat, file = paste(resdir, patient, tissue, "_feature_scatter_bf_filter.png", sep = ""),  width = 9, height = 6)


	sc_obj <- subset(sc_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20) ###check feature plot
	featplot_filt <- VlnPlot(sc_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	ggsave(plot = featplot_filt, file = paste(resdir, patient, tissue, "_feature_plot_aft_filter.png", sep = ""),  width = 9, height = 6)

	####Normalizing the data
	sc_obj <- NormalizeData(sc_obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
	#if (proFlag) {sc_obj = tryCatch(NormalizeData(sc_obj, assay = "Protein", normalization.method = "CLR"), error=function(e) sc_obj)}
	if (proFlag) {sc_obj = NormalizeData(sc_obj, assay = "Protein", normalization.method = "CLR")}

	sc_obj <- FindVariableFeatures(sc_obj, selection.method = "vst", nfeatures = 2000)	
	top10 <- head(VariableFeatures(sc_obj), 10)
	# plot variable features with and without labels
	plot1 <- VariableFeaturePlot(sc_obj)
	plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
	ggsave(plot = plot2, filename = paste(resdir, patient, tissue, "_vari_feature.png", sep = ""),  width = 7, height = 6)
	return (sc_obj)

}

sr_clst = function(sc_obj, resdir, proFlag) {

	#sc_obj <- FindVariableFeatures(sc_obj, selection.method = "vst", nfeatures = 2000)	
	#top10 <- head(VariableFeatures(sc_obj), 10)
	 #plot variable features with and without labels
	#plot1 <- VariableFeaturePlot(sc_obj)
	#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
	#ggsave(plot = plot2, filename = paste(res_dir, "vari_feature.png", sep = ""),  width = 7, heigh = 6)

	

	#allGenes = rownames(sc_obj[["RNA"]])
	sc_obj = ScaleData(sc_obj) #, features = allGenes, assay = "RNA")
	if (proFlag) {sc_obj = ScaleData(sc_obj, assay = "Protein")}
	
	sc_obj <- RunPCA(sc_obj) #, features = VariableFeatures(object = sc_obj))
	pca_vizdim = VizDimLoadings(sc_obj, dim = 1:6, reduction = "pca")
	ggsave(plot = pca_vizdim, filename = paste(resdir, "pca_vizdim_plot.png", sep = ""), width = 10, height = 15)
	
	pca_dim = DimPlot(sc_obj, reduction = "pca")
	ggsave(plot = pca_dim, filename = paste(resdir, "pca_dim_plot.png", sep = ""), width = 9, height = 6)

	png(paste(resdir, "pca_heatmap.png", sep = ""), res = 180, width = 12, height = 16, units = "in")
	DimHeatmap(sc_obj, dims = 1:15, cells = 500, balanced = TRUE)
	gar = dev.off()

	#sc_obj = JackStraw(sc_obj, num.replicate = 100)
	#sc_obj = ScoreJackStraw(sc_obj, dims = 1:20)
	#js_plot = JackStrawPlot(sc_obj, dims = 1:20)
	#ggsave(plot = js_plot, filename = paste(resdir, "JackStraw_plot.png", sep = ""), width = 9, height = 6)

	elb_plot = ElbowPlot(sc_obj, ndims = 50)
	ggsave(plot = elb_plot, filename = paste(resdir, "elbow_plot.png", sep = ""), width = 9, height = 6)

	saveRDS(sc_obj, paste(resdir, "seurat_obj_cor_bf_pc.RDS", sep = ""))

	return (sc_obj)
}

sr_plot = function (sc_obj, ndim, resdir, proFlag) {
	sc_obj = FindNeighbors(sc_obj, dims = 1:ndim)
	sc_obj = FindClusters(sc_obj, resolution = 0.5, algorithm = 4, random.seed = 256)

	sc_obj = RunUMAP(sc_obj, dims = 1:ndim)
	cls_umap = DimPlot(sc_obj, reduction = "umap", label=TRUE, pt.size=0.5)

	ggsave(plot = cls_umap, filename = paste(resdir, ndim, "_cluster_umap_plot.png", sep = ""), width = 9, height = 6)

	sc_obj = RunTSNE(sc_obj, dims = 1:ndim)
	cls_tsne = DimPlot(sc_obj, reduction = "tsne")
	ggsave(plot = cls_tsne, filename = paste(resdir, ndim, "_cluster_tsne_plot.png", sep = ""), width = 9, heigh = 6)

	btc = DimPlot(sc_obj, reduction = "umap", group.by = "orig.ident")
	ggsave(plot = btc, filename = paste(resdir, "pca", ndim, "_batch_check_aft_inte.png", sep = ""), width = 9, height = 6)

	#tisSource = DimPlot(sc_obj, reduction = "umap", group.by = "mtis")
	#ggsave(plot = tisSource, filename = paste(resdir, "pca", ndim, "_tissue_source_aft_inte.png", sep = ""), width = 9, height = 6)

	RNAexpr = c("PTPRC", "EPCAM", "ATL1", "CD3D", "CD8A", "rna_CD4", "PDCD1", "TCF7", "NKG7", "MS4A1", "CD27", "CD19", 
	             "CD68", "FCER1A", "CD14", "FCGR3A", "ITGAX", "IL3RA", "NCAM1", "TPSAB1", "GNLY", "FOXP3", "TBX21", "GATA3")

	#expr = c("PTPRC", "CD3E", "CD4", "CD8A", "TCF7", "IL7R", "ITGAE", "PDCD1", "ENTPD1", "DPP4", "FOXP3")
	expr_ol = FeaturePlot(sc_obj, features = RNAexpr, reduction = "umap", repel = TRUE,
					    min.cutoff = "q05", max.cutoff = "q95", ncol = 3, cols = c("orange1", "grey", "purple"))
	ggsave(plot = expr_ol, filename = paste(resdir, "pca", ndim ,"_expression_overlay.png", sep = ""),  width = 12, height = 22)


	#expr_ol = FeaturePlot(sc_obj, features = c("ACTA2", "MCAM","THY1", "NT5E", "CD44", "PROM1", "HLA-DRA", "CD34"), reduction = "umap", 
	#min.cutoff = "q05", max.cutoff = "q95",  cols = c("orange1", "grey", "purple"))
	#ggsave(plot = expr_ol, filename = paste(resdir, "pca", ndim ,"stromal_expression_overlay.png", sep = ""),  width = 16, height = 12)

	if (resdir == reclst_res_dir) {
		expr = c( "CD3E","rna_CD4", "CD8A", "TCF7", "IL7R", "CCR7", "FOXP3","DPP4","ITGAE", "KLRG1","CXCL13","CD69","PRF1","GZMB", "PDCD1", "ENTPD1", "HAVCR2","MKI67")
		expr_ol = FeaturePlot(sc_obj, features = expr, reduction = "umap",
		                      min.cutoff = "q05", max.cutoff = "q95", ncol = 3, cols = c("orange1", "grey", "purple"))
		ggsave(plot = expr_ol, filename = paste(resdir, "tcell_expression_overlay.png", sep = ""),  width = 9, heigh = 12)
		expr = c( "CD3E","rna_CD4", "CD8A", "TCF7", "IL7R", "CCR7","KLRG1","CD69","GZMB", "PDCD1", "ENTPD1", "HAVCR2" "CXCL13","SLAMF6",'CX3CR1',"MKI67")
		expr_ol = FeaturePlot(sc_obj, features = expr, reduction = "umap",
		                      min.cutoff = "q05", max.cutoff = "q95", ncol = 3, cols = c("orange1", "grey", "purple"))
		#ggsave(plot = expr_ol, filename = paste(resdir, "cd4_tcell_expression_overlay.png", sep = ""),  width = 9, heigh = 12)
	}

	RNAexpr = c("PTPRC", "EPCAM", "ATL1", "CD3D", "CD8A", "rna_CD4", "PDCD1", "TCF7", "NKG7", "MS4A1", "CD27", "CD19", 
		             "CD68", "FCER1A", "CD14",  "FCGR3A","ITGAX", "IL3RA", "NCAM1", "TPSAB1", "GNLY", "FOXP3", "TBX21", "GATA3")

	vln_feat = VlnPlot(sc_obj, features = RNAexpr, slot = "counts", log = TRUE, ncol = 3, pt.size = 0.5)
	ggsave(plot = vln_feat, filename = paste(resdir, "pca", ndim, "_violin_gene_feature.png", sep = ""), width = 15, height = 20)

 	#bcell_marker<-c("PTPRC","MS4A1","CD19","CD27","CD38","SDC1","CXCR4,"LRRC32","NT5E","ENTPD1","IL10","CD274","IGHD")

 	#bcell_marker<-c("PTPRC","MS4A1","CD19","CD27","CD38","SDC1","CXCR4","IGHD")
 	#dc_macro_marker = c("PTPRC","CD4", "HLA-DRA", "CD68","CCR1","CD14","LAMP3","BATF3", "CD1C", "PTCRA","CLEC4C","FCER1A", "KIT","TPSAB1","ITGAX","IL3RA", "THBD","CLEC9A","AXL")
 	#expr_ol = FeaturePlot(sc_obj, features = dc_macro_marker, reduction = "umap", min.cutoff = "q05", max.cutoff = "q95",  ncol=3 ,cols = c("orange1", "grey", "purple"))
	#ggsave(plot = expr_ol, filename = paste(resdir, "pca", ndim ,"macro_expression_overlay.png", sep = ""),  width = 12, height = 18)

 	#plasma_marker = c('SDC1','CD27','CXCR4','CD38')

 	cd26_marker=c("TBX21","GATA3","EOMES","STAT1","GATA1","RORA","RORC","IRF4","RORB","STAT3","CEBPA","ELK3","RUNX2")

	if (proFlag) {
		pro_expr = paste0("Protein_", rownames(sc_obj[["Protein"]]))
		prot_ol = FeaturePlot(sc_obj, features = pro_expr, reduction = "umap",
					    min.cutoff = "q05", max.cutoff = "q95", ncol = 3, cols = c("orange1", "grey", "purple"))
		ggsave(plot = prot_ol, filename = paste(resdir, "pca", ndim, "_TotalSeq_overlay.png", sep = ""), width = 10, height = 16)

		ridge_plot = RidgePlot(sc_obj, features = pro_expr, ncol = 4)
		ggsave(plot = ridge_plot, filename = paste(resdir, "pca", ndim, "_RidgePlot_clusters.png", sep = ""), width = 16, height = 16)
			}
	cluster.markers <- FindAllMarkers(sc_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
	write.csv(cluster.markers, paste(resdir, "pca", ndim, "_cluster_marker.csv", sep = ""))

	top.markers = cluster.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
	write.csv(top.markers, file = paste(resdir, "top_pos_cluster_markers.csv", sep = ""))

	png(paste(resdir, "marker_heatmap.png", sep = ""),  res= 180, width = 16, height = 12, units = "in")
	print(DoHeatmap(object = sc_obj, features = top.markers$gene) + NoLegend())
	gar = dev.off()
	saveRDS(sc_obj, paste(resdir, "seurat_obj_cor_aft_pca", ndim, ".RDS", sep = ""))
	
	return (sc_obj)
}


read10XSelf = function(data_dir) {
	allFiles = list.files(data_dir)
	#if (length(allFiles) != 3) {stop("Too many/few input 10x files!")}
	if (all(grepl(".gz" ,allFiles))) {
		cat("CellRange v3 output is detected!\n")
		barcodeDf = read.table(gzfile(paste(data_dir, "/barcodes.tsv.gz", sep = "")), stringsAsFactors = FALSE)
		colnames(barcodeDf) = c("CellID")
		barcodeDf$CellID = str_replace(barcodeDf$CellID, "-1", "")
		featureDf = read.table(gzfile(paste(data_dir, "/features.tsv.gz", sep = "")), 
				       #sep = "\t", stringsAsFactors = FALSE)
		#featureDf = read.table(gzfile(paste(data_dir, "/genes.tsv.gz", sep = "")), 
				       sep = "\t", stringsAsFactors = FALSE)
		colnames(featureDf) = c("ENSG", "GeneSymbol", "Type")
		matrixDf = read.table(gzfile(paste(data_dir, "/matrix.mtx.gz", sep = "")), skip = 2)
		colnames(matrixDf) = c("GeneID", "CellID", "Count")

	} else if (any(grepl("gene", allFiles))) {
		cat("CellRange v2 output is detected!\n")
		barcodeDf = read.table(paste(data_dir, "/barcodes.tsv", sep = ""), stringsAsFactors = FALSE)
		colnames(barcodeDf) = c("CellID")
		barcodeDf$CellID = str_replace(barcodeDf$CellID, "-1", "")
		featureDf = read.table(paste(data_dir, "/genes.tsv", sep = ""), sep = "\t", stringsAsFactors = FALSE)
		colnames(featureDf) = c("ENSG", "GeneSymbol")
		matrixDf = read.table(paste(data_dir, "/matrix.mtx", sep = ""), skip = 2)
		colnames(matrixDf) = c("GeneID", "CellID", "Count")
	}
	cat("\tCell number:", dim(barcodeDf)[1], "\n")
	cat("\tGene number:", dim(featureDf)[1], "\n")
	cat("\tTotal available counts:", dim(matrixDf)[1], "\n")
	matDf = as.data.frame(matrixDf)
	matDf = matDf[complete.cases(matDf[,"Count"]),]
	sprMatDf = spread(matDf, "CellID", "Count")
	rownames(sprMatDf) = sprMatDf$GeneID
	sprMatDf$GeneID = NULL

	colnames(sprMatDf) = barcodeDf[colnames(sprMatDf), "CellID"]
	mergeMatDf = merge(featureDf, sprMatDf, by = 0)
	print(dim(mergeMatDf))
	mergeMatDf = mergeMatDf[complete.cases(mergeMatDf[,c("ENSG", "GeneSymbol")]),]
	cleanMatDf = mergeMatDf
	rownames(cleanMatDf) = make.unique(names = cleanMatDf$GeneSymbol)
	cleanFeatDf = cleanMatDf[,1:(ncol(featureDf)+1)]
	cleanMatDf[,1:(ncol(featureDf)+1)] = NULL
	cleanMatDf[is.na(cleanMatDf)] = 0
	print(cleanMatDf[1:9,1:6])
	return(cleanMatDf)
}

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(stringr))
suppressMessages(library(scales))
suppressMessages(library(monocle))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))


work_dir = "..."
sample_type = "10X"#"10X" #"public_br_scRNAseq" #Public scRNA raw data
data_name = "..." #output folder
cluster = "immune_res" # all immune cells clustering
reclst = "tcell_res" #recluster all CD3+ T cells


data_dir = paste(work_dir, sample_type, "/", sep = "") 
res_dir = paste(work_dir, data_name, "/", cluster, "/", sep = "")
reclst_res_dir = paste(work_dir, data_name, "/", reclst, "/", sep = "")



#batch_check = TRUE #if batch check is needed
#batch_cor = TRUE #if batch correction is needed
recluster = TRUE
proteinflag = TRUE
datatype = "filtered" #filtered"

if (datatype == "raw") {
	inFiles = list.files(data_dir, pattern = "_raw_feature_bc_matrix")
	saName = str_split_fixed(inFiles, "_raw", n = 2)[,1]
}
if (datatype == "filtered") {
	#inFiles = list.files(data_dir, pattern = "_filtered_feature_bc_matrix")
	inFiles = list.files(data_dir, pattern = "_filtered_gene_bc_matrix")
	saName = str_split_fixed(inFiles, "_fil", n = 2)[,1]
}
print(saName)

srData = vector(mode = "list", length = length(saName))
srObjs = vector(mode = "list", length = length(saName))
names(srData) = saName
names(srObjs) = saName

for (isa in 1:length(inFiles)) {
	st = Sys.time()
	message(saName[isa], ":" , inFiles[isa])
	tmpDir = paste(data_dir, inFiles[isa], sep = "")
	tmpData = try(Read10X(tmpDir))
	if (class(tmpData) == "try-error") {
			cat("Using customized function to read 10X result!\n")
			tmpData = read10XSelf(tmpDir)}

	#srData[saName[isa]] = tmpData

	tisName = str_split_fixed(saName[isa], "_", n = 3)[,3]
	patName = str_split_fixed(saName[isa], "_", n = 3)[,1]
	primeName = str_split_fixed(saName[isa], "_", n = 3)[,2]
	
	if (proteinflag) { cat("Protein data is available!\n")
		if (primeName == "five") {
			if(patName == "BC394"){
				oldRowName = rownames(tmpData$`Antibody Capture`)
				newRowName = str_split_fixed(oldRowName, "_", n = 4)[,3]
				newRowName = str_split_fixed(newRowName, "-", n = 2)[,1]
			
			}

			if (patName == "BC401") {
				oldRowName = rownames(tmpData$`Antibody Capture`)
				newRowName = str_split_fixed(oldRowName, " ", n = 4)[,3]			
				newRowName[18]<-paste(newRowName[18], "3", sep = "")
				newRowName[19]<-paste(newRowName[19], "1", sep = "")

			}
			rownames(tmpData$`Antibody Capture`) = newRowName
			tmpObj = CreateSeuratObject(counts = tmpData$`Gene Expression`, project = saName[isa], min.cells = 3, min.features = 5)
			tmpObj[["Protein"]] = CreateAssayObject(counts = tmpData$`Antibody Capture`)
		} else {
		tmpObj = CreateSeuratObject(counts = tmpData, project = saName[isa],  min.cells = 3, min.features = 5)
		}
	}
	cat(tisName, patName, "\n")
	tmpObj@meta.data$tissue = tisName
	tmpObj@meta.data$patient = patName
	tmpObj@meta.data$prime = primeName


	#srObjs[saName[isa]] = tmpObj
	#srObjs[isa] = tmpObj

	#print(tmpObj)
	tmpObj = sr_QC(tmpObj, resdir = res_dir, proFlag = proteinflag, patient = patName, tissue = tisName)
	srObjs[[isa]] = tmpObj
	assign(saName[isa], tmpObj)
	
}

saveRDS(srObjs, paste(res_dir, "seurat_obj_bf_integrate.RDS"))

#rPCA integration
features <- SelectIntegrationFeatures(object.list = srObjs)
srObjs <- lapply(X = srObjs, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

st = Sys.time()
anchors <- FindIntegrationAnchors(object.list = srObjs,  reduction = "rpca",
    dims = 1:50)
inteObjs <- IntegrateData(anchorset = anchors, dims = 1:50)
#inteObjs <- ScaleData(inteObjs, verbose = FALSE)
#inteObjs <- RunPCA(inteObjs, verbose = FALSE)
#inteObjs <- RunUMAP(inteObjs, dims = 1:50)
#DimPlot(inteObjs, group.by = "orig.ident")
saveRDS(inteObjs, 'seurat_obj_allsamp_aft_integrate_rpca_correct.RDS')

print(Sys.time()-st)


#classic integration
integrate.anchors <- FindIntegrationAnchors(object.list = srObjs, dims = 1:30)
inteObjs <- IntegrateData(anchorset = integrate.anchors, dims = 1:30)
DefaultAssay(inteObjs) <- "integrated"
saveRDS(inteObjs, 'seurat_obj_alltumor_aft_integrate.RDS')


inteObjs = sr_clst(inteObjs, resdir = res_dir, proFlag = proteinflag)

inteObjs = sr_plot(inteObjs, ndim = 15, resdir = res_dir, proFlag = proteinflag)

#saveRDS(inteObjs, paste(res_dir, "seurat_obj_af_pc.RDS", sep=""))

new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
    "NK", "DC", "Platelet")
clst_id <- new.cluster.ids
names(clst_id) <- levels(inteObjs)
annot_obj <- RenameIdents(inteObjs, clst_id)
annot_plot<-DimPlot(annot_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(plot = annot_plot, filename = paste(res_dir, "pca", ndim, "_cluster_annot.png", sep = ""),  width = 9, heigh = 6)

ic <- c(1,2,3,4)  #interested clusters for reclustering

if (recluster) {

	sub_clst <- subset(annot_obj, ident = new.cluster.ids[ic]) ###recluster T cell based on CD3 and CD45
	sub_check <- DimPlot(sub_clst)
	ggsave(plot = sub_check, filename = paste(reclst_res_dir, "subset_check.png", sep = ""),  width = 9, heigh = 6)
	bt_check <- DimPlot(sub_clst, group.by = "orig.ident")
	ggsave(plot = bt_check, filename = paste(reclst_res_dir, "subset_batch_check.png", sep = ""),  width = 9, heigh = 6)
	saveRDS(sub_clst, paste(reclst_res_dir, "seurat_obj_subset_cluster.RDS", sep = ""))

	subObj <- CreateSeuratObject(counts = sub_clst[["RNA"]]@counts, meta.data = sub_clst@meta.data, project = "Tcell") # Reset the Seurat Objects, tmpSubSrsc[["RNA"]]@counts is the merged matrix
	if (proteinflag) {
		subObj[["Protein"]] = CreateAssayObject(counts = sub_clst[["Protein"]]@counts)
	}

	subObj <- SplitObject(subObj, split.by = "orig.ident") # Split by patients or samples

	#check cell number for each sample to determine the k.weight for data integration and filter out the sample(s) with low cell number
	cellN=as.data.frame(matrix(nrow = length(subObj), ncol=1))
	rownames(cellN)=names(subObj)
	colnames(cellN)='#cell'

	for (i in 1:length(subObj)) {
		sub_obj=subObj[[i]]
		cellN[i, 1]=dim(sub_obj)[2]

	}
	write.csv(cellN, paste(reclst_res_dir,'sample_cellN_table.csv', sep=''))
	rm_obj=rownames(cellN)[cellN$'#cell' <= 30]
	subObj[rm_obj]=NULL

	for (i in 1:length(subObj)) {
	    subObj[[i]] <- NormalizeData(subObj[[i]], verbose = FALSE)
	    if (proteinflag) {subObj[[i]] = NormalizeData(subObj[[i]], assay = "Protein", normalization.method = "CLR")}
	    subObj[[i]] <- FindVariableFeatures(subObj[[i]], selection.method = "vst", 
	        nfeatures = 2000, verbose = FALSE)
	}

	st=Sys.time()
	integrate.anchors <- FindIntegrationAnchors(object.list = subObj, dims = 1:30)
	sub_inteObj <- IntegrateData(anchorset = integrate.anchors, dims = 1:30, k.weight = 66)
	DefaultAssay(sub_inteObj) <- "integrated"
	saveRDS(sub_inteObj,  "seurat_obj_sub_clst_af_inte_cd4_reclst.RDS")
	print(Sys.time()-st)

	sub_inteObj = sr_clst(sub_inteObj, resdir = reclst_res_dir, proFlag = proteinflag)
	sub_inteObj = sr_plot(sub_inteObj, ndim = 19, resdir = reclst_res_dir, proFlag = proteinflag)