library(maftools)
library(TCGAbiolinks)

#data_dir = "//isi-dcnl/user_data/plee/Group/public_data/tcga/"
data_dir = "Y:/Joyce/TCGA_pancancer/mutation/"
cancertype = "read"
res_folder = paste(data_dir, cancertype, "_varscan/", sep = "")
TCGAtype = "READ"

pan_type="CRC_MSI"

#dir.create(file.path(res_folder), showWarnings = FALSE)

maf.file = paste(res_folder, "TCGA.", TCGAtype, ".varscan.maf.gz", sep = "")

# TCGAbiolinks mutation data download
query <- GDCquery(
    project = "TCGA-", #download specific cancer type 
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    access = "open"
)
GDCdownload(query)
maf.file <- GDCprepare(query)


mafObj = read.maf(maf = maf.file)

sample_total = getSampleSummary(mafObj)

sample_reid=sample_total
sample_reid$reid = substr(sample_reid$Tumor_Sample_Barcode, 1,12)
sample_reid$reid=gsub("-",".",sample_reid$reid)

all_mut_gene=mafObj@data
tex_sig_mut=subset(sub_cluster, cancertype %in% pan_type)

sub_sample=subset(sample_reid, reid %in% tex_sig_mut$pid)
sub_sample$cluster= tex_sig_mut$cluster[match(sub_sample$reid, tex_sig_mut$pid)] 
#tex_sig_mut=subset(tex_sig_mut, X %in% sample_reid$reid)


clst1=subset(sub_sample, cluster == 1)
clst2=subset(sub_sample, cluster == 2)
clst3=subset(sub_sample, cluster == 3)
clst4=subset(sub_sample, cluster == 4)
clst5=subset(sub_sample, cluster == 5)
clst0=subset(sub_sample, cluster == 0)


clst0_mafobj=subsetMaf(mafObj, tsb = clst0$Tumor_Sample_Barcode)

clst1_mafobj=subsetMaf(mafObj, tsb = clst1$Tumor_Sample_Barcode)
clst2_mafobj=subsetMaf(mafObj, tsb = clst2$Tumor_Sample_Barcode)
clst3_mafobj=subsetMaf(mafObj, tsb = clst3$Tumor_Sample_Barcode)
clst4_mafobj=subsetMaf(mafObj, tsb = clst4$Tumor_Sample_Barcode)
clst5_mafobj=subsetMaf(mafObj, tsb = clst5$Tumor_Sample_Barcode)

#clst2.vs.5 <- mafCompare(m1 = clst2_mafobj, m2 = clst5_mafobj, m1Name = 'cluster2', m2Name = 'cluster5', minMut = 5)
#print(clst2.vs.5)
#forestPlot(mafCompareRes = clst2.vs.5, pVal = 0.05, color = c('maroon','royalblue'))

#coOncoplot(m1 = clst3_mafobj, m2 = clst5_mafobj, m1Name = 'cluster3', m2Name = 'cluster5', genes = genes, removeNonMutated = TRUE)
#coBarplot(m1 = clst2_mafobj, m2 = clst5_mafobj, m1Name = "cluster2", m2Name = "cluster5")

#drugInteractions(maf = clst4_mafobj, fontSize = 0.75)
#drugInteractions(maf = clst5_mafobj, fontSize = 0.75)

#OncogenicPathways(maf = clst4_mafobj)
#OncogenicPathways(maf = clst5_mafobj)

#PlotOncogenicPathways(maf = clst4_mafobj, pathways = "RTK-RAS")
#PlotOncogenicPathways(maf = clst5_mafobj, pathways = "RTK-RAS")


#sample_reid$cluster= tex_sig_mut$cluster[match(sample_reid$reid, tex_sig_mut$pid)]
mafobjs=NULL

mafobjs = vector(mode = "list", length = 6)
names(mafobjs) = c('cluster0','cluster1','cluster2','cluster3','cluster4','cluster5')
names(mafobjs) = c('cluster0','cluster1','cluster2','cluster3','cluster5')

mafobjs[1]=clst0_mafobj
mafobjs[2]=clst1_mafobj
mafobjs[3]=clst2_mafobj
mafobjs[4]=clst3_mafobj
mafobjs[5]=clst4_mafobj
mafobjs[6]=clst5_mafobj


#pathdb <- system.file("extdata", "oncogenic_sig_patwhays.tsv", package = "maftools")
#pathdb = data.table::fread(input = pathdb)
#pathdb_size = pathdb[,.N,Pathway]
#pathdb = split(pathdb, as.factor(pathdb$Pathway))


pwy_aff = pathdb_size[,1]
samp_aff = pathdb_size[,1]

for (icmaf in 1:length(mafobjs)) {
  cluster = names(mafobjs)[icmaf]
  maf = mafobjs[[icmaf]]
  altered_pws = lapply(pathdb, function(pw){
    x = suppressMessages(try(genesToBarcodes(maf = maf, genes = pw$Gene)))
    if(class(x) == "try-error"){
      pw_genes = NULL
    }else{
      pw_genes = names(genesToBarcodes(maf = maf, genes = pw$Gene, justNames = TRUE, verbose = FALSE))
    }
    pw_genes
  })
  
  mut_load = lapply(altered_pws, function(x){
    if(is.null(x)){
      nsamps =  0
    }else{
      nsamps = length(unique(as.character(unlist(
        genesToBarcodes(
          maf = maf,
          genes = x,
          justNames = TRUE
        )
      ))))
    }
    nsamps
  })
  
  altered_pws = as.data.frame(t(data.frame(lapply(altered_pws, length))))
  data.table::setDT(x = altered_pws, keep.rownames = TRUE)
  colnames(altered_pws) = c("Pathway", "n_affected_genes")
  altered_pws$Pathway = gsub(pattern = "\\.", replacement = "-", x = altered_pws$Pathway)
  altered_pws = merge(pathdb_size, altered_pws, all.x = TRUE)
  
  altered_pws[, fraction_affected := n_affected_genes/N]
  altered_pws$Mutated_samples = unlist(mut_load)
  nsamps = as.numeric(maf@summary[ID == "Samples", summary])
  altered_pws[,Fraction_mutated_samples := Mutated_samples/nsamps]
  altered_pws = altered_pws[order(n_affected_genes, fraction_affected, decreasing = FALSE)]
  
  altered_pws = altered_pws[!n_affected_genes %in% 0]
  
  #altered_pws=as.data.frame(altered_pws)
  #rownames(altered_pws)=altered_pws$Pathway
  #altered_pws$Pathway=NULL
  
  
  
  #colnames(altered_pws)=paste(colnames(altered_pws), "_cluster",icmaf-1, sep = '')
  pwy_aff[,cluster]=altered_pws$fraction_affected[match(pwy_aff$Pathway, altered_pws$Pathway)]
  samp_aff[,cluster]=altered_pws$Fraction_mutated_samples[match(samp_aff$Pathway, altered_pws$Pathway)]
  
  
}

write.csv(samp_aff, paste(res_dir, pan_type,'_sample_aff_across_cluster.csv', sep = ''))

write.csv(pwy_aff, paste(res_dir, pan_type, '_pwy_aff_across_cluster.csv', sep = ''))


plot_pwy<-gather(pwy_aff, "cluster", "Fraction", -Pathway)
plot_samp<-gather(samp_aff, "cluster", "Fraction", -Pathway)

g=ggplot(plot_samp, aes(fill=cluster, y=Fraction, x=Pathway)) + 
  geom_bar(position= position_fill(reverse = TRUE), stat="identity")+coord_flip()+scale_fill_manual(values = c("royalblue4","cornflowerblue","#9FCAE6","lightsalmon1", "orangered2","firebrick4"))
ggsave(plot = g, filename = paste(res_dir, pan_type,"crcall_sample_affected_pwy_cluster.png", sep = ""),  dpi=300,width = 7, heigh = 6)

c=ggplot(plot_pwy, aes(fill=cluster, y=Fraction, x=Pathway)) + 
  geom_bar(position= position_fill(reverse = TRUE),stat="identity")+coord_flip()+scale_fill_manual(values = c("royalblue4","cornflowerblue","#9FCAE6","lightsalmon1", "orangered2","firebrick4"))
ggsave(plot = c, filename = paste(res_dir, pan_type, "crcall_pwy_affected_pwy_cluster.png", sep = ""),  dpi=300,width = 7, heigh = 6)


OncogenicPathways = function(maf, pathways = NULL, fontSize = 1, panelWidths = c(2, 4, 4)){
  if(is.null(pathways)){
    pathdb <- system.file("extdata", "oncogenic_sig_patwhays.tsv", package = "maftools")
    pathdb = data.table::fread(input = pathdb)
  }else{
    pathdb = data.table::copy(x = pathways)
    colnames(pathdb) = c("Gene", "Pathway")
    data.table::setDT(x = pathdb)
  }
  pathdb_size = pathdb[,.N,Pathway]
  pathdb = split(pathdb, as.factor(pathdb$Pathway))
  
  
  altered_pws = lapply(pathdb, function(pw){
    x = suppressMessages(try(genesToBarcodes(maf = maf, genes = pw$Gene)))
    if(class(x) == "try-error"){
      pw_genes = NULL
    }else{
      pw_genes = names(genesToBarcodes(maf = maf, genes = pw$Gene, justNames = TRUE, verbose = FALSE))
    }
    pw_genes
  })
  
  mut_load = lapply(altered_pws, function(x){
    if(is.null(x)){
      nsamps =  0
    }else{
      nsamps = length(unique(as.character(unlist(
        genesToBarcodes(
          maf = maf,
          genes = x,
          justNames = TRUE
        )
      ))))
    }
    nsamps
  })
  
  altered_pws = as.data.frame(t(data.frame(lapply(altered_pws, length))))
  data.table::setDT(x = altered_pws, keep.rownames = TRUE)
  colnames(altered_pws) = c("Pathway", "n_affected_genes")
  altered_pws$Pathway = gsub(pattern = "\\.", replacement = "-", x = altered_pws$Pathway)
  altered_pws = merge(pathdb_size, altered_pws, all.x = TRUE)
  
  altered_pws[, fraction_affected := n_affected_genes/N]
  altered_pws$Mutated_samples = unlist(mut_load)
  nsamps = as.numeric(maf@summary[ID == "Samples", summary])
  altered_pws[,Fraction_mutated_samples := Mutated_samples/nsamps]
  altered_pws = altered_pws[order(n_affected_genes, fraction_affected, decreasing = FALSE)]
  
  altered_pws = altered_pws[!n_affected_genes %in% 0]
  
  if(nrow(altered_pws) == 0){
    stop("None of the provided pathways are altered!")
  }
  
  lo = layout(mat = matrix(c(1, 2, 3), ncol = 3, nrow = 1), widths = panelWidths)
  par(mar = c(2, 2, 2, 0))
  plot(NA, xlim = c(0, 1), ylim = c(0, nrow(altered_pws)), axes = FALSE)
  text(x = 1, y = seq(0.5, nrow(altered_pws), by = 1), labels = altered_pws$Pathway, adj = 1, xpd = TRUE, cex = fontSize)
  title(main = "Pathway", adj = 0)
  
  par(mar = c(2, 0, 2, 1), xpd = TRUE)
  plot(NA, xlim = c(0, 1.2), ylim = c(0, nrow(altered_pws)), axes = FALSE)
  rect(xleft = 0, ybottom = seq(0.1, nrow(altered_pws), by = 1), xright = 1, ytop = seq(0.2, nrow(altered_pws), by = 1)+0.7, col = '#ecf0f1', border = "white")
  rect(xleft = 0, ybottom = seq(0.1, nrow(altered_pws), by = 1), xright = altered_pws$fraction_affected, ytop = seq(0.2, nrow(altered_pws), by = 1)+0.7, col = "#c0392b", border = "white")
  text(x = 1.05, y = seq(0.5, nrow(altered_pws), by = 1), labels = paste0(altered_pws$n_affected_genes, "/", altered_pws$N), adj = 0, cex = fontSize)
  axis(side = 1, at = seq(0, 1, 0.25), line = -1, cex.axis = fontSize)
  title(main = "Fraction of pathway affected", adj = 0)
  
  par(mar = c(2, 0, 2, 1), xpd = TRUE)
  plot(NA, xlim = c(0, 1.2), ylim = c(0, nrow(altered_pws)), axes = FALSE)
  rect(xleft = 0, ybottom = seq(0.1, nrow(altered_pws), by = 1), xright = 1, ytop = seq(0.2, nrow(altered_pws), by = 1)+0.7, col = '#ecf0f1', border = "white")
  rect(xleft = 0, ybottom = seq(0.1, nrow(altered_pws), by = 1), xright = altered_pws$Fraction_mutated_samples, ytop = seq(0.2, nrow(altered_pws), by = 1)+0.7, col = "#c0392b", border = "white")
  text(x = 1.05, y = seq(0.5, nrow(altered_pws), by = 1), labels = paste0(altered_pws$Mutated_samples, "/", nsamps), adj = 0, cex = fontSize)
  axis(side = 1, at = seq(0, 1, 0.25), line = -1, cex.axis = fontSize)
  title(main = "Fraction of samples affected")
  
  altered_pws
}




rownames(expr) = expr[,1]
ssdata = expr[, 3:dim(expr)[2]]
ssdata = as.data.frame(t(ssdata))
cat("Expression data dimension, sample: probe ", dim(ssdata), "\n")
print(ssdata[1:9,1:6])

ssannot = expr[,1:2]
cat("Annotation data dimension, probe: annot ", dim(ssannot), "\n")
print(head(ssannot))

sign_gene = toupper(rownames(sign))  
cat("Signature gene number: ", length(sign_gene),"\n")
anno_gene = rownames(expr)
ol_gene = intersect(sign_gene, anno_gene)
unava_gene = setdiff(sign_gene, ol_gene)
cat("Unavailable",length(unava_gene) ," genes: ", unava_gene, "\n")
sub_annot = ssannot[ssannot$Hugo_Symbol %in% sign_gene,]
sub_annot=as.data.frame(sub_annot)
cat("Available probe: ", dim(sub_annot)[1], "\n")
sssign = as.data.frame(matrix(ncol = 3, nrow = dim(sub_annot)[1])) # Initialize sig.score x
colnames(sssign) = c("probe", "EntrezGene.ID", "coefficient")
rownames(sssign) = rownames(sub_annot)
sssign$probe = rownames(sub_annot)
sssign$EntrezGene.ID = 1:18#sub_annot$Entrez_Gene_Id

sssign$coefficient = 1#sign[sssign$probe,"logfc"]
cat("Run sig.score in all patient samples\n")
sig_score = sig.score(x=sssign, data=ssdata, annot=ssannot, do.mapping=FALSE, signed=TRUE, verbose=TRUE)
# print(head(sssign))
}
ssin = list("data" = ssdata, "annot" = ssannot, "x" = sssign)
if (sig_save) {
  st = Sys.time()
  cat("Saving sig.score input and output to ", res_dir, "\n")
  saveRDS(ssin, file = paste(res_dir, "sig_score_inputs.RDS", sep = ""))
  saveRDS(sig_score, file = paste(res_dir, "sig_score_outputs.RDS", sep = ""))
  print(Sys.time()-st)
}








