#conda activate vitessce_jiwoo4

# Load data
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
options(future.globals.maxSize = 100000 * 1024^2)
library(gghighlight)
library(RColorBrewer)
library(reticulate)
library(magrittr)
library(EnvStats)
library(matrixStats)
library(InSituType)
library(UCell)
library(stringr)
library(SCINA)
library(hoodscanR)
library(SpatialExperiment)
library(CellChat)
library(viridis)
library(Matrix)
source("/home/singlecell/Rfunction/splitDimPlot2.R")
source("/home/singlecell/Rfunction/getPCs.R")
source("/st07/sccDB/Spatial/Anchor_R/anchor_method2.r")
source("/home/singlecell/Rfunction/computeCommunProb_jiwoo.r")
source("/st07/kobic/data/code/cellchat_input.r")


meta_df <- read.csv("/st07/sccDB/Spatial/HISSTA/pipeline/db_kobic.csv")

# -----------------------------------------------------------
# meta filtering
# -----------------------------------------------------------
meta_imsi <- meta_df %>% filter(tissue=="Ovary")
sample_list <- meta_imsi$sample_id


# -----------------------------------------------------------
# celltype list
# -----------------------------------------------------------
celltype_list <- c("Malignant","Ciliated","FTSEC",
"TNKILC","B","Plasma", "Macrophage","Monocyte","DC","pDC","Mast","Neutrophil",
"Endothelial","Fibroblast","unknown") 
names(celltype_list) <- c("#bf4040",hue_pal()(length(celltype_list)- 2),"gray88")



############################################
# 1. load data & QC & normalize
############################################

for( ID in sample_list){
    
    tmp_meta <- meta_imsi%>%filter(sample_id==ID)
    
    # set path
    if(tmp_meta$platform == "Xenium"){
        obj_path = paste0("/st07/sccDB/Spatial/01.Xenium/14.human_ovary/01.public_data/",ID)
    } else if(tmp_meta$platform == "MERFISH"){
        obj_path = paste0("/st07/sccDB/Spatial/03.MERFISH/06.human_ovary/01.public_data/",ID)
    } else if(tmp_meta$platform == "CosMX"){
        obj_path = paste0("/st07/sccDB/Spatial/02.CosMx/09.human_ovary/",ID)
    } 

    # set kobic_path
    path =paste0("/st07/kobic/data/dataset/",ID)

    if (!dir.exists(paste0(path))) {
        dir.create(paste0(path))
    }

    if (!dir.exists(paste0(path, "/rds/"))) {
        dir.create(paste0(path,"/rds/"))
    }


    #-------------------------------------------------------------
    # 1.processing
    #-------------------------------------------------------------
    print(paste0("1. Processing ",ID))

    if (!dir.exists(paste0(obj_path, "/processed/"))) {
        dir.create(paste0(obj_path,"/processed/"))
    }

    # Load data
    if(tmp_meta$platform=="Xenium"){

        data <- ReadXenium(
            data.dir = paste0(obj_path,"/rawdata/"),
            type = c("centroids", "segmentations"),
            )

        segmentations.data <- list(
            "centroids" = CreateCentroids(data$centroids),
            "segmentation" = CreateSegmentation(data$segmentations))

        coords <- CreateFOV(
            coords = segmentations.data,
            type = c("segmentation", "centroids"),
            molecules = data$microns,
            assay = "Xenium"
        )

        raw_obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = "Xenium")

        if("Blank Codeword" %in% names(data$matrix)){
            raw_obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Blank Codeword"]])
            raw_obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
            raw_obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
        }else{
            raw_obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
            raw_obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
        }

        raw_obj[["fov"]] <- coords  


        cells <- read.csv(paste0(obj_path,"/rawdata/cells.csv.gz"))

        raw_obj$cell_id <- colnames(raw_obj)
        raw_obj$nucleus_area <- cells$nucleus_area
        raw_obj$cell_area <- cells$cell_area

        coords <- raw_obj@images$fov@boundaries$centroids@coords %>%as.data.frame
        raw_obj$x_global_um <- coords$x
        raw_obj$y_global_um <- coords$y

        tmp_assay="Xenium"
        tmp_ctrl_assay="ControlProbe"
        tmp_area_column ="nucleus_area"
        tmp_x ="x_global_um" 
        tmp_y ="y_global_um"

    } else if(tmp_meta$platform=="MERFISH") {
        raw_obj <- LoadVizgen(paste0(obj_path,"/rawdata"), fov = "fov")
        cells <- read.csv(paste0(obj_path,"/rawdata/cell_metadata.csv"))

        raw_obj$cell_id <- colnames(raw_obj)
        raw_obj$cell_area <- cells$volume

        # Get blank matrix
        mat <- read.csv(paste0(obj_path,"/rawdata/cell_by_gene.csv"))
        blank <- colnames(mat)[grep("Blank", colnames(mat))]
        blank_mat <- t(mat[, blank])
        colnames(blank_mat)<-colnames(raw_obj)

        # Add BlankAssay
        raw_obj[["BlankAssay"]] <- CreateAssayObject(counts = blank_mat)

        coords <- raw_obj@images$fov@boundaries$centroids@coords %>%as.data.frame
        raw_obj$x_global_um <- coords$x
        raw_obj$y_global_um <- coords$y

        tmp_assay="Vizgen"
        tmp_ctrl_assay="BlankAssay"
        tmp_area_column ="cell_area"
        tmp_x ="x_global_um" 
        tmp_y ="y_global_um" 

    } else if(tmp_meta$platform=="CosMX") {
        ## count matrix --------------------------------------------------
        count.file <- list.files(paste0(obj_path, "/rawdata"))
        count.file <- count.file[grepl('[_a-zA-Z0-9]*_exprMat_file.csv', count.file)]
        count.file <- paste0(obj_path, "/rawdata/", count.file)

        count_matrix <- read.csv(count.file)
        rownames(count_matrix) <- paste0(as.character(count_matrix$cell_ID), "_", count_matrix$fov)
        count_matrix <- count_matrix[count_matrix$cell_ID != 0, -c(1,2)] #cell_ID가 0인건 cell에 assign안된 transcript들을 모아놓은것
        count_matrix <- t(as.matrix(count_matrix))


        genes <- rownames(count_matrix)
        genes <- gsub("\\.", "-", genes)
        rownames(count_matrix) <- genes
        
        count_matrix <- as(count_matrix, "dgCMatrix")

        negprb <- count_matrix[grepl("^NegPrb", rownames(count_matrix)), ]
        count_matrix <-  count_matrix[!grepl("^NegPrb", rownames(count_matrix)), ]

        ## metadata --------------------------------------------------
        metadata.file <- list.files(paste0(obj_path, "/rawdata"))
        metadata.file <- metadata.file[grepl('[_a-zA-Z0-9]*_metadata_file.csv', metadata.file)]
        metadata.file <- paste0(obj_path, "/rawdata/", metadata.file)

        metadata <- read.csv(metadata.file)
        metadata$cell_id <- paste0(as.character(metadata$cell_ID), "_", metadata$fov)
        rownames(metadata) <- metadata$cell_id

        print(head(metadata))

        ## fov --------------------------------------------------
        cents <- metadata[, c("CenterX_global_px", "CenterY_global_px", "cell_id")]
        colnames(cents) <- c("x", "y", "cell")
        cents <- CreateCentroids(cents)

        ## polygon --------------------------------------------------
        polygon.file <- list.files(paste0(obj_path, "/rawdata"))
        polygon.file <- polygon.file[grepl('[_a-zA-Z0-9]*-polygons.csv', polygon.file)]
        polygon.file <- paste0(obj_path, "/rawdata/", polygon.file)

        polygons <- read.csv(polygon.file)
        polygons$cell <- paste0(as.character(polygons$cellID), "_", polygons$fov)

        polygons <- polygons %>% dplyr::select(x_global_px, y_global_px, cell)
        colnames(polygons) <- c("x", "y", "cell")
        segs <- CreateSegmentation(polygons)

        ## molecules --------------------------------------------------
        molecules.file <- list.files(paste0(obj_path, "/rawdata"))
        molecules.file <- molecules.file[grepl('[_a-zA-Z0-9]*_tx_file.csv', molecules.file)]
        molecules.file <- paste0(obj_path, "/rawdata/", molecules.file)

        molecules_xy <- read.csv(molecules.file)
        colnames(molecules_xy)[colnames(molecules_xy) == "target"] <- "gene"
        colnames(molecules_xy)[colnames(molecules_xy) == "x_global_px"] <- "x"
        colnames(molecules_xy)[colnames(molecules_xy) == "y_global_px"] <- "y"
        
        molecules_xy <- molecules_xy[, c("x", "y", "gene")]
        molecules_xy <- CreateMolecules(molecules_xy)

        # creating fov --------------------------------------------------
        segmentations.data <- list(
            "centroids" = cents,
            "segmentation" = segs
        )

        coords <- CreateFOV(coords = segmentations.data, 
                            type = c("segmentation", "centroids"), 
                            molecules = molecules_xy, 
                            assay = "Nanostring") #cosmx
        
        raw_obj <- CreateSeuratObject(counts = count_matrix, meta.data = metadata, assay = "Nanostring")
        negprb <- CreateAssayObject(counts = negprb)

        raw_obj[["negprobes"]] <- negprb
        raw_obj[["fov"]] <- subset(coords, cell=Cells(raw_obj))

        colnames(raw_obj@meta.data)[grepl("Area",colnames(raw_obj@meta.data))] <- "cell_area"
        colnames(raw_obj@meta.data)[grepl("CenterX_global_px",colnames(raw_obj@meta.data))] <- "x_global_px"
        colnames(raw_obj@meta.data)[grepl("CenterY_global_px",colnames(raw_obj@meta.data))] <- "y_global_px"

        colnames(raw_obj@meta.data)[grepl("CenterX_local_px",colnames(raw_obj@meta.data))] <- "x_fov_px"
        colnames(raw_obj@meta.data)[grepl("CenterY_local_px",colnames(raw_obj@meta.data))] <- "y_fov_px"
        colnames(raw_obj@meta.data)[grepl("fov",colnames(raw_obj@meta.data))] <- "FOV"

        tmp_assay="Nanostring"
        tmp_ctrl_assay="negprobes"
        tmp_area_column ="cell_area"
        tmp_x ="x_global_px" 
        tmp_y ="y_global_px" 
    }
 
    print(dim(raw_obj))


    # Add meta data
    raw_obj$sample_id <- tmp_meta$sample_id
    raw_obj$patient_id <- tmp_meta$patient_id
    raw_obj$project_id <- tmp_meta$project_id
    raw_obj$platform <- tmp_meta$platform
    raw_obj$data_resource <- tmp_meta$data_resource
    raw_obj$tissue <- tmp_meta$tissue
    raw_obj$condition <- tmp_meta$condition
    raw_obj$disease <- tmp_meta$disease
    raw_obj$panel <- tmp_meta$panel
    raw_obj$panel_plex <- tmp_meta$panel_plex
    raw_obj$preservation_method <- tmp_meta$preservation_method
    raw_obj$staining_method <- tmp_meta$staining_method
    raw_obj$sample_name <- tmp_meta$sample_name

    raw_obj@misc$celltype_list <- celltype_list



    # QC
    if(tmp_meta$platform=="Xenium"){
        obj <- subset(raw_obj,
                    subset = nCount_Xenium >= 20 & nFeature_Xenium >= 10)
        raw_obj$QC <-"pass"
        raw_obj$QC[raw_obj$nCount_Xenium<20|raw_obj$nFeature_Xenium < 10] <-"low_quality_cells"

        cols <-c("gray","red")
        FeatureScatter(object = raw_obj, feature1 = "nCount_Xenium", feature2 = "nFeature_Xenium", 
            raster=FALSE,group.by="QC",col=cols )+theme_bw()
        ggsave(paste0(obj_path,"/processed/QC_Check.png"),width=8)


    } else if (tmp_meta$platform=="MERFISH") {
        obj <- subset(raw_obj,
                subset = nCount_Vizgen >= 20 & nFeature_Vizgen >= 10)
        raw_obj$QC <-"pass"
        raw_obj$QC[raw_obj$nCount_Vizgen<20|raw_obj$nFeature_Vizgen < 10] <-"low_quality_cells"

        cols <-c("gray","red")
        FeatureScatter(object = raw_obj, feature1 = "nCount_Vizgen", feature2 = "nFeature_Vizgen",
            raster=FALSE,group.by="QC",col=cols )+theme_bw()
        ggsave(paste0(obj_path,"/processed/QC_Check.png"),width=8)

    } else if (tmp_meta$platform=="CosMX") {
        obj <- subset(raw_obj,
                subset = nCount_Nanostring >= 20 & nFeature_Nanostring >= 10)
        raw_obj$QC <-"pass"
        raw_obj$QC[raw_obj$nCount_Nanostring<20|raw_obj$nFeature_Nanostring < 10] <-"low_quality_cells"

        cols <-c("gray","red")
        FeatureScatter(object = raw_obj, feature1 = "nCount_Nanostring", feature2 = "nFeature_Nanostring",
            raster=FALSE,group.by="QC",col=cols )+theme_bw()
        ggsave(paste0(obj_path,"/processed/QC_Check.png"),width=8)

    }
    print(dim(obj))

    removed <- 100 - (ncol(obj)/ncol(raw_obj)*100)
    removed <- round(removed,2)

    print(paste0("before QC : ",ncol(raw_obj),"      after QC : ",ncol(obj),"      ",removed," % filtered "))

    print("QC result")
    print(table(raw_obj$QC))

    saveRDS(raw_obj,paste0(path, "/rds/obj_raw.rds"))


    # SCTransform
    obj <- SCTransform(obj, assay = tmp_assay)
    obj <- RunPCA(obj,npcs = 50, features = rownames(obj))

    # Get PCs and make elbow plot
    npcs <- min(getPCs(obj)) 
    npcs 
    ElbowPlot(obj, ndims = 50) + geom_vline(xintercept = npcs, lty = "dashed",color = "red")+theme_bw()
    ggsave(paste0(obj_path,"/processed/elbowplot_pc.png"),width=14)

    obj <- RunUMAP(obj, dims = 1:npcs)
    obj <- FindNeighbors(obj, reduction = "pca", dims = 1:npcs)
    obj <- FindClusters(obj,resolution = c(0.3,0.5,0.7,1))

    Res=c("SCT_snn_res.0.3","SCT_snn_res.0.5","SCT_snn_res.0.7","SCT_snn_res.1")

    for(res in Res){
        obj@meta.data[,res] <- paste("Cluster", obj@meta.data[,res])
        obj@meta.data[,res]<- droplevels(factor(obj@meta.data[,res], levels = paste("Cluster", c(0:100))))
        print(table(obj@meta.data[,res]))

        DimPlot(obj,raster=FALSE,label = TRUE, group.by=res,label.size=5)
        ggsave(paste0(obj_path,"/processed/Dimplot_pc",npcs,"_",res,".png"), width = 8, height = 7)

        splitDimPlot(obj, res, dim_method = "pc", npcs, sub_title = NULL, width = 20, dim_reduc = "umap", path = paste0(obj_path,"/processed/"))
    }

    saveRDS(obj,paste0(obj_path,"/processed/obj_processed.rds"))
}


############################################
# 2. cell type annotation
############################################

# -----------------------------------------------------------
# marker list
# -----------------------------------------------------------
l2_markers <- read.csv("/st07/sccDB/Spatial/04.scRNA_reference/09.ovary/findmarkers/ovary_markers_celltype_l2.csv",row.names=1)
add_genes <- read.csv("/st07/sccDB/Spatial/04.scRNA_reference/09.ovary/add_genes.csv")

# -----------------------------------------------------------
# malignant marker definition
# -----------------------------------------------------------
malignant_list<-c("SOX17","ESR1","EPCAM","ERBB3","KRAS","BBOX1","S100A1","ERBB2","C1R")

# -----------------------------------------------------------
# celltype list
# -----------------------------------------------------------
celltype_list <- c("Malignant","Ciliated","FTSEC",
"TNKILC","B","Plasma", "Macrophage","Monocyte","DC","pDC","Mast","Neutrophil",
"Endothelial","Fibroblast","unknown") 
names(celltype_list) <- c("#bf4040",hue_pal()(length(celltype_list)- 2),"gray88")

# -----------------------------------------------------------
# sc reference profile
# -----------------------------------------------------------
profiles <- read.csv("/st07/sccDB/Spatial/04.scRNA_reference/09.ovary/ovary_insitutype_profiles.csv",row.names=1)
profiles <- as.matrix(profiles)

for( ID in sample_list[4:6]){

    tmp_meta <- meta_imsi%>%filter(sample_id==ID)

    # set path
    if(tmp_meta$platform == "Xenium"){
        obj_path = paste0("/st07/sccDB/Spatial/01.Xenium/14.human_ovary/01.public_data/",ID)
        tmp_assay="Xenium"
        tmp_ctrl_assay="ControlProbe"
        tmp_area_column ="nucleus_area"
        tmp_x ="x_global_um" 
        tmp_y ="y_global_um"
    } else if(tmp_meta$platform == "MERFISH"){
        obj_path = paste0("/st07/sccDB/Spatial/03.MERFISH/06.human_ovary/01.public_data/",ID)
        tmp_assay="Vizgen"
        tmp_ctrl_assay="BlankAssay"
        tmp_area_column ="cell_area"
        tmp_x ="x_global_um" 
        tmp_y ="y_global_um" 
    } else if(tmp_meta$platform == "CosMX"){
        obj_path = paste0("/st07/sccDB/Spatial/02.CosMx/09.human_ovary/",ID)
        tmp_assay="Nanostring"
        tmp_ctrl_assay="negprobes"
        tmp_area_column ="cell_area"
        tmp_x ="x_global_px" 
        tmp_y ="y_global_px" 
    } 

    # set kobic_path
    path =paste0("/st07/kobic/data/dataset/",ID)

    if (!dir.exists(paste0(path))) {
        dir.create(paste0(path))
    }

    if (!dir.exists(paste0(path, "/rds/"))) {
        dir.create(paste0(path,"/rds/"))
    }
    
    
    obj <- readRDS(paste0(obj_path,"/processed/obj_processed.rds"))
    obj@misc$celltype_list <- celltype_list

    #-------------------------------------------------------------
    # cell type annotation
    #-------------------------------------------------------------
    print(paste0("2. Assign cell types ",ID))

    ## anchor-based
    if (!dir.exists(paste0(obj_path, "/processed/anchor_based/"))) {
        dir.create(paste0(obj_path,"/processed/anchor_based"))
    }

    scref_marker_list <-SCMarkerSelection(obj,
                                        marker_list=l2_markers,
                                        n_top=20,
                                        logfc.threshold =1,
                                        rm_overlap = TRUE,
                                        add_genes = add_genes)
    head(scref_marker_list)
    print(table(scref_marker_list$cluster))


    core_marker_list<-CoreMarkerSelection(obj,
                                        marker_list=l2_markers,
                                        n_top=5,
                                        logfc.threshold =1,
                                        rm_overlap = TRUE,
                                        add_genes = add_genes)
    head(core_marker_list)
    print(table(core_marker_list$cluster))


    if( tmp_meta$condition == "Cancer"){
        obj <- FindMalignantCells_bimodal(obj,
                                assay=tmp_assay,
                                ctrl_assay = tmp_ctrl_assay,
                                marker_list=malignant_list,
                                add_area=TRUE,
                                area_column =tmp_area_column,
                                smooth=TRUE,
                                cutoff=0.95,
                                path=paste0(obj_path,"/processed/anchor_based/"))

        obj <- FindNonMalignantCells(obj,
                                assay=tmp_assay,
                                ctrl_assay = tmp_ctrl_assay,
                                malignant_column = "sup.malignant_km.clust",
                                core_marker=core_marker_list,
                                merged_batch_name = NULL,
                                smooth=TRUE,
                                scref_marker=scref_marker_list,
                                anchor_percent=20,
                                cutoff=0.95,
                                path=paste0(obj_path,"/processed/anchor_based/"))

    } else if (tmp_meta$condition %in% c("Healthy","Other disease")){
        obj <- FindNonMalignantCells(obj,
                                assay=tmp_assay,
                                ctrl_assay = tmp_ctrl_assay,
                                malignant_column = NULL,
                                core_marker=core_marker_list,
                                merged_batch_name = NULL,
                                smooth=TRUE,
                                scref_marker=scref_marker_list,
                                anchor_percent=20,
                                cutoff=0.95,
                                path=paste0(obj_path,"/processed/anchor_based/"))
    }

    obj$celltype_anchor <- obj$anchor_final_clust
    print("Anchor result : ")
    print(table(obj$celltype_anchor ))

    obj$celltype_anchor <- droplevels(factor(obj$celltype_anchor, levels = celltype_list))
    colz <- names(celltype_list[celltype_list %in% names(table(obj$celltype_anchor))])
    DimPlot(obj,raster=FALSE, label = TRUE, group.by= "celltype_anchor",
        label.size=5, pt.size = 0.001, cols = colz, repel = TRUE)
    ggsave(paste0(obj_path,"/processed/dimplot_celltype_anchor.png"), width = 12, height = 8)


    # 불필요한 columns 지우기
    patterns <- c("malignant", "signature", "clust","cut")
    cols_to_remove <- sapply(patterns, function(p) grepl(p, names(obj@meta.data)))
    remove_indices <- which(rowSums(cols_to_remove) > 0)
    obj@meta.data <- obj@meta.data[ , -remove_indices]


    ## marker-based
    # marker list 만들기
    markers <- data.frame(matrix("", nrow = 50, ncol = length(unique(scref_marker_list$cluster))))

    if(tmp_meta$condition == "Cancer"){ # cancer의 경우
        marker_list <- rbind(scref_marker_list, 
                            data.frame("cluster" = "Malignant", "gene" = malignant_list))
        markers <- data.frame(matrix("", nrow = max(table(marker_list$cluster)), ncol = length(table(marker_list$cluster))))

        types <- unique(marker_list$cluster)
        for(i in 1:length(types)){
                tmp <- marker_list %>% filter(cluster==types[i])
                markers[1:nrow(tmp),i] <- tmp$gene
                names(markers)[i] <- types[i]
        }
        # write.csv(markers,paste0(path, path2, "scina_marker_list.csv"), row.names=FALSE)

    } else if (tmp_meta$condition %in% c("Healthy","Other disease")){
        marker_list <- scref_marker_list
        markers <- data.frame(matrix("", nrow = max(table(marker_list$cluster)), ncol = length(table(marker_list$cluster))))

        types <- unique(marker_list$cluster)
        for(i in 1:length(types)){
                tmp <- marker_list %>% filter(cluster==types[i])
                markers[1:nrow(tmp),i] <- tmp$gene
                names(markers)[i] <- types[i]
        }
        # write.csv(markers,paste0(path, path2, "scina_marker_list.csv"), row.names=FALSE)
    }

    mat <-GetAssayData(obj, assay = "SCT", layer = "data") %>% as.matrix

    results = SCINA(mat, 
                markers,
                max_iter = 100,
                convergence_n = 10,
                convergence_rate = 0.999,
                sensitivity_cutoff = 0.9,
                allow_unknown = 0,
                rm_overlap=FALSE,
                log_file = paste0(obj_path,"/processed/SCINA.log")) #allow_unknown=0 : False, 1 : Ture

    print("SCINA result : ")
    print(table(results$cell_labels))
    obj$celltype_marker <- results$cell_labels

    obj$celltype_marker <- droplevels(factor(obj$celltype_marker, levels = celltype_list))
    colz <- names(celltype_list[celltype_list %in% names(table(obj$celltype_marker))])
    DimPlot(obj,raster=FALSE, label = TRUE, group.by= "celltype_marker",
        label.size=5, pt.size = 0.001, cols = colz, repel = TRUE)
    ggsave(paste0(obj_path,"/processed/dimplot_celltype_marker.png"), width = 12, height = 8)



    ## profile-based
    counts <-GetAssayData(obj, assay = tmp_assay, layer = "counts")%>%as.matrix
    neg <- GetAssayData(obj, assay = tmp_ctrl_assay, layer = "counts")%>%as.matrix
    obj$negmean <- Matrix::colMeans(neg)

    ist.res = insitutypeML(x = t(counts),
						neg = obj@meta.data[, "negmean"],
						reference_profiles = as.matrix(profiles),
						cohort = NULL)

    print("InSituType result : ")
    print(table(ist.res$clust))


    ist.res2 <- cbind(ist.res$clust, ist.res$prob, ist.res$clust)
    colnames(ist.res2) <- c("ist.res.clust", "ist.res.prob", "ist.res.cut")
    ist.res2 <- as.data.frame(ist.res2)
    ist.res2$ist.res.cut[ist.res2$ist.res.prob < 0.95] <- "noise"
    print(table(ist.res2$ist.res.cut))

    write.csv(ist.res2, paste0(obj_path,"/processed/IST_res.csv"))

    obj$celltype_profile <- ist.res$clust

    obj$celltype_profile <- droplevels(factor(obj$celltype_profile, levels = celltype_list))
    colz <- names(celltype_list[celltype_list %in% names(table(obj$celltype_profile))])
    DimPlot(obj,raster=FALSE, label = TRUE, group.by= "celltype_profile",
        label.size=5, pt.size = 0.001, cols = colz, repel = TRUE)
    ggsave(paste0(obj_path,"/processed/dimplot_celltype_profile.png"), width = 12, height = 8)


    saveRDS(obj,paste0(obj_path,"/processed/obj_celltype.rds"))
}


############################################
# 3. downstream analysis
############################################
library(Seurat)
library(GSVA)
library(qusage)
library(ggplot2)
library(Biobase)


for( ID in sample_list){

    tmp_meta <- meta_imsi%>%filter(sample_id==ID)

    # set path
    if(tmp_meta$platform == "Xenium"){
        obj_path = paste0("/st07/sccDB/Spatial/01.Xenium/14.human_ovary/01.public_data/",ID)
        tmp_assay="Xenium"
        tmp_ctrl_assay="ControlProbe"
        tmp_area_column ="nucleus_area"
        tmp_x ="x_global_um" 
        tmp_y ="y_global_um"
    } else if(tmp_meta$platform == "MERFISH"){
        obj_path = paste0("/st07/sccDB/Spatial/03.MERFISH/06.human_ovary/01.public_data/",ID)
        tmp_assay="Vizgen"
        tmp_ctrl_assay="BlankAssay"
        tmp_area_column ="cell_area"
        tmp_x ="x_global_um" 
        tmp_y ="y_global_um" 
    } else if(tmp_meta$platform == "CosMX"){
        obj_path = paste0("/st07/sccDB/Spatial/02.CosMx/09.human_ovary/",ID)
        tmp_assay="Nanostring"
        tmp_ctrl_assay="negprobes"
        tmp_area_column ="cell_area"
        tmp_x ="x_global_px" 
        tmp_y ="y_global_px" 
    } 
    # set kobic_path
    path =paste0("/st07/kobic/data/dataset/",ID)

    if (!dir.exists(paste0(path))) {
        dir.create(paste0(path))
    }

    if (!dir.exists(paste0(path, "/rds/"))) {
        dir.create(paste0(path,"/rds/"))
    }
    
    
    obj <- readRDS(paste0(obj_path,"/processed/obj_celltype.rds"))


    ##findmarkers
    tmp_list<- list()
    cols =c("SCT_snn_res.0.3","SCT_snn_res.0.5","SCT_snn_res.0.7",
    "SCT_snn_res.1","celltype_anchor","celltype_marker","celltype_profile")
    for( col in cols){
        Idents(obj) <- obj@meta.data[,col]
        findmarkers <- FindAllMarkers(obj, only.pos = TRUE)
        tmp_list[[col]] <- findmarkers
    }
    obj@misc$DEG <- tmp_list


    #-------------------------------------------------------------
    # colocalization
    #-------------------------------------------------------------
    print(paste0("3. Calculate colocalization ",ID))

    sce <- as.SingleCellExperiment(obj, assay="SCT")

    tmp_list<- list()
    methods=c("celltype_anchor","celltype_marker","celltype_profile")
    
    for( method in methods){
        spe <- SpatialExperiment(
            assays = list(counts = assay(sce)),
            rowData = rowData(sce), 
            colData = colData(sce),
            spatialCoordsNames = c(tmp_x, tmp_y))

        spe@colData[,method]<- as.character(spe@colData[,method])
        spe <- readHoodData(spe, anno_col = method)

        fnc <- findNearCells(spe, k = 500)

        pm <- scanHoods(fnc$distance)
        hoods <- mergeByGroup(pm, fnc$cells)
        spe <- mergeHoodSpe(spe, hoods)
        spe <- calcMetrics(spe, pm_cols = colnames(hoods))

        dat <- as.data.frame(colData(spe), optional = TRUE)
        pm_cols = colnames(hoods)
        corr <- round(cor(dat[, pm_cols]), 1)

        tmp_list[[method]] <- corr
    }
    obj@misc$spe <- tmp_list


    #-------------------------------------------------------------
    # cellchat
    #-------------------------------------------------------------
    
    spatial.locs = GetTissueCoordinates(obj, scale = NULL)
    rownames(spatial.locs) <- spatial.locs$cell
    spatial.locs <- spatial.locs[,-3]


    if(tmp_meta$platform=="Xenium"){
        conversion.factor = 1 ##um인 경우
    }else if(tmp_meta$platform=="MERFISH"){
        conversion.factor = 1 ##um인 경우
    }else if(tmp_meta$platform=="CosMX"){
        conversion.factor = 0.18 ##mm인 경우
    }

    spot.size = 10 # use the typical human cell size
    spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)
    
    ## use all CellChatDB for cell-cell communication analysis ++++++++
    CellChatDB.use <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

    ptm = Sys.time()

    tmp_list<- list()
    for(method in methods){
        cellchat <- createCellChat(object = obj, group.by = method, 
            assay = "SCT", datatype = "spatial", coordinates = spatial.locs, 
            spatial.factors = spatial.factors)

        cellchat@DB <- CellChatDB.use

        cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
        future::plan("multisession", workers = 4) # do parallel
        cellchat <- identifyOverExpressedGenes(cellchat, do.DE = FALSE, min.cells = 10) #panel 수가 적은 경우
        cellchat <- identifyOverExpressedInteractions(cellchat)

        cellchat <- computeCommunProb_jiwoo(cellchat, type = "truncatedMean", contact.range = 10)
        cellchat <- filterCommunication(cellchat, min.cells = 10)

        ## Infer the cell-cell communication at a signaling pathway level
        cellchat <- computeCommunProbPathway(cellchat)

        ## Calculate the aggregated cell-cell communication network
        cellchat <- aggregateNet(cellchat)
        cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

        tmp_list[[method]] <- cellchat
    }
    obj@misc$cellchat <- tmp_list



    #-------------------------------------------------------------
    # niche
    #-------------------------------------------------------------
    
    coord <- GetTissueCoordinates(obj)
    ratio <- c(range(coord$x)[2] - range(coord$x)[1], range(coord$y)[2] - range(coord$y)[1])
    max_value <- 14; 
    scale_factor <- max_value / max(ratio); scaled_ratio <- ratio * scale_factor


    methods=c("celltype_anchor","celltype_marker","celltype_profile")

    
    niche_list <- list()
    ssgsea_list <- list()
    heatmap_list <- list()
    for( method in methods){
        
        obj <- BuildNicheAssay(object = obj, fov = "fov", group.by = method,niches.k = 5, neighbors.k = 30)

        obj$niches  <- paste("Niche", obj$niches )
        obj$niches<- factor(obj$niches ,levels = paste("Niche", c(1:5)))

        mycols <- names(obj@misc$celltype_list[obj@misc$celltype_list %in% names(table(obj@meta.data[,method]))])

        celltype.plot <- ImageDimPlot(obj, group.by = method, size = 0.5, cols = mycols,
            dark.background = F) + ggtitle(method)+coord_flip()+
            theme(legend.text=element_text(size=14) ,legend.title=element_text(size=20),
            plot.title = element_text(size=25,hjust = 0.5))

        niche.plot <- ImageDimPlot(obj, group.by = "niches", size = 0.5, dark.background = F) + 
            ggtitle(paste0(method,"_niches")) +
            scale_fill_manual(values = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B"))+coord_flip()+
            theme(legend.text=element_text(size=14) ,legend.title=element_text(size=20),
            plot.title = element_text(size=25,hjust = 0.5))

        celltype.plot | niche.plot
        ggsave(paste0(obj_path,"/processed/",method,"_niche_com_imgdim.png"),
        width = (scaled_ratio[1]+1)*2, height = scaled_ratio[2])

        niche.plot+theme(plot.title = element_text(size=0))
        ggsave(paste0(obj_path,"/processed/",method,"_niche_imgdim.png"),
                    width = (scaled_ratio[1]+2), height = scaled_ratio[2])



        df <- obj@meta.data %>%
        dplyr::select(niches, !!sym(method)) %>%
        group_by(niches, !!sym(method)) %>%
        summarise(count = n(), .groups = 'drop') %>%
        group_by(niches) %>%
        mutate(proportion = count / sum(count)) %>%
        ungroup()

        colnames(df)[1] <-"niches"
        colnames(df)[2] <-"celltype"

        niche_list[[method]] <- df
        obj@meta.data[,paste0(method,"_niches")] <- obj$niches
    }

    ##ssGSEA
    ssgsea_list <- list()
    heatmap_list <- list()

    Niches=c("celltype_anchor","celltype_marker","celltype_profile")

    for( niche in Niches){

        niche = paste0(niche,"_niches")
        mat <-GetAssayData(obj, assay = "SCT", layer = "data") %>% as.matrix
        colnames(mat) <- obj@meta.data[,niche]
    
        niches <- names(table(obj@meta.data[,niche]))

        df <- do.call(cbind, lapply(niches, function(type) {
        celltype_index <- which(colnames(mat) == type)
        celltype_mat <- mat[, celltype_index, drop = FALSE]
        avg_mat <- rowMeans(celltype_mat, na.rm = TRUE)
        setNames(data.frame(avg_mat), type)
        }))
        df <- as.matrix(df)

        sc.pheno <- data.frame(check.names=F, check.rows=F, stringsAsFactors=F,
        row.names=colnames(df),cellType=colnames(df))

        sc.meta <- data.frame(labelDescription="cellType",row.names="cellType")
        sc.pdata <- new("AnnotatedDataFrame",data=sc.pheno,varMetadata=sc.meta)
        sc.eset <- Biobase::ExpressionSet(assayData=data.matrix(df),
            phenoData=sc.pdata)

        hallmark <- read.gmt("/st07/sccDB/Spatial/ssGSEA/h.all.v7.5.1.symbols.gmt")

        test_es <- gsva(expr=sc.eset, gset.idx.list=hallmark, mx.diff=F,
            method="ssgsea", verbose=F, ssgsea.norm=T, min.sz=1, tau=0.25)

        res <- test_es@assayData$exprs %>%as.data.frame

        rownames(res) <- sub("^HALLMARK_", "", rownames(res))

        res_df <- lapply(niches, function(celltype) {
            tmp <- res[,celltype] %>% as.data.frame
            colnames(tmp) <- "Enrichment_score"
            tmp <- mutate(tmp, Pathway = rownames(res)) %>%
                    mutate(tmp, Niches = celltype) %>%
                dplyr::select(Pathway, Niches,Enrichment_score)
            return(tmp)
        }) %>%
                bind_rows()

 
        ssgsea_list[[niche]] <- res
        heatmap_list[[niche]] <- res_df
    }
    obj@misc$niche <- niche_list
    obj@misc$ssgsea <- ssgsea_list
    obj@misc$heatmap <- heatmap_list  

    obj$niches <- NULL

    saveRDS(obj,paste0(path,"/rds/obj_final.rds"))
}


#################################################################################
# 4. Make plotting data
#################################################################################

for( ID in sample_list){

    print(paste0("Start  ", ID))
    

    # set path
    path=paste0("/st07/kobic/data/dataset/",ID)
    
    obj<- readRDS(paste0(path,"/rds/obj_raw.rds"))


    if (!dir.exists(paste0(path, "/meta/"))) {
        dir.create(paste0(path,"/meta/"))
    }

    # -----------------------------
    ## raw_cells_info.csv 
    #------------------------------
    feat <- colnames(obj@meta.data)
    feat <- feat[grepl("(nCount|nFeature)", feat) & (grepl("(Xenium|Nanostring|Vizgen)", feat))]

    raw_cells_info <- obj@meta.data[, c("cell_id", feat,"QC")]


    obj <- readRDS(paste0(path, "/rds/obj_final.rds"))


    # -----------------------------
    ## processed_cells_info.csv 
    #------------------------------
    feat <- colnames(obj@meta.data)
    feat <- feat[grepl("SCT_snn", feat) | (grepl("(celltype_anchor|marker|profile)", feat))] #clustering, celltype annotation column명 추출

    cells_info <- obj@meta.data[, c("cell_id",feat)]

    cells_info_merged <- left_join(raw_cells_info,cells_info)
    cells_info_merged%>% write.csv(paste0(path,"/meta/",ID,"_cells_info.csv"), row.names=FALSE)



    # -----------------------------
    ## deg.csv 
    #------------------------------
    feat <- names(obj@misc$DEG)
    df <- lapply(feat, function(featcol) {
        findmarkers <- obj@misc$DEG[[featcol]]
        findmarkers <- mutate(findmarkers, type = featcol) %>%
            dplyr::select(type, everything())
        return(findmarkers)
        }) %>%
            bind_rows()
    rownames(df) <- NULL
    write.csv(df,paste0(path,"/meta/",ID,"_deg.csv"), row.names=FALSE)


    # -----------------------------
    ## dot.csv 파일
    #------------------------------
    feat <- names(obj@misc$DEG)
    dot_list<-list()
    for( featcol in feat){
        Idents(obj) <- obj@meta.data[, featcol]
        findmarkers <- obj@misc$DEG[[featcol]]

        unique_genes <- findmarkers %>%
            group_by(gene) %>%
            slice_max(order_by = avg_log2FC, n = 1) %>%
            ungroup()

        top_markers <- unique_genes %>% 
            group_by(cluster) %>% # 각 cluster별로 
            slice_max(order_by = avg_log2FC, n = 3) %>% # avg_log2FC가 제일 큰 3개의 행을 추출
            arrange(cluster, desc(avg_log2FC))

        top_markers$gene

        # idents를 기준으로 cell 꺼내기
        cells <- unlist(x = CellsByIdentities(object = obj, cells = colnames(obj[["SCT"]]), idents = obj@meta.data[, featcol]))

        #feature*cell 발현값꺼내기 (df)
        DefaultAssay(obj) <- "SCT"

        data.features <- FetchData(object = obj, vars = top_markers$gene, cells = cells)
        data.features$id <- obj[[featcol, drop = TRUE]][cells, drop = TRUE] # ident를 id에 넣어줌
        id.levels <- levels(x = data.features$id)
        data.features$id <- as.vector(x = data.features$id) # id를 vector로 변환하고 level 설정

        #특징 평균 및 백분율 계산
        data.plot <- lapply(
        X = unique(x = data.features$id),
        FUN = function(ident) {
            data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
            avg.exp <- apply(
            X = data.use,
            MARGIN = 2,
            FUN = function(x) {
                return(mean(x = expm1(x = x)))
            }
            )
            pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
            return(list(avg.exp = avg.exp, pct.exp = pct.exp))
            }
        )
        names(x = data.plot) <- unique(x = data.features$id)

        #데이터 클러스터링 및 변환
        data.plot <- lapply(
        X = names(x = data.plot),
        FUN = function(x) {
            data.use <- as.data.frame(x = data.plot[[x]])
            data.use$features.plot <- rownames(x = data.use)
            data.use$id <- x
            return(data.use)
            }
        )
        data.plot <- do.call(what = 'rbind', args = data.plot)


        #평균 표현 및 스케일링 계산
        avg.exp.scaled <- sapply(
        X = unique(x = data.plot$features.plot),
        FUN = function(x) {
            data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
            data.use <- scale(x = log1p(data.use))
            data.use <- MinMax(data = data.use, min = -2.5, max = 2.5)
            return(data.use)
            }
        )
        avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))

        #색상 분할 처리
        data.plot$avg.exp.scaled <- avg.exp.scaled
        data.plot$features.plot <- factor(
        x = data.plot$features.plot,
        levels = top_markers$gene
        )

        #백분율 계산
        data.plot$pct.exp[data.plot$pct.exp < 0] <- NA
        data.plot$pct.exp <- data.plot$pct.exp * 100

        if(featcol %in%c("celltype_anchor","celltype_marker","celltype_profile")){
            data.plot$id <- droplevels(factor(data.plot$id, levels = obj@misc$celltype_list))
        }else{
            data.plot$id<- droplevels(factor(data.plot$id, levels = paste("Cluster", c(0:100))))
        }
        data.plot <- data.plot %>%arrange(id)
        data.plot <- data.plot %>% dplyr::rename(celltype=id)

        dot_list[[featcol]] <- data.plot
    }

    # merged dot plot data
    df <- lapply(names(dot_list), function(name) {
        dot_list[[name]] %>%
        mutate(type = name) %>%
        dplyr::select(type, everything())
    }) %>%
    bind_rows() 
    write.csv(df,paste0(path,"/meta/",ID,"_dot.csv"),row.names=FALSE)

    
    
    # -----------------------------
    # colocalization
    #------------------------------
    for( method in names(obj@misc$spe)){
        cor_df <- obj@misc$spe[[method]]
        write.csv(cor_df,paste0(path,"/meta/",ID,"_colocal_",method,".csv"))
    }


    # -----------------------------
    # circle  & heatmap   
    #------------------------------
    if("cellchat" %in% names(obj@misc)){
        for(method in names(obj@misc$cellchat)){
            
            cellchat <- obj@misc$cellchat[[method]]
            groupSize <- as.numeric(table(cellchat@idents))

            tmp <- netVisual_circle_input(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F)
            write.csv(tmp, paste0(path, "/meta/",ID,"_num_of_interactions_",method,".csv"))

            tmp <- netVisual_circle_input(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F)
            write.csv(tmp, paste0(path, "/meta/",ID,"_interaction_weights_strength_",method,".csv"))

            tmp <- netAnalysis_signalingRole_heatmap_input(cellchat, pattern = "incoming")
            write.csv(tmp, paste0(path, "/meta/",ID,"_incoming_",method,".csv"), row.names=FALSE)

            tmp <- netAnalysis_signalingRole_heatmap_input(cellchat, pattern = "outgoing")
            write.csv(tmp, paste0(path,"/meta/",ID,"_outgoing_",method,".csv"), row.names=FALSE)
            }
    }
    

    # -----------------------------
    # niche bar   
    #------------------------------
    niche_list <- obj@misc$niche
    df <- lapply(names(niche_list), function(name) {
        niche_list[[name]] %>%
        mutate(type = name) %>%
        dplyr::select(type, everything())
        }) %>%
        bind_rows()
    write.csv(df,paste0(path,"/meta/",ID,"_niche_bar.csv"),row.names=FALSE)


    # merged heatmap
    heatmap_list <- obj@misc$heatmap
    df <- lapply(names(heatmap_list), function(name) {
        heatmap_list[[name]] %>%
        mutate(type = name) %>%
        dplyr::select(type, everything())
        }) %>%
        bind_rows()
    write.csv(df,paste0(path,"/meta/",ID,"_niche_heatmap.csv"),row.names=FALSE)

}


#################################################################################
# 5. Vitessce segmentation json file
#################################################################################
library(Seurat)
library(SeuratObject)
library(SeuratData)
library(ggplot2)
library(dplyr)
library(scales)
library(R.utils)
library(jsonlite)

source("/st07/sccDB/Spatial/spatialDB/convert_to_json_format.r")

meta_df <- read.csv("/st07/sccDB/Spatial/HISSTA/pipeline/db_kobic.csv")

# -----------------------------------------------------------
# meta filtering
# -----------------------------------------------------------
meta_df <- meta_df %>% filter(tissue=="Ovary")
sample_list <- meta_df$sample_id

for( ID in sample_list){

    print(paste0("Start segmentation ",ID))

    tmp_meta <- meta_df%>%filter(sample_id==ID)
    
    # set path
    if(tmp_meta$platform == "Xenium"){
        obj_path = paste0("/st07/sccDB/Spatial/01.Xenium/14.human_ovary/01.public_data/",ID)
    } else if(tmp_meta$platform == "MERFISH"){
        obj_path = paste0("/st07/sccDB/Spatial/03.MERFISH/06.human_ovary/01.public_data/",ID)
    } else if(tmp_meta$platform == "CosMX"){
        obj_path = paste0("/st07/sccDB/Spatial/02.CosMx/09.human_ovary/",ID)
    } 

    obj <- readRDS(paste0(obj_path, "/processed/obj_processed.rds"))

    # set path
    path=paste0("/st07/kobic/data/dataset/",ID)


    if (!dir.exists(paste0(path, "/vitessce/"))) {
        dir.create(paste0(path,"/vitessce/"))
    }


    # -----------------------------------------------------------
    # cells.segmentations.json
    # -----------------------------------------------------------
    if(sum(grepl("segmentation", Boundaries(obj@images$fov)))){
        json_data <- obj@images$fov@boundaries$segmentation@polygons %>% convert_to_json_format
        json_string <- toJSON(json_data, pretty = TRUE)
        write(json_string,paste0(path, "/vitessce/",ID, "_cells.segmentations.json"))
    }
}


#################################################################################
# 6. Convert to h5ad
#################################################################################

for( ID in sample_list){
    print(paste0("Start ",ID))


    use_condaenv(condaenv = "vitessce-tutorial-env", required = TRUE)
    ad <- import("anndata", convert = FALSE)
    np <- import("numpy", convert = FALSE)

    path =paste0("/st07/kobic/data/dataset/",ID)

    obj <- readRDS(paste0(path,"/rds/obj_final.rds"))


    if (!dir.exists(paste0(path, "/vitessce/"))) {
        dir.create(paste0(path,"/vitessce/"))
    }

    meta <- cbind(obj@meta.data, obj@reductions$umap@cell.embeddings, GetTissueCoordinates(obj))
    meta <- meta %>% dplyr::select(cell_id, umap_1, umap_2, x, y, SCT_snn_res.0.3, SCT_snn_res.0.5, 
    SCT_snn_res.0.7, SCT_snn_res.1, celltype_anchor, celltype_marker, celltype_profile,
    celltype_anchor_niches, celltype_marker_niches, celltype_profile_niches)
    colnames(meta)[2:3] <- c("UMAP_1", "UMAP_2")
    colnames(meta)[4:5] <- c("X", "Y")
    meta <- reticulate::r_to_py(meta)

    spatial <- GetTissueCoordinates(obj); rownames(spatial) <- spatial$cell
    spatial <- spatial[, -3]
    spatial <- np$array(as.matrix(spatial), dtype = "float32")

    umap <- obj@reductions$umap@cell.embeddings
    umap <- np$array(as.matrix(umap), dtype = "float32")

    gene_symbol <- data.frame(gene_symbol = colnames(Matrix::t(GetAssayData(obj, layer = "data", assay ="SCT"))))
    rownames(gene_symbol) <- gene_symbol$gene_symbol
    gene_symbol <- reticulate::r_to_py(gene_symbol)

    tmp = Matrix::t(GetAssayData(obj, layer = "data", assay ="SCT"))
    tmp = as(tmp, "RsparseMatrix")
    adata = ad$AnnData(X = tmp)

    adata$obs_names <- rownames(Matrix::t(GetAssayData(obj, layer = "data", assay ="SCT")))
    adata$var_names <- colnames(Matrix::t(GetAssayData(obj, layer = "data", assay ="SCT")))

    adata$obs <- meta
    adata$var <- gene_symbol
    adata$obsm['spatial'] <- spatial
    adata$obsm['X_umap'] <- umap

    adata$write_h5ad(paste0(path, "/vitessce/obj_final.h5ad"))
}


# segmentation이 없는 경우는 polygon.py로 가공






############################################
# 7. Gene panel
############################################
library(dplyr)

for( ID in sample_list){
    
    tmp_meta <- meta_imsi%>%filter(sample_id==ID)
    
    # set path
    if(tmp_meta$platform == "Xenium"){
        obj_path = paste0("/st07/sccDB/Spatial/01.Xenium/14.human_ovary/01.public_data/",ID)
    } else if(tmp_meta$platform == "MERFISH"){
        obj_path = paste0("/st07/sccDB/Spatial/03.MERFISH/06.human_ovary/01.public_data/",ID)
    } else if(tmp_meta$platform == "CosMX"){
        obj_path = paste0("/st07/sccDB/Spatial/02.CosMx/09.human_ovary/",ID)
    } 


    # set kobic_path
    path =paste0("/st07/kobic/data/dataset/",ID)


    if (!dir.exists(paste0(path, "/panel/"))) {
        dir.create(paste0(path,"/panel/"))
    }


    #-------------------------------------------------------------
    # 1.processing
    #-------------------------------------------------------------
    print(paste0(ID))

    # Load data
    if(tmp_meta$platform=="Xenium"){
        gene <- read.delim(paste0(obj_path,"/rawdata/cell_feature_matrix/features.tsv.gz"),header=F)
        gene <- gene$V2 %>%as.data.frame
        colnames(gene) <-"gene"
        write.csv(gene,paste0(path,"/panel/gene_panel.csv"))
    } else if(tmp_meta$platform == "MERFISH"){
        gene <- read.csv(paste0(obj_path,"/rawdata/cell_by_gene.csv"),nrow=1)
        gene <- t(gene)
        gene <- rownames(gene)
        gene2 <- gene[!grepl("^Blank\\.\\d+$", gene) ]

        gene2 <- as.data.frame(gene2[-1])
        colnames(gene2) <-"gene"

        write.csv(gene2,paste0(path,"/panel/gene_panel.csv"))

    } else if(tmp_meta$platform == "CosMX"){
        count.file <- list.files(paste0(obj_path, "/rawdata"))
        count.file <- count.file[grepl('[_a-zA-Z0-9]*_exprMat_file.csv', count.file)]
        count.file <- paste0(obj_path, "/rawdata/", count.file)

        gene <- read.csv(count.file,nrow=1)
        gene <- t(gene)
        gene <- rownames(gene)
        gene2 <- gene[!grepl("^NegPrb\\d+$", gene) ]

        gene2 <- as.data.frame(gene2[c(-1,-2)])
        colnames(gene2) <-"gene"

        write.csv(gene2,paste0(path,"/panel/gene_panel.csv"))


    }
}

#################################################################################
# 8. Draw plot
#################################################################################
library(ggplot2)
library(magick)
library(cowplot)
mycols <- read.csv("/st07/sccDB/Spatial/HISSTA/pipeline/color_list.csv")
meta_df <- read.csv("/st07/sccDB/Spatial/HISSTA/pipeline/db_kobic.csv")
sample_list <- meta_df$sample_id

ds=0.7

for( ID in sample_list){
    
    print(paste0("Start ",ID))

    tmp_meta <- meta_df%>%filter(sample_id==ID)

    # set path
    path =paste0("/st07/kobic/data/dataset/",ID)

    if (!dir.exists(paste0(path, "/plot_image/"))) {
        dir.create(paste0(path,"/plot_image/"))
    }


    if (!dir.exists(paste0(path, "/rds/"))) {
        dir.create(paste0(path,"/rds/"))
    }


    obj <- readRDS(paste0(path, "/rds/obj_raw.rds"))

    
    # Imgdim plot의 비율 조정
    coord <- GetTissueCoordinates(obj)
    ratio <- c(range(coord$x)[2] - range(coord$x)[1], range(coord$y)[2] - range(coord$y)[1])
    max_value <- 14; 
    scale_factor <- max_value / max(ratio); scaled_ratio <- ratio * scale_factor


    # n_Count ImageFeaturePlot
    obj$imsi <- " "
    Idents(obj)<- obj$imsi
    feat <- colnames(obj@meta.data)
    feat <- feat[grepl("nCount", feat) & (grepl("(Xenium|Nanostring|Vizgen)", feat))] #nCount column명 추출

    ImageFeaturePlot(obj, features = feat, size = ds, fov = "fov",scale="none",
            dark.background = FALSE)+
            scale_fill_viridis()+
            theme(legend.position = "none",
            plot.title = element_text(size=0),
            plot.margin = unit(c(0, 1, 1, 1), "cm"))
    ggsave(filename = paste0(path,"/plot_image/",ID,"_Imgfeat_nCount_nolegend.png"),
        width = scaled_ratio[2], height = scaled_ratio[1])
    
    image <- image_read(paste0(path,"/plot_image/",ID,"_Imgfeat_nCount_nolegend.png"))
    image_rotated <- image_rotate(image, 90)
    image_write(image_rotated, path = paste0(path,"/plot_image/",ID,"_Imgfeat_nCount_rotated.png"))

    if(tmp_meta$platform=="Xenium"){
        wd=2.9
    } else if(tmp_meta$platform=="CosMX"){
        wd=3.43
    } else{
        wd=2.75
    }

    plot_legend <- ImageFeaturePlot(obj, features = feat, size = ds, fov = "fov", scale="none",
                                    dark.background = FALSE) +
                scale_fill_viridis() +
                theme(legend.position = "right",
                    legend.text=element_text(size=18) ,legend.title=element_text(size=24),
                legend.key.height = unit(1.5, 'cm'),legend.key.width = unit(1.2, 'cm'),
                plot.title = element_text(size=0))
    ggsave(filename = paste0(path,"/plot_image/",ID,"_Imgfeat_nCount_legend.png"),
        plot = plot_legend, width = wd, height = scaled_ratio[2])  


    image_rotated <- image_read(paste0(path,"/plot_image/",ID,"_Imgfeat_nCount_rotated.png"))
    image_legend <- image_read(paste0(path,"/plot_image/",ID,"_Imgfeat_nCount_legend.png"))


    # 회전된 이미지에 범례 추가
    combined_image <- image_append(c(image_rotated, image_legend), stack = FALSE)
    image_write(combined_image, path = paste0(path,"/plot_image/",ID,"_Imgfeat_nCount.png"))
    

    
    # n_Feature vlnplot
    feat <- colnames(obj@meta.data)
    feat <- feat[grepl("nFeature", feat) & (grepl("(Xenium|Nanostring|Vizgen)", feat))]
 
    ImageFeaturePlot(obj, features = feat, size = ds, fov = "fov",scale="none",
            dark.background = FALSE)+
            scale_fill_viridis()+
            theme(legend.position = "none",
            plot.title = element_text(size=0),
            plot.margin = unit(c(0, 1, 1, 1), "cm"))
    ggsave(filename = paste0(path,"/plot_image/",ID,"_Imgfeat_nFeature_nolegend.png"),
        width = scaled_ratio[2], height = scaled_ratio[1])
    
    image <- image_read(paste0(path,"/plot_image/",ID,"_Imgfeat_nFeature_nolegend.png"))
    image_rotated <- image_rotate(image, 90)
    image_write(image_rotated, path = paste0(path,"/plot_image/",ID,"_Imgfeat_nFeature_rotated.png"))

    if(tmp_meta$platform=="Xenium"){
        wd=3.15
    } else if(tmp_meta$platform=="CosMX"){
        wd=3.5
    } else{
        wd=3
    }

    plot_legend <- ImageFeaturePlot(obj, features = feat, size = ds, fov = "fov", scale="none",
                                    dark.background = FALSE) +
                scale_fill_viridis() +
                theme(legend.position = "right",
                    legend.text=element_text(size=18) ,legend.title=element_text(size=24),
                legend.key.height = unit(1.5, 'cm'),legend.key.width = unit(1.2, 'cm'),
                plot.title = element_text(size=0))
    ggsave(filename = paste0(path,"/plot_image/",ID,"_Imgfeat_nFeature_legend.png"),
        plot = plot_legend, width = wd, height = scaled_ratio[2])  # 여기서 width는 범례의 너비로 조정


    image_rotated <- image_read(paste0(path,"/plot_image/",ID,"_Imgfeat_nFeature_rotated.png"))
    image_legend <- image_read(paste0(path,"/plot_image/",ID,"_Imgfeat_nFeature_legend.png"))


    # 회전된 이미지에 범례 추가
    combined_image <- image_append(c(image_rotated, image_legend), stack = FALSE)
    image_write(combined_image, path = paste0(path,"/plot_image/",ID,"_Imgfeat_nFeature.png"))
}


#### cluster, cell type, niche

for( ID in sample_list){
    
    print(paste0("Start ",ID))

    tmp_meta <- meta_df%>%filter(sample_id==ID)

    # set path
    path =paste0("/st07/kobic/data/dataset/",ID)

    if (!dir.exists(paste0(path, "/plot_image/"))) {
        dir.create(paste0(path,"/plot_image/"))
    }


    if (!dir.exists(paste0(path, "/rds/"))) {
        dir.create(paste0(path,"/rds/"))
    }


    obj <- readRDS(paste0(path, "/rds/obj_final.rds"))

    coord <- GetTissueCoordinates(obj)
    ratio <- c(range(coord$x)[2] - range(coord$x)[1], range(coord$y)[2] - range(coord$y)[1])
    max_value <- 14; 
    scale_factor <- max_value / max(ratio); scaled_ratio <- ratio * scale_factor
   

    feat <- colnames(obj@meta.data)
    feat <- c(feat[grepl("SCT_snn_", feat)])

    # --------------------------------------
    # clustering image dimplot
    #---------------------------------------

    for( featcol in feat){
        tmp_cols <- mycols %>% filter(type %in% unique(obj@meta.data[, featcol]))

        ImageDimPlot(obj, fov = "fov", group.by=featcol,
                size = ds,cols = tmp_cols$color, dark.background = FALSE)+ 
                theme(legend.position = "none",
                plot.title = element_text(size=0),
                plot.margin = unit(c(0, 1, 1, 1), "cm"))
        ggsave(filename = paste0(path,"/plot_image/",ID,"_Imgdim_", featcol,"_nolegend.png"),width = scaled_ratio[2]+4, height = scaled_ratio[1])



        image <- image_read(paste0(path,"/plot_image/",ID,"_Imgdim_", featcol,"_nolegend.png"))
        image_rotated <- image_rotate(image, 90)
        image_write(image_rotated, path = paste0(path,"/plot_image/",ID,"_Imgdim_", featcol,"_rotated.png"))


        if(length(table(obj@meta.data[, featcol]))<21){
            wd=3
        } else{
            wd=4.3
        }
        # wd=6.3 #40개 이상
        # wd=2.7 res 1

        plot_legend <- ImageDimPlot(obj, fov = "fov", group.by=featcol,
                    size = ds,cols = tmp_cols$color, dark.background = FALSE) +
                    theme(legend.position = "right",
                        legend.text=element_text(size=20) ,legend.title=element_text(size=24),
                    legend.key.height = unit(1, 'cm'),legend.key.width = unit(1.2, 'cm'),
                    plot.title = element_text(size=0))
        ggsave(filename = paste0(path,"/plot_image/",ID,"_Imgdim_", featcol,"_legend.png"),
            plot = plot_legend, width = wd, height = scaled_ratio[2]+4)


        image_rotated <- image_read(paste0(path,"/plot_image/",ID,"_Imgdim_", featcol,"_rotated.png"))
        image_legend <- image_read(paste0(path,"/plot_image/",ID,"_Imgdim_", featcol,"_legend.png"))


        # 회전된 이미지에 범례 추가
        combined_image <- image_append(c(image_rotated, image_legend), stack = FALSE)
        image_write(combined_image, path = paste0(path,"/plot_image/",ID,"_Imgdim_",featcol,".png"))
    }


    # --------------------------------------
    # celltype image dimplot
    #---------------------------------------
    feat <-  c("celltype_anchor", "celltype_marker", "celltype_profile") #clustering, celltype annotation column명 추출

    for( featcol in feat){
        tmp_cols <- names(obj@misc$celltype_list[obj@misc$celltype_list %in% names(table(obj@meta.data[,featcol]))])

        ImageDimPlot(obj, fov = "fov", group.by=featcol,
                size = ds,cols = tmp_cols, dark.background = FALSE) +
                theme(legend.position = "none",
                plot.title = element_text(size=0),
                plot.margin = unit(c(0, 1, 1, 1), "cm"))
        ggsave(filename = paste0(path,"/plot_image/",ID,"_Imgdim_", featcol,"_nolegend.png"),width = scaled_ratio[2]+4, height = scaled_ratio[1])


        image <- image_read(paste0(path,"/plot_image/",ID,"_Imgdim_", featcol,"_nolegend.png"))
        image_rotated <- image_rotate(image, 90)
        image_write(image_rotated, path = paste0(path,"/plot_image/",ID,"_Imgdim_", featcol,"_rotated.png"))

        if(featcol=="celltype_anchor"){
            title="Anchor-based Cell Type"
            wd=4.15
        } else if(featcol=="celltype_marker"){
            title="Marker-based Cell Type"
            wd=4.15
        } else{
            title="Profile-based Cell Type"
            wd=4
        }

        #wd=5.8 # Uterine
        #wd=6.8 #2줄

        plot_legend <- ImageDimPlot(obj, fov = "fov", group.by=featcol,
                    size = ds,cols = tmp_cols, dark.background = FALSE) +
                    theme(legend.position = "right",
                        legend.text=element_text(size=20) ,legend.title=element_text(size=24),
                    legend.key.height = unit(1, 'cm'),legend.key.width = unit(1.2, 'cm'),
                    plot.title = element_text(size=0))+
                    labs(fill=paste(title))
        ggsave(filename = paste0(path,"/plot_image/",ID,"_Imgdim_", featcol,"_legend.png"),
            plot = plot_legend, width = wd, height = scaled_ratio[2]+4)


        image_rotated <- image_read(paste0(path,"/plot_image/",ID,"_Imgdim_", featcol,"_rotated.png"))
        image_legend <- image_read(paste0(path,"/plot_image/",ID,"_Imgdim_", featcol,"_legend.png"))


        # 회전된 이미지에 범례 추가
        combined_image <- image_append(c(image_rotated, image_legend), stack = FALSE)
        image_write(combined_image, path = paste0(path,"/plot_image/",ID,"_Imgdim_",featcol,".png"))
        
    }


    # --------------------------------------
    # kobic_plot 4. niche image dimplot
    #---------------------------------------
    feat <-  names(obj@misc$niche)
    coord <- GetTissueCoordinates(obj)
    ratio <- c(range(coord$x)[2] - range(coord$x)[1], range(coord$y)[2] - range(coord$y)[1])

    ratio <- ratio/ratio[2]

    for( featcol in feat){

        feat_niches <- paste0(featcol,"_niches")
        obj$niches <- obj@meta.data[,feat_niches]
        obj$niches<- factor(obj$niches ,levels = paste("Niche", c(1:5)))

        tmp_cols <- names(obj@misc$celltype_list[obj@misc$celltype_list %in% names(table(obj@meta.data[,featcol]))])

        niche.plot <- ImageDimPlot(obj, group.by = "niches", size = ds, dark.background = F) + 
            ggtitle(paste0(feat_niches)) +
            scale_fill_manual(values = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B"))+
            theme(legend.position = "none",
                plot.title = element_text(size=0),
            plot.margin = margin(10, 10, 10, 10))
        ggsave(filename = paste0(path,"/plot_image/",ID,"_Imgdim_", feat_niches,"_nolegend.png"),width = scaled_ratio[2]+4, height = scaled_ratio[1])


        image <- image_read(paste0(path,"/plot_image/",ID,"_Imgdim_", feat_niches,"_nolegend.png"))
        image_rotated <- image_rotate(image, 90)
        image_write(image_rotated, path = paste0(path,"/plot_image/",ID,"_Imgdim_", feat_niches,"_rotated.png"))


        plot_legend <- ImageDimPlot(obj, group.by = "niches", size = ds, dark.background = F) + 
            ggtitle(paste0(feat_niches)) +
            scale_fill_manual(values = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B")) +
                    theme(legend.position = "right",
                        legend.text=element_text(size=20) ,legend.title=element_text(size=24),
                    legend.key.height = unit(1, 'cm'),legend.key.width = unit(1.2, 'cm'),
                    plot.title = element_text(size=0))
        ggsave(filename = paste0(path,"/plot_image/",ID,"_Imgdim_", feat_niches,"_legend.png"),
            plot = plot_legend, width = 1.85, height = scaled_ratio[2]+4)


        image_rotated <- image_read(paste0(path,"/plot_image/",ID,"_Imgdim_", feat_niches,"_rotated.png"))
        image_legend <- image_read(paste0(path,"/plot_image/",ID,"_Imgdim_", feat_niches,"_legend.png"))


        # 회전된 이미지에 범례 추가
        combined_image <- image_append(c(image_rotated, image_legend), stack = FALSE)
        image_write(combined_image, path = paste0(path,"/plot_image/",ID,"_Imgdim_",feat_niches,".png"))
        
        
        # cell type + niche
        image_ct <- image_read(paste0(path,"/plot_image/",ID,"_Imgdim_", featcol,".png"))
        image_niche <- image_read(paste0(path,"/plot_image/",ID,"_Imgdim_", feat_niches,".png"))

        combined_image <- image_append(c(image_ct, image_niche), stack = FALSE)
        image_write(combined_image, path = paste0(path,"/plot_image/",ID,"_Imgdim_",feat_niches,".png"))
    

        if(ratio[1] > 3){
            combined_image <- image_append(c(image_ct, image_niche), stack = TRUE)
            image_write(combined_image, path = paste0(path,"/plot_image/",ID,"_Imgdim_",feat_niches,".png"))
        } else{
            combined_image <- image_append(c(image_ct, image_niche), stack = FALSE)
            image_write(combined_image, path = paste0(path,"/plot_image/",ID,"_Imgdim_",feat_niches,".png"))
        }
    }
}









# -----------------------------------------------------------
# polygon O
# -----------------------------------------------------------
# conda activate vitessce-tutorial-env

import pandas as pd
import anndata as ad
import numpy as np
import scipy.cluster
from scipy import sparse
from vitessce.data_utils import (
    to_diamond,
    rgb_img_to_ome_tiff,
    optimize_adata,
    to_dense,
    VAR_CHUNK_SIZE,
    sort_var_axis,
)
from os.path import join


df = pd.read_csv("/st07/sccDB/Spatial/HISSTA/pipeline/db_kobic.csv")
df.columns

filtered_df = df[df['tissue'].isin(['Ovary'])]
datasets = filtered_df.sample_id


for dataset in datasets:
print(dataset)
input_path = f"/st07/kobic/data/dataset/{dataset}/vitessce/obj_final.h5ad"
output_path = f"/st07/kobic/data/dataset/{dataset}/vitessce/obj_final.zarr"
adata = ad.read_h5ad(input_path)
adata = optimize_adata(
    adata,
    obs_cols=['SCT_snn_res.0.3', 'SCT_snn_res.0.5', 'SCT_snn_res.0.7', 'SCT_snn_res.1', 
    'celltype_anchor', 'celltype_marker', 'celltype_profile', 
    'celltype_anchor_niches', 'celltype_marker_niches', 'celltype_profile_niches'],
    obsm_keys=["X_umap", "spatial"],
    optimize_X=True,
    to_dense_X=True
    )
adata.write_zarr(output_path, chunks=[adata.shape[0], VAR_CHUNK_SIZE])
















# -----------------------------------------------------------
# polygon X
# -----------------------------------------------------------
# conda activate vitessce-tutorial-env

import pandas as pd
import anndata as ad
import numpy as np
import scipy.cluster
from scipy import sparse
from vitessce.data_utils import (
    to_diamond,
    rgb_img_to_ome_tiff,
    optimize_adata,
    to_dense,
    VAR_CHUNK_SIZE,
    sort_var_axis,
)
from os.path import join

datasets=['MER089']

for dataset in datasets:
    print(dataset)
    input_path = f"/st07/kobic/data/dataset/{dataset}/vitessce/obj_final2.h5ad"
    output_path = f"/st07/kobic/data/dataset/{dataset}/vitessce/obj_final2.zarr"
    adata = ad.read_h5ad(input_path)
    num_cells = adata.obs.shape[0]
    scale_factor = 1 #um로 변환
    adata.obsm['spatial'] = (adata.obsm['spatial'] * scale_factor)
    adata.obsm['segmentations'] = np.zeros((num_cells, 4, 2))
    radius = 10
    for i in range(num_cells):
        adata.obsm['segmentations'][i, :, :] = to_diamond(adata.obsm['spatial'][i, 0], adata.obsm['spatial'][i, 1], radius)
    adata = optimize_adata(
        adata,
        obs_cols=['SCT_snn_res.0.3', 'SCT_snn_res.0.5', 'SCT_snn_res.0.7', 'SCT_snn_res.1', 
        'celltype_anchor', 'celltype_marker', 'celltype_reference', 
        'celltype_anchor_niches', 'celltype_marker_niches', 'celltype_reference_niches'],
        obsm_keys=["X_umap", "spatial", "segmentations"],
        optimize_X=True,
        to_dense_X=True
        )
    adata.write_zarr(output_path, chunks=[adata.shape[0], VAR_CHUNK_SIZE])