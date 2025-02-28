# HISSTA
HISSTA: a Human In Situ Single-cell Transcriptome Atlas

https://kbds.re.kr/hissta/

Article: [] (논문 링크 수정하기)

HISSTA provides several advanced tools for spatial context analyses to aid users in gaining scientific insights:
1. Cell type annotation. Three ways of cell type annotation are provided: 1) marker-based using SCINA, 2) profile-based using Insitutype, and 3) a hybrid method combining both approaches. Results are used in the subsequent neighborhood analyses and visualized interactively with Vitessce.
2. Spatial colocalization analysis. hoodscanR is used to identify neighborhood cells and calculate colocalization correlations between different cell types.
3. Spatial cellular communication. CellChat v2 is employed to identify spatially proximal cell-cell communication, visualized with circular plots.
4. Spatial niche analysis. The BuildNicheAssay function of Seurat identifies cell niches based on 30 neighborhood cells. We further show the cell type composition and niche-specific pathways.


![Graphical Abstract](images/0_graphical_abstract.png)
