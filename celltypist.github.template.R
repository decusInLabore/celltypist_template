####################################################################################################
### Run CellTypist on a Seurat object
####################################################################################################


## Environment for this task ##

You can either use your personal R/python3 setup (no gurantees) or you can use the singularity container
with R/Seurat5 and python 3 in it.

The environment consists of two parts: a singularity container setting up the LINUX environment and
a renv part for the R-environment.

## Cloning the github repo ##
Start by cloning the github repo into a place on nemo, e.g. somewhere in your project.

git clone git@github.com:FrancisCrickInstitute/celltypist.git

## cd into the celltypist directory
cd celltypist

On Nemo, run at the command line the following:

ml Singularity/3.6.4
singularity shell --cleanenv --bind /nemo:/nemo,/camp:/camp /flask/apps/containers/all-singularity-images/r431.ubuntu.22.04.sif;

# Start R by typing
R

####################################################################################
## Now we continue in the singulariy container R environment or your R environment##
####################################################################################

# If you're using the containerized R version you can restore the R-environment by running the following
if (!require("renv")){
  install.packages("renv")
}

renv::restore("renv.lock")

## Load your Seurat object. In this example the Seurat object is called OsC ##

Once you have loaded your Seurat object, rename it as follows:
OsC <- [yourSeuratObject]

## First the remove all unnecessary items from the Seurat object ##

OsCsub <- Seurat::DietSeurat(
  OsC,
  # counts = TRUE, # so, raw counts save to adata.layers['counts']
  layers = c("counts", "data"), # so, log1p counts save to adata.X when scale.data = False, else adata.layers['data']
  scale.data = FALSE, # if only scaled highly variable gene, the export to h5ad would fail. set to false
  features = rownames(OsC), # export all genes, not just top highly variable genes
  assays = "RNA",
  dimreducs = c("pca","umap"),
  graphs = c("RNA_nn", "RNA_snn"), # to RNA_nn -> distances, RNA_snn -> connectivities
  misc = TRUE
)

# sce <- Seurat::as.SingleCellExperiment(OsCsub, assay = "RNA")
sce <- Seurat::as.SingleCellExperiment(OsCsub)

#unlink("../../../../data/MYC.h5ad")

# Option A (faster)
outfile <- "conversted.seurat.object.h5ad"

if (!require("sceasy")){
  renv::install("cellgeni/sceasy")
}

## Create a h5ad file to read into python
sceasy::convertFormat(
  sce, 
  from = "sce",
  to = "anndata",
  outFile = outfile
)

## Now we can exit R and continue in python
q()

####################################################################################
## Running celltypist in python3                                                  ##
####################################################################################

## Before we start python we need to make sure the celltypist package is installed
pip install celltypist
pip install scanpy

## To start python in the singularity environment, type python3 at the command line.
python3


## Now in python3 do:
import scanpy as sc
import os
import celltypist
from celltypist import models


# Load data
adata =  sc.read_h5ad("conversted.seurat.object.h5ad")

adata.layers['counts'] = adata.X

## Some basic transformation
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

## Download the model you'd like to use from https://www.celltypist.org/models
## Alternatively you can create and save your own model

modSel = 'Healthy_Mouse_Liver.pkl'
os.path.isfile(modSel)


## Now we run celltypist

predictions = celltypist.annotate(adata, model = modSel , majority_voting = True)

## If pca dims are not sufficient ##
n_comps_required = 50  # Adjust this value as needed

# Rerun PCA with adjusted n_comps
sc.pp.pca(adata, n_comps=n_comps_required)
adata.obsm['X_pca'].shape
predictions = celltypist.annotate(adata, model = modSel , majority_voting = True)

predictions.predicted_labels

## Optional dotplot with original clusters
dotplot = celltypist.dotplot(predictions, use_as_reference = 'seurat_clusters', use_as_prediction = 'predicted_labels')

import matplotlib.pyplot as plt
plt.savefig('report_figures/celltypist_dotplot.pdf', format='pdf')


# Get an `AnnData` with predicted labels embedded into the cell metadata columns.
adata = predictions.to_adata()

adata .obs["clusterName"]= adata.obs["majority_voting"]

## optional plot
sc.pl.umap(adata, color = ['clusterName', 'majority_voting'], legend_loc = 'on data')

fn = "test1.pdf"
import matplotlib.pyplot as plt
plt.savefig(fn, bbox_inches='tight', dpi = 50)

## Now we make an outputfile so we can add this information to the original Seurat object
colVec = ['cellID', 'clusterName','predicted_labels', 'conf_score']

df = adata.obs[colVec]

## Get HEX colors
cluster_colors = adata.uns['clusterName_colors']


# Convert RGB colors to hex
import matplotlib.colors as mcolors
import pandas as pd

cluster_hex_colors = [mcolors.rgb2hex(color) for color in cluster_colors]


df_clusters = pd.DataFrame({
'clusterName': adata.obs['clusterName'].cat.categories,
'clusterColor': cluster_hex_colors
})

merged_df = pd.merge(df, df_clusters, on='clusterName')

fnOut = "celltypist.annotation.csv"
merged_df.to_csv(fnOut)

## now we are done in python and can return to the command line
exit()

#####################################################################################
## Now we return to R to add the celltypist annotation to the Seurat object        ##
#####################################################################################

# At the command line start R
R

## Back to R

## Load the original Seurat object and rename it to OsC if you want to run this script directly
OsC <- [your Seurat object]

df <- read.csv("celltypist.annotation.csv")
df$X <- NULL
row.names(df) <- df$cellID
df$cellID <- NULL
df$over_clustering = NULL

## Check for special characters.
# Load necessary libraries
library(dplyr)
library(stringr)

# Define a function to clean column values
clean_column <- function(column) {
column %>%
str_replace_all(" ", "_") %>%  # Replace spaces with underscores
str_replace_all("[^[:alnum:]_]", "")  # Remove non-standard characters
}

# Example usage with a data frame `df` and a column `columnName`
df <- df %>%
mutate(clusterName = clean_column(clusterName)) %>%
mutate(predicted_labels = clean_column(predicted_labels))

df$clusterName <- gsub("__", "_", df$clusterName)

unique(df$clusterName)
df$conf_score <- round(df$conf_score, 3)


OsC <- biologicToolsSC::addDf2seuratMetaData(obj = OsC, dfAdd = df)

## Save annotated Seurat object
save(OsC,
  file = paste0(
  Obio@parameterList$localWorkDir,
  Obio@parameterList$project_id,
  ".Seurat.Robj"
  )
)

file = paste0(
  Obio@parameterList$localWorkDir,
  Obio@parameterList$project_id,
  ".SeuratV5.obj.Rds"
)

saveRDS(
  object = OsC,
  file = file
)



######### Optional ##################
## Training

adata = sc.read('celltypist_demo_folder/gut_cell_atlas_James.h5ad', backup_url = 'https://cellgeni.cog.sanger.ac.uk/gutcellatlas/Colon_cell_atlas.h5ad')
data_James.obs.cell_type.unique()
sampled_cell_index = celltypist.samples.downsample_adata(adata_Elmentaite, mode = 'each', n_cells = 500, by = 'Integrated_05', return_index = True)

t_start = time.time()
model_fs = celltypist.train(adata_Elmentaite[sampled_cell_index], 'Integrated_05', n_jobs = 10, max_iter = 5, use_SGD = True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")
gene_index = np.argpartition(np.abs(model_fs.classifier.coef_), -100, axis = 1)[:, -100:]
gene_index = np.unique(gene_index)
# Add `check_expression = False` to bypass expression check with only a subset of genes.
t_start = time.time()
model = celltypist.train(adata_Elmentaite[sampled_cell_index, gene_index], 'Integrated_05', check_expression = False, n_jobs = 10, max_iter = 100)
t_end = time.time()
print(f"Time elapsed: {(t_end - t_start)/60} minutes")
model.write('celltypist_demo_folder/model_from_Elmentaite_2021.pkl')
t_start = time.time()
predictions = celltypist.annotate(adata_James, model = 'celltypist_demo_folder/model_from_Elmentaite_2021.pkl')
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")
t_start = time.time()
predictions = celltypist.annotate(adata_James, model = 'celltypist_demo_folder/model_from_Elmentaite_2021.pkl')
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")

## Create a model 
max_iter = 300 
150 genes

############# done optional ##############################








```{r, eval=TRUE, echo=F, results=FALSE, message = FALSE, warning = FALSE}
tempDir <- "../../../../workdir/temp/"

if (!dir.exists(tempDir)){
  dir.create(tempDir, recursive = T)
}

if (exists("whiteListWorkspace")){
  rm(list = setdiff(ls(), whiteListWorkspace))
}

FN <- "../../../../workdir/temp/temp.workspace.RData"
save.image(FN)

```

## Create loop file

if (!require("loupeR")){
  renv::install("10XGenomics/loupeR")
}

# Gene Expression RNA assay
assay <- OsC[["RNA"]]

# get counts matrix from either the old or newer formats of assay
counts <- loupeR::counts_matrix_from_assay(assay)

# convert the count matrix, clusters, and projections into a Loupe file

## Preselect cluster slot
OsC@meta.data$clusterName <- gsub("_", "", OsC@meta.data$clusterName)

Seurat::Idents(OsC) <- "clusterName"

L1 <- loupeR::select_clusters(OsC)
L1 <- L1[c("clusterName", "sampleName")]

P1 <- loupeR::select_projections(OsC)
P1 <- P1["umap"]

loupeR::create_loupe(
  count_mat = counts,
  clusters = L1,
  projections = P1,
  output_name = "loupe_object"
)


###################################
## Trying cellhint


import cellhint

cellhint.integrate(
  adata = adata, 
  cell_type = 'clusterName'
)

sc.pl.umap(adata, color = 'clusterName')

cellhint.integrate(
  adata = adata, 'Dataset', 'Curated_annotation')

sc.pl.umap(adata, color = 'clusterName')

import matplotlib.pyplot as plt
fn = "../../../../html_local/report_figures/test1.pdf"
plt.savefig(fn, dpi = 50)

sc.pl.umap(adata, color = 'sampleName')


fn = "../../../../html_local/report_figures/test2.pdf"
plt.savefig(fn, dpi = 50)

