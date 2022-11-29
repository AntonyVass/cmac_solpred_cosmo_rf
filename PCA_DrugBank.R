library(dplyr)
library(plotly)

# ----- files in
f_DrugBank <- "path/to/datadir/drugbank2022_moe.txt"
f_MOEsolu <- "path/to/datadir/dataset_solutes_moe.txt"
f_MOEsolv <- "path/to/datadir/dataset_solvents_moe.txt"
data_DrugBank <- read.csv(f_DrugBank, header = TRUE, sep = "\t")
data_MOEsolu <- read.csv(f_MOEsolu, header = TRUE, sep = "\t")
data_MOEsolv <- read.csv(f_MOEsolv, header = TRUE, sep = "\t")

# ----- prep DrugBank data
data_DrugBank_clean <- data_DrugBank[,-c(1:4)]
data_DrugBank_clean_noNA <- data_DrugBank_clean[ , colSums(is.na(data_DrugBank_clean)) == 0]
data_DrugBank_clean_noZV <- data_DrugBank_clean_noNA[,apply(data_DrugBank_clean_noNA, 2, var, na.rm=TRUE) != 0]

# ----- prep our solute and solvent datasets
solutes <- 
  data_MOEsolu[match(unique(data_MOEsolu$solute_name), data_MOEsolu$solute_name),] %>% 
  rename(compound = solute_name) %>% 
  select(-c(1:3,5:20))
solvents <- 
  data_MOEsolv[match(unique(data_MOEsolv$solvent_name), data_MOEsolv$solvent_name),] %>%
  rename(compound = solvent_name) %>% 
  select(-c(1:2,4:20))
all_compounds <- rbind(solutes,solvents)
all_compounds <- cbind(Tag = c(rep("solute",75),rep("solvent",49)), all_compounds)
all_compounds_clean <- all_compounds[,-c(1:2)]

# ----- match features to those in DrugBank data
all_compounds_clean_matched <- all_compounds_clean %>% select(colnames(data_DrugBank_clean_noZV))

# ----- PCA on DrugBank 
pca <- prcomp(data_DrugBank_clean_noZV, center = T, scale = T) 

# ----- project our dataset into the PC space
pca_dataset <- scale(all_compounds_clean_matched, pca$center, pca$scale) %*% pca$rotation 

# ----- get loadings
sdev <- summary(pca)[["sdev"]]
pca_loadings <- pca$rotation # init size
for (i in seq(sdev)) {
  pca_loadings[, i] <- pca$rotation[, i] * sdev[i]
}

# ----- plot
pca__p1 <- 
  plot_ly(
    data = data.frame(pca$x),
    type = 'scatter',
    #type = 'scatter3d',
    mode = 'markers',
    x = ~ PC1,
    y = ~ PC2,
    #z = ~ PC3,
    opacity = 0.6,
    marker = list(size = 6),
    name = "DrugBank (2022)",
    text = data_DrugBank$DATABASE_ID,
    hovertemplate = paste('<i>%{text}</i>:',
                          '<br><b>PC1</b>: %{x}<br>',
                          '<b>PC2</b>: %{y}')
  ) %>%
  
  add_trace(
    data = data.frame(pca_dataset)[all_compounds$Tag == "solute",],
    x = ~ PC1,
    y = ~ PC2,
    #z = ~ PC3,
    mode = 'markers',
    opacity = 0.8,
    marker = list(size = 14),
    name = "this study - solutes",
    text = all_compounds$compound[all_compounds$Tag == "solute"],
    hovertemplate = paste('<i>%{text}</i>:',
                          '<br><b>PC1</b>: %{x}<br>',
                          '<b>PC2</b>: %{y}')
  ) %>%
  
  # add_trace(
  #   data = data.frame(pca_dataset)[all_compounds$Tag == "solvent",],
  #   x = ~ PC1,
  #   y = ~ PC2,
  #   #z = ~ PC3,
  #   mode = 'markers',
  #   opacity = 0.8,
  #   marker = list(size = 14),
  #   name = "this study - solvents",
  #   text = all_compounds$compound[all_compounds$Tag == "solvent"],
  #   hovertemplate = paste('<i>%{text}</i>:',
  #                         '<br><b>PC1</b>: %{x}<br>',
  #                         '<b>PC2</b>: %{y}')
  # ) %>%
  
  layout(
    xaxis = list(title = "PC1"),
    yaxis = list(title = "PC2"),
    legend = list(x = 100, y = 0.5)
  ) 
pca__p1

# ----- add feature vectors for bi-plot
# (we manually picked these two - select any feature(s) you want in the same way)
my_feature1 <- "logP.o.w." 
my_feature2 <- "Weight"
pca__p2 <- pca__p1 %>%
  add_segments(
    x = 0,
    xend = pca_loadings[my_feature1, 1] * 10,
    y = 0,
    yend = pca_loadings[my_feature1, 2] * 10,
    line = list(color = 'black'),
    inherit = FALSE,
    showlegend = FALSE
  ) %>%
  add_annotations(
    x = pca_loadings[my_feature1, 1] * 12,
    y = pca_loadings[my_feature1, 2] * 12,
    ax = 0,
    ay = 0,
    text = paste('<b>', "log P", '</b>'),
    font = list(size = 18),
    xanchor = 'center',
    yanchor = 'bottom'
  ) %>%
  add_segments(
    x = 0,
    xend = pca_loadings[my_feature2, 1] * 10,
    y = 0,
    yend = pca_loadings[my_feature2, 2] * 10,
    line = list(color = 'black'),
    inherit = FALSE,
    showlegend = FALSE
  ) %>%
  add_annotations(
    x = pca_loadings[my_feature2, 1] * 12,
    y = pca_loadings[my_feature2, 2] * 12,
    ax = 0,
    ay = 0,
    text = paste('<b>', "molecular weight", '</b>'),
    font = list(size = 18),
    xanchor = 'center',
    yanchor = 'bottom'
  )
pca__p2

# ----- scree plot
pca__explainedVar <- summary(pca)$importance[2,1:80]
pca__explainedVar_cumul <- summary(pca)$importance[3,1:80]
pca__p3 <- plot_ly(
  name = "% Variance Explained",
  y = pca__explainedVar * 100,
  x = paste0("PC", seq(pca__explainedVar)),
  type = 'bar'
) %>%
  add_trace(
    name = "Cumulative % Variance Explained",
    yaxis = "y2",
    y = pca__explainedVar_cumul * 100,
    x = paste0("PC", seq(pca__explainedVar)),
    type = 'scatter',
    mode = 'lines+markers',
    marker = list(size = 10)
  ) %>%
  layout(
    legend = list(
      xanchor = "center",
      x = 0.5,
      y = 100,
      orientation = 'h'
    ),
    plot_bgcolor = 'rgba(0,0,0,0)',
    paper_bgcolor = 'rgba(0,0,0,0)',
    xaxis = list(
      title = "Principal Component Axis",
      categoryorder = "array",
      categoryarray = paste0("PC", seq(pca__explainedVar))
    ),
    yaxis = list(title = "% Explained Variance"),
    yaxis2 = list(
      title = "% Explained Variance (cumulative)",
      tickfont = list(color = "orange"),
      showgrid = F,
      side = "right",
      automargin = T,
      standoff = 40L
    )
  )
pca__p3
