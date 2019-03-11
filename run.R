#!/usr/local/bin/Rscript

task <- dyncli::main()

library(dplyr, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)
library(dynwrap, warn.conflicts = FALSE)
library(dyncli, warn.conflicts = FALSE)
library(CellTrails, warn.conflicts = FALSE)

#####################################
###           LOAD DATA           ###
#####################################

# load data
expression <- task$expression %>% as.matrix()
params <- task$params
priors <- task$priors

# TIMING: done with preproc
timings <- list(method_afterpreproc = Sys.time())

#####################################
###        INFER TRAJECTORY       ###
#####################################
# steps from the vignette https://dcellwanger.github.io/CellTrails/

sce <- SingleCellExperiment(assays = list(logcounts = t(expression)))

# filter features
if (isTRUE(params$filter_features)) {
  trajFeatureNames(sce) <- filterTrajFeaturesByDL(sce, threshold = params$threshold_dl, show_plot = FALSE)
  trajFeatureNames(sce) <- filterTrajFeaturesByCOV(sce, threshold = params$threshold_cov, show_plot = FALSE)
  trajFeatureNames(sce) <- filterTrajFeaturesByFF(sce, threshold = params$threshold_ff, min_expr = params$min_expr, show_plot = FALSE)
}

# filter cells based on the features
sce <- sce[,apply(logcounts(sce[trajFeatureNames(sce), ]), 2, sd) > 0]

# dimensionality reduction
se <- CellTrails::embedSamples(sce)
d <- CellTrails::findSpectrum(se$eigenvalues, frac = params$frac)
CellTrails::latentSpace(sce) <- se$components[, d]

# find states
CellTrails::states(sce) <- sce %>% CellTrails::findStates(
  min_size = params$min_size,
  min_feat = params$min_feat,
  max_pval = params$max_pval,
  min_fc = params$min_fc
)

# construct tree
sce <- CellTrails::connectStates(sce, l = params$l)

# fit trajectory
# this object can contain multiple trajectories (= "components"), so we have to extract information for every one of them and combine afterwards
components <- CellTrails::trajComponents(sce)

trajectories <- map(
  seq_along(components),
  function(ix) {
    if (length(components[[ix]]) > 1) {
      traj <- CellTrails::selectTrajectory(sce, ix)
      CellTrails::fitTrajectory(traj)
    } else {
      components[[ix]]
    }
  }
)


timings$method_aftermethod <- Sys.time()

#   ____________________________________________________________________________
#   Process cell graph                                                      ####

cell_ids <- CellTrails::sampleNames(sce)
grouping <- CellTrails::states(sce) %>% as.character() %>% set_names(cell_ids)
dimred <- SingleCellExperiment::reducedDim(sce, type = "CellTrails")

cell_graph <- map_dfr(
  trajectories,
  function(traj) {
    if (is.character(traj)) {
      cell_ids <- colnames(sce)[which(states(sce) == traj)]
      data_frame(
        from = cell_ids[-length(cell_ids)],
        to = cell_ids[-1],
        length = 0,
        directed = FALSE
      )
    } else {
      graph <- CellTrails:::.trajGraph(traj)
      cell_ids_graph <- igraph::vertex.attributes(graph)$sampleName
      cell_graph <- graph %>%
        igraph::as_data_frame() %>%
        mutate(
          from = cell_ids_graph[as.numeric(from)],
          to = cell_ids_graph[as.numeric(to)],
          directed = FALSE
        ) %>%
        dplyr::rename(
          length = weight
        )
    }
  }
)

to_keep <- unique(c(cell_graph$from, cell_graph$to))

# save output
output <-
  wrap_data(
    cell_ids = cell_ids
  ) %>%
  add_dimred(
    dimred = dimred
  ) %>%
  add_grouping(
    grouping = grouping
  ) %>%
  add_cell_graph(
    cell_graph = cell_graph,
    to_keep = to_keep
  )  %>%
  add_timings(
    timings = timings
  )

dyncli::write_output(output, task$output)
