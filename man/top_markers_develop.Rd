% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/top_markers.R
\name{top_markers_develop}
\alias{top_markers_develop}
\title{Identify the genes most specifically expressed in groups of cells}
\usage{
top_markers_develop(
  cds,
  group_cells_by = "cluster",
  genes_to_test_per_group = 25,
  reduction_method = "UMAP",
  marker_sig_test = TRUE,
  reference_cells = NULL,
  speedglm.maxiter = 25,
  cores = 1,
  verbose = FALSE
)
}
\arguments{
\item{cds}{A cell_data_set object to calculate top markers for.}

\item{group_cells_by}{String indicating what to group cells by for
comparison. Default is "cluster".}

\item{genes_to_test_per_group}{Numeric, how many genes of the top ranked
specific genes by Jenson-Shannon to do the more expensive regression test
on.}

\item{reduction_method}{String indicating the method used for dimensionality
reduction. Currently only "UMAP" is supported.}

\item{marker_sig_test}{A flag indicating whether to assess the discriminative
power of each marker through logistic regression. Can be slow, consider
disabling to speed up top_markers().}

\item{reference_cells}{If provided, top_markers will perform the marker
significance test against a "reference set" of cells. Must be either a list
of cell ids from colnames(cds), or a positive integer. If the latter, top_markers()
will randomly select the specified number of reference cells. Accelerates
the marker significance test at some cost in sensitivity.}

\item{speedglm.maxiter}{Maximum number of iterations allowed for fitting GLM
models when testing markers for cell group.}

\item{cores}{Number of cores to use.}

\item{verbose}{Whether to print verbose progress output.}
}
\value{
a data.frame where the rows are genes and the columns are
\itemize{
\item gene_id vector of gene names
\item gene_short_name vector of gene short names
\item cell_group character vector of the cell group to which the cell belongs
\item marker_score numeric vector of marker scores as the fraction expressing scaled by the specificity. The value ranges from 0 to 1.
\item mean_expression numeric vector of mean normalized expression of the gene in the cell group
\item fraction_expressing numeric vector of fraction of cells expressing the gene within the cell group
\item specificity numeric vector of a measure of how specific the gene's expression is to the cell group based on the Jensen-Shannon divergence. The value ranges from 0 to 1.
\item pseudo_R2 numeric vector of pseudo R-squared values, a measure of how well the gene expression model fits the categorical data relative to the null model. The value ranges from 0 to 1.
\item marker_test_p_value numeric vector of likelihood ratio p-values
\item marker_test_q_value numeric vector of likelihood ratio q-values
}
}
\description{
Identify the genes most specifically expressed in groups of cells
}
\examples{
  \donttest{
    library(dplyr)

    cell_metadata <- readRDS(system.file('extdata',
                                         'worm_embryo/worm_embryo_SummarizedExperiment::colData.rds',
                                         package='monocle3'))
    gene_metadata <- readRDS(system.file('extdata',
                                         'worm_embryo/worm_embryo_SummarizedExperiment::rowData.rds',
                                         package='monocle3'))
    expression_matrix <- readRDS(system.file('extdata',
                                             'worm_embryo/worm_embryo_expression_matrix.rds',
                                             package='monocle3'))

    cds <- new_cell_data_set(expression_data=expression_matrix,
                             cell_metadata=cell_metadata,
                             gene_metadata=gene_metadata)

    cds <- preprocess_cds(cds)
    cds <- reduce_dimension(cds)
    cds <- cluster_cells(cds)
    marker_test_res <- top_markers(cds, group_cells_by="partition", reference_cells=1000)
    top_specific_markers <- marker_test_res \%>\%
                             filter(fraction_expressing >= 0.10) \%>\%
                             group_by(cell_group) \%>\%
                             top_n(1, pseudo_R2)
    top_specific_marker_ids <- unique(top_specific_markers \%>\% pull(gene_id))
  }

}
