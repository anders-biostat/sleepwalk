# What is it for?

People, who work with lots of samples, generally want to have some nice visual represantation of their data.
One intuitive and commonly used approach to that is to imagine each sample as a dot and its features (gene expressions, 
protein abundances, drug responses or any other measured values) as its coordinates. Yet, since we usually work with tens,
hundreds or thousands features, this alone is not enough. Have you ever tried to imagine 1000-dimensional space?

What we want is to get an **embedding** - a relatively low-dimensional space, into which we can project our points/samples so that
the distances between them are preserved. We want similar samples to be represented by the points that are close to each other,
while big differences between them should correspond to large distances. For visualisation we usually want a 2D embedding.

Dozens of [dimensionality reduction](https://en.wikipedia.org/wiki/Dimensionality_reduction) approaches have been developed
by now to help you find an optimal embedding for your data: [PCA](https://en.wikipedia.org/wiki/Principal_component_analysis),
[Diffusion Map](https://en.wikipedia.org/wiki/Diffusion_map), [t-SNE](https://en.wikipedia.org/wiki/T-distributed_stochastic_neighbor_embedding),
[UMAP](https://umap-learn.readthedocs.io/en/latest/) and many others. But can you be sure that the visualisation you get out is really 
an optimal one? Are the points that lie almost on top of each other really so much similar? Does the large distance on your 2D representation
always mean lots of dissimilarities? Our `sleepwalk` package for **R** statistical programming environment can help you answer these questions.

## Explore an embedding

Bellow, there is an example of visualisation, created by the `sleepwalk` package.

Here, we are looking at single-cell transcriptomics data from the "CiteSeq" paper 
(Stoeckius et al., *Simultaneous epitope and transcriptome measurement in single cells*, [Nature Methods, 14:865](https://doi.org/10.1038/nmeth.4380), 2017).
Each point is a cell from a human cord blood sample, normalised gene expression values were used as coordinates, and the 2D embedding was created using 
[t-SNE](https://en.wikipedia.org/wiki/T-distributed_stochastic_neighbor_embedding).

Try to move your mouse over the points. At any moment the colour will show you the *real* (e.g. based on the original coordinates) distance from
the cells, at which your mouse is pointing right now, to all other cells. Therefore, simply by moving the mouse one can explore the structure of 
the data and decide wether the current representation is good enough.

<div class="aspect-ratio">
	<iframe src="examples/single_emb.html"></iframe>
</div>

For instance, we immediately can see that the big upper cluster (lets call it **A**) is not as dense as the big lower one (this one we are going to 
call **B**)

## Compare two embeddings

## Compare two sets of samples

# Installation

# Tutorial

# Usage

