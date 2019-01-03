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

For instance, we immediately can see that the big upper cluster (lets call it **A**) is not as dense as the big lower one (this one we are going to call **B**).
When you move your mouse over any of the cells of cluster **B** all others immediately light up, but that's not the case for cluster **A**. You can also
observe the gradual change from **B** to the smaller cluster to its upper left. A very dense cluster at the very top of the plot is in reality the least dense
one. Actually, to see any similarity between neighbours in this cluster, you will have to change the colour scale (click the `+` sign two times). Try to 
look closely at all those smaller clusters between **A** and **B**. `sleepwalk` can help you to figure out, how justified is their positioning on the plot.
You can also notice that at leas some cells at the edges of the clusters are drastically different from everything else around them. It's a good chance that they
are some sort of outliers.

## Compare two embeddings

`sleepwalk` can also help you to compare two different embeddings. For example, let's look at the same t-SNE plot from above and put next to it a visualisation
obtained with another dimensionality reduction approach. In this case, we've used [UMAP](https://umap-learn.readthedocs.io/en/latest/).

<div class="aspect-ratio">
	<iframe src="examples/comp_emb.html"></iframe>
</div>

Main idea is the same: you move the mouse over an embedding and colour shows you the real distance from the current point to all others. But now you 
see distances from the same point on the other embedding as well. This allows you immediately identify the same clusters from both visualisations 
and gives you an easy way to compare two (or more if you want) dimensionality reduction techniques applied to your data.

## Compare two sets of samples

You can also use `sleepwalk` to compare not only dimensionality reduction approaches, but also two different samples. Now, as example data we are going
to use glioma samples from two patients (out of six, presented in the paper) 
from (Filbin et al., *Developmental and oncogenic programs in H3K27M gliomas dissected by single-cell RNA-seq*, [Science, 360:6386](https://doi.org/10.1126/science.aao4750), 2018).
Once again, single-cell transcriptomics was used to study each sample, and once again we will use a t-SNE plot to visualise them.

<div class="aspect-ratio">
	<iframe src="examples/comp_samp.html"></iframe>
</div>

Now, our goal is to figure out if there are correponding groups of cells in the two samples and what are those groups. And `sleepwalk` can help us here as well. The colour 
now shows the distance from one cell to all other cells in all the embeddings, allowing us to find the most similar clusters. We can immediately see that the two biggest
clusters not only correspond to each other in both samples, but also are aranged in a similar manner: Move the mouse along the biggest cluster in one of the plots to see it.
Small cluster in the bottom of the plot to the right, probably, corresponds to a group of clusters on the other plot. We also can clearly see that each sample has a population
of cells not present in the other one.

# Installation

To install `sleepwalk`, start R and type

```r
install.packages( "devtools" )
devtools::install_github( "anders-biostat/JsRCom" )
devtools::install_github( "anders-biostat/sleepwalk" )
```

# Usage

`sleepwalk` package has only two functions. `sleepwalk()` generates a single plot as in the first example, `sleepwalkMulti()`
generates several plots.

Both functions take a list of embeddings as their first argument (each embedding must be an `n_i x 2` matrix, where `n_i` is
a number of points in the sample `i`). Second argument is the list of features to calculate the real distances (each feature 
matrix must have points as rows and features as columns). Third required argument is maximal distance for colour scale.
An optional parameter `pointSize` allows one to change the size of the points. `same` is an optional parameter for 
`sleepwalkMulti()` and defines whether the samples share the same `objects` (as in the second example with two embeddings compared)
or the same `features` (as in the third examples, where we compared two samples).


