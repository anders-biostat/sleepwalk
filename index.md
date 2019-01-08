<h2><center>Sleepwalk</center></h2>

# What is it?

When working with single-cell data, e.g., from single-cell RNA-Seq, we need an intuitive visual representation of our data. To do so, [dimensionality reduction](https://en.wikipedia.org/wiki/Dimensionality_reduction) approaches are useful, and dozens of these have been developed
by now: [PCA](https://en.wikipedia.org/wiki/Principal_component_analysis),
[Diffusion Map](https://en.wikipedia.org/wiki/Diffusion_map), [t-SNE](https://en.wikipedia.org/wiki/T-distributed_stochastic_neighbor_embedding),
[UMAP](https://umap-learn.readthedocs.io/en/latest/) and many others.

The aim of these dimension-reduced embeddings is simple: Each cell is represented by a point on two-dimensional plot, and the points are arranged such that similar cells appear close to each other and cells that are different farer apart. 

The same ideas are useful, if, insetad of many cells, we have many samples that have been assessed with some bulk omics (e.g., ordinary RNA-Seq), and we want to see which samples are similar.

*What do we mean by similar and different?* We think of each cell (or sample) as a point in "feature space" -- an imagined high-dimensional space whose axes are all the features (i.e., gene expressions, protein abundances, drug responses 
or any other measured values) and the cell's measured values for all the features are its coordinates. Since we usually work with tens, hundreds or thousands features, this is not very intuitive. Have you ever tried to imagine 1000-dimensional space?

However, we can imagine to measure **distances** in this high-dimensional space. The Euclidean distance between two cells, for example, is to simply take for each feature the difference between the two cells, add up the squares of all this distances and take the square root -- exactly as one does in oridnary 3D space using Pythagoras' theorem.

A dimension-reduced **embedding** is now an arrangement of point representing the cells in less dimensions -- only 2, if we want to show it on a two-deimensional computer screen -- that, to some extent, preserves these distances. 

This can never work perfectly -- we will always introduce **distortions** when reducing dimensions.

**So, can you be sure that the visualisation you get by using t-SNE, UMAP, MDS or the like really give you a faithful representation of your data?** Are the points that lie almost on top of each other really all similar? Does the large distance on your 2D representation
always mean lots of dissimilarities? **Our `sleepwalk` package for the *R* statistical programming environment can help you answer these questions.**

## Explore an embedding

Below is an example of a t-SNE visualisation that you can explore with `sleepwalk`.

This is single-cell transcriptomics data from the "CiteSeq" paper 
(Stoeckius et al., *Simultaneous epitope and transcriptome measurement in single cells*, [Nature Methods, 14:865](https://doi.org/10.1038/nmeth.4380), 2017).
Each point is a cell from a human cord blood sample. Normalised gene expression values were used as coordinates for the feature space, and we have used the [Rtsne package](https://cran.r-project.org/web/packages/Rtsne/index.html) to get a two-dimensional [t-SNE](https://en.wikipedia.org/wiki/T-distributed_stochastic_neighbor_embedding) embedding

<div class="aspect-ratio">
	<iframe src="examples/single_emb.html"></iframe>
</div>

**Move your mouse over the points.** At any moment the colours will show you the *real* (e.g. based on the original, feature-space coordinates) distance of all the cells to the one cell right under your mouse cursor. Thus, simply by moving the mouse you can explore the structure of the data and see whether the embedding gives a faithful representation and where it is distorted.

For instance, we immediately can see that the big upper cluster (let's call it *A*) is not as 
dense as the big lower one (this one we are going to call *B*).
When you move your mouse over any of the cells of cluster *B* all other cells in *B* immediately become dark green, showing that they are all really close, but that is not the case for cluster *A*: you can wander along its length and see the green colour "following" you and so see that *B* is extended over a gradient. 

You can also
observe the gradual change from *B* to the smaller cluster to its upper left, when moving through their thin connection. The small cluster at the very top of the plot is depcited really dense but you cann notice that, in reality, it is the least dense
one. Actually, to see any similarity between neighbours in this cluster, you will have to change the colour scale (click the `+` button twice). 

Try to 
look closely at all those smaller clusters between *A* and *B*: `sleepwalk` can help you to figure out how justified  their positioning on the plot is, which of them really are "between" *A* and *B*.
You can also notice that some of the cells at the edges of the clusters are drastically different from everything else around them. They likely are some sort of outliers.

## Compare two embeddings

`sleepwalk` can also help you compare two different embeddings. For example, let's look at the same t-SNE plot as before and put next to it a visualisation of the same data 
obtained with another dimensionality reduction approach. This time, we have used [UMAP](https://umap-learn.readthedocs.io/en/latest/).

<div class="aspect-ratio">
	<iframe src="examples/comp_emb.html"></iframe>
</div>

The main idea is the same: you move the mouse over an embedding and colour shows you the real distance from the current point to all others. But now you 
see distances from the same point on the other embedding as well. This allows you immediately identify the same clusters from both visualisations 
and gives you an easy way to compare two (or more if you want) dimensionality reduction techniques applied to your data.

## Compare two sets of samples

You can also use `sleepwalk` to compare not only different embeddings for the same data, but also to compare different samples. As example data we use here glioma samples from two patients (out of six, presented in the paper) 
by Filbin et al. (*Developmental and oncogenic programs in H3K27M gliomas dissected by single-cell RNA-seq*, [Science, 360:6386](https://doi.org/10.1126/science.aao4750), 2018).
Again, single-cell transcriptomics has been used to study each sample, and again we have used t-SNE to visualise them.

<div class="aspect-ratio">
	<iframe src="examples/comp_samp.html"></iframe>
</div>

Here, our goal is to figure out if there are correponding groups of cells in the two samples and what are those groups. `sleepwalk` can help here, too. The colour 
now shows the distance from the cell under the mouse cursor to all other cells in <i>all</i> the embeddings, allowing us to find the most cells not only in the current but also in the other samples. We can immediately see that the two largest
clusters not only correspond to each other in both samples, but also are aranged in a similar manner: Move the mouse along the large cluster in one of the plots to see it.
The small cluster at the bottom of the right-hand plot probably corresponds to a group of clusters on the other plot. We also can clearly see that each sample has a population of cells not present in the other one.

# Installation

To install `sleepwalk`, start R and type

```r
install.packages( "devtools" )
devtools::install_github( "anders-biostat/JsRCom" )
devtools::install_github( "anders-biostat/sleepwalk" )
```

# Usage

The `sleepwalk` package is easy to use. It has only two functions: `sleepwalk()` generates a single plot as in the first example, `sleepwalkMulti()` generates several plots as in the later examples.

Both functions take the following arguments

- First, pass the embedding, i.e., a *n* x 2 matrix of 2D coordinates for the *n* points. For `sleepwalkMulti`, pass a list with several such matrices, one for each of the plots.

- As second argument, pass the feature matrix: this is the *n* x *m* matrix with one row for each of the *n* points (e.g., cells) and one columns for each of the *m* features (e.g., genes) that should be used to calculate the actual distances. Sleepwalk determines the colour of each point from the Euclidean distance between the row of the point to be coloured and the row of the point under the mouse cursor. For `sleepwalkMulti`, pass a list of several such matrices, one for each plot. If *m* is large, it can be a good idea to only pass the first 30 or so principal components of the sample. Use the `prcomp_irlba` function from the [`irlba` package](https://bwlewis.github.io/irlba/), to efficiently calculate them.

- The third argument is the maximum value for the colour scale. You can adjust this afterwards with the "+" and "-" buttons next to the colour scale. For `sleepwalkMulti`, pass a vector, with one value for each plot.

- The optional parameter `pointsize` allwos you to change the size of the points.

- `same` is an optional parameter, only for `sleepwalkMulti`. It specifies whether the samples share the same objects (`same="objects"`), as in the second example with two embeddings compared, or the same features (`same="features"`), as in the third examples, where we compared two samples). The default is `same="objects"`. Note that for same `same="objects"`, all the matrices in the first and second argument must have the same number of rows. For `same="features"`, the feature matrices in the second argument must have the same number of columns.


# Authors

Sleepwalk is beeing developed by [Svetlana Ovchinnikova](mailto:s.ovchinnikova@zmbh.uni-heidelberg.de) and [Simon Anders](mailto:s.anders@zmbh.uni-heidelberg.de) at the [Center for Molecular Biology of the University of Heidelberg](https://www.zmbh.uni-heidelberg.de/anders). Please contact us for questions or feedback.


