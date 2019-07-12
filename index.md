<h2><center>Sleepwalk</center></h2>

# What is it?

When working with single-cell data, e.g., from single-cell RNA-Seq, we need an intuitive visual representation of our data. To do so, [dimensionality reduction](https://en.wikipedia.org/wiki/Dimensionality_reduction) approaches are useful, and dozens of these have been developed
by now: [PCA](https://en.wikipedia.org/wiki/Principal_component_analysis),
[Diffusion Map](https://en.wikipedia.org/wiki/Diffusion_map), [t-SNE](https://en.wikipedia.org/wiki/T-distributed_stochastic_neighbor_embedding),
[UMAP](https://umap-learn.readthedocs.io/en/latest/) and many others.

The aim of these dimension-reduced embeddings is simple: Each cell is represented by a point on two-dimensional plot, and the points are arranged such that similar cells appear close to each other and cells that are different farer apart. 

The same ideas are useful, if, instead of many cells, we have many samples that have been assessed with some bulk omics (e.g., ordinary RNA-Seq), and we want to see which samples are similar.

*What do we mean by similar and different?* We think of each cell (or sample) as a point in "feature space" -- an imagined high-dimensional space whose axes are all the features (i.e., gene expressions, protein abundances, drug responses 
or any other measured values) and the cell's measured values for all the features are its coordinates. Since we usually work with tens, hundreds or thousands features, this is not very intuitive. Have you ever tried to imagine 1000-dimensional space?

However, we can imagine to measure **distances** in this high-dimensional space. The Euclidean distance between two cells, for example, is to simply take for each feature the difference between the two cells, add up the squares of all this distances and take the square root -- exactly as one does in ordinary 3D space using Pythagoras' theorem.

A dimension-reduced **embedding** now is an arrangement of point representing the cells in less dimensions -- only 2, if we want to show it on a two-dimensional computer screen -- that, to some extent, preserves these distances. 

This can never work perfectly -- we will always introduce **distortions** when reducing dimensions.

**So, can you be sure that the visualisation you get by using t-SNE, UMAP, MDS or the like really give you a faithful representation of your data?** Are the points that lie almost on top of each other really all similar? Does the large distance on your 2D representation
always mean lots of dissimilarities? **Our `sleepwalk` package for the *R* statistical programming environment can help you answer these questions.**

## Explore an embedding

Below is an example of a t-SNE visualisation that you can explore with `sleepwalk`.

This is single-cell transcriptomics data from the "CiteSeq" paper 
(Stoeckius et al., *Simultaneous epitope and transcriptome measurement in single cells*, [Nature Methods, 14:865](https://doi.org/10.1038/nmeth.4380), 2017).
Each point is a cell from a human cord blood sample. A two-dimensional [t-SNE](https://en.wikipedia.org/wiki/T-distributed_stochastic_neighbor_embedding) embedding was obtained by running a [Seurat 
workflow](https://satijalab.org/seurat/multimodal_vignette.html). Gene expression values were normalised and scaled, [PCA](https://en.wikipedia.org/wiki/Principal_component_analysis) was performed
and all but the first 13 principal components were omitted to reduce noise levels. After that, tSNE was used to further reduce the number of dimensions and get a 2D visualisation.

<div class="aspect-ratio">
	<iframe src="examples/single_emb.html"></iframe>
</div>

To the left you can see cell types assigned to each cluster in the [Seurat vignette](https://satijalab.org/seurat/multimodal_vignette.html).

**Move your mouse over the gray points.** At any moment the colours will show you the "real" distances of all the cells to the one cell right under your mouse cursor. By "real distances", we mean those based on the original data, before the dimension reduction, and, by default, we calculate them as Euclidean distances. Thus, simply by moving the mouse you can explore the structure of the data and see where the embedding gives a faithful representation and where it is distorted.

For instance, we immediately can see that the cluster of T cells (red and green), which is the largest on the plot, is in fact in fact among the most dense ones. If you put the mouse over one of the T cells, almost the entire 
cluster immediately shows up in black or dark-green. Compare this, for example, with the seemingly small and compact cluster of B cells (grass green), or, even more so, with the megakaryocytes (lilac, close to the T cells) or erythrocytes (pink, upper-right corner). These small clusters contain cells that are further away from each other than some entire clusters are. To see any similarity between megakaryocytes you need to increase the range of the colour scale. (Press the “+” button).

We can also see the reason why tSNE failed to separate CD4+ (red) and CD8+ T cells (green). The differences between the two are so subtle that we can hardly see them. Moreover, there is no distinct border between them.
Try to move the mouse to the very right tip of the CD8+ T cells cluster to see a continuous change from the CD8+ T cells to CD4+ T cells (from right to the left). Yet we know that CD4+ T cells do not mature to become CD8+ nor vice verse. 

Notice how, when you explore the T cells cluster, the upper tip of the B cells cluster lights up. After moving the mouse there we may notice that these cells, even though they are clustered as B cells, resemble other cells in their cluster as much as T cells. This may be an indication of doublets and can be worth further investigation. If you had a running R session instead of exploring this web page, you would be able to select some cells with your mouse. For this you need to press the left mouse button and enclose the points in a selection contour. When you finish the selection, a variable with indices of all selected points will appear in your R session and you can go on with further analysis. 

You can also notice that some of the cells at the edges of the clusters are drastically different from everything else around them. They likely are some sort of outliers.

## Compare two embeddings

`sleepwalk` can also help you compare two different embeddings. For example, let's look at the same t-SNE plot as before and put next to it a visualisation of the same data 
obtained with another dimensionality reduction approach. This time, we have used [UMAP](https://umap-learn.readthedocs.io/en/latest/).

<div class="aspect-ratio">
	<iframe src="examples/comp_emb.html"></iframe>
</div>

The main idea is the same: you move the mouse over an embedding and colour shows you the real distance from the current point to all others. But now you see distances from the same point on the other embedding as well. This allows you immediately identify the same clusters from both visualisations 
and gives you an easy way to compare two (or more if you want) dimensionality reduction techniques applied to your data.

## Compare several sets of samples

You can use `sleepwalk` to compare not only different embeddings for the same data, but also to compare different samples. As example data we use here samples of the developing murine cerebellum at different time points by Carter et al., (*A Single-Cell Transcriptional Atlas of the Developing Murine Cerebellum*, [Current Biology, 28:18](https://doi.org/10.1016/j.cub.2018.07.062), 2018).
Again, single-cell transcriptomics has been used to study each sample, but here we have used UMAP to visualise them. Two of the samples are biological replicates both taken at time point E13.5. The third one has been taken at time point E14.5.

<div class="aspect-ratio">
	<iframe src="examples/comp_samp.html"></iframe>
</div>

Here, our goal is to figure out if there are corresponding groups of cells in the two samples and what those groups are. `sleepwalk` can help here, too. The colour now shows the distance from the cell under the mouse cursor to all other cells in <i>all</i> the embeddings, allowing us to find the most similar cells not only in the current but also in the other samples. Exploring the data with the mouse shows the two replicas of E13.5 samples are almost identical. The two branches (GABAergic and glutamatergic neurons) can be easily followed from the early progenitor cells up to the most differentiated ones. Comparing between the two E13.5 replicates reveals which aspects of the peculiar two-pronged shape of the glutamatergic branch is simply due to random variation and what seems reproducible. In the later E14.5 sample, the branches have disconnected from the progenitor cells, but Sleepwalk allows us to still identify corresponding cells. Interestingly, Sleepwalk can show that the GABAergic lineage is differentiated further in the E14.5 than in E13.5 samples, as the endpoint of the branch in E14.5 corresponds to an intermediate point in E13.5. Sleepwalk allows one to discover such details immediately, with minimal effort. Of course, such a visual exploration cannot replace a tailored detailed analysis but it is does provide a starting point and a first overview.

# Installation

Sleepwalk is available on [CRAN](https://cran.r-project.org/) and can be installed like any other CRAN package, with

```r
install.packages( "sleepwalk" )
```

# Usage

The `sleepwalk` package is easy to use. It has only one function, also called `sleepwalk()`.

As its first argument (`embeddings`), `sleepwalk` expects an embedding, i.e., a *n* x 2 matrix of 2D coordinates for the *n* points. To get multiple plot (as in the second and third example above), pass a list with several such matrices, one for each of the plots.

As second argument (`featureMatrices`), pass the feature matrix: this is the *n* x *m* matrix with one row for each of the *n* points (e.g., cells) and one columns for each of the *m* features (e.g., genes) that should be used to calculate the actual distances. Sleepwalk determines the colour of each point from the Euclidean distance between the row of the point to be coloured and the row of the point under the mouse cursor. For multiple plots, pass a list of several such matrices, one for each plot. 

If *m* is large, it can be a good idea to only pass the first 30 or so principal components of the sample. Use the `prcomp_irlba` function from the [`irlba` package](https://bwlewis.github.io/irlba/), to efficiently calculate them.

Instead of passing the *n* x *m* feature matrix, you can use the positional argument `distances` and pass a square *n* x *n* distance matrix, giving the distances between all pairs of points that should be used to determine the point colours. Use this if you think that simple Euclidean distances are not appropriate for your data and you have something more suitable.

You also have to specify `maxdists`: the maximum value for the colour scale. You can adjust this afterwards with the "+" and "-" buttons next to the colour scale. For mutiple plots, you can pass a vector, with one value for each plot.

If you ask for mroe than one plot, you can specify, with the parameter `same`, whether the plots share the same objects (`same="objects"`), as above in the second example with two embeddings compared, or the same features (`same="features"`), as in the third examples, where we compared two samples. The default is `same="objects"`. Note that for same `same="objects"`, all the matrices in the first and second argument must have the same number of rows. For `same="features"`, the feature matrices in the second argument must have the same number of columns.

The optional parameter `pointsize` allows you to change the size of the points.

Finally, instead of opening the app in a browser, you can ask to **save everything to an HTML file**, by passing a file name in the optional `saveToFile` parameter. This file will then contaiin all your data and the sleepwalk app code (in JavaScript), so that you can view it in any web browser without needing to have an R session running. This can be useful to share your sleepwalk visualization with colleagues unfamiliar with R or to provide your embedding as an interactive supplement in a publication.

# For Seurat users

If you have a Seurat data object, and have already run Seurat's `RunPCA` and `RunTSNE` functions on it, then you can use Sleepwalk simply by executing the following command.

for Seurat version 2.x:

```
sleepwalk( seu@dr$tsne@cell.embeddings, seu@dr$pca@cell.embeddings )
```

for Seurat version 3.x:

```
sleepwalk( seu$tsne@cell.embeddings, seu$pca@cell.embeddings )
```

This takes the t-SNE embedding stored in the Seurat data object `seu` and displays it with t-SNE. Replace `tsne` with `umap` in case you have used `RunUMAP` instead of `RunTSNE`.

(Note: The feature-space distances are calculated from the PCA coordinates, as this is what the t-SNE function gets by default. It is, of course, more accurate to pass as second parameter not the principal components buts the full data, and to do so, replace `seu@dr$pca@cell.embeddings` with `seu@data`. However, `data` is a huge matrix, and Sleepwalk might become quite slow.)

# Authors

Sleepwalk is being developed by [Svetlana Ovchinnikova](mailto:s.ovchinnikova@zmbh.uni-heidelberg.de) and [Simon Anders](mailto:s.anders@zmbh.uni-heidelberg.de) at the [Center for Molecular Biology of the University of Heidelberg](https://www.zmbh.uni-heidelberg.de/anders). Please contact us for questions or feedback, or file an [issue](https://github.com/anders-biostat/sleepwalk/issues) on GitHub if you find a bug.

# Preprint

A preprint describing Sleepwalk in more detail is available at BioRxiv:

S. Ovchinnikova and S. Anders: *Exploring dimension-reduced embeddings with Sleepwalk.* BioRvix 603589 (2019). [doi:10.1101/603589](https://doi.org/10.1101/603589).

Please cite this paper if you use Sleepwalk in your research work.

<!-- Default Statcounter code for anders-biostat.github.io
https://anders-biostat.github.io -->
<script type="text/javascript">
var sc_project=11967766; 
var sc_invisible=1; 
var sc_security="4f345b9a"; 
</script>
<script type="text/javascript"
src="https://www.statcounter.com/counter/counter.js"
async></script>
<noscript><div class="statcounter"><a title="Web Analytics"
href="https://statcounter.com/" target="_blank"><img
class="statcounter"
src="https://c.statcounter.com/11967766/0/4f345b9a/1/"
alt="Web Analytics"></a></div></noscript>
<!-- End of Statcounter Code -->