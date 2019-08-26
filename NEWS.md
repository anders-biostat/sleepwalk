# sleepwalk 0.2.0

* Now HTML Canvas is used to plot the embedding. It makes Sleepwalk faster and allows to simultaneously display more points.

* New parameter `mode = c("canvas", "svg")` is added, that allows user to go back to the old SVG-based version of Sleepwalk app.

* Bug in `slw_snapshot` is fixed. The function no longer returns a list of identical plots, when used with several different embeddings.