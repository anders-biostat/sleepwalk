# sleepwalk 0.3.1

* broken path to the start page, caused by `jrc` update fixed

# sleepwalk 0.3.0

* New argument `metric` allows to use angular distance (`metric = "cosine"`) as an alternative to default Euclidean distance
(`meric = "euclid"`).

* If `compare = "distances"`, it is no longer required to provide several embeddings. If only one embedding is given, it
will be used for all the distances.

# sleepwalk 0.2.1

* Changes due to an update of the `jrc` package.

* Indices of selected points are no longer stored in a variable and can be accessed only via the callback function. 
Thus, no changes to the global environment are made, unless user specifies them his- or herself.

* Added the possibility to pass arguments to `jrc::openPage` (such as port number or browser in which to open the app.)

# sleepwalk 0.2.0

* Now HTML Canvas is used to plot the embedding. It makes Sleepwalk faster and allows to simultaneously display more points.

* New parameter `mode = c("canvas", "svg")` is added, that allows user to go back to the old SVG-based version of Sleepwalk app.

* Bug in `slw_snapshot` is fixed. The function no longer returns a list of identical plots, when used with several different embeddings.