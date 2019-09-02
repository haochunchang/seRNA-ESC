# seRNA-ESC

* The raw data is stored using [git-lfs](https://git-lfs.github.com/)

```bash
git pull https://github.com/haochunchang/seRNA-ESC.git
```

------------------------------------------------------------------
## A starting docker image

* Installed most of the packages needed.

```bash
docker run -it -v /path/to/seRNA-ESC:/seRNA-ESC/ \
           -p 8888:8888 \ # for displaying jupyter notebooks
           --name container_name \
           jr55662003/serna-esc:v1 bash
```

* In this docker, use ```python3.5``` to call scripts instead of ```python3```.
* _Since some computation takes time, this image did not test all of the scripts._
* Tested:
    * Preprocessing
    * NMF
    * Co-expression analysis



----------------------------------------------------------------
## Requirements

The following are the packages we used for developing the scripts:
**Please make sure you have those packages before running**

```bash
apt-get install bedtools
```

### Python 3.5
* [numpy](https://docs.scipy.org/doc/numpy/) 1.12.1
* [pandas](https://pandas.pydata.org/pandas-docs/stable/) 0.19.2
* [scikit-learn](http://scikit-learn.org/stable/) 0.18.1
* [networkx](https://networkx.github.io/) 1.11 (For stitching enhancers and BiCoxNet)
* [fastcluster](http://danifold.net/fastcluster.html) 1.1.20 (For NMF & WGCNA)
* [matplotlib](https://matplotlib.org/)
* [seaborn](https://seaborn.pydata.org/)

For checking python packages and install them, simply go into main directory and run:
```bash
pip3 install -r requirements.txt
```

### R 3.3.2
* [WGCNA](https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/)
* [limma](https://bioconductor.org/packages/release/bioc/html/limma.html)
* [impute](http://bioconductor.org/packages/release/bioc/html/impute.html)
* MASS, class, cluster, data.table can all be installed from [CRAN](https://cran.r-project.org/)
* [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html) (Enable WGCNA uses all CPUs in your machine)

For checking R packages and install them, you can run the R script in main directory:
```bash
Rscript check_requirement.R
```

