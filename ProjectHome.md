# Authors #
  * Biao Li, Baylor College of Medicine
  * Gao Wang, Baylor College of Medicine
  * Suzanne M. Leal, Baylor College of Medicine

---


# Acknowledgment #
Development of PhenoMan was supported by the [NHLBI Exome Sequencing Project](https://esp.gs.washington.edu/drupal/), the [Minority Health GRID project](https://rcenterportal.msm.edu/node/146), and the [Centers for Mendelian Genomics](http://www.mendelian.org).

---


# Introduction #
Recently, the next generation sequencing and other high-throughput technology advances to promote great interest in detecting associations between complex phenotypic traits and genetic variants. Phenotype quality control procedure is crucial and could be a complicated issue that can largely impact association analysis results. Although various decisions are likely to be made on different traits, there is lacking a simple, effective and uniformed way to perform and record phenotype data cleaning steps. Arbitrariness and ambiguity often arise in managing phenotype data quality. This problem is normally neglected, which could cause biased meta-analysis or comparison of association results between studies on same traits. PhenoMan is an interactive program that integrates data exploration, management and quality control using a unified platform. It is featured by dissecting ambiguous data cleaning steps into a series of simple commands with dynamic arguments, which are clearly recording the entire data cleaning procedure. PhenoMan provides approaches in efficient exploration and management of phenotype data. We should perform them before the association analysis so as to ensure effective estimates from association results and consistent comparisons between different projects. PhenoMan can be used on new and existing rich sets of phenotypic data for association analysis of both quantitative and qualitative traits.

Please go to [\*Downloads\*](http://code.google.com/p/phenoman/downloads/list) area to download PhenoMan [source code](http://phenoman.googlecode.com/files/phenoman-src-0.1.0.tar.gz), [tutorial/documentation](http://phenoman.googlecode.com/files/doc_phenoman.pdf) and [example data](http://phenoman.googlecode.com/files/example_data.tar.gz).


---

# Installation #
  * PhenoMan requires installing Python2.7+ (Python 3.2+ preferred), R and ggplot2 0.9+ (a R package) as dependencies.
  1. **install [python](http://www.python.org/download/)** (need version 2.7 or above)
  1. **check if command 'python' (or 'python3') is in system path**
> > Open a terminal (Linux or Mac) or run a command prompt (Windows) and type in 'python' (or python3). If an error is
> > returned, it need to add python to system path.
    * For Linux or Mac, find where python has been installed and add the directory to ~/.bashrc (or ~/.bash\_profile).
    * For windows, refer to [this](http://stackoverflow.com/questions/3701646/how-to-add-to-the-pythonpath-in-windows-7)
  1. **install [R](http://www.r-project.org/)**
    * For Linux (Ubuntu-like system): run command "sudo apt-get install r-base r-base-dev" or compile from source code.
    * For Mac and Windows: download precompiled packages and add R to system path after installation.
  1. **install ggplot2**
    * Launch R and run command 'install.packages("ggplot2").
  1. **install PhenoMan**
    * Download [source code](http://phenoman.googlecode.com/files/phenoman-src-xxx.tar.gz)
    * Unarchive and use Python to install: "tar -xf phenoman-src-xxx.tar.gz, sudo python3 setup.py install"
    * Or to install locally:
      * python3 setup.py install --prefix=/path/to/your/local/folder
      * Set environmental variables by opening ~/.bash\_profile and add the following:
```
    PATH=/path/to/your/local/folder/bin:${PATH}
    export PATH
    PYTHONPATH=/path/to/your/local/folder/lib/pythonx.x/site-packages:$PYTHONPATH 
    export PYTHONPATH
```


---

# Tutorial #
  * Please download PhenoMan [documentation](http://phenoman.googlecode.com/files/doc_phenoman.pdf) and refer to Chapter 2, which includes tutorial on how to use PhenoMan and a variety of examples. Please also download and unzip [example data](http://phenoman.googlecode.com/files/example_data.tar.gz), for more details see Chapter 2.1 in the documentation.