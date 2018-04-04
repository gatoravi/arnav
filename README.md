## bcall
Call variants using a site-specific binomial model.

bcall uses the [Rmath](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#The-standalone-Rmath-library)
libraries for the pbinom function
and uses [cereal](https://uscilab.github.io/cereal/)
for serializing C++ objects and storing them to disk.

### Install
```
git clone --recursive https://github.com/gatoravi/bcall
mkdir build && cd build && cmake .. && make
```

### Usage

#### List all options available with bcall
```
./bcall
```

#### Estimate site-specific priors for a list of sites specified in a bed file
```
./bcall prior-dump-fixed samples_readcountfile_list.tsv merged_priors.dump list_of_regions.bed.gz
```

#### To call mutations using the merged priors
```
./bcall call-using-merged sample_readcountfile_list.tsv merged_priors.dump
```

For details on the file-formats used in the steps above please see the `Workflow` section below.

### Workflow
A workflow describing how bcall can be used to detect putative somatic mutations in RNA sequencing
data is described [here.](http://gatoravi.github.io/genetics/2017/10/12/mutations-rnaseq.html)
