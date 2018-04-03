##bcall

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

#### Get all options
```
./bcall
```

#### Estimate priors on fixed sites
```
./bcall prior-dump-fixed samples_rclist.tsv dump_output SeqCap_EZ_Exome_v3_primary.bed.gz
```

#### To call mutations using the merged priors
```
./bcall call-using-merged sample_readcountfile_list.tsv merged_priors.dump
```
