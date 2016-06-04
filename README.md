##bcall
Call variants using a site-specific binomial model.
bcall uses the [Rmath](https://github.com/gatoravi/Rmath)
libraries for the pbinom function
and uses [cereal](https://uscilab.github.io/cereal/)
for serializing C++ objects and storing them to disk.

###Install
```
mkdir build && cd build && cmake .. && make
```

###Usage

####Get all options
```
./bcall
```

####Estimate priors on fixed sites
```
./bcall prior-dump-fixed samples_rclist.tsv dump_output SeqCap_EZ_Exome_v3_primary.bed.gz
```
