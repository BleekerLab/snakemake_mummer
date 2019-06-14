# Snakemake pipeline for genome to genome alignment using MUMMER
MUMMER4
One Paragraph of project description goes here


## Table of contents
* [Getting started](#getting-started)
* [Prerequisites](#prerequisites)
* [Installing](#installing)
* [Tests](#running-the-tests)
* [Deployment](#deployment)
* [Status](#status)
* [Inspiration](#inspiration)
* [Contact](#contact)

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them

All packages were installed using conda
Get miniconda from https://docs.conda.io/en/latest/miniconda.html which contains Python 3.7 and the conda package manager.

The script also requires R


### Installing

A step by step series of examples that tell you how to get a development env running

The alignments and some processing are performed using functions from the mummer package.
```
conda install -c bioconda mummer
```

The R scripts require the following packages:
Biostrings
```
conda install -c bioconda bioconductor-biostrings
```
argparse
```
conda install -c bioconda r-argparse
```
tidyr
```
conda install -c r r-tidyr
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

To run some simple tests we have provided a set of publicly available sequences:  
10 reference sequences, Solanum_lycopersicum.SL3.0.dna.chromosome.8 was cut-up in 10 pieces H8.1-H8.10.  
4 query sequences, 4 BAC sequences from Matsuba et al. (2013) which belong to a region of chromosome 8.  
These sequences are located in test/refs/ and test/queries/ respectively.

The config.yaml file contains 4 important variables:  
queries_dir should be the directory in which your query sequences are stored.  
refs_dir should be the directory in which your reference sequences are stored.  
length_threshold is the minimum length you want your alignments to be.
identity_threshold is the minimum identity score you want your alignments to have.

The BAC sequences for this test are relatively small compared to the scaffold or super-scaffold to chromosome alignments that the pipeline is made for. Therefore we need to use a low length threshold, 100.
To ensure some matches of our sequence against the less likely to be matched references we set the identity threshold to 85.

Once the config.yaml file is set up we perform a dryrun to ensure we get the correct number of jobs.
```
snakemake -np --use-conda
```
This should yield a total of 162 jobs.  
40 jobs for each of nucmer, delta_to_coords, sort_coords, and calculate_alignment_percentage.  
1 job for create_results_matrix.  
And finally 1 job for all.  

If all seems right we run the real deal.
```
snakemake --use-conda
```
This should yield 41 files in total.  
40 query_vs_ref.sorted.coords  
And 1 results.tsv

The results.tsv file contains the % of aligned bases for each query-reference combination in a convenient format.

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags).

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
