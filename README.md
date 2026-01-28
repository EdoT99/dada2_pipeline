## Dockerized-DADA2-pipeline 
### THis is an example of DADA2 pipeline for the analysis of 16S amplicon seqeunces using the
### DADA2 package from (Github)[https://github.com/benjjneb/dada2]

This is a test

### Creates an R working environment using a preexisitng image, and install required packages
```bash

    sudo docker build -f dockerfile -t dada2_pipeline

```

### Downloads Latest Silva taxonomy database and perform the pipeline

```bash

    sudo docker compose up

```