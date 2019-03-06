## iMKT data repository
Pipelines to extract unfolded site frequency spectrum from 1000GP VCF and DGN alignments.
Each pipeline returns a tab-delimited file including information unfolded site frequency spectrum, analyzable sites by largest transcrip and divergence by genes, populations and MKT functional classes (0-fold: selected class; 4-fold: neutral class). Both allows the *Drosophila melanogaster* and humans proteins analysis through [iMKT web-service](https://imkt.uab.cat) and [iMKT R-package](https://github.com/BGD-UAB/iMKT).  
This repository only include raw code to get main results. notebooks/ folder include two main Jupyter Notebooks running on Python 3.6 kernel to execute step by step the pipeline. src/ folder contain raw scripts to needed to execute the pipelin. Please note that multiple step could be parallelized, in this case create yourself customs bash scripts or run it on your server manually.   
Pipeline were developed in the conda enviroment *imktData.yml* in local server: 100GB RAM and 16 Intel(R) Xeon(R) CPU.  
In addition *structure.sh* deposited in scr/ create the folders we used to complete the whole process. If you decided execute it, ovewrite notebook/ and src/ with the same folders deposited at this repository.

## Data retrieve
To execute the pipelines requiere download the following files. Paths would need to be changed too. 
### *D. melanogaster* population genomic data.
Variation data generated by the Drosophila Genome Nexus, together with divergence data between D. melanogaster and D. simulans, was retrieved from PopFly ([Hervás et al. 2017](https://academic.oup.com/bioinformatics/article/33/17/2779/3796397)) in FASTA format (also available in [DGN web site](http://www.johnpool.net/genomes.html)). Recomb data from [Comeron et al. 2012](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002905)  
### Human population genomic data. 
Genome variation data and information of the ancestral state of the variants generated by the 1000GP Phase III (1000 Genomes Project Consortium 2015), together with divergence between humans and chimpanzees, were retrieved from PopHuman ([Casillas et al. 2018](https://academic.oup.com/nar/article/46/D1/D1003/4559406)) in Variant Call Format (VCF). Recomb data from [Bhèrer et al. 2017](https://www.nature.com/articles/ncomms14994). Pilot mask to exclude low quality variants download from [1000GP ftp](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks)

## Jupyter notebooks
### [01_*D. melanogaster* pipeline (python/bash)](https://github.com/jmurga/iMKTData/blob/master/notebooks/dmelProteins.ipynb)
### [02_Human pipeline (python/bash)](https://github.com/jmurga/iMKTData/blob/master/notebooks/humanProteins.ipynb)
