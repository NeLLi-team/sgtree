# Simple Genome Tree (SGTree)

Modified June 2022 - Juan C. Villada (JGI)
Created (v0.1_) 05/19/2019

Simple Genome Tree (SGTree) is a computational pipeline for fast and easy construction of phylogenetic trees from a set of user provided genomes and a set of phylogenetic markers, in a taxonomic framework of de-replicated reference genomes. SGTree identifies conserved phylogenetic marker proteins and evaluates additional copies of markers derived from either duplications, horizontal gene transfer or contamination, to build a phylogenetic tree based on the concatenated alignment of selected marker proteins. 

## Setup

### Create conda environment and run sgtree

1. First clone the git repository.

```bash
git clone https://github.com/NeLLi-team/sgtree.git
```

2. Make sure you have anaconda installed, you can install it here. https://www.anaconda.com/distribution/#download-section

3. Depending on your OS, choose the osx_env.txt for mac or the linux_env.txt file for Ubuntu or linux. 

4. Next run (where `<spec-file>` is either `linux_env.txt` or `osx_env.txt`): 

```bash
conda create --name sgtree --file <spec-file>
conda activate sgtree
```  


5. Make `sgtree` executable:

```bash
cd sgtree/
chmod u+x sgtree.py
```


## Run SGTree

5. Run sgtree with the provided set of query genomes and models for testing, user can control the number of CPUs used by the computer, the minimum percentage of models with hits for a genome to be considered as part of the dataset as well as a directory with reference genomes. Genomes from the query genomes directory and the reference genomes directory will be colored red and grey respectively. 

```bash
# test example
./sgtree.py testgenomes/Chloroflexi hmms/UNI56 --num_cpus 8

# general example
./sgtree.py <genomes directory> <models directory> --num_cpus <integer> --percent_models <integer> --ref <reference genomes directory>
```
	
For running sgtree_final.py, please make a ref_concat folder wherever you want to store reference db runs. 


```bash
python3 sgtree_final.py testgenomes/Chloroflexi hmms/UNI56 sgtree/references_concat --num_cpus 10 --save_dir sgtree/test --ref sgtree/testgenomes/chlorref

python3 sgtree_final.py <genomes directory> <models directory> <references directory>--num_cpus <integer> --percent_models <integer> --ref <reference genomes directory>
```


## Important note
1. Genomes must have the header format as follows: 

``` bash
>IMG2684622718|2685462912
MLCAFAEEEAKIAETVGKVATELKVKKLLSDFATKEGEEHISTYNKIAMTAKAEGYADIEAMLCAFAEEEAKLQKL
```
where the first field before the pipe contains the genome identifier which matches the filename base, and the second field a unique protein identifier.
