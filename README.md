# Simple Genome Tree (SGTree)

Simple Genome Tree (SGTree) is a computational pipeline for fast and easy construction of phylogenetic trees from a set of user provided genomes and a set of phylogenetic markers, in a taxonomic framework of de-replicated reference genomes. SGTree identifies conserved phylogenetic marker proteins and evaluates additional copies of markers derived from either duplications, horizontal gene transfer or contamination, to build a phylogenetic tree based on the concatenated alignment of selected marker proteins. 

## ⚙️ Setup

### Create conda environment and run sgtree

1. Clone the git repository.

```bash
git clone https://github.com/NeLLi-team/sgtree.git
```

2. Make sure you have [anaconda](https://www.anaconda.com/distribution/#download-section) installed.

3. Depending on your OS, choose the osx_env.txt for macOS or the linux_env.txt file for Linux. 

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


## 🚀 Run SGTree

6. Run sgtree with the provided set of query genomes and models for testing, user can control the number of CPUs used by the computer, the minimum percentage of models with hits for a genome to be considered as part of the dataset as well as a directory with reference genomes. Genomes from the query genomes directory and the reference genomes directory will be colored red and grey respectively. 

```bash
# test example
./sgtree.py testgenomes/Chloroflexi hmms/UNI56 --num_cpus 8

# general example
./sgtree.py \
	<genomes_dir> \
	<models_dir> \
	--num_cpus 8 \
	--percent_models 50 \
	--marker_selection yes \
	--aln mafft \
	--ref <reference_genomes_dir> \
	--save_dir <output_dir>
```

## 🚨 Important note 
Genomes must have the header format as follows: 

```bash
>IMG2684622718|2685462912
MLCAFAEEEAKIAETVGKVATELKVKKLLSDFATKEGEEHISTYNKIAMTAKAEGYADIEAMLCAFAEEEAKLQKL
```
where the first field before the pipe contains the genome identifier which matches the filename base, and the second field a unique protein identifier.


## Authors

|Authors		            | Email			            | Date		    |
|-----------------------|-----------------------|-------------|
|Ewan Whittaker-Walker 	| ewanww@berkeley.edu	  | 05/19/2019	|
|Frederik Schulz		    | fschulz@lbl.gov 	    | Since 2019	|
|Juan C. Villada 	      | jvillada@lbl.gov 	    | Since 2021	|
|Marianne Buscaglia	    | mbuscaglia@lbl.gov 	  | Since 2022	|
