# GMEmbeddings


## 1. Install package

### Prerequisites

Before installing this package, make sure that you have the required prerequisite packages already installed. These packages include `plyr` and `seqinr`.
If not already installed, use:
```
install.packages(c("plyr", "seqinr"))
```

Also, before installing `GMEmbeddings`, load the `devtools` package:
```
library(devtools)
```

### Installing

Download and install the `GMEmbeddings` package from GitHub:

```
install_git("https://github.com/MaudeDavidLab/GMEmbeddings")
```

Now load the package:

```
library(GMEmbeddings)
```


## 2. Read in the embedding matrix you'd like to use. 
Standard in .07_100

```
embedding_matrix <- getEmbeddingMatrix("path to embedding matrix file")
```

This package has an example embedding matrix for you to test if you like. It can be obtained using the following command:
```
embedding_matrix <- getEmbeddingMatrix(system.file("extdata", "embed_.07_100dim.txt", package = "GMEmbeddings"))
```

## 3. Read in your sequence table. It should look like this:
```
seqtab <- getExampleSeqtab("path to sequence table file")
```

An example sequence table can be obtained using the following command:
```
seqtab <- getSeqtab(system.file("extdata", "seqtab_test.csv", package = "GMEmbeddings"))
```

## 4. Read in your fasta file. ID's should be ASV ids from the column names of your sequence table, and sequences should be the full length ASV sequence. It should look like this:
```
fasta_file <- getExampleFasta("path to fasta file")
```


An example fasta file can be obtained using the following command:
```
fasta_file <- getFasta(system.file("extdata", "fasta_test.fasta", package = "GMEmbeddings"))
```

## 5.Run blast to align your sequences to the sequences in our embedding database. Here's how:
### 5a: install blast
To install BLAST, follow instructions at the following link

https://www.ncbi.nlm.nih.gov/books/NBK279671/ 
[](https://www.ncbi.nlm.nih.gov/books/NBK279671/)

### 5b: Download blastdb embedding sequences from CGRB: 
The files can be found at: 

http://files.cgrb.oregonstate.edu/David_Lab/microbiome_embeddings/data/blastdb/
[](http://files.cgrb.oregonstate.edu/David_Lab/microbiome_embeddings/data/blastdb/)

### 5c: Align your sequences to those sequences from 5b. Here is how:

Run BLAST to rename your sequence with the nearest sequence available in the embedding matrix.
  Command should be similar to:
 
  ```
  blast_software_dir/blastn -db path_to_blastdb_dir/embedding_db_.07 -query path_to_fasta_file -out output_file_name -outfmt "6 qseqid sseqid qseq sseq evalue bitscore length pident"
  ```

Here is an example that I would use on my own machine:
```
ncbi-blast-2.11.0+/bin/blastn -db blastdb/embedding_db_.07 -query fasta_test.fasta -out this_is_me_running_BLAST.tsv -outfmt "6 qseqid sseqid qseq sseq evalue bitscore length pident"
```


## 6. Read in the hits from running blast. 
```
blast_hits <- read.delim("path to blast hits file", header = FALSE)
```

We now need to add column names to our blast_hits file. To do this, use the following command:
```
colnames(blast_hits) <- c("qseqid", "sseqid", "qseq", "sseq", "evalue", "bitscore", "length", "pident")
```

If you would like to use the blast hits file that is preinstalled with the package, then use the following command. Please note that this file already has the column names added and no other changes need to be done.
```
blast_hits <- read.delim(system.file("extdata", "blast_hits.tsv", package = "GMEmbeddings"))
```

## 7. Embed your sequence table
```
EmbedAsvTable(seqtab, fasta_file, blast_hits, embedding_matrix)
```
Please keep in mind that the arguments passed in to the `EmbedAsvTable` function must be named: **seqtab**, **fasta_file**, **blast_hits**, and **embedding_file**! If the names of these objects are saved as anything else the `EmbedAsvTable` function will not execute properly.


## Authors

## License
GPL-3.0-or-later

## Acknowledgments
