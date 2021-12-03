[BLASThits.pdf](https://github.com/MaudeDavidLab/GMEmbeddings/files/7047752/BLASThits.pdf)
[SequenceTable.pdf](https://github.com/MaudeDavidLab/GMEmbeddings/files/7047698/SequenceTable.pdf)
# GMEmbeddings
# See embed_example.Rmd as well

## 1. Run BLAST to align your sequences against embeddings sequences.
### 1a: install blast
To install BLAST, follow instructions at the following link

https://www.ncbi.nlm.nih.gov/books/NBK279671/ 
[](https://www.ncbi.nlm.nih.gov/books/NBK279671/)

### 1b: Download blastdb embedding sequences: 
The files can be found at: 

http://files.cgrb.oregonstate.edu/David_Lab/microbiome_embeddings/blastdb_fullseq//
[](http://files.cgrb.oregonstate.edu/David_Lab/microbiome_embeddings/blastdb_fullseq//)

### 1c: Align your sequences to the database sequences. Here is how:

Run BLAST to rename your sequence with the nearest sequence available in the embedding matrix.
  Command should be similar to:
 
  ```
  blast_software_dir/blastn -db path_to_blastdb_dir/embedding_db_.07 -query path_to_fasta_file -out blast_hits.tsv -outfmt "6 qseqid sseqid qseq sseq evalue bitscore length pident"
  ```

Here is an example that I would use on my own machine:
```
ncbi-blast-2.11.0+/bin/blastn -db blastdb/embedding_db_.07 -query fasta_test.fasta -out blast_hits.tsv -outfmt "6 qseqid sseqid qseq sseq evalue bitscore length pident"
```

## 2. To increase speed of the next steps, filter blast hits outside of R. Output will be called best_hits.tsv and will be in the data_dir folder you provide.
```
GMEmbeddings/R/making_embedding_transformation_matrix/filter_blast_hits.sh data_dir/with/blast_hits/
```

## 3. Install package

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


## 4. Read in your sequence table. 
Read in your sequence table file. Different methods may be used depending on what type of file format you have.
After being read in the sequence table should look like this, with ids in the columns and sample ids in the rows. The ids of the columns must match the ids in the fasta file used above:

<img width="500" alt="Screen Shot 2021-08-25 at 10 32 06 AM" src="https://user-images.githubusercontent.com/68047298/130810454-6852a55a-5e1b-469f-b1d3-2ce774ce76ab.png">


An example sequence table can be obtained using the following command:
```
seqtab <- read.csv(system.file("extdata", "test_dataset_1/asv_table.csv", package = "GMEmbeddings"), row.names = 1)
seqtab <- t(seqtab)
```

## 5. Read in the hits from running blast. 
```
best_hits <- read.delim("path to best hits file", header = FALSE, sep = " ")
```
An example file can be read in
```
best_hits <- read.delim(system.file("extdata", "test_dataset_1/best_hits.tsv", package = "GMEmbeddings"), header = FALSE, sep = " ")
colnames(best_hits) <- c("qseqid", "sseqid", "qseq", "sseq", "evalue", "bitscore", "length", "pident")
```

We now need to add column names to our blast_hits file. To do this, use the following command:
```
colnames(best_hits) <- c("qseqid", "sseqid", "qseq", "sseq", "evalue", "bitscore", "length", "pident")
```

<img width="500" alt="Screen Shot 2021-08-25 at 10 39 19 AM" src="https://user-images.githubusercontent.com/68047298/130811280-88875daa-bb60-4b39-aa1b-837532558855.png">

## 6. Read in your chosen embedding transformation matrix. Options are available in glove_transformation_matrices and pca_transformation_matrices. Use one of the files with "id" in the filename. We recommend using 50 dimensions of either GloVe or PCA matrices.
```
embedding_filepath <- system.file("extdata/glove_transformation_matrices/", "glove_emb_id_50.txt", package = "GMEmbeddings")
embedding_matrix <- read.delim(embedding_filepath, row.names = 1, sep = "\t")
embedding_matrix <- embedding_matrix[rownames(embedding_matrix) != "<unk>", ]
```

## 7. Embed your sequence table
```
results = EmbedAsvTable(seqtab, best_hits, embedding_matrix)
embedded <- result$embedded
num_seqs_aligned <- result$num_seqs_aligned
percent_sequences_aligned <- result$percent_sequences_aligned
```
Please keep in mind that the column names of the seqtab MUST match the qseqid in the blast_hits file! If they do not match patterns, the `EmbedAsvTable` function will throw an error.

## Authors
Christine Tataru

Austin Eaton

## License
GPL-3.0-or-later

## Acknowledgments
