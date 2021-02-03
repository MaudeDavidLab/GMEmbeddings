# GMEmbeddings

## 0: Readme: will download a bunch of stuff

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

## 3. Read in your sequence table. It should look like this:
```
seqtab <- getExampleSeqtab("path to sequence table file")
```

## 4. Read in your fasta file. ID's should be ASV ids from the column names of your sequence table, and sequences should be the full length ASV sequence. It should look like this:
```
fasta_file <- getExampleFasta("path to fasta file")
```

## 5.Run blast to align your sequences to the sequences in our embedding database. Here's how:
### 5a: install blast
To install BLAST, follow instructions at the following link

https://www.ncbi.nlm.nih.gov/books/NBK279671/ - automatic!
[NCBI](https://www.ncbi.nlm.nih.gov/books/NBK279671/)

### 5b: download blastdb embdding sequences from cgrb: http://files.cgrb.oregonstate.edu/David_Lab/microbiome_embeddings/
### 5c: Align your sequences to those sequences ^. Here's how: *line from other comments*. It should look like this:
```
#Instructions:
#1. Download all files in the folder http://files.cgrb.oregonstate.edu/David_Lab/microbiome_embeddings/data/blastdb/. Set variable below equal to path name of the download on your machine
blast_db = ""


#2. Run blast to rename your sequences with the nearest sequence available in the embedding matrix
blast_software_dir/blastn -db path_to_blast_db -query fasta.fasta -out output_file_name  -outfmt "6 qseqid sseqid qseq sseq evalue bitscore length pident"

best_hits = getBestHits(blast_hits = output_file_name, id_thresh = 99)


#3. Download the embedding transformation matrix from here: http://files.cgrb.oregonstate.edu/David_Lab/microbiome_embeddings/data/embed/embed_.07_100dim.txt
#set the variable below to the file name on your local machine
embedding_file_name = "C:/Users/ctata/Documents/Lab/GME/data/embed_matrices/embed_.07_100dim.txt"


#4. Embed your data
embedded = embedSeqtab(seqtab, best_hits, embedding_file_name)


#Usage Example:
#1. seqtab <- embedSeqtab(getExampleSeqtab(), fasta_file = getExampleFasta())

#2. seqtab <- getExampleSeqtab_Halfvarson()
#   seqtab_embed <- embedSeqtab(seqtab, "data/halfvarson/sequences.fasta")

#3.

#blast hits should have column names "qseqid"   "sseqid"   "qseq"     "sseq"     "evalue"   "bitscore"  "length"   "pident"
#blast_hits <- read.delim("C:/Users/ctata/Documents/Lab/GME/data/halfvarson/best_hits.tsv", header=FALSE)
#colnames(blast_hits) <- c("qseqid", "sseqid", "qseq", "sseq", "evalue", "bitscore", "length", "pident" )
#seqtab <- getExampleSeqtab_Halfvarson()
#seqtab_embedded <- embedSeqtab(seqtab = seqtab, fasta_file = "data/halfvarson/sequences.fasta", blast_hits = blast_hits)

```


## 6. Read in the hits from running blast. 
```
blast_hits <- getBlastHits("path to blast hits file")
```

## 7. Embed your sequence table
```
EmbedAsvTable(seqtab, fasta_file, blast_hits, embedding_matrix)
```
Please keep in mind that the arguments passed in to the `EmbedAsvTable` function must be named: **seqtab**, **fasta_file**, **blast_hits**, and **embedding_file**! If the names of these objects are saved as anything else the `EmbedAsvTable` function will not execute.

## Built in Example
This package comes with example **seqtab**, **fasta_file**, **blast_hits**, and **embedding_file** for you to test.
To access and read these files follow the instructions below:

a. Read in blast_hits file
```
blast_hits <- getBlastHits(system.file("extdata", "blast_hits.tsv", package = "GMEmbeddings"))
```

b. Read in embedding file
```
embedding_matrix <- getEmbeddingMatrix(system.file("extdata", "embed_.07_100dim.txt", package = "GMEmbeddings"))
```

c. Read in sequence table file
```
seqtab <- getSeqtab(system.file("extdata", "seqtab_test.csv", package = "GMEmbeddings"))
```

d. Read in fasta file
```
fasta_file <- getFasta(system.file("extdata", "fasta_test.fasta", package = "GMEmbeddings"))
```

e. Embed the sequence table
```
EmbedAsvTable(seqtab, fasta_file, blast_hits, embedding_matrix)
```
