# GMEmbeddings

## 0: Readme: will download a bunch of stuff

## 1. Install package

## 2. Select the embedding matrix you'd like to use. Standard is .07_100: 
```
embedding_file_name <- system.file("extdata", "embed_.07_100dim.txt", package = "GMEmbeddings")
```

## 3. Read in your sequence table. It should look like this:
```
seqtab <- getExampleSeqtab()
```

## 4. Read in your fasta file. ID's should be ASV ids from the column names of your sequence table, and sequences should be the full length ASV sequence. It should look like this:
```
fasta_file_name <- getExampleFasta()
```

## 5.Run blast to align your sequences to the sequences in our embedding database. Here's how:
## 5a: install blast
## 5b: download blastdb embdding sequences from cgrb: http://files.cgrb.oregonstate.edu/David_Lab/microbiome_embeddings/
## 5c: Align your sequences to those sequences ^. Here's how: *line from other comments*. It should look like this:


## 6. Read in the hits from running blast. Ultimately, it should look like this: 
```
blast_hits <- read.delim(system.file("extdata", "blast_hits.tsv", package = "GMEmbeddings"))
```

## 7. Embed your sequence table
```
EmbedAsvTable(seqtab, fasta_file_name, blast_hits, embedding_file_name)
```
