# Complex-SV-code

Download GRCh38 genome TRF bed file from code: 
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.trf.bed.gz
gunzip -c hg38.trf.bed.gz

Download chr20 and chr21 fasta files from code:
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/assembly_GRCh38/chr20.fa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/assembly_GRCh38/chr21.fa.gz
gunzip chr20.fa.gz chr21.fa.gz
cat chr20.fa chr21.fa > total_chr.fa
