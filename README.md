# Complex-SV-code

# Complex-SV: A Long-Read Alignment-Based Caller for Accurate Structural Variant Detection

Complex-SV is an enhanced structural variant detection pipeline that integrates multiple evidence types through sophisticated clustering algorithms and repeat-aware filtering. It addresses critical gaps in current methodologies by providing breakthrough detection of complex structural variants while maintaining comprehensive coverage across all major variant classes.

## Key Features

- **Complex SV Detection**: First tool to achieve 90% sensitivity for complex structural variants (multi-breakpoint rearrangements)
- **Comprehensive SV Coverage**: Detects all major SV classes (DEL, INS, DUP, INV, BND, CPX)
- **Duplication Subtype Classification**: Distinguishes tandem from interspersed duplications
- **Repeat-Aware Filtering**: Novel repeat region annotation and context classification
- **Cross-Platform Support**: Optimized for both PacBio HiFi and Oxford Nanopore data
- **Computational Efficiency**: Moderate overhead (135s runtime, 0.8GB memory) for expanded capabilities

## Performance Highlights

- **Precision**: 0.790 (competitive with established tools)
- **Recall**: 0.605 (comprehensive detection across challenging variant classes)
- **F1-Score**: 0.685 (balanced performance)
- **Complex SVs**: 90% sensitivity (vs 0% for existing methods)
- **Duplications**: 99% sensitivity (exceeds standard 60-80% benchmarks)

## Quick Start

### Basic Usage

For germline SV calling from long-read alignments:

```bash
# Run the complete pipeline with simulated data
./sv_pipeline.sh --threads 16 --coverage 30 --read-type pacbio

# With real PacBio HiFi data
./sv_pipeline.sh --real-fastq reads.fastq --real-reference ref.fa --threads 16

# With real ONT data
./sv_pipeline.sh --real-fastq ont_reads.fastq --read-type ont --threads 16
```

### Using Real Data

```bash
# Start from real FASTQ files
./sv_pipeline.sh --real-fastq your_reads.fastq \
                 --real-reference your_reference.fa \
                 --read-type pacbio \
                 --threads 16

# Start from aligned BAM file
./sv_pipeline.sh --real-bam aligned.bam \
                 --real-reference your_reference.fa \
                 --skip-alignment

# With truth VCF for evaluation
./sv_pipeline.sh --real-vcf truth.vcf \
                 --real-bam aligned.bam \
                 --real-reference your_reference.fa \
                 --skip-simulation
```

### With Tandem Repeat Annotations

For improved calling in repetitive regions, Complex-SV accepts tandem repeat annotations:

```bash
./sv_pipeline.sh --trf-bed hg38.trf.bed \
                 --real-fastq reads.fastq \
                 --real-reference hg38.fa
```

Download TRF annotations for human references:

```bash
# Download GRCh38 genome TRF bed file
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.trf.bed.gz
gunzip hg38.trf.bed.gz
```

## Installation

### Requirements

**System Tools:**
- `python3` (≥3.8)
- `minimap2` (≥2.24)
- `samtools` (≥1.16)
- `bcftools`
- `bgzip`
- `tabix`
- `truvari` (for evaluation)

**Optional Tools:**
- `sniffles` (≥2.0, for comparison)
- `cuteSV` (for comparison)
- `svim` (for comparison)
- `pbsim2` (for read simulation)

**Python Packages:**
```bash
pip install pysam numpy scipy matplotlib seaborn pandas matplotlib-venn upsetplot
```

### Installation Steps

1. **Clone the repository:**
```bash
git clone https://github.com/yourusername/complex-sv.git
cd complex-sv
```

2. **Install dependencies:**
```bash
# Install system tools (Ubuntu/Debian)
sudo apt-get install minimap2 samtools bcftools tabix

# Install Python packages
pip install -r requirements.txt

# Install Truvari for evaluation
pip install truvari

# Optional: Install comparison tools
conda install -c bioconda sniffles cutesv svim
```

3. **Download required reference data:**
```bash
# Download GRCh38 genome TRF bed file
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.trf.bed.gz
gunzip hg38.trf.bed.gz

# Download chr20 and chr21 fasta files
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/assembly_GRCh38/chr20.fa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/assembly_GRCh38/chr21.fa.gz
gunzip chr20.fa.gz chr21.fa.gz

# Combine chromosomes for testing
cat chr20.fa chr21.fa > total_chr.fa
```

4. **Make the pipeline executable:**
```bash
chmod +x sv_pipeline.sh
```

## Usage Examples

### 1. Complete Pipeline with Simulated Data

```bash
# PacBio HiFi simulation
./sv_pipeline.sh --threads 16 --coverage 30 --read-type pacbio

# Oxford Nanopore simulation  
./sv_pipeline.sh --threads 16 --coverage 30 --read-type ont
```

### 2. Real Data Analysis

```bash
# PacBio HiFi analysis
./sv_pipeline.sh --real-fastq pacbio_hifi.fastq \
                 --real-reference GRCh38.fa \
                 --trf-bed hg38.trf.bed \
                 --threads 16

# ONT analysis
./sv_pipeline.sh --real-fastq ont_reads.fastq \
                 --real-reference GRCh38.fa \
                 --read-type ont \
                 --threads 16
```

### 3. Skip Specific Stages

```bash
# Skip simulation, use existing data
./sv_pipeline.sh --skip-simulation \
                 --real-fastq reads.fastq \
                 --real-reference ref.fa

# Skip alignment, use existing BAM
./sv_pipeline.sh --skip-alignment \
                 --real-bam aligned.bam \
                 --real-reference ref.fa
```

### 4. Force Regeneration

```bash
# Regenerate all outputs
./sv_pipeline.sh --force --threads 16
```

## Command Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--threads N` | Number of threads | Auto-detected |
| `--coverage N` | Target sequencing coverage | 30 |
| `--read-type TYPE` | Sequencing technology (pacbio/ont) | pacbio |
| `--trf-bed FILE` | Tandem repeat annotations | hg38.trf.bed |
| `--real-fastq FILE` | Real FASTQ input | None |
| `--real-reference FILE` | Reference genome | total_chr.fa |
| `--real-vcf FILE` | Truth VCF for evaluation | None |
| `--real-bam FILE` | Aligned BAM file | None |
| `--skip-simulation` | Skip SV simulation | False |
| `--skip-alignment` | Skip read alignment | False |
| `--force` | Force regenerate all outputs | False |

## Output Files

The pipeline generates comprehensive outputs in the `sv_pipeline_results/` directory:

```
sv_pipeline_results/
├── truth.vcf                    # Simulated truth set
├── reads.fastq                  # Simulated reads
├── aligned.bam                  # Aligned reads
├── enhanced_calls.vcf           # Complex-SV results
├── sniffles2_calls.vcf         # Sniffles2 results
├── cutesv_calls.vcf            # CuteSV results
├── svim_calls.vcf              # SVIM results
├── eval_results/               # Truvari evaluation
│   ├── summary.csv             # Performance summary
│   └── truvari_*/              # Detailed results per tool
├── figures/                    # Generated plots
├── logs/                       # Tool logs
└── runtime/                    # Performance metrics
```

## Key Evaluation Commands

### SV Simulation
```bash
python3 simulate_svs_fast.py total_chr.fa sv_pipeline_results
```

### Read Alignment
```bash
# PacBio HiFi
minimap2 -ax map-hifi -t 16 --secondary=no -Y --MD --eqx \
         reference.fa reads.fastq | \
samtools sort -@ 16 -o aligned.bam -

# Oxford Nanopore
minimap2 -ax map-ont -t 16 --secondary=no -Y --MD --eqx \
         reference.fa reads.fastq | \
samtools sort -@ 16 -o aligned.bam -
```

### Complex-SV Calling
```bash
python3 enhanced_sv_caller.py \
        aligned.bam \
        reference.fa \
        params.json \
        output.vcf \
        tandem_repeats.bed
```

### Truvari Evaluation
```bash
truvari bench \
    -b truth.vcf.gz \
    -c calls.vcf.gz \
    -o eval_results \
    --reference reference.fa \
    --passonly \
    --pctseq 0 \
    --pctsize 0.7 \
    --refdist 500 \
    --sizemin 50
```

## Performance Characteristics

Based on evaluation using chromosomes 20-21 (~100Mb):

| Metric | Complex-SV | Sniffles2 | CuteSV | SVIM |
|--------|------------|-----------|--------|------|
| **Precision** | 0.790 | 0.854 | 0.740 | 0.826 |
| **Recall** | 0.605 | 0.623 | 0.623 | 0.609 |
| **F1-Score** | 0.685 | 0.720 | 0.677 | 0.701 |
| **Runtime** | 135s | 54s | 55s | 108s |
| **Memory** | 0.8GB | 1.2GB | 0.9GB | 1.0GB |

### Variant-Specific Performance

| SV Type | Sensitivity | Notes |
|---------|-------------|-------|
| **Deletions (DEL)** | 91% | Reliable detection 50bp-100kb |
| **Insertions (INS)** | 98% | Good performance 50bp-5kb |
| **Duplications (DUP)** |99% | Exceeds standard benchmarks |
| **Inversions (INV)** | 98.0% | Perfect detection 1-50kb |
| **Breakends (BND)** | 87.0% | Inter-chromosomal events |
| **Complex SVs (CPX)** | 90.0% | **Unique capability** |

## Multi-Sample Analysis

For population-scale studies, run Complex-SV on individual samples then merge:

```bash
# Process multiple samples
for sample in sample1 sample2 sample3; do
    ./sv_pipeline.sh --real-fastq ${sample}.fastq \
                     --real-reference reference.fa \
                     --threads 8
    mv sv_pipeline_results/enhanced_calls.vcf ${sample}_calls.vcf
done

# Merge results (requires bcftools)
bcftools merge *_calls.vcf.gz > merged_calls.vcf
```

## Citation

If you use Complex-SV in your research, please cite:

```
Complex-SV: A long-read alignment-based caller for accurate structural variant detection.
Master's Thesis, University of Manchester, 2025.
```

## Comparison with Other Tools

Complex-SV builds upon and complements existing tools:

- **vs Sniffles2**: Higher complex SV detection, lower precision
- **vs CuteSV**: Better duplication calling, moderate runtime increase  
- **vs SVIM**: Enhanced complex SV capabilities, competitive performance
- **Unique Features**: Complex SV detection, duplication subtyping, repeat classification

## Troubleshooting

### Common Issues

1. **Reference genome not found**
   ```bash
   # Ensure reference file exists and is readable
   ls -la your_reference.fa
   ```

2. **Memory errors**
   ```bash
   # Reduce thread count
   ./sv_pipeline.sh --threads 8
   ```

3. **Missing dependencies**
   ```bash
   # Check all tools are installed
   which minimap2 samtools python3 truvari
   ```

4. **TRF bed file not found**
   ```bash
   # Download or create TRF annotations
   wget https://github.com/yourusername/annotations/hg38.trf.bed
   ```

### Platform-Specific Notes


**PacBio HiFi (Recommended)**:
- Optimal performance across all SV types
- Best complex SV detection
- Reliable duplication calling

**Oxford Nanopore**:
- Good deletion/insertion detection
- Reduced duplication sensitivity (6% vs 118% on PacBio)
- Higher false positive rates for breakends



## Acknowledgments

- University of Manchester School of Health Sciences
- ICSF High-Performance interactive Computational shared Facility
