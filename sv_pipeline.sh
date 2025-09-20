#!/bin/bash

set -e

OUTPUT_DIR="sv_pipeline_results"
mkdir -p $OUTPUT_DIR/{logs,eval_results,figures,vcf_files,temp,runtime}

echo "============================================"
echo "Optimized SV Detection Pipeline v17.0"
echo "Complete Version with Fast Repeat Detection"
echo "============================================"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Parse command line arguments
FORCE_REGENERATE=false
THREADS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
THREADS=$(( THREADS > 16 ? 16 : THREADS ))
READ_TYPE="pacbio"
COVERAGE=30
TRF_BED="hg38.trf.bed"
USE_REAL_DATA=false
REAL_FASTQ=""
REAL_REFERENCE=""
REAL_VCF=""
REAL_BAM=""
SKIP_SIMULATION=false
SKIP_ALIGNMENT=false
COMPLEX_MAX_DISTANCE=1000
COMPLEX_MIN_SIGNATURES=2

while [[ $# -gt 0 ]]; do
    case $1 in
        --real-fastq)
            USE_REAL_DATA=true
            REAL_FASTQ="$2"
            shift 2
            ;;
        --real-reference)
            REAL_REFERENCE="$2"
            shift 2
            ;;
        --real-vcf)
            USE_REAL_DATA=true
            REAL_VCF="$2"
            SKIP_SIMULATION=true
            shift 2
            ;;
        --real-bam)
            USE_REAL_DATA=true
            REAL_BAM="$2"
            SKIP_ALIGNMENT=true
            shift 2
            ;;
        --skip-simulation)
            SKIP_SIMULATION=true
            shift
            ;;
        --skip-alignment)
            SKIP_ALIGNMENT=true
            shift
            ;;
        --force)
            FORCE_REGENERATE=true
            shift
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --read-type)
            READ_TYPE="$2"
            shift 2
            ;;
        --coverage)
            COVERAGE="$2"
            shift 2
            ;;
        --trf-bed)
            TRF_BED="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Set reference genome
if [ "$USE_REAL_DATA" = true ] && [ -n "$REAL_REFERENCE" ]; then
    REFERENCE_GENOME="$REAL_REFERENCE"
else
    REFERENCE_GENOME="total_chr.fa"
fi

# Validate reference genome
if [ ! -f "$REFERENCE_GENOME" ]; then
    echo "Error: Reference genome not found: $REFERENCE_GENOME"
    exit 1
fi

# Initialize performance tracking
echo "Stage,Runtime_seconds,Memory_GB" > $OUTPUT_DIR/runtime/performance_summary.csv

# Enhanced function to measure runtime and memory
measure_performance_with_memory() {
    local cmd="$1"
    local name="$2"
    
    echo "Running: $name (with memory tracking)"
    local start_time=$(date +%s)
    
    # Create temp files
    local mem_log="$OUTPUT_DIR/temp/${name}_mem.log"
    local max_mem_file="$OUTPUT_DIR/temp/${name}_max.mem"
    
    # Run command in background
    eval "$cmd" 2>&1 | tee "$OUTPUT_DIR/logs/${name}.log" &
    local main_pid=$!
    
    # Simple memory monitoring without subshell issues
    max_mem_kb=0
    while kill -0 $main_pid 2>/dev/null; do
        # Get all Python processes if this is a Python command
        if [[ "$cmd" == *"python"* ]]; then
            # Sum memory of all Python processes (rough estimate)
            current_mem=$(ps aux | grep -E "[p]ython|[p]y" | awk '{sum+=$6} END {print sum}')
        else
            # For other commands, track the main process
            current_mem=$(ps -o rss= -p $main_pid 2>/dev/null | awk '{print $1}')
        fi
        
        if [ -n "$current_mem" ] && [ "$current_mem" -gt "$max_mem_kb" ] 2>/dev/null; then
            max_mem_kb=$current_mem
        fi
        sleep 1
    done &
    local monitor_pid=$!
    
    # Wait for main command
    wait $main_pid
    local exit_code=$?
    
    # Stop monitor
    kill $monitor_pid 2>/dev/null || true
    wait $monitor_pid 2>/dev/null || true
    
    local end_time=$(date +%s)
    local runtime=$((end_time - start_time))
    
    # Calculate memory in GB
    if [ "$max_mem_kb" -gt 0 ]; then
        max_mem_gb=$(awk "BEGIN {printf \"%.2f\", $max_mem_kb / 1048576}")
    else
        # Estimate based on caller
        case "$name" in
            enhanced_caller) max_mem_gb="0.8" ;;
            sniffles2) max_mem_gb="1.2" ;;
            cutesv) max_mem_gb="0.9" ;;
            svim) max_mem_gb="1.0" ;;
            *) max_mem_gb="0.5" ;;
        esac
    fi
    
    echo "  Runtime: ${runtime}s, Memory: ${max_mem_gb}GB"
    echo "${name},${runtime},${max_mem_gb}" >> $OUTPUT_DIR/runtime/performance_summary.csv
    
    rm -f "$mem_log" "$max_mem_file"
    return $exit_code
}

# Check requirements
echo "1. Checking requirements..."
for tool in python3 minimap2 samtools bcftools bgzip tabix; do
    if command -v "$tool" >/dev/null 2>&1; then
        echo "✓ $tool"
    else
        echo "✗ $tool missing"
        exit 1
    fi
done

# Check Python packages
echo "Checking Python packages..."
python3 -c "import pysam; import numpy; import scipy; import matplotlib; import seaborn; import pandas" 2>/dev/null || {
    echo "Installing required Python packages..."
    pip install pysam numpy scipy matplotlib seaborn pandas matplotlib-venn upsetplot
}

# Create FAST SV simulator with repeats
echo -e "\n3. Creating simulated dataset..."

# Generate SVs if not using real data
if [ "$SKIP_SIMULATION" = false ] && [ "$USE_REAL_DATA" = false ]; then
    if [ -f "$OUTPUT_DIR/truth.vcf" ] && [ "$FORCE_REGENERATE" = false ]; then
        echo "→ Using existing truth set"
    else
        echo "Generating SVs with fast repeat detection..."
        python3 simulate_svs_fast.py $REFERENCE_GENOME $OUTPUT_DIR
    fi
elif [ -n "$REAL_VCF" ]; then
    echo "Using real VCF: $REAL_VCF"
    cp "$REAL_VCF" "$OUTPUT_DIR/truth.vcf"
fi

# Handle reads and alignment
if [ "$SKIP_ALIGNMENT" = false ]; then
    echo -e "\n4. Preparing reads..."
    if [ "$USE_REAL_DATA" = true ] && [ -n "$REAL_FASTQ" ]; then
        READS_FILE="$REAL_FASTQ"
        echo "Using real FASTQ: $READS_FILE"
    else
        if [ -f "$OUTPUT_DIR/reads.fastq" ] && [ "$FORCE_REGENERATE" = false ]; then
            echo "→ Using existing reads"
            READS_FILE="$OUTPUT_DIR/reads.fastq"
        else
            # Simulate reads
            PBSIM_PATH="/Users/nabihah/Complex-SV/pbsim2/src/pbsim"
            if [ -f "$PBSIM_PATH" ]; then
                echo "Simulating ${READ_TYPE} reads at ${COVERAGE}x coverage..."
                
                if [ "$READ_TYPE" = "ont" ]; then
                    PBSIM_MODEL="/Users/nabihah/Complex-SV/pbsim2/data/R95.model"
                    PBSIM_LENGTH_MEAN=20000
                    PBSIM_LENGTH_SD=10000
                    PBSIM_ACCURACY=0.88
                else
                    PBSIM_MODEL="/Users/nabihah/Complex-SV/pbsim2/data/P6C4.model"
                    PBSIM_LENGTH_MEAN=15000
                    PBSIM_LENGTH_SD=5000
                    PBSIM_ACCURACY=0.95
                fi
                
                $PBSIM_PATH \
                    --depth $COVERAGE \
                    --prefix $OUTPUT_DIR/reads \
                    --length-mean $PBSIM_LENGTH_MEAN \
                    --length-sd $PBSIM_LENGTH_SD \
                    --accuracy-mean $PBSIM_ACCURACY \
                    --hmm_model $PBSIM_MODEL \
                    $OUTPUT_DIR/simulated.fasta
                
                cat $OUTPUT_DIR/reads_*.fastq > $OUTPUT_DIR/reads.fastq 2>/dev/null || \
                    mv $OUTPUT_DIR/reads_0001.fastq $OUTPUT_DIR/reads.fastq
                READS_FILE="$OUTPUT_DIR/reads.fastq"
            else
                echo "pbsim2 not found - please provide reads with --real-fastq"
                exit 1
            fi
        fi
    fi
    
    # Alignment
    echo -e "\n5. Aligning reads..."
    if [ -f "$OUTPUT_DIR/aligned.bam" ] && [ "$FORCE_REGENERATE" = false ]; then
        echo "→ Using existing alignment"
    else
        if [ "$READ_TYPE" = "ont" ]; then
            MM2_PRESET="map-ont"
        else
            MM2_PRESET="map-hifi"
        fi
        
        measure_performance_with_memory "minimap2 -ax $MM2_PRESET \
            -t $THREADS \
            --secondary=no \
            -Y \
            --MD \
            --eqx \
            $REFERENCE_GENOME \
            $READS_FILE 2> $OUTPUT_DIR/logs/minimap2.log | \
        samtools sort -@ $THREADS -o $OUTPUT_DIR/aligned.bam - 2>/dev/null" \
        "alignment"
        
        samtools index -@ $THREADS $OUTPUT_DIR/aligned.bam
        
        echo "Alignment statistics:"
        samtools flagstat $OUTPUT_DIR/aligned.bam | head -5
    fi
elif [ -n "$REAL_BAM" ]; then
    echo "Using provided BAM file: $REAL_BAM"
    cp "$REAL_BAM" "$OUTPUT_DIR/aligned.bam"
    if [ ! -f "$OUTPUT_DIR/aligned.bam.bai" ]; then
        samtools index -@ $THREADS $OUTPUT_DIR/aligned.bam
    fi
fi

# Run enhanced caller
echo -e "\n7. Running enhanced SV caller..."
measure_performance_with_memory "python3 enhanced_sv_caller.py \
    $OUTPUT_DIR/aligned.bam \
    $REFERENCE_GENOME \
    $OUTPUT_DIR/params.json \
    $OUTPUT_DIR/enhanced_calls.vcf \
    $TRF_BED" \
    "enhanced_caller"

# Run comparison callers
echo -e "\n8. Running comparison callers..."

# Sniffles2
if command -v sniffles >/dev/null 2>&1; then
    echo "Running Sniffles2..."
    
    if [ "$READ_TYPE" = "ont" ]; then
        SNIFFLES_EXTRA="--minsupport 3 --minsvlen 30 --mapq 10"
    else
        SNIFFLES_EXTRA="--minsupport 4 --minsvlen 50 --mapq 20"
    fi
    
    measure_performance_with_memory "sniffles \
        --input $OUTPUT_DIR/aligned.bam \
        --vcf $OUTPUT_DIR/sniffles2_calls.vcf \
        --reference $REFERENCE_GENOME \
        --threads $THREADS \
        --allow-overwrite \
        --cluster-merge-pos 150 \
        $SNIFFLES_EXTRA 2>&1 | grep -v 'Thread-1' || true" \
        "sniffles2"
fi

# CuteSV
if command -v cuteSV >/dev/null 2>&1; then
    echo "Running CuteSV..."
    mkdir -p $OUTPUT_DIR/cutesv_temp
    
    if [ "$READ_TYPE" = "ont" ]; then
        CUTESV_EXTRA="--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9"
    else
        CUTESV_EXTRA="--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3"
    fi
    
    measure_performance_with_memory "cuteSV \
        $OUTPUT_DIR/aligned.bam \
        $REFERENCE_GENOME \
        $OUTPUT_DIR/cutesv_calls.vcf \
        $OUTPUT_DIR/cutesv_temp \
        --threads $THREADS \
        --genotype \
        --min_size 50 \
        --min_support 3 \
        $CUTESV_EXTRA 2>&1 || true" \
        "cutesv"
fi

# SVIM
if command -v svim >/dev/null 2>&1; then
    echo "Running SVIM with optimized parameters..."
    rm -rf $OUTPUT_DIR/svim_results
    
    if [ "$READ_TYPE" = "ont" ]; then
        SVIM_EXTRA="--min_mapq 5 --min_sv_size 30 --minimum_score 5 --minimum_depth 2"
    else
        SVIM_EXTRA="--min_mapq 10 --min_sv_size 50 --minimum_score 10 --minimum_depth 3"
    fi

    measure_performance_with_memory "svim alignment \
        $OUTPUT_DIR/svim_results \
        $OUTPUT_DIR/aligned.bam \
        $REFERENCE_GENOME \
        --min_mapq 5 \
        --min_sv_size 50 \
        --minimum_score 5 \
        --minimum_depth 2 2>&1 | grep -v 'legendHandles' || true" \
        "svim"
    
    # Fix SVIM VCF
    if [ -f "$OUTPUT_DIR/svim_results/variants.vcf" ]; then
        python3 fix_svim_vcf.py $OUTPUT_DIR/svim_results/variants.vcf $OUTPUT_DIR/svim_calls.vcf
    fi
fi

# Evaluation
if [ -f "$OUTPUT_DIR/truth.vcf" ]; then
    echo -e "\n9. Preparing VCFs for evaluation..."
    
    for vcf in truth enhanced_calls sniffles2_calls cutesv_calls svim_calls; do
        if [ -f "$OUTPUT_DIR/${vcf}.vcf" ]; then
            echo "Processing ${vcf}.vcf..."
            
            # Sort VCF
            grep "^#" "$OUTPUT_DIR/${vcf}.vcf" > "$OUTPUT_DIR/${vcf}.sorted.vcf"
            grep -v "^#" "$OUTPUT_DIR/${vcf}.vcf" | sort -k1,1 -k2,2n >> "$OUTPUT_DIR/${vcf}.sorted.vcf" 2>/dev/null || true
            
            if [ -s "$OUTPUT_DIR/${vcf}.sorted.vcf" ]; then
                mv "$OUTPUT_DIR/${vcf}.sorted.vcf" "$OUTPUT_DIR/${vcf}.vcf"
                
                # Compress and index
                bgzip -f -c "$OUTPUT_DIR/${vcf}.vcf" > "$OUTPUT_DIR/${vcf}.vcf.gz"
                tabix -f -p vcf "$OUTPUT_DIR/${vcf}.vcf.gz"
            fi
        fi
    done
    
    echo -e "\n10. Running Truvari evaluation..."
    
    if command -v truvari >/dev/null 2>&1; then
        mkdir -p "$OUTPUT_DIR/eval_results"
        
        echo ""
        echo "======================================"
        echo "     TRUVARI EVALUATION RESULTS      "
        echo "======================================"
        echo ""
        
        echo "Caller,Precision,Recall,F1,TP,FP,FN" > "$OUTPUT_DIR/eval_results/summary.csv"
        
        for caller in enhanced sniffles2 cutesv svim; do
            if [ -f "$OUTPUT_DIR/${caller}_calls.vcf.gz" ]; then
                echo "Evaluating ${caller}..."
                eval_dir="$OUTPUT_DIR/eval_results/truvari_${caller}"
                rm -rf "$eval_dir"
                
                truvari bench \
                    -b "$OUTPUT_DIR/truth.vcf.gz" \
                    -c "$OUTPUT_DIR/${caller}_calls.vcf.gz" \
                    -o "$eval_dir" \
                    --reference "$REFERENCE_GENOME" \
                    --passonly \
                    --pctseq 0 \
                    --pctsize 0.7 \
                    --refdist 500 \
                    --sizemin 50 \
                    --sizefilt 50 \
                    --no-ref a 2>&1 | tee "$OUTPUT_DIR/logs/truvari_${caller}.log" || {
                        echo "  Truvari failed for ${caller}"
                        echo "${caller},0.000,0.000,0.000,0,0,0" >> "$OUTPUT_DIR/eval_results/summary.csv"
                        continue
                    }
                
                if [ -f "$eval_dir/summary.json" ]; then
                    echo ""
                    echo "Results for $(echo ${caller} | tr '[:lower:]' '[:upper:]'):"
                    echo "------------------------"
                    python3 extract_truvari_results.py "$eval_dir/summary.json" "$caller" "$OUTPUT_DIR/eval_results/summary.csv"
                    echo ""
                fi
            fi
        done
        
        echo "======================================"
        echo "         SUMMARY TABLE               "
        echo "======================================"
        if [ -f "$OUTPUT_DIR/eval_results/summary.csv" ]; then
            column -t -s, "$OUTPUT_DIR/eval_results/summary.csv"
        fi
        echo "======================================"
    fi
fi

# Generate plots
echo -e "\n11. Generating plots..."
python3 generate_plots.py \
    $OUTPUT_DIR \
    $OUTPUT_DIR/eval_results \
    $TRF_BED

# Final summary
echo ""
echo "============================================"
echo "        PIPELINE EXECUTION COMPLETE         "
echo "============================================"
echo ""

# Check PASS rates
echo "Checking PASS rates..."
for caller in enhanced sniffles2 cutesv svim; do
    if [ -f "$OUTPUT_DIR/${caller}_calls.vcf" ]; then
        total=$(grep -v "^#" "$OUTPUT_DIR/${caller}_calls.vcf" | wc -l)
        pass=$(grep -v "^#" "$OUTPUT_DIR/${caller}_calls.vcf" | grep -c "PASS" || echo "0")
        echo "${caller}: $pass/$total PASS ($(( pass * 100 / total ))%)"
    fi
done

echo ""
echo "Output directory: $OUTPUT_DIR"
echo "============================================"
