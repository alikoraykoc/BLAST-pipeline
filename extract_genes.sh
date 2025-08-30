#!/bin/bash

# Gene extraction pipeline from genome assemblies
# Usage: ./extract_genes.sh REFERENCE_GENES ASSEMBLIES_DIR OUTPUT_DIR GENE_NAME [EMAIL]

REFERENCE_GENES=$1
ASSEMBLIES_DIR=$2
OUTPUT_DIR=$3
GENE_NAME=$4
EMAIL=$5

# Check arguments
if [ $# -lt 4 ]; then
    echo "Usage: $0 REFERENCE_GENES ASSEMBLIES_DIR OUTPUT_DIR GENE_NAME [EMAIL]"
    echo "Example: $0 rbcL_refs.fasta ./assemblies/ ./results/ rbcL user@email.com"
    exit 1
fi

# Set default email if not provided
if [ -z "$EMAIL" ]; then
    EMAIL="dummy@example.com"
fi

# Check if required tools are available
command -v makeblastdb >/dev/null 2>&1 || { echo "❌ makeblastdb not found. Install BLAST+"; exit 1; }
command -v blastn >/dev/null 2>&1 || { echo "❌ blastn not found. Install BLAST+"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "❌ samtools not found. Install samtools"; exit 1; }
command -v seqtk >/dev/null 2>&1 || { echo "❌ seqtk not found. Install seqtk"; exit 1; }

# Check input files
if [ ! -f "$REFERENCE_GENES" ]; then
    echo "❌ Reference genes file not found: $REFERENCE_GENES"
    exit 1
fi

if [ ! -d "$ASSEMBLIES_DIR" ]; then
    echo "❌ Assemblies directory not found: $ASSEMBLIES_DIR"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "🧬 Gene Extraction Pipeline Starting..."
echo "🎯 Target gene: $GENE_NAME"
echo "📚 Reference file: $REFERENCE_GENES"
echo "📁 Assemblies directory: $ASSEMBLIES_DIR"
echo "💾 Output directory: $OUTPUT_DIR"
echo ""

# Count reference sequences
ref_count=$(grep -c ">" "$REFERENCE_GENES" 2>/dev/null || echo 0)
echo "📊 Found $ref_count reference sequences"

# Count assembly files
assembly_count=$(find "$ASSEMBLIES_DIR" -name "*.fasta" -o -name "*.fa" -o -name "*.fna" | wc -l)
echo "📊 Found $assembly_count assembly files"
echo ""

# Initialize counters
success_count=0
failed_count=0
log_file="$OUTPUT_DIR/extraction_log.txt"

# Clear log file
echo "Gene Extraction Log - $(date)" > "$log_file"
echo "=================================" >> "$log_file"

# Process each assembly
for assembly in "$ASSEMBLIES_DIR"/*.{fasta,fa,fna}; do
    # Skip if no files match the pattern
    [ ! -f "$assembly" ] && continue

    # Extract taxon name from filename
    TAXON=$(basename "$assembly" | sed 's/\.\(fasta\|fa\|fna\)$//')

    echo "🔍 Processing: $TAXON"

    # Create BLAST database for this assembly
    db_path="$OUTPUT_DIR/${TAXON}_db"
    makeblastdb -in "$assembly" -dbtype nucl -out "$db_path" -logfile /dev/null 2>/dev/null

    if [ $? -ne 0 ]; then
        echo "  ❌ Failed to create BLAST database"
        echo "$TAXON: Failed to create BLAST database" >> "$log_file"
        failed_count=$((failed_count + 1))
        continue
    fi

    # Run BLAST search
    blast_output="$OUTPUT_DIR/${TAXON}_blast.tsv"
    blastn -query "$REFERENCE_GENES" \
           -db "$db_path" \
           -out "$blast_output" \
           -outfmt "6 qseqid sseqid pident qcovs evalue bitscore sstart send sstrand" \
           -evalue 1e-10 \
           -max_target_seqs 5 \
           2>/dev/null

    if [ $? -ne 0 ]; then
        echo "  ❌ BLAST search failed"
        echo "$TAXON: BLAST search failed" >> "$log_file"
        failed_count=$((failed_count + 1))
        rm -f "${db_path}".*
        continue
    fi

    # Check if we have any hits
    if [ ! -s "$blast_output" ]; then
        echo "  ❌ No BLAST hits found"
        echo "$TAXON: No BLAST hits found" >> "$log_file"
        failed_count=$((failed_count + 1))
        rm -f "${db_path}".*
        rm -f "$blast_output"
        continue
    fi

    # Get best hit (highest bitscore)
    best_hit=$(sort -k6,6nr "$blast_output" | head -1)

    if [ -z "$best_hit" ]; then
        echo "  ❌ No valid hits found"
        echo "$TAXON: No valid hits found" >> "$log_file"
        failed_count=$((failed_count + 1))
        rm -f "${db_path}".*
        rm -f "$blast_output"
        continue
    fi

    # Parse best hit information
    IFS=$'\t' read -r qseqid sseqid pident qcovs evalue bitscore sstart send sstrand <<< "$best_hit"

    # Quality control checks
    identity_ok=$(echo "$pident >= 70" | bc -l 2>/dev/null || echo 0)
    coverage_ok=$(echo "$qcovs >= 70" | bc -l 2>/dev/null || echo 0)

    if [ "$identity_ok" -eq 0 ] || [ "$coverage_ok" -eq 0 ]; then
        printf "  ❌ Poor quality hit (%.1f%% identity, %.1f%% coverage)\n" "$pident" "$qcovs"
        echo "$TAXON: Poor quality hit ($pident% identity, $qcovs% coverage)" >> "$log_file"
        failed_count=$((failed_count + 1))
        rm -f "${db_path}".*
        rm -f "$blast_output"
        continue
    fi

    # Index the assembly for samtools
    samtools faidx "$assembly" 2>/dev/null

    if [ $? -ne 0 ]; then
        echo "  ❌ Failed to index assembly"
        echo "$TAXON: Failed to index assembly" >> "$log_file"
        failed_count=$((failed_count + 1))
        rm -f "${db_path}".*
        rm -f "$blast_output"
        continue
    fi

    # Handle coordinate order (start might be > end for reverse hits)
    if [ "$sstart" -gt "$send" ]; then
        coord_start=$send
        coord_end=$sstart
    else
        coord_start=$sstart
        coord_end=$send
    fi

    # Extract sequence using samtools
    output_fasta="$OUTPUT_DIR/${DISPLAY_NAME}_${GENE_NAME}.fasta"
    samtools faidx "$assembly" "${sseqid}:${coord_start}-${coord_end}" 2>/dev/null | \
    sed "s/>.*/>${DISPLAY_NAME}_${GENE_NAME}/" > "$output_fasta"

    if [ $? -ne 0 ] || [ ! -s "$output_fasta" ]; then
        echo "  ❌ Failed to extract sequence"
        echo "$TAXON: Failed to extract sequence" >> "$log_file"
        failed_count=$((failed_count + 1))
        rm -f "${db_path}".*
        rm -f "$blast_output"
        rm -f "${assembly}.fai"
        continue
    fi

    # Handle reverse complement if needed
    if [[ "$sstrand" == "minus" ]] || [[ "$sstart" -gt "$send" ]]; then
        seqtk seq -r "$output_fasta" > "${output_fasta}.tmp"
        mv "${output_fasta}.tmp" "$output_fasta"

        if [ $? -ne 0 ]; then
            echo "  ❌ Failed to reverse complement sequence"
            echo "$TAXON: Failed to reverse complement sequence" >> "$log_file"
            failed_count=$((failed_count + 1))
            rm -f "${db_path}".*
            rm -f "$blast_output"
            rm -f "${assembly}.fai"
            rm -f "$output_fasta"
            continue
        fi
    fi

    # Verify extracted sequence
    seq_length=$(grep -v ">" "$output_fasta" | tr -d '\n' | wc -c)

    if [ "$seq_length" -lt 50 ]; then
        echo "  ❌ Extracted sequence too short ($seq_length bp)"
        echo "$TAXON: Extracted sequence too short ($seq_length bp)" >> "$log_file"
        failed_count=$((failed_count + 1))
        rm -f "$output_fasta"
    else
        printf "  ✅ Success (%.1f%% ID, %.1f%% cov, %d bp)\n" "$pident" "$qcovs" "$seq_length"
        echo "$TAXON: Success ($pident% ID, $qcovs% cov, $seq_length bp) - $sseqid:$coord_start-$coord_end" >> "$log_file"
        success_count=$((success_count + 1))
    fi

    # Clean up temporary files
    rm -f "${db_path}".*
    rm -f "$blast_output"
    rm -f "${assembly}.fai"

done

# Combine all extracted sequences
echo ""
echo "🔗 Combining extracted sequences..."
combined_file="$OUTPUT_DIR/all_${GENE_NAME}_extracted.fasta"

if ls "$OUTPUT_DIR"/*_${GENE_NAME}.fasta 1> /dev/null 2>&1; then
    cat "$OUTPUT_DIR"/*_${GENE_NAME}.fasta > "$combined_file" 2>/dev/null
    combined_count=$(grep -c ">" "$combined_file" 2>/dev/null || echo 0)
    echo "✅ Combined file created: $combined_file ($combined_count sequences)"
    echo "" >> "$log_file"
    echo "SUMMARY:" >> "$log_file"
    echo "Combined $combined_count sequences into $combined_file" >> "$log_file"
else
    echo "⚠️  No sequences extracted to combine"
    echo "" >> "$log_file"
    echo "SUMMARY: No sequences extracted" >> "$log_file"
    combined_count=0
fi

# Final summary
echo ""
echo "🎉 Extraction completed!"
echo "✅ Successful extractions: $success_count"
echo "❌ Failed extractions: $failed_count"
echo "📊 Success rate: $(echo "scale=1; $success_count * 100 / ($success_count + $failed_count)" | bc -l 2>/dev/null || echo "N/A")%"
echo "📂 Results saved in: $OUTPUT_DIR"
echo "📋 Log file: $log_file"

# Final log entry
echo "SUCCESS: $success_count, FAILED: $failed_count" >> "$log_file"
echo "Completed at: $(date)" >> "$log_file"

exit 0