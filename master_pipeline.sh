#!/bin/bash

set -e

# Configuration
GENES=("CLOCK" "ARNTL" "PER1" "PER2" "CRY1" "CRY2" "NR1D1" "RORA")
INPUT_FILE=""  # Can be species or accessions
EMAIL=""
TAXONOMIC_GROUPS=("mammals" "birds" "vertebrates")  # Default groups

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

print_status() {
    echo -e "${BLUE}[$(date +'%H:%M:%S')]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[$(date +'%H:%M:%S')] ✅${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[$(date +'%H:%M:%S')] ⚠️${NC} $1"
}

print_error() {
    echo -e "${RED}[$(date +'%H:%M:%S')] ❌${NC} $1"
}

show_help() {
    echo "🧬 BLAST PIPELINE v1.0"
    echo "=============================="
    echo ""
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --email EMAIL           Email for NCBI API (required)"
    echo "  --input FILE            Input file: species TSV or accessions list"
    echo "  --genes GENE1,GENE2     Gene names, accession lists, or file paths"
    echo "  --taxa GROUP1,GROUP2    Comma-separated taxonomic groups (default: mammals,birds,vertebrates)"
    echo "  --list-taxa            Show available taxonomic groups"
    echo "  --references-only      Only collect reference sequences"
    echo "  --assemblies-only      Only download assemblies"
    echo "  --extract-only         Only extract genes (refs and assemblies must exist)"
    echo "  --help                 Show this help"
    echo ""
    echo "Input file formats:"
    echo "  Species: Genus\\tspecies\\ttag"
    echo "  Accessions: GCF_000001405.40 or GCA_000001405.29"
    echo "  Mixed format also supported"
    echo ""
    echo "Taxonomic Groups:"
    echo "  mammals, birds, reptiles, fish, arthropods, insects,"
    echo "  plants, fungi, bacteria, vertebrates, primates, rodents"
    echo ""
    echo "Examples:"
    echo "  $0 --email user@example.com --input accessions.txt"
    echo "  $0 --email user@example.com --genes CLOCK,PER1 --taxa mammals,birds"
    echo "  $0 --email user@example.com --genes \"MT536114.1 MT536113.1\" --input species.tsv"
    echo "  $0 --email user@example.com --genes 18S --taxa 6993,Orthoptera"
}


get_gene_name_interactively() {
    local gene_input="$1"
    local suggested_name=""

    # Try to suggest a name based on filename
    if [ -f "$gene_input" ]; then
        local basename_file=$(basename "$gene_input" | sed 's/\.[^.]*$//')

        # Look for common gene patterns
        if [[ "$basename_file" =~ ([0-9]+[sS]) ]]; then
            suggested_name="${BASH_REMATCH[1]}"
        elif [[ "$basename_file" =~ (28[sS]|18[sS]|16[sS]|ITS|COI|rbcL|matK|CLOCK|PER1|PER2|CRY1|CRY2|ARNTL|NR1D1|RORA) ]]; then
            suggested_name="${BASH_REMATCH[1]}"
        elif [[ "$basename_file" =~ [a-zA-Z]+$ ]]; then
            suggested_name=$(echo "$basename_file" | sed 's/.*[_-]//')
        else
            suggested_name="$basename_file"
        fi

        # Clean up suggestion
        suggested_name=$(echo "$suggested_name" | tr '[:lower:]' '[:upper:]')

        # FIXED: Use simple echo instead of print_status
        echo "📁 Processing file: $gene_input" >&2
        echo "💡 Suggested gene name: $suggested_name" >&2
        echo "" >&2

        while true; do
            if [ -n "$suggested_name" ]; then
                read -p "🧬 Enter gene name (or press Enter for '$suggested_name'): " user_input
                if [ -z "$user_input" ]; then
                    # FIXED: Return ONLY the gene name, no extra text
                    echo "$suggested_name"
                    return 0
                else
                    # FIXED: Return ONLY the user input, clean it
                    clean_input=$(echo "$user_input" | tr -d '[:space:]' | tr -d '[:cntrl:]')
                    echo "$clean_input"
                    return 0
                fi
            else
                read -p "🧬 Enter gene name for $gene_input: " user_input
                if [ -n "$user_input" ]; then
                    clean_input=$(echo "$user_input" | tr -d '[:space:]' | tr -d '[:cntrl:]')
                    echo "$clean_input"
                    return 0
                else
                    echo "Gene name cannot be empty. Please enter a valid gene name." >&2
                fi
            fi
        done
    else
        # For accession lists or direct gene names
        read -p "🧬 Enter gene name for '$gene_input': " user_input
        if [ -n "$user_input" ]; then
            clean_input=$(echo "$user_input" | tr -d '[:space:]' | tr -d '[:cntrl:]')
            echo "$clean_input"
        else
            echo "custom_gene"  # Fallback
        fi
    fi
}


list_taxonomic_groups() {
    print_status "Available taxonomic groups:"
    python3 collect_gene_by_taxon.py --list_taxa --email "$EMAIL" 2>/dev/null || {
        echo "mammals - Mammals"
        echo "birds - Birds"
        echo "vertebrates - Vertebrates"
        echo "arthropods - Arthropods"
        echo "insects - Insects"
        echo "plants - Plants"
        echo "fungi - Fungi"
        echo "or use GenBank taxon code"
    }
}

check_dependencies() {
    print_status "Checking dependencies..."

    local missing=0

    # Python packages
    if ! python3 -c "import Bio" 2>/dev/null; then
        print_error "BioPython not found. Install with: pip install biopython"
        missing=1
    fi

    # BLAST tools
    if ! command -v makeblastdb >/dev/null 2>&1; then
        print_warning "makeblastdb not found. Install BLAST+ for gene extraction"
    fi

    if [ "$missing" -eq 1 ]; then
        exit 1
    fi

    print_success "Dependencies checked"
}

setup_directories() {
    print_status "Setting up directories..."

    mkdir -p assemblies references results logs

    for gene in "${GENES[@]}"; do
        mkdir -p "results/$gene"
    done

    print_success "Directory structure ready"
}

collect_reference_genes() {
    echo "🎯 Gene name detection: Interactive mode"
    echo ""

    local total_collected=0

    for gene in "${GENES[@]}"; do
        local gene_name

        # Check if gene is a file path
        if [ -f "$gene" ]; then
            print_status "📋 Processing accession file: $gene"

            # FIXED: Capture gene name cleanly
            gene_name=$(get_gene_name_interactively "$gene")

            # FIXED: Clean any remaining control characters
            gene_name=$(echo "$gene_name" | tr -d '[:cntrl:]' | tr -d '\033' | sed 's/\[[0-9;]*[mGK]//g')

            print_success "Using gene name: $gene_name"
            echo ""

            # FIXED: Pass clean gene name to Python script
            python3 accession_reference_downloader.py \
                --accessions "$gene" \
                --gene "$gene_name" \
                --email "$EMAIL" \
                --output "./references" \
                > "logs/${gene_name}_accession_collection.log" 2>&1

            if [ $? -eq 0 ]; then
                local ref_file="references/${gene_name}_accession_references.fasta"
                if [ -f "$ref_file" ]; then
                    local count
                    count=$(grep -c ">" "$ref_file" 2>/dev/null || echo 0)
                    print_success "Downloaded $count accession references for $gene_name"
                    total_collected=$((total_collected + count))
                else
                    print_warning "No accession references collected for $gene_name"
                fi
            else
                print_error "Failed to collect accession references from $gene"
            fi

        # Rest of function stays same but with similar fixes...
        elif [[ "$gene" =~ ^[A-Z]+[0-9]+\.[0-9]+ ]] || [[ "$gene" == *" "* ]]; then
            print_status "📋 Processing accession list for custom gene"

            gene_name=$(get_gene_name_interactively "$gene")
            gene_name=$(echo "$gene_name" | tr -d '[:cntrl:]' | tr -d '\033' | sed 's/\[[0-9;]*[mGK]//g')

            print_success "Using gene name: $gene_name"
            echo ""

            python3 accession_reference_downloader.py \
                --accessions "$gene" \
                --gene "$gene_name" \
                --email "$EMAIL" \
                --output "./references" \
                > "logs/${gene_name}_accession_collection.log" 2>&1

            if [ $? -eq 0 ]; then
                local ref_file="references/${gene_name}_accession_references.fasta"
                if [ -f "$ref_file" ]; then
                    local count
                    count=$(grep -c ">" "$ref_file" 2>/dev/null || echo 0)
                    print_success "Downloaded $count accession references for $gene_name"
                    total_collected=$((total_collected + count))
                else
                    print_warning "No accession references collected for $gene_name"
                fi
            else
                print_error "Failed to collect accession references"
            fi
        else
            # Regular gene name case
            print_status "📋 Processing gene name: $gene"
            read -p "🧬 Confirm gene name '$gene' or enter new name: " user_input

            if [ -n "$user_input" ]; then
                gene_name=$(echo "$user_input" | tr -d '[:cntrl:]' | tr -d '\033' | sed 's/\[[0-9;]*[mGK]//g')
            else
                gene_name="$gene"
            fi

            print_success "Using gene name: $gene_name"
            echo ""

            # Continue with taxonomic collection...
            print_status "Collecting $gene_name from groups: ${TAXONOMIC_GROUPS[*]}"

            python3 collect_gene_by_taxon.py \
                --gene "$gene_name" \
                --taxa "${TAXONOMIC_GROUPS[@]}" \
                --email "$EMAIL" \
                --output "./references" \
                --max_per_species 3 \
                --min_length 200 \
                > "logs/${gene_name}_reference_collection.log" 2>&1

            if [ $? -eq 0 ]; then
                local count=0
                for taxa in "${TAXONOMIC_GROUPS[@]}"; do
                    local ref_file="references/${gene_name}_${taxa}_references.fasta"
                    if [ -f "$ref_file" ]; then
                        local gene_count
                        gene_count=$(grep -c ">" "$ref_file" 2>/dev/null || echo 0)
                        count=$((count + gene_count))
                    fi
                done

                if [ "$count" -gt 0 ]; then
                    print_success "Collected $count $gene_name references"
                    total_collected=$((total_collected + count))

                    local combined_file="references/${gene_name}_all_references.fasta"
                    cat references/${gene_name}_*_references.fasta > "$combined_file" 2>/dev/null
                    print_success "Combined references: $combined_file"
                else
                    print_warning "No $gene_name references collected"
                fi
            else
                print_error "Failed to collect $gene_name references"
            fi
        fi

        echo "---"
    done

    print_success "Total reference sequences collected: $total_collected"
}

download_assemblies() {
    if [ -z "$INPUT_FILE" ] || [ ! -f "$INPUT_FILE" ]; then
        print_warning "No input file specified, skipping assembly download"
        return 0
    fi

    print_status "Downloading genome assemblies..."

    python3 datasets_assembly_download.py \
        --input "$INPUT_FILE" \
        --output "./assemblies" \
        --max_size 4000 \
        > logs/assembly_download.log 2>&1

    if [ $? -eq 0 ]; then
        local count
        count=$(find assemblies -name "*.fasta" | wc -l 2>/dev/null || echo 0)
        print_success "Downloaded $count assemblies"
    else
        print_error "Assembly download failed. Check logs/assembly_download.log"
        return 1
    fi
}

extract_genes() {
    print_status "Extracting genes from assemblies..."

    local assembly_count
    assembly_count=$(find assemblies -name "*.fasta" | wc -l 2>/dev/null || echo 0)
    if [ "$assembly_count" -eq 0 ]; then
        print_error "No assembly files found for gene extraction"
        return 1
    fi

    local total_extractions=0

    for gene in "${GENES[@]}"; do
        local ref_file=""
        local gene_name=""

        if [ -f "$gene" ]; then
            # For file inputs, we need to find the actual gene name used
            # Look for *_accession_references.fasta files
            local possible_files=(references/*_accession_references.fasta)
            if [ ${#possible_files[@]} -eq 1 ] && [ -f "${possible_files[0]}" ]; then
                ref_file="${possible_files[0]}"
                gene_name=$(basename "$ref_file" | sed 's/_accession_references.fasta$//')
            else
                # Multiple files or not found - ask user
                echo "Available reference files:"
                ls -la references/*_accession_references.fasta 2>/dev/null || echo "No accession reference files found"
                read -p "Enter gene name for extraction: " gene_name
                ref_file="references/${gene_name}_accession_references.fasta"
            fi
        elif [[ "$gene" =~ ^[A-Z]+[0-9]+\.[0-9]+ ]] || [[ "$gene" == *" "* ]]; then
            # Similar logic for accession lists
            local possible_files=(references/*_accession_references.fasta)
            if [ ${#possible_files[@]} -eq 1 ] && [ -f "${possible_files[0]}" ]; then
                ref_file="${possible_files[0]}"
                gene_name=$(basename "$ref_file" | sed 's/_accession_references.fasta$//')
            else
                read -p "Enter gene name for extraction: " gene_name
                ref_file="references/${gene_name}_accession_references.fasta"
            fi
        else
            # Regular gene case
            ref_file="references/${gene}_all_references.fasta"
            gene_name="$gene"
        fi

        if [ -f "$ref_file" ]; then
            print_status "Extracting $gene_name from $assembly_count assemblies..."

            ./extract_genes.sh \
                "$ref_file" \
                "./assemblies" \
                "./results/$gene_name" \
                "$gene_name" \
                > "logs/${gene_name}_extraction.log" 2>&1

            if [ $? -eq 0 ]; then
                local extracted_file="results/$gene_name/all_${gene_name}_extracted.fasta"
                if [ -f "$extracted_file" ]; then
                    local count
                    count=$(grep -c ">" "$extracted_file" 2>/dev/null || echo 0)
                    print_success "Extracted $count $gene_name sequences"
                    total_extractions=$((total_extractions + count))
                else
                    print_warning "No $gene_name sequences extracted"
                fi
            else
                print_error "Failed to extract $gene_name"
            fi
        else
            print_warning "No reference file found for $gene_name, skipping extraction"
        fi
    done

    print_success "Total gene sequences extracted: $total_extractions"
}
generate_summary() {
    print_status "Generating pipeline summary..."

    local summary_file="logs/pipeline_summary.txt"

    {
        echo "BLAST PIPELINE v1.0 SUMMARY"
        echo "=================================="
        echo "Run date: $(date)"
        echo "Email: $EMAIL"

        # FIXED: Only show taxonomic groups if actually used
        local used_taxonomic_collection=false
        for gene in "${GENES[@]}"; do
            # Check if this gene used taxonomic collection (not file-based)
            if [ ! -f "$gene" ] && [[ ! "$gene" =~ ^[A-Z]+[0-9]+\.[0-9]+ ]] && [[ "$gene" != *" "* ]]; then
                used_taxonomic_collection=true
                break
            fi
        done

        if [ "$used_taxonomic_collection" = true ]; then
            echo "Taxonomic groups: ${TAXONOMIC_GROUPS[*]}"
        else
            echo "Input method: Custom reference files/accessions"
        fi

        echo "Genes: ${GENES[*]}"
        echo ""

        echo "ASSEMBLIES:"
        local assembly_count
        assembly_count=$(find assemblies -name "*.fasta" 2>/dev/null | wc -l || echo 0)
        echo "Total downloaded: $assembly_count"
        if [ "$assembly_count" -gt 0 ]; then
            find assemblies -name "*.fasta" -exec ls -lh {} \; 2>/dev/null | head -5 | awk '{print "  " $9 " (" $5 ")"}'
            if [ "$assembly_count" -gt 5 ]; then
                echo "  ... and $((assembly_count - 5)) more"
            fi
        fi
        echo ""

        echo "REFERENCE GENES:"
        # FIXED: Accurate reference counting
        for gene in "${GENES[@]}"; do
            local total=0
            local gene_display=""
            local actual_gene_name=""

            if [ -f "$gene" ]; then
                # File path case - find the actual gene name used
                actual_gene_name=$(find references -name "*_accession_references.fasta" -exec basename {} \; | sed 's/_accession_references.fasta$//' | tail -1)

                if [ -z "$actual_gene_name" ]; then
                    # Fallback: try to extract from filename
                    actual_gene_name=$(basename "$gene" | sed 's/\.[^.]*$//')
                fi

                gene_display="$gene (extracted as $actual_gene_name)"

                # Count from accession references file
                local accession_ref_file="references/${actual_gene_name}_accession_references.fasta"
                if [ -f "$accession_ref_file" ]; then
                    total=$(grep -c ">" "$accession_ref_file" 2>/dev/null || echo 0)
                fi

            elif [[ "$gene" =~ ^[A-Z]+[0-9]+\.[0-9]+ ]] || [[ "$gene" == *" "* ]]; then
                # Accession list case
                gene_display="$gene (custom accessions)"

                # Find any accession references file
                local custom_ref_file=$(find references -name "*_accession_references.fasta" | head -1)
                if [ -f "$custom_ref_file" ]; then
                    total=$(grep -c ">" "$custom_ref_file" 2>/dev/null || echo 0)
                fi

            else
                # Regular gene case - taxonomic collection
                gene_display="$gene (taxonomic collection)"
                actual_gene_name="$gene"

                for taxa in "${TAXONOMIC_GROUPS[@]}"; do
                    local ref_file="references/${gene}_${taxa}_references.fasta"
                    if [ -f "$ref_file" ]; then
                        local count
                        count=$(grep -c ">" "$ref_file" 2>/dev/null || echo 0)
                        total=$((total + count))
                    fi
                done

                # Also check combined file
                local combined_file="references/${gene}_all_references.fasta"
                if [ -f "$combined_file" ]; then
                    total=$(grep -c ">" "$combined_file" 2>/dev/null || echo 0)
                fi
            fi

            echo "$gene_display: $total reference sequences"
        done
        echo ""

        echo "GENE EXTRACTION RESULTS:"
        # FIXED: Accurate extraction counting
        local total_extracted=0

        for gene in "${GENES[@]}"; do
            local gene_display=""
            local actual_gene_name=""
            local extracted_count=0

            if [ -f "$gene" ]; then
                # File path case - find actual gene name used
                actual_gene_name=$(find references -name "*_accession_references.fasta" -exec basename {} \; | sed 's/_accession_references.fasta$//' | tail -1)

                if [ -z "$actual_gene_name" ]; then
                    actual_gene_name=$(basename "$gene" | sed 's/\.[^.]*$//')
                fi

                gene_display="$gene (extracted as $actual_gene_name)"

            elif [[ "$gene" =~ ^[A-Z]+[0-9]+\.[0-9]+ ]] || [[ "$gene" == *" "* ]]; then
                # Accession list case - find the gene name that was used
                actual_gene_name=$(find references -name "*_accession_references.fasta" -exec basename {} \; | sed 's/_accession_references.fasta$//' | tail -1)

                if [ -z "$actual_gene_name" ]; then
                    actual_gene_name="custom_gene"
                fi

                gene_display="$gene (extracted as $actual_gene_name)"

            else
                # Regular gene case
                actual_gene_name="$gene"
                gene_display="$gene"
            fi

            # Count extracted sequences
            local extracted_file="results/$actual_gene_name/all_${actual_gene_name}_extracted.fasta"
            if [ -f "$extracted_file" ]; then
                extracted_count=$(grep -c ">" "$extracted_file" 2>/dev/null || echo 0)
            else
                # Try individual files in results directory
                if [ -d "results/$actual_gene_name" ]; then
                    extracted_count=$(find "results/$actual_gene_name" -name "*_${actual_gene_name}.fasta" -exec grep -c ">" {} \; 2>/dev/null | awk '{sum+=$1} END {print sum+0}')
                fi
            fi

            echo "$gene_display: $extracted_count sequences extracted"
            total_extracted=$((total_extracted + extracted_count))
        done

        echo ""
        echo "TOTAL EXTRACTED: $total_extracted sequences"

        # FIXED: Add extraction quality summary if available
        echo ""
        echo "EXTRACTION QUALITY:"
        for gene in "${GENES[@]}"; do
            local actual_gene_name=""

            # Determine actual gene name used
            if [ -f "$gene" ]; then
                actual_gene_name=$(find references -name "*_accession_references.fasta" -exec basename {} \; | sed 's/_accession_references.fasta$//' | tail -1)
                if [ -z "$actual_gene_name" ]; then
                    actual_gene_name=$(basename "$gene" | sed 's/\.[^.]*$//')
                fi
            elif [[ "$gene" =~ ^[A-Z]+[0-9]+\.[0-9]+ ]] || [[ "$gene" == *" "* ]]; then
                actual_gene_name=$(find references -name "*_accession_references.fasta" -exec basename {} \; | sed 's/_accession_references.fasta$//' | tail -1)
                if [ -z "$actual_gene_name" ]; then
                    actual_gene_name="custom_gene"
                fi
            else
                actual_gene_name="$gene"
            fi

            # Check extraction log for quality metrics
            local log_file="logs/${actual_gene_name}_extraction.log"
            if [ -f "$log_file" ]; then
                local success_entries
                success_entries=$(grep "Success" "$log_file" | wc -l 2>/dev/null || echo 0)
                if [ "$success_entries" -gt 0 ]; then
                    echo "$actual_gene_name:"
                    grep "Success" "$log_file" | head -3 | sed 's/^/  /'
                    if [ "$success_entries" -gt 3 ]; then
                        echo "  ... and $((success_entries - 3)) more successful extractions"
                    fi
                fi
            fi
        done

    } > "$summary_file"

    cat "$summary_file"
    print_success "Summary saved: $summary_file"
}


main() {
    # Parse arguments FIRST - before any validation
    while [[ $# -gt 0 ]]; do
        case $1 in
            --email)
                EMAIL="$2"
                shift 2
                ;;
            --species)
                INPUT_FILE="$2"
                shift 2
                ;;
            --input)
                INPUT_FILE="$2"
                shift 2
                ;;
            --genes)
                IFS=',' read -r -a GENES <<< "$2"
                shift 2
                ;;
            --taxa)
                IFS=',' read -r -a TAXONOMIC_GROUPS <<< "$2"
                shift 2
                ;;
            --list-taxa)
                # Set dummy email for list-taxa command
                if [ -z "$EMAIL" ]; then
                    EMAIL="dummy@example.com"
                fi
                list_taxonomic_groups
                exit 0
                ;;
            --references-only)
                REFERENCES_ONLY=1
                shift
                ;;
            --assemblies-only)
                ASSEMBLIES_ONLY=1
                shift
                ;;
            --extract-only)
                EXTRACT_ONLY=1
                shift
                ;;
            --help)
                show_help
                exit 0
                ;;
            *)
                print_error "Unknown option: $1"
                echo "Use --help for usage information"
                exit 1
                ;;
        esac
    done

    # NOW validate email - AFTER parsing
    if [ -z "$EMAIL" ]; then
        print_error "Email is required for NCBI API access"
        echo "Use: $0 --email your@email.com"
        exit 1
    fi

    # Show startup messages
    print_status "🧬 PHYLOGENETIC PIPELINE v2.0 STARTING"
    print_status "Genes: ${GENES[*]}"

    # FIXED: Check gene types before showing taxonomic groups
    local uses_taxonomic_collection=false
    for gene in "${GENES[@]}"; do
        # Check if this gene uses taxonomic collection (not file-based)
        if [ ! -f "$gene" ] && [[ ! "$gene" =~ ^[A-Z]+[0-9]+\.[0-9]+ ]] && [[ "$gene" != *" "* ]]; then
            uses_taxonomic_collection=true
            break
        fi
    done

    if [ "$uses_taxonomic_collection" = true ]; then
        print_status "Taxonomic groups: ${TAXONOMIC_GROUPS[*]}"
    else
        print_status "Input method: Custom reference files/accessions"
    fi

    # Check dependencies and setup
    check_dependencies
    setup_directories

    # Log start time
    local start_time
    start_time=$(date)
    echo "Started at: $start_time" > logs/pipeline.log

    # Run pipeline steps based on options
    if [ "$ASSEMBLIES_ONLY" == "1" ]; then
        print_status "Running assemblies-only mode"
        download_assemblies

    elif [ "$REFERENCES_ONLY" == "1" ]; then
        print_status "Running references-only mode"
        collect_reference_genes

    elif [ "$EXTRACT_ONLY" == "1" ]; then
        print_status "Running extract-only mode"
        extract_genes
        generate_summary

    else
        # Full pipeline
        print_status "Running full pipeline"
        collect_reference_genes

        if [ -n "$INPUT_FILE" ] && [ -f "$INPUT_FILE" ]; then
            download_assemblies
            extract_genes
        else
            print_warning "No input file provided or file not found, skipping assembly download and extraction"
            print_warning "Use --input accessions.txt or --input species_list.tsv"
        fi

        generate_summary
    fi

    # Log completion time
    local end_time
    end_time=$(date)
    echo "Completed at: $end_time" >> logs/pipeline.log

    # Final success message
    print_success "🎉 Pipeline completed!"
    echo ""
    print_status "📂 Check results in: ./results/"
    print_status "📊 Summary: ./logs/pipeline_summary.txt"
}

main "$@"