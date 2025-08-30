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
    print_status "Collecting reference genes..."

    local total_collected=0

    for gene in "${GENES[@]}"; do
        # Check if gene is a file path
        if [ -f "$gene" ]; then
            print_status "Processing accession file: $gene"

            # Extract gene name from filename (remove path and extension)
            local gene_name
            gene_name=$(basename "$gene" | sed 's/\.[^.]*$//' | sed 's/.*_//')
            if [[ "$gene_name" == *"accession"* ]] || [[ "$gene_name" == *"ref"* ]]; then
                gene_name="unknown"  # Default fallback
            fi

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

        # Check if gene looks like accession list (contains accession patterns or spaces)
        elif [[ "$gene" =~ ^[A-Z]+[0-9]+\.[0-9]+ ]] || [[ "$gene" == *" "* ]]; then
            print_status "Processing accession list for custom gene"

            # Treat as accession list
            python3 accession_reference_downloader.py \
                --accessions "$gene" \
                --gene "custom_gene" \
                --email "$EMAIL" \
                --output "./references" \
                > "logs/custom_accession_collection.log" 2>&1

            if [ $? -eq 0 ]; then
                local ref_file="references/custom_gene_accession_references.fasta"
                if [ -f "$ref_file" ]; then
                    local count
                    count=$(grep -c ">" "$ref_file" 2>/dev/null || echo 0)
                    print_success "Downloaded $count accession references"
                    total_collected=$((total_collected + count))
                else
                    print_warning "No accession references collected"
                fi
            else
                print_error "Failed to collect accession references"
            fi
        else
            # Regular gene name - collect from taxonomic groups
            print_status "Collecting $gene from groups: ${TAXONOMIC_GROUPS[*]}"

            python3 collect_gene_by_taxon.py \
                --gene "$gene" \
                --taxa "${TAXONOMIC_GROUPS[@]}" \
                --email "$EMAIL" \
                --output "./references" \
                --max_per_species 3 \
                --min_length 200 \
                > "logs/${gene}_reference_collection.log" 2>&1

            if [ $? -eq 0 ]; then
                # Count collected sequences
                local count=0
                for taxa in "${TAXONOMIC_GROUPS[@]}"; do
                    local ref_file="references/${gene}_${taxa}_references.fasta"
                    if [ -f "$ref_file" ]; then
                        local gene_count
                        gene_count=$(grep -c ">" "$ref_file" 2>/dev/null || echo 0)
                        count=$((count + gene_count))
                    fi
                done

                if [ "$count" -gt 0 ]; then
                    print_success "Collected $count $gene references"
                    total_collected=$((total_collected + count))

                    # Combine all taxonomic groups for this gene
                    local combined_file="references/${gene}_all_references.fasta"
                    cat references/${gene}_*_references.fasta > "$combined_file" 2>/dev/null
                    print_success "Combined references: $combined_file"
                else
                    print_warning "No $gene references collected"
                fi
            else
                print_error "Failed to collect $gene references"
            fi
        fi
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
        # Determine reference file based on gene type
        local ref_file=""
        local gene_name=""

        if [ -f "$gene" ]; then
            # File path case
            gene_name=$(basename "$gene" | sed 's/\.[^.]*$//' | sed 's/.*_//')
            if [[ "$gene_name" == *"accession"* ]] || [[ "$gene_name" == *"ref"* ]]; then
                gene_name="18S"  # Default fallback
            fi
            ref_file="references/${gene_name}_accession_references.fasta"
        elif [[ "$gene" =~ ^[A-Z]+[0-9]+\.[0-9]+ ]] || [[ "$gene" == *" "* ]]; then
            # Accession list case
            ref_file="references/custom_gene_accession_references.fasta"
            gene_name="custom_gene"
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
            print_warning "No reference file found for $gene, skipping extraction"
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
        echo "Taxonomic groups: ${TAXONOMIC_GROUPS[*]}"
        echo "Genes: ${GENES[*]}"
        echo ""

        echo "ASSEMBLIES:"
        local assembly_count
        assembly_count=$(find assemblies -name "*.fasta" | wc -l 2>/dev/null || echo 0)
        echo "Total downloaded: $assembly_count"
        if [ "$assembly_count" -gt 0 ]; then
            find assemblies -name "*.fasta" -exec ls -lh {} \; | head -5 | awk '{print "  " $9 " (" $5 ")"}'
            if [ "$assembly_count" -gt 5 ]; then
                echo "  ... and $((assembly_count - 5)) more"
            fi
        fi
        echo ""

        echo "REFERENCE GENES:"
        for gene in "${GENES[@]}"; do
            local total=0
            local gene_display=""

            if [ -f "$gene" ]; then
                # File path case - extract gene name
                local gene_name
                gene_name=$(basename "$gene" | sed 's/\.[^.]*$//' | sed 's/.*_//')
                if [[ "$gene_name" == *"accession"* ]] || [[ "$gene_name" == *"ref"* ]]; then
                    gene_name="18S"
                fi
                gene_display="$gene (extracted as $gene_name)"

                # Count from accession references file
                local accession_ref_file="references/${gene_name}_accession_references.fasta"
                if [ -f "$accession_ref_file" ]; then
                    total=$(grep -c ">" "$accession_ref_file" 2>/dev/null || echo 0)
                fi
            elif [[ "$gene" =~ ^[A-Z]+[0-9]+\.[0-9]+ ]] || [[ "$gene" == *" "* ]]; then
                # Accession list case
                gene_display="$gene (custom accessions)"
                local custom_ref_file="references/custom_gene_accession_references.fasta"
                if [ -f "$custom_ref_file" ]; then
                    total=$(grep -c ">" "$custom_ref_file" 2>/dev/null || echo 0)
                fi
            else
                # Regular gene case
                gene_display="$gene"
                for taxa in "${TAXONOMIC_GROUPS[@]}"; do
                    local ref_file="references/${gene}_${taxa}_references.fasta"
                    if [ -f "$ref_file" ]; then
                        local count
                        count=$(grep -c ">" "$ref_file" 2>/dev/null || echo 0)
                        total=$((total + count))
                    fi
                done
            fi

            echo "$gene_display: $total reference sequences"
        done
        echo ""

        echo "GENE EXTRACTION RESULTS:"
        local total_extracted=0
        for gene in "${GENES[@]}"; do
            local gene_name=""
            local gene_display=""

            if [ -f "$gene" ]; then
                # File path case
                gene_name=$(basename "$gene" | sed 's/\.[^.]*$//' | sed 's/.*_//')
                if [[ "$gene_name" == *"accession"* ]] || [[ "$gene_name" == *"ref"* ]]; then
                    gene_name="18S"
                fi
                gene_display="$gene (extracted as $gene_name)"
            elif [[ "$gene" =~ ^[A-Z]+[0-9]+\.[0-9]+ ]] || [[ "$gene" == *" "* ]]; then
                # Accession list case
                gene_name="custom_gene"
                gene_display="$gene (custom accessions)"
            else
                # Regular gene case
                gene_name="$gene"
                gene_display="$gene"
            fi

            local extracted_file="results/$gene_name/all_${gene_name}_extracted.fasta"
            if [ -f "$extracted_file" ]; then
                local count
                count=$(grep -c ">" "$extracted_file" 2>/dev/null || echo 0)
                echo "$gene_display: $count sequences extracted"
                total_extracted=$((total_extracted + count))
            else
                echo "$gene_display: 0 sequences extracted"
            fi
        done
        echo ""
        echo "TOTAL EXTRACTED: $total_extracted sequences"

    } > "$summary_file"

    cat "$summary_file"
    print_success "Summary saved: $summary_file"
}

main() {
    # Parse arguments
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
                get_email() { EMAIL="dummy@example.com"; }
                get_email
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

    # Validate email
    if [ -z "$EMAIL" ]; then
        print_error "Email is required for NCBI API access"
        echo "Use: $0 --email your@email.com"
        exit 1
    fi

    print_status "🧬 PHYLOGENETIC PIPELINE v2.0 STARTING"
    print_status "Genes: ${GENES[*]}"
    print_status "Taxonomic groups: ${TAXONOMIC_GROUPS[*]}"

    check_dependencies
    setup_directories

    local start_time
    start_time=$(date)
    echo "Started at: $start_time" > logs/pipeline.log

    # Run pipeline steps
    if [ "$ASSEMBLIES_ONLY" == "1" ]; then
        download_assemblies

    elif [ "$REFERENCES_ONLY" == "1" ]; then
        collect_reference_genes

    elif [ "$EXTRACT_ONLY" == "1" ]; then
        extract_genes
        generate_summary

    else
        # Full pipeline
        collect_reference_genes

        if [ -n "$INPUT_FILE" ]; then
            download_assemblies
            extract_genes
        else
            print_warning "No input file provided, skipping assembly download and extraction"
            print_warning "Use --input accessions.txt or --input species_list.tsv"
        fi

        generate_summary
    fi

    local end_time
    end_time=$(date)
    echo "Completed at: $end_time" >> logs/pipeline.log

    print_success "🎉 Pipeline completed!"
    echo ""
    print_status "📂 Check results in: ./results/"
    print_status "📊 Summary: ./logs/pipeline_summary.txt"
}

# Run main function
main "$@"