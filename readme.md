# BLAST Pipeline

A comprehensive bioinformatics pipeline for automated gene extraction from genome assemblies using taxonomically-filtered reference sequences or custom accession lists.

## Overview

This pipeline automates the process of:
1. **Reference Collection**: Gathering high-quality gene references from taxonomic groups or accession lists
2. **Assembly Download**: Retrieving genome assemblies from NCBI
3. **Gene Extraction**: BLAST-based extraction of target genes from assemblies
4. **Organism Annotation**: Automatic species identification and metadata enrichment

## Features

- **Flexible Reference Sources**: Use taxonomic groups, custom accession lists, or file paths
- **Taxonomic Filtering**: Collect references from specific clades (mammals, insects, etc.)
- **Custom Accession Support**: Use manually curated reference sequences
- **Automatic Assembly Download**: NCBI datasets integration for genome retrieval
- **Quality Control**: Identity and coverage thresholds for reliable extractions
- **Organism Annotation**: Automatic species name lookup and file labeling
- **Comprehensive Logging**: Detailed logs for troubleshooting and analysis
- **Resume Capability**: Skip previously processed files

## Installation

### Prerequisites

**System Dependencies:**
```bash
# BLAST+ tools
sudo apt-get install ncbi-blast+

# Samtools for sequence extraction
sudo apt-get install samtools

# Seqtk for sequence manipulation
sudo apt-get install seqtk

# Calculator for quality checks
sudo apt-get install bc

# NCBI datasets CLI
conda install -c conda-forge ncbi-datasets-cli
```

**Python Dependencies:**
```bash
pip install biopython requests
```

### Setup

1. Clone or download the pipeline files
2. Make shell scripts executable:
```bash
chmod +x master_pipeline.sh extract_genes.sh
```
3. Test dependencies:
```bash
./master_pipeline.sh --help
```

## Usage

### Quick Start

```bash
# Create species list (if using assemblies)
cat > species_list.tsv << 'EOF'
Homo	sapiens	human
Mus	musculus	mouse
EOF

# Run full pipeline with taxonomic gene collection
./master_pipeline.sh \
    --email your@email.com \
    --input species_list.tsv \
    --genes 18S \
    --taxa mammals

# Or use custom accession references
./master_pipeline.sh \
    --email your@email.com \
    --input assemblies.txt \
    --genes accession_refs.txt
```

### Input File Formats

**Species List (TSV):**
```
Genus	species	tag
Homo	sapiens	human
Mus	musculus	mouse
```

**Assembly Accessions:**
```
GCA_034510155.1
GCF_000001405.40
GCA_050613505.1
```

**Reference Accessions:**
```
MT536114.1
KM853220.1
AF423791.1
```

### Advanced Options

**Taxonomic Gene Collection:**
```bash
# Use predefined taxonomic groups
./master_pipeline.sh \
    --email your@email.com \
    --genes CLOCK,PER1 \
    --taxa mammals,birds

# Use TaxIDs directly
./master_pipeline.sh \
    --email your@email.com \
    --genes 18S \
    --taxa 6993,40674

# Mix taxonomic groups and TaxIDs
./master_pipeline.sh \
    --email your@email.com \
    --genes rRNA \
    --taxa mammals,6993,Orthoptera
```

**Step-by-step Execution:**
```bash
# Only collect references
./master_pipeline.sh \
    --email your@email.com \
    --genes 18S \
    --taxa mammals \
    --references-only

# Only download assemblies
./master_pipeline.sh \
    --email your@email.com \
    --input assemblies.txt \
    --assemblies-only

# Only extract genes (requires existing references and assemblies)
./master_pipeline.sh \
    --email your@email.com \
    --extract-only
```

## Available Taxonomic Groups

- **mammals** - Mammalia (40,674 species)
- **birds** - Aves  
- **reptiles** - Reptilia
- **amphibians** - Amphibia
- **fish** - Actinopterygii (ray-finned fishes)
- **arthropods** - Arthropoda
- **insects** - Insecta
- **nematodes** - Nematoda
- **plants** - Viridiplantae
- **fungi** - Fungi
- **bacteria** - Bacteria
- **archaea** - Archaea
- **vertebrates** - Vertebrata
- **primates** - Primates
- **rodents** - Rodentia

View complete list:
```bash
./master_pipeline.sh --list-taxa
```

## Output Structure

```
blast_pipeline/
├── assemblies/           # Downloaded genome assemblies
├── references/           # Reference gene sequences
│   ├── combined/        # Combined references per gene
│   └── individual files # Species-specific references
├── results/             # Extracted gene sequences
│   └── GENE_NAME/       # Per-gene results
│       ├── individual/  # Per-assembly extractions
│       └── all_GENE_extracted.fasta  # Combined results
└── logs/               # Execution logs
    ├── pipeline_summary.txt
    └── detailed logs
```

## File Naming Convention

**Assemblies:**
- `Species_name_tag_GCA123456.1.fasta`

**Extracted Genes:**
- `Species_name_GCA123456.1_GENE.fasta`

**FASTA Headers:**
- `>Species_name_GCA123456.1_GENE`

## Quality Control

**Reference Collection:**
- Prioritizes RefSeq over GenBank sequences
- Filters by sequence length (default: 100-50,000 bp)
- Excludes partial, predicted, and low-quality sequences
- Limits sequences per species (default: 3)

**Gene Extraction:**
- Minimum identity threshold: 70%
- Minimum query coverage: 70%
- Minimum sequence length: 50 bp
- Automatic reverse complement handling

**Assembly Selection:**
- Prioritizes RefSeq over GenBank
- Prefers Complete Genome > Chromosome > Scaffold
- Quality scoring based on N50 values
- Size limits (default: 4GB maximum)

## Examples

### Example 1: 18S rRNA from Multiple Taxa
```bash
# Collect 18S references from vertebrates and extract from insect assemblies
./master_pipeline.sh \
    --email user@example.com \
    --genes 18S \
    --taxa vertebrates,insects \
    --input insect_assemblies.txt
```

### Example 2: Custom Gene References
```bash
# Use manually curated gene references
echo -e "MT536114.1\nKM853220.1\nAF423791.1" > my_refs.txt

./master_pipeline.sh \
    --email user@example.com \
    --genes my_refs.txt \
    --input target_assemblies.txt
```

### Example 3: Multiple Genes
```bash
# Extract multiple circadian rhythm genes
./master_pipeline.sh \
    --email user@example.com \
    --genes CLOCK,ARNTL,PER1,PER2 \
    --taxa mammals,birds \
    --input vertebrate_assemblies.txt
```

## Component Scripts

**Individual scripts can be used standalone:**

```bash
# Collect taxonomic references
python3 collect_gene_by_taxon.py \
    --gene 18S \
    --taxa mammals \
    --email user@example.com

# Download assemblies
python3 datasets_assembly_download.py \
    --input assemblies.txt \
    --output ./assemblies \
    --email user@example.com

# Download from accession list
python3 accession_reference_downloader.py \
    --accessions refs.txt \
    --gene 18S \
    --email user@example.com

# Extract genes
./extract_genes.sh \
    references.fasta \
    ./assemblies/ \
    ./results/ \
    18S \
    user@example.com
```

## Troubleshooting

**Common Issues:**

1. **BLAST tools not found**
   ```bash
   conda install -c bioconda blast
   ```

2. **Python import errors**
   ```bash
   pip install biopython requests
   ```

3. **No sequences extracted**
   - Check reference quality
   - Verify assembly integrity
   - Review logs in `logs/` directory

4. **NCBI API errors**
   - Ensure valid email address
   - Check internet connection
   - Review rate limiting in logs

**Log Files:**
Always check relevant logs for detailed error information:
- `logs/pipeline_summary.txt` - Overall results
- `logs/GENE_extraction.log` - Gene-specific extraction
- `logs/assembly_download.log` - Assembly download status
- `logs/reference_collection.log` - Reference collection status

## Performance Notes

- Assembly downloads: ~1-5 minutes per genome
- Reference collection: ~30 seconds per gene per taxonomic group
- Gene extraction: ~1-2 minutes per assembly per gene (depends on size)
- Memory usage: Primarily limited by BLAST database size

## Citation

If you use this pipeline in your research, please cite the relevant databases:

- NCBI Nucleotide Database
- NCBI Assembly Database
- BLAST+ Tools
- NCBI Datasets

## License

This pipeline is provided for research purposes. Please respect NCBI usage policies and rate limits.

## Contributing

Report issues and suggestions through the project repository. Contributions welcome for:
- Additional taxonomic groups
- New gene extraction methods
- Performance optimizations
- Documentation improvements
