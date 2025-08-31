#!/usr/bin/env python3

import argparse
import os
import sys
import time
import json
import subprocess
import shutil
import zipfile
import re
from pathlib import Path
from Bio import SeqIO


def clean_species_name(species_name):
    """Clean species name for filename usage"""
    if not species_name:
        return "Unknown_species"
    
    # Remove common prefixes and clean
    cleaned = species_name.replace("PREDICTED: ", "")
    cleaned = re.sub(r'\[.*?\]', '', cleaned)  # Remove [brackets]
    cleaned = re.sub(r'\(.*?\)', '', cleaned)  # Remove (parentheses)
    cleaned = cleaned.strip()
    
    # Extract genus + species (first two words)
    words = cleaned.split()
    if len(words) >= 2:
        genus_species = f"{words[0]}_{words[1]}"
    elif len(words) == 1:
        genus_species = words[0]
    else:
        genus_species = "Unknown_species"
    
    # Clean special characters
    genus_species = re.sub(r'[^\w\-_]', '_', genus_species)
    genus_species = re.sub(r'_+', '_', genus_species)
    
    return genus_species


def detect_input_type(input_file):
    """Detect if input file contains accessions or species names"""
    with open(input_file, 'r') as f:
        first_line = f.readline().strip()

    # Check if it looks like an accession
    if first_line.startswith(('GCF_', 'GCA_')):
        return 'accessions'

    # Check if it has tab-separated fields (species format)
    if '\t' in first_line:
        parts = first_line.split('\t')
        if len(parts) >= 2:
            # Check if first part looks like genus name (starts with capital)
            if parts[0].istitle() and parts[1].islower():
                return 'species'

    # Default fallback - try to determine from content
    return 'auto'


def read_accession_file(accession_file):
    """Read accession list from file"""
    accessions = []
    with open(accession_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            # Handle different formats
            if '\t' in line:
                # Tab-separated: accession\ttag or accession\tspecies_name\ttag
                parts = line.split('\t')
                accession = parts[0]
                tag = parts[-1] if len(parts) > 1 else accession
            else:
                # Simple list: one accession per line
                accession = line
                tag = accession

            accessions.append((accession, tag))

    return accessions


def read_species_file(species_file):
    """Read species list from TSV file"""
    species_list = []
    with open(species_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split('\t')
            if len(parts) >= 2:
                species = f"{parts[0]} {parts[1]}"
                tag = parts[2] if len(parts) > 2 else parts[1].lower()
                species_list.append((species, tag))

    return species_list


def check_datasets_installed():
    """Check if NCBI datasets CLI is installed"""
    try:
        result = subprocess.run(['datasets', '--version'], capture_output=True, text=True)
        if result.returncode == 0:
            print(f"✅ NCBI datasets found: {result.stdout.strip()}")
            return True
    except FileNotFoundError:
        pass

    print("❌ NCBI datasets CLI not found!")
    print("📦 Install with: conda install -c conda-forge ncbi-datasets-cli")
    return False


def enhance_assembly_headers(fasta_file, species_name, accession):
    """Enhance FASTA headers in assembly file with species name"""
    try:
        records = []
        with open(fasta_file, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                # Create enhanced header: >Species_name_AccessionVersion_scaffold/chromosome_ID
                scaffold_id = record.id
                new_id = f"{species_name}_{accession}_{scaffold_id}"
                record.id = new_id
                record.description = f"{species_name.replace('_', ' ')} {accession} {scaffold_id}"
                records.append(record)
        
        # Write back with enhanced headers
        with open(fasta_file, 'w') as f:
            SeqIO.write(records, f, "fasta")
        
        return True
    except Exception as e:
        print(f"⚠️ Could not enhance headers: {e}")
        return False


def download_by_accession(accession, output_dir, tag, max_size_mb=5000):
    """Download assembly by direct accession with species name enhancement"""
    print(f"📥 Downloading {accession} ({tag})...")

    temp_dir = Path(output_dir) / "temp_download"
    temp_dir.mkdir(exist_ok=True)

    try:
        # Download using datasets
        cmd = [
            'datasets', 'download', 'genome', 'accession', accession,
            '--filename', str(temp_dir / f"{accession}.zip"),
            '--include', 'genome'
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

        if result.returncode != 0:
            print(f"❌ Download failed for {accession}: {result.stderr}")
            return False

        zip_file = temp_dir / f"{accession}.zip"
        if not zip_file.exists():
            print(f"❌ Downloaded file not found for {accession}")
            return False

        # Extract the zip file
        with zipfile.ZipFile(zip_file, 'r') as zip_ref:
            zip_ref.extractall(temp_dir)

        # Find genomic.fna file
        genomic_files = list(temp_dir.glob("ncbi_dataset/data/*/*.fna"))

        if not genomic_files:
            print(f"❌ No genomic.fna file found for {accession}")
            return False

        genomic_file = genomic_files[0]

        # Check size
        size_mb = genomic_file.stat().st_size / (1024 * 1024)
        if size_mb > max_size_mb:
            print(f"⚠️ Assembly too large: {size_mb:.1f} MB > {max_size_mb} MB")
            return False

        # Get species name from assembly metadata
        species_name = get_species_from_assembly_metadata(temp_dir, accession)
        
        # Move to final location with enhanced filename
        final_filename = f"{species_name}_{accession}.fasta"
        final_path = Path(output_dir) / final_filename

        shutil.copy(str(genomic_file), str(final_path))

        # Enhance FASTA headers
        enhance_assembly_headers(str(final_path), species_name, accession)

        print(f"✅ Downloaded: {final_filename} ({size_mb:.1f} MB)")

        # Cleanup
        shutil.rmtree(temp_dir)
        return final_filename

    except subprocess.TimeoutExpired:
        print(f"⏰ Download timeout for {accession}")
        return False
    except Exception as e:
        print(f"❌ Error downloading {accession}: {e}")
        return False
    finally:
        if temp_dir.exists():
            shutil.rmtree(temp_dir, ignore_errors=True)


def get_species_from_assembly_metadata(temp_dir, accession):
    """Extract species name from downloaded assembly metadata"""
    try:
        # Look for assembly_data_report.jsonl
        report_files = list(temp_dir.glob("ncbi_dataset/data/assembly_data_report.jsonl"))
        
        if report_files:
            with open(report_files[0], 'r') as f:
                for line in f:
                    try:
                        data = json.loads(line.strip())
                        if data.get('accession') == accession:
                            organism = data.get('organism', {}).get('organismName', '')
                            if organism:
                                return clean_species_name(organism)
                    except json.JSONDecodeError:
                        continue
        
        # Fallback: try to extract from directory structure
        data_dirs = list(temp_dir.glob("ncbi_dataset/data/*"))
        if data_dirs:
            dir_name = data_dirs[0].name
            if dir_name != accession:
                # Sometimes directory name contains species info
                return clean_species_name(dir_name)
        
    except Exception as e:
        print(f"⚠️ Could not extract species from metadata: {e}")
    
    # Final fallback
    return f"Unknown_{accession.replace('.', '_')}"


def search_and_download_species(species, tag, output_dir, max_size_mb=5000):
    """Search for species and download best assembly with enhanced naming"""
    print(f"🔍 Searching assemblies for: {species}")

    try:
        # Search for assemblies
        cmd = [
            'datasets', 'summary', 'genome', 'taxon', species,
            '--limit', '10'
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)

        if result.returncode != 0:
            print(f"❌ Search failed for {species}: {result.stderr}")
            return False

        # Parse results
        try:
            data = json.loads(result.stdout)
            if 'reports' not in data or not data['reports']:
                print(f"❌ No assembl
