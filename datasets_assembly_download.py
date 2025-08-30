#!/usr/bin/env python3

import argparse
import os
import sys
import time
import json
import subprocess
import shutil
import zipfile
from pathlib import Path


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


def download_by_accession(accession, output_dir, tag, max_size_mb=5000):
    """Download assembly by direct accession"""
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

        # Move to final location
        final_filename = f"{tag}_{accession}.fasta"
        final_path = Path(output_dir) / final_filename

        shutil.move(str(genomic_file), str(final_path))

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


def search_and_download_species(species, tag, output_dir, max_size_mb=5000):
    """Search for species and download best assembly"""
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
                print(f"❌ No assemblies found for {species}")
                return False

            assemblies = data['reports']

            # Filter and rank assemblies
            valid_assemblies = []
            for assembly in assemblies:
                try:
                    accession = assembly.get('accession', '')
                    assembly_info = assembly.get('assembly_info', {})
                    assembly_stats = assembly.get('assembly_stats', {})

                    total_length = int(assembly_stats.get('total_sequence_length', 0))
                    size_mb = total_length / (1024 * 1024) if total_length else 0

                    if size_mb <= max_size_mb:
                        valid_assemblies.append({
                            'accession': accession,
                            'name': assembly_info.get('assembly_name', ''),
                            'level': assembly_info.get('assembly_level', 'Unknown'),
                            'size_mb': size_mb,
                            'is_refseq': accession.startswith('GCF_'),
                            'status': assembly_info.get('assembly_status', 'unknown')
                        })
                except:
                    continue

            if not valid_assemblies:
                print(f"❌ No suitable assemblies for {species}")
                return False

            # Sort by quality
            def ranking_key(asm):
                refseq_score = 1 if asm['is_refseq'] else 0
                status_score = 1 if asm['status'] == 'current' else 0
                level_score = {'Complete Genome': 4, 'Chromosome': 3, 'Scaffold': 2, 'Contig': 1}.get(asm['level'], 0)
                return (refseq_score, status_score, level_score)

            valid_assemblies.sort(key=ranking_key, reverse=True)

            # Download best assembly
            best = valid_assemblies[0]
            print(f"🎯 Selected: {best['accession']} ({best['level']}, {best['size_mb']:.1f} MB)")

            return download_by_accession(best['accession'], output_dir, tag, max_size_mb)

        except json.JSONDecodeError as e:
            print(f"❌ Failed to parse search results for {species}")
            return False

    except Exception as e:
        print(f"❌ Error processing {species}: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(description="Download genome assemblies - flexible input support")
    parser.add_argument("--input", required=True,
                        help="Input file: TSV with species (Genus\\tspecies\\ttag) or accessions list")
    parser.add_argument("--output", required=True,
                        help="Output directory for assemblies")
    parser.add_argument("--max_size", type=int, default=5000,
                        help="Maximum assembly size in MB (default: 5000)")
    parser.add_argument("--input_type", choices=["auto", "species", "accessions"],
                        default="auto", help="Input file type (default: auto-detect)")
    parser.add_argument("--dry_run", action="store_true",
                        help="Show what would be downloaded without downloading")

    args = parser.parse_args()

    # Check datasets CLI
    if not check_datasets_installed():
        sys.exit(1)

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Detect input type
    if args.input_type == "auto":
        input_type = detect_input_type(args.input)
        if input_type == 'auto':
            print("⚠️ Could not auto-detect input format. Assuming species format.")
            input_type = 'species'
        print(f"🔍 Detected input type: {input_type}")
    else:
        input_type = args.input_type

    # Setup logging
    log_file = output_dir / "download_log.tsv"
    failed_log = output_dir / "failed_downloads.tsv"

    # Process input file
    if input_type == 'accessions':
        items = read_accession_file(args.input)
        print(f"📋 Processing {len(items)} accessions")
    else:
        items = read_species_file(args.input)
        print(f"📋 Processing {len(items)} species")

    # Track downloads
    downloaded_species = set()
    if log_file.exists():
        with open(log_file, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith("Input"):
                    parts = line.strip().split('\t')
                    if parts:
                        downloaded_species.add(parts[0])

    # Initialize log files
    if not log_file.exists():
        with open(log_file, 'w') as f:
            f.write("Input\tAccession\tAssemblyName\tLevel\tSize_MB\tFilename\tInputType\n")

    if not failed_log.exists():
        with open(failed_log, 'w') as f:
            f.write("Input\tReason\tInputType\n")

    print(f"🚀 Starting downloads with max size: {args.max_size} MB")
    print(f"📊 Input type: {input_type}")
    print()

    success_count = 0
    failed_count = 0

    for item, tag in items:
        item_key = f"{tag}_{item}" if input_type == 'accessions' else item.replace(" ", "_")

        # Check if already downloaded
        if item_key in downloaded_species:
            print(f"⏭️ Already downloaded: {item}")
            success_count += 1
            continue

        print(f"\n🔍 Processing: {item}")

        if args.dry_run:
            print(f"🧪 [DRY RUN] Would process {item}")
            continue

        try:
            if input_type == 'accessions':
                filename = download_by_accession(item, args.output, tag, args.max_size)
            else:
                filename = search_and_download_species(item, tag, args.output, args.max_size)

            if filename:
                # Log successful download
                with open(log_file, 'a') as f:
                    f.write(f"{item}\t{item}\t-\t-\t-\t{filename}\t{input_type}\n")
                success_count += 1
                print(f"🎉 Success: {filename}")
            else:
                with open(failed_log, 'a') as f:
                    f.write(f"{item}\tDownload failed\t{input_type}\n")
                failed_count += 1

        except Exception as e:
            print(f"❌ Error processing {item}: {e}")
            with open(failed_log, 'a') as f:
                f.write(f"{item}\tProcessing error: {e}\t{input_type}\n")
            failed_count += 1

        # Rate limiting
        time.sleep(1)

    print(f"\n🎉 Download summary:")
    print(f"✅ Successful: {success_count}")
    print(f"❌ Failed: {failed_count}")
    print(f"📂 Files saved in: {args.output}")


if __name__ == "__main__":
    main()