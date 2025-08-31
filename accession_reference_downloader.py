#!/usr/bin/env python3

import argparse
import os
import time
from Bio import Entrez, SeqIO
import re

class AccessionReferenceDownloader:
    def __init__(self, email, output_dir="references"):
        self.email = email
        self.output_dir = output_dir
        Entrez.email = email
        os.makedirs(output_dir, exist_ok=True)
        
        # Species name cache
        self.species_cache = {}
        
    def log(self, message):
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{timestamp}] {message}")
        
    def read_accession_list(self, input_string):
        """Parse accession list from string or file"""
        accessions = []
        
        # Check if it's a file
        if os.path.isfile(input_string):
            with open(input_string, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        accessions.append(line)
        else:
            # Treat as space/comma separated list
            accessions = input_string.replace(',', ' ').split()
        
        return accessions
    
    def clean_species_name(self, species_name):
        """Clean species name for filename/header usage"""
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
    
    def get_species_from_accession(self, accession):
        """Get species name from GenBank/RefSeq accession"""
        # Check cache first
        if accession in self.species_cache:
            return self.species_cache[accession]
        
        try:
            # Search for the accession
            handle = Entrez.esearch(db="nucleotide", term=f"{accession}[Accession]")
            search_results = Entrez.read(handle)
            handle.close()
            
            if not search_results["IdList"]:
                self.species_cache[accession] = "Unknown_species"
                return "Unknown_species"
            
            seq_id = search_results["IdList"][0]
            
            # Get sequence summary to find organism
            handle = Entrez.esummary(db="nucleotide", id=seq_id)
            summary = Entrez.read(handle)[0]
            handle.close()
            
            organism = summary.get("Organism", "Unknown organism")
            cleaned_species = self.clean_species_name(organism)
            
            # Cache the result
            self.species_cache[accession] = cleaned_species
            
            time.sleep(0.3)  # Rate limiting
            return cleaned_species
            
        except Exception as e:
            self.log(f"Warning: Could not fetch species for {accession}: {e}")
            fallback = f"Unknown_{accession.replace('.', '_')}"
            self.species_cache[accession] = fallback
            return fallback
    
    def download_accession_sequences(self, accessions, gene_name, batch_size=20):
        """Download sequences for given accessions with species names"""
        self.log(f"Downloading {len(accessions)} accession references for {gene_name}")
        
        output_file = os.path.join(self.output_dir, f"{gene_name}_accession_references.fasta")
        
        downloaded = 0
        failed = 0
        
        # First, resolve species names for all accessions
        self.log("🔍 Resolving species names...")
        species_map = {}
        for i, acc in enumerate(accessions):
            if i % 10 == 0:
                self.log(f"Resolving species names: {i+1}/{len(accessions)}")
            species_map[acc] = self.get_species_from_accession(acc)
        
        # Process in batches to be efficient
        with open(output_file, "w") as outf:
            for i in range(0, len(accessions), batch_size):
                batch = accessions[i:i + batch_size]
                
                try:
                    # Batch download
                    handle = Entrez.efetch(
                        db="nucleotide",
                        id=",".join(batch),
                        rettype="fasta",
                        retmode="text"
                    )
                    
                    records = SeqIO.parse(handle, "fasta")
                    
                    for record in records:
                        # Extract accession from record ID
                        original_id = record.id
                        # Try to find matching accession
                        matched_acc = None
                        for acc in batch:
                            if acc in original_id or original_id.startswith(acc.split('.')[0]):
                                matched_acc = acc
                                break
                        
                        if not matched_acc:
                            # Fallback - use first accession in batch
                            matched_acc = batch[0]
                        
                        species_name = species_map.get(matched_acc, "Unknown_species")
                        
                        # Create enhanced header with species name
                        record.id = f"{species_name}_{matched_acc}_{gene_name}"
                        record.description = f"{gene_name} from {species_name.replace('_', ' ')} ({matched_acc})"
                        
                        SeqIO.write(record, outf, "fasta")
                        downloaded += 1
                        
                        if downloaded % 10 == 0:
                            self.log(f"Downloaded {downloaded}/{len(accessions)} sequences")
                    
                    handle.close()
                    time.sleep(0.5)  # Rate limiting
                    
                except Exception as e:
                    self.log(f"Error downloading batch starting with {batch[0]}: {e}")
                    # Try individual downloads for failed batch
                    for acc in batch:
                        try:
                            handle = Entrez.efetch(
                                db="nucleotide",
                                id=acc,
                                rettype="fasta",
                                retmode="text"
                            )
                            record = SeqIO.read(handle, "fasta")
                            handle.close()
                            
                            species_name = species_map.get(acc, "Unknown_species")
                            
                            # Create enhanced header with species name
                            record.id = f"{species_name}_{acc}_{gene_name}"
                            record.description = f"{gene_name} from {species_name.replace('_', ' ')} ({acc})"
                            
                            SeqIO.write(record, outf, "fasta")
                            downloaded += 1
                            
                            time.sleep(0.3)  # Rate limiting for individual requests
                            
                        except Exception as individual_error:
                            self.log(f"Failed to download {acc}: {individual_error}")
                            failed += 1
                    
                    continue
        
        self.log(f"Download complete: {downloaded} success, {failed} failed")
        
        if downloaded > 0:
            self.log(f"References saved to: {output_file}")
            return output_file
        else:
            self.log("No sequences downloaded successfully")
            return None

def main():
    parser = argparse.ArgumentParser(description="Download sequences from accession list as references with species names")
    parser.add_argument("--accessions", required=True,
                      help="Accession list (space/comma separated) or file path")
    parser.add_argument("--gene", required=True,
                      help="Gene name for output file")
    parser.add_argument("--email", required=True,
                      help="Email for NCBI API")
    parser.add_argument("--output", default="references",
                      help="Output directory")
    
    args = parser.parse_args()
    
    downloader = AccessionReferenceDownloader(args.email, args.output)
    
    # Parse accession list
    accessions = downloader.read_accession_list(args.accessions)
    
    if not accessions:
        print("❌ No accessions found in input")
        return 1
    
    print(f"📋 Processing {len(accessions)} accessions:")
    for acc in accessions[:5]:  # Show first 5
        print(f"  {acc}")
    if len(accessions) > 5:
        print(f"  ... and {len(accessions) - 5} more")
    
    # Download sequences
    result_file = downloader.download_accession_sequences(accessions, args.gene)
    
    if result_file:
        print(f"✅ Success! Reference file created: {result_file}")
        print("🧬 Headers now include species names: Species_name_Accession_Gene")
        return 0
    else:
        print("❌ Failed to create reference file")
        return 1

if __name__ == "__main__":
    exit(main())
