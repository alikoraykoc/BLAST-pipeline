#!/usr/bin/env python3
"""
Enhanced species naming utility for FASTA headers
Fetches organism names from NCBI and creates clean, informative headers
"""

import re
import time
import sys
from Bio import Entrez
from collections import defaultdict
import pickle
import os

class SpeciesNameResolver:
    def __init__(self, email, cache_file="species_cache.pkl"):
        self.email = email
        self.cache_file = cache_file
        Entrez.email = email
        
        # Load cached species names
        self.species_cache = self._load_cache()
        
    def _load_cache(self):
        """Load species name cache from file"""
        if os.path.exists(self.cache_file):
            try:
                with open(self.cache_file, 'rb') as f:
                    return pickle.load(f)
            except:
                pass
        return {}
    
    def _save_cache(self):
        """Save species name cache to file"""
        try:
            with open(self.cache_file, 'wb') as f:
                pickle.dump(self.species_cache, f)
        except Exception as e:
            print(f"Warning: Could not save cache: {e}")
    
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
            # Try nucleotide database first
            species_name = self._fetch_from_nucleotide(accession)
            
            # If not found, try assembly database
            if not species_name:
                species_name = self._fetch_from_assembly(accession)
            
            # Cache the result
            cleaned_name = self.clean_species_name(species_name) if species_name else "Unknown_species"
            self.species_cache[accession] = cleaned_name
            
            return cleaned_name
            
        except Exception as e:
            print(f"Warning: Could not fetch species for {accession}: {e}")
            fallback = f"Unknown_{accession.replace('.', '_')}"
            self.species_cache[accession] = fallback
            return fallback
    
    def _fetch_from_nucleotide(self, accession):
        """Fetch species name from nucleotide database"""
        try:
            handle = Entrez.esearch(db="nucleotide", term=f"{accession}[Accession]")
            search_results = Entrez.read(handle)
            handle.close()
            
            if not search_results["IdList"]:
                return None
            
            seq_id = search_results["IdList"][0]
            
            handle = Entrez.esummary(db="nucleotide", id=seq_id)
            summary = Entrez.read(handle)[0]
            handle.close()
            
            time.sleep(0.3)  # Rate limiting
            
            return summary.get("Organism", None)
            
        except:
            return None
    
    def _fetch_from_assembly(self, accession):
        """Fetch species name from assembly database"""
        try:
            # Extract base accession (remove version if present)
            base_accession = accession.split('.')[0]
            
            handle = Entrez.esearch(db="assembly", term=f"{base_accession}[Assembly Accession]")
            search_results = Entrez.read(handle)
            handle.close()
            
            if not search_results["IdList"]:
                return None
            
            assembly_id = search_results["IdList"][0]
            
            handle = Entrez.esummary(db="assembly", id=assembly_id)
            summaries = Entrez.read(handle, validate=False)
            handle.close()
            
            time.sleep(0.3)  # Rate limiting
            
            if summaries and "DocumentSummarySet" in summaries:
                summary = summaries["DocumentSummarySet"]["DocumentSummary"][0]
                return summary.get("Organism", None)
                
        except:
            return None
        
        return None
    
    def batch_resolve_species(self, accessions, max_batch=20):
        """Resolve species names for multiple accessions efficiently"""
        results = {}
        uncached = []
        
        # Check cache first
        for acc in accessions:
            if acc in self.species_cache:
                results[acc] = self.species_cache[acc]
            else:
                uncached.append(acc)
        
        print(f"🔍 Resolving species names: {len(results)} cached, {len(uncached)} new")
        
        # Process uncached accessions in batches
        for i in range(0, len(uncached), max_batch):
            batch = uncached[i:i + max_batch]
            print(f"📡 Processing batch {i//max_batch + 1}/{(len(uncached)-1)//max_batch + 1}")
            
            for acc in batch:
                species = self.get_species_from_accession(acc)
                results[acc] = species
                
            time.sleep(1)  # Rate limiting between batches
        
        # Save updated cache
        self._save_cache()
        
        return results
    
    def extract_accession_from_text(self, text):
        """Extract accession numbers from text (headers, filenames, etc.)"""
        patterns = [
            r'\b(GC[AF]_\d+\.\d+)\b',  # Assembly accessions
            r'\b([A-Z]{1,2}\d{5,6}\.\d+)\b',  # GenBank accessions
            r'\b(NM_\d+\.\d+)\b',  # RefSeq mRNA
            r'\b(NR_\d+\.\d+)\b',  # RefSeq ncRNA
            r'\b(XM_\d+\.\d+)\b',  # RefSeq predicted mRNA
            r'\b(XR_\d+\.\d+)\b',  # RefSeq predicted ncRNA
        ]
        
        for pattern in patterns:
            match = re.search(pattern, text)
            if match:
                return match.group(1)
        
        return None


def enhance_fasta_headers(input_fasta, output_fasta, email, gene_name=None):
    """Enhance FASTA headers with species names"""
    print(f"🧬 Enhancing FASTA headers: {input_fasta}")
    
    resolver = SpeciesNameResolver(email)
    
    # First pass: collect all accessions
    accessions = []
    with open(input_fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                acc = resolver.extract_accession_from_text(line)
                if acc:
                    accessions.append(acc)
    
    if not accessions:
        print("⚠️ No accessions found in FASTA file")
        return False
    
    # Resolve species names
    species_map = resolver.batch_resolve_species(accessions)
    
    # Second pass: rewrite with enhanced headers
    enhanced_count = 0
    with open(input_fasta, 'r') as inf, open(output_fasta, 'w') as outf:
        for line in inf:
            if line.startswith('>'):
                acc = resolver.extract_accession_from_text(line.strip())
                if acc and acc in species_map:
                    species = species_map[acc]
                    
                    # Create enhanced header
                    if gene_name:
                        new_header = f">{species}_{acc}_{gene_name}"
                    else:
                        new_header = f">{species}_{acc}"
                    
                    outf.write(new_header + "\n")
                    enhanced_count += 1
                else:
                    outf.write(line)  # Keep original if no accession found
            else:
                outf.write(line)
    
    print(f"✅ Enhanced {enhanced_count} headers in {output_fasta}")
    return True


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Enhance FASTA headers with species names")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output FASTA file")
    parser.add_argument("--email", required=True, help="Email for NCBI API")
    parser.add_argument("--gene", help="Gene name to add to headers")
    
    args = parser.parse_args()
    
    success = enhance_fasta_headers(args.input, args.output, args.email, args.gene)
    
    if success:
        print("🎉 Header enhancement completed!")
    else:
        print("❌ Header enhancement failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()