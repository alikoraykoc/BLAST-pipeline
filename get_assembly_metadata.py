#!/usr/bin/env python3

import sys
import re
from Bio import Entrez
import time

def get_organism_from_assembly(accession, email):
    """Get organism name from assembly accession"""
    Entrez.email = email
    
    try:
        # Search assembly database
        handle = Entrez.esearch(db="assembly", term=f"{accession}[Assembly Accession]")
        search_results = Entrez.read(handle)
        handle.close()
        
        if not search_results["IdList"]:
            return None
            
        assembly_id = search_results["IdList"][0]
        
        # Get assembly summary
        handle = Entrez.esummary(db="assembly", id=assembly_id)
        summaries = Entrez.read(handle, validate=False)
        handle.close()
        
        if summaries and "DocumentSummarySet" in summaries:
            summary = summaries["DocumentSummarySet"]["DocumentSummary"][0]
            organism = summary.get("Organism", "")
            if organism:
                # Clean organism name for filename usage
                clean_organism = organism.replace(" ", "_").replace("/", "_").replace("(", "").replace(")", "")
                return clean_organism
    
    except Exception as e:
        print(f"Error getting metadata for {accession}: {e}", file=sys.stderr)
    
    return None

def extract_accession_from_filename(filename):
    """Extract assembly accession from filename"""
    # Pattern: GCA_123456789.1 or GCF_123456789.1
    pattern = r'(GC[AF]_\d+\.\d+)'
    match = re.search(pattern, filename)
    return match.group(1) if match else None

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 get_assembly_metadata.py <filename> <email>")
        sys.exit(1)
    
    filename = sys.argv[1]
    email = sys.argv[2]
    
    accession = extract_accession_from_filename(filename)
    if not accession:
        print("Unknown", end="")
        sys.exit(0)
    
    organism = get_organism_from_assembly(accession, email)
    if organism:
        print(organism, end="")
    else:
        print("Unknown", end="")
