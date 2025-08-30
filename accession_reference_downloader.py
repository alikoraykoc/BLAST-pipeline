#!/usr/bin/env python3

import argparse
import os
import time
from Bio import Entrez, SeqIO

class AccessionReferenceDownloader:
    def __init__(self, email, output_dir="references"):
        self.email = email
        self.output_dir = output_dir
        Entrez.email = email
        os.makedirs(output_dir, exist_ok=True)
        
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
    
    def download_accession_sequences(self, accessions, gene_name, batch_size=20):
        """Download sequences for given accessions"""
        self.log(f"Downloading {len(accessions)} accession references for {gene_name}")
        
        output_file = os.path.join(self.output_dir, f"{gene_name}_accession_references.fasta")
        
        downloaded = 0
        failed = 0
        
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
                        # Clean up the header
                        original_id = record.id
                        record.id = f"{original_id}_{gene_name}"
                        record.description = f"{gene_name} reference from {original_id}"
                        
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
                            
                            record.id = f"{acc}_{gene_name}"
                            record.description = f"{gene_name} reference from {acc}"
                            
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
    parser = argparse.ArgumentParser(description="Download sequences from accession list as references")
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
        return 0
    else:
        print("❌ Failed to create reference file")
        return 1

if __name__ == "__main__":
    exit(main())