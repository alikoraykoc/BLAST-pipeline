#!/usr/bin/env python3

import argparse
import os
import time
import re
from Bio import Entrez, SeqIO
from collections import defaultdict

# Predefined taxonomic groups
TAXONOMIC_GROUPS = {
    "mammals": {
        "taxid": "40674",
        "name": "Mammalia",
        "description": "Mammals (Mammalia)"
    },
    "birds": {
        "taxid": "8782",
        "name": "Aves",
        "description": "Birds (Aves)"
    },
    "reptiles": {
        "taxid": "8504",
        "name": "Reptilia",
        "description": "Reptiles (Reptilia)"
    },
    "amphibians": {
        "taxid": "8292",
        "name": "Amphibia",
        "description": "Amphibians (Amphibia)"
    },
    "fish": {
        "taxid": "7898",
        "name": "Actinopterygii",
        "description": "Ray-finned fishes"
    },
    "arthropods": {
        "taxid": "6656",
        "name": "Arthropoda",
        "description": "Arthropods (insects, crustaceans, arachnids)"
    },
    "insects": {
        "taxid": "50557",
        "name": "Insecta",
        "description": "Insects"
    },
    "nematodes": {
        "taxid": "6231",
        "name": "Nematoda",
        "description": "Roundworms"
    },
    "plants": {
        "taxid": "33090",
        "name": "Viridiplantae",
        "description": "Green plants"
    },
    "fungi": {
        "taxid": "4751",
        "name": "Fungi",
        "description": "Fungi"
    },
    "bacteria": {
        "taxid": "2",
        "name": "Bacteria",
        "description": "Bacteria"
    },
    "archaea": {
        "taxid": "2157",
        "name": "Archaea",
        "description": "Archaea"
    },
    "vertebrates": {
        "taxid": "7742",
        "name": "Vertebrata",
        "description": "Vertebrates"
    },
    "primates": {
        "taxid": "9443",
        "name": "Primates",
        "description": "Primates"
    },
    "rodents": {
        "taxid": "9989",
        "name": "Rodentia",
        "description": "Rodents"
    }
}


class TaxonomicGeneCollector:
    def __init__(self, email, output_dir="gene_references"):
        self.email = email
        self.output_dir = output_dir
        Entrez.email = email
        os.makedirs(output_dir, exist_ok=True)

    def log(self, message):
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{timestamp}] {message}")

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

    def list_taxonomic_groups(self):
        """List available taxonomic groups"""
        print("🦠 Available taxonomic groups:")
        print("=" * 50)
        for key, info in TAXONOMIC_GROUPS.items():
            print(f"📊 {key:12} - {info['description']} (TaxID: {info['taxid']})")
        print()

    def parse_taxon_input(self, taxon_input):
        """Parse and validate taxonomic input - handles groups, TaxIDs, and scientific names"""
        valid_taxons = []
        
        for taxon in taxon_input:
            taxon = taxon.strip()
            
            # Check if it's a predefined group
            if taxon in TAXONOMIC_GROUPS:
                valid_taxons.append({
                    'type': 'group',
                    'input': taxon,
                    'taxid': TAXONOMIC_GROUPS[taxon]['taxid'],
                    'name': TAXONOMIC_GROUPS[taxon]['name'],
                    'description': TAXONOMIC_GROUPS[taxon]['description']
                })
            
            # Check if it's a numeric TaxID
            elif taxon.isdigit():
                valid_taxons.append({
                    'type': 'taxid',
                    'input': taxon,
                    'taxid': taxon,
                    'name': f"TaxID_{taxon}",
                    'description': f"Taxonomic ID {taxon}"
                })
            
            # Assume it's a scientific name
            elif any(char.isalpha() for char in taxon):
                # For scientific names, we'll use them directly in search
                valid_taxons.append({
                    'type': 'name',
                    'input': taxon,
                    'taxid': None,
                    'name': taxon.replace(' ', '_'),
                    'description': f"Scientific name: {taxon}"
                })
            
            else:
                self.log(f"⚠️ Invalid taxonomic input: {taxon}")
                continue
        
        return valid_taxons

    def search_gene_in_taxon(self, gene_name, taxon_info, max_sequences=100):
        """Search for a gene within a taxonomic group"""
        
        taxon_type = taxon_info['type']
        taxon_input = taxon_info['input']
        
        self.log(f"🔍 Searching '{gene_name}' in {taxon_info['description']}")

        # Construct search queries based on taxon type
        if taxon_type == 'group' or taxon_type == 'taxid':
            taxid = taxon_info['taxid']
            search_queries = [
                f"{gene_name}[Gene Name] AND txid{taxid}[Organism:exp]",
                f"{gene_name}[All Fields] AND txid{taxid}[Organism:exp] AND (mRNA[Title] OR CDS[Title] OR gene[Title])",
                f"{gene_name}[Title] AND txid{taxid}[Organism:exp]"
            ]
        else:  # scientific name
            scientific_name = taxon_input
            search_queries = [
                f"{gene_name}[Gene Name] AND {scientific_name}[Organism]",
                f"{gene_name}[All Fields] AND {scientific_name}[Organism] AND (mRNA[Title] OR CDS[Title] OR gene[Title])",
                f"{gene_name}[Title] AND {scientific_name}[Organism]"
            ]

        all_results = []

        for i, query in enumerate(search_queries):
            try:
                self.log(f"🎯 Query {i + 1}: {query}")

                handle = Entrez.esearch(
                    db="nucleotide",
                    term=query,
                    retmax=max_sequences,
                    sort="relevance"
                )
                search_result = Entrez.read(handle)
                handle.close()

                ids = search_result["IdList"]
                if ids:
                    self.log(f"✅ Found {len(ids)} sequences with query {i + 1}")
                    all_results.extend(ids)
                else:
                    self.log(f"⚠️  No results for query {i + 1}")

                time.sleep(0.5)  # Rate limiting

            except Exception as e:
                self.log(f"❌ Error with query {i + 1}: {e}")
                continue

        # Remove duplicates while preserving order
        unique_ids = []
        seen = set()
        for seq_id in all_results:
            if seq_id not in seen:
                unique_ids.append(seq_id)
                seen.add(seq_id)

        self.log(f"📊 Total unique sequences found: {len(unique_ids)}")
        return unique_ids

    def get_sequence_info(self, seq_ids, batch_size=50):
        """Get detailed information about sequences"""
        if not seq_ids:
            return []

        self.log(f"📋 Getting sequence information for {len(seq_ids)} sequences")

        sequence_info = []

        # Process in batches
        for i in range(0, len(seq_ids), batch_size):
            batch = seq_ids[i:i + batch_size]

            try:
                handle = Entrez.esummary(db="nucleotide", id=",".join(batch))
                summaries = Entrez.read(handle)
                handle.close()

                for summary in summaries:
                    try:
                        info = {
                            "accession": summary["AccessionVersion"],
                            "title": summary["Title"],
                            "organism": summary.get("Organism", "Unknown"),
                            "length": int(summary.get("Length", 0)),
                            "create_date": summary.get("CreateDate", ""),
                            "taxonomy": summary.get("TaxId", "")
                        }
                        sequence_info.append(info)
                    except Exception as e:
                        self.log(f"⚠️  Error processing summary: {e}")
                        continue

                time.sleep(0.5)  # Rate limiting

            except Exception as e:
                self.log(f"❌ Error getting batch info: {e}")
                continue

        return sequence_info

    def filter_sequences(self, sequence_info, min_length=100, max_length=50000):
        """Filter sequences by quality criteria"""
        filtered = []

        for seq in sequence_info:
            # Length filter
            if seq["length"] < min_length or seq["length"] > max_length:
                continue

            # Title quality check
            title_lower = seq["title"].lower()

            # Skip poor quality indicators
            skip_terms = [
                "partial", "fragment", "hypothetical", "predicted",
                "uncharacterized", "unknown", "synthetic", "artificial"
            ]

            if any(term in title_lower for term in skip_terms):
                continue

            # Prefer complete sequences
            good_terms = ["complete", "full", "mrna", "cds", "gene"]
            if any(term in title_lower for term in good_terms):
                seq["quality_score"] = 2
            else:
                seq["quality_score"] = 1

            # Prefer RefSeq
            if seq["accession"].startswith(("NM_", "NR_", "XM_", "XR_")):
                seq["quality_score"] += 1

            filtered.append(seq)

        # Sort by quality score, then by length
        filtered.sort(key=lambda x: (-x["quality_score"], -x["length"]))

        self.log(f"✅ Filtered to {len(filtered)} high-quality sequences")
        return filtered

    def download_sequences(self, sequence_info, gene_name, taxon_info, max_per_species=3):
        """Download the actual sequences"""
        if not sequence_info:
            return 0

        taxon_name = taxon_info['name']

        # Group by organism to limit per species
        by_organism = defaultdict(list)
        for seq in sequence_info:
            organism = seq["organism"].replace(" ", "_")
            by_organism[organism].append(seq)

        # Select best sequences per organism
        selected_sequences = []
        for organism, seqs in by_organism.items():
            # Take top sequences per organism
            selected_sequences.extend(seqs[:max_per_species])

        self.log(f"📥 Downloading {len(selected_sequences)} sequences from {len(by_organism)} organisms")

        # Create output file
        output_file = os.path.join(
            self.output_dir,
            f"{gene_name}_{taxon_name}_references.fasta"
        )

        downloaded = 0
        failed = 0

        with open(output_file, "w") as outf:
            for seq_info in selected_sequences:
                try:
                    accession = seq_info["accession"]
                    organism = seq_info["organism"]

                    # Download sequence
                    handle = Entrez.efetch(
                        db="nucleotide",
                        id=accession,
                        rettype="fasta",
                        retmode="text"
                    )

                    record = SeqIO.read(handle, "fasta")
                    handle.close()

                    # Create enhanced header with species name
                    clean_organism = self.clean_species_name(organism)
                    record.id = f"{clean_organism}_{accession}_{gene_name}"
                    record.description = f"{gene_name} from {organism} ({accession})"

                    SeqIO.write(record, outf, "fasta")
                    downloaded += 1

                    if downloaded % 10 == 0:
                        self.log(f"📥 Downloaded {downloaded} sequences...")

                    time.sleep(0.3)  # Rate limiting

                except Exception as e:
                    self.log(f"⚠️  Failed to download {accession}: {e}")
                    failed += 1
                    continue

        self.log(f"✅ Downloaded {downloaded} sequences to {output_file}")
        if failed > 0:
            self.log(f"⚠️  Failed to download {failed} sequences")

        return downloaded

    def collect_gene_references(self, gene_name, taxon_input, max_sequences=100,
                                min_length=100, max_length=50000, max_per_species=3):
        """Main function to collect gene references from taxonomic input"""

        self.log(f"🧬 Collecting '{gene_name}' references from taxonomic input: {taxon_input}")

        # Parse and validate taxonomic input
        valid_taxons = self.parse_taxon_input(taxon_input)
        
        if not valid_taxons:
            self.log("❌ No valid taxonomic groups provided")
            return 0

        total_downloaded = 0

        for taxon_info in valid_taxons:
            self.log(f"\n🦠 Processing: {taxon_info['description']}")

            try:
                # Search for sequences
                seq_ids = self.search_gene_in_taxon(gene_name, taxon_info, max_sequences)

                if not seq_ids:
                    self.log(f"❌ No sequences found for {gene_name} in {taxon_info['input']}")
                    continue

                # Get sequence information
                seq_info = self.get_sequence_info(seq_ids)

                if not seq_info:
                    self.log(f"❌ Could not retrieve sequence info for {taxon_info['input']}")
                    continue

                # Filter sequences
                filtered_info = self.filter_sequences(seq_info, min_length, max_length)

                if not filtered_info:
                    self.log(f"❌ No sequences passed quality filter for {taxon_info['input']}")
                    continue

                # Download sequences
                downloaded = self.download_sequences(
                    filtered_info, gene_name, taxon_info, max_per_species
                )

                total_downloaded += downloaded

                # Summary for this taxon
                organisms = set(seq["organism"] for seq in filtered_info)
                self.log(f"📊 {taxon_info['input']}: {downloaded} sequences from {len(organisms)} organisms")

            except Exception as e:
                self.log(f"❌ Error processing {taxon_info['input']}: {e}")
                continue

        self.log(f"\n🎉 Collection complete!")
        self.log(f"📊 Total sequences downloaded: {total_downloaded}")
        self.log(f"📂 Files saved in: {self.output_dir}")

        return total_downloaded


def main():
    parser = argparse.ArgumentParser(description="Collect gene references from taxonomic groups")
    parser.add_argument("--gene", required=True, help="Gene name to search for")
    parser.add_argument("--email", required=True, help="Email for NCBI API")
    parser.add_argument("--taxa", nargs="+", help="Taxonomic groups to search",
                        metavar="GROUP")
    parser.add_argument("--output", default="gene_references", help="Output directory")
    parser.add_argument("--max_sequences", type=int, default=100,
                        help="Maximum sequences to search per taxon")
    parser.add_argument("--min_length", type=int, default=100,
                        help="Minimum sequence length")
    parser.add_argument("--max_length", type=int, default=50000,
                        help="Maximum sequence length")
    parser.add_argument("--max_per_species", type=int, default=3,
                        help="Maximum sequences per species")
    parser.add_argument("--list_taxa", action="store_true",
                        help="List available taxonomic groups")

    args = parser.parse_args()

    collector = TaxonomicGeneCollector(args.email, args.output)

    if args.list_taxa:
        collector.list_taxonomic_groups()
        return

    if not args.taxa:
        print("❌ Please specify taxonomic groups with --taxa")
        print("💡 Use --list_taxa to see available groups")
        return

    # Run collection - DÜZELTME: taxon_input parametresi kullanılıyor
    total = collector.collect_gene_references(
        gene_name=args.gene,
        taxon_input=args.taxa,  # ✅ Doğru parametre adı
        max_sequences=args.max_sequences,
        min_length=args.min_length,
        max_length=args.max_length,
        max_per_species=args.max_per_species
    )

    if total > 0:
        print(f"\n🎯 Success! Use the generated files as references for BLAST searches.")
    else:
        print(f"\nNo sequences collected. Try different taxonomic groups or gene names.")


if __name__ == "__main__":
    main()
