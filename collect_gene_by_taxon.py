#!/usr/bin/env python3

import argparse
import os
import time
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

    def list_taxonomic_groups(self):
        """List available taxonomic groups"""
        print("🦠 Available taxonomic groups:")
        print("=" * 50)
        for key, info in TAXONOMIC_GROUPS.items():
            print(f"📊 {key:12} - {info['description']} (TaxID: {info['taxid']})")
        print()

    def search_gene_in_taxon(self, gene_name, taxon_group, max_sequences=100):
        """Search for a gene within a taxonomic group"""

        if taxon_group not in TAXONOMIC_GROUPS:
            raise ValueError(f"Unknown taxonomic group: {taxon_group}")

        taxon_info = TAXONOMIC_GROUPS[taxon_group]
        taxid = taxon_info["taxid"]

        self.log(f"🔍 Searching '{gene_name}' in {taxon_info['description']} (TaxID: {taxid})")

        # Construct search query
        search_queries = [
            f"{gene_name}[Gene Name] AND txid{taxid}[Organism:exp]",
            f"{gene_name}[All Fields] AND txid{taxid}[Organism:exp] AND (mRNA[Title] OR CDS[Title] OR gene[Title])",
            f"{gene_name}[Title] AND txid{taxid}[Organism:exp]"
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

    def download_sequences(self, sequence_info, gene_name, taxon_group, max_per_species=3):
        """Download the actual sequences"""
        if not sequence_info:
            return 0

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
            f"{gene_name}_{taxon_group}_references.fasta"
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

                    # Create clean header
                    clean_organism = organism.replace(" ", "_")
                    record.id = f"{clean_organism}_{accession}"
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

    def collect_gene_references(self, gene_name, taxon_groups, max_sequences=100,
                                min_length=100, max_length=50000, max_per_species=3):
        """Main function to collect gene references from taxonomic groups"""

        self.log(f"🧬 Collecting '{gene_name}' references from taxonomic groups: {taxon_groups}")

        total_downloaded = 0

        for taxon_group in taxon_groups:
            self.log(f"\n🦠 Processing taxonomic group: {taxon_group}")

            try:
                # Search for sequences
                seq_ids = self.search_gene_in_taxon(gene_name, taxon_group, max_sequences)

                if not seq_ids:
                    self.log(f"❌ No sequences found for {gene_name} in {taxon_group}")
                    continue

                # Get sequence information
                seq_info = self.get_sequence_info(seq_ids)

                if not seq_info:
                    self.log(f"❌ Could not retrieve sequence info for {taxon_group}")
                    continue

                # Filter sequences
                filtered_info = self.filter_sequences(seq_info, min_length, max_length)

                if not filtered_info:
                    self.log(f"❌ No sequences passed quality filter for {taxon_group}")
                    continue

                # Download sequences
                downloaded = self.download_sequences(
                    filtered_info, gene_name, taxon_group, max_per_species
                )

                total_downloaded += downloaded

                # Summary for this taxon
                organisms = set(seq["organism"] for seq in filtered_info)
                self.log(f"📊 {taxon_group}: {downloaded} sequences from {len(organisms)} organisms")

            except Exception as e:
                self.log(f"❌ Error processing {taxon_group}: {e}")
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
                        choices=list(TAXONOMIC_GROUPS.keys()))
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

    # Validate taxonomic groups
    invalid_taxa = []
    for t in args.taxa:
        if not (t in TAXONOMIC_GROUPS or t.isdigit() or any(char.isalpha() for char in t)):
            invalid_taxa.append(t)

    if invalid_taxa:
        print(f"❌ Invalid taxonomic input: {invalid_taxa}")
        print("💡 Use:")
        print("  - Predefined groups (--list_taxa to see options)")
        print("  - Numeric TaxIDs (e.g., 6656 for Arthropoda)")
        print("  - Scientific names (e.g., 'Orthoptera')")
        collector.list_taxonomic_groups()
        return

    # Run collection
    total = collector.collect_gene_references(
        gene_name=args.gene,
        taxon_input=args.taxa,
        max_sequences=args.max_sequences,
        min_length=args.min_length,
        max_length=args.max_length,
        max_per_species=args.max_per_species
    )

    if total > 0:
        print(f"\n🎯 Success! Use the generated files as references for BLAST searches.")
    else:
        print(f"\n😞 No sequences collected. Try different taxonomic groups or gene names.")


if __name__ == "__main__":
    main()