#!/usr/bin/env python3
"""
TaxID Resolver - Retrieve NCBI taxids using hierarchical taxonomy.

Performs progressive search from species -> phylum, validates matches against
provided higher taxonomy, and handles multiple matches by scoring them.

Written by Dan Parsons @ Natural History Museum, London
Version: 1.0.0
"""

import argparse
import csv
import logging
import sys
import time
from functools import wraps
from http.client import IncompleteRead
from pathlib import Path
from random import uniform
from time import sleep
from typing import Any, Dict, List, Optional, Tuple
from urllib.error import HTTPError

from Bio import Entrez
from ratelimit import limits, sleep_and_retry

# =============================================================================
# Logging setup
# =============================================================================

logger = logging.getLogger("taxid_resolver")


def setup_logging(log_file: Optional[Path] = None):
    """Configure logging to console and optionally to file."""
    logger.setLevel(logging.DEBUG)
    
    # Console handler - INFO level
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_format = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    console_handler.setFormatter(console_format)
    logger.addHandler(console_handler)
    
    # File handler - DEBUG level
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_format = logging.Formatter(
            "%(asctime)s - %(levelname)s - %(funcName)s - %(message)s"
        )
        file_handler.setFormatter(file_format)
        logger.addHandler(file_handler)


# =============================================================================
# Retry decorator
# =============================================================================

def enhanced_retry(
    exceptions: tuple,
    tries: int = 4,
    initial_delay: int = 10,
    backoff: int = 2,
    max_delay: int = 240,
):
    """Retry decorator with exponential backoff."""
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            mdelay = initial_delay
            for i in range(tries):
                try:
                    return func(*args, **kwargs)
                except exceptions as e:
                    if i == tries - 1:
                        logger.error(f"Final attempt failed: {str(e)}")
                        return None

                    # Add jitter to avoid thundering herd
                    delay_with_jitter = mdelay + uniform(-0.1 * mdelay, 0.1 * mdelay)
                    logger.warning(
                        f"{str(e)}, Retrying in {delay_with_jitter:.2f} seconds..."
                    )
                    sleep(delay_with_jitter)

                    # Progressive backoff
                    mdelay = min(mdelay * backoff, max_delay)

            return None

        return wrapper

    return decorator


# =============================================================================
# TaxID Resolver class
# =============================================================================

class TaxIDResolver:
    """Handles NCBI Entrez API calls for taxonomy lookups."""
    
    def __init__(self, email: str, api_key: str):
        self.email = email
        self.api_key = api_key
        Entrez.email = email
        Entrez.api_key = api_key
        
        # In-memory taxonomy cache
        self.taxonomy_cache: Dict[str, Tuple[List[str], Dict[str, str], str, Dict[str, str]]] = {}
    
    @enhanced_retry((HTTPError, RuntimeError, IOError, IncompleteRead))
    @sleep_and_retry
    @limits(calls=10, period=1.1)
    def _fetch(self, **kwargs) -> Optional[Any]:
        """Rate-limited Entrez efetch wrapper."""
        result = Entrez.efetch(**kwargs)
        return result
    
    def _search(self, **kwargs) -> Optional[Dict]:
        """Entrez esearch with batching for large result sets."""
        @enhanced_retry((HTTPError, RuntimeError, IOError, IncompleteRead))
        @sleep_and_retry
        @limits(calls=10, period=1.1)
        def _do_search(**kwargs):
            handle = Entrez.esearch(**kwargs)
            result = Entrez.read(handle)
            handle.close()
            return result
        
        # First search to get total count
        initial_result = _do_search(**kwargs)
        if not initial_result:
            return None
        
        total_count = int(initial_result["Count"])
        all_ids = initial_result["IdList"]
        
        # If there are more results, fetch them in batches
        batch_size = 100
        if total_count > len(all_ids):
            for start in range(len(all_ids), total_count, batch_size):
                kwargs["retstart"] = start
                kwargs["retmax"] = batch_size
                result = _do_search(**kwargs)
                if result and result.get("IdList"):
                    all_ids.extend(result["IdList"])
                sleep(uniform(1, 2))
        
        initial_result["IdList"] = all_ids
        return initial_result
    
    @enhanced_retry(
        (HTTPError, RuntimeError, IOError, IncompleteRead),
        tries=5,
        initial_delay=15,
    )
    def fetch_taxonomy(
        self, taxid: str
    ) -> Tuple[List[str], Dict[str, str], str, Dict[str, str]]:
        """Fetch taxonomy information for a given taxid."""
        # Check cache first
        if taxid in self.taxonomy_cache:
            logger.debug(f"Using cached taxonomy for taxID: {taxid}")
            return self.taxonomy_cache[taxid]
        
        logger.debug(f"Fetching taxonomy from NCBI for taxID: {taxid}")
        
        # Verify taxid format
        taxid = taxid.strip()
        if not taxid.isdigit():
            logger.error(f"Invalid taxID format: {taxid} (must be numerical)")
            return [], {}, "", {}
        
        sleep(uniform(0.5, 1.0))
        
        try:
            params = {
                "db": "taxonomy",
                "id": taxid,
                "email": self.email,
                "api_key": self.api_key,
                "tool": "taxid_resolver",
            }
            
            max_retries = 3
            records = None
            for attempt in range(max_retries):
                try:
                    handle = self._fetch(**params)
                    if handle is None:
                        return [], {}, "", {}
                    records = Entrez.read(handle)
                    handle.close()
                    break
                except HTTPError as e:
                    if e.code == 400:
                        if attempt < max_retries - 1:
                            delay = (attempt + 1) * 2
                            logger.warning(
                                f"HTTP 400 error for taxID {taxid}, attempt "
                                f"{attempt + 1}/{max_retries}. Retrying in {delay}s..."
                            )
                            sleep(delay)
                            continue
                        else:
                            logger.error(
                                f"Failed to fetch taxonomy after {max_retries} "
                                f"attempts for taxID {taxid}"
                            )
                            return [], {}, "", {}
                    else:
                        raise
            
            if not records:
                logger.error(f"No taxonomy records found for taxID {taxid}")
                return [], {}, "", {}
            
            record = records[0]
            rank_info = {}
            taxid_info = {}
            lineage_nodes = record.get("LineageEx", [])
            lineage = []
            
            for node in lineage_nodes:
                name = node.get("ScientificName", "")
                rank = node.get("Rank", "no rank")
                node_taxid = str(node.get("TaxId", ""))
                
                lineage.append(name)
                
                if name and rank != "no rank":
                    rank_info[name] = rank
                    taxid_info[name] = node_taxid
            
            current_name = record.get("ScientificName", "")
            current_rank = record.get("Rank", "no rank")
            complete_lineage = lineage + [current_name]
            
            if current_rank != "no rank":
                rank_info[current_name] = current_rank
                taxid_info[current_name] = taxid
            
            logger.debug(f"Retrieved taxonomy for taxID {taxid}: {current_name} ({current_rank})")
            
            # Cache the result
            self.taxonomy_cache[taxid] = (
                complete_lineage,
                rank_info,
                current_rank,
                taxid_info,
            )
            
            return complete_lineage, rank_info, current_rank, taxid_info
        
        except IncompleteRead as e:
            logger.warning(f"IncompleteRead error for taxID {taxid}: {e}")
            raise
        except Exception as e:
            if isinstance(e, HTTPError) and e.code == 400:
                logger.error(f"HTTP 400 error for taxID {taxid}, skipping")
            else:
                logger.error(f"Error fetching taxonomy for taxID {taxid}: {e}")
            return [], {}, "", {}
    
    def validate_taxonomy_consistency(
        self, taxid: str, provided_taxonomy: Dict[str, str]
    ) -> Tuple[bool, Dict[str, Any]]:
        """Validate that a taxid's NCBI lineage is consistent with provided taxonomy."""
        logger.debug(f"Validating taxid {taxid} against provided taxonomy")
        
        try:
            complete_lineage, rank_info, current_rank, taxid_info = self.fetch_taxonomy(taxid)
            
            if not complete_lineage:
                logger.warning(f"No taxonomy data available for taxid {taxid}")
                return False, {"match_score": 0, "has_higher_conflict": False}
            
            # Build dictionary of fetched taxonomy by rank
            fetched_taxonomy = {}
            for taxon in complete_lineage:
                if taxon in rank_info:
                    rank = rank_info[taxon]
                    fetched_taxonomy[rank.lower()] = taxon
            
            if current_rank != "no rank":
                current_taxon = complete_lineage[-1] if complete_lineage else ""
                fetched_taxonomy[current_rank.lower()] = current_taxon
            
            logger.debug(f"NCBI taxonomy: {fetched_taxonomy}")
            logger.debug(f"Provided taxonomy: {provided_taxonomy}")
            
            match_score = 0
            match_details = {}
            
            # Check each rank
            ranks_to_check = ["phylum", "class", "order", "family", "genus"]
            for rank in ranks_to_check:
                if provided_taxonomy.get(rank) and rank in fetched_taxonomy:
                    matches = (
                        provided_taxonomy[rank].lower() == fetched_taxonomy[rank].lower()
                    )
                    match_details[rank] = "match" if matches else "mismatch"
                    if matches:
                        match_score += 1
                    else:
                        logger.debug(
                            f"{rank.capitalize()} mismatch for taxid {taxid}: "
                            f"provided '{provided_taxonomy[rank]}' vs "
                            f"fetched '{fetched_taxonomy[rank]}'"
                        )
            
            # Check for higher taxonomy conflicts (phylum/class)
            higher_taxonomy_conflict = False
            
            if (
                provided_taxonomy.get("phylum")
                and "phylum" in fetched_taxonomy
                and match_details.get("phylum") == "mismatch"
            ):
                higher_taxonomy_conflict = True
            
            if (
                provided_taxonomy.get("class")
                and "class" in fetched_taxonomy
                and match_details.get("class") == "mismatch"
            ):
                higher_taxonomy_conflict = True
            
            is_valid = not higher_taxonomy_conflict
            
            if is_valid:
                logger.debug(f"Taxid {taxid} passed validation (score: {match_score})")
            else:
                logger.debug(f"Taxid {taxid} failed validation (higher taxonomy conflict)")
            
            return is_valid, {
                "match_score": match_score,
                "details": match_details,
                "lineage": complete_lineage,
                "fetched_taxonomy": fetched_taxonomy,
                "has_higher_conflict": higher_taxonomy_conflict,
            }
        
        except Exception as e:
            logger.error(f"Error validating taxonomy for taxid {taxid}: {e}")
            return False, {"match_score": 0, "has_higher_conflict": False}
    
    def _search_and_validate_taxids(
        self,
        search_term: str,
        rank_name: str,
        provided_taxonomy: Dict[str, str],
    ) -> Tuple[Optional[str], str, bool]:
        """Search for taxids and validate them against provided taxonomy."""
        try:
            result = self._search(db="taxonomy", term=search_term)
            if not result or not result.get("IdList"):
                return None, "", False
            
            taxids = result["IdList"]
            
            if len(taxids) > 1:
                logger.info(f"Multiple taxids found for {rank_name}: {taxids}")
                
                # Validate each and find best match
                valid_taxids = []
                for taxid in taxids:
                    is_valid, details = self.validate_taxonomy_consistency(
                        taxid, provided_taxonomy
                    )
                    if is_valid:
                        valid_taxids.append((taxid, details))
                
                if not valid_taxids:
                    logger.warning(
                        f"None of the taxids for {rank_name} match provided higher taxonomy"
                    )
                    # Return best match anyway but flag the mismatch
                    best_taxid = taxids[0]
                    _, details = self.validate_taxonomy_consistency(
                        best_taxid, provided_taxonomy
                    )
                    return best_taxid, rank_name, True
                elif len(valid_taxids) == 1:
                    return valid_taxids[0][0], rank_name, False
                else:
                    # Multiple valid - use best match score
                    best = max(valid_taxids, key=lambda x: x[1]["match_score"])
                    return best[0], rank_name, False
            else:
                # Single match - validate it
                taxid = taxids[0]
                is_valid, details = self.validate_taxonomy_consistency(
                    taxid, provided_taxonomy
                )
                has_conflict = details.get("has_higher_conflict", False)
                
                if is_valid:
                    logger.info(f"Found valid taxid {taxid} for {rank_name}")
                    return taxid, rank_name, False
                else:
                    logger.warning(
                        f"Taxid {taxid} for {rank_name} has higher taxonomy conflict"
                    )
                    # Return it anyway but flag the mismatch
                    return taxid, rank_name, True
        
        except Exception as e:
            logger.error(f"Error searching for {rank_name} taxid: {e}")
            return None, "", False
    
    def resolve_taxid_from_taxonomy(
        self,
        phylum: str,
        class_name: str,
        order: str,
        family: str,
        genus: str,
        species: str,
    ) -> Tuple[Optional[str], str, bool]:
        """Resolve NCBI taxID from hierarchical taxonomic information."""
        logger.info(
            f"Resolving taxid for: {genus} {species} "
            f"(Family: {family}, Order: {order}, Class: {class_name}, Phylum: {phylum})"
        )
        
        # Store provided taxonomy for validation
        provided_taxonomy = {
            "phylum": phylum.strip() if phylum else "",
            "class": class_name.strip() if class_name else "",
            "order": order.strip() if order else "",
            "family": family.strip() if family else "",
            "genus": genus.strip() if genus else "",
            "species": species.strip() if species else "",
        }
        
        # Try species level first
        if genus and species:
            # Handle species that already contains genus name
            if species.startswith(genus):
                full_species = species
            else:
                full_species = f"{genus} {species}"
            
            search_term = f"{full_species}[Scientific Name]"
            logger.info(f"Searching species: {search_term}")
            
            taxid, rank, mismatch = self._search_and_validate_taxids(
                search_term, "species", provided_taxonomy
            )
            if taxid:
                return taxid, rank, mismatch
        
        # Try genus level
        if genus:
            search_term = f"{genus}[Scientific Name] AND Metazoa[Organism]"
            logger.info(f"Searching genus: {search_term}")
            
            taxid, rank, mismatch = self._search_and_validate_taxids(
                search_term, "genus", provided_taxonomy
            )
            if taxid:
                return taxid, rank, mismatch
        
        # Try family level
        if family:
            search_term = f"{family}[Scientific Name] AND Metazoa[Organism]"
            logger.info(f"Searching family: {search_term}")
            
            taxid, rank, mismatch = self._search_and_validate_taxids(
                search_term, "family", provided_taxonomy
            )
            if taxid:
                return taxid, rank, mismatch
        
        # Try order level
        if order:
            search_term = f"{order}[Scientific Name] AND Metazoa[Organism]"
            logger.info(f"Searching order: {search_term}")
            
            taxid, rank, mismatch = self._search_and_validate_taxids(
                search_term, "order", provided_taxonomy
            )
            if taxid:
                return taxid, rank, mismatch
        
        # Try class level
        if class_name:
            search_term = f"{class_name}[Scientific Name] AND Metazoa[Organism]"
            logger.info(f"Searching class: {search_term}")
            
            taxid, rank, mismatch = self._search_and_validate_taxids(
                search_term, "class", provided_taxonomy
            )
            if taxid:
                return taxid, rank, mismatch
        
        # Try phylum level (minimal validation at this level)
        if phylum:
            search_term = f"{phylum}[Scientific Name] AND Metazoa[Organism]"
            logger.info(f"Searching phylum: {search_term}")
            
            try:
                result = self._search(db="taxonomy", term=search_term)
                if result and result.get("IdList"):
                    taxids = result["IdList"]
                    if len(taxids) > 1:
                        logger.warning(f"Multiple taxids found for phylum {phylum}: {taxids}")
                    taxid = taxids[0]
                    logger.info(f"Found taxid {taxid} at phylum level")
                    return taxid, "phylum", False
            except Exception as e:
                logger.error(f"Error searching for phylum taxid: {e}")
        
        logger.warning(f"Could not find any valid taxid for {genus} {species}")
        return None, "", False


# =============================================================================
# CSV processing
# =============================================================================

def process_taxonomy_csv(
    input_path: Path,
    output_path: Path,
    resolver: TaxIDResolver,
) -> None:
    """Process input CSV and add taxid, matched_rank, and lineage_mismatch columns."""
    logger.info(f"Reading input CSV: {input_path}")
    
    # Read input CSV
    with open(input_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        
        if not fieldnames:
            raise ValueError("Input CSV has no headers")
        
        # Validate required columns (case-insensitive check)
        fieldnames_lower = [fn.lower() for fn in fieldnames]
        required_columns = ["id", "phylum", "class", "order", "family", "genus", "species"]
        
        missing = [col for col in required_columns if col not in fieldnames_lower]
        if missing:
            raise ValueError(f"Input CSV missing required columns: {missing}")
        
        # Create mapping from lowercase to actual column names
        column_map = {fn.lower(): fn for fn in fieldnames}
        
        rows = list(reader)
    
    logger.info(f"Processing {len(rows)} samples...")
    
    # Process each row
    output_rows = []
    for i, row in enumerate(rows, 1):
        sample_id = row.get(column_map["id"], "")
        phylum = row.get(column_map["phylum"], "")
        class_name = row.get(column_map["class"], "")
        order = row.get(column_map["order"], "")
        family = row.get(column_map["family"], "")
        genus = row.get(column_map["genus"], "")
        species = row.get(column_map["species"], "")
        
        logger.info(f"[{i}/{len(rows)}] Processing sample: {sample_id}")
        
        taxid, matched_rank, lineage_mismatch = resolver.resolve_taxid_from_taxonomy(
            phylum=phylum,
            class_name=class_name,
            order=order,
            family=family,
            genus=genus,
            species=species,
        )
        
        # Add new columns
        row["taxid"] = taxid if taxid else ""
        row["matched_rank"] = matched_rank if matched_rank else ""
        row["lineage_mismatch"] = "YES" if lineage_mismatch else "NO"
        
        output_rows.append(row)
        
        if taxid:
            logger.info(
                f"  -> taxid: {taxid}, matched_rank: {matched_rank}, "
                f"lineage_mismatch: {'YES' if lineage_mismatch else 'NO'}"
            )
        else:
            logger.warning(f"  -> No taxid found for sample {sample_id}")
    
    # Write output CSV
    output_fieldnames = list(fieldnames) + ["taxid", "matched_rank", "lineage_mismatch"]
    
    logger.info(f"Writing output CSV: {output_path}")
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=output_fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)
    
    logger.info(f"Done! Processed {len(output_rows)} samples.")


# =============================================================================
# CLI
# =============================================================================

def setup_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Resolve NCBI taxids from hierarchical taxonomic information.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python taxid_resolver.py -i samples.csv -e your@email.com -k YOUR_API_KEY

Input CSV must have columns: ID, phylum, class, order, family, genus, species
Output overwrites input CSV with additional columns: taxid, matched_rank, lineage_mismatch
        """,
    )
    
    parser.add_argument(
        "--input", "-i",
        required=True,
        type=Path,
        help="Path to input CSV file with taxonomy columns (will be overwritten with results)",
    )
    
    parser.add_argument(
        "--email", "-e",
        required=True,
        type=str,
        help="Email for NCBI API requests (required by NCBI)",
    )
    
    parser.add_argument(
        "--api-key", "-k",
        required=True,
        type=str,
        help="NCBI API key for increased rate limits",
    )
    
    parser.add_argument(
        "--log-file", "-l",
        type=Path,
        default=None,
        help="Optional log file path (logs DEBUG level messages)",
    )
    
    return parser


def validate_credentials(email: str, api_key: str) -> bool:
    """Validate NCBI credentials by making a test request."""
    if not email or "@" not in email:
        raise ValueError("Invalid email address")
    
    if not api_key or len(api_key) < 10:
        raise ValueError("Invalid API key (too short)")
    
    Entrez.email = email
    Entrez.api_key = api_key
    
    try:
        handle = Entrez.einfo()
        _ = Entrez.read(handle)
        handle.close()
        return True
    except Exception as e:
        raise ValueError(f"NCBI credential validation failed: {e}")


def main():
    print("=" * 60)
    print("         TaxID Resolver - NCBI Taxonomy Lookup")
    print("=" * 60)
    print()
    
    parser = setup_argument_parser()
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.log_file)
    
    # Validate input file exists
    if not args.input.exists():
        logger.error(f"Input file not found: {args.input}")
        sys.exit(1)
    
    # Validate credentials
    logger.info("Validating NCBI credentials...")
    try:
        validate_credentials(args.email, args.api_key)
        logger.info("Credential validation passed")
    except ValueError as e:
        logger.error(f"Credential validation failed: {e}")
        sys.exit(1)
    
    # Create resolver
    resolver = TaxIDResolver(email=args.email, api_key=args.api_key)
    
    # Process CSV (overwrites input file)
    try:
        process_taxonomy_csv(
            input_path=args.input,
            output_path=args.input,
            resolver=resolver,
        )
    except ValueError as e:
        logger.error(f"Error processing CSV: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        logger.exception("Full traceback:")
        sys.exit(1)
    
    print()
    print("=" * 60)
    print("              TaxID Resolver complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
