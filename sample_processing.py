#!/usr/bin/env python3

"""
Sample Processing Script:

This script processes raw DNA sequencing read files, and matches them with sample metadata, 
supporting flexible directory structures and file naming patterns. 
It features progress tracking and detailed logging of the processing steps.

Features:
-------------
- Can process files directly in the input directory (flat)
- Can recursively search subdirectories (nested)
- Identifies IDs in both filenames and directory paths
- Supports various delimiter patterns (_, -, /)
- Handles multiple fastq file extensions (.fastq, .fq, w/wo .gz)
- Automatically excludes 'Undetermined' and 'NC' (negative control) files
- Logs both to console and file
- Use -op argument to specify custom output file prefix (Default to directory name
 if no prefix specified)
- Use --merge flag to prefer *_merged directories when available (Falls back to non-suffixed 
 directories if no merged version exists & warns about duplicate IDs 
 when --merge is not specified)

Output Files:
-------------
1. <prefix>.csv: All matched samples
2. <prefix>_types.csv: Samples marked as types

Each output file contains:
- ID: Unique ID of the sample
- forward: Absolute path to R1 (forward) read file
- reverse: Absolute path to R2 (reverse) read file
- phylum: Phylum from metadata
- class: Class from metadata
- order: Order from metadata
- family: Family from metadata
- genus: Genus from metadata
- species: Species from metadata
- type_status: Type specimen status from metadata

Usage:
------
python 4_sample_sheet.py [parent_directory] [sample_metadata_file] [output_directory] [-op output_prefix] [--merge] [--identifier column_name]

Arguments:
    parent_directory: Directory containing FASTQ/FQ files (directly or in subdirectories)
    sample_metadata_file: CSV file with sample metadata
    output_directory: Directory where output files will be saved (will be created if it doesn't exist)
    -op output_prefix: Optional prefix for output files (default: uses directory name)
    --merge: Prefer *_merged directories when available, fall back to non-suffixed directories
    --identifier/--id: Column name in metadata CSV containing sample IDs (default: "Process ID")

Required Metadata CSV Columns:
    - <identifier column>: Unique identifier for each sample (default column name: "Process ID")
    - phylum: Phylum classification
    - class: Class classification
    - order: Order classification
    - family: Family classification
    - genus: Genus classification
    - species: Species classification
    - type_status: Type specimen status

Example Directory Structures Supported:
1. Flat structure:
   /parent_dir/
   ├── BSNHM593-24_R1.fq
   ├── BSNHM593-24_R2.fq
   └── ...

2. Nested structure:
   /parent_dir/XE-4013/
   └── 20240906_LH00179_0123_A22CKGHLT4/
       ├── Sample_XE-4013-BGSNL096-23/
       │   ├── BGSNL096-23_R1.fastq.gz
       │   └── BGSNL096-23_R2.fastq.gz
       └── ...

Dependencies:
    - Python 3.x
    - pandas
    - tqdm (for progress bars)
    - Standard library modules: os, csv, sys, logging, argparse
"""



import os
import csv
import sys
import pandas as pd
import logging
import argparse
from tqdm import tqdm

# Taxonomic columns to extract from metadata
TAXONOMY_COLUMNS = ['phylum', 'class', 'order', 'family', 'genus', 'species']


def setup_logging(output_dir, log_filename):
    """Set up logging to both file and console"""
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    file_handler = logging.FileHandler(os.path.join(output_dir, log_filename))
    file_handler.setLevel(logging.DEBUG)
    
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    
    formatter = logging.Formatter('%(message)s')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)
    
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    return logger


def get_directory_name(path):
    path = path.rstrip(os.path.sep)
    return os.path.basename(path)


def is_fastq_file(filename):
    """Check if a filename corresponds to a fastq file."""
    fastq_extensions = ('.fastq.gz', '.fq.gz', '.fastq', '.fq')
    filename_lower = filename.lower()
    return any(filename_lower.endswith(ext) for ext in fastq_extensions)


def get_directory_suffix(dirpath):
    """
    Extract suffix from directory name.
    Returns:
        - 'merged' if directory ends with '_merged'
        - 'numbered' if directory ends with '_<number>'
        - 'base' if no suffix
    """
    dirname = os.path.basename(dirpath)
    if dirname.endswith('_merged'):
        return 'merged'
    elif '_' in dirname:
        last_part = dirname.split('_')[-1]
        if last_part.isdigit():
            return 'numbered'
    return 'base'


def find_id_in_string(string, ids, logger=None):
    """
    Find an ID in a string.
    Handles IDs in filenames, directory names, and compound IDs.
    """
    if logger:
        logger.debug("Trying to find ID in string: %s", string)
    
    # Clean up the string & try each ID
    base_string = os.path.splitext(os.path.splitext(string)[0])[0]  
    
    for sid in ids:
        sid_str = str(sid)  # Ensure ID is a string
        
        if not sid_str.strip():
            continue
            
        if sid_str in base_string:
            # Get the index where the ID was found
            idx = base_string.index(sid_str)
            
            before = base_string[idx-1] if idx > 0 else None
            after = base_string[idx+len(sid_str)] if idx+len(sid_str) < len(base_string) else None
            
            # ID is valid if it's bounded by separators or string boundaries
            if (before is None or before in ['_', '-', '/', ' ']) and \
               (after is None or after in ['_', '-', '/', ' ']):
                if logger:
                    logger.debug(f"Found ID: {sid_str}")
                    logger.debug(f"Position: {idx}")
                    logger.debug(f"Before: {before}, After: {after}")
                return sid_str
    
    if logger:
        logger.debug("No ID found")
    return None


def extract_id_reads_taxonomy_type(parent_dir, metadata_df, logger, use_merge=False, id_column='Process ID'):
    results = {}
    ids = set(metadata_df[id_column].astype(str))
    
    logger.debug("\nDEBUG: Number of IDs in metadata: %d", len(ids))
    logger.debug("DEBUG: First few IDs from metadata: %s", list(ids)[:5])
    logger.debug("DEBUG: Merge mode enabled: %s", use_merge)
    
    # Build dictionaries for each taxonomy column
    taxonomy_dicts = {}
    for col in TAXONOMY_COLUMNS:
        if col in metadata_df.columns:
            taxonomy_dicts[col] = dict(zip(metadata_df[id_column].astype(str), metadata_df[col].fillna('')))
        else:
            taxonomy_dicts[col] = {}
            logger.warning("WARNING: Column '%s' not found in metadata, will use empty values", col)
    
    type_status_dict = dict(zip(metadata_df[id_column].astype(str), metadata_df['type_status']))
    
    # First check if there are fastq files in the input directory
    input_dir_files = os.listdir(parent_dir)
    fastq_files_in_input = [f for f in input_dir_files if is_fastq_file(f)]
    using_subdirs = len(fastq_files_in_input) == 0
    
    logger.debug("\nDEBUG: FASTQ files found in input directory: %d", len(fastq_files_in_input))
    logger.debug("DEBUG: Will search in subdirectories: %s", using_subdirs)
    
    if not using_subdirs:
        logger.debug("\nDEBUG: Processing files in input directory")
        r1_files = {}
        r2_files = {}
        
        for filename in tqdm(input_dir_files, desc="Processing input files"):
            if not is_fastq_file(filename) or "Undetermined" in filename or "NC" in filename:
                continue
            
            logger.debug("\nDEBUG: Processing file: %s", filename)
            sample_id = find_id_in_string(filename, ids, logger)
            
            if not sample_id:
                logger.debug("DEBUG: No ID found in filename: %s", filename)
                continue
            
            logger.debug("DEBUG: Found ID: %s", sample_id)
                
            filepath = os.path.abspath(os.path.join(parent_dir, filename))
            logger.debug("File path: %s", filepath)
            if "R1" in filename:
                logger.debug("Found R1 file - adding to r1_files with ID: %s", sample_id)
                r1_files[sample_id] = filepath
            elif "R2" in filename:
                logger.debug("Found R2 file - adding to r2_files with ID: %s", sample_id)
                r2_files[sample_id] = filepath
            else:
                logger.debug("No R1/R2 pattern found in filename: %s", filename)
    else:
        logger.debug("\nDEBUG: Searching in subdirectories")
        
        # Store all candidate files with their directory info
        # Structure: {sample_id: {dir_suffix: {'R1': path, 'R2': path, 'dir': dirpath}}}
        candidate_files = {}
        
        # First, count total files for progress bar
        total_files = sum([len([f for f in files if is_fastq_file(f)])
                         for _, _, files in os.walk(parent_dir)])
        
        processed_files = 0
        pbar = tqdm(total=total_files, desc="Processing FASTQ files")
        
        for root, _, files in os.walk(parent_dir):
            fastq_files = [f for f in files if is_fastq_file(f)]
            if not fastq_files:
                continue
                
            logger.debug("\nDEBUG: Found %d FASTQ files in %s", len(fastq_files), root)
            
            # Get directory suffix type
            dir_suffix = get_directory_suffix(root)
            
            for filename in fastq_files:
                if "Undetermined" in filename or "NC" in filename:
                    continue
                
                # Try to find ID in both filename and directory path
                sample_id = find_id_in_string(filename, ids)
                if not sample_id:
                    dir_path = os.path.relpath(root, parent_dir)
                    sample_id = find_id_in_string(dir_path, ids)
                
                if sample_id:
                    logger.debug("DEBUG: Found ID %s in %s (dir_suffix: %s)", 
                               sample_id, os.path.join(root, filename), dir_suffix)
                    
                    # Initialise structure for this sample_id if needed
                    if sample_id not in candidate_files:
                        candidate_files[sample_id] = {}
                    if dir_suffix not in candidate_files[sample_id]:
                        candidate_files[sample_id][dir_suffix] = {'dir': root}
                    
                    filepath = os.path.abspath(os.path.join(root, filename))
                    if "_R1_" in filename:
                        candidate_files[sample_id][dir_suffix]['R1'] = filepath
                    elif "_R2_" in filename:
                        candidate_files[sample_id][dir_suffix]['R2'] = filepath
                
                processed_files += 1
                pbar.update(1)
        
        pbar.close()
        
        # Now select the appropriate files based on merge mode and directory suffixes
        r1_files = {}
        r2_files = {}
        
        for sample_id, dir_variants in candidate_files.items():
            logger.debug("\nDEBUG: Processing variants for ID: %s", sample_id)
            logger.debug("DEBUG: Available directory types: %s", list(dir_variants.keys()))
            
            selected_variant = None
            
            if use_merge:
                # Merge mode: prefer 'merged', fall back to 'base'
                if 'merged' in dir_variants:
                    selected_variant = 'merged'
                    logger.debug("DEBUG: Selected 'merged' variant")
                elif 'base' in dir_variants:
                    selected_variant = 'base'
                    logger.debug("DEBUG: Selected 'base' variant (no merged available)")
                else:
                    # Only numbered variants available - pick the first one but log a warning
                    selected_variant = list(dir_variants.keys())[0]
                    logger.warning("WARNING: ID %s only has numbered variants, using '%s'", 
                                 sample_id, selected_variant)
            else:
                # Non-merge mode: check for duplicates and warn
                if len(dir_variants) > 1:
                    logger.warning("WARNING: ID %s found in multiple directories: %s", 
                                 sample_id, [v['dir'] for v in dir_variants.values()])
                    logger.warning("         Consider using --merge flag to handle duplicates")
                # Use whichever variant we have (will use last one found if multiple)
                selected_variant = list(dir_variants.keys())[0]
            
            # Extract the R1 and R2 files from the selected variant
            if selected_variant and 'R1' in dir_variants[selected_variant]:
                r1_files[sample_id] = dir_variants[selected_variant]['R1']
                logger.debug("DEBUG: Selected R1: %s", r1_files[sample_id])
            if selected_variant and 'R2' in dir_variants[selected_variant]:
                r2_files[sample_id] = dir_variants[selected_variant]['R2']
                logger.debug("DEBUG: Selected R2: %s", r2_files[sample_id])
    
    logger.debug("\nDEBUG: Found R1 files for IDs: %s", list(r1_files.keys()))
    logger.debug("DEBUG: Number of R1 files found: %d", len(r1_files))
    logger.debug("DEBUG: Found R2 files for IDs: %s", list(r2_files.keys()))
    logger.debug("DEBUG: Number of R2 files found: %d", len(r2_files))
    
    matched_pairs = set(r1_files.keys()) & set(r2_files.keys())
    logger.debug("\nDEBUG: Number of matched pairs found: %d", len(matched_pairs))
    if len(matched_pairs) > 0:
        logger.debug("DEBUG: First few matched pairs: %s", list(matched_pairs)[:5])
    
    for sample_id in set(r1_files.keys()) & set(r2_files.keys()):
        # Get taxonomy values for this sample
        taxonomy_values = tuple(
            taxonomy_dicts[col].get(str(sample_id), '') for col in TAXONOMY_COLUMNS
        )
        results[sample_id] = (
            r1_files[sample_id], 
            r2_files[sample_id],
            *taxonomy_values,
            type_status_dict.get(str(sample_id), '')
        )
    
    logger.debug("\nDEBUG: Final number of matched pairs: %d", len(results))
    return results


def write_results(results, output_prefix, output_dir, logger):
    """Write all samples and type samples to output files."""
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    samples_filename = os.path.join(output_dir, f'{output_prefix}.csv')
    types_filename = os.path.join(output_dir, f'{output_prefix}_types.csv')
    
    type_samples = {}
    
    # Results tuple structure: (r1_path, r2_path, phylum, class, order, family, genus, species, type_status)
    # Identify type samples
    for sample_id, values in results.items():
        type_status = values[-1]  # Last element is type_status
        if isinstance(type_status, str) and 'type' in type_status.lower():
            type_samples[sample_id] = values
    
    # Define output header
    header = ['ID', 'forward', 'reverse'] + TAXONOMY_COLUMNS + ['type_status']
    
    # Write all samples
    with open(samples_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        for sample_id, values in results.items():
            writer.writerow([sample_id] + list(values))
    
    # Write type samples
    with open(types_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        for sample_id, values in type_samples.items():
            writer.writerow([sample_id] + list(values))
    
    return (len(results), len(type_samples), samples_filename, types_filename)


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Process DNA sequencing read files and match with sample metadata',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Using default prefix (directory name)
  python %(prog)s /path/to/reads metadata.csv /output
  
  # Using custom prefix
  python %(prog)s /path/to/reads metadata.csv /output -op samples_batch1
  
  # Using merge mode to prefer merged directories
  python %(prog)s /path/to/reads metadata.csv /output --merge
  
  # Using a custom ID column name
  python %(prog)s /path/to/reads metadata.csv /output --id "sample ID"
        """
    )
    
    parser.add_argument('parent_directory', 
                       help='Directory containing FASTQ/FQ files (directly or in subdirectories)')
    parser.add_argument('sample_metadata_file', 
                       help='CSV file with sample metadata')
    parser.add_argument('output_directory', 
                       help='Directory where output files will be saved')
    parser.add_argument('-op', '--output-prefix', 
                       help='Prefix for output files (default: uses directory name)',
                       default=None)
    parser.add_argument('--merge', 
                       action='store_true',
                       help='Prefer *_merged directories when available, fall back to non-suffixed directories')
    parser.add_argument('--identifier', '--id',
                       dest='id_column',
                       default='Process ID',
                       help='Column name in metadata CSV containing sample IDs (default: "Process ID")')
    
    return parser.parse_args()


def main():
    args = parse_arguments()
    
    parent_dir = args.parent_directory
    sample_metadata_file = args.sample_metadata_file
    output_dir = args.output_directory
    use_merge = args.merge
    id_column = args.id_column
    
    # Determine output prefix
    if args.output_prefix:
        output_prefix = args.output_prefix
        log_filename = f'{output_prefix}.log'
    else:
        dir_name = get_directory_name(parent_dir)
        output_prefix = f'samples_{dir_name}'
        log_filename = f'samples_{dir_name}.log'
    
    # Set up logging
    logger = setup_logging(output_dir, log_filename)
    
    logger.debug("DEBUG: Using output prefix: %s", output_prefix)
    logger.debug("DEBUG: Using ID column: %s", id_column)
    logger.debug("DEBUG: Reading metadata file: %s", sample_metadata_file)
    
    metadata_df = pd.read_csv(sample_metadata_file, low_memory=False)
    logger.debug("DEBUG: Metadata shape: %s", metadata_df.shape)
    
    # Validate that the ID column exists
    if id_column not in metadata_df.columns:
        logger.error("ERROR: Column '%s' not found in metadata file.", id_column)
        logger.error("Available columns: %s", list(metadata_df.columns))
        sys.exit(1)
    
    results = extract_id_reads_taxonomy_type(parent_dir, metadata_df, logger, use_merge, id_column)
    
    (sample_count, type_count,
     samples_file, types_file) = write_results(
         results, output_prefix, output_dir, logger)
    
    logger.debug("\nProcessing complete:")
    logger.debug("- %d samples written to '%s'", sample_count, samples_file)
    logger.debug("- %d type specimens written to '%s'", type_count, types_file)


if __name__ == "__main__":
    # For backward compatibility with old usage pattern
    if len(sys.argv) == 4 and not any(arg.startswith('-') for arg in sys.argv[1:]):
        # Old usage: python script.py parent_dir metadata_file output_dir
        parent_dir = sys.argv[1]
        sample_metadata_file = sys.argv[2]
        output_dir = sys.argv[3]
        id_column = 'Process ID'  # Default for backward compatibility
        
        # Use old behavior with directory name as prefix
        dir_name = get_directory_name(parent_dir)
        output_prefix = f'samples_{dir_name}'
        log_filename = f'samples_{dir_name}.log'
        
        logger = setup_logging(output_dir, log_filename)
        logger.debug("DEBUG: Using backward compatibility mode")
        logger.debug("DEBUG: Using output prefix: %s", output_prefix)
        logger.debug("DEBUG: Using ID column: %s", id_column)
        logger.debug("DEBUG: Reading metadata file: %s", sample_metadata_file)
        
        metadata_df = pd.read_csv(sample_metadata_file, low_memory=False)
        logger.debug("DEBUG: Metadata shape: %s", metadata_df.shape)
        
        # Validate that the ID column exists
        if id_column not in metadata_df.columns:
            logger.error("ERROR: Column '%s' not found in metadata file.", id_column)
            logger.error("Available columns: %s", list(metadata_df.columns))
            sys.exit(1)
        
        results = extract_id_reads_taxonomy_type(parent_dir, metadata_df, logger, use_merge=False, id_column=id_column)
        
        (sample_count, type_count,
         samples_file, types_file) = write_results(
             results, output_prefix, output_dir, logger)
        
        logger.debug("\nProcessing complete:")
        logger.debug("- %d samples written to '%s'", sample_count, samples_file)
        logger.debug("- %d type specimens written to '%s'", type_count, types_file)
    else:
        main()
