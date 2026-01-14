"""Command-line interface for PDBAbRenum."""

import argparse
import sys
from pathlib import Path
from typing import Optional

from .renumber import renumber, RenumberingError, SUPPORTED_SCHEMES
from .utils import find_pdb_files, generate_output_filename


def process_single_file(
    input_file: str,
    output_file: Optional[str],
    scheme: str,
    heavy_chains: Optional[list],
    light_chains: Optional[list],
    ab_only: bool,
    verbose: bool,
    output_abnumber: Optional[str],
    output_fasta: Optional[str]
) -> bool:
    """
    Process a single PDB file.

    Returns:
        True if successful, False otherwise
    """
    try:
        if verbose:
            print(f"\n{'='*60}")
            print(f"Processing: {input_file}")
            print(f"{'='*60}")

        # Generate output filename if not provided
        if output_file is None:
            output_file = generate_output_filename(
                input_file,
                suffix=f"_{scheme}",
                extension=".pdb"
            )

        result = renumber(
            in_pdb=input_file,
            out_pdb=output_file,
            scheme=scheme,
            heavy_chains=heavy_chains,
            light_chains=light_chains,
            ab_only=ab_only,
            verbose=verbose,
            output_abnumber=output_abnumber,
            output_fasta=output_fasta
        )

        print(f"\n[SUCCESS] Renumbered PDB saved to: {output_file}")
        print(f"  Heavy chains: {', '.join(result['heavy']) if result['heavy'] else 'None'}")
        print(f"  Light chains: {', '.join(result['light']) if result['light'] else 'None'}")
        if result['other']:
            print(f"  Other chains: {', '.join(result['other'])}")

        return True

    except RenumberingError as e:
        print(f"[ERROR] Failed to process {input_file}: {str(e)}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"[ERROR] Unexpected error processing {input_file}: {str(e)}", file=sys.stderr)
        if verbose:
            import traceback
            traceback.print_exc()
        return False


def process_batch(
    input_pattern: str,
    output_dir: Optional[str],
    scheme: str,
    heavy_chains: Optional[list],
    light_chains: Optional[list],
    ab_only: bool,
    verbose: bool,
    output_abnumber_dir: Optional[str],
    output_fasta_dir: Optional[str]
) -> None:
    """Process multiple PDB files matching a pattern."""

    # Find input files
    input_path = Path(input_pattern)

    if '*' in input_pattern or '?' in input_pattern:
        # Glob pattern
        parent_dir = input_path.parent if input_path.parent.exists() else Path.cwd()
        pattern = input_path.name
        pdb_files = find_pdb_files(str(parent_dir), pattern)
    elif input_path.is_dir():
        # Directory
        pdb_files = find_pdb_files(input_pattern)
    else:
        print(f"[ERROR] Invalid batch input: {input_pattern}", file=sys.stderr)
        sys.exit(1)

    if not pdb_files:
        print(f"[WARNING] No PDB files found matching: {input_pattern}")
        return

    print(f"\n[INFO] Found {len(pdb_files)} PDB file(s) to process")

    # Create output directory if specified
    if output_dir:
        Path(output_dir).mkdir(parents=True, exist_ok=True)

    if output_abnumber_dir:
        Path(output_abnumber_dir).mkdir(parents=True, exist_ok=True)

    if output_fasta_dir:
        Path(output_fasta_dir).mkdir(parents=True, exist_ok=True)

    # Process each file
    success_count = 0
    for pdb_file in pdb_files:
        output_file = generate_output_filename(
            str(pdb_file),
            output_dir=output_dir,
            suffix=f"_{scheme}",
            extension=".pdb"
        )

        # Generate output paths for abnumber and FASTA
        abnumber_file = None
        fasta_file = None

        if output_abnumber_dir:
            abnumber_file = generate_output_filename(
                str(pdb_file),
                output_dir=output_abnumber_dir,
                suffix=f"_{scheme}",
                extension=".json"
            )

        if output_fasta_dir:
            fasta_file = generate_output_filename(
                str(pdb_file),
                output_dir=output_fasta_dir,
                suffix=f"_{scheme}",
                extension=".fasta"
            )

        if process_single_file(
            str(pdb_file),
            output_file,
            scheme,
            heavy_chains,
            light_chains,
            ab_only,
            verbose,
            abnumber_file,
            fasta_file
        ):
            success_count += 1

    print(f"\n{'='*60}")
    print(f"[SUMMARY] Processed {success_count}/{len(pdb_files)} files successfully")
    print(f"{'='*60}")


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description='Renumber antibody structures in PDB files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with IMGT numbering
  pdbabrenum input.pdb output.pdb --scheme imgt

  # Specify heavy and light chains
  pdbabrenum input.pdb output.pdb --heavy H --light L

  # Extract only antibody chains
  pdbabrenum input.pdb output.pdb --ab-only

  # Batch processing
  pdbabrenum --batch "structures/*.pdb" --output-dir renumbered/

  # Generate additional outputs
  pdbabrenum input.pdb output.pdb --abnumber results.json --fasta sequences.fasta

  # Batch with all outputs
  pdbabrenum --batch "*.pdb" --output-dir out/ --abnumber-dir analysis/ --fasta-dir seqs/
        """
    )

    # Input/Output
    parser.add_argument(
        'input',
        nargs='?',
        help='Input PDB file (or pattern for batch mode)'
    )
    parser.add_argument(
        'output',
        nargs='?',
        help='Output PDB file (not used in batch mode)'
    )

    # Numbering scheme
    parser.add_argument(
        '--scheme', '-s',
        choices=SUPPORTED_SCHEMES,
        default='imgt',
        help='Numbering scheme to use (default: imgt)'
    )

    # Chain specification
    parser.add_argument(
        '--heavy', '-H',
        nargs='+',
        metavar='CHAIN_ID',
        help='Chain ID(s) to treat as heavy chains'
    )
    parser.add_argument(
        '--light', '-L',
        nargs='+',
        metavar='CHAIN_ID',
        help='Chain ID(s) to treat as light chains'
    )

    # Options
    parser.add_argument(
        '--ab-only',
        action='store_true',
        help='Only output antibody chains (exclude other chains)'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Print detailed information'
    )

    # Additional outputs
    parser.add_argument(
        '--abnumber',
        metavar='FILE',
        help='Save abnumber analysis to JSON/TSV file (.json or .tsv extension)'
    )
    parser.add_argument(
        '--fasta',
        metavar='FILE',
        help='Save sequences to FASTA file'
    )

    # Batch processing
    parser.add_argument(
        '--batch', '-b',
        metavar='PATTERN',
        help='Process multiple PDB files (glob pattern or directory)'
    )
    parser.add_argument(
        '--output-dir', '-o',
        metavar='DIR',
        help='Output directory for batch processing'
    )
    parser.add_argument(
        '--abnumber-dir',
        metavar='DIR',
        help='Output directory for abnumber results in batch mode'
    )
    parser.add_argument(
        '--fasta-dir',
        metavar='DIR',
        help='Output directory for FASTA files in batch mode'
    )

    # Version
    parser.add_argument(
        '--version',
        action='version',
        version='PDBAbRenum 0.1.0'
    )

    args = parser.parse_args()

    # Validate arguments
    if args.batch:
        # Batch mode
        process_batch(
            input_pattern=args.batch,
            output_dir=args.output_dir,
            scheme=args.scheme,
            heavy_chains=args.heavy,
            light_chains=args.light,
            ab_only=args.ab_only,
            verbose=args.verbose,
            output_abnumber_dir=args.abnumber_dir,
            output_fasta_dir=args.fasta_dir
        )
    else:
        # Single file mode
        if not args.input:
            parser.error("Input PDB file is required (or use --batch for batch processing)")

        success = process_single_file(
            input_file=args.input,
            output_file=args.output,
            scheme=args.scheme,
            heavy_chains=args.heavy,
            light_chains=args.light,
            ab_only=args.ab_only,
            verbose=args.verbose,
            output_abnumber=args.abnumber,
            output_fasta=args.fasta
        )

        sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
