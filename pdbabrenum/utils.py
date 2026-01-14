"""Utility functions for PDBAbRenum."""

import json
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass, asdict

import abnumber
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


@dataclass
class AbnumberResult:
    """Result from abnumber analysis."""
    chain_id: str
    chain_type: str
    scheme: str
    sequence: str
    numbered_sequence: List[Tuple[str, str]]  # [(position, aa), ...]
    cdr1_seq: Optional[str] = None
    cdr2_seq: Optional[str] = None
    cdr3_seq: Optional[str] = None


def extract_abnumber_info(
    chain_id: str,
    abchain: abnumber.Chain,
    scheme: str
) -> AbnumberResult:
    """
    Extract detailed information from abnumber Chain object.

    Args:
        chain_id: Chain identifier
        abchain: abnumber Chain object
        scheme: Numbering scheme used

    Returns:
        AbnumberResult object with detailed information
    """
    # Get numbered sequence
    numbered_seq = []
    for pos, aa in abchain:
        pos_str = f"{pos.number}{pos.letter}" if pos.letter else str(pos.number)
        numbered_seq.append((pos_str, aa))

    # Extract CDR sequences
    cdr1 = abchain.cdr1_seq if hasattr(abchain, 'cdr1_seq') else None
    cdr2 = abchain.cdr2_seq if hasattr(abchain, 'cdr2_seq') else None
    cdr3 = abchain.cdr3_seq if hasattr(abchain, 'cdr3_seq') else None

    return AbnumberResult(
        chain_id=chain_id,
        chain_type=abchain.chain_type,
        scheme=scheme,
        sequence=str(abchain.seq),
        numbered_sequence=numbered_seq,
        cdr1_seq=cdr1,
        cdr2_seq=cdr2,
        cdr3_seq=cdr3
    )


def save_abnumber_results(
    results: List[AbnumberResult],
    output_file: str,
    format: str = 'json'
) -> None:
    """
    Save abnumber results to file.

    Args:
        results: List of AbnumberResult objects
        output_file: Output file path
        format: Output format ('json' or 'tsv')
    """
    output_path = Path(output_file)

    if format == 'json':
        # Convert to JSON
        data = [asdict(r) for r in results]
        with open(output_path, 'w') as f:
            json.dump(data, f, indent=2)

    elif format == 'tsv':
        # Convert to TSV
        with open(output_path, 'w') as f:
            # Header
            f.write("chain_id\tchain_type\tscheme\tsequence\tcdr1\tcdr2\tcdr3\n")
            # Data
            for r in results:
                f.write(
                    f"{r.chain_id}\t{r.chain_type}\t{r.scheme}\t{r.sequence}\t"
                    f"{r.cdr1_seq or 'N/A'}\t{r.cdr2_seq or 'N/A'}\t{r.cdr3_seq or 'N/A'}\n"
                )

    else:
        raise ValueError(f"Unsupported format: {format}")


def save_fasta(
    sequences: Dict[str, str],
    output_file: str,
    description: str = ""
) -> None:
    """
    Save sequences to FASTA file.

    Args:
        sequences: Dictionary of {chain_id: sequence}
        output_file: Output FASTA file path
        description: Optional description for sequences
    """
    records = []
    for chain_id, seq in sequences.items():
        desc = f"{description} chain {chain_id}" if description else f"chain {chain_id}"
        record = SeqRecord(
            Seq(seq),
            id=chain_id,
            description=desc
        )
        records.append(record)

    with open(output_file, 'w') as f:
        SeqIO.write(records, f, 'fasta')


def find_pdb_files(directory: str, pattern: str = "*.pdb") -> List[Path]:
    """
    Find PDB files in a directory.

    Args:
        directory: Directory to search
        pattern: File pattern (default: "*.pdb")

    Returns:
        List of Path objects
    """
    dir_path = Path(directory)
    if not dir_path.is_dir():
        raise ValueError(f"Not a directory: {directory}")

    return sorted(dir_path.glob(pattern))


def generate_output_filename(
    input_file: str,
    output_dir: Optional[str] = None,
    suffix: str = "_renum",
    extension: str = ".pdb"
) -> str:
    """
    Generate output filename from input filename.

    Args:
        input_file: Input file path
        output_dir: Output directory (use input dir if None)
        suffix: Suffix to add before extension
        extension: File extension

    Returns:
        Output file path
    """
    input_path = Path(input_file)
    stem = input_path.stem

    if output_dir:
        output_path = Path(output_dir) / f"{stem}{suffix}{extension}"
    else:
        output_path = input_path.parent / f"{stem}{suffix}{extension}"

    return str(output_path)
