"""Core renumbering functionality for antibody PDB structures."""

import warnings
from typing import List, Tuple, Optional, Dict
from dataclasses import dataclass

import abnumber
from Bio import PDB
from Bio.PDB import Model, Chain, Residue, Selection
from Bio.PDB.Polypeptide import protein_letters_3to1

from .utils import AbnumberResult, extract_abnumber_info, save_abnumber_results, save_fasta


# Supported numbering schemes
SUPPORTED_SCHEMES = ['imgt', 'chothia', 'kabat', 'aho']


@dataclass
class ChainInfo:
    """Information about a renumbered chain."""
    chain_id: str
    chain_type: str  # 'H', 'L', 'K', or None
    scheme: str
    num_residues: int
    sequence: str
    cdr_regions: Optional[Dict[str, List[Tuple[int, str]]]] = None


class RenumberingError(Exception):
    """Exception raised during renumbering."""
    pass


def biopython_chain_to_sequence(chain: Chain.Chain) -> Tuple[str, List[Residue.Residue]]:
    """
    Convert a BioPython chain to its amino acid sequence.

    Args:
        chain: BioPython chain object

    Returns:
        Tuple of (sequence string, list of residues)
    """
    residue_list = Selection.unfold_entities(chain, 'R')
    seq = ''.join([
        protein_letters_3to1.get(r.resname, 'X')
        for r in residue_list
    ])
    return seq, residue_list


def assign_number_to_sequence(seq: str, scheme: str = 'imgt') -> Tuple[List, abnumber.Chain]:
    """
    Assign numbering to an antibody sequence.

    Args:
        seq: Amino acid sequence
        scheme: Numbering scheme ('imgt', 'chothia', 'kabat', 'aho')

    Returns:
        Tuple of (numbering list, abnumber Chain object)

    Raises:
        RenumberingError: If sequence cannot be numbered
    """
    if scheme not in SUPPORTED_SCHEMES:
        raise RenumberingError(
            f"Unsupported scheme '{scheme}'. Supported schemes: {SUPPORTED_SCHEMES}"
        )

    try:
        abchain = abnumber.Chain(seq, scheme=scheme)
    except abnumber.ChainParseError as e:
        raise RenumberingError(f"Failed to parse sequence: {str(e)}")

    offset = seq.index(abchain.seq)
    if offset < 0:
        raise RenumberingError(
            'The identified Fv sequence is not a subsequence of the original sequence.'
        )

    numbers = [None for _ in range(len(seq))]
    for i, (pos, aa) in enumerate(abchain):
        resseq = pos.number
        icode = pos.letter if pos.letter else ' '
        numbers[i + offset] = (resseq, icode)

    return numbers, abchain


def renumber_biopython_chain(
    chain_id: str,
    residue_list: List[Residue.Residue],
    numbers: List[Optional[Tuple[int, str]]]
) -> Chain.Chain:
    """
    Create a renumbered BioPython chain.

    Args:
        chain_id: New chain identifier
        residue_list: List of residues
        numbers: List of (resseq, icode) tuples

    Returns:
        Renumbered BioPython chain
    """
    chain = Chain.Chain(chain_id)
    for residue, number in zip(residue_list, numbers):
        if number is None:
            continue
        residue = residue.copy()
        new_id = (residue.id[0], number[0], number[1])
        residue.id = new_id
        chain.add(residue)
    return chain


def get_cdr_regions(abchain: abnumber.Chain) -> Dict[str, List[Tuple[int, str]]]:
    """
    Extract CDR region positions from an abnumber Chain.

    Args:
        abchain: abnumber Chain object

    Returns:
        Dictionary with CDR regions (cdr1, cdr2, cdr3) and their positions
    """
    cdr_regions = {}

    try:
        for cdr_name in ['cdr1', 'cdr2', 'cdr3']:
            cdr_seq = getattr(abchain, f'{cdr_name}_seq', None)
            if cdr_seq:
                cdr_positions = []
                for pos, aa in abchain:
                    # Check if position is in this CDR
                    if hasattr(abchain, f'is_{cdr_name}'):
                        # This is a simplified version - abnumber provides CDR definitions
                        pass
                cdr_regions[cdr_name] = cdr_seq
    except AttributeError:
        # CDR information not available
        pass

    return cdr_regions


def renumber(
    in_pdb: str,
    out_pdb: str,
    scheme: str = 'imgt',
    heavy_chains: Optional[List[str]] = None,
    light_chains: Optional[List[str]] = None,
    ab_only: bool = False,
    verbose: bool = False,
    output_abnumber: Optional[str] = None,
    output_fasta: Optional[str] = None
) -> Dict[str, any]:
    """
    Renumber antibody chains in a PDB file.

    Args:
        in_pdb: Input PDB file path
        out_pdb: Output PDB file path
        scheme: Numbering scheme ('imgt', 'chothia', 'kabat', 'aho')
        heavy_chains: List of chain IDs to treat as heavy chains (auto-detect if None)
        light_chains: List of chain IDs to treat as light chains (auto-detect if None)
        ab_only: If True, only output antibody chains
        verbose: If True, print detailed information
        output_abnumber: Path to save abnumber analysis results (JSON/TSV)
        output_fasta: Path to save sequences in FASTA format

    Returns:
        Dictionary with chain information and results

    Raises:
        RenumberingError: If renumbering fails
    """
    parser = PDB.PDBParser(QUIET=not verbose)

    try:
        structure = parser.get_structure(None, in_pdb)
    except Exception as e:
        raise RenumberingError(f"Failed to parse PDB file: {str(e)}")

    model = structure[0]
    model_new = Model.Model(0)

    result = {
        'heavy': [],
        'light': [],
        'other': []
    }

    chain_infos = []
    abnumber_results = []
    sequences = {}

    for chain in model:
        chain_id = chain.id

        # Check if this chain should be processed
        is_specified_heavy = heavy_chains and chain_id in heavy_chains
        is_specified_light = light_chains and chain_id in light_chains

        try:
            seq, reslist = biopython_chain_to_sequence(chain)
            numbers, abchain = assign_number_to_sequence(seq, scheme=scheme)
            chain_new = renumber_biopython_chain(chain_id, reslist, numbers)

            # Determine chain type
            if is_specified_heavy:
                chain_type = 'H'
            elif is_specified_light:
                chain_type = abchain.chain_type if abchain.chain_type in ('K', 'L') else 'L'
            else:
                chain_type = abchain.chain_type

            # Categorize chain
            if chain_type == 'H':
                result['heavy'].append(chain_id)
            elif chain_type in ('K', 'L'):
                result['light'].append(chain_id)
            else:
                result['other'].append(chain_id)

            # Create chain info
            chain_info = ChainInfo(
                chain_id=chain_id,
                chain_type=chain_type,
                scheme=scheme,
                num_residues=len([n for n in numbers if n is not None]),
                sequence=seq
            )
            chain_infos.append(chain_info)

            # Store abnumber result
            if output_abnumber or output_fasta:
                abnumber_result = extract_abnumber_info(chain_id, abchain, scheme)
                abnumber_results.append(abnumber_result)
                sequences[chain_id] = str(abchain.seq)

            if verbose:
                print(f'[INFO] Renumbered chain {chain_id} ({chain_type}) using {scheme.upper()} scheme: {chain_info.num_residues} residues')

            model_new.add(chain_new)

        except (abnumber.ChainParseError, RenumberingError) as e:
            if verbose:
                print(f'[INFO] Chain {chain_id} does not contain valid antibody Fv: {str(e)}')

            # Keep non-antibody chain if not ab_only
            if not ab_only:
                chain_new = chain.copy()
                model_new.add(chain_new)
                result['other'].append(chain_id)

    # Save the renumbered structure
    pdb_io = PDB.PDBIO()
    pdb_io.set_structure(model_new)

    # Add header comment
    class SelectWithHeader(PDB.Select):
        def __init__(self, header_lines):
            self.header_lines = header_lines
            self.header_written = False

    try:
        pdb_io.save(out_pdb)

        # Add REMARK header
        with open(out_pdb, 'r') as f:
            content = f.read()

        header = (
            f"REMARK   1 RENUMBERED BY PDBABRENUM\n"
            f"REMARK   1 NUMBERING SCHEME: {scheme.upper()}\n"
            f"REMARK   1 HEAVY CHAINS: {', '.join(result['heavy']) if result['heavy'] else 'NONE'}\n"
            f"REMARK   1 LIGHT CHAINS: {', '.join(result['light']) if result['light'] else 'NONE'}\n"
        )

        with open(out_pdb, 'w') as f:
            f.write(header)
            f.write(content)

    except Exception as e:
        raise RenumberingError(f"Failed to save PDB file: {str(e)}")

    # Save abnumber results if requested
    if output_abnumber and abnumber_results:
        format = 'json' if output_abnumber.endswith('.json') else 'tsv'
        save_abnumber_results(abnumber_results, output_abnumber, format=format)
        if verbose:
            print(f'[INFO] Saved abnumber results to {output_abnumber}')

    # Save FASTA if requested
    if output_fasta and sequences:
        save_fasta(sequences, output_fasta, description=f"{scheme.upper()} numbered")
        if verbose:
            print(f'[INFO] Saved sequences to {output_fasta}')

    result['chain_infos'] = chain_infos
    result['abnumber_results'] = abnumber_results if abnumber_results else None
    result['sequences'] = sequences if sequences else None

    return result
