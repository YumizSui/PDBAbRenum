# PDBAbRenum

Renumber antibody structures in PDB files using IMGT and other numbering schemes.

Extracts and renumbers only the variable regions (VH/VL) of antibody chains.

## Installation

```bash
# Using pixi
git clone https://github.com/Yumizsui/PDBAbRenum.git
cd PDBAbRenum
pixi install && pixi run install

# Using pip
pip install git+https://github.com/Yumizsui/PDBAbRenum.git
```

## Usage

```bash
# Basic
pdbabrenum input.pdb output.pdb --scheme imgt

# Specify chains
pdbabrenum input.pdb output.pdb --heavy H --light L

# Extract antibody chains only
pdbabrenum input.pdb output.pdb --ab-only

# Additional outputs
pdbabrenum input.pdb output.pdb --abnumber results.json --fasta sequences.fasta

# Batch processing
pdbabrenum --batch "*.pdb" --output-dir renumbered/ --scheme chothia
```

## Python API

```python
from pdbabrenum import renumber

result = renumber(
    in_pdb='input.pdb',
    out_pdb='output.pdb',
    scheme='imgt',
    output_abnumber='analysis.json',
    output_fasta='sequences.fasta'
)
```
