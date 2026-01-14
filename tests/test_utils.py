"""Tests for utility functions."""

import os
import tempfile
from pathlib import Path

import pytest
from pdbabrenum.utils import (
    generate_output_filename,
    find_pdb_files,
    save_fasta,
)


class TestUtils:
    """Test utility functions."""

    def test_generate_output_filename(self):
        """Test output filename generation."""
        input_file = "/path/to/input.pdb"
        output = generate_output_filename(
            input_file,
            suffix="_renum",
            extension=".pdb"
        )

        assert output.endswith("input_renum.pdb")
        assert "/path/to/" in output

    def test_generate_output_filename_with_dir(self):
        """Test output filename generation with custom directory."""
        input_file = "/path/to/input.pdb"
        output = generate_output_filename(
            input_file,
            output_dir="/custom/dir",
            suffix="_numbered",
            extension=".pdb"
        )

        assert output == "/custom/dir/input_numbered.pdb"

    def test_find_pdb_files(self):
        """Test finding PDB files in directory."""
        # Create temporary directory with test files
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create some test files
            Path(tmpdir, "test1.pdb").touch()
            Path(tmpdir, "test2.pdb").touch()
            Path(tmpdir, "test.txt").touch()

            pdb_files = find_pdb_files(tmpdir)

            assert len(pdb_files) == 2
            assert all(f.suffix == '.pdb' for f in pdb_files)

    def test_save_fasta(self):
        """Test saving sequences to FASTA."""
        sequences = {
            'H': 'EVQLQQSGAELVKPGASVKM',
            'L': 'DIVMTQSPSSLSASVGDRVT'
        }

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            fasta_path = f.name

        try:
            save_fasta(sequences, fasta_path, description="Test antibody")

            assert os.path.exists(fasta_path)

            # Read and verify
            with open(fasta_path, 'r') as f:
                content = f.read()
                assert '>H' in content
                assert '>L' in content
                assert 'EVQLQQSGAELVKPGASVKM' in content
                assert 'DIVMTQSPSSLSASVGDRVT' in content

        finally:
            if os.path.exists(fasta_path):
                os.remove(fasta_path)
