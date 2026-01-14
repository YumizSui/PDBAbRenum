"""Tests for renumbering functionality."""

import os
import tempfile
from pathlib import Path

import pytest
from pdbabrenum import renumber, RenumberingError, SUPPORTED_SCHEMES


# Get path to examples directory
EXAMPLES_DIR = Path(__file__).parent.parent / "examples"
TEST_PDB_1A0Q = EXAMPLES_DIR / "1a0q.pdb"


@pytest.fixture
def temp_output():
    """Create a temporary output file."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        temp_path = f.name
    yield temp_path
    # Cleanup
    if os.path.exists(temp_path):
        os.remove(temp_path)


class TestRenumbering:
    """Test renumbering functionality."""

    def test_supported_schemes(self):
        """Test that all expected schemes are supported."""
        assert 'imgt' in SUPPORTED_SCHEMES
        assert 'chothia' in SUPPORTED_SCHEMES
        assert 'kabat' in SUPPORTED_SCHEMES
        assert 'aho' in SUPPORTED_SCHEMES

    @pytest.mark.skipif(not TEST_PDB_1A0Q.exists(), reason="Example PDB not found")
    def test_basic_renumbering_imgt(self, temp_output):
        """Test basic renumbering with IMGT scheme."""
        result = renumber(
            in_pdb=str(TEST_PDB_1A0Q),
            out_pdb=temp_output,
            scheme='imgt',
            verbose=False
        )

        assert os.path.exists(temp_output)
        assert 'heavy' in result
        assert 'light' in result
        assert len(result['heavy']) > 0 or len(result['light']) > 0

    @pytest.mark.skipif(not TEST_PDB_1A0Q.exists(), reason="Example PDB not found")
    def test_renumbering_chothia(self, temp_output):
        """Test renumbering with Chothia scheme."""
        result = renumber(
            in_pdb=str(TEST_PDB_1A0Q),
            out_pdb=temp_output,
            scheme='chothia',
            verbose=False
        )

        assert os.path.exists(temp_output)
        assert 'heavy' in result
        assert 'light' in result

    @pytest.mark.skipif(not TEST_PDB_1A0Q.exists(), reason="Example PDB not found")
    def test_renumbering_with_chain_specification(self, temp_output):
        """Test renumbering with specified chains."""
        result = renumber(
            in_pdb=str(TEST_PDB_1A0Q),
            out_pdb=temp_output,
            scheme='imgt',
            heavy_chains=['H'],
            light_chains=['L'],
            verbose=False
        )

        assert os.path.exists(temp_output)

    @pytest.mark.skipif(not TEST_PDB_1A0Q.exists(), reason="Example PDB not found")
    def test_renumbering_ab_only(self, temp_output):
        """Test renumbering with ab_only option."""
        result = renumber(
            in_pdb=str(TEST_PDB_1A0Q),
            out_pdb=temp_output,
            scheme='imgt',
            ab_only=True,
            verbose=False
        )

        assert os.path.exists(temp_output)

    @pytest.mark.skipif(not TEST_PDB_1A0Q.exists(), reason="Example PDB not found")
    def test_renumbering_with_fasta_output(self, temp_output):
        """Test renumbering with FASTA output."""
        fasta_output = temp_output.replace('.pdb', '.fasta')

        result = renumber(
            in_pdb=str(TEST_PDB_1A0Q),
            out_pdb=temp_output,
            scheme='imgt',
            output_fasta=fasta_output,
            verbose=False
        )

        assert os.path.exists(temp_output)
        assert os.path.exists(fasta_output)

        # Cleanup
        if os.path.exists(fasta_output):
            os.remove(fasta_output)

    @pytest.mark.skipif(not TEST_PDB_1A0Q.exists(), reason="Example PDB not found")
    def test_renumbering_with_abnumber_output(self, temp_output):
        """Test renumbering with abnumber results output."""
        abnumber_output = temp_output.replace('.pdb', '.json')

        result = renumber(
            in_pdb=str(TEST_PDB_1A0Q),
            out_pdb=temp_output,
            scheme='imgt',
            output_abnumber=abnumber_output,
            verbose=False
        )

        assert os.path.exists(temp_output)
        assert os.path.exists(abnumber_output)

        # Cleanup
        if os.path.exists(abnumber_output):
            os.remove(abnumber_output)

    def test_invalid_scheme(self, temp_output):
        """Test that invalid scheme raises error."""
        with pytest.raises(RenumberingError):
            renumber(
                in_pdb=str(TEST_PDB_1A0Q),
                out_pdb=temp_output,
                scheme='invalid_scheme',
                verbose=False
            )

    def test_invalid_pdb_file(self, temp_output):
        """Test that invalid PDB file raises error."""
        with pytest.raises(RenumberingError):
            renumber(
                in_pdb='/nonexistent/file.pdb',
                out_pdb=temp_output,
                scheme='imgt',
                verbose=False
            )
