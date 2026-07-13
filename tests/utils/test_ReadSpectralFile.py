"""Tests for src/janus/utils/ReadSpectralFile.py.

Exercises the spectral-file band-edge parser: metre-to-nanometre conversion
on a synthetic band block, edge continuity, the no-band-block fallback, and
the missing-file error contract. See docs/How-to/test.md.
"""

import numpy as np
import pytest

from janus.utils.ReadSpectralFile import ReadBandEdges

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]

BAND_BLOCK = """*BLOCK: TYPE =    1
Band        Lower limit         Upper limit
    1        1.000000000E-06     2.000000000E-06
    2        2.000000000E-06     4.000000000E-06
    3        4.000000000E-06     1.000000000E-05
*END
trailing content that must not be parsed
"""


def test_band_edges_parsed_and_converted_to_nm(tmp_path):
    """Three bands give four contiguous edges converted from m to nm.

    The parser stores the lower limit of the first band and the upper limit
    of every band; a dropped 1e9 conversion factor would leave values nine
    orders of magnitude too small, and any discontinuity between bands
    would break the shared-edge structure asserted here.
    """
    sf = tmp_path / 'toy.sf'
    sf.write_text(BAND_BLOCK)
    edges = ReadBandEdges(str(sf))

    np.testing.assert_allclose(edges, [1000.0, 2000.0, 4000.0, 10000.0], rtol=1e-12)
    # Unit guard: edges in nm for micron-scale bands are O(1e3), not O(1e-6).
    assert min(edges) > 1.0
    # Monotone increasing edges: bands tile the interval without overlap.
    assert all(b > a for a, b in zip(edges[:-1], edges[1:]))


def test_file_without_band_block_and_missing_file(tmp_path):
    """A file lacking the band header yields no edges; a missing path raises.

    The header sentinel is the only trigger for parsing, so an arbitrary
    text file must return an empty list rather than misread its lines; the
    missing-file branch is the documented error contract.
    """
    plain = tmp_path / 'plain.sf'
    plain.write_text('no bands here\njust text\n*END\n')
    assert ReadBandEdges(str(plain)) == []

    with pytest.raises(Exception, match='spectral file'):
        ReadBandEdges(str(tmp_path / 'absent.sf'))
