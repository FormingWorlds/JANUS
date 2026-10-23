"""Tests for src/janus/utils/data.py and the nightly data cache steps.

Exercises the fwl-io data plumbing without the network: the spectral-file and
stellar-spectra keys and directories, the downloads each entry point asks
fwl-io for, the unknown-name error contract, and the fwl-io cache steps of the nightly.
See docs/How-to/test.md.
"""

import pytest

from janus.utils import data as jdata

pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]


class _Fetcher:
    """Records the datasets the data module asks fwl-io to fetch."""

    fetched = []

    def __init__(self, **kwargs):
        self.kwargs = kwargs
        self.target_dir = kwargs['data_root'] / kwargs['subdir']

    def fetch_all(self):
        self.fetched.append(self.kwargs['subdir'])


@pytest.fixture
def fetches(monkeypatch, tmp_path):
    """Route janus data fetches to a recorder below tmp_path and return the log."""
    import fwl_io

    _Fetcher.fetched = []
    monkeypatch.setattr(jdata, 'FWL_DATA_DIR', tmp_path)
    monkeypatch.setattr(fwl_io, 'create_fetcher', lambda **kw: _Fetcher(**kw))
    return _Fetcher.fetched


@pytest.mark.parametrize(
    ('fname', 'nband', 'expected'),
    [
        ('Dayspring', 48, ['atmos_clim/spectral_files/dayspring/48']),
        ('Frostflow', 4096, ['atmos_clim/spectral_files/frostflow/4096']),
        ('Oak', 256, ['atmos_clim/spectral_files/oak/318']),
        ('Dayspring', None, ['atmos_clim/spectral_files/dayspring/256']),
        (
            '',
            16,
            [
                'atmos_clim/spectral_files/dayspring/256',
                'atmos_clim/spectral_files/frostflow/256',
                'atmos_clim/spectral_files/oak/318',
            ],
        ),
    ],
)
def test_spectral_download_asks_fwl_io_for_the_named_dataset(fetches, fname, nband, expected):
    """A group with several band counts fetches the one asked for; Oak has only 318.

    The empty name fetches the default list whatever nband says, and nothing else.
    """
    jdata.DownloadSpectralFiles(fname=fname, nband=nband)
    assert fetches == expected
    assert len(set(fetches)) == len(fetches)


@pytest.mark.parametrize(
    ('fname', 'nband', 'message'),
    [
        (
            'Dayspring',
            100,
            r"Dayspring has no 100-band spectral file; declared: \['16', '48', '256', '4096'\]",
        ),
        ('NotADataset', 256, "No spectral file group 'NotADataset'"),
    ],
)
def test_unknown_spectral_file_raises_before_any_fetch(fetches, fname, nband, message):
    """A band count or group the manifest does not declare raises its own error and
    fetches nothing."""
    with pytest.raises(ValueError, match=message):
        jdata.DownloadSpectralFiles(fname=fname, nband=nband)
    assert fetches == []


def test_cli_reports_an_unknown_group_without_a_traceback():
    """`janus download spectral` turns the ValueError into a one-line CLI error."""
    from click.testing import CliRunner

    from janus.cli import cli

    result = CliRunner().invoke(cli, ['download', 'spectral', '-n', 'NotADataset'])
    assert result.exit_code == 1
    assert "Error: No spectral file group 'NotADataset'" in result.output
    assert 'Traceback' not in result.output


def test_fetches_carry_the_manifest_pins_of_the_dataset(fetches, monkeypatch, tmp_path):
    """Each fetch gets the Zenodo record, the DataverseNL mirror and the registry of
    its own manifest entry, for a spectral file and for the named spectra."""
    import fwl_io

    seen = []
    monkeypatch.setattr(
        fwl_io, 'create_fetcher', lambda **kw: seen.append(kw) or _Fetcher(**kw)
    )
    jdata.DownloadSpectralFiles('Oak')
    jdata.DownloadStellarSpectra()
    keys = ['atmos_clim.spectral_files.oak.318', jdata.STELLAR_SPECTRA_NAMED]
    for kw, key in zip(seen, keys, strict=True):
        ds = jdata._shared_datasets()[key]
        for field in ('subdir', 'zenodo', 'dataverse', 'extract'):
            assert kw[field] == getattr(ds, field)
        assert kw['registry'] == ds.registry() and kw['data_root'] == tmp_path.absolute()
    assert 'Oak.sf' in seen[0]['registry'] and 'sun.txt' in seen[1]['registry']
    zenodo = ['10.5281/zenodo.15743843', '10.5281/zenodo.23197931']
    assert [kw['zenodo'] for kw in seen] == zenodo
    assert [kw['dataverse'] for kw in seen] == ['10.34894/K3UKBX', '10.34894/URD02X']


def test_a_single_band_group_warns_only_on_an_explicit_other_count(fetches, caplog):
    """Oak resolves to 318; a different count given explicitly is ignored with a
    warning at the default log level, and no count gives no warning."""
    jdata.DownloadSpectralFiles('Oak')
    assert jdata.spectral_file_key('Oak', 318) == jdata.spectral_file_key('Oak', '318')
    assert caplog.text == ''
    assert jdata.spectral_file_key('Oak', 4096).endswith('oak.318')
    assert 'Oak has only 318 bands; ignoring the requested 4096' in caplog.text


def test_a_band_count_without_a_group_is_ignored_with_a_warning(fetches, caplog):
    """A count given with no group still fetches the default list, and says so."""
    jdata.DownloadSpectralFiles(nband=4096)
    assert 'ignoring the band count 4096' in caplog.text
    assert fetches == [
        'atmos_clim/spectral_files/dayspring/256',
        'atmos_clim/spectral_files/frostflow/256',
        'atmos_clim/spectral_files/oak/318',
    ]


def test_a_multi_band_group_is_located_where_it_was_downloaded(fetches, tmp_path):
    """Without a count, download and lookup both use DEFAULT_BANDS for Dayspring."""
    jdata.DownloadSpectralFiles('Dayspring')
    assert fetches == ['atmos_clim/spectral_files/dayspring/256']
    assert jdata.spectral_file_dir('Dayspring') == tmp_path / fetches[0]
    assert jdata.DEFAULT_BANDS == 256


def test_a_key_with_a_non_numeric_suffix_does_not_break_its_group(monkeypatch):
    """A manifest key under a group prefix that is not a band count is ignored."""
    declared = dict(jdata._shared_datasets())
    declared['atmos_clim.spectral_files.oak.318_k'] = declared[
        'atmos_clim.spectral_files.oak.318'
    ]
    monkeypatch.setattr(jdata, '_shared_datasets', lambda: declared)
    assert jdata.spectral_file_key('Oak') == 'atmos_clim.spectral_files.oak.318'
    assert jdata.spectral_file_key('Dayspring', 48).endswith('dayspring.48')


def test_stellar_spectra_and_directories_resolve_through_fwl_io(fetches, tmp_path):
    """The named spectra are fetched by key, and the helpers return the dataset
    directories that hold sun.txt and Oak.sf."""
    jdata.DownloadStellarSpectra()
    assert fetches == ['star/spectra/named']
    assert jdata.stellar_spectra_dir() == tmp_path / 'star/spectra/named'
    assert jdata.spectral_file_dir('Oak') == tmp_path / 'atmos_clim/spectral_files/oak/318'
    assert jdata.GetFWLData() == tmp_path.absolute()


def test_directories_are_the_version_directories_fwl_io_fills(monkeypatch, tmp_path):
    """With the real fwl-io fetcher, the helpers return the r<record-id> directory and
    create the data root, but not the dataset directory."""
    from fwl_io.manifest import zenodo_record_id

    root = tmp_path / 'absent_root'
    monkeypatch.setattr(jdata, 'FWL_DATA_DIR', root)
    oak = jdata.spectral_file_dir('Oak')
    named = jdata.stellar_spectra_dir()
    records = {k: zenodo_record_id(ds.zenodo) for k, ds in jdata._shared_datasets().items()}
    oak_record = records['atmos_clim.spectral_files.oak.318']
    assert oak == root / f'atmos_clim/spectral_files/oak/318/r{oak_record}'
    assert named == root / f'star/spectra/named/r{records[jdata.STELLAR_SPECTRA_NAMED]}'
    assert root.is_dir() and not oak.exists()


def test_selectable_groups_are_the_manifest_spectral_sets(fetches, monkeypatch):
    """A group the manifest declares is accepted, and a name it does not declare is
    refused before any fetch, whatever the installed fwl-io declares."""
    keys = ['atmos_clim.spectral_files.honeyside.48', 'atmos_clim.spectral_files.honeyside.256']
    keys.append('atmos_clim.spectral_files.oak.318')
    monkeypatch.setattr(jdata, '_shared_datasets', lambda: dict.fromkeys(keys))
    assert jdata.spectral_file_key('Honeyside', 48) == keys[0]
    assert jdata.spectral_file_key('Oak') == keys[2]
    with pytest.raises(ValueError, match="No spectral file group 'Legacy'"):
        jdata.DownloadSpectralFiles('Legacy')
    assert fetches == []
