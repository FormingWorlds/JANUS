# Testing

JANUS uses [pytest](https://docs.pytest.org/) for its test suite. The tests
verify physical correctness by comparing model output against known reference
values, and check utility functions and constants.

## Prerequisites

Install the test dependencies:

```bash
pip install pytest
pip install -r examples/requirements.txt
```

Make sure `RAD_DIR` and `FWL_DATA` are set before running:

```bash
echo $RAD_DIR    # should point to your SOCRATES installation
echo $FWL_DATA   # should point to your data directory
```

## Running the tests

From the root of the JANUS repository:

```bash
pytest tests/
```

To run a specific test file:

```bash
pytest tests/test_runaway_greenhouse.py
pytest tests/test_instellation.py
pytest tests/test_constants.py
pytest tests/test_code.py
```

To run with verbose output:

```bash
pytest tests/ -v
```

## What the tests check

### `test_runaway_greenhouse.py`

Runs the full JANUS pipeline at three surface temperatures for a pure H₂O
atmosphere at 0.3 AU and checks the OLR against reference values to a relative
tolerance of 1 × 10⁻⁵:

| T_surf (K) | Expected OLR (W m⁻²) | Regime |
|---|---|---|
| 200 | 90.73 | Sub-runaway |
| 1705 | 278.88 | Simpson–Nakajima plateau |
| 2800 | 6581.74 | Post-runaway |

Uses `RadConvEqm()` with `cp_dry=False`, `trppD=False`, `rscatter=False`.

### `test_instellation.py`

Runs `MCPA_CBL()` — the coupled boundary layer solver — at two orbital
distances (0.3 AU and 1.4 AU) using `config_instellation.toml` and checks
five output quantities against reference values:

```
(SW_flux_down[0], LW_flux_up[0], net_flux[0], ts, trppT)
```

Unlike the runaway greenhouse test this uses `rscatter=True` and
`setTropopauseTemperature()` (dynamic skin-temperature tropopause).

### `test_constants.py`

Checks physical constants and utility functions:

- `gravity()` decreases with radius as expected
- `Earth` planet constants match database values (`g = 9.798`, `albedo = 0.306`, `Tsbar = 288.0`)
- `cpv()` returns finite, positive heat capacities for H₂O, CO₂, N₂, H₂, CH₄, O₂

### `test_code.py`

Checks utility code and the CLI:

- `setup_logger()` sets the correct log level and handler count; raises `ValueError` for invalid levels
- `natural_sort()` correctly sorts filenames with embedded numbers
- `janus env` CLI command exits successfully and prints `RAD_DIR`

