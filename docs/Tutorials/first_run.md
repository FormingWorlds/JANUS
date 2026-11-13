# First run tutorial

This tutorial runs the runaway greenhouse example that ships with JANUS. By the
end you will have a plot of outgoing longwave radiation (OLR) versus surface
temperature for a pure steam atmosphere, validated against four published
datasets.

!!! info "Prerequisites"
    Complete the [installation](../How-to/installation.md) steps first. You will also need
    the `mors` stellar evolution package and `toml`:

    ```bash
    pip install -r examples/requirements.txt
    ```

## 1. What the example does

`examples/demo_runaway_greenhouse.py` increases the surface temperature of a
pure H₂O atmosphere from 200 K to 2800 K in 20 steps, computing the OLR at
each step using the full JANUS pipeline (moist pseudoadiabat + SOCRATES
radiative transfer). The result is the classic **runaway greenhouse curve**,
plotted against reference data from Kopparapu et al. (2013)[^cite-kopparapu2013], Goldblatt et al. (2013)[^cite-goldblatt2013],
Hamano et al. (2015)[^cite-hamano2015], and Selsis et al. (2023)[^cite-selsis2023].

## 2. Run it

From the root of the JANUS repository:

```bash
python examples/demo_runaway_greenhouse.py
```

The run takes a few minutes, since each of the 20 temperature steps calls SOCRATES
twice. You should see progress logged to the terminal:

```
INFO - Start
INFO - Running JANUS...
INFO - T_surf = 200 K
...
INFO - T_surf = 2800 K
INFO - Making plot
INFO - Done!
```

Output files are written to `output/`:

```
output/runaway_demo.pdf
output/runaway_demo.png
```

---

## 3. Understanding the output

The plot shows OLR as a function of surface temperature.

**The plateau**: OLR barely changes between ~350 K and ~1600 K, staying near
277 W m⁻². This is the **Simpson–Nakajima radiation limit**: the photosphere
is trapped within the saturated water vapour region, so increasing the surface
temperature does not increase the emitting temperature. Any planet receiving
more stellar radiation than this limit cannot radiate fast enough and enters a
runaway greenhouse; its oceans evaporate entirely.

**The steep rise**: above ~1800 K the deep atmosphere becomes supercritical
(water has no condensed phase above 647 K). The condensing region disappears,
the photosphere couples back to the surface, and OLR rises steeply into the
post-runaway regime.

**The JANUS curve** should sit close to Hamano et al. (2015)[^cite-hamano2015] and Selsis et al. (2023)[^cite-selsis2023]. Small
differences between models come from different spectral line databases and
lapse rate assumptions.

## 4. The other example

`examples/demo_instellation.py` increases the **stellar instellation** instead of
the surface temperature, showing how the energy balance changes with orbital
distance. Run it the same way:

```console
python examples/demo_instellation.py
```

## 5. Next steps

- Modify `vol_mixing` in the script to add a CO₂ or N₂ background gas and see
  how it shifts the radiation limit
- Change `mean_distance` in the config located at [src/janus/data/tests/config_runaway](https://github.com/FormingWorlds/JANUS/blob/main/src/janus/data/tests/config_runaway.toml) to move the planet closer or further
  from its star
- See the [physical model overview](../Explanations/model.md) for the equations
  behind the pseudoadiabat and radiative transfer

[^cite-kopparapu2013]: R. K. Kopparapu, R. Ramirez, J. F. Kasting, V. Eymet, T. D. Robinson, S. Mahadevan, R. C. Terrien, S. Domagal-Goldman, V. Meadows, R. Deshpande, *[Habitable Zones around Main-sequence Stars: New Estimates](https://doi.org/10.1088/0004-637X/765/2/131)*, The Astrophysical Journal, 765, 131, 2013. [SciX](https://scixplorer.org/abs/2013ApJ...765..131K/abstract).

[^cite-goldblatt2013]: C. Goldblatt, T. D. Robinson, K. J. Zahnle, D. Crisp, *[Low simulated radiation limit for runaway greenhouse climates](https://doi.org/10.1038/ngeo1892)*, Nature Geoscience, 6, 661-667, 2013. [SciX](https://scixplorer.org/abs/2013NatGe...6..661G/abstract).

[^cite-hamano2015]: K. Hamano, H. Kawahara, Y. Abe, M. Onishi, G. L. Hashimoto, *[Lifetime and Spectral Evolution of a Magma Ocean with a Steam Atmosphere: Its Detectability by Future Direct Imaging](https://doi.org/10.1088/0004-637X/806/2/216)*, The Astrophysical Journal, 806, 216, 2015. [SciX](https://scixplorer.org/abs/2015ApJ...806..216H/abstract).

[^cite-selsis2023]: F. Selsis, J. Leconte, M. Turbet, G. Chaverot, É. Bolmont, *[A cool runaway greenhouse without surface magma ocean](https://doi.org/10.1038/s41586-023-06258-3)*, Nature, 620, 287-291, 2023. [SciX](https://scixplorer.org/abs/2023Natur.620..287S/abstract).
