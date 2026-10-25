# Flood2Resource: SWMM model of a street-scale flood-to-irrigation unit, MBZ City, Abu Dhabi

Simulation behind the Flood2Resource submission (team CMAT-00067, NUS) to the Khalifa University and UNESCO Global Water Hackathon 2026, Theme 5.

The model uses EPA SWMM 5.2, run from Python with PySWMM. It represents one MBZ City residential street (2,000 m², 90 % paved) under three designs:

- **none**: no unit, all runoff goes to the municipal drain.
- **singapore**: Singapore ABC Waters-style rain garden (unlined, underdrain to the drain).
- **f2r**: Flood2Resource unit. An 80 m² bio-retention LID (300 mm berm, 600 mm filter, 300 mm gravel, lined) whose underdrain and overflow fill a 70 m³ storage tank. The tank supplies 1.5 m³/day of drip irrigation and overflows to the drain when full.
- **f2r_halffull**: the Flood2Resource unit with the tank half full at the start of the storm.

Storms: March 2026 MBZ City (78.7 mm) and Al Wathba (88.2 mm), an April 2024-scale storm (143 mm), a 10-year 30-minute burst (29.6 mm), ordinary 10, 20 and 30 mm storms, and a continuous 1991 to 2025 run on ERA5 hourly rainfall. Evaporation comes from NASA POWER temperatures via the Hargreaves equation.

**Results are simulated, not field measurements.**

## Key results (March 2026, MBZ City, 78.7 mm)

| Design | Peak flow to drain (L/s) | Peak cut |
|---|---|---|
| No unit | 51.1 | |
| Singapore-style rain garden | 50.8 | 0.6 % |
| Flood2Resource | 9.7 | 81 % |
| Flood2Resource, tank half full | 39.9 | 22 % |

Full tables: `results/events/event_results.csv`, `results/sensitivity/sensitivity_mar2026_mbz.csv`, `results/longterm/longterm_results.json`.

## Requirements

- Python 3.10 or later
- `pip install pyswmm pandas numpy matplotlib` (PySWMM bundles the EPA SWMM 5.2 engine)
- Optional: the [EPA SWMM 5.2 desktop program](https://www.epa.gov/water-research/storm-water-management-model-swmm) to open the `.inp` files directly

## How to run

From the repository root:

```
python scripts/flood2resource_swmm.py events results/events
python scripts/flood2resource_swmm.py sensitivity results/sensitivity
python scripts/flood2resource_swmm.py longterm results/longterm
python scripts/make_figures.py
```

`events` and `sensitivity` take a few minutes; `longterm` (35 years of hourly rain) takes longer. Each command rewrites the `.inp`, rainfall `.dat` and results files in its folder. The binary `.out` files are not included (the long-term ones are about 20 MB each); they are recreated on a run.

`make_figures.py` writes the report's results figures to `figures/` (hydrograph, peak reduction by storm, and the combined two-panel figure). It uses the Calibri font and falls back to the default font if Calibri is missing.

## Files

| Path | Contents |
|---|---|
| `scripts/flood2resource_swmm.py` | Builds every SWMM input file, runs it with PySWMM and writes the results tables. All design parameters are at the top. |
| `scripts/make_figures.py` | Report figures from the results CSVs. |
| `data/era5_hourly_mbz.csv` | ERA5 hourly precipitation at MBZ City, 1991 to 2026. |
| `data/power_mbz.csv` | NASA POWER daily temperature and precipitation at MBZ City, 1991 to March 2026. |
| `results/events/` | SWMM inputs (`<storm>_<design>.inp`), rainfall series (`rain_*.dat`), SWMM reports (`.rpt`), `event_results.csv` (peak and volume per storm and design), `event_hydrographs.csv` (flow to drain, 1-minute), `event_tank_depth.csv`. |
| `results/sensitivity/` | One-at-a-time sensitivity for the March 2026 MBZ City storm: tank 50/70/100/150 m³, paved share 80/90/100 %, filter infiltration rate 100/200/300/400 mm/h. |
| `results/longterm/` | 1991 to 2025 continuous run: inputs, annual drain and irrigation volumes per design, daily tank depth, summary JSON. |
| `swmm_gui/flood2resource_mbz_detailed.inp` | Display version of the Flood2Resource March 2026 model for the SWMM desktop program, drawn on the real street with a drain line and manholes (dynamic-wave routing). Open it in SWMM 5.2 together with `rain_mar2026_mbz.dat` in the same folder. The satellite backdrop is not included. |

## Data credits

- **ERA5**: Hersbach, H. et al. (2020), The ERA5 global reanalysis, *Quarterly Journal of the Royal Meteorological Society* 146, 1999 to 2049. Copernicus Climate Change Service (C3S). Hourly series retrieved through the Open-Meteo historical weather API.
- **NASA POWER**: data obtained from the NASA Langley Research Center (LaRC) POWER Project, funded through the NASA Earth Science Applied Science Program.
- **Storm totals** for March 2026 and April 2024: UAE National Centre of Meteorology (NCM), as reported in the press.
- Storm shape: NRCS (TR-55) Type II 24-hour distribution. 30-minute burst depth from a Sharjah IDF curve.
