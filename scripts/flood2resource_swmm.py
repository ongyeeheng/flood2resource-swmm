"""Flood2Resource: EPA SWMM 5.2 model of one street-scale unit in MBZ City, Abu Dhabi.

Three configurations of the same 2,000 m2 street catchment:
  none       - no unit, all runoff goes to the municipal drain
  singapore  - Singapore ABC Waters-style rain garden (unlined, sandy loam, underdrain to the drain)
  f2r        - Flood2Resource: desert filter media, lined, underdrain and burst bypass to a 70 m3 tank,
               tank overflow to the drain, drip irrigation drawn from the tank between storms

Usage:
  python flood2resource_swmm.py events <out_dir>       design/observed storm runs
  python flood2resource_swmm.py longterm <out_dir>     1991-2025 continuous run on ERA5 hourly rain
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from pyswmm import Simulation, Nodes, Links

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"

# ---------------------------------------------------------------- design parameters (report Table 3)
CATCH_HA = 0.20          # 2,000 m2
IMPERV = 90              # % of non-LID area
GARDEN_M2 = 80
TANK_AREA = 35.0         # m2 plan area -> 70 m3 at 2.0 m
TANK_CREST = 2.0         # m, overflow weir crest
IRRIG_M3_DAY = 1.5       # drip to garden (80 m2) + ~300 m2 of verge at ~4 L/m2/day  [VERIFY local norm]
LAT = 24.34

LID = {
    # SURFACE: berm(mm) vegfrac rough slope xslope | SOIL: thick por fc wp ksat kcoef suction
    # STORAGE: thick voidratio seepage(mm/h) clog | DRAIN: coeff expon offset delay
    "singapore": dict(berm=250, soil=(600, 0.45, 0.20, 0.10, 150, 40, 80),
                      store=(200, 0.54, 5.0, 0), drain=(12.0, 0.5, 0, 6)),
    "f2r": dict(berm=300, soil=(600, 0.40, 0.10, 0.05, 300, 30, 50),
                store=(300, 0.54, 0.0, 0), drain=(17.0, 0.5, 0, 6)),
}

# NRCS (TR-55) Type II 24-hour cumulative rainfall distribution
TYPE_II = [(0, 0), (2, .022), (4, .048), (6, .080), (7, .098), (8, .120), (8.5, .133), (9, .147),
           (9.5, .163), (9.75, .172), (10, .181), (10.5, .204), (11, .235), (11.5, .283), (11.75, .357),
           (12, .663), (12.5, .735), (13, .772), (13.5, .799), (14, .820), (16, .880), (20, .952), (24, 1.0)]


def type_ii(depth_mm, step_min=5):
    t = np.arange(0, 24 * 60 + step_min, step_min) / 60
    cum = np.interp(t, *zip(*TYPE_II)) * depth_mm
    inc = np.diff(cum)
    idx = pd.to_timedelta(np.arange(len(inc)) * step_min, unit="min")  # integer minutes, no float drift
    return pd.Series(inc * 60 / step_min, index=idx)  # mm/h


def uniform(depth_mm, minutes, step_min=5):
    n = minutes // step_min
    return pd.Series([depth_mm / minutes * 60] * n, index=pd.to_timedelta(np.arange(n) * step_min, unit="min"))


def monthly_et0():
    """FAO-56 Hargreaves reference ET (mm/day) by calendar month from NASA POWER daily temperatures."""
    p = pd.read_csv(DATA / "power_mbz.csv", skiprows=_power_header_rows())
    p = p.replace(-999, np.nan)
    p["date"] = pd.to_datetime(p.YEAR.astype(str) + p.DOY.astype(str).str.zfill(3), format="%Y%j")
    j = p.date.dt.dayofyear.values
    phi = np.radians(LAT)
    dr = 1 + 0.033 * np.cos(2 * np.pi * j / 365)
    dec = 0.409 * np.sin(2 * np.pi * j / 365 - 1.39)
    ws = np.arccos(-np.tan(phi) * np.tan(dec))
    ra = 24 * 60 / np.pi * 0.0820 * dr * (ws * np.sin(phi) * np.sin(dec) + np.cos(phi) * np.cos(dec) * np.sin(ws))
    ra_mm = 0.408 * ra
    et0 = 0.0023 * (p.T2M + 17.8) * np.sqrt((p.T2M_MAX - p.T2M_MIN).clip(lower=0)) * ra_mm
    return et0.groupby(p.date.dt.month).mean().round(1).tolist()


def _power_header_rows():
    with open(DATA / "power_mbz.csv") as f:
        for i, line in enumerate(f):
            if line.startswith("-END HEADER-"):
                return i + 1


def build_inp(config, rain_file, start, end, init_tank=0.0, interval="0:05", evap=None, report_step="00:05:00"):
    evap = evap or [5.0] * 12
    has_lid = config in LID
    outlet = "TANK" if config == "f2r" else "DRAIN"
    s = [
        "[TITLE]", f"Flood2Resource {config}",
        "[OPTIONS]",
        "FLOW_UNITS LPS", "INFILTRATION HORTON", "FLOW_ROUTING KINWAVE",
        f"START_DATE {start:%m/%d/%Y}", f"START_TIME {start:%H:%M:%S}",
        f"REPORT_START_DATE {start:%m/%d/%Y}", f"REPORT_START_TIME {start:%H:%M:%S}",
        f"END_DATE {end:%m/%d/%Y}", f"END_TIME {end:%H:%M:%S}",
        f"REPORT_STEP {report_step}", "WET_STEP 00:01:00", "DRY_STEP 01:00:00", "ROUTING_STEP 0:00:30",
        "ALLOW_PONDING NO",
        "[EVAPORATION]", "MONTHLY " + " ".join(str(v) for v in evap), "DRY_ONLY NO",
        "[RAINGAGES]", f"RG1 INTENSITY {interval} 1.0 FILE \"{rain_file}\" RG1 MM",
        "[SUBCATCHMENTS]", f"S1 RG1 {outlet} {CATCH_HA} {IMPERV} 40 1.0 0",
        "[SUBAREAS]", "S1 0.012 0.15 1.5 5 0 OUTLET",
        "[INFILTRATION]", "S1 100 30 4 7 0",
    ]
    if has_lid:
        p = LID[config]
        s += ["[LID_CONTROLS]", "RG BC",
              f"RG SURFACE {p['berm']} 0.1 0.1 1.0 5",
              "RG SOIL " + " ".join(map(str, p["soil"])),
              "RG STORAGE " + " ".join(map(str, p["store"])),
              "RG DRAIN " + " ".join(map(str, p["drain"])),
              "[LID_USAGE]",
              f"S1 RG 1 {GARDEN_M2} 0 10 100 0 * {outlet} 0"]
    s += ["[OUTFALLS]", "DRAIN 0 FREE NO", "IRRIG 0 FREE NO"]
    if config == "f2r":
        q = IRRIG_M3_DAY * 1000 / 86400  # L/s
        s += ["[STORAGE]", f"TANK 0 {TANK_CREST + 0.3} {init_tank} FUNCTIONAL 0 0 {TANK_AREA} 0 0",
              "[WEIRS]", f"OVERFLOW TANK DRAIN TRANSVERSE {TANK_CREST} 1.84 NO 0 0",
              "[OUTLETS]", "DRIP TANK IRRIG 0.0 TABULAR/DEPTH irr NO",
              "[CURVES]", f"irr Rating 0 0", f"irr 0.02 {q:.5f}", f"irr 5 {q:.5f}",
              "[XSECTIONS]", "OVERFLOW RECT_OPEN 0.3 2.0 0 0"]
    else:
        s += ["[JUNCTIONS]", "J0 0 1 0 0 0"]  # unused placeholder keeps sections valid
    return "\n".join(s) + "\n"


def write_rain(series, start, path):
    """series: mm/h indexed by timedelta from start (or by timestamps)."""
    idx = series.index if isinstance(series.index, pd.DatetimeIndex) else start + series.index
    with open(path, "w") as f:
        for t, v in zip(idx, series.values):
            f.write(f"RG1 {t.year} {t.month} {t.day} {t.hour} {t.minute} {v:.3f}\n")


def run(inp_path, track=True, step_s=60):
    """Run SWMM; return time series of drain inflow (L/s), tank depth (m), irrigation flow (L/s)."""
    rows = []
    with Simulation(str(inp_path)) as sim:
        sim.step_advance(step_s)
        nodes = Nodes(sim)
        drain = nodes["DRAIN"]
        tank = nodes["TANK"] if "TANK" in [n.nodeid for n in nodes] else None
        irr = nodes["IRRIG"]
        for _ in sim:
            if track:
                rows.append((sim.current_time, drain.total_inflow, tank.depth if tank else 0.0, irr.total_inflow))
        stats = {"tank_end_m3": (tank.depth * TANK_AREA) if tank else 0.0,
                 "runoff_err_pct": sim.runoff_error, "routing_err_pct": sim.flow_routing_error}
    ts = pd.DataFrame(rows, columns=["time", "drain_lps", "tank_m", "irrig_lps"]).set_index("time")
    stats["drain_m3"] = ts.drain_lps.sum() * step_s / 1000
    stats["irrig_m3"] = ts.irrig_lps.sum() * step_s / 1000
    return ts, stats


# ---------------------------------------------------------------- experiments
EVENTS = {
    "mar2026_mbz": ("March 2026, MBZ City (78.7 mm, Type II 24 h)", lambda: type_ii(78.7)),
    "mar2026_wathba": ("March 2026, Al Wathba (88.2 mm, Type II 24 h)", lambda: type_ii(88.2)),
    "apr2024_coastal": ("April 2024, coastal proxy (143 mm, Type II 24 h)", lambda: type_ii(143.0)),
    "burst_30min": ("10-year 30-min burst, Sharjah IDF (29.6 mm in 30 min)", lambda: uniform(29.6, 30)),
    "typical_10": ("Typical storm 10 mm (Type II 24 h)", lambda: type_ii(10.0)),
    "typical_20": ("Typical storm 20 mm (Type II 24 h)", lambda: type_ii(20.0)),
    "typical_30": ("Typical storm 30 mm (Type II 24 h)", lambda: type_ii(30.0)),
}


def events(out):
    out.mkdir(parents=True, exist_ok=True)
    start = pd.Timestamp("2026-03-23 00:00")
    end = start + pd.Timedelta(hours=48)
    results, series = [], {}
    for key, (label, fn) in EVENTS.items():
        rain = fn()
        rain_mm = float((rain * 5 / 60).sum())
        rf = out / f"rain_{key}.dat"
        write_rain(rain, start, rf)
        for config, init in [("none", 0), ("singapore", 0), ("f2r", 0), ("f2r_halffull", TANK_CREST / 2)]:
            cfg = "f2r" if config.startswith("f2r") else config
            inp = out / f"{key}_{config}.inp"
            inp.write_text(build_inp(cfg, rf.name, start, end, init_tank=init))
            ts, st = run(inp)
            series[(key, config)] = ts
            peak = ts.drain_lps.max()
            results.append(dict(event=key, label=label, rain_mm=round(rain_mm, 1), config=config,
                                peak_drain_lps=round(peak, 1),
                                t_peak_h=round((ts.drain_lps.idxmax() - start).total_seconds() / 3600, 2) if peak > 0.01 else None,
                                drain_m3=round(st["drain_m3"], 1),
                                kept_m3=round(st["tank_end_m3"] - init * TANK_AREA + st["irrig_m3"], 1),
                                runoff_err=round(st["runoff_err_pct"], 2), routing_err=round(st["routing_err_pct"], 2)))
    df = pd.DataFrame(results)
    base = df[df.config == "none"].set_index("event")
    df["peak_reduction_pct"] = df.apply(lambda r: round(100 * (1 - r.peak_drain_lps / base.loc[r.event, "peak_drain_lps"]), 1), axis=1)
    df["volume_reduction_pct"] = df.apply(lambda r: round(100 * (1 - r.drain_m3 / base.loc[r.event, "drain_m3"]), 1), axis=1)
    df.to_csv(out / "event_results.csv", index=False)
    pd.concat({f"{k}|{c}": v.drain_lps for (k, c), v in series.items()}, axis=1).to_csv(out / "event_hydrographs.csv")
    pd.concat({f"{k}|{c}": v.tank_m for (k, c), v in series.items() if c.startswith("f2r")}, axis=1).to_csv(out / "event_tank_depth.csv")
    print(df.to_string(index=False))


def longterm(out):
    """Continuous 1991-2025 run, hourly ERA5 rainfall at MBZ City."""
    out.mkdir(parents=True, exist_ok=True)
    e = pd.read_csv(DATA / "era5_hourly_mbz.csv", skiprows=3, parse_dates=["time"]).set_index("time")["precipitation (mm)"]
    e = e["1991-01-01":"2025-12-31 23:00"]
    keep = (e > 0) | (e.shift(1) > 0)          # rain hours plus the zero that ends each burst
    rf = out / "rain_era5_1991_2025.dat"
    write_rain(e[keep], None, rf)
    evap = monthly_et0()
    start, end = pd.Timestamp("1991-01-01"), pd.Timestamp("2026-01-01")
    res = {"et0_monthly_mm_day": evap}
    for config in ["none", "singapore", "f2r"]:
        inp = out / f"longterm_{config}.inp"
        inp.write_text(build_inp(config, rf.name, start, end, interval="1:00", evap=evap, report_step="01:00:00"))
        ts, st = run(inp, track=True, step_s=3600)
        ts = ts.resample("1h").mean().fillna(0)
        annual = pd.DataFrame({"drain_m3": ts.drain_lps.resample("YE").sum() * 3.6,
                               "irrig_m3": ts.irrig_lps.resample("YE").sum() * 3.6,
                               "rain_mm": e.resample("YE").sum().values})
        res[config] = {"drain_m3_total": round(st["drain_m3"], 1), "irrig_m3_total": round(st["irrig_m3"], 1),
                       "runoff_err": st["runoff_err_pct"], "routing_err": st["routing_err_pct"],
                       "irrig_m3_per_year_mean": round(annual.irrig_m3.mean(), 1)}
        if config == "f2r":
            days_full = (ts.tank_m.resample("D").max() >= TANK_CREST - 0.01).sum()
            res[config]["days_tank_at_overflow"] = int(days_full)
            ts.tank_m.resample("D").mean().to_csv(out / "longterm_tank_depth_daily.csv")
        annual.to_csv(out / f"longterm_annual_{config}.csv")
    res["rain_mm_per_year_mean"] = round(e.resample("YE").sum().mean(), 1)
    (out / "longterm_results.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))


def sensitivity(out):
    """One-at-a-time sensitivity for the March 2026 MBZ City storm."""
    global TANK_AREA, IMPERV
    out.mkdir(parents=True, exist_ok=True)
    start = pd.Timestamp("2026-03-23 00:00")
    end = start + pd.Timedelta(hours=48)
    rf = out / "rain_mar2026_mbz.dat"
    write_rain(type_ii(78.7), start, rf)
    base_inp = out / "base_none.inp"
    base_inp.write_text(build_inp("none", rf.name, start, end))
    base_ts, base_st = run(base_inp)
    rows = []
    defaults = (TANK_AREA, IMPERV, LID["f2r"]["soil"])
    cases = [("tank_m3", v) for v in (50, 70, 100, 150)] + [("imperv_pct", v) for v in (80, 90, 100)] + \
            [("ksat_mm_h", v) for v in (100, 200, 300, 400)]
    for param, val in cases:
        TANK_AREA, IMPERV, soil = defaults
        if param == "tank_m3":
            TANK_AREA = val / TANK_CREST
        elif param == "imperv_pct":
            IMPERV = val
        else:
            soil = soil[:4] + (val,) + soil[5:]
        LID["f2r"]["soil"] = soil
        if param == "imperv_pct":  # baseline changes with imperviousness
            b = out / f"none_{val}.inp"; b.write_text(build_inp("none", rf.name, start, end))
            bts, bst = run(b)
        else:
            bts, bst = base_ts, base_st
        inp = out / f"f2r_{param}_{val}.inp"
        inp.write_text(build_inp("f2r", rf.name, start, end))
        ts, st = run(inp)
        rows.append(dict(parameter=param, value=val,
                         peak_reduction_pct=round(100 * (1 - ts.drain_lps.max() / bts.drain_lps.max()), 1),
                         volume_kept_pct=round(100 * (1 - st["drain_m3"] / bst["drain_m3"]), 1)))
    TANK_AREA, IMPERV, LID["f2r"]["soil"] = defaults
    df = pd.DataFrame(rows)
    df.to_csv(out / "sensitivity_mar2026_mbz.csv", index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    {"events": events, "longterm": longterm, "sensitivity": sensitivity}[sys.argv[1]](Path(sys.argv[2]))
