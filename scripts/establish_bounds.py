"""Diagnostic for the validity rules applied in validate_panel.py.

Three checks, the first being the most reliable:

 1. Internal consistency (independent of any manufacturer's manual).
    MCV, MCH, MCHC and PCT are calculated by the analyser from measured
    values, so Mindray gives no linearity range for them. Each can be
    recomputed from its inputs and compared with the stored value:
        MCV  (fL)   = HCT(%) * 10 / RBC(10^12/L)
        MCH  (pg)   = HGB(g/L) / RBC(10^12/L)
        MCHC (g/L)  = HGB(g/L) * 100 / HCT(%)
        PCT  (%)    = PLT(10^9/L) * MPV(fL) / 10000
    A record whose stored value disagrees with its own inputs is
    unreliable, even if the value looks plausible on its own.

 2. Constraints on the differential. Percentages must lie in 0-100 and
    the five main fractions must sum to about 100.

 3. Physiological limits from the literature, for measured parameters,
    where consistency cannot be checked.

This script applies nothing. It reports how many records fail each check
and by how much, as the basis for the limits in validate_panel.py.
"""
import csv
import os
import statistics

DATA_DIR   = os.environ.get("DATA_DIR", "data")
with open(f"{DATA_DIR}/cbc_cohort.csv", encoding="utf-8") as f:
    rows = list(csv.DictReader(f))

print(f"records: {len(rows)}\n")


def val(row, key):
    v = row.get(key, "")
    if v == "" or v is None:
        return None
    try:
        return float(v)
    except ValueError:
        return None


def report_consistency(name, stored_key, fn, tol_pct):
    """Compare stored vs recomputed value; report disagreement."""
    diffs = []
    missing = 0
    for r in rows:
        stored = val(r, stored_key)
        calc = fn(r)
        if stored is None or calc is None or calc == 0:
            missing += 1
            continue
        pct_diff = 100 * abs(stored - calc) / calc
        diffs.append((pct_diff, stored, calc))

    if not diffs:
        print(f"{name}: no comparable records")
        return

    diffs.sort(reverse=True)
    over = [d for d in diffs if d[0] > tol_pct]
    pcts = [d[0] for d in diffs]
    print(f"{name} (stored vs recomputed, tolerance {tol_pct}%)")
    print(f"   comparable records : {len(diffs)}  (skipped {missing})")
    print(f"   median discrepancy : {statistics.median(pcts):.2f}%")
    # Interpolated 95th percentile, consistent with the median above.
    p95 =statistics.quantiles(pcts, n=20)[-1] if len(pcts) > 1 else pcts[0]
    print(f"   95th percentile    : {p95:.2f}%")
    print(f"   EXCEEDING tolerance: {len(over)}  ({100*len(over)/len(diffs):.2f}%)")
    for pct, stored, calc in diffs[:5]:
        print(f"      stored={stored:9.2f}  recomputed={calc:8.2f}  diff={pct:8.1f}%")
    print()


print("=" * 68)
print("1. INTERNAL CONSISTENCY OF CALCULATED PARAMETERS")
print("=" * 68)

report_consistency(
    "MCV", "MCV",
    lambda r: (val(r, "HCT") * 10 / val(r, "RBC"))
    if val(r, "HCT") and val(r, "RBC") else None,
    tol_pct=5,
)
report_consistency(
    "MCH", "MCH",
    lambda r: (val(r, "HGB") / val(r, "RBC"))
    if val(r, "HGB") and val(r, "RBC") else None,
    tol_pct=5,
)
report_consistency(
    "MCHC", "MCHC",
    lambda r: (val(r, "HGB") * 100 / val(r, "HCT"))
    if val(r, "HGB") and val(r, "HCT") else None,
    tol_pct=7,
)
report_consistency(
    "PCT", "PCT",
    lambda r: (val(r, "PLT") * val(r, "MPV") / 10000)
    if val(r, "PLT") and val(r, "MPV") else None,
    tol_pct=15,
)

print("=" * 68)
print("2. MATHEMATICAL CONSTRAINTS ON THE DIFFERENTIAL")
print("=" * 68)

pct_fields = ["NEU_PCT", "LYMPH_PCT", "MON_PCT", "EOS_PCT", "BAS_PCT"]
out_of_range = {f: 0 for f in pct_fields}
sums = []
bad_sum = 0
for r in rows:
    vals = [val(r, f) for f in pct_fields]
    for f, v in zip(pct_fields, vals):
        if v is not None and not (0 <= v <= 100):
            out_of_range[f] += 1
    if all(v is not None for v in vals):
        s = sum(vals)
        sums.append(s)
        if not (95 <= s <= 105):
            bad_sum += 1

print("values outside 0-100:")
for f, n in out_of_range.items():
    print(f"   {f:10s} {n}")
if sums:
    print(f"\ndifferential sum (NEU+LYMPH+MON+EOS+BAS):")
    print(f"   records checked : {len(sums)}")
    print(f"   median sum      : {statistics.median(sums):.2f}%")
    print(f"   outside 95-105% : {bad_sum}  ({100*bad_sum/len(sums):.2f}%)")

print()
print("=" * 68)
print("3. DISTRIBUTION OF EACH PARAMETER (for setting literature bounds)")
print("=" * 68)
params = ["HGB", "RBC", "HCT", "MCV", "MCH", "MCHC", "RDW_CV", "RDW_SD",
          "PLT", "MPV", "PCT", "PDW", "WBC", "NEU_ABS", "LYMPH_ABS",
          "MON_ABS", "EOS_ABS", "BAS_ABS", "IMG_ABS", "NRBC_ABS"]
print(f"{'param':10s} {'n':>7s} {'min':>9s} {'p1':>8s} {'median':>8s} {'p99':>9s} {'max':>10s}")
for p in params:
    vals = sorted(v for v in (val(r, p) for r in rows) if v is not None)
    if not vals:
        continue
    n = len(vals)
    print(f"{p:10s} {n:7d} {vals[0]:9.2f} {vals[int(0.01*n)]:8.2f} "
          f"{statistics.median(vals):8.2f} {vals[int(0.99*n)]:9.2f} {vals[-1]:10.2f}")
