"""Apply validity rules to the full CBC panel and write a validated dataset
and a report.

The rationale for each rule is in RULES below, which is also written to the
report.

This is a separate step from prepare_cohort.py because the consistency
checks need several parameters of the same sample at once, which exist only
after the data are pivoted to one row per sample.
"""
import csv
import os
import statistics

DATA_DIR   = os.environ.get("DATA_DIR", "data")
RESULT_DIR = os.environ.get("RESULT_DIR", "results")

IN_FILE  = f"{DATA_DIR}/cbc_cohort.csv"
OUT_FILE = f"{DATA_DIR}/cbc_validated.csv"
REPORT   = f"{RESULT_DIR}/validation_report.txt"

# --- tolerances for calculated parameters -------------------------------
# From the manufacturer's acceptance criteria (the Mindray BC-6800 operator
# documentation gives MCH +/-5%, MCHC +/-7%, RDW-CV +/-10%, RDW-SD +/-15%,
# MPV +/-15%, PLT +/-10%), rather than thresholds of our own. Observed
# agreement is much closer (median discrepancy 0.08-0.27%, 95th percentile
# <0.75%), so these tolerances catch only corrupt records.
TOL = {"MCV": 5.0, "MCH": 5.0, "MCHC": 7.0, "PCT": 15.0}

# --- physiological limits for measured parameters -----------------------
# Measured values cannot be checked by recomputation, so they get explicit
# limits. The limits are wide: they exclude values incompatible with life,
# not unusual pathology (e.g. WBC above 300 x10^9/L occurs in leukaemia and
# is kept).
MEASURED_BOUNDS = {
    # HGB 14 g/L: at the lower edge of values reported in any surviving
    # patient (see prepare_cohort.py).
    "HGB": (14, 250),
    # RBC below ~1.0 x10^12/L is incompatible with life; the upper limit is
    # the analyser's stated linearity for this class of instrument.
    "RBC": (1.0, 8.5),
    # HCT below 10% corresponds to a haemoglobin far below the HGB limit
    # above; 70% exceeds the highest recorded polycythaemia.
    "HCT": (10, 70),
    # PLT: values below 1 x10^9/L are another parameter in the wrong field;
    # the upper limit is the highest reported reactive thrombocytosis.
    "PLT": (1, 1703),
    # WBC: the lower limit allows profound aplasia, the upper limit
    # hyperleukocytosis.
    "WBC": (0.1, 500),
}

# --- limits for calculated parameters ------------------------------------
# Passing the consistency check is not enough: two wrong inputs can give a
# consistent but impossible index (seen here as MCHC 862 g/L and MCH 108 pg,
# both of which passed). These limits are applied as well.
CALCULATED_BOUNDS = {
    # MCHC is limited by how much haemoglobin a red cell can hold. Values
    # above about 400 g/L indicate interference (cold agglutinins, lipaemia,
    # in-vitro haemolysis or spherocytosis), not a true result.
    "MCHC": (200, 400),
    # MCH (haemoglobin per cell) does not approach 50 pg even in severe
    # macrocytic anaemia.
    "MCH": (10, 50),
    # Red cell distribution width: wide limits that keep severe
    # anisocytosis and exclude impossible values.
    "RDW_CV": (5, 60),
    "RDW_SD": (20, 120),
    # Platelet indices.
    "MPV": (5, 20),
    "PDW": (5, 30),
    "PCT": (0.001, 2.0),
    # Nucleated red cells: a count this high would exceed the total
    # leukocyte count.
    "NRBC_ABS": (0, 50),
}

RULES = """
VALIDATION RULES APPLIED
========================
1. Internal consistency of calculated parameters. MCV, MCH, MCHC and PCT
   are computed by the analyser from measured inputs and were recomputed
   independently:
       MCV  = HCT * 10 / RBC
       MCH  = HGB / RBC
       MCHC = HGB * 100 / HCT
       PCT  = PLT * MPV / 10000
   An encounter was rejected if any stored value disagreed with its
   recomputed value by more than the manufacturer's stated acceptance
   tolerance for that parameter. Because the discrepancy cannot be
   attributed to one specific parameter (either the stored index or one
   of its inputs may be at fault), the whole encounter is rejected rather
   than nulling a single field.

2. Physiological limits on directly measured parameters, which cannot be
   verified by recomputation. Limits are wide by design, chosen to exclude
   values incompatible with life rather than unusual pathology.

3. Mathematical constraints on the leukocyte differential: each
   percentage must lie in 0-100, and NEU+LYMPH+MON+EOS+BAS must sum to
   100 +/- 5.
"""


def num(row, key):
    v = row.get(key, "")
    if v in ("", None):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def main():
    with open(IN_FILE, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)

    reasons = {}
    kept = []

    for r in rows:
        why = []

        # --- rule 1: consistency of calculated parameters ---
        checks = [
            ("MCV", num(r, "MCV"),
             (num(r, "HCT") * 10 / num(r, "RBC"))
             if num(r, "HCT") and num(r, "RBC") else None),
            ("MCH", num(r, "MCH"),
             (num(r, "HGB") / num(r, "RBC"))
             if num(r, "HGB") and num(r, "RBC") else None),
            ("MCHC", num(r, "MCHC"),
             (num(r, "HGB") * 100 / num(r, "HCT"))
             if num(r, "HGB") and num(r, "HCT") else None),
            ("PCT", num(r, "PCT"),
             (num(r, "PLT") * num(r, "MPV") / 10000)
             if num(r, "PLT") and num(r, "MPV") else None),
        ]
        for name, stored, calc in checks:
            if stored is None or calc is None or calc == 0:
                continue
            if 100 * abs(stored - calc) / calc > TOL[name]:
                why.append(f"inconsistent_{name}")

        # --- rule 2: physiological limits on measured parameters ---
        for p, (lo, hi) in MEASURED_BOUNDS.items():
            v = num(r, p)
            if v is not None and not (lo <= v <= hi):
                why.append(f"out_of_range_{p}")

        # --- rule 2b: absolute limits on calculated parameters ---
        for p, (lo, hi) in CALCULATED_BOUNDS.items():
            v = num(r, p)
            if v is not None and not (lo <= v <= hi):
                why.append(f"out_of_range_{p}")

        # --- rule 2c: the red cell channel must have reported ---
        # The rules above apply only to values that are present, so a record
        # missing a whole analyser channel passes them all. Six of 30,774
        # records have a differential but no haemoglobin, red cell count or
        # platelet count (or the reverse). These are partial results or failed
        # transmissions; a real WBC-only test order would occur far more often.
        if num(r, "HGB") is None:
            why.append("incomplete_panel_no_HGB")

        # --- rule 3: differential constraints ---
        pcts = [num(r, f) for f in
                ("NEU_PCT", "LYMPH_PCT", "MON_PCT", "EOS_PCT", "BAS_PCT")]
        if all(v is not None for v in pcts):
            if any(not (0 <= v <= 100) for v in pcts):
                why.append("differential_out_of_0_100")
            elif not (95 <= sum(pcts) <= 105):
                why.append("differential_sum")

        if why:
            for w in why:
                reasons[w] = reasons.get(w, 0) + 1
            continue
        kept.append(r)

    with open(OUT_FILE, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(kept)

    lines = [RULES, ""]
    lines.append(f"input encounters   : {len(rows)}")
    lines.append(f"rejected           : {len(rows) - len(kept)}"
                 f"  ({100*(len(rows)-len(kept))/len(rows):.2f}%)")
    lines.append(f"validated output   : {len(kept)}")
    lines.append("")
    lines.append("rejections by reason (an encounter may fail several):")
    for k, v in sorted(reasons.items(), key=lambda kv: -kv[1]):
        lines.append(f"   {k:28s} {v}")

    lines.append("")
    lines.append("PARAMETER DISTRIBUTIONS AFTER VALIDATION")
    lines.append(f"{'param':10s} {'n':>7s} {'min':>9s} {'p1':>8s} "
                 f"{'median':>8s} {'p99':>9s} {'max':>10s}")
    params = ["HGB", "RBC", "HCT", "MCV", "MCH", "MCHC", "RDW_CV", "RDW_SD",
              "PLT", "MPV", "PCT", "PDW", "WBC", "NEU_ABS", "LYMPH_ABS",
              "MON_ABS", "EOS_ABS", "BAS_ABS", "IMG_ABS", "NRBC_ABS"]
    for p in params:
        vals = sorted(v for v in (num(r, p) for r in kept) if v is not None)
        if not vals:
            continue
        n = len(vals)
        lines.append(f"{p:10s} {n:7d} {vals[0]:9.2f} {vals[int(0.01*n)]:8.2f} "
                     f"{statistics.median(vals):8.2f} {vals[int(0.99*n)]:9.2f} "
                     f"{vals[-1]:10.2f}")

    with open(REPORT, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    print("\n".join(lines))
    print(f"\nSaved: {OUT_FILE}")
    print(f"Saved: {REPORT}")


if __name__ == "__main__":
    main()
