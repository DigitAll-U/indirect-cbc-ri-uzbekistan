"""Compare the derived intervals with published reference values.

Two comparators, for different purposes.

Primary - a direct-method study in a Turkish population:

    Ozarda Y, Ichihara K, Bakan E, Polat H, Ozturk N, Baygutalp NK, et al.
    A nationwide multicentre study in Turkey for establishing reference
    intervals of haematological parameters with novel use of a panel of whole
    blood. Biochem Med (Zagreb). 2017;27(2):350-377. doi:10.11613/BM.2017.038

    Values verified against the open full text on PubMed Central
    (PMC5493180, PMID 28694726). Direct estimation under the IFCC C-RIDL
    protocol, 3,363 healthy volunteers, 12 laboratories across 7 regions.
    It tests whether the indirect estimates agree with direct estimates
    from a comparable population.

Secondary - the textbook values that laboratories commonly use:

    Lewis SM, Bain BJ, Bates I, eds. Dacie and Lewis Practical Haematology.
    10th ed. Philadelphia: Churchill Livingstone/Elsevier; 2006.
    Chapter 2, "Reference ranges and normal values", Table 2.2, p. 14.

    This is not a validation standard. Laboratories without their own
    intervals often adopt textbook values, so the comparison shows the gap
    between those values and the population actually tested.

Reading Table 2.2. Its header reads "expressed as a mean +/-2SD (95% range)",
so the quoted spread is already 2SD: "5.0 +/- 0.5" means 4.5-5.5, not
4.0-6.0. This agrees with well-known ranges: haemoglobin 150 +/- 20 g/l gives
the male 130-170, and PCV 0.45 +/- 0.05 gives 0.40-0.50.

Paediatric comparison. Ozarda recruited only adults aged 18-79. Dacie & Lewis
Table 2.4, p. 15, gives values for normal children in three age bands (1 year,
2-6 years, 6-12 years), against our single 1-14 year partition, of which 3.7%
are aged 1, 40.3% 2-6, 43.2% 7-12 and 12.7% 13-14. No single band matches, so
the envelope of all three is used. The 1-year band is included because the
partition contains one-year-olds, and that band has the lowest red cell limits.

The 12.7% aged 13-14 lie beyond the oldest band, where values approach adult
ones, so the envelope is both wider than the interval for any single age and
short at the top of our age range. The paediatric comparison is therefore
weaker than the adult one, and is reported as containment (whether a derived
limit falls inside the envelope) rather than as a percentage difference.

Table 2.4 covers 12 of the 17 parameters; RDW, MPV, plateletcrit, platelet
distribution width and basophils are absent.

Neither source covers plateletcrit or platelet distribution width, in any
partition.
"""

import csv
import io
import os
import statistics
import sys

RESULT_DIR = os.environ.get("RESULT_DIR", "results")
IN_FILE = os.path.join(RESULT_DIR, "reference_intervals.csv")

ADULTS = ["Men >=15", "Women >=15"]

# Ozarda 2017. HCT is published in L/L and converted to % to match our units;
# that factor of 100 is the only unit conversion applied in this file.
OZARDA_BY_SEX = {
    "HGB": {"Men >=15": (131, 175), "Women >=15": (110, 152)},
    "RBC": {"Men >=15": (4.43, 6.07), "Women >=15": (3.96, 5.31)},
    "HCT": {"Men >=15": (39.2, 52.2), "Women >=15": (33.7, 46.1)},
}
OZARDA_COMBINED = {
    "MCV": (77.2, 95.7),
    "MCH": (25.2, 32.2),
    "MCHC": (319, 350),
    "RDW_CV": (12.2, 16.3),
    "PLT": (152, 383),
    "MPV": (6.2, 11.8),
    "WBC": (4.39, 11.59),
    "NEU_ABS": (2.04, 7.54),
    "LYMPH_ABS": (1.21, 3.77),
    "MON_ABS": (0.26, 0.94),
    "EOS_ABS": (0.02, 0.50),
    "BAS_ABS": (0.01, 0.12),
}

# Dacie & Lewis Table 2.2, converted from mean +/-2SD to the 95% range.
# The arithmetic is written out so each entry can be checked against the page.
DACIE_BY_SEX = {
    "HGB": {"Men >=15": (150 - 20, 150 + 20),          # 150 +/- 20 g/l
            "Women >=15": (135 - 15, 135 + 15)},       # 135 +/- 15 g/l
    "RBC": {"Men >=15": (5.0 - 0.5, 5.0 + 0.5),        # 5.0 +/- 0.5 x10^12/l
            "Women >=15": (4.3 - 0.5, 4.3 + 0.5)},     # 4.3 +/- 0.5
    "HCT": {"Men >=15": (100 * (0.45 - 0.05), 100 * (0.45 + 0.05)),
            "Women >=15": (100 * (0.41 - 0.05), 100 * (0.41 + 0.05))},
}
DACIE_COMBINED = {
    "MCV": (92 - 9, 92 + 9),                 # 92 +/- 9 fl, men and women
    "MCH": (29.5 - 2.5, 29.5 + 2.5),         # 29.5 +/- 2.5 pg
    "MCHC": (330 - 15, 330 + 15),            # 330 +/- 15 g/l
    "RDW_CV": (12.8 - 1.2, 12.8 + 1.2),      # 12.8 +/- 1.2 %
    "PLT": (280 - 130, 280 + 130),           # 280 +/- 130 x10^9/l
    "WBC": (4.0, 10.0),                      # printed as a range, not mean+/-SD
    "NEU_ABS": (2.0, 7.0),
    "LYMPH_ABS": (1.0, 3.0),
    "MON_ABS": (0.2, 1.0),
    "EOS_ABS": (0.02, 0.5),
    "BAS_ABS": (0.02, 0.1),
}

SOURCES = [
    ("Ozarda 2017", OZARDA_BY_SEX, OZARDA_COMBINED),
    ("Dacie & Lewis", DACIE_BY_SEX, DACIE_COMBINED),
]

# Dacie & Lewis Table 2.4, p. 15, "Haematological values for normal children".
# Same mean +/-2SD convention as Table 2.2 where a spread is given; the
# leucocyte differential, platelets and reticulocytes are printed as ranges.
DACIE_CHILD_1 = {
    "RBC": (4.5 - 0.6, 4.5 + 0.6),
    "HGB": (126 - 15, 126 + 15),
    "HCT": (100 * (0.34 - 0.04), 100 * (0.34 + 0.04)),
    "MCV": (78 - 6, 78 + 6),
    "MCH": (27 - 2, 27 + 2),
    "MCHC": (340 - 20, 340 + 20),
    "WBC": (11 - 5, 11 + 5),
    "NEU_ABS": (1, 7),
    "LYMPH_ABS": (3.5, 11),
    "MON_ABS": (0.2, 1.0),
    "EOS_ABS": (0.1, 1.0),
    "PLT": (200, 550),
}
DACIE_CHILD_2_6 = {
    "RBC": (4.6 - 0.6, 4.6 + 0.6),
    "HGB": (125 - 15, 125 + 15),
    "HCT": (100 * (0.37 - 0.03), 100 * (0.37 + 0.03)),
    "MCV": (81 - 6, 81 + 6),
    "MCH": (27 - 3, 27 + 3),
    "MCHC": (340 - 30, 340 + 30),
    "WBC": (10 - 5, 10 + 5),
    "NEU_ABS": (1.5, 8),
    "LYMPH_ABS": (6, 9),
    "MON_ABS": (0.2, 1.0),
    "EOS_ABS": (0.1, 1.0),
    "PLT": (200, 490),
}
DACIE_CHILD_6_12 = {
    "RBC": (4.6 - 0.6, 4.6 + 0.6),
    "HGB": (135 - 20, 135 + 20),
    "HCT": (100 * (0.40 - 0.05), 100 * (0.40 + 0.05)),
    "MCV": (86 - 9, 86 + 9),
    "MCH": (29 - 4, 29 + 4),
    "MCHC": (340 - 30, 340 + 30),
    "WBC": (9 - 4, 9 + 4),
    "NEU_ABS": (2, 8),
    "LYMPH_ABS": (1, 5),
    "MON_ABS": (0.2, 1.0),
    "EOS_ABS": (0.1, 1.0),
    "PLT": (170, 450),
}

# The lymphocyte bands are inconsistent: 3.5-11 at one year, 6-9 at 2-6 years
# and 1-5 at 6-12. Lymphocyte counts fall steadily through childhood, so the
# narrow 2-6 band is probably a misprint. The envelope is reported, but the
# parameter is marked and no conclusion should rest on it.
CHILD_SUSPECT = {"LYMPH_ABS"}

CHILD_PARTITION = "Children 1-14"

# Ozarda reports that these required partitioning by analyser manufacturer.
# We run a Mindray BC-6800Plus, which was not among the three platforms in the
# Turkish study (Abbott, Beckman Coulter, Sysmex), so differences on these
# parameters are expected from the platform alone.
MANUFACTURER_DEPENDENT = {"MCHC", "RDW_CV", "MPV", "BAS_ABS"}

NOT_REPORTED = ["PCT", "PDW"]

# The 17 parameters the paper reports, in the order used in the tables.
ALL_PARAMS = ["HGB", "RBC", "HCT", "MCV", "MCH", "MCHC", "RDW_CV",
              "PLT", "MPV", "PCT", "PDW",
              "WBC", "NEU_ABS", "LYMPH_ABS", "MON_ABS", "EOS_ABS", "BAS_ABS"]


def load():
    with io.open(IN_FILE, encoding="utf-8") as fh:
        return {(r["partition"], r["parameter"]):
                (float(r["ri_lower"]), float(r["ri_upper"]))
                for r in csv.DictReader(fh)}


def load_ci():
    """Bootstrap 90% confidence bounds on each limit, keyed as in load()."""
    out = {}
    with io.open(IN_FILE, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["lower_ci_lo"] in ("NA", ""):
                continue
            out[(r["partition"], r["parameter"])] = (
                (float(r["lower_ci_lo"]), float(r["lower_ci_hi"])),
                (float(r["upper_ci_lo"]), float(r["upper_ci_hi"])))
    return out


def ci_containment(ours, cis):
    """Does each published limit fall inside our 90% CI for that limit?

    This is the strict form of the comparison. A percentage difference says
    how far apart two limits are, but not whether they are distinguishable
    given the uncertainty of our estimate. Containment tests that.
    """
    print("\n" + "=" * 94)
    print("Published limits against our 90% confidence intervals")
    print("=" * 94)
    header = "%-11s %-11s %-13s %-22s %-9s" % (
        "parameter", "partition", "source", "our 90% CI", "contains?")
    print(header)
    print("-" * len(header))

    tally = {}
    for name, by_sex, combined in SOURCES:
        hit = tot = 0
        hits = []
        for p in list(by_sex) + list(combined):
            for prt in ADULTS:
                ref = by_sex[p][prt] if p in by_sex else combined[p]
                ci = cis.get((prt, p))
                if ci is None:
                    continue
                for i, lim in enumerate(("lower", "upper")):
                    tot += 1
                    lo, hi = ci[i]
                    if lo <= ref[i] <= hi:
                        hit += 1
                        hits.append("%-11s %-11s %-5s published %g in [%.3f, %.3f]"
                                    % (p, prt, lim, ref[i], lo, hi))
        tally[name] = (hit, tot)
        print("\n%s - the %d limits our CI does contain:" % (name, hit))
        for h in hits:
            print("  " + h)

    print()
    for name, (hit, tot) in tally.items():
        print("  %-15s %2d of %d published limits fall inside our 90%% CI (%.0f%%)"
              % (name, hit, tot, 100.0 * hit / tot))

    # The containment rate is low because, with this many observations, the
    # bootstrap CI is far narrower than any real difference between
    # populations or platforms. The CI widths below show this.
    widths = []
    for (prt, p), ((la, lb), (ua, ub)) in cis.items():
        est = ours.get((prt, p))
        if not est:
            continue
        for w, e in (((lb - la), est[0]), ((ub - ua), est[1])):
            if e:
                widths.append(100.0 * w / abs(e))
    widths.sort()
    n = len(widths)
    print("\n90%% CI width as %% of the limit it surrounds, across all %d limits:" % n)
    print("  median %.2f%%   p10 %.2f%%   p90 %.2f%%   max %.2f%%"
          % (statistics.median(widths), widths[n // 10], widths[(9 * n) // 10],
             widths[-1]))
    return tally


def compare(ours, name, by_sex, combined):
    print("\n" + "=" * 94)
    print("Against %s" % name)
    print("=" * 94)
    header = "%-11s %-11s %-19s %-19s %8s %8s %9s %9s" % (
        "parameter", "partition", "derived", name,
        "d_lower", "d_upper", "abs_low", "abs_up")
    print(header)
    print("-" * len(header))

    flagged = []
    compared = 0

    for p in list(by_sex) + list(combined):
        for prt in ADULTS:
            ref = by_sex[p][prt] if p in by_sex else combined[p]
            got = ours.get((prt, p))
            if got is None:
                print("%-11s %-11s  MISSING from results CSV" % (p, prt))
                continue
            compared += 1
            # Both relative and absolute differences are printed. The relative
            # difference supports statements such as "within 5%", but for a
            # limit near zero it overstates the gap: the eosinophil lower limit
            # can differ by 95% and by only 0.019 x10^9/L. The absolute columns
            # exist mainly for eosinophils and basophils.
            dl = 100.0 * (got[0] - ref[0]) / ref[0]
            du = 100.0 * (got[1] - ref[1]) / ref[1]
            al, au = got[0] - ref[0], got[1] - ref[1]
            # The 10% mark only highlights rows for reading; it is not a
            # clinical or statistical criterion.
            mark = " *" if max(abs(dl), abs(du)) >= 10 else ""
            print("%-11s %-11s %-19s %-19s %+7.1f%% %+7.1f%% %+9.3f %+9.3f%s" % (
                p, prt,
                "%g-%g" % (round(got[0], 3), round(got[1], 3)),
                "%g-%g" % (round(ref[0], 4), round(ref[1], 4)),
                dl, du, al, au, mark))
            if mark:
                flagged.append((p, prt, dl, du, al, au))

    print("\n%d intervals compared, %d differ by 10%% or more on a limit."
          % (compared, len(flagged)))
    for p, prt, dl, du, al, au in sorted(
            flagged, key=lambda f: -max(abs(f[2]), abs(f[3]))):
        cause = ("  <- manufacturer-dependent in Ozarda"
                 if p in MANUFACTURER_DEPENDENT else "")
        print("  %-11s %-11s lower %+6.1f%% (%+.3f)  upper %+6.1f%% (%+.3f)%s"
              % (p, prt, dl, al, du, au, cause))
    return compared, {(f[0], f[1]) for f in flagged}


def compare_children(ours):
    print("\n" + "=" * 94)
    print("Children 1-14 against Dacie & Lewis Table 2.4 "
          "(envelope of the 1, 2-6 and 6-12 year bands)")
    print("=" * 94)
    header = "%-11s %-19s %-19s %9s %9s" % (
        "parameter", "derived", "envelope", "lower_in", "upper_in")
    print(header)
    print("-" * len(header))

    bands = (DACIE_CHILD_1, DACIE_CHILD_2_6, DACIE_CHILD_6_12)
    inside = 0
    for p in DACIE_CHILD_2_6:
        lo = min(b[p][0] for b in bands)
        hi = max(b[p][1] for b in bands)
        got = ours.get((CHILD_PARTITION, p))
        if got is None:
            print("%-11s  MISSING from results CSV" % p)
            continue
        # Containment, not percentage difference. The envelope spans three age
        # bands, so whether a derived limit falls inside it is all the
        # comparison can show; how far inside carries no meaning.
        li = lo <= got[0] <= hi
        ui = lo <= got[1] <= hi
        inside += li and ui
        mark = "  <- see CHILD_SUSPECT" if p in CHILD_SUSPECT else ""
        print("%-11s %-19s %-19s %9s %9s%s" % (
            p, "%g-%g" % (round(got[0], 3), round(got[1], 3)),
            "%g-%g" % (round(lo, 4), round(hi, 4)),
            "yes" if li else "NO", "yes" if ui else "NO", mark))

    print("\n%d of %d paediatric intervals fall entirely within the envelope."
          % (inside, len(DACIE_CHILD_2_6)))
    print("Not covered by Table 2.4: %s"
          % ", ".join(sorted(set(ALL_PARAMS) - set(DACIE_CHILD_2_6))))


def main():
    ours = load()
    results = {}
    for name, by_sex, combined in SOURCES:
        results[name] = compare(ours, name, by_sex, combined)
    compare_children(ours)
    cis = load_ci()
    if cis:
        ci_containment(ours, cis)
    else:
        print("\n(no confidence bounds in the results file; "
              "run the estimation with NBootstrap > 0)")

    print("\n" + "=" * 94)
    print("Coverage")
    print("=" * 94)
    for name, (n, _) in results.items():
        print("  %-16s %2d of 51 derived intervals have a comparator"
              % (name, n))
    print("  %-16s %d of 17 paediatric intervals have a comparator "
          "(Dacie & Lewis Table 2.4 only)"
          % ("Children:", len(DACIE_CHILD_2_6)))
    print("  %-16s %s in both adult partitions"
          % ("no comparator:", " and ".join(NOT_REPORTED)))

    # A limit that disagrees with both sources more likely reflects our data;
    # one that disagrees with only one more likely reflects that source.
    #
    # This holds only where both sources give a value, so the comparison is
    # restricted to those parameters. MPV, for example, is in Ozarda but not
    # in Dacie & Lewis.
    both_cover = ((set(OZARDA_BY_SEX) | set(OZARDA_COMBINED)) &
                  (set(DACIE_BY_SEX) | set(DACIE_COMBINED)))
    o = {f for f in results["Ozarda 2017"][1] if f[0] in both_cover}
    d = {f for f in results["Dacie & Lewis"][1] if f[0] in both_cover}
    only_one = sorted((set(OZARDA_BY_SEX) | set(OZARDA_COMBINED)) ^
                      (set(DACIE_BY_SEX) | set(DACIE_COMBINED)))
    print("\nCompared on the %d parameters both sources publish; %s excluded "
          "from this breakdown (one source only)."
          % (len(both_cover), ", ".join(only_one)))
    for label, s in (("BOTH sources", o & d),
                     ("Ozarda only", o - d),
                     ("Dacie & Lewis only", d - o)):
        print("\nDisagree with %s (%d):" % (label, len(s)))
        for p, prt in sorted(s):
            print("  %-11s %s" % (p, prt))


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
