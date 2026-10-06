"""
Build one de-identified annual dataset from the monthly raw lab exports
(see README for the raw data directory).

Standard library only (zipfile + re); openpyxl and pandas are not needed.

Extraction:
 - finds each file's columns by header name, not position, because the
   files use at least three layouts (a Turkish LIS export and two Russian
   variants, with VOZRAST or YEAR as the age header)
 - does not load FIO/ADRES/TEL/ADRES1/TELEFON into the working data
 - drops rows with corrupted sex values (header text such as "POL" or
   "CINSIYET" appearing as data)
 - excludes non-numeric results (haemolysis, clotted sample, insufficient
   volume, instrument error, bare "*"/"****") and implausible zeros, and
   counts each reason
 - pivots to one row per patient-test-event, one column per CBC parameter
 - assigns a new sequential ID unrelated to the original REG_NO
"""
import zipfile
import re
import csv
import os
from collections import Counter, defaultdict

RAW_DIR    = os.environ.get("RAW_DIR", "raw_data")  # source LIS exports (contain identifiers)
DATA_DIR   = os.environ.get("DATA_DIR", "data")                   # de-identified output
RESULT_DIR = os.environ.get("RESULT_DIR", "results")                # exclusion table and logs

FILES = [
    ("ЯНВАРЬ.xlsx", "2023-01"),
    ("ФЕВРАЛЬ.xlsx", "2023-02"),
    ("МАРТ.xlsx", "2023-03"),
    ("АПРЕЛЬ.xlsx", "2023-04"),
    ("МАЙ.xlsx", "2023-05"),
    ("ИЮНЬ.xlsx", "2023-06"),
    ("ИЮЛЬ.xlsx", "2023-07"),
    ("АВГУСТ.xlsx", "2023-08"),
    ("СЕНТЯБРЬ.xlsx", "2023-09"),
    # "СЕНТЯБРЬ (2).xlsx" is excluded: it duplicates part of the September
    # extract (identical HGB/MCV/PLT values in both files).
    ("Октябрь.xlsx", "2023-10"),
    ("Ноябрь.xlsx", "2023-11"),
    ("Декабрь.xlsx", "2023-12"),
    ("Январь 2024.xlsx", "2024-01"),
]

ALIASES = {
    "reg_id": ["REG_NO", "PROTOKOL_NO"],
    "date":   ["DATA_REG", "KAYITTARIHI"],
    "sex":    ["POL", "CINSIYET"],
    "age":    ["VOZRAST", "YEAR", "YAS"],
    "test":   ["ANALIZ", "TEST_ADI"],
    "result": ["REZULTAT", "SONUC"],
}
PII_PATTERN = re.compile(r"FIO|ADRES|TEL|SOYAD|NAME", re.IGNORECASE)

SEX_MAP = {
    "Man": "Man", "M": "Man", "Erkek": "Man",
    "Woman": "Woman", "F": "Woman", "Kadın": "Woman", "Kadin": "Woman",
}

NUM_RE = re.compile(r"^-?[0-9]+([.,][0-9]+)?$")
STAR_RE = re.compile(r"^\*+$")

# ---- District classification (administrative geography, not PII) ----
# Matches word stems (e.g. "Чиланзар" also catches "Чиланзарский") against
# the raw address text. This is the only use of the address field: the raw
# text is not stored, printed or passed on; only the matched category is.
#
# Limitation: do not use AREA/DISTRICT as a research variable. About 86% of
# addresses are entered in an operator-dependent shorthand (e.g. "Ю 2 3 5")
# that, as the data owner confirmed, is used inconsistently for both city
# and region, so not even a city/region split can be recovered for them.
# Only the ~14% written out in full are classified, and that subset is not
# a random sample. AREA/DISTRICT is exploratory metadata only.
TASHKENT_CITY_DISTRICTS = {
    "Bektemir": ["Бектемир"],
    "Chilanzar": ["Чиланзар"],
    "Mirzo-Ulugbek": ["Мирзо.?Улугбек", "Мирзо.?Улуг.?бек"],
    "Mirabad": ["Мирабад"],
    "Almazar/Olmazor": ["Алмазар"],
    "Sergeli": ["Сергели"],
    "Shaykhantakhur": ["Шайхантах", "Шайхантаур"],
    "Uchtepa": ["Учтепа"],
    "Yakkasaray": ["Яккасарай"],
    "Yangihayot": ["Янгихаёт", "Янгихает"],
    "Yashnobod": ["Яшнобод", "Хамзинск"],  # Хамзинский = name before renaming
    # Yunusabad is matched separately by _YUNUSABAD_RE below, which also
    # catches abbreviated forms (ю.обод/юобод).
}
TASHKENT_REGION_DISTRICTS = {
    "Bekabad": ["Бекабад"],
    "Bostanliq": ["Бостанлык", "Бустонлик"],
    "Buka": ["Букинск", r"\bБука\b"],
    "Chinaz": ["Чиназ"],
    "Qibray": ["Кибрай"],
    "Ohangaron": ["Ахангаран"],
    "Oqqorgon": ["Аккурган"],
    "Parkent": ["Паркент"],
    "Piskent": ["Пскент"],
    "Tashkent_district": ["Ташкентский р", "Ташкентский район"],
    "Yuqori_Chirchiq": ["Верхнечирчик", "Верхний Чирчик"],
    "Quyi_Chirchiq": ["Нижнечирчик", "Нижний Чирчик"],
    "Zangiota": ["Зангиат"],
    "Yangiyol": ["Янгиюль"],
    "Angren": ["Ангрен"],
    "Almalyk": ["Алмалык"],
    "Chirchiq_city": [r"\bЧирчик\b"],
}


def _translit_vowel(text):
    """Collapse о to а (Юнусобод/Юнусабад). Safe for both address text and
    regex patterns, since 'о' has no regex meaning."""
    return text.lower().replace("о", "а")


def _normalize_translit(text):
    """Normalise raw address text (not patterns; see _translit_vowel):
    collapse о to а and remove dots and hyphens ("ю.обод" -> "юобод").
    Spaces are kept because some patterns use \\b word boundaries (e.g.
    Чирчик, so it does not match inside a longer word). Abbreviations
    without spaces are handled by _strip_spaces() where needed.
    """
    t = _translit_vowel(text)
    t = re.sub(r"[.\-]", "", t)
    return t


def _strip_spaces(text):
    return re.sub(r"\s", "", text)


# Patterns get only the vowel collapse: removing dots and hyphens, as
# _normalize_translit does for address text, would break regex syntax
# such as ".?".
_CITY_RE = {
    k: re.compile("|".join(_translit_vowel(p) for p in v))
    for k, v in TASHKENT_CITY_DISTRICTS.items()
}
_REGION_RE = {
    k: re.compile("|".join(_translit_vowel(p) for p in v))
    for k, v in TASHKENT_REGION_DISTRICTS.items()
}
_TASHKENT_CITY_GENERIC = re.compile(_translit_vowel("г. Ташкент") + "|" + _translit_vowel("город Ташкент"))
# Matches full forms ("юнусабад", "юнусобод", "юнус обод") and abbreviated
# ones ("ю.обод", "юобод"). Checked against the text with spaces removed,
# so "ю обод" and "ю.обод" both become "юобод".
_YUNUSABAD_RE = re.compile(r"^ю.*абад")
# Ц1-Ц5 (Yunusabad's Ц-quarter blocks) and Чинобод/Чинабад (Chinobod massif)
# are sub-areas of Yunusabad, as confirmed by the data owner.
_YUNUSABAD_EXTRA_RE = re.compile(r"^ц[1-5]\b|чинабад")


def classify_district(address_raw):
    """Returns (area, district) - area is 'City'/'Region'/'Unclassified'."""
    if not address_raw:
        return "Unclassified", None
    normalized = _normalize_translit(address_raw)
    if _YUNUSABAD_RE.search(_strip_spaces(normalized)) or _YUNUSABAD_EXTRA_RE.search(normalized):
        return "City", "Yunusabad"
    for name, pattern in _CITY_RE.items():
        if pattern.search(normalized):
            return "City", name
    for name, pattern in _REGION_RE.items():
        if pattern.search(normalized):
            return "Region", name
    if _TASHKENT_CITY_GENERIC.search(normalized):
        return "City", "Unspecified_city_district"
    return "Unclassified", None


def classify_result(raw, test_std=None):
    """Classify a raw result cell.

    test_std is the standardised parameter name (HGB, MCV, PLT, RBC, ...);
    it is used only to decide whether a zero result is implausible.

    The flag patterns below are copied from the distinct values in the
    source files, read as UTF-8. Console output of this data is often
    garbled, so do not write patterns from memory or from console output:
    list the distinct values from the files and copy the exact strings.
    """
    if raw is None:
        return "Missing", None
    x = raw.strip()
    if x == "":
        return "Missing", None
    # "Сгусток - L05" = clotted sample. The most common failure (216 readings
    # = 72 draws x 3 parameters). Clotting consumes platelets, so including
    # these would bias PLT downward.
    if re.search("Сгусток", x, re.IGNORECASE):
        return "Clotted_sample", None
    # "Недостаточно биоматериала - L06" = insufficient biomaterial
    if re.search("Недостаточно биоматериала", x, re.IGNORECASE):
        return "Insufficient_biomaterial", None
    # "Не правильно зобор" (sic; "забор" is misspelt in the source) = improperly
    # collected
    if re.search("правильно зобор", x, re.IGNORECASE):
        return "Improper_collection", None
    # "Отказ клиента - L07" = patient declined, no sample taken
    if re.search("Отказ клиента", x, re.IGNORECASE):
        return "Client_refusal", None
    # "Не доставлен материал - L01" = sample not delivered to the lab
    if re.search("доставлен материал", x, re.IGNORECASE):
        return "Sample_not_delivered", None
    # "анализатор не смог подсчитать лейкоцитарную форму"
    if re.search("анализатор", x, re.IGNORECASE):
        return "Analyzer_error", None
    # bare "*" / "****" - the analyser could not calculate a value
    if STAR_RE.match(x):
        return "Analyzer_error", None
    if NUM_RE.match(x):
        val = float(x.replace(",", "."))
        # A zero is a failed measurement for red cell, platelet and total
        # counts, but a real finding for the differential (see
        # ZERO_IS_CLINICALLY_VALID).
        if val == 0 and (test_std is None or test_std not in ZERO_IS_CLINICALLY_VALID):
            return "Zero_implausible", None
        return "Valid", val
    return "Other_invalid", None


# Cyrillic letters that look identical to Latin ones but are different
# characters. This data's "(МСV)" label uses Cyrillic М and С with a Latin V,
# which a plain ASCII regex does not match. Test names are normalised before
# matching so the letters are recognised in either script.
_HOMOGLYPHS = str.maketrans({
    "Н": "H", "В": "B", "М": "M", "С": "C", "Р": "P", "Т": "T",
})


# The CBC panel reported by the Mindray BC-6800Plus in this data. Each test
# is matched on the Latin abbreviation in parentheses in the test name
# (e.g. "Гемоглобин (HGB)"), after homoglyph normalisation.
#
# Some abbreviations contain others (MCHC contains MCH; RDW-CV and RDW-SD
# contain RDW; each "#" differential has a "%" partner), so every pattern
# matches the full parenthesised token.
PANEL = [
    ("HGB",      r"\(HGB\)"),
    ("RBC",      r"\(RBC\)"),
    ("HCT",      r"\(HCT\)"),
    ("MCV",      r"\(MCV\)"),
    ("MCH",      r"\(MCH\)"),
    ("MCHC",     r"\(MCHC\)"),
    ("RDW_CV",   r"\(RDW-CV\)"),
    ("RDW_SD",   r"\(RDW-SD\)"),
    ("PLT",      r"\(PLT\)"),
    ("MPV",      r"\(MPV\)"),
    ("PCT",      r"\(PCT\)"),
    ("PDW",      r"\(PDW\)"),
    ("WBC",      r"\(WBC\)"),
    ("NEU_ABS",  r"\(NEU#\)"),
    ("NEU_PCT",  r"\(NEU%\)"),
    ("LYMPH_ABS", r"\(LYMPH#\)"),
    ("LYMPH_PCT", r"\(LYMPH%\)"),
    ("MON_ABS",  r"\(MON#\)"),
    ("MON_PCT",  r"\(MON%\)"),
    ("EOS_ABS",  r"\(EOS#\)"),
    ("EOS_PCT",  r"\(EOS%\)"),
    ("BAS_ABS",  r"\(BAS#\)"),
    ("BAS_PCT",  r"\(BAS%\)"),
    ("IMG_ABS",  r"\(IMG#\)"),
    ("IMG_PCT",  r"\(IMG%\)"),
    ("NRBC_ABS", r"\(NRBC#\)"),
]
_PANEL_RE = [(name, re.compile(pat)) for name, pat in PANEL]

# Only these three parameters have plausibility bounds in this script (see
# below). The other 23 have pre-analytical flags and implausible zeros
# removed but no range limits here, so implausible values remain in this
# script's output, for example:
#   MCH      up to 241 pg     (physiological range roughly 27-33)
#   MCHC     up to 887 g/L    (physiological range roughly 320-360)
#   PCT      up to 325        (physiological range roughly 0.1-0.5;
#                              these look like platelet counts written
#                              into the plateletcrit field)
#   RBC      down to 0.01     (incompatible with life)
#   RDW_SD   up to 149.5
# validate_panel.py removes these in the next step, by internal
# consistency checks and physiological limits for every parameter.
CORE_PARAMS = {"HGB", "MCV", "PLT"}

# Parameters for which zero is a real result, not a failed measurement: the
# differential. Eosinophils, basophils, immature granulocytes and nucleated
# red cells are often zero in healthy people; neutrophils, lymphocytes and
# monocytes can reach zero in agranulocytosis, severe lymphopenia and
# monocytopenia. No other parameter (HGB, RBC, HCT, MCV, MCH, MCHC, RDW,
# PLT, MPV, PCT, PDW, WBC) can be zero in a successfully analysed sample,
# so a zero there means a failed run.
ZERO_IS_CLINICALLY_VALID = {
    "NEU_ABS", "NEU_PCT", "LYMPH_ABS", "LYMPH_PCT", "MON_ABS", "MON_PCT",
    "EOS_ABS", "EOS_PCT", "BAS_ABS", "BAS_PCT", "IMG_ABS", "IMG_PCT",
    "NRBC_ABS",
}


def match_test(name):
    if not name:
        return None
    normalized = name.upper().translate(_HOMOGLYPHS)
    for std_name, pattern in _PANEL_RE:
        if pattern.search(normalized):
            return std_name
    return None


class XlsxReader:
    """Minimal xlsx reader: shared strings + first worksheet, streamed by row."""

    def __init__(self, path):
        self.zf = zipfile.ZipFile(path)
        sst_raw = self.zf.read("xl/sharedStrings.xml").decode("utf-8", errors="ignore")
        self.strings = re.findall(r"<t[^>]*>(.*?)</t>", sst_raw, re.S)
        self.sheet_xml = self.zf.read("xl/worksheets/sheet1.xml").decode("utf-8", errors="ignore")

    def _resolve(self, cell_xml):
        if cell_xml is None:
            return None
        is_shared = 't="s"' in cell_xml
        v = re.search(r"<v>(.*?)</v>", cell_xml)
        if not v:
            return None
        val = v.group(1)
        if is_shared:
            idx = int(val)
            val = self.strings[idx] if idx < len(self.strings) else val
        return val

    def header_colmap(self):
        m = re.search(r'<row r="1"[^>]*>(.*?)</row>', self.sheet_xml, re.S)
        header_xml = m.group(1)
        colmap = {}
        for col_m in re.finditer(r'<c r="([A-Z]+)1"[^>]*(?:/>|>.*?</c>)', header_xml):
            col = col_m.group(1)
            full = re.search(
                r'<c r="' + col + r'1"[^/]*/>|<c r="' + col + r'1"[^>]*>.*?</c>',
                header_xml,
            )
            val = self._resolve(full.group(0)) if full else None
            if val:
                colmap[val] = col
        return colmap

    def rows(self):
        for row_m in re.finditer(r'<row r="\d+"[^>]*>(.*?)</row>', self.sheet_xml, re.S):
            yield row_m.group(1)

    def cell(self, row_xml, col):
        m = re.search(r'<c r="' + col + r'\d+"[^>]*>(.*?)</c>', row_xml, re.S)
        return self._resolve(m.group(0)) if m else None


def find_col(colmap, field):
    safe_headers = {h: c for h, c in colmap.items() if not PII_PATTERN.search(h)}
    for alias in ALIASES[field]:
        for header, col in safe_headers.items():
            if header.upper() == alias.upper():
                return col
    return None


ADDRESS_ALIASES = ["ADRES", "ADRES1"]


def find_address_col(colmap):
    # The one exception to the PII filter: the address is needed briefly to
    # derive a district category. The raw value is read, classified and
    # discarded in the same loop iteration; it is not stored, printed or
    # written to any file.
    for alias in ADDRESS_ALIASES:
        for header, col in colmap.items():
            if header.upper() == alias.upper():
                return col
    return None


def load_month(path, month_label, flag_counter, sex_drop_counter, district_counter,
               age_drop_counter):
    reader = XlsxReader(path)
    colmap = reader.header_colmap()

    cols = {f: find_col(colmap, f) for f in ALIASES}
    missing = [f for f, c in cols.items() if c is None]
    if missing:
        print(f"  WARNING {month_label}: missing columns for {missing} -- colmap={colmap}")
        return []

    addr_col = find_address_col(colmap)
    if addr_col is None:
        print(f"  WARNING {month_label}: no address column found - district will be Unclassified")

    long_rows = []
    first_row = True
    for row_xml in reader.rows():
        if first_row:
            first_row = False
            continue  # header row already parsed separately

        reg_id = reader.cell(row_xml, cols["reg_id"])
        date = reader.cell(row_xml, cols["date"])
        sex_raw = reader.cell(row_xml, cols["sex"])
        age_raw = reader.cell(row_xml, cols["age"])
        test_raw = reader.cell(row_xml, cols["test"])
        result_raw = reader.cell(row_xml, cols["result"])

        test_std = match_test(test_raw)
        if test_std is None:
            continue

        # The address is classified immediately and only `area`/`district`
        # are kept; address_raw goes out of scope at the end of this
        # iteration. Done only for rows that matched a CBC parameter.
        address_raw = reader.cell(row_xml, addr_col) if addr_col else None
        area, district = classify_district(address_raw)
        district_counter[(area, district)] += 1

        sex = SEX_MAP.get((sex_raw or "").strip())
        if sex is None:
            sex_drop_counter[sex_raw] += 1
            continue

        try:
            age = float(age_raw)
        except (TypeError, ValueError):
            continue

        # The oldest verified human age is 122 (Jeanne Calment); higher values
        # are data-entry errors. "124" recurs several times in this data,
        # which suggests a placeholder date of birth.
        if age > 122:
            age_drop_counter[age] += 1
            continue

        # Infants under 1 year are excluded. Age is recorded only in whole
        # years, so "0" covers 0-11 months, and infants under 6 months (for
        # whom WHO haemoglobin thresholds are not defined) cannot be told
        # apart from older ones. This also excludes 6-11-month-olds. If date
        # of birth becomes available, compute age in months and lower the
        # cut-off to 6 months.
        if age < 1:
            age_drop_counter["<1y (WHO threshold undefined <6mo)"] += 1
            continue

        flag, val = classify_result(result_raw, test_std)

        # Plausibility bounds, set in advance from the extremes reported in
        # published case reports. result_exclusion_table.csv shows how many
        # readings each bound removes.
        #
        # MCV 44.5-136.9 fL, checked against the primary sources:
        #   lower 44.5 fL - Goretti L, Adiatmaja CO, Kahar H. Severe
        #     microcytosis in a hemoglobin E/beta-thalassemia patient with
        #     signs of iron deficiency: a case report. Ann Med Surg.
        #     2022;78:103826. (PMC9207008)
        #   upper 136.9 fL - Akinola OO, Mandapati A, Douthit N. Unveiling
        #     the pernicious truth: a case report on the rare presentation
        #     of severe vitamin B12 deficiency. Cureus. 2026;18(1):e101860.
        #     (PMC12916071)
        # Hundreds of readings in this data lie between 50 and 60 fL,
        # consistent with the high burden of microcytic anaemia, so a
        # narrower range such as 50-130 would exclude real results.
        #
        # PLT < 1 x10^9/L: published severe thrombocytopenia is reported in
        # whole numbers (5-8 x10^9/L), never below 1. Values below 1 here are
        # more likely another parameter, such as plateletcrit (~0.1-0.5),
        # in the wrong field.
        # PLT > 1703 x10^9/L, checked against the primary source, which
        # reports it as the highest platelet count in iron deficiency to
        # date: Voigt W, Jordan K, Sippel C, Amoury M, Schmoll HJ, Wolf HH.
        # Severe thrombocytosis and anemia associated with celiac disease in
        # a young female patient: a case report. J Med Case Rep. 2008;2:96.
        # (PMC2329657). Iron deficiency is also the relevant mechanism here.
        # The bound removes one value in this data (1945).
        #
        # HGB < 14 g/L: at the lower edge of values reported in any surviving
        # patient, so such a value in an outpatient sample is more likely an
        # error. The published cases at this level were in shock, not
        # ambulatory. No value in this data falls below 35 g/L, so the rule
        # removes nothing.
        #
        # No upper HGB bound: the maximum here (211 g/L) is near the threshold
        # for severe polycythaemia in men (>210 g/L) and is not implausible.
        if flag == "Valid" and test_std == "MCV" and not (44.5 <= val <= 136.9):
            flag = "Implausible_MCV"
        if flag == "Valid" and test_std == "PLT" and not (1 <= val <= 1703):
            flag = "Implausible_PLT"
        if flag == "Valid" and test_std == "HGB" and val < 14:
            flag = "Implausible_HGB"

        flag_counter[flag] += 1
        if flag != "Valid":
            continue

        long_rows.append((reg_id, date, sex, age, test_std, val, month_label, area, district))

    return long_rows


def main():
    flag_counter = Counter()
    sex_drop_counter = Counter()
    district_counter = Counter()
    age_drop_counter = Counter()
    all_long = []

    for fn, label in FILES:
        path = f"{RAW_DIR}/{fn}"
        print(f"Loading {fn} ...")
        rows = load_month(path, label, flag_counter, sex_drop_counter, district_counter,
                          age_drop_counter)
        print(f"  -> {len(rows)} valid CBC readings kept")
        all_long.extend(rows)

    print("\n=== District classification match rate (CBC readings) ===")
    total_district = sum(district_counter.values())
    for (area, district), n in district_counter.most_common():
        print(f"  {area}/{district}: {n} ({100*n/total_district:.1f}%)")

    print(f"\nTotal valid long-format readings: {len(all_long)}")

    print("\n=== Rows dropped due to unmapped/corrupted sex values ===")
    for k, v in sex_drop_counter.most_common():
        print(f"  {k!r}: {v}")

    print("\n=== Readings dropped on age grounds (implausible >122y; infants <1y) ===")
    for k, v in age_drop_counter.most_common():
        print(f"  age {k}: {v} readings")

    print("\n=== Result-quality flag breakdown (all CBC readings, incl. excluded) ===")
    for k, v in flag_counter.most_common():
        print(f"  {k}: {v}")
    with open(f"{RESULT_DIR}/result_exclusion_table.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["result_flag", "n"])
        for k, v in flag_counter.most_common():
            w.writerow([k, v])

    # ---- pivot to wide: one row per (reg_id, date) patient-test-event ----
    groups = defaultdict(dict)
    meta = {}
    for reg_id, date, sex, age, test_std, val, month, area, district in all_long:
        key = (reg_id, date, month)
        groups[key][test_std] = val
        meta[key] = (sex, age, month, area, district)

    print(f"\nPatient-test-events (unique reg_id+date+month): {len(groups)}")

    both_hgb_mcv_plt = sum(
        1 for k, v in groups.items() if "HGB" in v and "MCV" in v and "PLT" in v
    )
    print(f"Events with all of HGB/MCV/PLT present: {both_hgb_mcv_plt}")

    # ---- remove duplicate records ----
    # The same SEX+AGE+HGB+MCV+PLT combination recurs up to 16 times, on up
    # to 11 dates months apart. Repeat testing cannot reproduce MCV to 0.1 fL
    # months later, so these are copies of one measurement. Keeping them
    # would inflate N and give that measurement extra weight.
    #
    # Chance matches were ruled out: about 5 exact five-way matches would be
    # expected by chance among ~32,000 records, against ~890 groups observed.
    # The first occurrence of each combination is kept.
    deduped = {}
    dup_removed = 0
    for key, vals in groups.items():
        reg_id, date, month = key
        sex, age, _, area, district = meta[key]
        sig = (sex, age, vals.get("HGB"), vals.get("MCV"), vals.get("PLT"))
        if sig in deduped:
            dup_removed += 1
            continue
        deduped[sig] = (key, vals)

    print(f"\nDuplicate records removed: {dup_removed}")
    print(f"Final unique patient-test-events: {len(deduped)}")

    # ---- completeness of each parameter across the final records ----
    param_counts = Counter()
    for key, vals in deduped.values():
        for p in vals:
            param_counts[p] += 1
    print("\n=== Parameter completeness across final events ===")
    n_final = len(deduped)
    for name, _ in PANEL:
        n = param_counts.get(name, 0)
        print(f"  {name:10s} {n:6d}  ({100*n/n_final:5.1f}%)")

    param_names = [name for name, _ in PANEL]

    out_path = f"{DATA_DIR}/cbc_cohort.csv"
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["ID", "SEX", "AGE"] + param_names + ["month", "date", "AREA", "DISTRICT"])
        for i, (key, vals) in enumerate(deduped.values(), start=1):
            reg_id, date, month = key
            sex, age, _, area, district = meta[key]
            w.writerow([i, sex, age]
                       + [vals.get(p, "") for p in param_names]
                       + [month, date, area, district or ""])

    print(f"\nSaved: {out_path}")
    print(f"Saved: {RESULT_DIR}/result_exclusion_table.csv")


if __name__ == "__main__":
    main()
