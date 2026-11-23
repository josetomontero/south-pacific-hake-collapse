"""Extract the historical acoustic-survey biomass of Chilean hake (Table 10 of the 2024 direct-assessment
report, IFOP project P-682119, published 2025: biomass (t) of the winter surveys 2002 and 2004-2024).
Usage: python3 -I extract_table10_2024.py <report.pdf> <out.csv>
"""
import csv, re, subprocess, sys
pdf, out = sys.argv[1:3]
txt = subprocess.run(["pdftotext", "-layout", pdf, "-"], capture_output=True, text=True, check=True).stdout
i = txt.index("Tabla 10      Composición porcentual de las estimaciones de biomasa")
block = txt[i:i + 8000].splitlines()
years, bio = [], []
for line in block:
    if re.match(r"\s*Años\s", line):
        years += [int(y) for y in re.findall(r"\b(20\d\d)\b", line)]
    m = re.match(r"\s*B \(t\)\s+(.*)$", line)
    if m:
        bio += [int(v.replace(".", "")) for v in m.group(1).split()]
assert len(years) == len(bio) == 22, (years, bio)
with open(out, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["year", "biomass_t"]); w.writerows(zip(years, bio))
print("wrote", out, len(years), "years")
