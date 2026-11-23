"""Extract Tables 2-5 (Cases 0-3) of the IFOP 2022 hake status report (pdftotext -layout text).

Source: Tascheri R. (2022) Segundo Informe Tecnico, Convenio de Desempeno 2021, Estatus y posibilidades
de explotacion biologicamente sustentables de los principales recursos pesqueros nacionales, ano 2022:
Merluza comun. IFOP. Tables 2-5, printed pp. 54-57: total biomass 2+ (kt), spawning biomass (kt),
exploitable biomass (kt), recruits (10^6), fishing mortality F and exploitation rate u, 1997-2021.
Usage: python3 extract_case_tables.py <report.pdf> <out.csv>
"""
import csv, re, subprocess, sys

pdf, out = sys.argv[1], sys.argv[2]
txt = subprocess.run(["pdftotext", "-layout", pdf, "-"], capture_output=True, text=True, check=True).stdout
rows = []
for case, table in enumerate(["Tabla 2.", "Tabla 3.", "Tabla 4.", "Tabla 5."]):
    starts = [m.start() for m in re.finditer(re.escape(table), txt)]
    # the first occurrence after the list of tables whose caption ends with "Caso <case>."
    block = None
    for s in starts:
        head = txt[s:s + 600]
        if re.search(r"1997 a 2021\.\s*Caso %d\." % case, head):
            block = txt[s:s + 6000]
            break
    if block is None:
        sys.exit("table for case %d not found" % case)
    n = 0
    for line in block.splitlines():
        m = re.match(r"\s*(19\d\d|20\d\d)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s*$", line)
        if m:
            y, bt, bd, bv, rec, F, u = m.groups()
            rows.append(dict(case=case, year=int(y), total_biomass_kt=int(bt), ssb_kt=int(bd),
                             exploitable_biomass_kt=int(bv), recruits_1e6=int(rec), F=float(F), u=float(u)))
            n += 1
        if n == 25:
            break
    if n != 25:
        sys.exit("case %d: %d rows instead of 25" % (case, n))
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)
print("wrote", out, len(rows), "rows")
