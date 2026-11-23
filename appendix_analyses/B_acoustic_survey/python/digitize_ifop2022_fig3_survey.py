"""Digitise the acoustic-survey biomass of Chilean hake from Figure 3 of the IFOP 2022 status report
(Tascheri 2022, Segundo Informe, printed p. 16, PDF p. 21; stacked bars by survey zone, 1995-2021).

Only needed for the survey years without a printed value: 1995, 1997, 1999, 2000, 2001. Years 2002 and
2004-2024 are taken from Table 10 of the 2024 survey report (IFOP P-682119, published 2025), which is also
used to check the digitisation.

Method: the embedded figure (911 x 810 px) is extracted with pdfimages. The y axis is calibrated with
the tick marks of the left axis (0, 500, 1000, 1500 thousand t); each bar is located on a row just above
the baseline as a run of non-white pixels, and its height is the topmost non-white pixel of the run's
central columns, found by walking up from the baseline to the first white pixel (so the legend,
which sits above some bars, is not picked up). Usage:
  python3 -I digitize_ifop2022_fig3_survey.py <IFOP_2022_report.pdf> <table10_2024.csv> <out.csv>
"""
import csv, subprocess, sys, tempfile, os
import numpy as np
from PIL import Image

pdf, t10_path, out = sys.argv[1:4]
YEARS = [1995, 1997, 1999, 2000, 2001, 2002] + list(range(2004, 2022))
with tempfile.TemporaryDirectory() as tmp:
    subprocess.run(["pdfimages", "-f", "21", "-l", "21", "-png", pdf, os.path.join(tmp, "f")], check=True)
    im = np.array(Image.open(os.path.join(tmp, "f-001.png")).convert("RGB")).astype(int)
assert im.shape[:2] == (810, 911), im.shape

dark = im.sum(2) < 200
axis_x = int(np.argmax(dark[:, :150].sum(0)))
ticks = np.where(dark[:, axis_x - 3])[0]
ticks = [t for i, t in enumerate(ticks) if i == 0 or t - ticks[i - 1] > 2]
assert len(ticks) == 4, ticks                        # 1500, 1000, 500, 0
y_tick = np.array(ticks[::-1], float); v_tick = np.array([0, 500, 1000, 1500], float)
slope, icpt = np.polyfit(y_tick, v_tick, 1)          # thousand t per pixel (negative), linear

white = im.min(2) > 235
row = int(y_tick[0]) - 3                             # just above the zero line
cols = ~white[row, axis_x + 3:]
runs, start = [], None
for i, c in enumerate(cols):
    if c and start is None: start = i
    if not c and start is not None: runs.append((start + axis_x + 3, i - 1 + axis_x + 3)); start = None
runs = [r for r in runs if r[1] - r[0] > 15]
assert len(runs) == len(YEARS), (len(runs), runs)
vals = []
for (a, b) in runs:
    centre = slice(a + 5, b - 4)
    tops = []
    for c in range(centre.start, centre.stop):
        y = row
        while y > 0 and not white[y - 1, c]: y -= 1
        tops.append(y)
    top = min(tops)
    vals.append(round(icpt + slope * top, 1))

t10 = {int(r["year"]): float(r["biomass_t"]) / 1000 for r in csv.DictReader(open(t10_path))}
diff = [v - t10[y] for y, v in zip(YEARS, vals) if y in t10]
print("check vs Table 10 (2002, 2004-2021): n = %d, bias = %.1f kt, mean abs = %.1f kt, max abs = %.1f kt, r = %.4f" % (
    len(diff), np.mean(diff), np.mean(np.abs(diff)), np.max(np.abs(diff)),
    np.corrcoef([v for y, v in zip(YEARS, vals) if y in t10], [t10[y] for y in YEARS if y in t10])[0, 1]))
with open(out, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["year", "biomass_kt_digitised", "biomass_kt_table10_2024"])
    for y, v in zip(YEARS, vals): w.writerow([y, v, t10.get(y, "")])
print("wrote", out)
