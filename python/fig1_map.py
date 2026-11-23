#!/usr/bin/env python3
"""Figure 1: study area of the South Pacific hake fishery.

Cartography follows the style of the SDM_marino_extendido atlas maps (Albers equal-area
projection, paper and ink colours, 5-degree graticule). Graticule labels are placed outside the
frame as axis ticks; the South America inset sits inside the map frame.

Layers (all public or produced by the R pipeline):
  data/external/natural_earth_paises_10m.geojson   Natural Earth 1:10m countries (public domain)
  data/derived/LME_Humboldt.geojson                Humboldt Current LME (NOAA LME programme)
  data/derived/fishery_unit_polygon.geojson        hake fishery management unit (Tascheri 2022),
                                                   built by R/00_study_area_polygons.R
  figures/hake_illustration_Salinas_Nass.png       illustration by Pablo Salinas Nass (permission)

Environment: python >= 3.9 with matplotlib, numpy, pyproj, shapely (no GDAL needed).
Usage:  python python/fig1_map.py [--out results/figs]
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon as MplPolygon, Rectangle
import numpy as np
from pyproj import Transformer
from shapely.geometry import shape, box
from shapely.ops import transform as sh_transform

ROOT = Path(__file__).resolve().parents[1]
EXT, DER, FIG = ROOT / 'data/external', ROOT / 'data/derived', ROOT / 'figures'

CRS = ('+proj=aea +lat_1=-25 +lat_2=-42 +lat_0=-33 +lon_0=-74 '
       '+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
FWD = Transformer.from_crs('EPSG:4326', CRS, always_xy=True)
INK, MUTED, PAPER, SEA = '#19333d', '#53666d', '#ffffff', '#eef3f5'
LAND, LAND_EDGE = '#dddcd0', '#a8ada4'
LME_COL, UNIT_COL, UNIT_EDGE = '#c9dff0', '#2f5f9e', '#1d3f6e'
GRAT = '#819c9b'

LON0, LON1, LAT0, LAT1 = -84.0, -66.0, -45.0, -20.0          # map window (degrees)
WIDTH_CM, HEIGHT_CM = 14.0, 19.0


def proj_geom(g):
    return sh_transform(lambda x, y, z=None: FWD.transform(x, y), g)


def polys(g):
    return list(g.geoms) if g.geom_type == 'MultiPolygon' else [g]


def add_geojson(ax, path, **kw):
    feats = json.loads(Path(path).read_text())['features']
    for f in feats:
        g = shape(f['geometry'])
        if g.is_empty or not g.intersects(box(-95, -60, -30, 15)):
            continue
        for p in polys(proj_geom(g)):
            ax.add_patch(MplPolygon(np.asarray(p.exterior.coords), closed=True, **kw))
            for ring in p.interiors:
                ax.add_patch(MplPolygon(np.asarray(ring.coords), closed=True, facecolor=SEA,
                                        edgecolor=kw.get('edgecolor', 'none'), lw=kw.get('lw', 0),
                                        zorder=kw.get('zorder', 1)))


def graticule_with_ticks(ax, xmin, xmax, ymin, ymax, step=5):
    """Light graticule inside the frame; degree labels as outward ticks on the frame."""
    xt, xl, yt, yl = [], [], [], []
    for la in range(-60, 0, step):
        lons = np.linspace(-100, -40, 600)
        gx, gy = FWD.transform(lons, np.full_like(lons, la))
        ax.plot(gx, gy, color=GRAT, lw=.3, alpha=.5, zorder=5)
        yv = np.interp(xmin, gx, gy)                 # where the parallel meets the left edge
        if ymin < yv < ymax:
            yt.append(yv); yl.append(f'{abs(la)}°S')
    for lo in range(-100, -40, step):
        las = np.linspace(-70, 10, 600)
        gx, gy = FWD.transform(np.full_like(las, lo), las)
        ax.plot(gx, gy, color=GRAT, lw=.3, alpha=.5, zorder=5)
        xv = np.interp(ymin, gy, gx)                 # where the meridian meets the bottom edge
        if xmin < xv < xmax:
            xt.append(xv); xl.append(f'{abs(lo)}°W')
    ax.set_xticks(xt); ax.set_xticklabels(xl)
    ax.set_yticks(yt); ax.set_yticklabels(yl)
    ax.tick_params(direction='out', length=3, width=.6, colors=INK, labelsize=8, pad=2)


def geodetic_rect(lo0, lo1, la0, la1):
    ring = [(lo, la0) for lo in np.linspace(lo0, lo1, 20)] + [(lo1, la) for la in np.linspace(la0, la1, 40)] + \
           [(lo, la1) for lo in np.linspace(lo1, lo0, 20)] + [(lo0, la) for la in np.linspace(la1, la0, 40)]
    rx, ry = FWD.transform([p[0] for p in ring], [p[1] for p in ring])
    return np.column_stack([rx, ry])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=str(ROOT / 'results/figs'))
    args = ap.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'text.color': INK,
                         'axes.labelcolor': INK, 'pdf.fonttype': 42})
    fig = plt.figure(figsize=(WIDTH_CM / 2.54, HEIGHT_CM / 2.54), facecolor=PAPER)

    cx, cy = FWD.transform([LON0, LON1, LON0, LON1], [LAT0, LAT0, LAT1, LAT1])
    xmin, xmax, ymin, ymax = min(cx), max(cx), min(cy), max(cy)
    ax = fig.add_axes([.09, .05, .89, .93], facecolor=SEA)
    ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax); ax.set_aspect('equal')
    for sp in ax.spines.values():
        sp.set_edgecolor(INK); sp.set_linewidth(.7)

    add_geojson(ax, DER / 'LME_Humboldt.geojson', facecolor=LME_COL, edgecolor='none', lw=0, zorder=0.5)
    add_geojson(ax, EXT / 'natural_earth_paises_10m.geojson', facecolor=LAND, edgecolor=LAND_EDGE, lw=.35, zorder=1)
    add_geojson(ax, DER / 'fishery_unit_polygon.geojson', facecolor=UNIT_COL, edgecolor=UNIT_EDGE, lw=.5, alpha=.6, zorder=3)
    graticule_with_ticks(ax, xmin, xmax, ymin, ymax)

    # country labels: horizontal, centred on the east-west chord of the country at the chosen latitude,
    # and checked to lie entirely inside the country polygon
    countries = {f['properties'].get('ADMIN') or f['properties'].get('name') or f['properties'].get('NAME'): shape(f['geometry'])
                 for f in json.loads((EXT / 'natural_earth_paises_10m.geojson').read_text())['features']}
    from shapely.geometry import LineString, Polygon as ShPolygon
    for lab, country, la, fs in [('CHILE', 'Chile', -35.0, 7), ('ARGENTINA', 'Argentina', -34.5, 7.5)]:
        geom = countries[country]
        chord = geom.intersection(LineString([(LON0 + 1.0, la), (LON1 - 0.8, la)]))   # chord clipped to the map window
        xs = [c[0] for c in (chord.coords if chord.geom_type == 'LineString' else max(chord.geoms, key=lambda g: g.length).coords)]
        xv, yv = FWD.transform((min(xs) + max(xs)) / 2, la)
        txt = ax.text(xv, yv, lab, fontsize=fs, color=MUTED, ha='center', va='center', zorder=6)
        fig.canvas.draw()
        bb = txt.get_window_extent(fig.canvas.get_renderer()).transformed(ax.transData.inverted())
        inv = Transformer.from_crs(CRS, 'EPSG:4326', always_xy=True)
        corners = [inv.transform(x, y) for x, y in [(bb.x0, bb.y0), (bb.x1, bb.y0), (bb.x1, bb.y1), (bb.x0, bb.y1)]]
        assert geom.contains(ShPolygon(corners)), f'{lab} label is not fully inside {country}'
    xv, yv = FWD.transform(-81.0, -33.0)
    ax.text(xv, yv, 'PACIFIC OCEAN', rotation=90, ha='center', va='center', fontsize=8.5, color='#749499', zorder=6)
    xv, yv = FWD.transform(-76.0, -36.8)
    ax.text(xv, yv, 'Humboldt Current\nLarge Marine Ecosystem', rotation=90, ha='center', va='center',
            fontsize=7.5, color='#3f6e9c', zorder=6, linespacing=1.3)

    # legend in the upper-left corner (no frame), illustration right below it, left-aligned
    handles = [Patch(facecolor=LME_COL, edgecolor='none', label='Humboldt Current LME'),
               Patch(facecolor=UNIT_COL, edgecolor=UNIT_EDGE, alpha=.6, label='Hake fishery management unit')]
    leg = ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(.012, .992), frameon=False, fontsize=7.5,
                    handlelength=1.8, labelspacing=.6, borderpad=.2, borderaxespad=.2)
    leg.set_zorder(10)
    fig.canvas.draw()
    lb = leg.get_window_extent(fig.canvas.get_renderer()).transformed(ax.transData.inverted())

    img = plt.imread(FIG / 'hake_illustration_Salinas_Nass.png')
    w = lb.width; h = w * img.shape[0] / img.shape[1]
    ix0, iy1 = lb.x0, lb.y0 - .01 * (ymax - ymin)
    ax.imshow(img, extent=(ix0, ix0 + w, iy1 - h, iy1), zorder=6, interpolation='lanczos')
    ax.text(ix0 + w / 2, iy1 - h - .01 * (ymax - ymin), 'South Pacific hake', ha='center', va='top', fontsize=7.5, zorder=7)
    ax.text(ix0 + w / 2, iy1 - h - .034 * (ymax - ymin), 'Merluccius gayi gayi', ha='center', va='top', fontsize=7.5,
            fontstyle='italic', zorder=7)

    # South America inset, inside the frame (lower left, over the open ocean)
    ring = geodetic_rect(-82.5, -34.0, -56.5, 13.0)                      # South America only, little Atlantic
    ix_min, ix_max, iy_min, iy_max = ring[:, 0].min(), ring[:, 0].max(), ring[:, 1].min(), ring[:, 1].max()
    pos = ax.get_position()
    in_w = .25 * pos.width
    in_h = in_w * (iy_max - iy_min) / (ix_max - ix_min) * (WIDTH_CM / HEIGHT_CM)
    ins = fig.add_axes([pos.x0, pos.y0, in_w, in_h], facecolor=SEA)      # shares the frame's lower-left corner
    ins.set_xlim(ix_min, ix_max); ins.set_ylim(iy_min, iy_max)
    ins.set_aspect('equal', adjustable='box', anchor='SW')               # box shrinks to the data, anchored to the corner
    ins.set_xticks([]); ins.set_yticks([])
    for sp in ins.spines.values():
        sp.set_edgecolor(INK); sp.set_linewidth(.6)
    add_geojson(ins, DER / 'LME_Humboldt.geojson', facecolor=LME_COL, edgecolor='none', lw=0, zorder=0.5)
    add_geojson(ins, EXT / 'natural_earth_paises_10m.geojson', facecolor=LAND, edgecolor=LAND_EDGE, lw=.2, zorder=1)
    ins.add_patch(MplPolygon(geodetic_rect(LON0, LON1, LAT0, LAT1), closed=True, facecolor='none', edgecolor=INK, lw=.8, zorder=5))

    # north arrow from the projection (true north at the arrow's own longitude), upper-right corner
    nx_, ny_ = xmin + .95 * (xmax - xmin), ymin + .905 * (ymax - ymin)
    nlo, nla = Transformer.from_crs(CRS, 'EPSG:4326', always_xy=True).transform(nx_, ny_)
    exx, eyy = FWD.transform(nlo, nla + 1.6)
    ax.annotate('', xy=(exx, eyy), xytext=(nx_, ny_), arrowprops={'arrowstyle': '-|>', 'color': INK, 'lw': 1.1}, zorder=9)
    ax.text(exx, eyy + .008 * (ymax - ymin), 'N', fontsize=9, ha='center', va='bottom', weight='bold', zorder=9)

    # scale bar, two segments of 250 km, bottom of the frame to the right of the inset
    sx, sy = xmin + .31 * (xmax - xmin), ymin + .035 * (ymax - ymin)
    alto = .006 * (ymax - ymin)
    for j in range(2):
        ax.add_patch(Rectangle((sx + j * 250000, sy), 250000, alto, facecolor=INK if j == 0 else 'white',
                               edgecolor=INK, lw=.6, zorder=9))
    for j, lab in enumerate(['0', '250', '500 km']):
        ax.text(sx + j * 250000, sy + 1.6 * alto, lab, ha='center', va='bottom', fontsize=6.5, zorder=9)

    fig.canvas.draw()
    for t in fig.texts + ax.texts:
        b = t.get_window_extent(fig.canvas.get_renderer())
        assert b.x0 >= 0 and b.y0 >= 0 and b.x1 <= fig.bbox.width and b.y1 <= fig.bbox.height, t.get_text()

    fig.savefig(out / 'Fig_1.png', dpi=600, facecolor=PAPER)
    print('wrote', out / 'Fig_1.png')


if __name__ == '__main__':
    main()
