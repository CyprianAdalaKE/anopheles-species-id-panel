#!/usr/bin/env python3
"""Fig. 1: study-site map. Kenya county boundaries (geoBoundaries gbOpen KEN ADM1, RCMRD, public domain);
collection villages = mean GPS position per village from SPP_METADATA.xlsx (collection_sites_clean.tsv).
Village points are used for the site map only; they are not linked to individual sequenced specimens."""
import json, csv, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, ConnectionPatch
from matplotlib.path import Path
HERE = os.path.dirname(os.path.abspath(__file__))
G = json.load(open(f'{HERE}/ken_adm1.geojson'))
SITES = list(csv.DictReader(open(f'{HERE}/collection_sites_clean.tsv'), delimiter='\t'))
STUDY = {'Siaya', 'Turkana'}
FILL_S, EDGE_S, FILL_O, EDGE_O, PT = '#F2997A', '#9A3A1E', '#EAF2DF', '#3E6B1F', '#1F3B73'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9})

def rings(f):
    g = f['geometry']; polys = g['coordinates'] if g['type'] == 'MultiPolygon' else [g['coordinates']]
    return [p[0] for p in polys]

def draw(ax, lw_s=1.0, lw_o=0.5):
    for f in G['features']:
        s = f['properties']['shapeName'] in STUDY
        for r in rings(f):
            ax.add_patch(Polygon(r, closed=True, fc=FILL_S if s else FILL_O, ec=EDGE_S if s else EDGE_O, lw=lw_s if s else lw_o, zorder=2 if s else 1))

# sanity: each village lies inside its county polygon
for s in SITES:
    f = next(f for f in G['features'] if f['properties']['shapeName'] == s['county'])
    inside = any(Path(r).contains_point((float(s['longitude']), float(s['latitude']))) for r in rings(f))
    assert inside, s

fig = plt.figure(figsize=(7.2, 5.9))
ax = fig.add_axes([0.075, 0.08, 0.58, 0.9]); draw(ax)
for s in SITES: ax.plot(float(s['longitude']), float(s['latitude']), 'o', ms=3.5, mfc=PT, mec='white', mew=0.4, zorder=5)
for name, (x, y, dx, dy) in {'Turkana': (35.4, 3.6, 0.0, 0), 'Siaya': (34.3, 0.5, 0, 0)}.items():
    ax.text(x, y, name, fontsize=9, fontweight='bold', color='#3a2a22', ha='center', zorder=6)
ax.plot(36.817, -1.286, marker='*', ms=11, color='black', zorder=6); ax.text(36.95, -1.55, 'Nairobi', fontsize=8, style='italic', zorder=6)
ax.set_xlim(33.6, 42.2); ax.set_ylim(-5.0, 5.4); ax.set_aspect('equal')
ax.set_xlabel('Longitude (°E)'); ax.set_ylabel('Latitude (°N)'); ax.grid(ls=':', lw=0.4, color='#bbbbbb', zorder=0)
ax.annotate('', xy=(41.6, 4.85), xytext=(41.6, 4.2), arrowprops=dict(arrowstyle='-|>', lw=1.2, color='black')); ax.text(41.6, 5.12, 'N', ha='center', va='center', fontweight='bold', fontsize=9)
ax.plot([33.9, 33.9 + 100 / 111.3], [-4.6, -4.6], color='black', lw=2.2); ax.text(33.9 + 50 / 111.3, -4.45, '100 km', ha='center', fontsize=8)
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
ax.legend(handles=[Patch(fc=FILL_S, ec=EDGE_S, label='Study counties'), Patch(fc=FILL_O, ec=EDGE_O, label='Other counties'),
                   Line2D([], [], marker='o', ls='', mfc=PT, mec='white', label='Collection villages')], loc='lower right', fontsize=8, framealpha=0.95)

OFF = {'Arumrum': (-5, 2, 'right'), 'Nakwalele A': (5, -2, 'left'), 'Lopur': (5, 3, 'left'), 'Naperobei': (-5, -2, 'right'),
       'Tomokia': (5, 2, 'left'), "Ng'adirkonyen": (5, -1, 'left'), 'Bukhoba ABE': (-1, 9, 'center'), 'Bukhoba CDF': (5, -6, 'left'),
       'Gendro': (-5, 0, 'right'), 'Kanyango': (5, -4, 'left')}
def inset(rect, county, xlim, ylim, label, scale_km):
    a = fig.add_axes(rect); draw(a, 1.2, 0.6)
    for s in [s for s in SITES if s['county'] == county]:
        x, y = float(s['longitude']), float(s['latitude']); a.plot(x, y, 'o', ms=5, mfc=PT, mec='white', mew=0.6, zorder=5)
        dx, dy, ha = OFF.get(s['village'], (5, 2, 'left'))
        a.annotate(s['village'], (x, y), xytext=(dx, dy), textcoords='offset points', fontsize=6.8, zorder=6, ha=ha, va='center',
                   bbox=dict(fc='white', ec='none', alpha=0.7, pad=0.3))
    a.set_xlim(*xlim); a.set_ylim(*ylim); a.set_aspect('equal'); a.set_xticks([]); a.set_yticks([])
    a.set_title(label, fontsize=8.5, fontweight='bold', pad=3)
    x0 = xlim[0] + 0.06 * (xlim[1] - xlim[0]); y0 = ylim[0] + 0.07 * (ylim[1] - ylim[0])
    a.plot([x0, x0 + scale_km / 111.3], [y0, y0], color='black', lw=1.8); a.text(x0 + scale_km / 222.6, y0 + 0.04 * (ylim[1] - ylim[0]), f'{scale_km} km', ha='center', fontsize=7)
    for sp in a.spines.values(): sp.set_edgecolor('#555555')
    ax.add_patch(Rectangle((xlim[0], ylim[0]), xlim[1] - xlim[0], ylim[1] - ylim[0], fill=False, ec='#333333', lw=0.8, zorder=7))
    return a
a1 = inset([0.68, 0.55, 0.30, 0.40], 'Turkana', (35.33, 35.53), (2.20, 2.42), 'Turkana County', 5)
a2 = inset([0.68, 0.06, 0.30, 0.40], 'Siaya', (34.03, 34.21), (0.00, 0.15), 'Siaya County (Alego-Usonga)', 5)
fig.savefig(f'{HERE}/figure1_study_sites.png', dpi=300, facecolor='white'); fig.savefig(f'{HERE}/figure1_study_sites.pdf', facecolor='white')
print('saved')
