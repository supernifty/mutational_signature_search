#!/usr/bin/env python
'''
  plot signatures seen across two inputs
'''

import argparse
import csv
import logging
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import rcParams
# rcParams['figure.figsize'] = 16, 10

LABELS = {'DP': 'Minimum Tumour Depth at Variant', 'AF': 'Minimum Variant Allele Fraction'}
DPI = 300

def plot_heat(data_fh, samples, xlabel, ylabel, target, filters, title, highlight, y_multiples, max_error, width=10, height=10, fontsize=8):
  figsize = (width, height)
  logging.info('starting...')

  import matplotlib.style
  matplotlib.style.use('default')
  rcParams.update({'font.size': fontsize})

  # Sample  Tags    Caller  DP      AF      Error   Variants        Multiplier      Signature.1     ...    Signature.30
  included = total = 0
  results = {}
  tags = set()
  xvals = set()
  yvals = set()
  max_zval = 0.0
  for row in csv.DictReader(data_fh, delimiter='\t'):
    if row['Sample'] in samples:
      tags.add(row['Tags'])
      ok = True
      for f in filters:
        col, val = f.split('=')
        if val != row[col]:
          ok = False

      if y_multiples is not None and int(row['DP']) % y_multiples != 0:
        ok = False

      if ok:
        included += 1
        xval = float(row[xlabel]) # x axis value
        yval = float(row[ylabel]) # y axis value
        xvals.add(xval)
        yvals.add(yval)
        #zval = float(row['Signature.{}'.format(highlight)])
        zval = float(row[highlight])
        if row['Error'] == 'nan':
          results['{},{}'.format(xval, yval)] = (0.0, 1)
        elif float(row['Error']) > max_error:
          results['{},{}'.format(xval, yval)] = (0.0, float(row['Error']))
        else:
          results['{},{}'.format(xval, yval)] = (zval, float(row['Error']))
          max_zval = max(max_zval, zval)
        logging.debug(row)

    total += 1

  logging.info('finished reading %i of %i records', included, total)
  logging.debug('tags: %s', ' '.join(sorted(list(tags))))

  if len(results) == 0:
    logging.warn('No data to plot')
    return

  xvals = sorted(list(xvals))
  yvals = sorted(list(yvals))[::-1] # bottom left

  zvals = []
  tvals = []

  for y in yvals:
    zrow = []
    trow = []
    for x in xvals:
      key = '{},{}'.format(x, y)
      if key in results:
        zrow.append(results[key][0])
        if results[key][1] > max_error:
          trow.append('-\n{:.0f}%'.format(results[key][1] * 100))
        else:
          trow.append('{:.0f}%\n{:.0f}%'.format(results[key][0] * 100, results[key][1] * 100))
      else:
        zrow.append(0.0)
        trow.append('')
    zvals.append(zrow)
    tvals.append(trow)

  fig = plt.figure(figsize=(width, height))
  ax = fig.add_subplot(111)
  im = ax.imshow(zvals, cmap='viridis')

  cbar = ax.figure.colorbar(im, ax=ax, fraction=0.04, pad=0.01, shrink=0.9)
  cbar.ax.set_ylabel('Signature {} %'.format(highlight), rotation=-90, va="bottom")

  ax.set_xticks(range(len(xvals)))
  ax.set_yticks(range(len(yvals)))
  ax.set_xticklabels(xvals, rotation=90)
  ax.set_yticklabels(yvals)

  ax.set_ylabel(LABELS[ylabel])
  ax.set_xlabel(LABELS[xlabel])

  for y in range(len(yvals)):
    for x in range(len(xvals)):
      if zvals[y][x] > max_zval * 0.5:
        text = ax.text(x, y, tvals[y][x], ha="center", va="center", color="k")
      else:
        text = ax.text(x, y, tvals[y][x], ha="center", va="center", color="w")

  if title is None:
    ax.set_title('Signature {} given {} and {}'.format(highlight, xlabel, ylabel))
  else:
    ax.set_title(title)

  logging.info('done processing %i of %i', included, total)
  plt.tight_layout()
  plt.savefig(target, dpi=DPI)
  matplotlib.pyplot.close('all')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Plot changes in signature')
  parser.add_argument('--data', required=True, help='data file')
  parser.add_argument('--samples', nargs='+', required=True, help='sample filter')
  parser.add_argument('--filters', nargs='*', default=[], help='other filters')
  parser.add_argument('--highlight', required=True, help='signature(s) to highlight')
  parser.add_argument('--x', required=True, help='x column name')
  parser.add_argument('--y', required=True, help='y column name')
  parser.add_argument('--y_multiples', type=int, required=False, help='y column name')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--title', required=False, help='sample filter')
  parser.add_argument('--target', required=False, default='plot.png', help='plot filename')
  parser.add_argument('--max_error', type=float, default=0.5, help='x column name')
  parser.add_argument('--height', required=False, type=float, default=10, help='height of plot')
  parser.add_argument('--width', required=False, type=float, default=10, help='width of plot')
  parser.add_argument('--fontsize', required=False, type=float, width=8, help='width of plot')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  plot_heat(open(args.data, 'r'), args.samples, args.x, args.y, args.target, args.filters, args.title, args.highlight, args.y_multiples, args.max_error, args.width, args.height, args.fontsize)
