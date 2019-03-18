#!/usr/bin/env python
'''
  given tumour and normal vcf pairs, explore msi status
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
FIGSIZE = (16, 8)

# from http://tools.medialab.sciences-po.fr/iwanthue/
colors = [ 
  '#cc5ac5',
  '#61c350',
  '#9358ca',
  '#9bb932',
  '#626cd4',
  '#d8ab39',
  '#5e73b8',
  '#dc842f',
  '#60a4da',
  '#cd4b33',
  '#4ec584',
  '#db3666',
  '#3d9037',
  '#d5509a',
  '#69862a',
  '#cb91d4',
  '#b4ad4a',
  '#98518b',
  '#9db76f',
  '#a44657',
  '#42c3c2',
  '#e07a8a',
  '#2b7650',
  '#e38f6b',
  '#5aad86',
  '#995a2b',
  '#5d8547',
  '#c39e64',
  '#676828',
  '#977b1f'
]

def plot_hist(data_fh, samples, x, target, filters, title, logx, highlight, error_plot):
  logging.info('starting...')

  # Sample  Tags    Caller  DP      AF      Error   Variants        Multiplier      Signature.1     ...    Signature.30
  included = total = 0
  results = []
  tags = set()
  for row in csv.DictReader(data_fh, delimiter='\t'):
    if row['Sample'] in samples:
      tags.add(row['Tags'])
      ok = True
      for f in filters:
        col, val = f.split('=')
        if val != row[col]:
          ok = False

      if ok:
        included += 1
        xval = float(row[x]) # x axis value
        sigs = [float(row['Signature.{}'.format(n)]) for n in range(1, 31)] # array of sig values
        results.append((xval, sigs, float(row['Error'])))
        logging.debug(row)

    total += 1

  logging.info('finished reading %i of %i records', included, total)
  logging.debug('tags: %s', ' '.join(sorted(list(tags))))

  # generate plot
  # sort xvals
  results = sorted(results, key=lambda r: r[0])

  if len(results) == 0:
    logging.warn('No data to plot')
    sys.exit(0)

  fig = plt.figure(figsize=FIGSIZE)
  if error_plot:
    grid = plt.GridSpec(6, 6, hspace=0, wspace=0)
    ax = fig.add_subplot(grid[0:-1, :])
    ax_err = fig.add_subplot(grid[-1, :], sharex=ax)
  else:
    ax = fig.add_subplot(111)

  if logx:
    ax.set_xscale("log", nonposx='clip')
  #fig, ax = plt.subplots()

  xs = [r[0] for r in results]
  base = [0] * len(xs)
  for sig in range(0, 30):
    ys = [r[1][sig] for r in results]
    ys_cumulative = [ys[i] + base[i] for i in range(len(ys))]

    if highlight is None:
      ax.fill_between(xs, base, ys_cumulative, color=colors[sig], label='{}'.format(sig + 1))
    else:
      label='{}'.format(sig + 1)
      if label in highlight:
        ax.fill_between(xs, base, ys_cumulative, color=colors[sig], label=label)
      else:
        ax.fill_between(xs, base, ys_cumulative, color=colors[sig], label=label, alpha=0.2)
    
    # new base
    base = ys_cumulative

  ax.set_ylabel('Signature proportion')
  ax.set_xlabel(x)
  if title is None:
    ax.set_title('Somatic mutational signatures detected by {} {} ({})'.format(x, ' '.join(samples), ' '.join(filters)))
  else:
    ax.set_title('Somatic mutational signatures detected by {} {} ({})'.format(x, ' '.join(samples), title))
  #ax.legend(loc='best')
  #ax.legend(loc="upper right", bbox_to_anchor=(0.99,0.90), bbox_transform=plt.gcf().transFigure)
  ax.legend(loc="upper right", bbox_to_anchor=(0.99,0.90), bbox_transform=plt.gcf().transFigure)

  if error_plot:
    ys = [r[2] for r in results]
    ax_err.plot(xs, ys, 'k-', linewidth=0.5)
    ax_err.grid(True)
    ax_err.set_ylabel('Error')
    ax_err.set_ylim((0, 1.0))

  logging.info('done processing %i of %i', included, total)
  plt.savefig(target)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Plot changes in signature')
  parser.add_argument('--data', required=True, help='data file')
  parser.add_argument('--samples', nargs='+', required=True, help='sample filter')
  parser.add_argument('--filters', nargs='*', help='other filters')
  parser.add_argument('--highlight', nargs='*', help='signature(s) to highlight')
  parser.add_argument('--x', required=True, help='x column name')
  parser.add_argument('--logx', action='store_true', help='log x')
  parser.add_argument('--error_plot', action='store_true', help='include erorr plot')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--title', required=False, help='sample filter')
  parser.add_argument('--target', required=False, default='plot.png', help='plot filename')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  plot_hist(open(args.data, 'r'), args.samples, args.x, args.target, args.filters, args.title, args.logx, args.highlight, args.error_plot)
