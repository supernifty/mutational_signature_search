#!/usr/bin/env python
'''
  plot evolution of signatures across a parameter for a single sample
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
#FIGSIZE = (16, 8)
FIGSIZE = (36, 18)
LABELS = {'AF': 'Minimum Variant Allele Fraction', 'DP': 'Minimum Depth', 'Variants': 'Variant Count'}
COUNT_SPLIT = 0.05

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
  '#cd4b33', # 10
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

def plot_hist(data_fh, samples, x, target, filters, title, logx, highlight, error_plot, count_plot, split_count, signature_ids, min_signature_val):
  logging.info('starting...')
  # split count currently only for both error and count plots

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
        sigs = [float(row[signature_id]) for signature_id in signature_ids] # array of sig values
        results.append((xval, sigs, float(row['Error']), float(row['Variants'])))
        logging.debug(row)

    total += 1

  logging.info('finished reading %i of %i records', included, total)
  logging.debug('tags: %s', ' '.join(sorted(list(tags))))

  # generate plot
  # sort xvals
  results = sorted(results, key=lambda r: r[0])

  if len(results) == 0:
    logging.warn('No data to plot')
    return

  fig = plt.figure(figsize=FIGSIZE)
  if error_plot and not count_plot: # just error plot
    grid = plt.GridSpec(6, 1, hspace=0, wspace=0)
    ax = fig.add_subplot(grid[0:-1, :])
    ax_err = fig.add_subplot(grid[-1, :], sharex=ax)
  elif error_plot and count_plot: # both
    if split_count:
      grid = plt.GridSpec(100, 1, hspace=0, wspace=0) # total height
      ax = fig.add_subplot(grid[0:-40, :]) # give to main plot
      ax_err = fig.add_subplot(grid[-39:-21, :], sharex=ax) # give to error plot
      ax_count_top = fig.add_subplot(grid[-20:-11, :], sharex=ax) # give to half of variant count plot
      ax_count_bottom = fig.add_subplot(grid[-10:, :], sharex=ax) # give to half of variant count plot
    else:
      grid = plt.GridSpec(7, 1, hspace=0, wspace=0)
      ax = fig.add_subplot(grid[0:-2, :])
      ax_err = fig.add_subplot(grid[-2, :], sharex=ax)
      ax_count = fig.add_subplot(grid[-1, :], sharex=ax)
  elif not error_plot and count_plot: # just count plot
    grid = plt.GridSpec(6, 1, hspace=0, wspace=0)
    ax = fig.add_subplot(grid[0:-1, :])
    ax_count = fig.add_subplot(grid[-1, :], sharex=ax)
  else:
    ax = fig.add_subplot(111)

  if logx:
    ax.set_xscale("log", nonposx='clip')
  #fig, ax = plt.subplots()

  xs = [r[0] for r in results]
  base = [0] * len(xs)
  for sig in range(0, len(signature_ids)):
    ys = [r[1][sig] for r in results]
    ys_cumulative = [ys[i] + base[i] for i in range(len(ys))]

    if max(ys) < min_signature_val:
      logging.debug('excluding signature %s due to no value greater than %.2f', min_signature_val)
      ax.fill_between(xs, base, ys_cumulative, color=colors[hash(signature_ids[sig]) % len(colors)])
    else:
      label=signature_ids[sig]
      if highlight is None:
        ax.fill_between(xs, base, ys_cumulative, color=colors[hash(signature_ids[sig]) % len(colors)], label=label)
      else:
        if label in highlight:
          ax.fill_between(xs, base, ys_cumulative, color=colors[hash(signature_ids[sig]) % len(colors)], label=label)
        else:
          ax.fill_between(xs, base, ys_cumulative, color=colors[hash(signature_ids[sig]) % len(colors)], label=label, alpha=0.2)
    
    # new base
    base = ys_cumulative

  ax.set_ylabel('Signature proportion')
  if not error_plot and not count_plot:
    ax.set_xlabel(LABELS[x])
  else:
    ax.tick_params(labeltop=False, labelbottom=False)

  if title is None:
    ax.set_title('Somatic mutational signatures detected by {} for {} ({})'.format(x, ' '.join(samples), ' '.join(filters)))
  else:
    ax.set_title(title)
  #ax.legend(loc="upper right", bbox_to_anchor=(0.99,0.90), bbox_transform=plt.gcf().transFigure)

  if error_plot:
    ys = [r[2] for r in results]
    ax_err.plot(xs, ys, 'k-', linewidth=0.5)
    ax_err.fill_between(xs, 0, ys, color='#a0a0a0')
    ax_err.grid(True)
    ax_err.set_ylabel('Error')
    ax_err.set_ylim((0, 1.0))
    if count_plot:
      ax_err.tick_params(labeltop=False, labelbottom=False)
    else:
      ax_err.set_xlabel(LABELS[x])

  if count_plot:
    ys = [r[3] for r in results]
    if split_count:
      ymax = max(ys)
      ymid = ys[int(len(ys) / 2)] # roughly the median?
      # based on https://matplotlib.org/examples/pylab_examples/broken_axis.html
      ax_count_top.plot(xs, ys, 'k-', linewidth=0.5)
      ax_count_bottom.plot(xs, ys, 'k-', linewidth=0.5)
      ax_count_top.fill_between(xs, 0, ys, color='#a0a0a0')
      ax_count_bottom.fill_between(xs, 0, ys, color='#a0a0a0')

      #ax_count_top.set_ylim(ymax * COUNT_SPLIT, ymax)  # outliers only
      #ax_count_bottom.set_ylim(0, ymax * COUNT_SPLIT)  # most of the data

      ax_count_top.set_ylim(ymid, ymax)  # top "half" of the data
      ax_count_bottom.set_ylim(0, ymid)  # bottom "half" of the data

      # hide the spines between ax and ax2
      ax_count_top.spines['bottom'].set_visible(False)
      ax_count_bottom.spines['top'].set_visible(False)
      ax_count_top.xaxis.tick_top()
      ax_count_top.tick_params(labeltop=False)  # don't put tick labels at the top
      ax_count_bottom.xaxis.tick_bottom()

      ax_count_top.grid(True)
      ax_count_bottom.grid(True)

      dx = .01  # how big to make the diagonal lines in axes coordinates
      dy = .01  # how big to make the diagonal lines in axes coordinates
      # arguments to pass to plot, just so we don't keep repeating them
      kwargs = dict(transform=ax_count_top.transAxes, color='k', clip_on=False)
      ax_count_top.plot((-dx, +dy), (-dx, +dy), **kwargs)        # top-left diagonal
      ax_count_top.plot((1 - dx, 1 + dy), (-dx, +dy), **kwargs)  # top-right diagonal

      kwargs.update(transform=ax_count_bottom.transAxes)  # switch to the bottom axes
      ax_count_bottom.plot((-dx, +dy), (1 - dx, 1 + dy), **kwargs)  # bottom-left diagonal
      ax_count_bottom.plot((1 - dx, 1 + dy), (1 - dx, 1 + dy), **kwargs)  # bottom-right diagonal
      ax_count_bottom.set_xlabel(LABELS[x])
      ax_count_bottom.set_ylabel('Variants')
    else:
      ax_count.plot(xs, ys, 'k-', linewidth=0.5)
      ax_count.fill_between(xs, 0, ys, color='#a0a0a0')
      ax_count.grid(True)
      ax_count.set_ylabel('Variants')
      ax_count.set_xlabel(LABELS[x])

  # place legend at right based on https://stackoverflow.com/questions/10101700/moving-matplotlib-legend-outside-of-the-axis-makes-it-cutoff-by-the-figure-box/10154763#10154763
  handles, labels = ax.get_legend_handles_labels()
  lgd = ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(1.01,1.0), borderaxespad=0)
  lgd.get_frame().set_edgecolor('#000000')
  fig.savefig(target, bbox_extra_artists=(lgd,), bbox_inches='tight')

  logging.info('done processing %i of %i', included, total)
  matplotlib.pyplot.close('all')


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Plot changes in signature')
  parser.add_argument('--data', required=True, help='data file')
  parser.add_argument('--samples', nargs='+', required=True, help='sample filter')
  parser.add_argument('--filters', nargs='*', help='other filters')
  parser.add_argument('--highlight', nargs='*', help='signature(s) to highlight')
  parser.add_argument('--x', required=True, help='x column name')
  parser.add_argument('--logx', action='store_true', help='log x')
  parser.add_argument('--error_plot', action='store_true', help='include erorr plot')
  parser.add_argument('--count_plot', action='store_true', help='include count plot')
  parser.add_argument('--split_count', action='store_true', help='split count plot')
  parser.add_argument('--signature_ids', nargs='+', help='signature ids for plot')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--title', required=False, help='sample filter')
  parser.add_argument('--target', required=False, default='plot.png', help='plot filename')
  parser.add_argument('--min_signature_val', required=False, default=-1, type=float, help='minimum value of signature to include')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  plot_hist(open(args.data, 'r'), args.samples, args.x, args.target, args.filters, args.title, args.logx, args.highlight, args.error_plot, args.count_plot, args.split_count, args.signature_ids, args.min_signature_val)
