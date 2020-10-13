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
LABELS = {'AF': 'Minimum Variant Allele Fraction', 'DP': 'Minimum Tumour Depth at Variant', 'Variants': 'Variant Count'}
COUNT_SPLIT = 0.05
X_HIGHLIGHT_COLOUR='#ff8000'

colors = {"SBS1": "#c0c0c0", "SBS2": "#41ac2f", "SBS3": "#7951d0", "SBS4": "#73d053", "SBS5": "#b969e9", "SBS6": "#91ba2c", "SBS7a": "#b4b42f", "SBS7b": "#5276ec", "SBS7c": "#daae36", "SBS7d": "#9e40b5", "SBS8": "#43c673", "SBS9": "#dd4cb0", "SBS10a": "#3d9332", "SBS10b": "#de77dd", "SBS11": "#7bad47", "SBS12": "#9479e8", "SBS13": "#487b21", "SBS14": "#a83292", "SBS15": "#83c67d", "SBS16": "#664db1", "SBS17a": "#e18d28", "SBS17b": "#588de5", "SBS18": "#e2672a", "SBS19": "#34c7dd", "SBS20": "#cf402b", "SBS21": "#5acdaf", "SBS22": "#d74587", "SBS23": "#428647", "SBS24": "#7b51a7", "SBS25": "#b4ba64", "SBS26": "#646cc1", "SBS27": "#a27f1f", "SBS28": "#3b63ac", "SBS29": "#dca653", "SBS30": "#505099", "SBS31": "#7d8529", "SBS32": "#bf8ade", "SBS33": "#516615", "SBS34": "#b65da7", "SBS35": "#57a87a", "SBS36": "#c84249", "SBS37": "#37b5b1", "SBS38": "#a14622", "SBS39": "#58b5e1", "SBS40": "#ba6e2f", "SBS41": "#589ed8", "SBS42": "#e98261", "SBS43": "#3176ae", "SBS44": "#656413", "SBS45": "#a19fe2", "SBS46": "#756121", "SBS47": "#7e4a8d", "SBS48": "#326a38", "SBS49": "#dd8abf", "SBS50": "#1a6447", "SBS51": "#e78492", "SBS52": "#30876c", "SBS53": "#9d4d7c", "SBS54": "#919d5b", "SBS55": "#9d70ac", "SBS56": "#5b6f34", "SBS57": "#65659c", "SBS58": "#c9a865", "SBS59": "#a1455d", "SBS60": "#5e622c", "SBS84": "#b66057", "SBS85": "#dca173", "DBS1": "#855524", "DBS2": "#9f7846", "DBS3": "#7951d0", "DBS4": "#73d053", "DBS5": "#b969e9", "DBS6": "#91ba2c", "DBS7": "#3656ca", "DBS8": "#b4b42f", "DBS9": "#5276ec", "DBS10": "#daae36", "DBS11": "#9e40b5", "ID1": "#de3860", "ID2": "#41ac2f", "ID3": "#7951d0", "ID4": "#73d053", "ID5": "#b969e9", "ID6": "#91ba2c", "ID7": "#9e40b5", "ID8": "#43c673", "ID9": "#dd4cb0", "ID10": "#3d9332", "ID11": "#de77dd", "ID12": "#7bad47", "ID13": "#9479e8", "ID14": "#487b21", "ID15": "#a83292", "ID16": "#83c67d", "ID17": "#664db1", "1": "#b66057", "2": "#dca173", "3": "#855524", "4": "#9f7846", "5": "#7951d0", "6": "#73d053", "7": "#b969e9", "8": "#91ba2c", "9": "#3656ca", "10": "#b4b42f", "11": "#5276ec", "12": "#daae36", "13": "#9e40b5", "14": "#de3860", "15": "#41ac2f", "16": "#7951d0", "17": "#73d053", "18": "#b969e9", "19": "#91ba2c", "20": "#9e40b5", "21": "#43c673", "22": "#dd4cb0", "23": "#3d9332", "24": "#de77dd", "25": "#7bad47", "26": "#9479e8", "27": "#487b21", "28": "#a83292", "29": "#83c67d", "30": "#664db1"}

#plot_hist(open(args.data, 'r'), args.samples, args.x, args.target, args.filters, args.title, args.logx, args.highlight, args.error_plot, args.count_plot, args.split_count, args.signature_ids, args.min_signature_val, args.x_highlight, args.height, args.width, args.fontsize)
def plot_hist(data_fh, samples, x, target, filters, title, logx, highlight, error_plot, count_plot, split_count, signature_ids, min_signature_val, x_highlight=None, height=8, width=12, fontsize=12):
  logging.info('starting...')
  # split count currently only for both error and count plots

  import matplotlib.style
  matplotlib.style.use('seaborn')
  rcParams.update({'font.size': fontsize})

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

  fig = plt.figure(figsize=(width, height))
  if error_plot and not count_plot: # just error plot
    grid = plt.GridSpec(6, 1, hspace=0, wspace=0)
    ax = fig.add_subplot(grid[0:-1, :])
    ax_err = fig.add_subplot(grid[-1, :], sharex=ax)
  elif error_plot and count_plot: # both
    if split_count:
      grid = plt.GridSpec(100, 1, hspace=0, wspace=0) # total height
      ax = fig.add_subplot(grid[0:-30, :]) # give to main plot
      ax_err = fig.add_subplot(grid[-29:-21, :], sharex=ax) # give to error plot
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

    label=signature_ids[sig]
    color = colors[signature_ids[sig]]
    if (highlight is None or label not in highlight) and max(ys) < min_signature_val:
      logging.debug('excluding signature %s due to no value greater than %.2f', sig, min_signature_val)
      ax.fill_between(xs, base, ys_cumulative, color=color)
    else:
      if highlight is None:
        ax.fill_between(xs, base, ys_cumulative, color=colors[signature_ids[sig]], label=label)
      else:
        if label in highlight:
          logging.info('highlighting %s', label)
          ax.fill_between(xs, base, ys_cumulative, color=colors[signature_ids[sig]], label=label)
        else:
          ax.fill_between(xs, base, ys_cumulative, color=colors[signature_ids[sig]], label=label, alpha=0.2)
    
    # new base
    base = ys_cumulative

  if x_highlight is not None:
    ax.axvline(x_highlight, color=X_HIGHLIGHT_COLOUR, ymin=0, ymax=1)

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
  parser.add_argument('--x_highlight', required=False, type=float, help='vertical line at this value')
  parser.add_argument('--logx', action='store_true', help='log x')
  parser.add_argument('--error_plot', action='store_true', help='include erorr plot')
  parser.add_argument('--count_plot', action='store_true', help='include count plot')
  parser.add_argument('--split_count', action='store_true', help='split count plot')
  parser.add_argument('--signature_ids', nargs='+', help='signature ids for plot')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--title', required=False, help='sample filter')
  parser.add_argument('--target', required=False, default='plot.png', help='plot filename')
  parser.add_argument('--min_signature_val', required=False, default=-1, type=float, help='minimum value of signature to include')
  parser.add_argument('--height', required=False, type=float, default=8, help='height of plot')
  parser.add_argument('--width', required=False, type=float, width=12, help='width of plot')
  parser.add_argument('--fontsize', required=False, default=12, type=int, help='plot font size')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  plot_hist(open(args.data, 'r'), args.samples, args.x, args.target, args.filters, args.title, args.logx, args.highlight, args.error_plot, args.count_plot, args.split_count, args.signature_ids, args.min_signature_val, args.x_highlight, args.height, args.width, args.fontsize)
