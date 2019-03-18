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
import matplotlib.lines
# rcParams['figure.figsize'] = 16, 10
FIGSIZE = (12, 11)

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

def plot_corr(data_fh, target, panel, exome, xaxis, yaxis):
  logging.info('starting...')

  # Sample  Tags    Caller  DP      AF      Error   Variants        Multiplier      Signature.1     ...    Signature.30
  total = 0
  tags = set()
  for row in csv.DictReader(data_fh, delimiter='\t'):
    ok = True
    tags.add(row['Tags'])
    for f in panel:
      col, val = f.split('=')
      if val != row[col]:
        ok = False

    if ok:
      logging.debug(row)
      sigs_panel = [float(row['Signature.{}'.format(n)]) for n in range(1, 31)] # array of sig values
      total += 1

    ok = True
    for f in exome:
      col, val = f.split('=')
      if val != row[col]:
        ok = False

    if ok:
      logging.debug(row)
      sigs_exome = [float(row['Signature.{}'.format(n)]) for n in range(1, 31)] # array of sig values
      total += 1

  logging.info('found %i matching records', total)
  logging.debug('tags: %s', ' '.join(sorted(list(tags))))

  # generate plot
  if total != 2:
    logging.warn('expected two items')

  fig = plt.figure(figsize=FIGSIZE)
  ax = fig.add_subplot(111)
  #ax.set_xscale("log")#, nonposx='clip')
  #ax.set_yscale("log")#, nonposx='clip')
  
  maxs = 0.0
  for idx, (x, y) in enumerate(zip(sigs_panel, sigs_exome)):
    label = '{}'.format(idx + 1)
    ax.scatter([x], [y], c=[colors[idx]], label=label)
    ax.annotate(
        label,
        xy=(x, y), xytext=(1, 3),
        textcoords='offset points', ha='right', va='bottom')
    maxs = max(maxs, x, y)

  if yaxis is None:
    ax.set_ylabel('Exome')
  else:
    ax.set_ylabel(yaxis)
  if xaxis is None:
    ax.set_xlabel('Panel')
  else:
    ax.set_xlabel(xaxis)

  ax.set_title('Correlation for {}'.format(' '.join(exome)))
  ax.legend(loc="upper right", bbox_to_anchor=(0.99,0.90), bbox_transform=plt.gcf().transFigure)
  ax.set_xlim([0,maxs + 0.05])
  ax.set_ylim([0,maxs + 0.05])
  ax.plot([0, maxs], [0, maxs], color='#808080')

  logging.info('done processing %i', total)
  plt.savefig(target)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Plot changes in signature')
  parser.add_argument('--data', required=True, help='data file')
  parser.add_argument('--panel', nargs='*', help='other filters')
  parser.add_argument('--exome', nargs='*', help='other filters')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--target', required=False, default='correlation.png', help='plot filename')
  parser.add_argument('--xaxis', required=False, help='plot filename')
  parser.add_argument('--yaxis', required=False, help='plot filename')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  plot_corr(open(args.data, 'r'), args.target, args.panel, args.exome, args.xaxis, args.yaxis)
