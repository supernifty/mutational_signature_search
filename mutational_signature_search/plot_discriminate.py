#!/usr/bin/env python
'''
  plot separation of samples by signature
'''

import argparse
import collections
import csv
import logging
import math
import sys

import numpy as np
import scipy
import scipy.optimize
import scipy.stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import rcParams
DPI=300

LINE_WIDTH=2
ALPHA=0.7

HIGHLIGHT_LINE_WIDTH=3
HIGHLIGHT_MARKER_SIZE=8
HIGHLIGHT_ALPHA=1.0

MEAN_LINEWIDTH=3
STD_LINEWIDTH=2

#PHENOTYPE_COLOR = '#61c350'
#CONTROL_COLOR =  '#98518b'
#EXTRA_COLOR = '#903000'

#GROUP_COLOURS = [
#  '#6cbe45', 
#  '#0039a6', 
#  '#ff2222', 
#  '#808183', 
#  '#462759', 
#  '#ee352e', 
#  '#b933ad', 
#  '#ff6319', 
#  '#903000'
#]

GROUP_COLOURS = [
  '#75ab3d',
  '#c554b2',
  '#7e6acc',
  '#b49242',
  '#6797d0',
  '#ca633b',
  '#54a676',
  '#ce4560',
  '#be6b91'
]

MARKERS = ('.', 'v', 'o', 'p', '^', '<', '>', '1', '+', '2', '3', 'x', '4', '8', ',', 's', 'P', '*', 'h', 'H', 'X', 'D', 'd')

LABELS = {'DP': 'Minimum Tumour Depth at Variant', 'AF': 'Minimum Variant Allele Fraction'}

# blue red green
#CONFIDENCE_COLORS=['#0000cc', '#cc0000', '#00cc00']
#CONFIDENCE_COLORS=['#cc0000', '#f09000', '#0000cc', '#00cc00']

# red blue green
CONFIDENCE_COLORS=['#cc0000', '#0000cc', '#00cc00']
CONFIDENCE_ALPHA=0.8
CONFIDENCE_LINEWIDTH=4

CONFIDENCE_MAXERROR=0.20
CONFIDENCE_PRIOR_PHENOTYPE=0.5
CONFIDENCE_DELTA=0.01
CONFIDENCE_RESOLUTION=0.005
CONFIDENCE_MIN_DENOM=1e-100

def find_confidence(goal, pu, cu, psd, csd, distribution, max_yval, report=False, point_estimate=True):
  def calculate(x):
    if distribution == 'normal':
      logging.info('normal parameters for positive: u=%.2f sd=%.2f', pu, psd)
      logging.info('normal parameters for negative: u=%.2f sd=%.2f', cu, csd)

      if point_estimate:
        prob_p = abs(scipy.stats.norm.cdf(x - CONFIDENCE_DELTA, loc=pu, scale=psd) - scipy.stats.norm.cdf(x + CONFIDENCE_DELTA, loc=pu, scale=psd))
        prob_c = abs(scipy.stats.norm.cdf(x - CONFIDENCE_DELTA, loc=cu, scale=csd) - scipy.stats.norm.cdf(x + CONFIDENCE_DELTA, loc=cu, scale=csd))
      else:
        prob_p = 1 - scipy.stats.norm.cdf(x, loc=pu, scale=psd)
        prob_c = 1 - scipy.stats.norm.cdf(x, loc=cu, scale=csd)
    
      prob = CONFIDENCE_PRIOR_PHENOTYPE * prob_p / (CONFIDENCE_PRIOR_PHENOTYPE * prob_p + (1 - CONFIDENCE_PRIOR_PHENOTYPE) * prob_c)
      #logging.debug('prob at %f is %f for %f to %f', x, prob, cu, pu)
      return abs(prob - goal) # error
    elif distribution == 'beta':
      # scale by max_yval
      pu_scaled = pu / max_yval
      psd_scaled = psd / max_yval
      cu_scaled = cu / max_yval
      csd_scaled = csd / max_yval
      x_scaled = x / max_yval

      pn = pu_scaled * (1 - pu_scaled) / (psd_scaled ** 2)
      pa = pu_scaled * pn
      pb = (1 - pu_scaled) * pn

      cn = cu_scaled * (1 - cu_scaled) / (csd_scaled ** 2)
      ca = cu_scaled * cn
      cb = (1 - cu_scaled) * cn

      logging.info('beta parameters for positive: a=%.2f b=%.2f', pa, pb)
      logging.info('beta parameters for negative: a=%.2f b=%.2f', ca, cb)

      if point_estimate:
        prob_p = abs(scipy.stats.beta.cdf(x_scaled - CONFIDENCE_DELTA, pa, pb) - scipy.stats.beta.cdf(x_scaled + CONFIDENCE_DELTA, pa, pb))
        prob_c = abs(scipy.stats.beta.cdf(x_scaled - CONFIDENCE_DELTA, ca, cb) - scipy.stats.beta.cdf(x_scaled + CONFIDENCE_DELTA, ca, cb))
      else:
        prob_p = 1 - scipy.stats.beta.cdf(x_scaled, pa, pb)
        prob_c = 1 - scipy.stats.beta.cdf(x_scaled, ca, cb)
    
      denominator = (CONFIDENCE_PRIOR_PHENOTYPE * prob_p + (1 - CONFIDENCE_PRIOR_PHENOTYPE) * prob_c)
      if denominator < CONFIDENCE_MIN_DENOM:
        logging.info('goal %f is unstable with denominator %f prob_p %f prob_c %f pu %f cu %f psd %f csd %f', goal, denominator, prob_p, prob_c, pu, cu, psd, csd)
        # let's continue
        #return 1.0 # max error
        #raise RuntimeError('unstable denominator')

      prob = CONFIDENCE_PRIOR_PHENOTYPE * prob_p / denominator
      if report:
        logging.info('prob at %f is %f goal is %f scaled params x %f pus %f psds %f cus %f csds %f unscaled %f %f %f %f', x, prob, goal, x_scaled, pu_scaled, psd_scaled, cu_scaled, csd_scaled, pu, psd, cu, csd)
      return abs(prob - goal) # error

  return calculate

def solve(conf, pu, cu, psd, csd, distribution, max_yval, point_estimate=True, first=True):
  '''
    return best value for given conf
  '''
  pu_scaled = pu / max_yval
  psd_scaled = psd / max_yval
  cu_scaled = cu / max_yval
  csd_scaled = csd / max_yval

  pn = pu_scaled * (1 - pu_scaled) / (psd_scaled ** 2)
  pa = pu_scaled * pn
  pb = (1 - pu_scaled) * pn

  cn = cu_scaled * (1 - cu_scaled) / (csd_scaled ** 2)
  ca = cu_scaled * cn
  cb = (1 - cu_scaled) * cn

  best = (None, 1.0)
  last = 0.0

  for x in np.arange(-0.01, max_yval + 0.01, CONFIDENCE_RESOLUTION):
    x_scaled = x / max_yval
    if point_estimate:
      prob_p = abs(scipy.stats.beta.cdf(x_scaled - CONFIDENCE_DELTA, pa, pb) - scipy.stats.beta.cdf(x_scaled + CONFIDENCE_DELTA, pa, pb))
      prob_c = abs(scipy.stats.beta.cdf(x_scaled - CONFIDENCE_DELTA, ca, cb) - scipy.stats.beta.cdf(x_scaled + CONFIDENCE_DELTA, ca, cb))
    else:
      prob_p = 1 - scipy.stats.beta.cdf(x_scaled, pa, pb)
      prob_c = 1 - scipy.stats.beta.cdf(x_scaled, ca, cb)
    denominator = (CONFIDENCE_PRIOR_PHENOTYPE * prob_p + (1 - CONFIDENCE_PRIOR_PHENOTYPE) * prob_c)
    if denominator < CONFIDENCE_MIN_DENOM:
      # too unlikely
      continue
    # phenotype probability
    prob = CONFIDENCE_PRIOR_PHENOTYPE * prob_p / denominator
    logging.debug('x {} xscaled {} prob {}'.format(x, x_scaled, prob))
    if first:
      if prob >= conf and last < conf:
        best = (x, abs(prob - conf))
        return best
      elif prob < conf and last > conf:
        best = (x, abs(prob - conf))
        return best
    else:
      if abs(prob - conf) < best[1]:
        best = (x, abs(prob - conf))

    last = prob
    
  return best

def calculate_confidence(conf, phenotype, control, distribution, max_yval, point_estimate):
  logging.debug('phenotype: %s; control: %s', phenotype, control)
  # assume len(control) > 1
  if len(phenotype) > 1:
    # want to find the point where tail(phenotype) / (tail(phenotype) + tail(control)) = conf
    try:
      start = (np.mean(phenotype) + np.mean(control)) / 2
      pmean = np.mean(phenotype)
      cmean = np.mean(control)
      pstd = np.std(phenotype)
      cstd = np.std(control)
      #discriminant = scipy.optimize.newton(find_confidence(conf, np.mean(phenotype), np.mean(control), np.std(phenotype), np.std(control)), start)
      #discriminant = scipy.optimize.brute(find_confidence(conf, pmean, cmean, pstd, cstd, distribution, max_yval), (slice(min(pmean, cmean), max(pmean, cmean), 0.01),))
      #discriminant = scipy.optimize.brute(find_confidence(conf, pmean, cmean, pstd, cstd, distribution, max_yval), (slice(-CONFIDENCE_MAXERROR, max_yval + CONFIDENCE_MAXERROR, 0.01),))
      #discriminant = scipy.optimize.brute(find_confidence(conf, pmean, cmean, pstd, cstd, distribution, max_yval, point_estimate=point_estimate), (slice(-0.01, max_yval, CONFIDENCE_RESOLUTION),))
      #error = find_confidence(conf, pmean, cmean, pstd, cstd, distribution, max_yval, report=True, point_estimate=point_estimate)(discriminant)
      discriminant, error = solve(conf, pmean, cmean, pstd, cstd, distribution, max_yval, point_estimate=point_estimate, first=True)
      if discriminant is None:
        logging.warn('no discriminant for confidence %.3f: %s vs %s error %.3f', conf, phenotype, control, error)
        return None

      if error > CONFIDENCE_MAXERROR or discriminant < -1 or discriminant > max_yval + CONFIDENCE_MAXERROR:
        logging.warn('failed to find confidence within error for confidence %.3f: %s vs %s error %.3f discriminant %.3f', conf, phenotype, control, error, discriminant)
        # let's just return it
        #return discriminant
        return None
      logging.debug('discriminant %f for %s vs %s at %f', discriminant, phenotype, control, conf)
      return discriminant
    except ValueError:
      logging.warn('calculate_confidence: discriminant did not converge with %s vs %s', phenotype, control)
      return None
    except RuntimeError:
      logging.warn('calculate_confidence: something went wrong with %s vs %s', phenotype, control)
      return None
  else:
    # here we assume the phenotype has a distribution matching the control group and calculate confidence similar to above
    try:
      start = (np.mean(phenotype) + np.mean(control)) / 2
      pmean = np.mean(phenotype)
      cmean = np.mean(control)
      cstd = np.std(control)
      pstd = cstd # assume same sd
      #discriminant = scipy.optimize.newton(find_confidence(conf, np.mean(phenotype), np.mean(control), np.std(phenotype), np.std(control)), start)
      #discriminant = scipy.optimize.brute(find_confidence(conf, pmean, cmean, pstd, cstd, distribution, max_yval), (slice(min(pmean, cmean), max(pmean, cmean), 0.01),))
      #discriminant = scipy.optimize.brute(find_confidence(conf, pmean, cmean, pstd, cstd, distribution, max_yval), (slice(-CONFIDENCE_MAXERROR, max_yval + CONFIDENCE_MAXERROR, 0.01),))
      #discriminant = scipy.optimize.brute(find_confidence(conf, pmean, cmean, pstd, cstd, distribution, max_yval, point_estimate=point_estimate), (slice(-0.01, max_yval, CONFIDENCE_RESOLUTION),))
      #error = find_confidence(conf, pmean, cmean, pstd, cstd, distribution, max_yval, report=True, point_estimate=point_estimate)(discriminant)
      discriminant, error = solve(conf, pmean, cmean, pstd, cstd, distribution, max_yval, point_estimate=point_estimate, first=True)
      if discriminant is None:
        logging.warn('no discriminant for confidence %.3f: %s vs %s error %.3f', conf, phenotype, control, error)
        return None
      if error > CONFIDENCE_MAXERROR or discriminant < -1 or discriminant > max_yval + CONFIDENCE_MAXERROR:
        logging.warn('failed to find confidence within error for confidence %.3f: %s vs %s error %.3f discriminant %.3f', conf, phenotype, control, error, discriminant)
        # let's just return it
        return discriminant
        #return None
      logging.debug('discriminant %f for %s vs %s at %f', discriminant, phenotype, control, conf)
      return discriminant
    except ValueError:
      logging.warn('calculate_confidence: discriminant did not converge with %s vs %s', phenotype, control)
      return None
    except RuntimeError:
      logging.warn('calculate_confidence: something went wrong with %s vs %s', phenotype, control)
      return None
 
    # this method only deals with the control group distribution and does not consider how far away the single phenotype sample is
    # so not so helpful because 1% confidence will just include all of the controls
    # we can't do lda - just calculate distance from mean in terms of sd
    #mean = np.mean(control)
    #sd = np.std(control)
    #if distribution == 'normal':
    #  if phenotype[0] < mean:
    #    #return scipy.stats.norm.ppf(1 - conf, loc=mean, scale=sd) # 2-sided
    #    return scipy.stats.norm.ppf((1 - conf) / 2, loc=mean, scale=sd) # 1-sided
    #  else:
    #    #return scipy.stats.norm.ppf(conf, loc=mean, scale=sd) # 2-sided
    #    return scipy.stats.norm.ppf(1 - ((1-conf) / 2), loc=mean, scale=sd)
    #elif distribution == 'beta':
    #  mean = mean / max_yval
    #  sd = sd / max_yval
    #
    #  n = mean * (1 - mean) / (sd ** 2)
    #  a = mean * n
    #  b = (1 - mean) * n
    #  centre = scipy.stats.beta.ppf(0.5, a, b) * max_yval
    #  if phenotype[0] < centre:
    #    return scipy.stats.beta.ppf((1 - conf) / 2, a, b) * max_yval # 1-sided
    #  else:
    #    return scipy.stats.beta.ppf(1 - ((1-conf) / 2), a, b) * max_yval

def mean(vals, sd=None):
  if sd is None:
    return np.mean(vals)
  else:
    return np.mean(vals) + np.std(vals) * sd

def ci(xs, distribution, interval=0.95, max_yval=1.0):
  if len(xs) == 1:
    return [xs[0], xs[0], xs[0]] # no interval 

  if distribution == 'normal':
    #return (mean(xs, -2), mean(xs, 2)) # TODO ignores interval
    mean = np.mean(xs)
    sd = np.std(xs)

    if abs(sd) < 1e-6:
      return [xs[0], xs[0], xs[0]] # no interval since no sd

    return scipy.stats.norm.ppf([(1-interval)/2, 0.5, 1-(1-interval)/2], loc=mean, scale=sd)
  elif distribution == 'beta':
    mean_scaled = np.mean(xs) / max_yval
    sd_scaled = np.std(xs) / max_yval

    if abs(sd_scaled) < 1e-6:
      return [xs[0], xs[0], xs[0]] # no interval since no sd

    n = mean_scaled * (1 - mean_scaled) / (sd_scaled ** 2)
    a = mean_scaled * n
    b = (1 - mean_scaled) * n

    ppf = scipy.stats.beta.ppf([(1-interval)/2, 0.5, 1-(1-interval)/2], a, b) * max_yval
    if any([math.isnan(x) for x in ppf]):
      logging.warn('ppf for %s at interval %.2f max %f is %s with mean_scaled %f sd_scaled %f', xs, interval, max_yval, ppf, mean_scaled, sd_scaled)
    else:
      logging.debug('ppf for %s at interval %.2f max %f is %s with mean_scaled %f sd_scaled %f', xs, interval, max_yval, ppf, mean_scaled, sd_scaled)
    return ppf

def plot_discriminate(data_fh, groups, x, target, filters, title, logx, signature, error_plot, count_plot, anonymise, highlight_groups, no_legend=False, confidence=None, confidence_phenotypes=set(), summarise_groups=False, x_highlight=None, distribution='normal', max_yval=1.0, width=12, height=8, fontsize=8, confidence_out=sys.stdout, point_estimate=True, confidence_in=None, markersize=6, linewidth=2.0):
  logging.info('v2. plotting signature %s with filter %s and confidence phenotype %s max_yval %.2f...', signature, filters, confidence_phenotypes, max_yval)

  import matplotlib.style
  matplotlib.style.use('seaborn')
  rcParams.update({'font.size': fontsize})

  anon_map = None
  if anonymise is not None:
    anon_map = {}
    for pair in anonymise:
      real, fake = pair.split('=')
      anon_map[real] = fake

  # Sample  Tags    Caller  DP      AF      Error   Variants        Multiplier      Signature.1     ...    Signature.30
  included = total = 0
  results = collections.defaultdict(list)
  tags = set()
  all_samples = set()
  for group in groups:
    all_samples.update(groups[group])

  logging.info('all samples: %s', all_samples)

  for row in csv.DictReader(data_fh, delimiter='\t'): # each results line
    # is it in one of the groups?
    if row['Sample'] not in all_samples:
      continue

    tags.add(row['Tags'])
    ok = True
    for f in filters:
      col, val = f.split('=')
      if val != row[col]:
        ok = False

    if ok:
      included += 1
      xval = float(row[x]) # x axis value
      if isinstance(signature, list):
        #highlight = sum([float(row['Signature.{}'.format(s)]) for s in signature])
        try:
          yval = sum([float(row[s]) for s in signature])
        except:
          logging.warn('failed to read %s', row)
          raise
      else:
        #highlight = float(row['Signature.{}'.format(signature)])
        logging.debug('signature %s has value %s', signature, row[signature])
        yval = float(row[signature])
      if row['Error'] == '':
        row['Error'] = 0
      if row['Variants'] == '':
        row['Variants'] = 0
      results[row['Sample']].append((xval, yval, float(row.get('Error', 0)), float(row.get('Variants', 0)))) # xvalue, yvalue, error, variants

    total += 1

  logging.info('finished reading %i of %i records with max yval %.2f', included, total, max_yval)

  logging.debug('tags: %s', ' '.join(sorted(list(tags))))
  logging.info('samples/values: %s', ['{}/{}'.format(sample, len(results[sample])) for sample in sorted(list(results.keys()))])

  # generate plot
  fig = plt.figure(figsize=(width, height))
  plt.rcParams.update({'font.size': fontsize})
  plt.rc('legend',fontsize=fontsize)
  if error_plot and not count_plot: # just error plot
    grid = plt.GridSpec(6, 1, hspace=0, wspace=0)
    ax = fig.add_subplot(grid[0:-1, :])
    ax_err = fig.add_subplot(grid[-1, :], sharex=ax)
  elif error_plot and count_plot: # both
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

  ax.tick_params(axis='x', labelsize=fontsize)
  ax.tick_params(axis='y', labelsize=fontsize)

  if logx:
    ax.set_xscale("log", nonposx='clip')

  if len(results) == 0:
    logging.warn('No data to plot')
    return

  groups_low = {}
  groups_high = {}
  for group in groups:
    if len(groups[group]) == 0:
      continue
    groups_low[group] = [max_yval] * len(results[groups[group][0]]) # all ones
    groups_high[group] = [0.0] * len(results[groups[group][0]]) # all zeros
    logging.debug('group range %s: %s %s', group, groups_low[group], groups_high[group])

  confidence_phenotype = collections.defaultdict(list) # [x] = { val1, val2, ...}
  confidence_control = collections.defaultdict(list)
  group_summary = {} # {'group': {10: [1, 2, 3...]}}

  names = []
  for group_num, group in enumerate(groups):
    if len(groups[group]) == 0:
      continue

    color = GROUP_COLOURS[group_num % len(GROUP_COLOURS)]
    if highlight_groups is not None and group in highlight_groups:
      alpha = HIGHLIGHT_ALPHA
    else:
      alpha = ALPHA

    # each sample in this group
    for idx, sample in enumerate(groups[group]):
      logging.info('adding group %s', group)
      if sample not in results:
        logging.warn('skipping %s', sample)
        continue

      sample_results = sorted(results[sample], key=lambda r: r[0]) # for this sample, get all in order results
      xs = [r[0] for r in sample_results]
      ys = [r[1] for r in sample_results]

      # a hack to add our new samples with just a single value
      if len(ys) == 1:
        logging.warn('missing data but STILL adding %s to group. only %i values instead of %i', sample, len(ys), len(groups_low[group])) 
        ys = [ys[0]] * len(groups_low[group])

      if len(ys) != len(groups_low[group]):
        logging.warn('missing data: not adding %s to group. only %i values instead of %i', sample, len(ys), len(groups_low[group])) 
      else:
        logging.debug('group xs %i sample xs %i', len(groups_low[group]), len(xs))
        # which group is it in
        groups_low[group] = [min(x, y) for x, y in zip(groups_low[group], ys)]
        groups_high[group] = [max(x, y) for x, y in zip(groups_high[group], ys)]
        #logging.debug('group range updated %s: %s %s', group, groups_low[group], groups_high[group])

      # values for confidence calculation
      if group in confidence_phenotypes: # phenotype
        [confidence_phenotype[x].append(y) for x, y in zip(xs, ys)]
      elif highlight_groups is not None and group in highlight_groups:
        pass # don't put highlight groups in either phenotype or control
      else: # control
        [confidence_control[x].append(y) for x, y in zip(xs, ys)]

      if group not in group_summary:
        group_summary[group] = collections.defaultdict(list)
      [group_summary[group][x].append(y) for x, y in zip(xs, ys)]
  
      if no_legend:
        label = None
      elif anon_map is not None:
        #name, rest = sample.split('_', 1)
        #if name not in names:
        #  names.append(name)
        #label = '{}_{}'.format(chr(ord('A') + names.index(name)), rest)
        label = anon_map.get(sample, sample)
      else:
        label = sample

      if not summarise_groups:
        ax.plot(xs, ys, label=label, color=color, marker=MARKERS[idx % len(MARKERS)], markersize=markersize, alpha=alpha, linewidth=linewidth)
      if summarise_groups and highlight_groups is not None and group in highlight_groups:
        ax.plot(xs, ys, label='{} {}'.format(group, label), color=color, marker=MARKERS[idx % len(MARKERS)], markersize=markersize, alpha=HIGHLIGHT_ALPHA, linewidth=linewidth*0.2)
      logging.debug('added sample %s for group %s with %i xs %i ys', label, group, len(xs), len(ys))
    
    if summarise_groups:
      xs = [r[0] for r in sorted(results[groups[group][0]])]
      cis = [ci(group_summary[group][x], distribution, interval=0.95, max_yval=max_yval) for x in xs]
      lower = [c[0] for c in cis]
      mid = [c[1] for c in cis]
      upper = [c[2] for c in cis]
      logging.info('lower %s mid %s upper %s', lower, mid, upper)
      if len(xs) == 0:
        logging.warn('xs is empty')
        continue
      if group in highlight_groups: # skip it, already plotted
        pass #ax.plot(xs, [group_summary[group][x] for x in xs], label='{}'.format(group), color=color, alpha=alpha, linewidth=HIGHLIGHT_LINE_WIDTH, marker=MARKERS[(group_num * 3) % len(MARKERS)], markersize=HIGHLIGHT_MARKER_SIZE)
      elif len(group_summary[group][xs[0]]) == 1: # just one sample, or our highlight group
        logging.info('group %s has just one sample', group)
        ax.plot(xs, [group_summary[group][x] for x in xs], label='{}'.format(group), color=color, alpha=alpha, linewidth=linewidth, marker=MARKERS[(group_num * 3) % len(MARKERS)], markersize=markersize)
      else:
        logging.info('group %s has %i samples', group, len(group_summary[group][xs[0]]))
        #ax.plot(xs, mid, label='{} 50th pctl n={}'.format(group, len(group_summary[group][xs[0]])), color=color, alpha=alpha, linewidth=MEAN_LINEWIDTH, marker=MARKERS[(group_num * 3) % len(MARKERS)], markersize=marker_size)
        ax.plot(xs, mid, label='{} 50th pctl'.format(group), color=color, alpha=alpha, linewidth=linewidth, marker=MARKERS[(group_num * 3) % len(MARKERS)], markersize=markersize)
        ax.plot(xs, upper, label='{} 95th pctl'.format(group), color=color, alpha=alpha, linewidth=linewidth, marker=MARKERS[(group_num * 3 + 1) % len(MARKERS)], markersize=markersize)
        ax.plot(xs, lower, label='{} 5th pctl'.format(group), color=color, alpha=alpha, linewidth=linewidth, marker=MARKERS[(group_num * 3 + 2) % len(MARKERS)], markersize=markersize)
        ax.fill_between(xs, lower, upper, color=color, alpha=0.2)

    if x_highlight is not None:
      ax.axvline(x_highlight, color='#ff8000', ymin=0, ymax=1, linewidth=linewidth)

    logging.info('added group %s', group)

  # now add shading for groups
  if not summarise_groups:
    for group_num, group in enumerate(groups):
      logging.debug('getting xs for group %i: %s', group_num, group)
      if len(groups[group]) == 0:
        logging.info('skipping group %s due to being empty', group)
      else:
        xs = [r[0] for r in sorted(results[groups[group][0]])]
        logging.debug('xs %i groups_low %i groups_high %i', len(xs), len(groups_low[group]), len(groups_high[group]))
        ax.fill_between(xs, groups_low[group], groups_high[group], color=GROUP_COLOURS[group_num % len(GROUP_COLOURS)], label=group, alpha=0.2)

  # confidence provided as percentages
  if confidence is not None:
    confidence_out.write('x\t{}\n'.format('\t'.join([str(c) for c in confidence])))
    confidence_results = {}

    confidence_min = [max_yval] * len(confidence_phenotype.keys())
    confidence_max = [0.0] * len(confidence_phenotype.keys())
    for conf_ix, conf in enumerate(confidence): # each confidence level
      confidence_xs = []
      confidence_ys = []
      confidence_result = {}
      logging.debug('calculating confidence %s', confidence[conf_ix])
      for ix, xval in enumerate(sorted(confidence_phenotype.keys())): # each x value
        yval = calculate_confidence(conf, confidence_phenotype[xval], confidence_control[xval], distribution, max_yval, point_estimate)
        if yval is not None and (len(confidence_ys) == 0 or yval - confidence_ys[-1] < 0.5 * max_yval): # TODO hack for now
          confidence_xs.append(xval)
          confidence_ys.append(yval)
          logging.debug('including result %s', yval)
          if isinstance(yval, float): # TODO hack for now
            confidence_result[str(xval)] = yval
          else:
            confidence_result[str(xval)] = yval[0]
          confidence_min[ix] = min(confidence_min[ix], confidence_ys[-1])
          confidence_max[ix] = max(confidence_max[ix], confidence_ys[-1])
      ax.plot(confidence_xs, confidence_ys, label='{:.0f}% confidence'.format(conf * 100), color=CONFIDENCE_COLORS[conf_ix % len(CONFIDENCE_COLORS)], alpha=CONFIDENCE_ALPHA, linewidth=linewidth * 1.5)
      confidence_results[str(conf)] = confidence_result
    # and the range
    #ax.fill_between(sorted(confidence_phenotype.keys()), confidence_min, confidence_max, color=CONFIDENCE_COLOR, alpha=0.2)
    
    # write confidence results
    logging.debug('confidence results {}'.format(confidence_results))
    for ix, xval in enumerate(sorted(confidence_phenotype.keys())): # each x value
      logging.debug('writing xval={} for confs={}'.format(xval, confidence))
      confidence_out.write('{}\t{}\n'.format(xval, '\t'.join(['{:.2f}'.format(confidence_results[str(conf)].get(str(xval), -1)) for conf in confidence])))

  # confidence provided as file
  if confidence_in is not None:
    rows = [x.strip('\n').split('\t') for x in open(confidence_in, 'r').readlines()]
    intervals = [float(x) for x in rows[0][1:]]
    for col, interval in enumerate(intervals):
      confidence_xs = [float(row[0]) for row in rows[1:] if float(row[col + 1]) >= -0.1]
      confidence_ys = [float(row[col + 1]) for row in rows[1:] if float(row[col + 1]) >= -0.1]
      logging.debug('confidence_xs %s confidence_ys %s', confidence_xs, confidence_ys)
      ax.plot(confidence_xs, confidence_ys, label='{:.0f}% confidence'.format(interval * 100), color=CONFIDENCE_COLORS[col % len(CONFIDENCE_COLORS)], alpha=CONFIDENCE_ALPHA, linewidth=linewidth * 1.5)

  # additional annotation
  if isinstance(signature, list):
    ax.set_ylabel('Sum of signatures {}'.format(' '.join(signature)), fontsize=fontsize)
  else:
    ax.set_ylabel('Signature {} proportion'.format(signature), fontsize=fontsize)

  ax.set_xlabel(LABELS[x], fontsize=fontsize)
  if title is None:
    ax.set_title('Separation of subtypes across {} using signature {} ({})'.format(x, signature, ' '.join(filters)), fontsize=fontsize)
  else:
    ax.set_title(title, fontsize=fontsize) #'Association of signature {} across {} ({})'.format(signature, x, title))

  # place legend at right based on https://stackoverflow.com/questions/10101700/moving-matplotlib-legend-outside-of-the-axis-makes-it-cutoff-by-the-figure-box/10154763#10154763
  handles, labels = ax.get_legend_handles_labels()
  lgd = ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(1.01,1.0), borderaxespad=0)
  lgd.get_frame().set_edgecolor('#000000')

  if error_plot:
    representative = list(all_samples)[0]
    sample_results = sorted(results[representative], key=lambda r: r[0]) # for this sample, get all in order results
    xs = [r[0] for r in sample_results]
    ys = [r[2] for r in sample_results] 
    ax_err.plot(xs, ys, 'k-', linewidth=linewidth * 0.25)
    ax_err.fill_between(xs, 0, ys, color='#a0a0a0')
    ax_err.grid(True)
    ax_err.set_ylabel('Error', fontsize=fontsize)
    ax_err.set_ylim((0, 1.0))
    ax_err.set_xlabel(x, fontsize=fontsize)

  if count_plot:
    representative = list(all_samples)[0]
    sample_results = sorted(results[representative], key=lambda r: r[0]) # for this sample, get all in order results
    xs = [r[0] for r in sample_results]
    ys = [r[3] for r in sample_results] 
    ax_count.plot(xs, ys, 'k-', linewidth=linewidth * 0.25)
    ax_count.fill_between(xs, 0, ys, color='#a0a0a0')
    ax_count.grid(True)
    ax_count.set_ylabel('Mutations', fontsize=fontsize)
    ax_count.set_xlabel(x, fontsize=fontsize)
    #ax_count.set_xscale("log", nonposx='clip')

  #plt.savefig(target)
  fig.savefig(target, bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=DPI)
  matplotlib.pyplot.close('all')
  logging.info('done processing %i of %i', included, total)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Plot changes in signature')
  parser.add_argument('--data', required=True, help='data file')
  parser.add_argument('--groups', nargs='+', required=True, help='group1=sample1,sample2 group2=sample3,sample4')
  parser.add_argument('--highlight_groups', nargs='*', required=False, help='groups to highlight')
  parser.add_argument('--summarise_groups', action='store_true', required=False, help='show mean +- 2sd')
  parser.add_argument('--filters', nargs='*', help='other filters')
  parser.add_argument('--highlight', required=True, help='signature(s) to highlight')
  parser.add_argument('--x_highlight', required=False, type=float, help='vertical line at this value')
  parser.add_argument('--x', required=True, help='x column name')
  parser.add_argument('--logx', action='store_true', help='log x')
  parser.add_argument('--error_plot', action='store_true', help='include erorr plot')
  parser.add_argument('--count_plot', action='store_true', help='include count plot')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--anonymise', required=False, nargs='+', help='map samples to new names')
  parser.add_argument('--no_legend', action='store_true', help='no samples in legend')
  parser.add_argument('--title', required=False, help='sample filter')
  parser.add_argument('--target', required=False, default='plot.png', help='plot filename')
  parser.add_argument('--confidence', required=False, nargs='*', type=float, help='show confidence intervals e.g. 0.01 0.99')
  parser.add_argument('--confidence_phenotypes', required=False, nargs='*', default=set(), help='which groups are phenotype')
  parser.add_argument('--confidence_in', required=False, help='file containing existing confidence calculations')
  parser.add_argument('--distribution', required=False, default='normal', help='assumption of distribution [normal, lognormal]')
  parser.add_argument('--max_yval', required=False, type=float, default=1.0, help='vertical line at this value')
  parser.add_argument('--markersize', required=False, type=float, default=6, help='size of markers')
  parser.add_argument('--linewidth', required=False, type=float, default=2, help='line width')
  parser.add_argument('--height', required=False, type=float, default=8, help='height of plot')
  parser.add_argument('--width', required=False, type=float, width=12, help='width of plot')
  parser.add_argument('--fontsize', required=False, default=12, type=int, help='plot font size')
  parser.add_argument('--point_estimate', action='store_true', help='calculate point estimate confidence levels')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  groups = collections.defaultdict(list)
  for item in groups:
    group, samples = item.split('=')
    groups[group] = samples.split(',')

  plot_discriminate(open(args.data, 'r'), groups, args.x, args.target, args.filters, args.title, args.logx, args.highlight, args.error_plot, args.count_plot, args.anonymise, args.highlight_groups, args.no_legend, args.confidence, args.confidence_phenotypes, args.summarise_groups, args.x_highlight, args.distribution, args.max_yval, args.height, args.width, args.fontsize, sys.stdout, args.point_estimate, args.confidence_in, args.markersize, args.linewidth)

