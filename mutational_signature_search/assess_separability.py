#!/usr/bin/env python

import argparse
import logging
import math
import sys

import numpy as np
import scipy.stats
import sklearn.discriminant_analysis
import sklearn.metrics

def t_test(a, b, one_sided=False):
  if len(a) == 1 or len(b) == 1:
    return -1 # cannot calculate
  if one_sided:
    if np.mean(a) < np.mean(b):
      return 1
    else:
      return scipy.stats.ttest_ind(a, b)[1] / 2
  else:
    return scipy.stats.ttest_ind(a, b)[1]

def z_test(a, b):
  '''
    a is a single item, b is a list
    this is a one sided test as we assume the phenotype will show increased value of a signature
  '''
  z = (a - np.mean(b)) / np.std(b)
  #return 1 - scipy.stats.norm.cdf(z)
  return scipy.stats.norm.sf(z) # more accurate

def linear_discriminant_ratio(a, b):
  if np.mean(a) == np.mean(b):
    return 0 # no good

  # if all equal, add change
  if len(a) > 1 and all(a[0] == x for x in a):
    a = a.copy()
    a[0] += 0.001
  if len(b) > 1 and all(b[0] == x for x in b):
    b = b.copy()
    b[0] += 0.001
  
  # do grubb instead
  if len(a) == 1 or len(b) == 1:
    if len(a) == 1:
      mean_diff = (a - np.mean(b)) ** 2
      var = 2 * np.var(b)
    elif len(b) == 1:
      mean_diff = (b - np.mean(a)) ** 2
      var = 2 * np.var(a)

    # enforce a minimum difference
    if mean_diff < 0.001:
      return -1

    return float(mean_diff / var)

  # difference in means
  delta_mean = (np.mean(a) - np.mean(b)) ** 2
  var_sum = np.var(a) + np.var(b)

  # enforce a minimum difference
  if delta_mean < 0.001:
    return -1
  
  return float(delta_mean / var_sum) 

def separable(a, b):
  return min(a) > max(b) or min(b) > max(a)

def linear_discriminant_accuracy(a, b):
  if len(a) == 1 or len(b) == 1:
    return (-1, 1.0, 1.0) # cannot calculate
  xs = [[x] for x in a] + [[x] for x in b]
  ys = ['a'] * len(a) + ['b'] * len(b)
  logging.debug('a %s b %s', a, b)

  # sanity checks
  if all(xs[0] == x for x in xs):
    logging.info('linear_discriminant_accuracy is 0.0: all values identical')
    return (0.0, 1.0, 1.0)

  classifier = sklearn.discriminant_analysis.LinearDiscriminantAnalysis()
  classifier.fit(xs, ys)
  discriminant = -classifier.intercept_  / classifier.coef_
  # predicted % error
  if np.mean(a) > np.mean(b):
    a_error = scipy.stats.norm(np.mean(a), np.std(a)).cdf(discriminant)
    b_error = 1 - scipy.stats.norm(np.mean(b), np.std(b)).cdf(discriminant)
  else:
    a_error = 1 - scipy.stats.norm(np.mean(a), np.std(a)).cdf(discriminant)
    b_error = scipy.stats.norm(np.mean(b), np.std(b)).cdf(discriminant)

  error = (a_error * len(a) + b_error * len(b)) / (len(a) + len(b))
  accuracy = float(1 - error)
  logging.debug('linear discriminant %.2f a_error %.2f b_error %.2f accuracy %.2f', discriminant, a_error, b_error, accuracy)

  return (accuracy, float(a_error), float(b_error))

def find_discriminant(classifier):
  def f(x):
    return classifier.decision_function([[x]])[0]
  return f

def quadratic_discriminant_accuracy(a, b):
  if len(a) == 1 or len(b) == 1:
    return (-1, 1.0, 1.0)  # cannot calculate

  xs = [[x] for x in a] + [[x] for x in b]
  ys = ['a'] * len(a) + ['b'] * len(b)
  logging.debug('a %s b %s', a, b)

  # sanity checks
  if all(xs[0] == x for x in xs):
    logging.info('quadratic_discriminant_accuracy is 0.0: all values identical')
    return (0.0, 1.0, 1.0)

  classifier = sklearn.discriminant_analysis.QuadraticDiscriminantAnalysis()
  classifier.fit(xs, ys)
  try:
    discriminant = scipy.optimize.newton(lambda x: classifier.decision_function([[x]])[0], (np.mean(a) + np.mean(b)) / 2)
  except ValueError:
    logging.info('quadratic_discriminant_accuracy: discriminant did not converge')
    return (-1, 1.0, 1.0)
  except RuntimeError:
    logging.info('quadratic_discriminant_accuracy: something went wrong')
    return (-1, 1.0, 1.0)
  #discriminant = scipy.optimize.newton(find_discriminant(classifier), (np.mean(a) + np.mean(b)) / 2)
  # predicted % error
  if np.mean(a) > np.mean(b):
    a_error = scipy.stats.norm(np.mean(a), np.std(a)).cdf(discriminant)
    b_error = 1 - scipy.stats.norm(np.mean(b), np.std(b)).cdf(discriminant)
  else:
    a_error = 1 - scipy.stats.norm(np.mean(a), np.std(a)).cdf(discriminant)
    b_error = scipy.stats.norm(np.mean(b), np.std(b)).cdf(discriminant)

  error = (a_error * len(a) + b_error * len(b)) / (len(a) + len(b))
  accuracy = float(1 - error)
  logging.debug('quadratic discriminant %.2f a_error %.2f b_error %.2f accuracy %.2f', discriminant, a_error, b_error, accuracy)

  return (accuracy, float(a_error), float(b_error))

def variation_discriminant_accuracy(a, b):
  # choose our own separation point
  if len(a) == 1 or len(b) == 1:
    return (-1, 1.0, 1.0)  # cannot calculate

  delta_mean = (np.mean(a) - np.mean(b)) ** 2
  var_sum = np.var(a) + np.var(b)

  discriminant = np.mean(a) + delta_mean / var_sum

  if np.mean(a) > np.mean(b):
    a_error = scipy.stats.norm(np.mean(a), np.std(a)).cdf(discriminant)
    b_error = 1 - scipy.stats.norm(np.mean(b), np.std(b)).cdf(discriminant)
  else:
    a_error = 1 - scipy.stats.norm(np.mean(a), np.std(a)).cdf(discriminant)
    b_error = scipy.stats.norm(np.mean(b), np.std(b)).cdf(discriminant)

  error = (a_error * len(a) + b_error * len(b)) / (len(a) + len(b))
  accuracy = float(1 - error)
  
  return (accuracy, float(a_error), float(b_error))

def measure_accuracy(a, b, v):
  # determine tp, tn, fp, fn at point v, given our set of a and b (assuming a to be the case)
  tp = sum([1 for x in a if x >= v]) # good case
  fp = sum([1 for x in b if x >= v]) # incorrect control
  tn = sum([1 for x in b if x < v]) # good control
  fn = sum([1 for x in a if x < v]) # incorrect case

  if tp + fn > 0:
    sensitivity = tp / (tp + fn)
  else:
    senitivity = 0
  if tn + fp > 0:
    specificity = tn / (tn + fp)
  else:
    specificity = 0

  if tp + fp > 0:
    ppv = tp / (tp + fp)
  else:
    ppv = 0
  if tn + fn > 0:
    npv = tn / (tn + fn)
  else:
    npv = 0

  accuracy = (tp + tn) / (tp + tn + fp + fn)

  return {'tp': tp, 'fp': fp, 'tn': tn, 'fn': fn, 'sensitivity': sensitivity, 'specificity': specificity, 'accuracy': accuracy, 'ppv': ppv, 'npv': npv}

def auc_ci(a, b):
  '''
    Hanley and McNeil 1982
    https://homes.cs.washington.edu/~jfogarty/publications/gi2005.pdf
  '''
  truth = [1] * len(a) + [0] * len(b)
  values = a + b
  auc = sklearn.metrics.roc_auc_score(truth, values)

  if auc == 1:
    # based on tables 1 and 2 from Obuchowski, N. A., & Lieber, M. L. (2002). Confidence Bounds When the Estimated ROC Area is 1.0. Academic Radiology, 9(5), 526–530. doi:10.1016/s1076-6332(03)80329-x 
    if len(a) >= 30 and len(b) >= 30:
      return (0.98, 1.0)
    if len(a) >= 15 and len(b) >= 15:
      # n,m n=cases
      table1 = {
        '15,15': 0.93, '20,15': 0.94, '25,15': 0.95, '30,15': 0.96,
        '15,20': 0.95, '20,20': 0.96, '25,20': 0.97, '30,20': 0.97,
        '15,25': 0.95, '20,25': 0.96, '25,25': 0.97, '30,25': 0.98,
        '15,30': 0.96, '20,30': 0.97, '25,30': 0.98, '30,30': 0.98
      }
      n = min(len(a) - len(a) % 5, 30)
      m = min(len(b) - len(b) % 5, 30)
      return (table1['{},{}'.format(n,m)], 1.0)
    else:
      n = min(len(a), 30)
      m = min(len(b), 30)
      n -= n % 5
      m -= m % 5
      if n > 0 and m > 0:
        table2 = {
          '5,5': 0.72, '5,10': 0.82, '5,15': 0.87, '5,20': 0.89, '5,25': 0.90, '5,30': 0.91,
          '10,5': 0.82, '10,10': 0.90, '10,15': 0.93, '10,20': 0.94, '10,25': 0.95, '10,30': 0.96,
          '15,5': 0.87, '15,10': 0.93,
          '20,5': 0.89, '20,10': 0.94,
          '25,5': 0.90, '25,10': 0.95,
          '30,5': 0.92, '30,10': 0.96
        }
        return (table2['{},{}'.format(n,m)], 1.0)
      else:      
        return (0.0, 1.0)

  else:
    dp = (len(a) - 1) * (auc/(2-auc)-auc*auc)
    dn = (len(b) - 1) * (2 * auc * auc / (1 + auc) - auc * auc)
    interim = (auc * (1-auc) + dp + dn) / (len(a) * len(b))
    if interim > 0:
      se = math.sqrt(interim)
    else:
      se = 0
    logging.debug('se of auc is %.3f', se)
    low = max(0, auc - 1.96 * se)
    high = min(1, auc + 1.96 * se)
    return (low, high)

def auc_sklearn(a, b):
  truth = [1] * len(a) + [0] * len(b)
  values = a + b
  auc = sklearn.metrics.roc_auc_score(truth, values)
  ci = auc_ci(a, b)
  return {'auc': auc, 'ci': ci}
  
def auc(a, b):
  points = sorted(a + b)
  fprs = []
  tprs = []
  area = 0
  for point in points:
    accuracy = measure_accuracy(a, b, point)
    if accuracy['tp'] + accuracy['fn'] > 0:
      tpr = accuracy['tp'] / (accuracy['tp'] + accuracy['fn'])
    else:
      tpr = 0
    if accuracy['tn'] + accuracy['fp'] > 0:
      fpr = accuracy['fp'] / (accuracy['tn'] + accuracy['fp'])
    else:
      fpr = 0
    if len(fprs) > 0 and fpr < fprs[-1]:
      area += (tprs[-1] + tpr) / 2 * (fprs[-1] - fpr)
    fprs.append(fpr)
    tprs.append(tpr)
  ci = auc_ci(a, b)
  return {'fprs': fprs, 'tprs': tprs, 'auc': area, 'ci': ci}

def write_accuracy(a, b, vs, out):
  out.write('threshold\taccuracy\tsensitivity\tspecificity\tppv\tnpv\tTP\tTN\tFP\tFN\n')
  for v in vs:
    r = measure_accuracy(a, b, v)
    out.write('{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{}\t{}\t{}\t{}\n'.format(v, r['accuracy'], r['sensitivity'], r['specificity'], r['ppv'], r['npv'], r['tp'], r['tn'], r['fp'], r['fn']))


#def calculate_auc(a, b):
if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Measure accuracy etc')
  parser.add_argument('--phenotype', type=float, nargs='+', required=True, help='phenotype')
  parser.add_argument('--control', type=float, nargs='+', required=True, help='control')
  parser.add_argument('--thresholds', type=float, nargs='+', required=True, help='thresholds')
  args = parser.parse_args()  
  write_accuracy(args.phenotype, args.control, args.thresholds, sys.stdout)
