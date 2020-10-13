#!/usr/bin/env python

import logging

import numpy as np
import scipy.stats
import sklearn.discriminant_analysis

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
  return {'tp': tp, 'fp': fp, 'tn': tn, 'fn': fn}

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
  return {'fprs': fprs, 'tprs': tprs, 'auc': area}

#def calculate_auc(a, b):
  
