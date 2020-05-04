#!/usr/bin/env python


import argparse
import logging
import sys

import numpy as np

import mutational_signature_search.plot_discriminate

def main(phenotype, control, confs):
  logging.info('starting...')
  start = (np.mean(phenotype) + np.mean(control)) / 2
  pmean = np.mean(phenotype)
  cmean = np.mean(control)
  pstd = np.std(phenotype)
  cstd = np.std(control)
  distribution = 'beta'
  max_yval = 1.0

  sys.stdout.write('Conf\tDiscriminant\tError\n')
  for conf in confs:
    discriminant, error = mutational_signature_search.plot_discriminate.solve(conf, pmean, cmean, pstd, cstd, distribution, max_yval, point_estimate=True, first=True)
    sys.stdout.write('{:.3f}\t{:.3f}\t{:.3f}\n'.format(conf, discriminant, error))
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--phenotype', required=True, nargs='+', type=float, help='phenotype values')
  parser.add_argument('--control', required=True, nargs='+', type=float, help='control values')
  parser.add_argument('--conf', required=True, nargs='+', type=float, help='confidences to calculate')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.phenotype, args.control, args.conf)
