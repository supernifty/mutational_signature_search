#!/usr/bin/env python
'''
  run multiple signature calculations
'''

import argparse
import logging
import io
import sys

import mutational_signature.count
import mutational_signature.decompose

def filter_mutect2(sample, dp, af):
  def filter_mutect2_instance(vcf, variant):
    sample_id = vcf.samples.index(sample)
    depths = variant.format('AD')[sample_id]
    total_depth = sum(depths)
    vcf_af = depths[1] / total_depth
    return total_depth >= dp and vcf_af > af

  return filter_mutect2_instance

def main(genome, signatures, vcfs, dps, afs):
  logging.info('starting...')

  chroms = {}
  first = True
  for vcf in vcfs:
    for dp in dps:
      for af in afs:
        # calculate counts
        logging.info('calculating counts for %s with dp %i af %.2f...', vcf, dp, af)
        sample = vcf.split('/')[-1].split('.')[0]

        if 'mutect2' in vcf:
          caller = 'mutect2'
          variant_filter = filter_mutect2(sample, dp, af)

        out = io.StringIO()
        counts = mutational_signature.count.count(open(genome, 'r'), vcf, out=out, chroms=chroms, variant_filter=variant_filter)
        chroms = counts['chroms'] # keep track of updates
        logging.info('calculating signature for %s with dp %i af %.2f with %i variants', vcf, dp, af, counts['total'])

        out.seek(0)
        result = mutational_signature.decompose.decompose(signatures, out, None, 'cosine', None, None, 'basin', None, 0.2)

        if first:
          first = False
          sys.stdout.write('Sample\tCaller\tDP\tAF\tError\tVariants\tNormalizer\t{}\n'.format('\t'.join(result['signature_names'])))

        sys.stdout.write('{}\t{}\t{}\t{:.2f}\t{:.3f}\t{}\t{}\t{}\n'.format(sample, caller, dp, af, result['error'], counts['total'], result['total'], '\t'.join(['{:.3f}'.format(x / result['total']) for x in result['signature_values']])))

        logging.info('calculating signature for %s with dp %i af %.2f: error %.2f total %.2f', vcf, dp, af, result['error'], result['total'])

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Signature search')

  parser.add_argument('--genome', required=True, help='reference genome')
  parser.add_argument('--signatures', required=True, help='reference genome')
  parser.add_argument('--vcfs', required=True, nargs='+', help='vcfs')
  parser.add_argument('--dps', required=True, nargs='+', type=int, help='depth settings')
  parser.add_argument('--afs', required=True, nargs='+', type=float, help='af settings')
  parser.add_argument('--verbose', action='store_true', help='more logging')

  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.genome, args.signatures, args.vcfs, args.dps, args.afs)
