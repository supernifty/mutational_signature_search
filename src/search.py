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
    #total_depth = sum(depths) # somatic depth
    total_depth = variant.INFO['DP'] # somatic + germline depth
    #vcf_af = depths[1] / sum(depths)
    vcf_af = variant.format('AF')[sample_id] # mutect2 af
    return total_depth >= dp and vcf_af > af

  return filter_mutect2_instance

def filter_strelka(sample, dp, af):
  def filter_strelka_instance(vcf, variant):
    sample_id = vcf.samples.index('TUMOR')
    total_depth = variant.INFO['DP'] # somatic + germline depth
    vcf_af = variant.INFO['AF']
    return total_depth >= dp and vcf_af > af

def main(genome, signatures, vcfs, dps, afs, context_cutoff, caller, tags):
  logging.info('starting...')

  chroms = {}
  first = True
  for vcf in vcfs:
    for dp in dps:
      for af in afs:
        # calculate counts
        logging.info('calculating counts for %s with dp %i af %.2f...', vcf, dp, af)
        sample = vcf.split('/')[-1].split('.')[0]

        if caller == 'mutect2':
          variant_filter = filter_mutect2(sample, dp, af)

        if caller == 'strelka':
          variant_filter = filter_strelka(sample, dp, af)

        out = io.StringIO()
        counts = mutational_signature.count.count(open(genome, 'r'), vcf, out=out, chroms=chroms, variant_filter=variant_filter)
        chroms = counts['chroms'] # keep track of updates
        logging.info('calculating signature for %s with dp %i af %.2f with %i variants', vcf, dp, af, counts['total'])

        out.seek(0)
        result = mutational_signature.decompose.decompose(signatures, out, None, 'cosine', None, None, 'basin', None, context_cutoff)

        if first:
          first = False
          sys.stdout.write('Sample\tTags\tCaller\tDP\tAF\tError\tVariants\tMultiplier\t{}\n'.format('\t'.join(result['signature_names'])))

        sys.stdout.write('{}\t{}\t{}\t{}\t{:.2f}\t{:.3f}\t{}\t{:.3f}\t{}\n'.format(sample, tags, caller, dp, af, result['error'], counts['total'], result['total'], '\t'.join(['{:.3f}'.format(x / result['total']) for x in result['signature_values']])))

        logging.info('calculating signature for %s with dp %i af %.2f: error %.2f total %.2f', vcf, dp, af, result['error'], result['total'])

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Signature search')

  parser.add_argument('--genome', required=True, help='reference genome')
  parser.add_argument('--signatures', required=True, help='reference genome')
  parser.add_argument('--vcfs', required=True, nargs='+', help='vcfs')
  parser.add_argument('--caller', required=True, help='variant caller mutect2 or strelka')
  parser.add_argument('--tags', required=False, default='', help='tags to include in output')
  parser.add_argument('--dps', required=True, nargs='+', type=int, help='depth settings')
  parser.add_argument('--afs', required=True, nargs='+', type=float, help='af settings')
  parser.add_argument('--context_cutoff', required=False, type=float, default=1e6, help='exclude signatures with contexts above this percent that are not represented in the sample') # deconstructSigs = 0.2
  parser.add_argument('--verbose', action='store_true', help='more logging')

  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.genome, args.signatures, args.vcfs, args.dps, args.afs, args.context_cutoff, args.caller, args.tags)
