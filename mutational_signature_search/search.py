#!/usr/bin/env python
'''
  run multiple signature calculations
'''

import argparse
import logging
import io
import sys

import cyvcf2

import mutational_signature.count
import mutational_signature.decompose

def filter_mutect2(sample, dp, af, use_bam_depth, pass_only, use_vaf):
  def filter_mutect2_instance(vcf, variant):
    if pass_only and variant.FILTER is not None:
      return False
    sample_id = vcf.samples.index(sample)
    if use_bam_depth is not None:
      #total_depth = variant.format('DP')[sample_id] # somatic
      total_depth = variant.INFO[use_bam_depth]
    else:
      total_depth = variant.INFO['DP'] # somatic + germline depth
    #vcf_af = depths[1] / sum(depths)
    if use_vaf is not None:
      if use_vaf == 'sample-AD':
        ad = variant.format('AD')[sample_id]
        ref = ad[0]
        alt = ad[1]
        if ref + alt > 0:
          vcf_af = alt / (ref + alt)
        else:
          vcf_af = 0
      else:
        vcf_af = variant.INFO[use_vaf]
    else:
      vcf_af = variant.format('AF')[sample_id] # mutect2 af
    return total_depth >= dp and vcf_af >= af

  return filter_mutect2_instance

def filter_strelka(sample, dp, af, use_bam_depth, pass_only, use_vaf):
  def filter_strelka_instance(vcf, variant):
    if pass_only and variant.FILTER is not None:
      return False
    sample_id = vcf.samples.index('TUMOR')
    if use_bam_depth is not None:
      total_depth = variant.INFO[use_bam_depth]
      #depths = variant.format('AD')[sample_id]
      #total_depth = sum(depths) # somatic depth
    else:
      total_depth = variant.INFO['DP'] # somatic + germline depth
    if use_vaf is not None:
      vcf_af = variant.INFO[use_vaf]
    else:
      vcf_af = variant.INFO['AF']
    return total_depth >= dp and vcf_af >= af
  return filter_strelka_instance

def filter_maf(sample, dp, af, use_bam_depth, pass_only, use_vaf, use_alt_depth):
  def filter_maf_instance(vcf, variant):
    if pass_only and variant.FILTER is not None and variant.FILTER != "PASS":
      return False
    # sample_id = vcf.samples.index('TUMOR') # maf is already filtered for sample of interest
    if use_bam_depth is not None:
      total_depth = float(variant.INFO[use_bam_depth])
      #depths = variant.format('AD')[sample_id]
      #total_depth = sum(depths) # somatic depth
    else:
      total_depth = float(variant.INFO['DP']) # somatic + germline depth
    if use_vaf is not None:
      vcf_af = float(variant.INFO[use_vaf])
    elif use_alt_depth is not None:
      vcf_af = float(variant.INFO[use_alt_depth]) / max(1, total_depth)
    else:
      vcf_af = float(variant.INFO['AF'])
    return total_depth >= dp and vcf_af >= af
  return filter_maf_instance

def main(genome, signatures, vcfs, dps, afs, context_cutoff, caller, tags, use_bam_depth, doublets, indels, just_indels, pass_only, sample, use_vaf, weight, maf_sample, maf_sample_column, maf_chrom_column, maf_pos_column, maf_ref_column, maf_alt_column, vcf_not_zipped, use_alt_depth, maf_filter_column):
  logging.info('starting...')

  chroms = {}

  first = True
  sample_index = 0
  for vcf in vcfs:
    outs = []
    variant_filters = []
    if maf_sample is not None and len(maf_sample) > 0:
      sample = maf_sample[sample_index]
    elif sample is None:
      sample = vcf.split('/')[-1].split('.')[0] # simplified sample name
    #sample = vcf
    for dp in dps:
      for af in afs:
        # calculate counts
        logging.debug('calculating counts for %s with dp %i af %.3f caller %s...', vcf, dp, af, caller)
        variant_filter = None

        if caller == 'mutect2':
          variant_filter = filter_mutect2(sample, dp, af, use_bam_depth, pass_only, use_vaf)

        elif caller == 'strelka':
          variant_filter = filter_strelka(sample, dp, af, use_bam_depth, pass_only, use_vaf)

        elif caller == 'maf':
          variant_filter = filter_maf(sample, dp, af, use_bam_depth, pass_only, use_vaf, use_alt_depth)

        else:
          logging.fatal('unrecognized caller %s', caller)

        out = io.StringIO()
        variant_filters.append(variant_filter)
        outs.append(out)

    # do all at once
    # def multi_count(genome_fh, vcf_in, outs=None, chroms=None, variant_filters=None, doublets=False, indels=False, just_indels=False)
    if maf_sample is not None and len(maf_sample) > 0:
      vcf_in = mutational_signature.count.maf_to_vcf(vcf, args.maf_sample[sample_index], args.maf_sample_column, args.maf_chrom_column, args.maf_pos_column, args.maf_ref_column, args.maf_alt_column, args.vcf_not_zipped, maf_filter_column)
    else:
      vcf_in = cyvcf2.VCF(vcf)

    all_counts = mutational_signature.count.multi_count(open(genome, 'r'), vcf_in, outs=outs, chroms=chroms, variant_filters=variant_filters, doublets=doublets, indels=indels, just_indels=just_indels, weight=weight)
    chroms = all_counts['chroms'] # keep track of updates

    idx = 0
    for dp in dps:
      for af in afs:
        counts = all_counts['all_counts'][idx]
        totals = all_counts['all_totals'][idx]
        logging.debug('calculating signature for %s with dp %i af %.3f with %i variants', vcf, dp, af, totals)
        out = outs[idx]
        idx += 1

        logging.debug('decomposing %s', out.getvalue())
        out.seek(0)
        result = mutational_signature.decompose.decompose(open(signatures, 'r'), out, None, 'cosine', None, None, 'basin', None, context_cutoff, False)

        if first:
          first = False
          sys.stdout.write('Sample\tTags\tCaller\tDP\tAF\tError\tVariants\tMultiplier\t{}\n'.format('\t'.join(result['signature_names'])))

        #sys.stdout.write('{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}\t{}\t{:.3f}\t{}\n'.format(sample, tags, caller, dp, af, result['error'], counts['total'], result['total'], '\t'.join(['{:.3f}'.format(x / result['total']) for x in result['signature_values']])))
        sys.stdout.write('{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}\t{}\t{:.3f}\t{}\n'.format(sample, tags, caller, dp, af, result['error'][0], totals, result['total_included'], '\t'.join(['{:.3f}'.format(x / result['total']) for x in result['signature_values']])))

        logging.info('calculating signature for %s with dp %i af %.3f: error %.3f total %.3f', vcf, dp, af, result['error'][0], result['total_included'])

    sample_index += 1
  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Signature search')

  parser.add_argument('--genome', required=True, help='reference genome')
  parser.add_argument('--signatures', required=True, help='context counts')
  parser.add_argument('--vcfs', required=True, nargs='+', help='vcfs')
  parser.add_argument('--caller', required=True, help='variant caller mutect2 or strelka or bam')
  parser.add_argument('--tags', required=False, default='', help='tags to include in output')
  parser.add_argument('--dps', required=True, nargs='+', type=int, help='depth settings')
  parser.add_argument('--afs', required=True, nargs='+', type=float, help='af settings')
  parser.add_argument('--use_bam_depth', required=False, help='use this info field for tumour depth e.g. BAM_DEPTH')
  parser.add_argument('--use_alt_depth', required=False, help='use this info field for alt depth')
  parser.add_argument('--use_vaf', required=False, help='use this info field for vaf e.g. VAF')
  parser.add_argument('--weight', required=False, help='field name to weight counts by')
  parser.add_argument('--pass_only', action='store_true', help='just pass variants')
  parser.add_argument('--context_cutoff', required=False, type=float, default=1e6, help='exclude signatures with contexts above this percent that are not represented in the sample') # deconstructSigs = 0.2
  parser.add_argument('--doublets', action='store_true', help='generate doublets')
  parser.add_argument('--indels', action='store_true', help='generate indels')
  parser.add_argument('--just_indels', action='store_true', help='generate only indels')
  parser.add_argument('--sample', required=False, help='sample name in vcf (else inferred from filename)')
  parser.add_argument('--maf_sample', required=False, nargs='+', help='vcf is actually a maf with these samples of interest')
  parser.add_argument('--maf_sample_column', required=False, default='Tumor_Sample_Barcode', help='maf chrom column name')
  parser.add_argument('--maf_chrom_column', required=False, default='Chromosome', help='maf chrom column name')
  parser.add_argument('--maf_pos_column', required=False, default='Start_Position', help='maf pos column name')
  parser.add_argument('--maf_ref_column', required=False, default='Reference_Allele', help='maf ref column name')
  parser.add_argument('--maf_alt_column', required=False, default='Tumor_Seq_Allele2', help='maf alt column name')
  parser.add_argument('--maf_filter_column', required=False, default='FILTER', help='maf alt column name')
  parser.add_argument('--vcf_not_zipped', action='store_true', help='for maf handling ?')
  parser.add_argument('--verbose', action='store_true', help='more logging')

  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.genome, args.signatures, args.vcfs, args.dps, args.afs, args.context_cutoff, args.caller, args.tags, args.use_bam_depth, args.doublets, args.indels, args.just_indels, args.pass_only, args.sample, args.use_vaf, args.weight, args.maf_sample, args.maf_sample_column, args.maf_chrom_column, args.maf_pos_column, args.maf_ref_column, args.maf_alt_column, args.vcf_not_zipped, args.use_alt_depth, args.maf_filter_column)
