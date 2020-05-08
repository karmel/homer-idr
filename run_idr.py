'''
Created on Aug 4, 2014

@author: karmel

Irreproducibility Discovery Rate is a method used by Encode and others
to analyze the concordance of chip-seq replicates and to select the peaks
that are replicable beyond a reasonable threshold.

Here we implement several processes for computing the IDR statistic over Homer
peak files and selecting peaks of interest to us. We leverage the R package 
for computing IDR that is referenced by Encode:

https://sites.google.com/site/anshulkundaje/projects/idr 

'''
from argparse import ArgumentParser
import math
import os
from random import randint

from idr.idr_caller import IdrCaller
from idr.utils import IdrUtilities
class IdrArgumentParser(ArgumentParser):
    def __init__(self):
        description = '''Functions for running Irreproducibility Discovery Rate
        (IDR) analysis on Homer peak files.'''
        
        super().__init__(description=description)
        
        self.add_argument('command', 
                help='Program to run; options are: idr, '
                + 'pseudoreplicate, pool-pseudoreplicates, '
                + 'homer2narrow, truncate.'),
        
        self.add_argument('-o','--output_dir', nargs='?', dest='output_dir',
                help='Directory name in which output files will be placed. ' +
                'Will be created if it does not exist.')
        
        
        self.add_argument('-d', '--tag_dirs', nargs='*', dest='tag_dirs',
                help='Space-separated list of input Homer tag directories.')
        self.add_argument('-p', '--peak_files', nargs='*', dest='peak_files',
                help='Space-separated list of input Homer peak files.')
        self.add_argument('-pr', '--pseudorep_files', nargs='*', dest='pseudorep_files',
                help='Space-separated list of pseudoreplicate peak files for IDR analysis.')
        self.add_argument('-ppr', '--pooled_pseudoreps', nargs='*', dest='pooled_pseudoreps',
                help='Space-separated list of pooled pseudoreplicate peak files for IDR analysis.')
        
        self.add_argument('--pooled_peaks', nargs='?', dest='pooled_peaks',
                help='Homer peak file for pooled tag directories. This will '
                + 'be used to generate the final set of peaks.')
        
        self.add_argument('--rep_narrowpeaks', nargs='*', dest='rep_narrowpeaks',
                help='Space-separated list of already-processed replicate '
                + 'narrowPeak files to be input directly into IDR analysis.')
        self.add_argument('--pseudorep_narrowpeaks', nargs='*', dest='pseudorep_narrowpeaks',
                help='Space-separated list of already-processed pseudoreplicate '
                + 'narrowPeak files to be input directly into IDR analysis.')
        self.add_argument('--pooled_narrowpeaks', nargs='*', dest='pooled_narrowpeaks',
                help='Space-separated list of already-processed pooled pseudoreplicate '
                + 'narrowPeak files to be input directly into IDR analysis.')

        self.add_argument('--rep_idr_peaks', nargs='*', dest='rep_idr_peaks',
                help='Space-separated list of already-processed replicate '
                + 'IDR peaks to be input directly into threshold analysis.')
        self.add_argument('--pseudorep_idr_peaks', nargs='*', dest='pseudorep_idr_peaks',
                help='Space-separated list of already-processed pseudoreplicate '
                + 'IDR peaks to be input directly into threshold analysis.')
        self.add_argument('--pooled_idr_peaks', nargs='*', dest='pooled_idr_peaks',
                help='Space-separated list of already-processed pooled pseudoreplicate '
                + 'IDR peaks to be input directly into threshold analysis.')
        
        self.add_argument('--pooled_dir_name', nargs='?', dest='pooled_dir_name',
                help='Base name for pooled pseudorep directories.')
        
        self.add_argument('--pseudorep_count', nargs='?', dest='pseudorep_count',
                type=int, default=2,
                help='Number of pseudoreplicates to create. Default: 2')
        
        self.add_argument('--ranking_measure', nargs='?', dest='ranking_measure',
                choices=['tag-count', 'p-value'], default='tag-count',
                help='Use tag-count or p-value for comparing replicates? ' +
                'Default: tag-count')
        self.add_argument('--number_of_peaks', nargs='?', dest='number_of_peaks',
                type=int, 
                help='If you are passing in already-processed IDR peak files, '
                + 'this is the original number of input peaks per Homer peak file '
                + 'that should be used for automatically determining the threshold.')
        self.add_argument('--threshold', nargs='?', dest='threshold',
                type=float, 
                help='Specificically set the IDR threshold if you do not want '
                + 'it to be auto-calculated based on the number of peaks.')
        self.add_argument('--pooled_threshold', nargs='?', dest='pooled_threshold',
                type=float, 
                help='Specificically set the IDR threshold for pooled pseudoreps '
                + 'if you do not want '
                + 'it to be auto-calculated based on the number of peaks.')
        
        
        
    def homer2narrow(self, options, peak_files, output_dir=None):
        '''
        Convert passed Homer peak files to narrowPeak files as specified by 
        the IdrUtilities object.
        
        Returns the set of filenames for generated narrowPeak files.
        '''
        output_dir = output_dir or options.output_dir
        self.check_output_dir(output_dir)
             
        idrutils = IdrUtilities()
        output_files = []
        for peak_file in peak_files:
            # Get extensionless name of file
            basename = os.path.splitext(os.path.basename(peak_file))[0]
            # Add a randint to avoid name collision
            basename = basename + '_' + str(randint(1,999))
            output_file = os.path.join(output_dir, basename + '.narrowPeak')
            
            data = idrutils.import_homer_peaks(peak_file)
            idrutils.homer_to_narrow_peaks(data, output_file)
            
            print('NarrowPeak file output to {}'.format(output_file))
            output_files.append(output_file)
        return output_files
    
    def pseudoreplicate(self, options, suffix='Pseudorep'):
        '''
        Generate pseudoreplicates for passed tag directory by splitting randomly.
        
        Returns sets of pseudoreps such that each numbered rep is grouped together:
        [(Sample1-Pseudorep1, Sample2-Pseudorep1, Sample3-Pseudorep1),
        (Sample1-Pseudorep2, Sample2-Pseudorep2, Sample3-Pseudorep2)...]
        '''
        self.check_output_dir(options.output_dir)
        
        idrutils = IdrUtilities()
        pseudoreps = []
        for tag_dir in options.tag_dirs:
            print('Generating {} pseudoreplicate tag directories for {}'.format(
                                options.pseudorep_count, tag_dir))
            pseudoreps.append(idrutils.create_pseudoreps(tag_dir, 
                                        options.output_dir, 
                                        count=options.pseudorep_count,
                                        suffix=suffix))
        
        return list(zip(*pseudoreps))
            
    def pool_pseudoreplicates(self, options):
        '''
        Generate pseudoreplicates for each directory, then pool the pseudoreps.
        '''
        if not options.pooled_dir_name:
            raise Exception('A name for the pooled directory is needed. '
                            + 'Please indicate one with the --pooled-dir-name option.')
            
        pseudorep_sets = self.pseudoreplicate(options, suffix='Pooling-Pseudorep')
        
        idrutils = IdrUtilities()
        for i, pseudorep_set in enumerate(pseudorep_sets):
            idrutils.clean_up_pseudoreps(os.path.join(options.output_dir,
                                            options.pooled_dir_name + 
                                            '-Pseudorep' + str(i + 1)), 
                                     pseudorep_set)
        
    def truncate(self, options, peak_files, output_dir=None):
        '''
        Truncate SORTED narrowPeak files so that they are all the same length.
        '''
        self.check_output_dir(output_dir or options.output_dir)
        
        idrutils = IdrUtilities()
        output_files = idrutils.standardize_peak_counts(peak_files, 
                                                        output_dir)
        
        return output_files
        
        
    def idr(self, options):
        '''
        Go through entire IDR pipeline, starting from replicate peak_files
        and pseudoreplicate peak files.
        
        '''
        self.check_output_dir(options.output_dir)
        
        peak_sets = [options.peak_files, 
                     options.pseudorep_files, 
                     options.pooled_pseudoreps]
        narrowpeak_sets = [options.rep_narrowpeaks,
                           options.pseudorep_narrowpeaks,
                           options.pooled_narrowpeaks]
        idr_peak_sets = [options.rep_idr_peaks,
                           options.pseudorep_idr_peaks,
                           options.pooled_idr_peaks]
        all_peaks = [item for sublist in peak_sets for item in sublist]
        all_narrowpeaks = [item for sublist in narrowpeak_sets for item in sublist]
        all_idr_peaks = [item for sublist in idr_peak_sets for item in sublist]
        
        # Set up our parameters
        if options.ranking_measure == 'p-value':
            ranking_measure = 'p.value'
        else: ranking_measure = 'signal.value'
        
        
        replicate_dir = os.path.join(options.output_dir, 'replicate_comparisons')
        pseudorep_dir = os.path.join(options.output_dir, 'pseudorep_comparisons')
        pooled_dir = os.path.join(options.output_dir, 'pooled_comparisons')
        for d in (replicate_dir, pseudorep_dir, pooled_dir):
            if not os.path.exists(d): os.mkdir(d)
            
        if len(all_idr_peaks) == 0:
            if len(all_narrowpeaks) == 0:
                narrowpeak_dir = os.path.join(options.output_dir, 'narrowpeaks')
                # First, convert all peak files and truncate to the same length
                narrows = self.homer2narrow(options, peak_files=all_peaks,
                                            output_dir=narrowpeak_dir)
                truncated = self.truncate(options, peak_files=narrows,
                                          output_dir=narrowpeak_dir)
                
                # Split back out into two separate groups
                rep_truncated = truncated[:len(options.peak_files)]
                pseudorep_truncated = truncated[len(options.peak_files):\
                                                (len(options.peak_files) + \
                                                len(options.pseudorep_files))]
                pooled_truncated = truncated[-len(options.pooled_pseudoreps):]
            else:
                rep_truncated = options.rep_narrowpeaks
                pseudorep_truncated = options.pseudorep_narrowpeaks
                pooled_truncated = options.pooled_narrowpeaks
            
            # Compare our replicates, pairwise.
            idrcaller = IdrCaller()
            rep_prefixes = idrcaller.compare_replicates(rep_truncated, 
                                                replicate_dir, ranking_measure)
            pseudorep_prefixes = idrcaller.compare_pseudoreps(pseudorep_truncated, 
                                                pseudorep_dir, ranking_measure)
            pooled_prefixes = idrcaller.compare_pseudoreps(pooled_truncated, 
                                                pooled_dir, ranking_measure)
            
            # Where did we output our files?
            suffix = '-overlapped-peaks.txt'
            rep_files = []
            for prefix in rep_prefixes:
                prefix = os.path.basename(prefix)
                rep_files.append(os.path.join(replicate_dir, prefix + suffix))
            pseudorep_files = []
            for prefix in pseudorep_prefixes:
                prefix = os.path.basename(prefix)
                pseudorep_files.append(os.path.join(pseudorep_dir, prefix + suffix))
            pooled_files = []
            for prefix in pooled_prefixes:
                prefix = os.path.basename(prefix)
                pooled_files.append(os.path.join(pooled_dir, prefix + suffix))
        
            # Plot all of our pairwise comparisons
            plot_dir = os.path.join(options.output_dir, 'plots')
            if not os.path.exists(plot_dir): os.mkdir(plot_dir)
            idrcaller.plot_comparisons(rep_prefixes, plot_dir, 
                                       output_prefix='Replicate_comparison')
            
            idrcaller.plot_comparisons(pseudorep_prefixes, plot_dir, 
                                       output_prefix='Pseudorep_comparison')
            
            idrcaller.plot_comparisons(pooled_prefixes, plot_dir, 
                                       output_prefix='Pooled_pseudorep_comparison')
         
        else:
            rep_files = options.rep_idr_peaks
            pseudorep_files = options.pseudorep_idr_peaks
            pooled_files = options.pooled_idr_peaks
        
        # We need the number of peaks input into analysis
        # to automatically determine a threshold
        try:
            number_of_peaks = len(open(rep_truncated[0], 'r').readlines())
        except Exception:
            # We skipped truncating our peaks; expect number of peaks
            # OR a threshold.
            if options.number_of_peaks is None \
                and (options.threshold is None or options.pooled_threshold is None):
                    raise Exception('You must pass the number_of_peaks '
                                    + 'or both a threshold and pooled_threshold '
                                    + 'to complete analysis.')
            number_of_peaks = options.number_of_peaks
            
        threshold = self.get_threshold(options, number_of_peaks, pooled=False)
        pooled_threshold = self.get_threshold(options, 
                                    number_of_peaks, pooled=True)
        
        if options.pooled_peaks:
            self.slice_pooled_peaks(threshold, pooled_threshold,
                                    rep_files, pseudorep_files, pooled_files,
                                    options.pooled_peaks, options.output_dir,
                                    ranking_measure=options.ranking_measure)
    
    def get_threshold(self, options, number_of_peaks, pooled=False):
        idrutil = IdrUtilities()
            
        # Determine our threshold
        if not pooled and options.threshold:
            threshold = options.threshold
        elif pooled and options.pooled_threshold:
            threshold = options.pooled_threshold
        else:
            threshold = idrutil.determine_threshold(number_of_peaks, 
                                                    pooled=pooled)
        
        return threshold
    
    def slice_pooled_peaks(self, threshold, pooled_threshold,
                           rep_files, pseudorep_files, pooled_files,
                           pooled_peaks, output_dir, ranking_measure='tag-count'):
        idrutil = IdrUtilities()
        # Determine how many peaks we want to keep.
        keep_count = idrutil.get_peaks_within_threshold(threshold, 
                                                          rep_files)
        idrutil.get_peaks_within_threshold(threshold, pseudorep_files)
        pooled_count = idrutil.get_peaks_within_threshold(pooled_threshold, 
                                                          pooled_files)
        
        # Pooled count should be within 2-fold of keep_count
        if abs(math.log(keep_count/pooled_count, 2)) > 1:
            print('!! Warning: The number of peaks within the replicate '
                  + 'threshold is not within two-fold of the number of '
                  + 'peaks within the pooled threshold. This could indicate '
                  + 'inconsistencies in the datasets.\n'
                  + 'Replicate count: {}, Pooled count: {}'.format(keep_count, 
                                                                   pooled_count))
        
        # Slice our pooled peak file accordingly.
        output_file = idrutil.slice_peaks(pooled_peaks, keep_count, 
                                          ranking_measure, output_dir)
        print('{} peaks output to {}'.format(keep_count, output_file))
    
    def sanitize_inputs(self, options):
        if not options.peak_files: options.peak_files = []
        if not options.tag_dirs: options.tag_dirs = []
        if not options.pseudorep_files: options.pseudorep_files = []
        if not options.pooled_pseudoreps: options.pooled_pseudoreps = []
        if not options.rep_narrowpeaks: options.rep_narrowpeaks = []
        if not options.pseudorep_narrowpeaks: options.pseudorep_narrowpeaks = []
        if not options.pooled_narrowpeaks: options.pooled_narrowpeaks = []
        if not options.rep_idr_peaks: options.rep_idr_peaks = []
        if not options.pseudorep_idr_peaks: options.pseudorep_idr_peaks = []
        if not options.pooled_idr_peaks: options.pooled_idr_peaks = []
        
        if options.output_dir: 
            options.output_dir = os.path.normpath(options.output_dir)
        
        for i, f in enumerate(options.peak_files):
            options.peak_files[i] = os.path.normpath(f)
        for i, f in enumerate(options.tag_dirs):
            options.tag_dirs[i] = os.path.normpath(f)
        for i, f in enumerate(options.pseudorep_files):
            options.pseudorep_files[i] = os.path.normpath(f)
        for i, f in enumerate(options.rep_narrowpeaks):
            options.rep_narrowpeaks[i] = os.path.normpath(f)
        for i, f in enumerate(options.pseudorep_narrowpeaks):
            options.pseudorep_narrowpeaks[i] = os.path.normpath(f)
        for i, f in enumerate(options.rep_idr_peaks):
            options.rep_idr_peaks[i] = os.path.normpath(f)
        for i, f in enumerate(options.pseudorep_idr_peaks):
            options.pseudorep_idr_peaks[i] = os.path.normpath(f)
        for i, f in enumerate(options.pooled_idr_peaks):
            options.pooled_idr_peaks[i] = os.path.normpath(f)
        
        return options

    def check_output_dir(self, output_dir):
        if not output_dir:
            raise Exception('An output directory is needed. '
                            + 'Please indicate one with the -o option.')
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    
if __name__ == '__main__':
    parser = IdrArgumentParser()
    options = parser.parse_args()

    options = parser.sanitize_inputs(options)
    
    if options.command == 'idr':
        parser.idr(options)
    elif options.command == 'pseudoreplicate':
        parser.pseudoreplicate(options)
    elif options.command == 'pool-pseudoreplicates':
        parser.pool_pseudoreplicates(options)
    elif options.command == 'homer2narrow':
        parser.homer2narrow(options, peak_files=options.peak_files)
    elif options.command == 'truncate':
        parser.truncate(options, peak_files=options.peak_files)
        
    
    else:
        print('Command {} not recognized.'.format(options.command))
        parser.print_help()