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
import os
from argparse import ArgumentParser
from idr.utils import IdrUtilities
from idr.idr_caller import IdrCaller

class IdrArgumentParser(ArgumentParser):
    def __init__(self):
        description = '''Functions for running Irreproducibility Discovery Rate
        (IDR) analysis on Homer peak files.'''
        
        super().__init__(description=description)
        
        self.add_argument('command', 
                help='Program to run; options are: idr, homer2narrow'),
        self.add_argument('-p', '--peak_files', nargs='*', dest='peak_files',
                help='Space-separated list of input Homer peak files.')
        self.add_argument('-d', '--tag_dirs', nargs='*', dest='tag_dirs',
                help='Space-separated list of input Homer tag directories.')
        self.add_argument('-pr', '--pseudorep_files', nargs='*', dest='pseudorep_files',
                help='Space-separated list of pseudoreplicate peak files for IDR analysis.')
        self.add_argument('--rep_narrowpeaks', nargs='*', dest='rep_narrowpeaks',
                help='Space-separated list of already-processed replicate '
                + 'narrowPeak files to be input directly into IDR analysis.')
        self.add_argument('--pseudorep_narrowpeaks', nargs='*', dest='pseudorep_narrowpeaks',
                help='Space-separated list of already-processed pseudoreplicate '
                + 'narrowPeak files to be input directly into IDR analysis.')
        self.add_argument('--pooled-dir-name', nargs='?', dest='pooled_dir_name',
                help='Base name for pooled pseudorep directories.')
        self.add_argument('--pseudorep_count', nargs='?', dest='pseudorep_count',
                type=int, default=2,
                help='Number of pseudoreplicates to create. Default: 2')
        self.add_argument('-o','--output_dir', nargs='?', dest='output_dir',
                help='Directory name in which output files will be placed. ' +
                'Will be created if it does not exist.')
        self.add_argument('--ranking-measure', nargs='?', dest='ranking_measure',
                choices=['tag-count', 'p-value'], default='tag-count',
                help='Use tag-count or p-value for comparing replicates? ' +
                'Default: tag-count')
        
        
        
    def homer2narrow(self, options, peak_files):
        '''
        Convert passed Homer peak files to narrowPeak files as specified by 
        the IdrUtilities object.
        
        Returns the set of filenames for generated narrowPeak files.
        '''
        self.check_output_dir(options)
             
        idrutils = IdrUtilities()
        output_files = []
        for peak_file in peak_files:
            # Get extensionless name of file
            basename = os.path.splitext(os.path.basename(peak_file))[0]
            output_file = os.path.join(options.output_dir, basename + '.narrowPeak')
            
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
        self.check_output_dir(options)
        
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
        
    def truncate(self, options, peak_files):
        '''
        Truncate SORTED narrowPeak files so that they are all the same length.
        '''
        self.check_output_dir(options)
        
        idrutils = IdrUtilities()
        output_files = idrutils.standardize_peak_counts(peak_files, 
                                                        options.output_dir)
        
        return output_files
        
        
    def idr(self, options):
        '''
        Go through entire IDR pipeline, starting from replicate peak_files
        and pseudoreplicate peak files.
        
        '''
        if len(options.rep_narrowpeaks) == 0 \
            and len(options.pseudorep_narrowpeaks) == 0:
                # First, convert all peak files and truncate to the same length
                all_peaks = options.peak_files + options.pseudorep_files
                
                narrows = self.homer2narrow(options, peak_files=all_peaks)
                truncated = self.truncate(options, peak_files=narrows)
                print('Files truncated to {} rows.'.format(len(truncated[0])))
                
                # Split back out into two separate groups
                rep_truncated = truncated[:len(options.peak_files)]
                pseudorep_truncated = truncated[-len(options.pseudorep_files):]
        else:
            rep_truncated = options.rep_narrowpeaks
            pseudorep_truncated = options.pseudorep_narrowpeaks
        # Set up our parameters
        if options.ranking_measure == 'p-value':
            ranking_measure = 'p.value'
        else: ranking_measure = 'signal.value'
        
        
        replicate_dir = os.path.join(options.output_dir, 'replicate_comparisons')
        pseudorep_dir = os.path.join(options.output_dir, 'pseudorep_comparisons')
        for d in (replicate_dir, pseudorep_dir):
            if not os.path.exists(d): os.mkdir(d)
            
        # Compare our replicates, pairwise.
        idrcaller = IdrCaller()
        rep_prefixes = idrcaller.compare_replicates(rep_truncated, 
                                            replicate_dir, ranking_measure)
        pseudorep_prefixes = idrcaller.compare_replicates(pseudorep_truncated, 
                                            pseudorep_dir, ranking_measure)
        
        # Plot all of our pairwise comparisons
        idrcaller.plot_comparisons(rep_prefixes, replicate_dir, 
                                   output_prefix='Replicate_comparison')
        
        idrcaller.plot_comparisons(pseudorep_prefixes, replicate_dir, 
                                   output_prefix='Pseudorep_comparison')
        
    def sanitize_inputs(self, options):
        if not options.peak_files: options.peak_files = []
        if not options.tag_dirs: options.tag_dirs = []
        if not options.pseudorep_files: options.pseudorep_files = []
        if not options.rep_narrowpeaks: options.rep_narrowpeaks = []
        if not options.pseudorep_narrowpeaks: options.pseudorep_narrowpeaks = []
        
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
        
        return options

    def check_output_dir(self, options):
        if not options.output_dir:
            raise Exception('An output directory is needed. '
                            + 'Please indicate one with the -o option.')
        if not os.path.exists(options.output_dir):
            os.makedirs(options.output_dir)
    
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