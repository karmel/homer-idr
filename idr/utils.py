'''
Created on Aug 4, 2014

@author: karmel

'''
import os
import re
import numpy as np
from pandas.io.parsers import read_csv
from collections import OrderedDict
from pandas import DataFrame, Series
import subprocess


class IdrUtilities(object):
    '''
    Various utilities for converting files, processing data, etc. that are 
    necessary to run the IDR R code.
    '''
    
    ######################################################
    # Creating pseudo-replicates
    ######################################################
    def create_pseudoreps(self, tag_dir, output_dir, count=2, suffix='Pseudorep'):
        '''
        Randomly split a Homer tag directory into two parts with approximately
        equal number of reads.
        
        @todo: this is super slow. Using Python IO is slower... maybe shuf 
        would be faster? But many systems, OS X included, don't come with shuf.
        '''
        tag_dir_name = os.path.basename(tag_dir)
        # Make destination directories
        pseudo_tag_dirs = []
        for i in range(1,count + 1):
            pseudo_tag_dirs.append(os.path.join(output_dir, 
                                    tag_dir_name + '-{}{}'.format(suffix,i)))
            # Make tmp directory
            os.mkdir(pseudo_tag_dirs[i - 1] + '-tmp')
        
        for f in os.listdir(tag_dir):
            # If this is a chr[X].txt file
            if re.match('chr[A-Za-z0-9]+\.tags\.tsv$',f):
                chr_file = os.path.join(tag_dir, f)
                shuffled_file = os.path.join(pseudo_tag_dirs[0] + '-tmp', 
                                             f + '.tmp')
                # Shuffle the source chr file:
                subprocess.check_call("awk 'BEGIN{srand()}{print rand(),$0}' " +
                                    "{} | sort -n | cut -d ' ' -f2- > {}".format(
                                    chr_file, shuffled_file), shell=True)
                
                # How many lines do we have to split up?
                line_count = subprocess.check_output('wc -l {}'.format(
                                    shuffled_file), shell=True)
                line_count = int(line_count.split()[0])
                per_file = line_count//count
                for i in range(1, count+1):
                    # Output into file.
                    subprocess.check_call('head -n {} {} | tail -n {} > {}'.format(
                                i*per_file, shuffled_file, per_file, 
                                os.path.join(pseudo_tag_dirs[i - 1] + '-tmp', f)),
                                shell=True)
        
                subprocess.check_call('rm {}'.format(shuffled_file), shell=True)
        
        # Finally, use Homer to clean up
        for i in range(0, count):
            self.clean_up_pseudoreps(pseudo_tag_dirs[i], 
                                     [pseudo_tag_dirs[i] + '-tmp'])
    
        return pseudo_tag_dirs
    
    def clean_up_pseudoreps(self, target_dir, source_dirs): 
        '''
        Use Homer to re-make tag directory, as we want an accurate tagInfo file.
        Make sure to use the /usr/bin/env clause to tell Python to look
        in the environment's path for Homer.
        '''
        cmd = '/usr/bin/env makeTagDirectory {} -d {}'.format(
                            target_dir,
                            ' '.join(source_dirs))
        
        subprocess.check_call(cmd.split())
        for source in source_dirs:
            subprocess.check_call('rm -rf {}'.format(source), shell=True)
        
        return target_dir
                
                        
    ######################################################
    # Converting Homer Peaks to narrowPeak files
    ######################################################
    def import_homer_peaks(self, filename):
        '''
        Takes filename, returns dataframe of Homer peaks after stripping
        out comment lines at the top.
        '''
        # Find header row
        header_row = None
        f = open(filename,'r')
        for i,line in enumerate(f):
            if line[:7] == '#PeakID': 
                header_row = i
                break
        if header_row is None: 
            raise Exception('There is no header in this Homer peak file!') 
        data = read_csv(filename, sep='\t', header=i)
        return data
        
    def homer_to_narrow_peaks(self, data, output_file):
        '''
        Given a Homer peak dataframe, extract necessary columns and convert
        to a broadPeak file. From the IDR package description:
        
            NarrowPeak files are in BED6+4 format. It consists of 10 tab-delimited columns
    
            1.chrom     string     Name of the chromosome
            2.chromStart     int     The starting position of the feature in the chromosome. The first base in a chromosome is numbered 0.
            3.chromEnd     int     The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the   feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
            4.name     string     Name given to a region (preferably unique). Use '.' if no name is assigned
            5.score     int     Indicates how dark the peak will be displayed in the browser (1-1000). If '0', the DCC will assign this based on signal value.         Ideally average signalValue per base spread between 100-1000.
            6.strand     char     +/- to denote strand or orientation (whenever applicable). Use '.' if no orientation is assigned.
            7.signalValue     float     Measurement of overall (usually, average) enrichment for the region.
            8.pValue     float     Measurement of statistical signficance (-log10). Use -1 if no pValue is assigned.
            9.qValue     float     Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.
            10.peak     int     Point-source called for this peak; 0-based offset from chromStart. Use -1 if no point-source called.
        
        '''

        columns = OrderedDict((
            ('chrom', self.get_first_column(data, ['chr','chrom', 'chromosome'])),
            ('chromStart', self.get_first_column(data, ['chromStart','start'])),
            ('chromEnd', self.get_first_column(data, ['chromEnd','end'])),
            ('name', self.get_first_column(data, ['#PeakID','PeakID','ID','name'])),
            ('score', Series([0]*data.shape[0])), # Leave zero so that signalValue column is used
            ('strand', self.get_first_column(data, ['strand'])),       
            ('signalValue', self.get_first_column(data, ['Normalized Tag Count', 'findPeaks Score', 'score', 'Total Tags'])),
            ('pValue', -np.log10(self.get_first_column(data, ['p-value vs Control', 'p-value vs Local', 'p-value']))),
            ('qValue', Series([-1]*data.shape[0])), # Leave -1 as no individual FDR is called for each peak
            ('peak', Series([-1]*data.shape[0])), # Leave -1 as no point-source is called for each peak
            ))
        df = DataFrame(columns)
        df = df.sort(['signalValue','pValue'], ascending=False)
        df.to_csv(output_file, sep='\t', header=False, index=False)
        
    def get_first_column(self, data, names):
        '''
        Given a dataset and a list of strings, return the first column that 
        exists from the passed names.
        '''
        for name in names:
            try: return data[name]
            except KeyError: pass
        raise Exception(
            'None of the columns "{}" were found.').format(', '.join(names))
        
    ######################################################
    # Standardizing peak counts for narrowPeak files
    ######################################################
    def standardize_peak_counts(self, peak_files, output_dir, max_count=None):
        '''
        Make sure all passed narrowPeak files have the same number of rows.
        
        Assumes peak files are sorted in descending order!
        '''
        # Find the peak file with the fewest rows
        min_rows = None
        for peak_file in peak_files:
            count = subprocess.check_output('wc -l {}'.format(peak_file), shell=True)
            count = int(count.split()[0])
            if min_rows is None or count < min_rows:
                min_rows = int(count)
        
        if max_count and min_rows > max_count: min_rows = max_count 
        
        # Truncate all to the min count
        output_files = []
        for peak_file in peak_files:
            basename, ext = os.path.splitext(os.path.basename(peak_file))
            output_file = os.path.join(output_dir, basename + '-truncated' + ext)
            subprocess.check_call('head -n {} {} > {}'.format(
                                        min_rows, peak_file, output_file), 
                                  shell=True)
            output_files.append(output_file)
        
        return output_files
            