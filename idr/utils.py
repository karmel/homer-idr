'''
Created on Aug 4, 2014

@author: karmel

'''
from collections import OrderedDict
import os
import re
import subprocess

from pandas import DataFrame, Series
from pandas.io.parsers import read_csv

import numpy as np
class IdrUtilities(object):
    '''
    Various utilities for converting files, processing data, etc. that are 
    necessary to run the IDR R code.
    '''
    
    p_value_columns = ['p-value vs Control', 'p-value vs Local', 'p-value']
    tag_count_columns = ['Normalized Tag Count', 'findPeaks Score', 
                         'score', 'Total Tags']
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

        # We don't want to require p-value, as Homer doesn't always output it.
        # Prep it here if it exists, or substitute tag count.
        pval_col = self.get_first_column(data,
            self.p_value_columns, required=False)
        if pval_col is not None:
            pvals = -np.log10(pval_col)
        else: 
            pvals = pvals = [-1]*data.shape[0]
            
        chrstart = Series(self.get_first_column(data, ['chromStart','start']))
        chrstart = chrstart.subtract(1) # use 0-based coordinate for narrowPeak (i.e. BED 6+4)
        
        columns = OrderedDict((
            ('chrom', self.get_first_column(data, ['chr','chrom', 'chromosome'])),
            ('chromStart', chrstart),
            ('chromEnd', self.get_first_column(data, ['chromEnd','end'])),
            ('name', self.get_first_column(data, ['#PeakID','PeakID','ID','name'])),
            ('score', Series([0]*data.shape[0])), # Leave zero so that signalValue column is used
            ('strand', self.get_first_column(data, ['strand'])),       
            ('signalValue', self.get_first_column(data, self.tag_count_columns)),
            ('pValue', pvals), # P-value if it exists, or tag count
            ('qValue', Series([-1]*data.shape[0])), # Leave -1 as no individual FDR is called for each peak
            ('peak', Series([-1]*data.shape[0])), # Leave -1 as no point-source is called for each peak
            ))
        df = DataFrame(columns)
        df = df.sort_values(['chrom', 'chromStart', 'chromEnd'], ascending=True) # sort by chromosome, start and end coordinates
        df['pValue'] = df['pValue'].round(4) # round pValue to 4 decimal places
        df.to_csv(output_file, sep='\t', header=False, index=False)
        
    def get_first_column(self, data, names, required=True):
        '''
        Given a dataset and a list of strings, return the first column that
        exists from the passed names.
        '''
        for name in names:
            try: return data[name]
            except KeyError: pass
        if required:
            raise Exception('None of the columns "{}" were found.'.format(
                ', '.join(names)))
        else:
            return None
        
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


    ######################################################
    # Post-processing
    ######################################################
    def determine_threshold(self, number_of_peaks, pooled=False):
        '''
        From the IDR documentation:
        
            - If you started with ~150 to 300K relaxed pre-IDR peaks for 
            large genomes (human/mouse), then threshold of 0.01 or 0.02 
            generally works well. 
            - If you started with < 100K pre-IDR peaks for large genomes 
            (human/mouse), then threshold of 0.05 is more appropriate.
            
            - If you started with ~150 to 300K relaxed pre-IDR peaks for large 
            genomes (human/mouse), then threshold of 0.0025 or 0.005 generally 
            works well. We use a tighter threshold for pooled-consistency 
            since pooling and subsampling equalizes the pseudo-replicates in 
            terms of data quality. So we err on the side of caution and 
            use more stringent thresholds. The equivalence between a 
            pooled-consistency threshold of 0.0025 and original replicate 
            consistency threshold of 0.01 was calibrated based on a 
            gold-standard pair of high quality replicate datasets for the CTCF 
            transcription factor in human. 
        
        So, we set 75K peaks and below at .05, then progress linearly to .01
        at 300K peaks.
        Similarly for pooled, but there we move from .0125 to .0025.
        '''
        print('Determining threshold based on {} peaks.'.format(number_of_peaks))
        
        if pooled: few_peaks, many_peaks = .0125, .0025
        else: few_peaks, many_peaks = .05, .01
        
        # Set up our linear equation, y = mx+b,
        # where x is the number of peaks and y is the desired threshold.
        m = (few_peaks - many_peaks)/(75000 - 300000)
        b = .01 - m*300000
        threshold = m*number_of_peaks + b
        print('Threshold: ' + str(threshold))
        
        return threshold
    
    def get_peaks_within_threshold(self, threshold, idr_files):
        '''
        For the generated IDR files, determine the greatest number
        that are within the given threshold.
        
        Meanwhile, compare counts across all the files. If any are
        substantially different (>2-fold), issue a warning.
        '''
        counts = []
        for filename in idr_files:
            data = read_csv(filename, sep=" ", header=0)
            within_thresh = len(data[data['IDR'] <= threshold])
            counts.append(within_thresh)
            
        if min(counts)*2 < max(counts):
            print('!! Warning: There is a large discrepancy between the number'
                  + ' of peaks within the threshold for each pair of compared'
                  + ' peak files. Please check individual output files and'
                  + ' determine whether one replicate is not like the others.')
            
        return max(counts)
    
    def slice_peaks(self, peak_file, number_of_peaks, 
                    ranking_measure, output_dir):
        '''
        Given a Homer tag file, import the file and output just the 
        specified number of peaks.
        '''
        data = self.import_homer_peaks(peak_file)
        
        sort_col = None
        if ranking_measure == 'p-value':
            ascending = True
            for col_name in self.p_value_columns: 
                if col_name in data.columns:
                    sort_col = col_name
                    break
        else:
            ascending = False
            for col_name in self.tag_count_columns: 
                if col_name in data.columns:
                    sort_col = col_name
                    break
        if not sort_col:
            raise Exception('Could not find column to sort final peaks by!')
            
        data = data.sort([sort_col], ascending=ascending)
        data = data[:number_of_peaks]
        
        # Output to file
        # Use the \r line ending because that is what Homer expects.
        basename, ext = os.path.splitext(os.path.basename(peak_file))
        output_file = os.path.join(output_dir, basename + '-top-set' + ext)
        data.to_csv(output_file, sep='\t', header=True, index=False, 
                    line_terminator='\r\n')
        return output_file