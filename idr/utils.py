'''
Created on Aug 4, 2014

@author: karmel

'''
import numpy as np
from pandas.io.parsers import read_csv
from collections import OrderedDict
from pandas import DataFrame, Series

class IdrUtilities(object):
    '''
    Various utilities for converting files, processing data, etc. that are 
    necessary to run the IDR R code.
    '''
    
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
        
        