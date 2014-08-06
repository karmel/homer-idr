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

class IdrArgumentParser(ArgumentParser):
    def __init__(self):
        description = '''Functions for running Irreproducibility Discovery Rate
        (IDR) analysis on Homer peak files.'''
        
        super().__init__(description=description)
        
        self.add_argument('command', 
                help='Program to run; options are: idr, homer2narrow'),
        self.add_argument('-p','--peak_files', nargs='*', dest='peak_files',
                help='Space-separated list of input Homer peak files.')
        self.add_argument('-o','--output_dir', nargs='?', dest='output_dir',
                help='Directory name in which output files will be placed. ' +
                'Will be created if it does not exist.')
        self.add_argument('--comparator', nargs='?', dest='comparator',
                choices=['tag-count', 'p-value'], default='tag-count',
                help='Use tag-count or p-value for comparing replicates? ' +
                'Default: tag-count')
        
        
        
    def homer2narrow(self, options):
        '''
        Convert passed Homer peak files to narrowPeak files as specified by 
        the IdrUtilities object.
        
        Returns the set of filenames for generated narrowPeak files.
        '''
        if not options.output_dir:
            raise Exception('An output directory is needed to run homer2narrow.')
        if not os.path.exists(options.output_dir):
            os.makedirs(options.output_dir)
             
        idrutils = IdrUtilities()
        output_files = []
        for peak_file in options.peak_files:
            # Get extensionless name of file
            basename = os.path.splitext(os.path.basename(peak_file))[0]
            output_file = os.path.join(options.output_dir, basename + '.narrow.bed')
            
            data = idrutils.import_homer_peaks(peak_file)
            idrutils.homer_to_narrow_peaks(data, output_file)
            
            print('NarrowPeak file output to {}'.format(output_file))
            output_files.append(output_file)
        return output_files
    
if __name__ == '__main__':
    parser = IdrArgumentParser()
    options = parser.parse_args()

    if options.command == 'idr':
        narrows = parser.homer2narrow(options)
        for narrow_peak in narrows:
            pass
    elif options.command == 'homer2narrow':
        parser.homer2narrow(options)
    else:
        print('Command {} not recognized.'.format(options.command))
        parser.print_help()