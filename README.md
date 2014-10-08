homer-idr
=========

A small package for applying Irreproducibility Discovery Rate (IDR) analysis for replicate chip-seq experiments analyzed using Homer.

Questions? Comments? Email me: <karmel@arcaio.com>

Last updated 2014-08.

Cite: 
[![DOI](https://zenodo.org/badge/5954/karmel/homer-idr.png)](http://dx.doi.org/10.5281/zenodo.11619)

## TOC

1. [Introduction](#introduction)
2. [Open Questions](#open-questions)
3. [Installation](#installation)
4. [Using homer-idr](#using-homer-idr)
5. [Other Usage Notes](#other-usage-notes)

## Introduction

The Irreproducibility Discovery Rate (IDR) statistic has been adopted by Encode in order to incorporate and interpret replicates in chip-sequencing experiments. [A procedure and an R package][IDR] have been developed to calculate the IDR statistic and call peaks accordingly; I highly suggest you read through the documentation there for a full understanding of what we are doing here.

The canonical IDR pipeline calls peaks with SPP or MACS. We here present some methods that (A) allow for the use of Homer peaks, and (B) make some of the initial data prep methods easier.

**This has not yet been extensively tested, and many important questions remain (see [below](#open-questions)).** Hopefully, this is enough to get you started with IDR analysis, and we can answer these questions together.

## Open Questions

There are many open questions about the best way to incorporate multiple replicates with IDR. I will highlight two below that seem especially pressing.

#### 1. How should we be calling peaks? 

The peaks that are input to IDR analysis are key to the output that is generated, so it is crucial to have high-quality peaks going in. The instructions for IDR indicate that peaks should be called permissively, such that a great proportion of the input peaks are just noise. The authors of the R package we are using suggest calling 150,000 to 300,000 peaks using SPP or about 100,000 using MACS.

Thus far, I have not found a great way to do this with Homer, especially with histone mark peaks. I have played with a number of parameters, and found a set that gives me more peaks, but not 150,000, and often substantially fewer than that, at least with histone marks. Further, relaxing thresholds does not just change the number of peaks, but the nature of the peaks, as more noise gets stitched into regions for the larger peaks. Finally, each replicate should have the same number of peaks called, and it does not seem possible with Homer to specify the returned number of peaks.

All that said, using the following parameters in initial tests seems to at least give me a wealth of peaks that are not too different on an individual basis from the peaks that are called with more stringent criteria:

	-P .1 -LP .1 -poisson .1

That is, p-value over input of up to .1, p-value over local tag count of up to .1, and overall Poisson p-value of up to .1. Note that because Homer has a number of different filtering mechanisms, it is not enough to just change the overall Poisson p-value or FDR.

Suggestions on how to best call peaks to all for lots of noise but not disturb the integrity of individual peaks much appreciated!

#### 2. How should we set the IDR threshold?

The IDR statistic is called algorithmically over a replicate set, but where to draw the line that separates noise from real peaks is determined by the user. According to the authors of the IDR package, the following guidelines apply:

- "If you started with ~150 to 300K relaxed pre-IDR peaks for 
large genomes (human/mouse), then threshold of 0.01 or 0.02 
generally works well."

- "If you started with < 100K pre-IDR peaks for large genomes 
(human/mouse), then threshold of 0.05 is more appropriate."

- "If you started with ~150 to 300K relaxed pre-IDR peaks for large 
genomes (human/mouse), then threshold of 0.0025 or 0.005 generally 
works well. We use a tighter threshold for pooled-consistency 
since pooling and subsampling equalizes the pseudoreplicates in 
terms of data quality. So we err on the side of caution and 
use more stringent thresholds. The equivalence between a 
pooled-consistency threshold of 0.0025 and original replicate 
consistency threshold of 0.01 was calibrated based on a 
gold-standard pair of high quality replicate datasets for the CTCF 
transcription factor in human."

We have interpreted this such that we establish a linear relationship between number of input peaks and the threshold, such that if you start with 75,000 peaks per sample, the threshold is .05, if you start with 300,000 peaks per sample, the threshold is .01, and everything between is scaled linearly. (And a comparable setup is employed for the pooled pseudoreplicates.) 

It is also possible to pass in explicit thresholds with the `--threshold` and `--pooled_threshold` options. In any case, make sure you review your returned peaks and consider the appropriateness of the selected threshold.

## Installation

Prerequisites:
- Homer, installed and available in your executable path
- Python 3, numpy, and pandas. We highly recommend installing the [Anaconda distribution](https://store.continuum.io/cshop/anaconda/) of Python 3.

To install homer-idr, clone the github repository into the directory from which you want to run it.

For example:

	cd ~/software
	mkdir homer-idr
	cd homer-idr
	git clone https://github.com/karmel/homer-idr.git

Then, add the Python package to your `PYTHONPATH`, either in your user's shell profile:

	nano ~/.bash_profile

And add the lines:

	PYTHONPATH=$PYTHONPATH:/home/me/software/homer-idr/homer-idr
	export PYTHONPATH

After saving your updated bash_profile, make sure to reload it with the command:

	source ~/.bash_profile

Alternatively, for less frequent use, you can just set the `PYTHONPATH` immediately before running, at the command line:

	PYTHONPATH=$PYTHONPATH:/home/me/software/homer-idr/homer-idr

You should now be able to run the run_idr.py script:

	python ~/software/homer-idr/homer-idr/idr/run_idr.py --help


## Using homer-idr

The procedure described here is a recapitulation of the procedure defined in the [IDR package][IDR] (as of 2014-08), but using Homer and homer-idr. For full explanations and more detail, please refer to the [IDR package][IDR] documentation.

#### 1. Create a combined input tag directory.

If you have separate inputs for each of your replicates, combine them all into one big input directory. For example:

	mkdir ~/CD4TCell-IDR
	cd ~/CD4TCell-IDR
	makeTagDirectory CD4TCell-Input-Combined -d /data/CD4TCell-Input-1 /data/CD4TCell-Input-2 /data/CD4TCell-Input-3

#### 2. Create a pooled replicate tag directory.

As with the inputs, you want one tag directory that combines all the separate replicate tag directories. For example:

	cd ~/CD4TCell-IDR
	makeTagDirectory CD4TCell-H3K4me2-Combined -d /data/CD4TCell-H3K4me2-1 /data/CD4TCell-H3K4me2-2 /data/CD4TCell-H3K4me2-3


#### 3. Permissively call peaks for each of your replicates using the combined input.

As discussed [above](#1-how-should-we-be-calling-peaks), IDR requires that we first call peaks permissively, such that a large proportion are actually noise. Read [above](#1-how-should-we-be-calling-peaks) for more detail on this, but I have been using the following Homer parameters to call peaks permissively:

	-P .1 -LP .1 -poisson .1

For example, to call peaks for a H3K4me2 histone chip-seq data set, I run:

	cd ~/CD4TCell-IDR
	mkdir -p peaks/replicates
	cd peaks/replicates
	findPeaks /data/CD4TCell-H3K4me2-1 -P .1 -LP .1 -poisson .1 -style histone -nfr -i ~/CD4TCell-IDR/CD4TCell-Input-Combined -o CD4TCell-H3K4me2-1_peaks.txt

Make sure to use the same parameters for all replicates.

#### 4. Call peaks on the pooled tag directory.

Make sure to use the **same permissive parameters** from above. For example:


	cd ~/CD4TCell-IDR
	mkdir -p peaks/pooled
	cd peaks/pooled
	findPeaks ~/CD4TCell-IDR/CD4TCell-H3K4me2-Combined -P .1 -LP .1 -poisson .1 -style histone -nfr -i ~/CD4TCell-IDR/CD4TCell-Input-Combined -o CD4TCell-H3K4me2-Combined_peaks.txt

Make sure to use the same parameters for all replicates.


#### 5. Create pseudoreplicates for each replicate directory.

We want to split the tags in each replicate randomly, so that we can analyze each replicate for internal consistency. For example:

	cd ~/CD4TCell-IDR
	mkdir -p pseudoreps/individual
	# python run_idr.py pseudoreplicate -d [tag_dirs to split] -o [output_file]
	python ~/software/homer-idr/homer-idr/idr/run_idr.py pseudoreplicate -d /data/CD4TCell-H3K4me2-1 /data/CD4TCell-H3K4me2-2 /data/CD4TCell-H3K4me2-3 -o pseudoreps/individual

Note: This process takes longer than expected. I am using `awk` and other command line tools to shuffle the tag order then split, so I feel like it should be relatively fast, but it's not. Let me know if you know of a better cross-platform way to complete this task.

#### 6. Create pseudoreplicates for the pooled directory.

We repeat the pseudoreplication process for our pooled tag directory. For example:

	cd ~/CD4TCell-IDR
	mkdir -p pseudoreps/pooled
	# python run_idr.py pseudoreplicate -d [tag_dirs to split] -o [output_file]
	python ~/software/homer-idr/homer-idr/idr/run_idr.py pseudoreplicate -d ~/CD4TCell-H3K4me2-Combined -o pseudoreps/individual

#### 7. Call peaks on each of the individual pseudoreplicate tag directories.

Make sure to use the **same permissive parameters** from above. For example:

	cd ~/CD4TCell-IDR
	mkdir -p peaks/pseudoreps
	cd peaks/pseudoreps 
	for f in ~/CD4TCell-IDR/pseudoreps/individual/*
		do
		findPeaks $f -P .1 -LP .1 -poisson .1 -style histone -nfr -i ~/CD4TCell-IDR/CD4TCell-Input-Combined -o ${f}_peaks.txt
		done

#### 8. Call peaks on each of the pooled pseudoreplicate tag directories.

Again, use the **same permissive parameters** from above. 

	cd ~/CD4TCell-IDR
	mkdir -p peaks/pooled-pseudoreps
	cd peaks/pooled-pseudoreps
	for f in ~/CD4TCell-IDR/pseudoreps/pooled-pseudoreps/*
		do
		findPeaks $f -P .1 -LP .1 -poisson .1 -style histone -nfr -i ~/CD4TCell-IDR/CD4TCell-Input-Combined -o ${f}_peaks.txt
		done

#### 9. Run IDR analysis.

Now that we have all of our peak files prepped, we can run the IDR analysis.

To sum the steps above, we need:

- Peaks from each of our replicates.
- Peaks from the pooled replicates.
- Peaks from each of the individual pseudoreplicates.
- Peaks from the pooled pseudoreplicates.

These sets of peaks then get fed into the run_idr.py program:

	# python run_idr.py idr -p [replicate peaks] \
	#	-pr [pseudorep peaks] \
	#	-ppr [pooled pseudorep peaks] \
	#	--pooled_peaks [pooled replicate peak file]
	#	-o [output directory to create]

Continuing our example, then:

	python ~/software/homer-idr/homer-idr/idr/run_idr.py idr -p ~/CD4TCell-IDR/peaks/replicates/* -pr ~/CD4TCell-IDR/peaks/pseudoreps/* -ppr ~/CD4TCell-IDR/peaks/pooled-pseudoreps --pooled_peaks ~/CD4TCell-IDR/peaks/pooled/CD4TCell-H3K4me2-Combined_peaks.txt -o ~/CD4TCell-IDR/idr-output

The same command, separated out into multiple lines for readability:

	python ~/software/homer-idr/homer-idr/idr/run_idr.py idr \
	-p ~/CD4TCell-IDR/peaks/replicates/* \
	-pr ~/CD4TCell-IDR/peaks/pseudoreps/* \
	-ppr ~/CD4TCell-IDR/peaks/pooled-pseudoreps \
	--pooled_peaks ~/CD4TCell-IDR/peaks/pooled/CD4TCell-H3K4me2-Combined_peaks.txt \
	-o ~/CD4TCell-IDR/idr-output

The homer-idr program will take the input peaks and do the following:

1. Convert the Homer peak files into the Encode narrowPeak format.
2. Truncate the narrowPeak files such that each replicate, pseudoreplicate, and pooled pseudoreplicate set has the same number of peaks.
3. Use the [IDR R package][IDR] to calculate the IDR statistic on the truncated narrowPeak files.
4. Calculate the IDR threshold to use for filtering (see [above](#2-how-should-we-set-the-idr-threshold) for more detail on how the threshold is selected).
5. Determine how many peaks should be kept in the final set based on the IDR threshold.
6. Sort the Homer peaks from the pooled replicate set and output a final Homer peak file with the chosen number of peaks.

In the specified output directory, you will find:

- The narrowPeak and truncated narrowPeak files that were created.
- The output from the [IDR R package][IDR], which includes
	- An -overlapped-peaks.txt file for each peak file, which lists peaks and their IDR statistics,
	- An -aboveIDR.txt file for each peak file, which lists how many peaks pass given IDR thresholds, and
	- A plots directory that contains generate plots comparing the replicates and pseudoreplicates.
- A **final Homer peak file**, named like the input `--pooled_peaks` file but suffixed with -top-set.txt, that has the peaks from the pooled replicate peak set cut off with only the top selected peaks. These are the peaks considered confident and likely real based on the IDR analysis.

The final Homer peak file can then be used for subsequent analysis with Homer or other programs to your replicated heart's content!


## Other Usage Notes

- When in doubt, ask for help:

		python ~/software/homer-idr/homer-idr/idr/run_idr.py --help

- If you want to try using different threshold values, you do not need to re-run the whole process. Just set the threshold manually, and pass in the already-processed IDR -overlapped-peaks.txt files like so:

		python ~/software/homer-idr/homer-idr/idr/run_idr.py idr \
		--rep_idr_peaks ~/CD4TCell-IDR/idr-output/replicate_comparisons/*overlapped-peaks.txt \
		--pseudorep_idr_peaks ~/CD4TCell-IDR/idr-output/pseudorep_comparisons/*overlapped-peaks.txt \
		--pooled_idr_peaks ~/CD4TCell-IDR/idr-output/pooled_comparisons/*overlapped-peaks.txt \
		--pooled_peaks ~/CD4TCell-IDR/peaks/pooled/CD4TCell-H3K4me2-Combined_peaks.txt \
		--threshold .04 --pooled_threshold .003 \
		-o ~/CD4TCell-IDR/idr-output-2

- By defauly, homer-idr uses the Normalized Tag Count to sort and compare peaks. If you would like to use p-value instead, use the parameter `--ranking_measure` when running the `idr` command:

		--ranking_measure p-value

- homer-idr can also be used to convert Homer peak files to narrowPeak files, or to truncate narrowPeak files to the same lenght:

		python ~/software/homer-idr/homer-idr/idr/run_idr.py homer2narrow -p ~/CD4TCell-Ets1_homer_peaks.txt -o ~/narrowPeak_files

		python ~/software/homer-idr/homer-idr/idr/run_idr.py truncate -p ~/narrowPeak_files/* -o ~/truncated_files
		
[IDR]: https://sites.google.com/site/anshulkundaje/projects/idr
