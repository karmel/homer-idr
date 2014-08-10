homer-idr
=========

A small package for applying Irreproducibility Discovery Rate (IDR) analysis for replicate chip-seq experiments analyzed using Homer.

Questions? Comments? Email me: <karmel@arcaio.com>

## TOC

1. [Introduction](#introduction)
2. [Open Questions](#open-questions)
3. [Installation](#installation)
4. [Using homer-idr](#using-homer-idr)

## Introduction

The Irreproducibility Discovery Rate (IDR) statistic has been adopted by Encode in order to incorporate and interpret replicates in chip-sequencing experiments. [A procedure and an R package](https://sites.google.com/site/anshulkundaje/projects/idr) have been developed to calculate the IDR statistic and call peaks accordingly; I highly suggest you read through the documentation there for a full understanding of what we are doing here.

The canonical IDR pipeline calls peaks with SPP or MACS. We here present some methods that (A) allow for the use of Homer peaks, and (B) make some of the initial data prep methods easier.

**This has not yet been extensively tested, and many important questions remain (see [below](#open-questions)).** Hopefully, this is enough to get you started with IDR analysis, and we can answer these questions together.

## Open Questions

There are many open questions about the best way to incorporate multiple replicates with IDR. I will highlight two below that seem especially pressing.

### 1. How should we be calling peaks? 

The peaks that are input to IDR analysis are key to the output that is generated, so it is crucial to have high-quality peaks going in. The instructions for IDR indicate that peaks should be called permissively, such that a great proportion of the input peaks are just noise. The authors of the R package we are using suggest calling 150,000 to 300,000 peaks using SPP or about 100,000 using MACS.

Thus far, I have not found a great way to do this with Homer, especially with histone mark peaks. I have played with a number of parameters, and found a set that gives me more peaks, but not 150,000, and often substantially fewer than that, at least with histone marks. Further, relaxing thresholds does not just change the number of peaks, but the nature of the peaks, as more noise gets stitched into regions for the larger peaks. Finally, each replicate should have the same number of peaks called, and it does not seem possible with Homer to specify the returned number of peaks.

All that said, using the following parameters in initial tests seems to at least give me a wealth of peaks that are not too different on an individual basis from the peaks that are called with more stringent criteria:

	-P .1 -LP .1 -poisson .1

That is, p-value over input of up to .1, p-value over local tag count of up to .1, and overall Poisson p-value of up to .1. Note that because Homer has a number of different filtering mechanisms, it is not enough to just change the overall Poisson p-value or FDR.

Suggestions on how to best call peaks to all for lots of noise but not disturb the integrity of individual peaks much appreciated!

### 2. How should we set the IDR threshold?

The IDR statistic is called algorithmically over a replicate set, but where to draw the line that separates noise from real peaks is determined by the user. According to the authors of the IDR package, the following guidelines apply:

- "If you started with ~150 to 300K relaxed pre-IDR peaks for 
large genomes (human/mouse), then threshold of 0.01 or 0.02 
generally works well."

- "If you started with < 100K pre-IDR peaks for large genomes 
(human/mouse), then threshold of 0.05 is more appropriate."

- "If you started with ~150 to 300K relaxed pre-IDR peaks for large 
genomes (human/mouse), then threshold of 0.0025 or 0.005 generally 
works well. We use a tighter threshold for pooled-consistency 
since pooling and subsampling equalizes the pseudo-replicates in 
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

	cd /home/me/software
	mkdir homer-idr
	cd homer-idr
	git clone https://github.com/karmel/homer-idr.git

Then, add the Python package to your `PYTHONPATH`, either in your user's shell profile:

	nano ~/.bash_profile

And add the lines:

	PYTHONPATH=$PYTHONPATH:/home/me/software/homer-idr/homer-idr
	export PYTHONPATH

Or, for less frequent use, you can just set the `PYTHONPATH` immediately before running, at the command line:

	PYTHONPATH=$PYTHONPATH:/home/me/software/homer-idr/homer-idr

You should now be able to run the run_idr.py script:

	python /home/me/software/homer-idr/homer-idr/run_idr.py --help


## Using homer-idr



