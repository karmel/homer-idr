homer-idr
=========

A small package for applying Irreproducibility Discovery Rate (IDR) analysis for replicate chip-seq experiments analyzed using Homer.

Questions? Comments? Email me: <karmel@arcaio.com>

# TOC

1. [Introduction](#introduction)
2. [Open Questions](#open-questions)
3. [Installation](#installation)

# Introduction

The Irreproducibility Discovery Rate (IDR) statistic has been adopted by Encode in order to incorporate and interpret replicates in chip-sequencing experiments. [A procedure and an R package](https://sites.google.com/site/anshulkundaje/projects/idr) have been developed to calculate the IDR statistic and call peaks accordingly; I highly suggest you read through the documentation there for a full understanding of what we are doing here.

The canonical IDR pipeline calls peaks with SPP or MACS. We here present some methods that (A) allow for the use of Homer peaks, and (B) make some of the initial data prep methods easier.

**This has not yet been extensively tested, and many important questions remain (see [below](#open-questions)).** Hopefully, this is enough to get you started with IDR analysis, and we can answer these questions together.

# Open Questions
