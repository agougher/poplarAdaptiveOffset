Poplar adaptive offset
================

Data and code for: Gougherty, A.V., Keller, S.R. and M.C. Fitzpatrick.
Future climate change promotes novel gene-climate associations in balsam
poplar (Populus balsamifera L.), a forest tree species. Link.

### Data

In the data folder are two files, a pairwise \(F_{ST}\) table, and a
table of minor allele frequencies (MAF), both with data from [Keller et
al. 2017](https://academic.oup.com/jhered/article/109/1/47/4605251), for
a subset of SNPs in the flowering time gene network in balsam poplar.
The \(F_{ST}\) table is used to fit a generalized dissimilarity model,
and the MAF table is used in a Gradient Forest model.

### Code

I have broken up the code into two scripts, one each for GDM and GF, as
it may not always be desired to calculate both.
