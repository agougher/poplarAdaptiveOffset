Poplar adaptive offset
================

Data and code for: Gougherty, A.V., Keller, S.R. and M.C. Fitzpatrick. 2021. 
Maladaptation, migration and extirpation fuel climate change risk in a forest 
tree species. Nature Climate Change 11:166–171 [Link.](https://www.nature.com/articles/s41558-020-00968-6)

### Data

In the data folder are several files:

  - popInfo.csv: names and locations (longitudes/latitudes) of 81 balsam
    poplar populations
  - fstTab.csv: pairwise *F<sub>ST</sub>* table for 81 populations
  - mafTab.csv: minor allele frequency (MAF) table for 75 SNPs for 81
    populations

Genetic data are from [Keller et
al. 2017](https://academic.oup.com/jhered/article/109/1/47/4605251),
for a subset of SNPs in the flowering time gene network in balsam
poplar. The *F<sub>ST</sub>* table is used to fit a generalized
dissimilarity model, and the MAF table is used in a Gradient Forest
model.

Climate data and the balsam poplar range file are available [here](https://drive.google.com/drive/folders/1NG-OPEmY9dmEWK_YoMCUF7YuFXJfOtvd?usp=sharing).

### Code

The code is broken up into two scripts, one each for GDM and GF, as it
may not always be desired to calculate both. The code was written to be
(relatively) easily modified for other applications. To recreate the
analyses in the article, download climate data from [WorldClim
v 1.4](https://www.worldclim.org/data/v1.4/worldclim14.html) and
Little’s range map for balsam poplar, now available
[here](https://github.com/wpetry/USTreeAtlas/).
