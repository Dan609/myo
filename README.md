# Danila E. Bobkov et al., Replicative senescence in MSCWJ-1 human umbilical cord mesenchymal stem cells is marked by characteristic changes in motility, cytoskeletal organization, and RhoA localization. Molecular Biology Reports, volume 47, pages 3867–3883 (2020). https://doi.org/10.1007/s11033-020-05476-6

Data and scripts used to generate the statistical analyses of protein colocalization, cell shape, F-actin fractality, and cell motility data.

## R scripts:

**CI.R** The proportion of cells with pronounced activity of *beta*-galactosidase during MSCWJ-1 cell line cultivation.

**coloc.R** Statistical analysis of colocalozation coefficients.

**traj_analysis_win.R** Cell movement tracks analysis.

**traj_analysis_mac.R** Mac version differs in path selection function.

**lmfrac.R** Actin cytoskeleton fractality analysis.

**shape.R** Cell shape analysis.


## Ipython notebook:

**significance_plot.ipynb** Multiple pairwise comparison of Myosin-9/F-actin colocalization data from long-term cultivation of WJMSC-1 cells.
Probe - passage number.
Rval,tM1,tM2,bTau,Rs - colocalization coefficients calculated in ImageJ Coloc2 plugin.


## Data:

**alltracks.csv** Tracks collected from intravital confocal microscopy using CQ1 cytometer.

**colocM9.csv** Myosin-9/F-actin colocalization data from long-term cultivation of WJMSC-1 cells. Probe - passage number.  Rval,tM1,tM2,bTau,Rs - colocalization coefficients calculated in ImageJ Coloc2 plugin.

**colocRhoAF.csv** RhoA/F-actin colocalization data from long-term cultivation of WJMSC-1 cells. Probe - passage number. Rval,tM1,tM2,bTau,Rs - colocalization coefficients calculated in ImageJ Coloc2 plugin.

**colocRhoAH.csv** RhoA/F-actin colocalization data from long-term cultivation of WJMSC-1 cells. Probe - passage number. Rval,tM1,tM2,bTau,Rs - colocalization coefficients calculated in ImageJ Coloc2 plugin.

**coloc_aA4-act.csv** Alpha-actinin-4/F-actin colocalization data from long-term cultivation of WJMSC-1 cells. Probe - passage number. Rval,tM1,tM2,bTau,Rs - colocalization coefficients calculated in ImageJ Coloc2 plugin.

**coloc_aA4-nuc.csv** Alpha-actinin-4/Hoeochst33444 colocalization data from long-term cultivation of WJMSC-1 cells. Probe - passage number. Rval,tM1,tM2,bTau,Rs - colocalization coefficients calculated in ImageJ Coloc2 plugin.

**wj1fracmod.csv** Data for Actin cytoskeleton fractality analysis. LCFD for F-actin was obtained from confocal images of rhodamine-phalloidin stained MSCWJ-1 cells using ImageJ FracLac plugin.

**wj1shape.csv** Data for cell shape analysis. Cell shape parameters were obtained from CQ1 confocal images of MSCWJ-1 cells at different passages.
