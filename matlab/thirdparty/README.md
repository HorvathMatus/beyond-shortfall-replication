# Third-party MATLAB code

The MATLAB code of this package builds on code written by others, which is not
redistributed here. Place it in this folder (any subfolder structure works,
`bs_paths.m` adds the folder recursively) or list its location in
`matlab/local_paths.m` (see `local_paths_example.m`).

## Replication files of Patton, Ziegel and Chen (2019)

From the replication package of Patton, A.J., J.F. Ziegel and R. Chen (2019),
"Dynamic semiparametric models for expected shortfall (and value-at-risk)",
*Journal of Econometrics* 211(2), 388–413, available from Andrew Patton's
website. The following files are used, unchanged:

    GAS_hybrid_LL3.m     GAS_onefactor_LL3.m   GAS_twofactor_LL3.m
    VEhatTAU.m           garch_FZ_LL.m         garchpq.m
    sample_VE.m          skewt_VE.m            skewtdis_LL.m
    skewtdis_cdf.m       skewtdis_inv.m

## Toolboxes

* **UCSD GARCH toolbox** (Kevin Sheppard): `garchlikelihood.m`, `garchcore.m`,
  `garchgrad.m`, called by `garchpq.m`.
* **Econometrics Toolbox** of James P. LeSage (http://www.spatial-econometrics.com),
  folder `distrib`: `tdis_inv.m`, `tdis_cdf.m`, `beta_cdf.m`, `beta_inv.m`,
  `beta_pdf.m`, `fdis_cdf.m` and `quantile.m`. Note that this `quantile.m`
  shadows the function of the same name in MATLAB's Statistics and Machine
  Learning Toolbox; both interpolate linearly so that F(x_(k)) = (k - 0.5)/n,
  so the empirical VaR in `sample_VE.m` is the same with either.
* **MFE Toolbox** (Kevin Sheppard, https://www.kevinsheppard.com/MFE_Toolbox):
  `normcdf.m`, `norminv.m`, `normpdf.m` (folder `duplication`) and
  `iscompatible.m`. These duplicate functions of the Statistics and Machine
  Learning Toolbox and were on the path when the results were produced.

The results in the paper were produced with all three toolboxes on the path in
the order MFE Toolbox, UCSD GARCH, Econometrics Toolbox (later folders take
precedence).
