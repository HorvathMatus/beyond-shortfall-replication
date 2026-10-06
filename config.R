# ============================================================
# Beyond Shortfall — paths used by the numbered scripts
#
# Every script sources this file; run the scripts from the folder that
# contains it (the root of the replication package).
#   MATLAB_DIR  folder with the MATLAB code of steps 02, 06, 10, 14 and 18.
#               The R scripts write the MATLAB input files there and read
#               the MATLAB output files from it.
#   TEX_DIR     LaTeX tables written by the R scripts
#   FIG_DIR     figures written by the R scripts
# Intermediate R results are written to data/.
# ============================================================
MATLAB_DIR <- normalizePath("matlab")
TEX_DIR    <- "output/tables"
FIG_DIR    <- "output/figures"
for (d_ in c("data", TEX_DIR, FIG_DIR)) dir.create(d_, showWarnings = FALSE, recursive = TRUE)
rm(d_)
