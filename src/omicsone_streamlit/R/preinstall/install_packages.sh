Rscript -e "if (!requireNamespace('dplyr', quietly = TRUE)) install.packages('dplyr', repos='http://cran.r-project.org')"
Rscript -e "if (!requireNamespace('plyr', quietly = TRUE)) install.packages('plyr', repos='http://cran.r-project.org')"
Rscript -e "if (!requireNamespace('tidyr', quietly = TRUE)) install.packages('tidyr', repos='http://cran.r-project.org')"
Rscript -e "if (!requireNamespace('tibble', quietly = TRUE)) install.packages('tibble', repos='http://cran.r-project.org')"
Rscript -e "if (!requireNamespace('data.table', quietly = TRUE)) install.packages('data.table', repos='http://cran.r-project.org')"
Rscript -e "if (!requireNamespace('NMF', quietly = TRUE)) install.packages('NMF', repos='http://cran.r-project.org')"
Rscript -e "if (!requireNamespace('jsonlite', quietly = TRUE)) install.packages('jsonlite', repos='http://cran.r-project.org')"
Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos='http://cran.r-project.org')"
Rscript -e "if (!requireNamespace('ComplexHeatmap', quietly = TRUE)) BiocManager::install('ComplexHeatmap')"
Rscript -e "if (!requireNamespace('Biobase', quietly = TRUE)) BiocManager::install('Biobase')"


