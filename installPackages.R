library('devtools')
# install fails without this options call:
# https://github.com/r-lib/devtools/issues/647
options(unzip = 'internal')
install_github("everettJK/gt23")
install.packages('backports', repos='http://cran.us.r-project.org')
install.packages('checkmate', repos='http://cran.us.r-project.org')
