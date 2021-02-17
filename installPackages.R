library('devtools')
# install fails without this options call:
# https://github.com/r-lib/devtools/issues/647
options(unzip = 'internal')
install_github("everettJK/gt23")
