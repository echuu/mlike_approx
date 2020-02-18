

# installing rstan package from source -----------------------------------------

# https://github.com/stan-dev/rstan/wiki/Installing-RStan-from-source-on-Windows


install.packages("fansi")
pkgbuild::has_build_tools(debug = TRUE)

dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, "Makevars.win")
if (!file.exists(M)) file.create(M)
cat("\nCXX14FLAGS=-O3 -march=corei7 -mtune=corei7",
    "CXX14 = $(BINPREF)g++ -m$(WIN) -std=c++1y",
    "CXX11FLAGS=-O3 -march=corei7 -mtune=corei7",
    file = M, sep = "\n", append = TRUE)


remove.packages("rstan")
if (file.exists(".RData")) file.remove(".RData")

install.packages("rstan", type = "source")

# enable option to run models in parallel
options(mc.cores = parallel::detectCores())

# save compiled stan program to hard disk so that don't need to recompile
rstan_options(auto_write = TRUE)

Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7')


# installing rstanarm (must be done after rstan installation) ------------------

install.packages("rstanarm")
#install.packages("bayesplot")
#install.packages("broom")











