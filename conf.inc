# %---------------------------------%
# |         CONFIGURATION           |
# %---------------------------------%
LPATH=./src
PPATH=./py_src
EPATH=./examples


# %---------------------------------%
# |  SECTION 1: PATHS AND LIBRARIES |
# %---------------------------------%

LIBN = soqcs
LIBB = .$(LPATH)
LIBINCLUDE = .$(LPATH)
LIBP= pthread

# %---------------------------------%
# |  SECTION 2: COMPILERS C++       |
# %---------------------------------%

GPP = g++
FLAGS= -O2 -std=c++17 -fopenmp -fPIC

#  %--------------------------------%
#  |  SECTION 3: THE ARCHIVER AR    |
#  %--------------------------------%

AR = ar
ARFLAGS = rcs


#  %--------------------------------%
#  |  SECTION 4: DON'T TOUCH        |
#  %--------------------------------%
CFLAGS += -I$(LIBINCLUDE)
LDFLAGS+= -L$(LIBB)
LDFLAGS += -l$(LIBN)
LDFLAGS += -l$(LIBP)


