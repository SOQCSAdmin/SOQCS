# %---------------------------------%
# |  SECTION 1: FLAGS C++           |
# %---------------------------------%

GPP = $(CXX)
UNAME := $(shell uname -s)
ifeq ($(UNAME),Linux)
    EIGEN3=/usr/include/eigen3/
    FLAGS= -O2 -Wall -pedantic -I $(EIGEN3)
	FLAGS += -std=c++17 -fopenmp -fPIC
endif
ifeq ($(UNAME),Darwin)
    EIGEN3=/usr/local/include/eigen3/
    FLAGS= -O2 -Wall -pedantic -I $(EIGEN3)
    FLAGS += -std=c++1z
endif

#  %--------------------------------%
#  |  SECTION 3: THE ARCHIVER AR    |
#  %--------------------------------%

AR = ar
ARFLAGS = rcs
