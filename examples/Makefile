# %---------------------------------%
# |  SECTION 1: SOQCS LIBRARY PATH  |
# %---------------------------------%

LIBN = soqcs
LPATH= ./src/soqcs

# %---------------------------------%
# |  SECTION 2: FLAGS C++           |
# %---------------------------------%

GPP = $(CXX)
FLAGS= -O2 -Wall -pedantic
UNAME := $(shell uname -s)
ifeq ($(UNAME),Linux)
	FLAGS += -std=c++17 -fopenmp -fPIC
endif
ifeq ($(UNAME),Darwin)
	arch := $(shell uname -p)
        FLAGS += -std=c++1z
endif

#  %--------------------------------%
#  |  SECTION 3: THE ARCHIVER AR    |
#  %--------------------------------%

AR = ar
ARFLAGS = rcs


#  %--------------------------------%
#  |  SECTION 4: DON'T TOUCH        |
#  %--------------------------------%
LIBB = .$(LPATH)
LIBINCLUDE = .$(LPATH)
LIBP= pthread

CFLAGS += -I$(LIBINCLUDE)
LDFLAGS+= -L$(LIBB)
LDFLAGS += -l$(LIBN)
LDFLAGS += -l$(LIBP)


#  %--------------------------------%
#  |  SECTION 5: COMPILATION RULES  |
#  %--------------------------------%

all:	     benchmark live1 live2 live3 live4 live5 live6 live7 live8 live9

benchlive:  
	     $(GPP) benchlive.cpp $(CFLAGS) $(LDFLAGS) $(FLAGS) -o benchlive.x

live1:  
	     $(GPP) live1.cpp $(CFLAGS) $(LDFLAGS) $(FLAGS) -o live1.x

live2:  
	     $(GPP) live2.cpp $(CFLAGS) $(LDFLAGS) $(FLAGS) -o live2.x

live3:  
	     $(GPP) live3.cpp $(CFLAGS) $(LDFLAGS) $(FLAGS) -o live3.x

live4:  
	     $(GPP) live4.cpp $(CFLAGS) $(LDFLAGS) $(FLAGS) -o live4.x
	     
live5:  
	     $(GPP) live5.cpp $(CFLAGS) $(LDFLAGS) $(FLAGS) -o live5.x

live6:  
	     $(GPP) live6.cpp $(CFLAGS) $(LDFLAGS) $(FLAGS) -o live6.x

live7:  
	     $(GPP) live7.cpp $(CFLAGS) $(LDFLAGS) $(FLAGS) -o live7.x

live8:  
	     $(GPP) live8.cpp $(CFLAGS) $(LDFLAGS) $(FLAGS) -o live8.x

live9:  
	     $(GPP) live9.cpp $(CFLAGS) $(LDFLAGS) $(FLAGS) -o live9.x	     
	     
clean:
	     rm *.x
