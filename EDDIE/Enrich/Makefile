

#Location to gsl headers
GSLCFLAGS     = -I/usr/include/gsl

#### local macros
# remove without fussing about it
RM = /bin/rm -f

# compiler name and flags
CCC = g++
CCFLAGS = -O3 -fomit-frame-pointer -funroll-loops -fforce-addr -fexpensive-optimizations -Wno-deprecated -g3 ${GSLCFLAGS} -std=c++11

### local program information
EXEC=run
SOURCES=fileReader.C buildSets.C NetworkEnrichment.C Driver.C

### intermediate objects
OBJECTS = $(SOURCES: .cc=.o)

### libraries
LIBS += -lgsl #-lgslcblas

### headers
HEADERS = Headers.h readFile.h buildSets.h NetworkEnrichment.h

### targets, dependencies and actions
$(EXEC): $(OBJECTS) Makefile
	$(LINK.cc) $(CCFLAGS) -o $(EXEC) $(OBJECTS) $(LIBS)

### sort out dependencies
depend:
	makedepend $(INCLUDES) $(HEADERS) $(SOURCES) $(LIBS)

clean: 
	$(RM) $(EXEC) *.o *~
