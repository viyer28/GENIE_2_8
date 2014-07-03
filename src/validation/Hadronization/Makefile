#
# Makefile for the GENIE hadronization validation package
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
include $(GENIE)/src/make/Make.include

CXXFLAGS += -g

GENIE_LIBS   = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
LIBRARIES   := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

HADRONZ_TEST_SOURCES   := RunConfig.cxx  Analyzer.cxx ExpData.cxx 
HADRONZ_TEST_SOURCES   += GenPlotsBase.cxx CommonCalc.cxx 
HADRONZ_TEST_SOURCES   += ChHadW2.cxx NegChHadD.cxx ChHadDoverMultW2.cxx 
HADRONZ_TEST_SOURCES   += Pi0W2.cxx Pi0D.cxx
HADRONZ_TEST_SOURCES   += gVldHadronzTest.cxx
HADRONZ_TEST_OBJECTS    = $(addsuffix .$(ObjSuf), $(basename $(HADRONZ_TEST_SOURCES)))

all:    gVldHadronzTest 

gVldHadronzTest: FORCE
	$(CXX) $(CXXFLAGS) -c $(HADRONZ_TEST_SOURCES) $(INCLUDES)
#	$(LD) $(LDFLAGS) $(HADRONZ_TEST_OBJECTS) $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_hadronz_test
	$(LD) $(LDFLAGS) $(HADRONZ_TEST_OBJECTS) $(LIBRARIES) -o gvld_hadronz_test

purge: FORCE
	$(RM) *.o *~ core

clean: FORCE
	$(RM) *.o *~ core
#	$(RM) $(GENIE_BIN_PATH)/gvld_hadronz_test
	$(RM) ./gvld_hadronz_test

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_hadronz_test

FORCE:

# DO NOT DELETE


