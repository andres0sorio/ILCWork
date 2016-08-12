SHELL = /bin/sh
UNAME = $(shell uname)

# Root variables
ROOTCFLAGS   = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS     = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --libs) 
ROOTGLIBS    = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --glibs)

################
##linux
CXX          = g++
RM           = rm -f 
AR	     = ar cru


##Flags
CXXFLAGS     = -O0 -g -Wall -fPIC -DGZSTREAM_NAMESPACE=gz
#CXXFLAGS     += -D_DEBUG
#CXXFLAGS     += -D_WARN
#CXXFLAGS     += -D_DEBUGKINFIT

#############
# Directories
INCDIR       = ../inc
SRCDIR       = ../src
EXEDIR       = ../exe
LIBDIR       = ../lib

#################
##Dependencies
#gzstream and zlib
GZSTREAM_DIR     = $(HOME)/opt/gzstream
CXXFLAGS    	+= -I$(GZSTREAM_DIR) 
GZSLIB       	 = -L$(GZSTREAM_DIR) -lgzstream
ZLIB         	 = -L/usr/lib -lz
###################################
### Linux
ifeq "$(UNAME)" "Linux"
RANLIB       = ranlib
CXXFLAGS    += -I$(KTJET_DIR)
CXXFLAGS    += -I$(CLHEP_DIR)/include -I$(INCDIR) $(ROOTCFLAGS)
KTLIB        = -L$(KTJET_DIR)/lib -lKtEvent 
CLHEPLIB     = -L$(CLHEP_DIR)/lib -lCLHEP-g++.1.8.0.0
LINKFLAGS    =
endif
##################################
###MAC OSX
ifeq "$(UNAME)" "Darwin"
RANLIB       = ranlib -s
#CXXFLAGS    += -Wunreachable-code
KTJET_DIR    = $(HOME)/HepSoft/ktjet_1.02
CXXFLAGS    += -I$(KTJET_DIR) -I$(CLHEP_DIR)/include -I$(ROOTSYS)/include -I$(INCDIR)
KTLIB        = -L$(KTJET_DIR)/lib -lKtEvent.1.02
CLHEPLIB     = -L$(CLHEP_DIR)/lib -lCLHEP-CC.1.8.2.0
LINKFLAGS    = -Xlinker -flat_namespace -bind_at_load
endif
###################################
##Libraries
LIBS       += $(ROOTLIBS) -lMinuit $(KTLIB) $(CLHEPLIB) $(GZSLIB) $(ZLIB)
#GLIBS      += $(ROOTGLIBS)
##############################################
## AnalysisTools RootTools GraphTools libaries

MYLIBS      = -L$(LIBDIR) -lAnalTools -lRootTools -lGraphTools
#######################################

LIB1OBJS    = SimdetStruct.o SimdetEvent.o SimdetLoader.o\
	       EFlowObject.o CutsLoader.o\
	       CombTable.o ParamLoader.o KfitsRoot.o\
               anal_utilities.o clhep_utilities.o

LIB2OBJS    = string_utilities.o root_utilities.o

LIB3OBJS    = BasicOutput.o SpecificOutput.o\
               PythiaSpecificOutput.o WhizardSpecificOutput.o\
	       HadronicOutput.o HadronicOutputOrganizer.o
#	       LeptonicOutput.o LeptonicOutputOrganizer.o

EXECOBJS1   = wwsHadronicA.o wwsHadronEventProc.o
#EXECOBJS2   = wwsLeptonicA.o wwsLeptonEventProc.o
EXECOBJS3   = etaphiEt.o
EXECOBJS4   = RunHadronicTest.o wwsHadronEventProc.o TestClass.o
EXECOBJS5   = RunMCanal.o wwsHadronEventProc.o MCData.o

GARBAGE     = $(EXEDIR)/doHadronAnalysis $(EXEDIR)/RunHadronicTest
GARBAGE    += $(EXEDIR)/RunMCanal #$(EXEDIR)/doLeptonAnalysis

########################################################################
.PHONY: libs exe

all:
	@echo make "[libs|exe|world|etrans|tests]  [clean:cleanall]"

exe: exec1 exec4 exec5 #exec2

etrans: exec3

tests: exec4

world: libs exe

########################################################################

%.o : %.cc 
	$(CXX) $(CXXFLAGS) -c $< 	

exec1 : $(EXECOBJS1)
	$(CXX) -o $(EXEDIR)/doHadronAnalysis $(LINKFLAGS) $(EXECOBJS1) $(MYLIBS) $(LIBS)

#
#exec2 : $(EXECOBJS2)
#	$(CXX) -o $(EXEDIR)/doLeptonAnalysis $(LINKFLAGS) $(EXECOBJS2) $(MYLIBS) $(LIBS)
#

exec3 : $(EXECOBJS3)
	 $(CXX) -o $(EXEDIR)/etaphiEt $(LINKFLAGS) $(EXECOBJS3) $(MYLIBS) $(LIBS)

exec4 : $(EXECOBJS4)
	 $(CXX) -o $(EXEDIR)/RunHadronicTest $(LINKFLAGS) $(EXECOBJS4) $(MYLIBS) $(LIBS)

exec5 : $(EXECOBJS5)
	 $(CXX) -o $(EXEDIR)/RunMCanal $(LINKFLAGS) $(EXECOBJS4) $(MYLIBS) $(LIBS)

libs : $(LIB1OBJS) $(LIB2OBJS) $(LIB3OBJS)
	$(AR) $(LIBDIR)/libAnalTools.a $(LIB1OBJS)
	$(RANLIB) $(LIBDIR)/libAnalTools.a
	$(AR) $(LIBDIR)/libRootTools.a $(LIB2OBJS)
	$(RANLIB) $(LIBDIR)/libRootTools.a
	$(AR) $(LIBDIR)/libGraphTools.a $(LIB3OBJS)
	$(RANLIB) $(LIBDIR)/libGraphTools.a

clean	:
	$(RM) $(GARBAGE) $(EXECOBJS1) $(EXECOBJS2)

cleanall:
	$(RM) $(GARBAGE) *.o $(LIBDIR)/libAnalTools.a $(LIBDIR)/libRootTools.a $(LIBDIR)/libGraphTools.a
