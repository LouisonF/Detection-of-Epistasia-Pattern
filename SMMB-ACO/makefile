#===============================================================================
# Makefile for Unix & Linux Systems using a GNU C++ compiler 		       #
#===============================================================================

#===============================================================================
# Compiler flags
#===============================================================================
# -g --Enable debugging
# -Wall --Turn on all warnings
# -fopenmp to use parallelized computation
# -Wall for all warnings
# -Wextra for extra warnings
# -DNDEBUG for using debugging tools such as dnd or valgrind and profiling of our program if those informations are don't wanted please use -s instead
# -I$(BOOST_FOLDER) include all libraries
# -O1 -O2 -O3  offers three levels of optimization, from least optimized to fastest. The program generated will generally be larger and, above all, the compilation time will increase as you move up levels.


# Add external library path
BOOST_FOLDER = /usr/include/boost/
# Give path of .hpp files
INC=-I$(BOOST_FOLDER) \
 -I./includes \

# Save compiler call
 CXX=g++ -o
 #Intial call
 # Flags for less verbose, include parallelized compute, optimization at level 2
 CXXFLAGS=-std=c++11 -s -O2 -fopenmp -g $(INC)
 LFLAGS=-std=c++11 -fopenmp -I$(BOOST_FOLDER) -lm
#  Flags with maximum verbose, include parallelized compute, optimization at level 2, debugger informations
 #CXXFLAGS=-std=c++11 -fopenmp -Wall -Wextra -DNDEBUG -O2 -g $(INC)
 #LFLAGS=-std=c++11 -fopenmp -Wall -Wextra -DNDEBUG -I$(BOOST_FOLDER) -lm
# Save directory names, store some shell commands
 SRCDIR=src
 OBJDIR=obj
 BINDIR=.
 LOGDIR=log
 rm=rm -f
 mkdir=mkdir -p

#Give .o files to create
OBJ=$(OBJDIR)/main.o \
$(OBJDIR)/Parametersfileparsing.o \
$(OBJDIR)/datainput.o \
$(OBJDIR)/Smmb_ACO.o \
$(OBJDIR)/Contingency.o \
$(OBJDIR)/G2_conditional_test_indep.o \
$(OBJDIR)/G2_test_indep.o \
$(OBJDIR)/Miscellaneous.o \


TARGET=$(BINDIR)/smmb_aco

all:$(TARGET)

$(TARGET):$(OBJ)
	@g++ $^ -o $@ $(LFLAGS)
	@echo "\033[32m"+"Compilation complete, executable produced.""\033[0m"

#===============================================================================
# Main
#===============================================================================
$(OBJDIR)/main.o: ./$(SRCDIR)/main.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "\033[32m"+"Compiled main.o""\033[0m"

#===============================================================================
# Parametersfileparsing
#===============================================================================
$(OBJDIR)/Parametersfileparsing.o: ./$(SRCDIR)/Parametersfileparsing.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "\033[32m"+"Compiled Parametersfileparsing.o""\033[0m"

$(OBJDIR)/datainput.o: ./$(SRCDIR)/datainput.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "\033[32m"+"Compiled datainput.o""\033[0m"
#===============================================================================
# Smmb_aco
#===============================================================================
$(OBJDIR)/Smmb_ACO.o: ./$(SRCDIR)/SmmbACO.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "\033[32m"+"Compiled Smmb_ACO.o""\033[0m"

#===============================================================================
# Statistics
#===============================================================================
$(OBJDIR)/Contingency.o: ./$(SRCDIR)/Contingency.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "\033[32m"+"Compiled Contingency.o""\033[0m"

$(OBJDIR)/G2_conditional_test_indep.o: ./$(SRCDIR)/G2_conditional_test_indep.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "\033[32m"+"Compiled G2_conditional_test_indep.o""\033[0m"

$(OBJDIR)/G2_test_indep.o: ./$(SRCDIR)/G2_test_indep.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "\033[32m"+"Compiled G2_test_indep.o""\033[0m"

#===============================================================================
# Global
#===============================================================================
$(OBJDIR)/Miscellaneous.o: ./$(SRCDIR)/Miscellaneous.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "\033[32m"+"Compiled Miscellaneous.o""\033[0m"



#===============================================================================
# PHONEY
#===============================================================================
.PHONEY: clean
clean:
	@$(rm) $(OBJ)
	@echo "\n\033[32m"+"Cleanup complete."+"\033[0m\n"

.PHONEY: remove
remove: clean
	@$(rm) $(TARGET)
	@echo "\n\033[32m"+"Executable removed."+"\033[0m\n"

.PHONEY: install
install:
	@$(mkdir) $(OBJDIR)
	@$(mkdir) $(SRCDIR)
	@$(mkdir) $(LOGDIR)
	@echo "\n\033[32m"+"Directories created"+"\033[0m\n"
