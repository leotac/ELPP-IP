#-------------------------
# Cplex and other libraries
# ------------------------

SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio126/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio126/concert

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CCLNFLAGS = -lconcert -lilocplex -lcplex -lm -pthread -lemon

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

CPLEXINC = -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)

#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	elpp
MAINOBJ		=	main.o elpp.o separation.o 

MAIN		=	$(MAINNAME)
MAINFILE	=	$(BINDIR)/$(MAIN)
MAINSHORTLINK	=	$(BINDIR)/$(MAINNAME)
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))
USRCXXFLAGS	=	-std=c++0x -DIL_STD
			#^C++11    ^Required by Cplex 

OBJDIR          =       obj
BINOBJDIR       =       $(OBJDIR)/bin
SRCDIR          =       src
BINDIR          =       bin
CXX		=	g++
CXX_o           =       -o # the trailing space is important
LINKCXX         =       g++
LINKCXX_o       =       -o # the trailing space is important
FLAGS 		= 	-I$(SRCDIR) -DNDEBUG -DROUNDING_FE $(USRCXXFLAGS)
OFLAGS 		= 	-O3
CXXFLAGS        +=      -m64
LDFLAGS         +=      -m64

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(MAINOBJFILES) $(MAINSHORTLINK)
endif

.PHONY: all
all:		$(MAINFILE) $(MAINSHORTLINK)

$(OBJDIR):
		@-mkdir -p $(OBJDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

.PHONY: clean
clean:		$(OBJDIR)
ifneq ($(OBJDIR),)
		-rm -f $(OBJDIR)/*.o
		-rmdir $(OBJDIR)
endif
		-rm -f $(MAINFILE)

$(MAINFILE):	$(BINDIR) $(OBJDIR)  $(LPILIBFILE) $(NLPILIBFILE) $(MAINOBJFILES)
		@echo "-> linking $@"
		$(LINKCXX) $(CPLEXINC) $(MAINOBJFILES)\
		$(CCLNDIRS) \
                $(OFLAGS) $(CCLNFLAGS) \
		$(LDFLAGS) $(LINKCXX_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(CPLEXINC) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c $< $(CXX_o)$@

#---- EOF --------------------------------------------------------------------
