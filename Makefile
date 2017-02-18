#@file    Makefile
#@brief   Makefile for ELPP tests 

#-------------------------
# Cplex
# ------------------------

SYSTEM     = x86-64_sles10_4.1
LIBFORMAT  = static_pic

CPLEXDIR      = /opt/ibm/cplex
CONCERTDIR    = /opt/ibm/concert

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
MAINOBJ		=	main.o elpp.o separation.o type.o 

MAIN		=	$(MAINNAME)
MAINFILE	=	$(BINDIR)/$(MAIN)
MAINSHORTLINK	=	$(BINDIR)/$(MAINNAME)
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))
USRCXXFLAGS	=	-std=c++0x -DIL_STD
			#^C++11    ^Required by Cplex 

OBJDIR          =       obj
BINOBJDIR       =       $(OBJDIR)/bin
LIBOBJDIR       =       $(OBJDIR)/lib
SRCDIR          =       src
BINDIR          =       bin
LIBDIR          =       lib
CXX		=	g++
CXX_c           =       -c # the trailing space is important
CXX_o           =       -o # the trailing space is important
LINKCXX         =       g++
LINKCXX_L       =       -L
LINKCXX_l       =       -l
LINKCXX_o       =       -o # the trailing space is important
FLAGS 		= 	-I$(SRCDIR) -DNDEBUG -DROUNDING_FE $(USRCXXFLAGS)
OFLAGS 		= 	-O3
CFLAGS          +=      -m64
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
                $(OFLAGS) $(LPSLDFLAGS) $(CCLNFLAGS) \
		$(LDFLAGS) $(LINKCXX_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(CPLEXINC) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c $< $(CXX_o)$@

#---- EOF --------------------------------------------------------------------
