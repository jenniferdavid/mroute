HOME          = .

CXX           = g++

INCLUDE       = -I/home/jendav/Videos/mroute/lib

#CXXFLAGS      = -g -w $(INCLUDE) 
CXXFLAGS      = -O -fno-implicit-templates $(INCLUDE) 
#CXXFLAGS      = -g -Wall -fno-implicit-templates $(INCLUDE) 

DEST	      = .

EXTHDRS	      =

HDRS	      =

INSTALL	      = cp

LD	      = g++

LDFLAGS	      = -L/usr/local/lib

LIBS	      = /home/jendav/Videos/mroute/lib/libmroute.a

MAKEFILE      = Makefile

OBJS	      = mroute.o

PRINT	      = maser -2 -b

PROGRAM       = mroute

SHELL	      = /bin/sh

SRCS	      = mroute.cc

SYSHDRS	      = /usr/include/getopt.h \
		/usr/include/kern/assert.h \
		/usr/include/kern/event.h \
		/usr/include/kern/lock.h \
		/usr/include/kern/macro_help.h \
		/usr/include/kern/queue.h \
		/usr/include/kern/spl.h \
		/usr/include/mach/boolean.h \
		/usr/include/mach/machine/boolean.h \
		/usr/include/mach/machine/vm_types.h \
		/usr/include/machine/spl.h \
		/usr/include/math.h \
		/usr/include/standards.h \
		/usr/include/stdlib.h \
		/usr/include/sys/select.h \
		/usr/include/sys/types.h

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) $(LIBS)
		@echo "Linking $(PROGRAM) ..."
		$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)
		@echo "done"

clean:;		@rm -f $(OBJS) core *~

clobber:;	@rm -f $(OBJS) $(PROGRAM) core tags *~

depend:;	@mkmf -f $(MAKEFILE)

echo:;		@echo $(HDRS) $(SRCS)

index:;		@ctags -wx $(HDRS) $(SRCS)

install:	$(PROGRAM)
		@echo Installing $(PROGRAM) in $(DEST)
		@-strip $(PROGRAM)
		@if [ $(DEST) != . ]; then \
		(rm -f $(DEST)/$(PROGRAM); $(INSTALL) $(PROGRAM) $(DEST)); fi

print:;		@$(PRINT) $(HDRS) $(SRCS)

tags:           $(HDRS) $(SRCS); @ctags $(HDRS) $(SRCS)
###
mroute.o: /usr/include/math.h /usr/include/stdlib.h \
