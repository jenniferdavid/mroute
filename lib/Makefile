HOME          = .

CXX           = g++

INCLUDE       = 

#CXXFLAGS      = -O -fno-implicit-templates $(INCLUDE) 
CXXFLAGS      = -O4 $(INCLUDE)
#CXXFLAGS      = -g1 -Wall -w $(INCLUDE) 
#CXXFLAGS      = -g -Wall -fno-implicit-templates $(INCLUDE) 
#CXXFLAGS      = -g  $(INCLUDE) 

DEST	      = .

EXTHDRS	      =

HDRS	      = NN_class.hh \
		Option.hh \
		PottsNeuron.hh \
		Propagator.hh \
		itoa.hh \
		link.hh \
		manswer.hh \
		mreq.hh \
		mroute.hh \
		net.hh \
		node.hh \
		random.hh \
		rmreq.hh \
		rnet.hh

INSTALL	      = cp

LIBRARY	      = libmroute.a

MAKEFILE      = Makefile

OBJS	      = NN_alg.o \
		NN_class.o \
		Option.o \
		PottsNeuron.o \
		Propagator.o \
		bellman.o \
		itoa.o \
		link.o \
		manswer.o \
		mreq.o \
		net.o \
		node.o \
		plot.o \
		random.o \
		rmreq.o \
		rnet.o \
		seq_bellman.o \
		solveexact.o

PRINT	      = maser -4

SHELL	      = /bin/sh

SRCS	      = NN_alg.cc \
		NN_class.cc \
		Option.cc \
		PottsNeuron.cc \
		Propagator.cc \
		bellman.cc \
		itoa.cc \
		link.cc \
		manswer.cc \
		mreq.cc \
		net.cc \
		node.cc \
		plot.cc \
		random.cc \
		rmreq.cc \
		rnet.cc \
		seq_bellman.cc \
		solveexact.cc

SYSHDRS	      = /usr/include/getopt.h \
		/usr/include/kern/assert.h \
		/usr/include/kern/event.h \
		/usr/include/kern/lock.h \
		/usr/include/kern/macro_help.h \
		/usr/include/kern/processor.h \
		/usr/include/kern/queue.h \
		/usr/include/kern/sched.h \
		/usr/include/kern/spl.h \
		/usr/include/mach/alpha/vm_param.h \
		/usr/include/mach/boolean.h \
		/usr/include/mach/machine/boolean.h \
		/usr/include/mach/machine/vm_param.h \
		/usr/include/mach/machine/vm_statistics.h \
		/usr/include/mach/machine/vm_types.h \
		/usr/include/mach/policy.h \
		/usr/include/mach/port.h \
		/usr/include/mach/vm_prot.h \
		/usr/include/mach/vm_statistics.h \
		/usr/include/machine/context.h \
		/usr/include/machine/counters.h \
		/usr/include/machine/cpu.h \
		/usr/include/machine/endian.h \
		/usr/include/machine/machlimits.h \
		/usr/include/machine/machparam.h \
		/usr/include/machine/machtime.h \
		/usr/include/machine/pmap.h \
		/usr/include/machine/signal.h \
		/usr/include/machine/spl.h \
		/usr/include/math.h \
		/usr/include/random.h \
		/usr/include/standards.h \
		/usr/include/stdio.h \
		/usr/include/stdlib.h \
		/usr/include/string.h \
		/usr/include/strings.h \
		/usr/include/sys/callout.h \
		/usr/include/sys/context_t.h \
		/usr/include/sys/limits.h \
		/usr/include/sys/lwc.h \
		/usr/include/sys/param.h \
		/usr/include/sys/procset.h \
		/usr/include/sys/resource.h \
		/usr/include/sys/rt_limits.h \
		/usr/include/sys/seek.h \
		/usr/include/sys/select.h \
		/usr/include/sys/siginfo.h \
		/usr/include/sys/signal.h \
		/usr/include/sys/syslimits.h \
		/usr/include/sys/sysmisc.h \
		/usr/include/sys/time.h \
		/usr/include/sys/timers.h \
		/usr/include/sys/types.h \
		/usr/include/sys/unix_defs.h \
		/usr/include/sys/wait.h \
		/usr/include/time.h \
		/usr/include/va_list.h

all:		$(LIBRARY)

$(LIBRARY):	$(OBJS) templates/mroute_templates.o
		@echo "Loading $(LIBRARY) ..."
		ar cru $(LIBRARY) templates/mroute_templates.o $(OBJS)
		@echo "done"

templates/mroute_templates.o: templates/mroute_templates.cc
		cd templates ; gmake mroute_templates.o

clean:;		@rm -f $(OBJS) core

clobber:;	@rm -f $(OBJS) $(LIBRARY) core tags *~
		@(cd templates ; rm -f *.o core *~)

depend:;	@mkmf -f $(MAKEFILE)

echo:;		@echo $(HDRS) $(SRCS)

extract:;	@ar x $(DEST)/$(LIBRARY)

index:;		@ctags -wx $(HDRS) $(SRCS)

install:	$(LIBRARY)
		@echo Installing $(LIBRARY) in $(DEST)
		@if [ $(DEST) != . ]; then \
		(rm -f $(DEST)/$(LIBRARY); $(INSTALL) $(LIBRARY) $(DEST)); fi

print:;		@$(PRINT) $(HDRS) $(SRCS)

tags:           $(HDRS) $(SRCS); @ctags $(HDRS) $(SRCS)
###
NN_alg.o: mroute.hh NN_class.hh /usr/include/stdio.h \
	PottsNeuron.hh /usr/include/stdlib.h Propagator.hh manswer.hh mreq.hh net.hh link.hh \
	node.hh Option.hh
NN_class.o: /usr/include/math.h \
	/usr/include/stdlib.h \
	PottsNeuron.hh Propagator.hh mroute.hh manswer.hh mreq.hh net.hh \
	link.hh node.hh Option.hh
Option.o: Option.hh /usr/include/stdlib.h random.hh
PottsNeuron.o: PottsNeuron.hh /usr/include/stdlib.h random.hh
Propagator.o: Propagator.hh /usr/include/stdlib.h PottsNeuron.hh
bellman.o: manswer.hh mroute.hh mreq.hh net.hh link.hh node.hh
itoa.o: /usr/include/stdio.h itoa.hh
link.o: mroute.hh link.hh node.hh
manswer.o: manswer.hh mroute.hh mreq.hh net.hh link.hh node.hh Option.hh \
	rmreq.hh rnet.hh
mreq.o: /usr/include/math.h /usr/include/stdlib.h \
        mreq.hh mroute.hh net.hh \
	link.hh node.hh Option.hh random.hh itoa.hh
net.o: itoa.hh mroute.hh net.hh link.hh node.hh Option.hh random.hh
node.o: node.hh mroute.hh link.hh
plot.o: /usr/include/stdio.h random.hh net.hh link.hh node.hh mroute.hh \
	manswer.hh mreq.hh
random.o: random.hh 
rmreq.o: rmreq.hh net.hh link.hh node.hh mroute.hh mreq.hh rnet.hh
rnet.o: net.hh link.hh node.hh mroute.hh itoa.hh rnet.hh
seq_bellman.o: manswer.hh mroute.hh mreq.hh net.hh link.hh node.hh
solveexact.o: net.hh link.hh node.hh mroute.hh manswer.hh mreq.hh Option.hh
