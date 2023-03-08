PRES_DIR = .
WORK_DIR = ./bi
console_DIR = ./console
DIRG = ${PRES_DIR}/general
DIRR = ${PRES_DIR}/randomc
DIRS = ${PRES_DIR}/special
norg  = ${WORK_DIR}/norg
VPATH = ${DIRG}:${DIRR}:${DIRS}

CC       = mpiicpc
CPPFLAGS = -std=c++17 -qmkl -O3
# CPPFLAGS = -std=c++17 -mkl -O3
CFLAGS   =
CXXFLAGS = $(CFLAGS)
COMPILE  = $(CC) $(CPPFLAGS) $(CXXFLAGS) -c

SRCS := $(wildcard ${DIRG}/*.cpp) $(wildcard ${DIRR}/*.cpp) $(wildcard ${DIRS}/*.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS))
DEPS := $(patsubst %.cpp,%.d,$(SRCS))

default: manpower

srun: ${EXE}
	@echo '______________________________________________________________________________________________ sbatch _______'
	@echo ''
	@qsub ${console_DIR}/slurm.bscc.t6.txt
	@echo '_____________________________________________________________________________________________________________'
	@echo ''

qsub: ${norg}
	@echo '______________________________________________________________ qsub ______________________________________________________________'
	@echo ''
	@qsub ${WORK_DIR}/qsub.DPC++CPU.mpi.txt > ${WORK_DIR}/jmwang.job.txt
	@chmod +x ${WORK_DIR}/job.process.DPC++CPU.txt
	@${WORK_DIR}/job.process.DPC++CPU.txt
	@rm -rf ${WORK_DIR}/jmwang.job.txt
	@echo '__________________________________________________________________________________________________________________________________'
	@echo ''

${norg}: compile $(DEPS) $(OBJS)
	@echo '_____________________________________________________________ link _____________________________________________________________'
	@echo ''
	$(CC) $(CPPFLAGS) $(CXXFLAGS) -o ${norg} $(OBJS) $(LIBS)

manpower: ${norg}
	@echo '_________________________________________ Now you can test by: mpirun -n 24 ./console/norg________________________________________'

start:
	@echo '______________________________________________________________ start ______________________________________________________________'
	@echo ''

explain:
	@echo '______________________________________________________________ explain ______________________________________________________________'
	@echo ''
	@echo "the following information represents your prgram"
	@echo "final norg name: $(norg)"
	@echo "source files: $(SRCS)"
	@echo "object files: $(OBJS)"
	@echo "dependency files: $(DEPS)"

compile:
	@echo '____________________________________________________________ compile ____________________________________________________________'
	@echo ''

%.d:%.cpp
	$(CC) -MM $(CPPFLAGS) $< > $@
	$(CC) -MM $(CPPFLAGS) $< | sed s/\\.o/\\.d/   >> $@

%.o:%.cpp
	$(COMPILE) -o $@ $<

clean:
	-rm -rf $(OBJS) $(DEPS) $(norg)

clear:

	-rm -rf bi/zic*	
	-rm -rf io/output.*
	-rm -rf io/*.txt
	# -rm -rf tso/*
	-rm -rf jmwang.*
	-rm -rf *jmw*
	-rm -rf log.*
deepclean:
	-rm -rf general/*.cpp
	-rm -rf general/*.h
	-rm -rf special/*.cpp
	-rm -rf special/*.h
	-rm -rf randomc/*.cpp
	-rm -rf randomc/*.h
	

depend:$(DEPS)
	@echo "dependencies are now up-to-date."

-include $(DEPS)
