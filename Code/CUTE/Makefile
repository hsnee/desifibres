####################################################
###          User definable stuff
#Monopole
NB_R = 64
R_MAX = 150.
#Angular CF
NB_THETA = 40
THETA_MAX = 10.
#3D CF (binning in sigma-pi)
NB_RL = 64
NB_RT = 64
RL_MAX = 150.
RT_MAX = 150.
#3D CF (binning in r-mu)
NB_R3D = 64
NB_CTH = 64
R3D_MAX = 150.
#Radial CF
NB_DZ = 64
DZ_MAX = 0.2
#Redshift range for full 2PCF
RED_0 = 0.43
RED_F = 0.70
NB_RED = 1
#Logarithmic binning
N_LOGINT = 10
#Behavior options
DEFINEOPTIONS = -D_VERBOSE
#DEFINEOPTIONS += -D_LOGBIN
#DEFINEOPTIONS += -D_DEBUG
#DEFINEOPTIONS += -D_TRUE_ACOS
DEFINEOPTIONS += -D_WITH_WEIGHTS
#GSL options
GSL_INC = /home/damonge/include
GSL_LIB = /home/damonge/lib
#CUDA options
NB_H2D = 64
OPT_PRECISION = -ftz=true -prec-div=false -prec-sqrt=false
CUDADIR = /usr/local/cuda
### End of user-definable stuff
####################################################

LGSL = -L$(GSL_LIB) -lgsl -lgslcblas

#Transform to inverse scales
I_R_MAX = $(shell echo "scale=5;1./$(R_MAX)" | bc)
I_THETA_MAX = $(shell echo "scale=5;57.29578/$(THETA_MAX)" | bc)
I_RL_MAX = $(shell echo "scale=5;1./$(RL_MAX)" | bc)
I_RT_MAX = $(shell echo "scale=5;1./$(RT_MAX)" | bc)
I_R3D_MAX = $(shell echo "scale=5;1./$(R3D_MAX)" | bc)
I_DZ_MAX = $(shell echo "scale=5;1./$(DZ_MAX)" | bc)
I_RED_INTERVAL = $(shell echo "scale=5;1./($(RED_F) - $(RED_0))" | bc)
LOG_R_MAX = $(shell echo "scale=9;l($(R_MAX))/l(10)" | bc -l)
LOG_TH_MAX = $(shell echo "scale=9;l($(THETA_MAX))/l(10)-1.75812263" | bc -l)
LOG_R3D_MAX = $(shell echo "scale=9;l($(R3D_MAX))/l(10)" | bc -l)

# DEFINES for the OpenMP version
DEFINEFLAGSCPU = $(DEFINEOPTIONS)
DEFINEFLAGSCPU += -DNB_R=$(NB_R) -DI_R_MAX=$(I_R_MAX) -DLOG_R_MAX=$(LOG_R_MAX)
DEFINEFLAGSCPU += -DNB_THETA=$(NB_THETA) -DI_THETA_MAX=$(I_THETA_MAX) -DLOG_TH_MAX=$(LOG_TH_MAX)
DEFINEFLAGSCPU += -DNB_RL=$(NB_RL) -DNB_RT=$(NB_RT) -DI_RL_MAX=$(I_RL_MAX) -DI_RT_MAX=$(I_RT_MAX)
DEFINEFLAGSCPU += -DNB_R3D=$(NB_R3D) -DNB_CTH=$(NB_CTH) -DI_R3D_MAX=$(I_R3D_MAX) -DLOG_R3D_MAX=$(LOG_R3D_MAX)
DEFINEFLAGSCPU += -DNB_DZ=$(NB_DZ) -DI_DZ_MAX=$(I_DZ_MAX)
DEFINEFLAGSCPU += -DI_RED_INTERVAL=$(I_RED_INTERVAL) -DRED_0=$(RED_0) -DNB_RED=$(NB_RED)
DEFINEFLAGSCPU += -D_HISTO_2D_$(NB_H2D) -DN_LOGINT=$(N_LOGINT)
#DEFINES for the CUDA version
DEFINEFLAGSGPU = $(DEFINEOPTIONS)
DEFINEFLAGSGPU += -DI_R_MAX=$(I_R_MAX) -DLOG_R_MAX=$(LOG_R_MAX)
DEFINEFLAGSGPU += -DI_THETA_MAX=$(I_THETA_MAX) -DLOG_TH_MAX=$(LOG_TH_MAX)
DEFINEFLAGSGPU += -DI_RL_MAX=$(I_RL_MAX) -DI_RT_MAX=$(I_RT_MAX)
DEFINEFLAGSGPU += -DI_R3D_MAX=$(I_R3D_MAX) -DLOG_R3D_MAX=$(LOG_R3D_MAX)
DEFINEFLAGSGPU += -D_HISTO_2D_$(NB_H2D) -DN_LOGINT=$(N_LOGINT)

# COMPILER AND OPTIONS
COMPCPU = gcc
COMPGPU = nvcc
OPTCPU = -Wall -O3 -fopenmp $(DEFINEFLAGSCPU)
OPTCPU_GPU = -Wall -O3 $(DEFINEFLAGSGPU)
OPTGPU = -O3 $(DEFINEFLAGSGPU) -arch compute_20 $(OPT_PRECISION) -Xcompiler -Wall

#INCLUDES AND LIBRARIES
INCLUDECOM = -I./src -I$(GSL_INC)
INCLUDECUDA = -I$(CUDADIR)/include
LIBCPU = $(LGSL) -lm -lhdf5
LIBGPU = $(LGSL) -L$(CUDADIR)/lib64 -lcudart -lpthread -lm

#.o FILES
#CUTE
DEF = src/define.o
COM = src/common.o
COSMO = src/cosmo.o
RANDOM = src/random.o
WOF = src/woftheta.o
CORR = src/correlator.o
BOX2D = src/boxes2D.o
BOX3D = src/boxes3D.o
IO = src/io.o
MAIN = src/main.c
OFILES = $(DEF) $(COM) $(COSMO) $(WOF) $(RANDOM) $(CORR) $(BOX2D) $(BOX3D) $(IO) $(MAIN)

#CU_CUTE
BOXCUDA = src/boxesCUDA.o
CORRCUDA = src/correlator_cuda.o
MAINCUDA = src/main_cuda.o
OFILESCUDA = $(DEF) $(COM) $(COSMO) $(MASK) $(CORRCUDA) $(BOXCUDA) $(RANDOM) $(IO) $(MAINCUDA)

#FINAL GOAL
EXE = CUTE
EXECUDA = CU_CUTE

#RULES
default : $(EXE) $(EXECUDA)
#RULE TO MAKE .o's FROM .c's
$(DEF) : src/define.c Makefile
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(COM) : src/common.c Makefile
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(COSMO) : src/cosmo.c Makefile
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(WOF) : src/woftheta.c Makefile
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(CORR) : src/correlator.c Makefile
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(RANDOM) : src/random.c Makefile
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(IO) : src/io.c Makefile
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(BOX2D) : src/boxes2D.c Makefile
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(BOX3D) : src/boxes3D.c Makefile
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(BOXCUDA) : src/boxesCUDA.c Makefile
	$(COMPCPU) $(OPTCPU) -c $< -o $@ $(INCLUDECOM)
$(CORRCUDA) : src/correlator_cuda.cu Makefile
	$(COMPGPU) $(OPTGPU) -c $< -o $@ $(INCLUDECOM) $(INCLUDECUDA)
$(MAINCUDA) : src/main_cuda.c Makefile
	$(COMPCPU) $(OPTCPU_GPU) -c $< -o $@ $(INCLUDECOM)

#RULES TO MAKE THE FINAL EXECUTABLES
$(EXE) : $(OFILES) Makefile
	$(COMPCPU) $(OPTCPU) $(OFILES) -o $(EXE) $(INCLUDECOM) $(LIBCPU)

$(EXECUDA) : $(OFILESCUDA) Makefile
	$(COMPCPU) $(OPTCPU_GPU) $(OFILESCUDA) -o $(EXECUDA) $(INCLUDECUDA) $(INCLUDECOM) $(LIBGPU)

#CLEANING RULES
clean :
	rm -f ./src/*.o

cleaner :
	rm -f ./src/*.o ./src/*~ *~ $(EXE) $(EXECUDA)
