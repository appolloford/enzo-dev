#default fortran compiler
fc = ifort
cc = icc
clib = -nofor_main

#executable name
exec = test

#default libraries
lib = -llapack

#guess if you have MKL
GREP_MKL = $(shell export | grep 'mkl')
ifneq ($(GREP_MKL),)
	lib = -mkl
endif

#if dgesv is not present non need for LAPACK
GREP = $(shell grep -i 'dgesv' krome_subs.f90)
ifeq ($(GREP),)
	lib =
endif

#flags
switchOPT = -O3 -ipo -ip -unroll -xHost -g -fp-model precise
switchDBG = -O0 -check all -warn all -fpe0 -u -traceback -warn nounused
switchDBG += -fp-model precise
switchPRO = $(switchOPT) -pg -traceback
switchOMP = $(switchOPT) -openmp

#default switch
switch = $(switchOPT)
switchc = $(switchOPT)

#no warning switch (for DLSODES)
nowarnIfort = -nowarn
nowarnGNU = -w

#objects
objs = opkda2.o
objs += opkda1.o
objs += opkdmain.o
objs += krome_commons.o
objs += krome_constants.o
objs += krome_fit.o
objs += krome_getphys.o
objs += krome_gadiab.o
objs += krome_grfuncs.o
objs += krome_phfuncs.o
objs += krome_user_commons.o
objs += krome_subs.o
objs += krome_stars.o
objs += krome_dust.o
objs += krome_photo.o
objs += krome_tabs.o
objs += krome_coolingGH.o
objs += krome_cooling.o
objs += krome_heating.o
objs += krome_ode.o
objs += krome_user.o
objs += krome_reduction.o
objs += krome.o

cobjs = krome_header.o

#default target
all: 	$(objs) test.o
	$(fc) $(objs) test.o -o $(exec) $(switch) $(lib)

#ifort full debug target
debug: switch = $(switchDBG)
debug: switchc = -traceback
debug: nowarn = $(nowarnIfort)
debug: all

#ifort profile target
profile: switch = $(switchPRO)
profile: nowarn = $(nowarnIfort)
profile: all

#gfortran target
gfortran: fc = gfortran
gfortran: switch = -ffree-line-length-none
gfortran: nowarn = $(nowarnGNU)
gfortran: all

#gfortran debug target
gfortran_dbg: fc = gfortran
gfortran_dbg: switch = -fbacktrace -g
gfortran_dbg: switch += -ffpe-trap=zero,overflow,invalid
gfortran_dbg: switch += -fbounds-check -ffree-line-length-none
gfortran_dbg: nowarn = $(nowarnGNU)
gfortran_dbg: all

#library accessible from python interface target
sharedlib:	$(objs)
	$(fc) $(objs) -o libkrome.so $(switch) -shared $(lib)

pyinterface: fc = gfortran
pyinterface: switch = -ffree-line-length-none -fPIC
pyinterface: nowarn = $(nowarnGNU)
pyinterface: sharedlib

pyinterface_opt: fc = gfortran
pyinterface_opt: switch = -O3 -ffree-line-length-none -fPIC
pyinterface_opt: nowarn = $(nowarnGNU)
pyinterface_opt: sharedlib

cinterface: $(cobjs)
cinterface: objs += $(cobjs)
cinterface: cc = gcc
cinterface: switchc =
cinterface: clib = -lm -lgfortran
cinterface: gfortran

cinterface_opt: $(cobjs)
cinterface_opt: fc = gfortran
cinterface_opt: cc = gcc
cinterface_opt: switch = -O3 -ffree-line-length-none
cinterface_opt: switchc = -O3
cinterface_opt: nowarn = $(nowarnGNU)
cinterface_opt: clib = -lm -lgfortran
cinterface_opt: all

#clean target
clean:
	rm -f *.o *.mod *__genmod.f90 *~ $(exec)

#rule for f90
%.o:%.f90
	$(fc) $(switch) -c $^ -o $@

#rule for f
%.o:%.f
	$(fc) $(switch) $(nowarn) -c $^ -o $@

#rule for c
%.o:%.c
	$(cc) $(switchc) $(nowarn) -c $^ -o $@
