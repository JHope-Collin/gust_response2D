##########################################


FSRC = main.f95 \
	steady.f95 \
	unsteady.f95 \
	poisson.f95 \
	lift.f95 \
	field.f95

FMOD = constantsmod.f95 \
	optionsmod.f95 \
	callsmod.f95


###########################################


FOBJ = $(FSRC:.f95=.o)

MOBJ = $(FMOD:.f95=.o)
MODS = $(FMOD:mod.f95=.mod)

FCMP = gfortran
FOPT = -g -Wall -Wextra -fbounds-check

LIBS = -llapack -lblas


###########################################


a.out : $(FOBJ) $(MOBJ)
	$(FCMP) $(FOPT) -o $@ $^ $(LIBS)


###########################################


%.o : %.f95
	$(FCMP) $(FOPT) -o $@ -c $<

%.mod : %mod.f95
	$(FCMP) $(FOPT) -c $<

$(FOBJ) : $(MOBJ)

unsteady.o : poisson.o


###########################################


.PHONY : clean

clean :
	rm -f a.out $(FOBJ) $(MODS) $(MOBJ)
