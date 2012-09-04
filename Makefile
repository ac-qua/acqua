#################################
#
# Debug Options 
#
#  -g 
#  -check all
#
# should be removed for speed tests
#
#

FC      = ifort -g
FFLAGS  = -O  -check all
CFLAGS  = -O -DLINUX
LINKER  = ifort -g 
LIBS    = 
PROG    = ./acqua

#####################################################################
# manually add one object file for each source file (.f or .f90)
OBJS= prmat.o matmath.o pointgroups.o symmetry.o readtheinput.o main.o
#####################################################################
# compile a .o from each .f
%.o: %.f
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@

# compile a .o from each .f90
%.o: %.f90
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@

# compile a .mod from each .f
#%.mod: %.f
#	@echo "making $@ from $<"
#	$(FC) $(FFLAGS) -c $< -o $@
#
## compile a .mod from each .f90
#%.mod: %.f90
#	@echo "making $@ from $<"
#	$(FC) $(FFLAGS) -c $< -o $@

# link
$(PROG): $(OBJS)
	$(LINKER) $(OBJS) $(LIBS) -o $(PROG)

clean:
	rm -f *.o $(PROG)
