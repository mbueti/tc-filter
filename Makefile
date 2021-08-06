FC = gfortran
#LDFLAGS = -L /usr/local/opt/openblas/lib -L /usr/local/lib
LDFLAGS = -L /usr/local/lib
#CPPFLAGS = -I /usr/local/opt/openblas/include -I /usr/local/include
CPPFLAGS = -I /usr/local/include
LOADLIBS = -lblas -lnetcdff -lnetcdf
#OPTIONS = -mmacosx-version-min=11.3 -static-libgcc -fno-range-check 
# OPTIONS = -fno-underscoring -mmacosx-version-min=11.3 -fno-range-check 
OPTIONS = -fcheck='all' -fno-range-check -ffpe-trap=zero,overflow,underflow -Og -g -fbacktrace

target = filter
object = filter.f90

filter: filter.f90

#default: $(target)

$(target): $(object)
	$(FC) $(OPTIONS) $(LDFLAGS) $(CPPFLAGS) $(LOADLIBS) -c filter_routines.f90 
	$(FC) $(OPTIONS) $(LDFLAGS) $(CPPFLAGS) $(LOADLIBS) -c datetime_module.f90
	$(FC) $(OPTIONS) $(LDFLAGS) $(CPPFLAGS) $(LOADLIBS) -o $@ $^ datetime_module.o filter_routines.o

clean:
	rm -rf datetime_module.mod
	rm -rf filter_routines.mod
	rm -rf datetime_module.o
	rm -rf filter_routines.o
	rm -f $(target)
