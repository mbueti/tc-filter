FC = gfortran
#LDFLAGS = -L /usr/local/opt/openblas/lib -L /usr/local/lib
LDFLAGS = -L /usr/local/lib
#CPPFLAGS = -I /usr/local/opt/openblas/include -I /usr/local/include
CPPFLAGS = -I /usr/local/include
LOADLIBS = -lblas -lnetcdff -lnetcdf
#OPTIONS = -mmacosx-version-min=11.3 -static-libgcc -fno-range-check 
# OPTIONS = -fno-underscoring -mmacosx-version-min=11.3 -fno-range-check 
OPTIONS = -fcheck='all' -fno-range-check -ffpe-trap=zero,overflow,underflow -Og -g -fbacktrace

target = filter-new
object = filter-new.f

filter-new: filter-new.f

#default: $(target)

$(target): $(object)
	$(FC) $(OPTIONS) $(LDFLAGS) $(CPPFLAGS) $(LOADLIBS) -o $@ datetime_module.f90 $^ 

clean:
	rm -rf datetime_module.mod
	rm -f $(target)
