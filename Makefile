FC = gfortran
LDFLAGS = -L /usr/local/opt/openblas/lib
CPPFLAGS = -I /usr/local/opt/openblas/include -I /usr/local/include
LOADLIBES = -lblas -lnetcdf -lnetcdff -ldatetime
OPTIONS = -fno-range-check

target = filter-new
object = filter-new.f

default: $(target)

$(target): $(object)
	$(FC) $(OPTIONS) $(LDFLAGS) $(CPPFLAGS) $(LOADLIBES) -o $@ $^

clean:
	rm -f $(target)