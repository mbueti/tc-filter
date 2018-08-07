# tc-filter

To build on mac, run:
 `gfortran -fno-range-check -v -L /usr/local/opt/openblas/lib -I /usr/local/opt/openblas/include -lblas -I /usr/local/include -lnetcdf -lnetcdff -o filter-new.o filter-new.`

After building, invoke the filter program against a particular JRA-55do dataset by running:
`./filter-new.o <u-wind-file>.nc <u-wind-field-name> <v-wind-file>.nc <v-wind-field-name>`

For example:
`./filter-new.o u_10.2012.18Oct2017.nc uas_10m v_10.2012.18Oct2017.nc vas_10m`