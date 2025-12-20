# Package of Numerical Weather Prediction model MOLOCH

MOLOCH model vers. July 2022 by P Malguzzi and implemented at LaMMA

0. Compile eccodes and radiation ECMWF (2012). Packages are not provided here.

1. Compile geo.F90 with informations from geo.inc and static files (to be downloaded separately)
gfortran -o geo.exe geo.F90

2. Set dimensions.inc accordingly

3. Then run:
   3.1 compila-premoloch
   3.2 compila-moloch
   3.3 compila-ppostmol
   3.4 compila-shf2grib2

For more information write to capecchi and pasi @ lamma toscana it
