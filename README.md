# CNMSSM_scan

This code interfaces FlexibleSUSY, lilith and MultiNest for a specific project.

The code depends on FlexibleSUSY, lilith and MultiNest.

************************** Initial setup *************************************

If you have not already done so must install FlexibleSUSY and SARAH yourself.

Then change the hardoceded path for FlexibleSUSY in the makefile.

To automatically download lilith and MultiNest

$ make init

Or if you want to build it at the same time instead do

$ make start

************************* Building code *************************************

To build the code

$ make all

To remove all files created during compilation do:

$ make clean

To also remove the downloaded pachkages MultiNest and lilith do:

$ make distclean


