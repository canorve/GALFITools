
## initgal.py

It creates files containing Galfit's initial
parameters to run it from a galfit file. 

Example of how to use it:

```
 ./initgal.py [GalfitFile] [-n 3] [-p4 1 8 ] 
```

The example above will create 3 files based on GalfitFile with different
values for the Sersic parameter in the range from 1 to 8.

The program will create the galfit input files and a file named
rungalfit.sh to run galfit with the galfit files from bash.


