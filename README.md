# RAPTOR

This is the public version of the RAPTOR code for radiative transfer in arbitrary spacetimes. 

<br />

# QUICK START

> These instructions guide the user in downloading and compiling the source code, and then producing an image-data file with RAPTOR, which can be plotted with the included Python plotting script.
Please note that these instructions expect the user to have a UNIX-like machine with a command terminal.

### Prerequisites / dependencies
-The gcc compiler 

-The OpenMP library

-The GNU Scientific computing Library (GSL)

-The CBLAS library

-Python 2.7 (only needed for plotting; a Python 3 compatible plotter script is coming soon)

### Downloading and compiling the RAPTOR code

Please download the repository, e.g., by typing:
```
git clone https://github.com/tbronzwaer/raptor.git
```
Next, it is necessary to compile the code using the included makefile. Please do so by typing:
```
make harm CPU=1
```
If everything goes correctly, RAPTOR will now have been compiled, and an executable named RAPTOR will have been created. 

### Running the code and plotting the result

Please create an 'output' directory for RAPTOR (from the directory in which the RAPTOR executable resides):
```
mkdir output
```
Finally, one can run the included run.sh batch script, which will instruct RAPTOR to create a single image-data file and a spectrum-data file in the output folder. In order to plot the output, please run the included plotter.py Python script (making sure to edit plotter.py first, so that it points to the correct location for the image-data file):
```
python plotter.py
```

<br />

# TROUBLESHOOTING & FAQ

**Q: RAPTOR failed to compile because the -fopenmp flag caused a compatibility issue.**

*A: A potential solution to this is to change the flag to -openmp instead. This ought to produce a fully functional executable, although I have not verified this. Thanks to user Saurabh for the feedback!*

<br />

# SUPPORT & CONTACT

Updates to RAPTOR, as well as the documentation, will soon be coming. If you have any questions, please contact me at t.bronzwaer[at]astro(dot)ru(dot)nl.
