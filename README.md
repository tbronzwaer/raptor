# RAPTOR

This is the public version of the RAPTOR code for radiative transfer in arbitrary spacetimes. 

# QUICK START

> These instructions guide the user in downloading and compiling the source code, and then producing an image-datafile with RAPTOR, which can be plotted with the included Python plotting script.
Please note that these instructions expect the user to have a UNIX-like machine with a command terminal, and that the gcc compiler, OpenMP, and Python are installed.

Please download the repository, e.g., by typing:
```
git clone https://github.com/tbronzwaer/raptor.git
```
Next, it is necessary to compile the code using the included makefile. Please do so by typing:
```
make harm CPU=1
```
If everything goes correctly, RAPTOR will now have been compiled, and an executable named RAPTOR will have been created. Please create an 'output' directory for RAPTOR (from the directory in which the RAPTOR executable resides):
```
mkdir output
```
Finally, one can run the included run.sh batch script, which will instruct RAPTOR to create a single image-data file and a spectrum-data file in the output folder. In order to plot the output, please run the included plotter.py Python script (making sure to edit plotter.py first, so that it points to the correct location for the image-data file):
```
python plotter.py
```

# QUESTIONS?

Updates to RAPTOR, as well as the documentation, will soon be coming. If you have any questions, please contact me at t.bronzwaer[at]astro(dot)ru(dot)nl.
