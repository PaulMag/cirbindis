# CirBinDis #

Version 0.2.1

Software for analysing simulated density maps of circumstellar/circumbinary disks of gas.

### What does it do? ###

* Create artificial lightcurves from a full rotation of the disk and stars for various inclinations.
* Create density profiles from various line-of-sights.
* Save results as csv-files and plots.

### What do I need? ###

* python
* numpy
* scipy
* matplotlib
* astropy

### How do I get it? ###

Clone this repository or download the zip-archive. Remember to update.

### How do I use it? ###

* Have you dataset with <x,y,ρ>- or <r,θ,ρ>-columns.
* Copy and configure 'xml/input.xml' with input parameters.
* Run 'src/cirbindis.py' (make an alias for this).
```
python cirbindis.py myinputfile.xml
```
Read the usermanual for more detailed instructions: 'doc/cirbindis_usermanual.pdf'

### Who do I talk to? ###

For questions, troubleshooting, suggestions, and/or feedback, contact the developer: paulmag91@gmail.com
