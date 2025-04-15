### Gnina:

Gnina has the [bin file](https://github.com/gnina/gnina/releases/tag/v1.0.3). Download it, `chmod +x gnina`, and then specify its path in the `./docking/gnina.py` 

### LeDock

LeDock has the [bin file](http://www.lephar.com/download.htm). Download the LeDock and Lepro, `chmod +x lepro_linux_x86`; `chmod +x ledock_linux_x86`, and then specify its path in the `./docking/ledock.py` 

### Qvina 

Download the qvina [binary](https://github.com/QVina/qvina/tree/master), or `conda install -c conda-forge qvina `. 

If you cannot access the internet, just download the conda .tar.gz file and `conda install --use-local ./filepath`. 

For the prepare_target and prepare_receptor command (used to generate pdbqt files)

Download [ADFR](https://ccsb.scripps.edu/adfr/downloads/), then `chmod a+x ADFRsuite_Linux-x86_64_1.0_install` `./ADFRsuite_Linux-x86_64_1.0_install`. 

Do not `conda install mgltools`, this command will upgrade your environment to the python 2.X version, it ruins everything!!!

### rDock

I recommend use `conda install -c bioconda rdock`, you can also compile it from the source. 

### Smina 

Down load the [smina.static](https://sourceforge.net/projects/smina/files/smina.static/download), and change the corresponding path in the `./docking/smina.py`. 

### SurfLex

Download the sybyl software, like [this link](https://sybyl-x.software.informer.com/2.1/). And then specify the TA_LICENCE and binary path in the `./docking/surflex.py`. 