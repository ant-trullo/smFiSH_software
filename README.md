This software was developed in Python™ 3.6.

- Packages needed for smFiSH_software to work:

	numpy 1.19.1 , pyqtgraph 0.11.0 , skimage 0.17.2 , sklearn 0.22.1 , scipy 1.5.2 , PyQt5 5.13.1 , xlwt 1.3.0, xlrd 1.2.0 , tifffile 2020.9.3 ,
        czifile 2019.7.2 , cython 0.26.1
	
	The number refers to the version we have used. smFiSH_software may work anyway 
	with newer version of the same packages, unless there are changes in syntax.
	All the dependecies of these packages must be fulfilled: you will be 
        required to install matplotlib, PIL and some others. Depending on the installation 
        technique, the Python installer can take care of this directly, or you will have to do this 
        manually. 

- Requirements
	
	In order to run the cython function, a C++ compiler is needed. For linux users, you need to install the gcc, 
	for Windows users instead you need to install the visualstudio compiler (you can follow this guide
	https://wiki.python.org/moin/WindowsCompilers )
	smFiSH_software was developed on the operative system Linux Ubuntu 18.04 64-bit 
        and tested on Linux Ubuntu 18.04 64-bit, Windows 8.1 Professional, Windows 10 Pro.

 
- Install smFiSH_software:

	Unzip the file smFiSH_software.zip and put all the files and the subfolder in a folder.
	In order to compile the cython files, open a terminal, or a cmd in case of Windows, move into the folder where .py and .pyx files are and then, one by         one, run the following commands:
	
	python3 setup_MNU.py build_ext --inplace
	
	python3 setup_NPU.py build_ext --inplace
	
	python3 setup_SBSUS.py build_ext --inplace

	
- Run smFiSH_software:	
	
	If the cython files are correctly compiled, open a terminal (cmd) in the folder with the .py, .pyx and the compiled files and run: 
	
	python3 smFish_GUI_v3_1.py
	
	The graphical user interface will pop up and let you work.
    	![me](https://github.com/ant-trullo/smFiSH_software/blob/main/smFiSH-software.gif)
    
    
           
For any question or issue send an email at:
    antonio.trullo@igmm.cnrs.fr            


