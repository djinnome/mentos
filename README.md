# mentos
Maximum Entropy production Of the Stoichiometric matrix


To install, first you need my forked version of pyOpt:

	pip install git+https://github.com/djinnome/pyopt.git

For this to successfully install, you will need to be in a Python2 environment. You may also need swig and a fortran compiler.  

If you are on a Mac, you can accomplish this with:

	brew install gfortran swig
	
Now you can install the rest of the dependencies this way:

	git clone https://github.com/djinnome/mentos.git
	cd mentos
	pip install --upgrade .
	
