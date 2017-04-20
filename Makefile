all:
	rm -f *.so
	LD_LIBRARY_PATH=./libs/_install/lib:/home/julien/tm/v2/libs/_install/lib:$LD_LIBRARY_PATH python setup.py build_ext --inplace
	rm -rf build/
