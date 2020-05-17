******** HOW TO COMPILE AND RUN EXAMPLES IN STOPT


----------------------------------------
		COMPILE IPOPT
-----------------------------------------

1-> cd $HOME/path/to/STOPT (The same level as this .txt)
2-> wget --continue http://www.coin-or.org/download/source/Ipopt/Ipopt-3.12.12.tgz
3-> tar xzvf Ipopt-3.12.12.tgz
4-> cd Ipopt-3.12.12/ThirdParty/Metis
5-> ./get.Metis
6-> cd ../Mumps
7-> ./get.Mumps
8-> cd ../Lapack
9-> ./get.Lapack
10-> cd ../Blas
11-> ./get.Blas
12-> cd ../ASL
13-> ./get.ASL
14-> cd $HOME/path/to/STOPT/Ipopt-3.12.12/
15->./configure
16->make
17->make install

-----------------------------------------
	ADD IPOPT TO STOPT
-----------------------------------------
18->cp -r include ../STOPT/ipopt/include
19->cp -r lib ../STOPT/ipopt/lib

-----------------------------------------
COMPILE STOPT AND EXECUTE A SIMPLE TEST
-----------------------------------------

20->cd ../STOPT/OCsolver

21->mkdir build
22->cd build
23->cmake .. && make

24->./test




