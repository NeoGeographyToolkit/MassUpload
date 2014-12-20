
mkdir data
mkdir data/google
mkdir programs
mkdir repo


# Install Python 2.7.1
wget http://www.python.org/ftp/python/2.7.1/Python-2.7.1.tgz
tar -zxvf Python-2.7.1.tgz
cd Python-2.7.1
make clean
./configure --prefix=/home/smcmich1/programs/Python-2.7.1-install
make -j 8
make install

# Download a virtual environment
cd ~/programs
wget https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.7.tar.gz --no-check-certificate
tar xvfz virtualenv-1.7.tar.gz
cd virtualenv-1.7
#python virtualenv.py ve2
python virtualenv.py ve2.7.1 --python=/home/smcmich1/programs/Python-2.7.1-install/bin/python2.7
. ve2.7.1/bin/activate

# Install the easy python modules
pip install beautifulsoup4
pip install argparse
pip install requests
pip install hashlib
pip install google-api-python-client
pip install simplekml


# Build sqlite and install python wrapper
cd ~/programs
wget http://www.sqlite.org/2014/sqlite-autoconf-3080600.tar.gz
tar -xvf sqlite-autoconf-3080600.tar.gz
cd sqlite-autoconf-3080600/
./configure --prefix=/home/smcmich1/programs/sqlite-install
make -j 8
make install
CPPFLAGS="-I/home/smcmich1/programs/sqlite-install/include  -L/home/smcmich1/programs/sqlite-install/lib"  pip install pysqlite

## Download CMAKE
#wget http://www.cmake.org/files/v2.8/cmake-2.8.12.2.tar.gz
#tar -xvf cmake-2.8.12.2.tar.gz
#cd cmake-2.8.12.2
#./bootstrap
#make
#make install


# Install openjpeg
cd ~/repo
#svn checkout http://openjpeg.googlecode.com/svn/trunk/
wget http://openjpeg.googlecode.com/files/openjpeg-2.0.0.tar.gz
tar -xvf openjpeg-2.0.0.tar.gz
cd openjpeg-2.0.0
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/home/smcmich1/programs/openjpeg2.0.0-install
make -j 8
make install

# Build GDAL
cd ~/repo
wget http://download.osgeo.org/gdal/1.11.0/gdal-1.11.0.tar.gz
tar xvfz gdal-1.11.0.tar.gz
cd gdal-1.11.0
./configure --with-hdf5=no --with-openjpeg=/home/smcmich1/programs/openjpeg2.0.0-install --with-python --prefix=/home/smcmich1/programs/gdal-1.11.0-install
make -j 8
make install

# TODO: Set up gdal_edit.py!

# Download repositories
git config --global user.email scott.t.mcmichael@nasa.gov
git config --global user.name Scott McMichael
cd ~/repo
git clone https://github.com/ScottMcMichael/MassUpload.git
git clone https://github.com/NeoGeographyToolkit/Tools.git
cd MassUpload
g++ -o fix_jp2 fix_jp2.cpp # HiRISE only


#firefox --no-remote & # Need ready to authorize!

# TODO: Copy keys.txt from local installation!
# TODO: Set up bash script to init paths!
# TODO: Copy appropriate database file!




