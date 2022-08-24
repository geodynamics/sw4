FROM ubuntu:latest

WORKDIR /home

# build environment
RUN apt-get update && \
  DEBIAN_FRONTEND='noninteractive' \
  DEBCONF_NONINTERACTIVE_SEEN='true' \
  apt-get install --yes wget build-essential git cmake gfortran openmpi-bin && \
\
\
# install hdf5 1.12.1, and remain packages \
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.1/src/hdf5-1.12.1.tar.gz && \
  tar xzf hdf5-1.12.1.tar.gz && rm hdf5-1.12.1.tar.gz && cd hdf5-1.12.1 && ./configure --prefix=/usr/local/hdf5 --enable-parallel --enable-build-mode=production && \
  make -j4 && make install && cd .. && \
  apt-get install --yes libproj-dev libhdf5-mpi-dev libfftw3-dev libfftw3-mpi-dev \
    liblapack-dev python3-h5py libcurl3-dev && \
\
\
# install zfp and H5Z-ZFP \
git -c http.sslVerify=false clone https://github.com/LLNL/zfp.git && cd zfp/ && \
  sed -i 's/# DEFS += -DBIT_STREAM_WORD_TYPE=uint8/DEFS += -DBIT_STREAM_WORD_TYPE=uint8/g' Config && \
  make && cd .. && \
git -c http.sslVerify=false clone https://github.com/LLNL/H5Z-ZFP.git && cd H5Z-ZFP/ && \
  make FC= CC=mpicc ZFP_HOME=/home/zfp HDF5_HOME=/usr/local/hdf5 all install && \
  cd ..

# build sw4 
RUN git -c http.sslVerify=false clone https://github.com/geodynamics/sw4.git && \
  cd sw4/ && git switch developer && rm -r .git/ && \
  ln -s /usr/lib/x86_64-linux-gnu/libssl.so.1.1  /usr/lib64/libssl.so && \
  ln -s /usr/lib/x86_64-linux-gnu/libcrypto.so.1.1  /usr/lib64/libcrypto.so && \
  make sw4 CXX=mpicxx FC=gfortran debug=no proj_6=yes hdf5=yes fftw=yes zfp=yes prec=double EXTRA_LINK_FLAGS="-L/usr/lib64 -lgfortran -lhdf5 -llapack"  HDF5ROOT=/usr/local/hdf5 FFTWHOME=/usr/lib/x86_64-linux-gnu ZFPROOT=/home/zfp H5ZROOT=/home/H5Z-ZFP/install -j4

ENV PATH="$PATH:/home/sw4/optimize_mp" \
  LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/hdf5/lib"

# run sw4
# I typically use:
# docker run --cap-add=SYS_PTRACE -it -v $(pwd):/home/model/ sw4:zfp
# to run the docker container in avoid of the CMA issue (https://github.com/open-mpi/ompi/issues/4948)
