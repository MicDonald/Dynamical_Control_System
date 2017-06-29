# mode 0 = uninstall
# mode 1 = install 

# type 0 = serial
# mode 1 = ompi
# mode 2 = mpi

mode=${1:-1}
runtype=${2:-2}

if (! test -e ../Makefile.package) then
  cp ../Makefile.package.empty ../Makefile.package
fi
rm ../fix_VRTransition*
sed -i -e 's/-I..\/USER-VRTRANSITION //' ../Makefile.package
sed -i -e 's/-L..\/USER-VRTRANSITION//' ../Makefile.package
sed -i -e 's/-lVRTransition //' ../Makefile.package

if (test $mode = 1) then
  cp fix_VRTransition.* ..
  sed -i -e 's/^PKG_INC =[ \t]*/&-I..\/USER-VRTRANSITION /' ../Makefile.package
  sed -i -e 's/^PKG_PATH =[ \t]*/&-L..\/USER-VRTRANSITION /' ../Makefile.package
  sed -i -e 's/^PKG_LIB =[ \t]*/&-lVRTransition /' ../Makefile.package

  if (test $runtype = 2) then
    cp fix_VRTransition.* ..
    make MPI_LIB
    cd ..
    make mpi -j 8  
#elif (test $runtype = 1) then
  #  cp fix_VRTransitionOMP.* ..
  #  make OpenMP_LIB
  if (test $runtype = 0) then
    make LIB
    cd ..
    make serial -j 8
  fi
fi
