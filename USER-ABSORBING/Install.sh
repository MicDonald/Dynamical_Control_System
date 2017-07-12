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
rm ../fix_absorbing*
sed -i -e 's/-I..\/USER-ABSORBING //' ../Makefile.package
sed -i -e 's/-L..\/USER-ABSORBING //' ../Makefile.package
sed -i -e 's/-lAbsorbing //' ../Makefile.package

if (test $mode = 1) then
  cp fix_absorbing.* ..
  sed -i -e 's/^PKG_INC =[ \t]*/&-I..\/USER-ABSORBING /' ../Makefile.package
  sed -i -e 's/^PKG_PATH =[ \t]*/&-L..\/USER-ABSORBING /' ../Makefile.package
  sed -i -e 's/^PKG_LIB =[ \t]*/&-lAbsorbing /' ../Makefile.package

  if (test $runtype = 2) then
    cp fix_absorbingMPI.* ..
    make MPI_LIB
  elif (test $runtype = 1) then
    cp fix_absorbingOMP.* ..
    make OpenMP_LIB
  elif (test $runtype = 0) then
    make LIB
  fi
fi

