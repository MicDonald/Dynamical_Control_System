# mode 0 = uninstall
# mode 1 = install 

# type 0 = serial
# type 1 = ompi
# type 2 = mpi

mode=1
runtype=1

if (! test -e ../Makefile.package) then
  cp ../Makefile.package.empty ../Makefile.package
fi
#rm ../fix_VRTransition*
rm ../fix_TransientResponse*
#sed -i -e 's/-I\/opt\/intel\/Compiler\/11.0\/074\/mkl\/include //' ../Makefile.package
sed -i -e 's/-I..\/USER-DYNAMICAL_CTRLSYSTEM //' ../Makefile.package
sed -i -e 's/-L..\/USER-DYNAMICAL_CTRLSYSTEM //' ../Makefile.package
sed -i -e 's/-lDynamicalCtrlSystem //' ../Makefile.package

if (test $mode = 1) then
#  sed -i -e 's/^PKG_INC =[ \t]*/&-I\/opt\/intel\/Compiler\/11.0\/074\/mkl\/include /' ../Makefile.package
  sed -i -e 's/^PKG_INC =[ \t]*/&-I..\/USER-DYNAMICAL_CTRLSYSTEM /' ../Makefile.package
  sed -i -e 's/^PKG_PATH =[ \t]*/&-L..\/USER-DYNAMICAL_CTRLSYSTEM /' ../Makefile.package
  sed -i -e 's/^PKG_LIB =[ \t]*/&-lDynamicalCtrlSystem /' ../Makefile.package

  #if (test $runtype = 2) then
   #sed -i -e 's/^PKG_INC =[ \t]*/&-I..\/opt/intel/Compiler/11.0/074/mkl/include /' ../Makefile.package
   #make MPI_LIB
   # cd ..
   # make mpi -j 8
   # cd USER-VRTRANSITION
  if (test $runtype = 1) then
    cp fix_* ..
    cd ..
    make omp -j 8
    cd _DCS
    #make OpenMP_LIB
  elif (test $runtype = 0) then 
    make LIB
    cd ..
    rm fix_DynamicalCtrlSystem.*
    make serial -j 8
    cd _DCS
  fi
fi

