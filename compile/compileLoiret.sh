# USER SET PARAMETERS - COMMENT/UNCOMMENT AS NECESSARY

PROC=4
REAL_TYPE="DOUBLE"
#REAL_TYPE="SINGLE"
#PLANAR="PLANAR"
#TIME_RUN="TIME_RUN"
EXECUTE_AFTER_COMPILE="EXECUTE"

# NO USER MODIFICATIONS NECESSARY BELOW THIS LINE

# REMOVE PRE-EXISTING EXECUTATBLES
rm -f ../loiret

# IF build DIRECTORY DOESN'T EXIST, CREATE IT
if [ ! -d build ]; then
    mkdir build
fi

# SWITCH TO build DIRECTORY
cd build

# RUN Cmake WITH NECESSARY FLAGS AS SET BY USER
if [ -z $PLANAR ]; then
    if [ -z $TIME_RUN ]; then
        if [ "$REAL_TYPE" == "DOUBLE" ]; then
            CC=mpicc CXX=mpicxx cmake ../../ -DREAL_DOUBLE=ON
        else
            CC=mpicc CXX=mpicxx cmake ../../ -DREAL_SINGLE=ON
        fi
    else
        CC=mpicc CXX=mpicxx cmake ../../ -DTIME_RUN=ON
    fi
else
    if [ -z $TIME_RUN ]; then
        CC=mpicc CXX=mpicxx cmake ../../ -DPLANAR=ON
    else
        CC=mpicc CXX=mpicxx cmake ../../ -DPLANAR=ON -DTIME_RUN=ON
    fi
fi

# COMPILE
make -j16

# SWITCH TO PARENT DIRECTORY
cd ../../

# RUN CODE IF REQUESTED BY USER
if ! [ -z $EXECUTE_AFTER_COMPILE ]; then
    mpirun -np $PROC ./loiret
fi
