#ifndef PARSER_H
#define PARSER_H

#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
#include <blitz/array.h>
#include <yaml-cpp/yaml.h>

#ifdef REAL_DOUBLE
#define H5T_NATIVE_REAL H5T_NATIVE_DOUBLE
#define MPI_FP_REAL MPI_DOUBLE
#define real double
#else
#define H5T_NATIVE_REAL H5T_NATIVE_FLOAT
#define MPI_FP_REAL MPI_FLOAT
#define real float
#endif

class parser {
    public:
        int nThreads;
        int npY, npX;
        int xInd, yInd, zInd;
        int vcDepth, vcCount;
        int preSmooth, postSmooth;

        int xGrid, yGrid, zGrid;

        bool xPer, yPer, zPer;

        real tolerance;
        real Lx, Ly, Lz;
        real betaX, betaY, betaZ;

        std::string meshType;

        std::vector<int> interSmooth;

        parser();

        void writeParams();

    private:
        std::string domainType;

        void parseYAML();
        void checkData();

        void setGrids();
        void setPeriodicity();
};

/**
 ********************************************************************************************************************************************
 *  \class parser parser.h "lib/io/parser.h"
 *  \brief  Contains all the global variables set by the user through the yaml file
 *
 *  The class parses the paramters.yaml file and stores all the simulation paramters in publicly accessible constants.
 *  The class also has a function to check the consistency of the user set paramters and throw exceptions.
 *  The class is best initialized as a constant to prevent inadvertent tampering of the global variables it contains.
 ********************************************************************************************************************************************
 */

#endif
