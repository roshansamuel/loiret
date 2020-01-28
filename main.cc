#include <iostream>

#include "parser.h"
#include "poisson.h"
#include "parallel.h"

/**
 ********************************************************************************************************************************************
 * \brief   Function to impose the Taylor Green Vortex initial condition on the given input velocity field
 *
 *          Depending on the preprocessor flag PLANAR, the function applies the Taylor Green vortex equation
 *          to set initial conditions in both 2D and 3D cases.
 *
 ********************************************************************************************************************************************
 */
void initializeField(plainsf &uField, grid &mesh) {
    if (mesh.rankData.rank == 0) std::cout << "Imposing Taylor-Green vortices initial condition" << std::endl << std::endl;

#ifdef PLANAR
    for (int i=uField.F.lbound(0); i <= uField.F.ubound(0); i++) {
        for (int k=uField.F.lbound(2); k <= uField.F.ubound(2); k++) {
            uField.F(i, 0, k) = sin(2.0*M_PI*mesh.xColloc(i)/mesh.xLen)*
                                cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
        }
    }
#else
    for (int i=uField.F.lbound(0); i <= uField.F.ubound(0); i++) {
        for (int j=uField.F.lbound(1); j <= uField.F.ubound(1); j++) {
            for (int k=uField.F.lbound(2); k <= uField.F.ubound(2); k++) {
                uField.F(i, j, k) = sin(2.0*M_PI*mesh.xColloc(i)/mesh.xLen)*
                                    cos(2.0*M_PI*mesh.yStaggr(j)/mesh.yLen)*
                                    cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
            }
        }
    }
#endif
}


int main() {
    // INITIALIZE MPI
    MPI_Init(NULL, NULL);

    // ALL PROCESSES READ THE INPUT PARAMETERS
    parser inputParams;

    // INITIALIZE PARALLELIZATION DATA
    parallel mpi(inputParams);

    // INITIALIZE GRID DATA
    grid gridData(inputParams, mpi);

    poisson *mgSolver;

    sfield P(gridData, "P");
    plainsf Pp(gridData, P);
    plainsf mgRHS(gridData, P);

    initializeField(Pp, gridData);

#ifdef PLANAR
    mgSolver = new multigrid_d2(gridData, inputParams);
#else
    mgSolver = new multigrid_d3(gridData, inputParams);
#endif

    mgSolver->mgSolve(Pp, mgRHS);

    // FINALIZE AND CLEAN-UP
    MPI_Finalize();

    return 0;
}
