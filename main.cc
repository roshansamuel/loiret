#include <iostream>

#include "parser.h"
#include "poisson.h"
#include "parallel.h"

/**
 ********************************************************************************************************************************************
 * \brief   Declaration of function to impose a sinusoidal variation for the input plain scalar field.
 *
 ********************************************************************************************************************************************
 */
void initializeField(plainsf &uField, grid &mesh);

int main() {
    // INITIALIZE MPI
    MPI_Init(NULL, NULL);

    // ALL PROCESSES READ THE INPUT PARAMETERS IN input/parameters.yaml FILE
    parser inputParams;

    // INITIALIZE PARALLELIZATION DATA
    parallel mpi(inputParams);

    // INITIALIZE GRID DATA
    grid gridData(inputParams, mpi);

    // POISSON SOLVER
    poisson *mgSolver;

    // TEMPLATE SCALAR FIELD (sfield) FOR PRESCRIBING PLAIN SCALAR FIELD (plainsf)
    sfield P(gridData, "P");

    // PLAIN SCALAR FIELDS THAT SERVE AS LHS AND RHS FOR THE POISSON SOLVER
    plainsf mgLHS(gridData, P);
    plainsf mgRHS(gridData, P);

    // INITIALIZE THE RHS WITH A SINUSOIDAL VARIATION
    initializeField(mgRHS, gridData);

#ifdef PLANAR
    mgSolver = new multigrid_d2(gridData, inputParams);
#else
    mgSolver = new multigrid_d3(gridData, inputParams);
#endif

    if (mpi.rank == 0 ) std::cout << "Value of LHS at (5, 5, 5) before solving: " << mgLHS.F(5, 5, 5) << std::endl;

    // SOLVE THE POISSON EQUATION
    mgSolver->mgSolve(mgLHS, mgRHS);

    if (mpi.rank == 0 ) std::cout << "Value of LHS at (5, 5, 5) after solving: " << mgLHS.F(5, 5, 5) << std::endl;

    // FINALIZE AND CLEAN-UP
    MPI_Finalize();

    return 0;
}


/**
 ********************************************************************************************************************************************
 * \brief   Definition of function to impose a sinusoidal variation for the input plain scalar field.
 *
 *          Depending on the preprocessor flag PLANAR, the function applies the sinusoidal variation
 *          to set initial conditions in both 2D and 3D cases.
 *
 ********************************************************************************************************************************************
 */
void initializeField(plainsf &uField, grid &mesh) {
    if (mesh.rankData.rank == 0) std::cout << "Imposing sinusoidal initial condition" << std::endl << std::endl;

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
