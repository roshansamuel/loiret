#ifndef PLAINVF_H
#define PLAINVF_H

#include "vfield.h"
#include "grid.h"

class plainvf {
    private:
        const grid &gridData;

    public:
        blitz::Array<real, 3> Vx, Vy, Vz;

        plainvf(const grid &gridData, const vfield &refV);

        mpidata *mpiVxData, *mpiVyData, *mpiVzData;

        plainvf& operator += (plainvf &a);
        plainvf& operator -= (plainvf &a);

        plainvf& operator += (vfield &a);
        plainvf& operator -= (vfield &a);

        plainvf& operator *= (real a);

        void operator = (plainvf &a);
        void operator = (vfield &a);
        void operator = (real a);

/**
 ********************************************************************************************************************************************
 * \brief   Function to synchronise data across all processors when performing parallel computations
 *
 *          Each of the individual field components have to send and receive data across its MPI decomposed sub-domains.
 *          This function calls the \ref mpidata#syncData "syncData" function for each component to update the sub-domain boundary pads.
 ********************************************************************************************************************************************
 */
        inline void syncData() {
            mpiVxData->syncData();
            mpiVyData->syncData();
            mpiVzData->syncData();
        }

/**
 ********************************************************************************************************************************************
 * \brief   Function to extract the maximum value from the Vx component of the plain vector field
 *
 *          The function uses the in-built blitz function to obtain the maximum value in an array.
 *          While performing parallel computation, the function performs an <B>MPI_Allreduce()</B> to get
 *          the global maximum from the entire computational domain.
 *
 * \return  The real value of the maximum is returned (it is implicitly assumed that only real values are used)
 ********************************************************************************************************************************************
 */
        inline real vxMax() {
            real localMax, globalMax;

            localMax = blitz::max(Vx);

            MPI_Allreduce(&localMax, &globalMax, 1, MPI_FP_REAL, MPI_MAX, MPI_COMM_WORLD);

            return globalMax;
        }

/**
 ********************************************************************************************************************************************
 * \brief   Function to extract the maximum value from the Vy component of the plain vector field
 *
 *          The function uses the in-built blitz function to obtain the maximum value in an array.
 *          While performing parallel computation, the function performs an <B>MPI_Allreduce()</B> to get
 *          the global maximum from the entire computational domain.
 *
 * \return  The real value of the maximum is returned (it is implicitly assumed that only real values are used)
 ********************************************************************************************************************************************
 */
        inline real vyMax() {
            real localMax, globalMax;

            localMax = blitz::max(Vy);

            MPI_Allreduce(&localMax, &globalMax, 1, MPI_FP_REAL, MPI_MAX, MPI_COMM_WORLD);

            return globalMax;
        }

/**
 ********************************************************************************************************************************************
 * \brief   Function to extract the maximum value from the Vz component of the plain vector field
 *
 *          The function uses the in-built blitz function to obtain the maximum value in an array.
 *          While performing parallel computation, the function performs an <B>MPI_Allreduce()</B> to get
 *          the global maximum from the entire computational domain.
 *
 * \return  The real value of the maximum is returned (it is implicitly assumed that only real values are used)
 ********************************************************************************************************************************************
 */
        inline real vzMax() {
            real localMax, globalMax;

            localMax = blitz::max(Vz);

            MPI_Allreduce(&localMax, &globalMax, 1, MPI_FP_REAL, MPI_MAX, MPI_COMM_WORLD);

            return globalMax;
        }

        ~plainvf() { };
};

/**
 ********************************************************************************************************************************************
 *  \class plainvf plainvf.h "lib/plainvf.h"
 *  \brief Plain vector field class to store simple vector fields with no differentiation or interpolation
 *
 *  The class stores vector fields in the form of three Blitz arrays
 *  The vector field is stored in such a way that the components are face-centered scalar fields, with:
 *      - x-component located at the face centers along the yz-plane
 *      - y-component located at the face centers along the zx-plane
 *      - z-component located at the face centers along the xy-plane
 *
 ********************************************************************************************************************************************
 */

#endif
