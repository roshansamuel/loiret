#ifndef PLAINSF_H
#define PLAINSF_H

#include "sfield.h"
#include "grid.h"

class plainsf {
    private:
        blitz::firstIndex i;
        blitz::secondIndex j;
        blitz::thirdIndex k;

        const grid &gridData;

    public:
        blitz::Array<real, 3> F;

        blitz::Range xColl, yColl, zColl;

        plainsf(const grid &gridData, const sfield &refF);

        mpidata *mpiHandle;

        plainsf& operator += (plainsf &a);
        plainsf& operator -= (plainsf &a);

        plainsf& operator += (sfield &a);
        plainsf& operator -= (sfield &a);

        plainsf& operator *= (real a);

        void operator = (plainsf &a);
        void operator = (sfield &a);
        void operator = (real a);

/**
 ********************************************************************************************************************************************
 * \brief   Function to synchronise data across all processors when performing parallel computations
 *
 *          This function calls the \ref mpidata#syncData "syncData" function of mpidata class to perform perform data-transfer and thus update
 *          the sub-domain boundary pads.
 ********************************************************************************************************************************************
 */
        inline void syncData() {
            mpiHandle->syncData();
        }

/**
 ********************************************************************************************************************************************
 * \brief   Function to extract the maximum value from the plain scalar field
 *
 *          The function uses the in-built blitz function to obtain the maximum value in an array.
 *          While performing parallel computation, the function performs an <B>MPI_Allreduce()</B> to get
 *          the global maximum from the entire computational domain.
 *
 * \return  The real value of the maximum is returned (it is implicitly assumed that only real values are used)
 ********************************************************************************************************************************************
 */
        inline real fxMax() {
            real localMax, globalMax;

            localMax = blitz::max(F);

            MPI_Allreduce(&localMax, &globalMax, 1, MPI_FP_REAL, MPI_MAX, MPI_COMM_WORLD);

            return globalMax;
        }

        ~plainsf() { };
};

/**
 ********************************************************************************************************************************************
 *  \class plainsf plainsf.h "lib/plainsf.h"
 *  \brief Plain scalar field class to store simple scalar fields with no differentiation or interpolation
 *
 *  The class stores scalar fields in the form of a Blitz array
 ********************************************************************************************************************************************
 */

#endif
