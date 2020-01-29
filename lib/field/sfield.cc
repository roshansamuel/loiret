#include "plainsf.h"
#include "sfield.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the sfield class
 *
 *          The instance of the field class to store the data of the scalar field is initialized, and the necessary grid
 *          transformation derivatives along each direction are chosen according to the grid staggering.
 *          The arrays to store the output from various operators like derivatives, convective derivatives, etc. are also
 *          allocated.
 *          Finally, an instance of the <B>mpidata</B> class is initialized to store the sub-arrays to be send/received
 *          across the processors during MPI communication.
 *
 * \param   gridData is a const reference to the global data contained in the grid class
 * \param   fieldName is a string value set by the user to name and identify the scalar field
 ********************************************************************************************************************************************
 */
sfield::sfield(const grid &gridData, std::string fieldName):
               gridData(gridData),
               F(gridData, fieldName, true, true, true),
               derS(gridData, F)
{
    this->fieldName = fieldName;

    interTempF.resize(F.fSize);
    interTempF.reindexSelf(F.flBound);

    derivTempF.resize(F.fSize);
    derivTempF.reindexSelf(F.flBound);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to synchronise data across all processors when performing parallel computations
 *
 *          This function calls the \ref mpidata#syncData "syncData" function of mpidata class to perform perform data-transfer and thus update
 *          the sub-domain boundary pads.
 ********************************************************************************************************************************************
 */
void sfield::syncData() {
    F.syncData();
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to add a given plain scalar field
 *
 *          The unary operator += adds a given plain scalar field to the entire field stored as sfield and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to a plainsf to be added to the member field
 *
 * \return  A pointer to itself is returned by the scalar field class to which the operator belongs
 ********************************************************************************************************************************************
 */
sfield& sfield::operator += (plainsf &a) {
    F.F += a.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to subtract a given plain scalar field
 *
 *          The unary operator -= subtracts a given plain scalar field from the entire field stored as sfield and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to a plainsf to be deducted from the member field
 *
 * \return  A pointer to itself is returned by the scalar field class to which the operator belongs
 ********************************************************************************************************************************************
 */
sfield& sfield::operator -= (plainsf &a) {
    F.F -= a.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to add a given scalar field
 *
 *          The unary operator += adds a given scalar field to the entire field stored as sfield and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another sfield to be added to the member field
 *
 * \return  A pointer to itself is returned by the scalar field class to which the operator belongs
 ********************************************************************************************************************************************
 */
sfield& sfield::operator += (sfield &a) {
    F.F += a.F.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to subtract a given scalar field
 *
 *          The unary operator -= subtracts a given scalar field from the entire field stored as sfield and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another sfield to be deducted from the member field
 *
 * \return  A pointer to itself is returned by the scalar field class to which the operator belongs
 ********************************************************************************************************************************************
 */
sfield& sfield::operator -= (sfield &a) {
    F.F -= a.F.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to multiply a scalar value to the scalar field
 *
 *          The unary operator *= multiplies a real value to the entire field stored as sfield and returns
 *          a pointer to itself.
 *
 * \param   a is a real number to be multiplied to the scalar field
 *
 * \return  A pointer to itself is returned by the scalar field class to which the operator belongs
 ********************************************************************************************************************************************
 */
sfield& sfield::operator *= (real a) {
    F.F *= a;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign a plain scalar field to the scalar field
 *
 *          The operator = copies the contents of the input plain scalar field to itself.
 *
 * \param   a is the plainsf to be assigned to the scalar field
 ********************************************************************************************************************************************
 */
void sfield::operator = (plainsf &a) {
    F.F = a.F;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign a scalar field to the scalar field
 *
 *          The operator = copies the contents of the input scalar field to itself.
 *
 * \param   a is the scalar field to be assigned to the scalar field
 ********************************************************************************************************************************************
 */
void sfield::operator = (sfield &a) {
    F.F = a.F.F;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign a scalar value to the scalar field
 *
 *          The operator = assigns a real value to all the scalar field.
 *
 * \param   a is a real number to be assigned to the scalar field
 ********************************************************************************************************************************************
 */
void sfield::operator = (real a) {
    F.F = a;
}
