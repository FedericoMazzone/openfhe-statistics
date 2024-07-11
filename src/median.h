#pragma once

#include "openfhe.h"


using namespace lbcrypto;


/**
 * @brief Computes the median of a vector.
 * 
 * This function computes the median of the elements in the input vector `vec`.
 * 
 * @param vec The input vector of doubles.
 * @return double The median of the input vector.
 */
double median(
    std::vector<double> vec
);


/**
 * @brief Computes the median of a ciphertext vector.
 * 
 * This function computes the median of a ciphertext vector `c`.
 * 
 * @param c The ciphertext vector for which to compute the median.
 * @param vectorLength The length of the vector.
 * @param leftBoundC The left bound for comparison's approximation.
 * @param rightBoundC The right bound for comparison's approximation.
 * @param degreeC The degree of the comparison's approximation.
 * @param degreeI The degree of the indicator function's approximation.
 * @return Ciphertext<DCRTPoly> A ciphertext bitmask representing the position
 * of the median in the input vector.
 */
Ciphertext<DCRTPoly> median(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    double leftBoundC,
    double rightBoundC,
    uint32_t degreeC,
    uint32_t degreeI
);


/**
 * @brief Computes the median of a vector stored in multiple
 * ciphertexts.
 * 
 * This function computes the median of a ciphertext vector `c`.
 * 
 * @param c A vector containing ciphertexts, each representing a portion of the
 * input vector.
 * @param subVectorLength The length of each sub-vector stored across the
 * ciphertexts.
 * @param leftBoundC The left bound for comparison's approximation.
 * @param rightBoundC The right bound for comparison's approximation.
 * @param degreeC The degree of the comparison's approximation.
 * @param degreeI The degree of the indicator function's approximation.
 * @return Ciphertext<DCRTPoly> A ciphertext bitmask representing the position
 * of the median in the input vector.
 */
Ciphertext<DCRTPoly> median(
    const std::vector<Ciphertext<DCRTPoly>> &c,
    const size_t subVectorLength,
    double leftBoundC,
    double rightBoundC,
    uint32_t degreeC,
    uint32_t degreeI
);
