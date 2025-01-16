#pragma once

#include "openfhe.h"


using namespace lbcrypto;


/**
 * @brief Finds the minimum element in a vector of doubles.
 * 
 * This function finds the minimum element in the input vector `vec`.
 * 
 * @param vec The input vector of doubles.
 * @return double The minimum element found in the input vector.
 */
double min(
    const std::vector<double> &vec
);


/**
 * @brief Computes the minimum value in a ciphertext vector.
 * 
 * This function computes the minimum value in a ciphertext vector `c`.
 * 
 * @param c The ciphertext vector for which to compute the minimum value.
 * @param vectorLength The length of the vector.
 * @param leftBoundC The left bound for comparison's approximation.
 * @param rightBoundC The right bound for comparison's approximation.
 * @param degreeC The degree of the comparison's approximation.
 * @param degreeI The degree of the indicator function's approximation.
 * @return Ciphertext<DCRTPoly> A ciphertext bitmask representing the position
 * of the minimum value in the input vector.
 */
Ciphertext<DCRTPoly> min(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    const double leftBoundC,
    const double rightBoundC,
    const uint32_t degreeC,
    const uint32_t degreeI
);


/**
 * @brief Computes the minimum value in a vector stored in multiple
 * ciphertexts.
 * 
 * This function computes the minimum value in a ciphertext vector `c`.
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
 * of the minimum value in the input vector.
 */
Ciphertext<DCRTPoly> min(
    const std::vector<Ciphertext<DCRTPoly>> &c,
    const size_t subVectorLength,
    const double leftBoundC,
    const double rightBoundC,
    const uint32_t degreeC,
    const uint32_t degreeI
);


/**
 * @brief Computes the minimum value in a vector stored in multiple
 * ciphertexts.
 * 
 * This function computes the minimum value in a ciphertext vector `c`. It uses
 * the fg approximation of the sign function.
 * 
 * @param c A vector containing ciphertexts, each representing a portion of the
 * input vector.
 * @param subVectorLength The length of each sub-vector stored across the
 * ciphertexts.
 * @param dg_c The composition degree of g for cmp (reduce the input gap).
 * @param df_c The composition degree of f for cmp (reduce the output error).
 * @param dg_i The composition degree of g for ind (reduce the input gap).
 * @param df_i The composition degree of f for ind (reduce the output error).
 * @return Ciphertext<DCRTPoly> A ciphertext bitmask representing the position
 * of the minimum value in the input vector.
 */
Ciphertext<DCRTPoly> minFG(
    const std::vector<Ciphertext<DCRTPoly>> &c,
    const size_t subVectorLength,
    const uint32_t dg_c,
    const uint32_t df_c,
    const uint32_t dg_i,
    const uint32_t df_i
);
