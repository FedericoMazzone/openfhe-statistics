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
    const double leftBoundC,
    const double rightBoundC,
    const uint32_t degreeC,
    const uint32_t degreeI
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
    const double leftBoundC,
    const double rightBoundC,
    const uint32_t degreeC,
    const uint32_t degreeI
);


/**
 * @brief Computes the median of a ciphertext vector.
 * 
 * This function computes the median of a ciphertext vector `c`. It uses the fg
 * approximation of the sign function.
 * 
 * @param c The ciphertext vector for which to compute the median.
 * @param vectorLength The length of the vector.
 * @param dg_c The composition degree of g for cmp (reduce the input gap).
 * @param df_c The composition degree of f for cmp (reduce the output error).
 * @param dg_i The composition degree of g for ind (reduce the input gap).
 * @param df_i The composition degree of f for ind (reduce the output error).
 * @return Ciphertext<DCRTPoly> A ciphertext bitmask representing the position
 * of the median in the input vector.
 */
Ciphertext<DCRTPoly> medianFG(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    const uint32_t dg_c,
    const uint32_t df_c,
    const uint32_t dg_i,
    const uint32_t df_i
);


/**
 * @brief Computes the median of a vector stored in multiple
 * ciphertexts.
 * 
 * This function computes the median of a ciphertext vector `c`. It uses the fg
 * approximation of the sign function.
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
 * of the median in the input vector.
 */
Ciphertext<DCRTPoly> medianFG(
    const std::vector<Ciphertext<DCRTPoly>> &c,
    const size_t subVectorLength,
    const uint32_t dg_c,
    const uint32_t df_c,
    const uint32_t dg_i,
    const uint32_t df_i
);
