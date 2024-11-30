#pragma once

#include "openfhe.h"


using namespace lbcrypto;


/**
 * @brief Sorts the elements of a vector in ascending order.
 * 
 * This function sorts the elements of the input vector `vec` in ascending
 * order.
 * 
 * @param vec The input vector of doubles to be sorted.
 * @return std::vector<double> The sorted form of `vec`.
 */
std::vector<double> sort(
    std::vector<double> vec
);


/**
 * @brief Sorts a ciphertext vector.
 * 
 * This function sorts a ciphertext vector `c`.
 * 
 * @param c The ciphertext vector to sort.
 * @param vectorLength The length of the vector.
 * @param leftBoundC The left bound for comparison's approximation.
 * @param rightBoundC The right bound for comparison's approximation.
 * @param degreeC The degree of the comparison's approximation.
 * @param degreeI The degree of the indicator function's approximation.
 * @return Ciphertext<DCRTPoly> A ciphertext representing the sorted vector.
 */
Ciphertext<DCRTPoly> sort(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    double leftBoundC,
    double rightBoundC,
    uint32_t degreeC,
    uint32_t degreeI
);


/**
 * @brief Sorts a vector stored in multiple ciphertexts.
 * 
 * This function sorts a ciphertext vector `c`.
 * 
 * @param c A vector containing ciphertexts, each representing a portion of the
 * input vector.
 * @param subVectorLength The length of each sub-vector stored across the
 * ciphertexts.
 * @param leftBoundC The left bound for comparison's approximation.
 * @param rightBoundC The right bound for comparison's approximation.
 * @param degreeC The degree of the comparison's approximation.
 * @param degreeI The degree of the indicator function's approximation.
 * @return Ciphertext<DCRTPoly> A ciphertext representing the sorted vector.
 */
std::vector<Ciphertext<DCRTPoly>> sort(
    const std::vector<Ciphertext<DCRTPoly>> &c,
    const size_t subVectorLength,
    double leftBoundC,
    double rightBoundC,
    uint32_t degreeC,
    uint32_t degreeI
);


/**
 * @brief Sorts a ciphertext vector.
 * 
 * This function sorts a ciphertext vector `c`, while handling ties at the same
 * time.
 * 
 * @param c The ciphertext vector to sort.
 * @param vectorLength The length of the vector.
 * @param leftBoundC The left bound for comparison's approximation.
 * @param rightBoundC The right bound for comparison's approximation.
 * @param degreeC The degree of the comparison's approximation.
 * @param degreeI The degree of the indicator function's approximation.
 * @return Ciphertext<DCRTPoly> A ciphertext representing the sorted vector.
 */
Ciphertext<DCRTPoly> sortWithCorrection(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    double leftBoundC,
    double rightBoundC,
    uint32_t degreeC,
    uint32_t degreeI
);


/**
 * @brief Sorts a vector stored in multiple ciphertexts.
 * 
 * This function sorts a ciphertext vector `c`, while handling ties at the same
 * time.
 * 
 * @param c A vector containing ciphertexts, each representing a portion of the
 * input vector.
 * @param subVectorLength The length of each sub-vector stored across the
 * ciphertexts.
 * @param leftBoundC The left bound for comparison's approximation.
 * @param rightBoundC The right bound for comparison's approximation.
 * @param degreeC The degree of the comparison's approximation.
 * @param degreeI The degree of the indicator function's approximation.
 * @param advanced Whether to use the higher-precision comparison function.
 * @return Ciphertext<DCRTPoly> A ciphertext representing the sorted vector.
 */
std::vector<Ciphertext<DCRTPoly>> sortWithCorrection(
    const std::vector<Ciphertext<DCRTPoly>> &c,
    const size_t subVectorLength,
    double leftBoundC,
    double rightBoundC,
    uint32_t degreeC,
    uint32_t degreeI,
    const bool advanced=false
);
