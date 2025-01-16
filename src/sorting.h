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
 * @return Ciphertext<DCRTPoly> A ciphertext representing the sorted vector.
 */
std::vector<Ciphertext<DCRTPoly>> sortWithCorrection(
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
 * This function sorts a ciphertext vector `c`. It uses the fg approximation of
 * the sign function.
 * 
 * @param c The ciphertext vector to sort.
 * @param vectorLength The length of the vector.
 * @param dg_c The composition degree of g for cmp (reduce the input gap).
 * @param df_c The composition degree of f for cmp (reduce the output error).
 * @param dg_i The composition degree of g for ind (reduce the input gap).
 * @param df_i The composition degree of f for ind (reduce the output error).
 * @return Ciphertext<DCRTPoly> A ciphertext representing the sorted vector.
 */
Ciphertext<DCRTPoly> sortFG(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    uint32_t dg_c,
    uint32_t df_c,
    uint32_t dg_i,
    uint32_t df_i
);


/**
 * @brief Sorts a ciphertext vector.
 * 
 * This function sorts a ciphertext vector `c`, while handling ties at the same
 * time. It uses the fg approximation of the sign function.
 * 
 * @param c The ciphertext vector to sort.
 * @param vectorLength The length of the vector.
 * @param dg_c The composition degree of g for cmp (reduce the input gap).
 * @param df_c The composition degree of f for cmp (reduce the output error).
 * @param dg_i The composition degree of g for ind (reduce the input gap).
 * @param df_i The composition degree of f for ind (reduce the output error).
 * @return Ciphertext<DCRTPoly> A ciphertext representing the sorted vector.
 */
Ciphertext<DCRTPoly> sortWithCorrectionFG(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    uint32_t dg_c,
    uint32_t df_c,
    uint32_t dg_i,
    uint32_t df_i
);


/**
 * @brief Sorts a vector stored in multiple ciphertexts.
 * 
 * This function sorts a ciphertext vector `c`, while handling ties at the same
 * time. It uses the fg approximation of the sign function.
 * 
 * @param c A vector containing ciphertexts, each representing a portion of the
 * input vector.
 * @param subVectorLength The length of each sub-vector stored across the
 * ciphertexts.
 * @param dg_c The composition degree of g for cmp (reduce the input gap).
 * @param df_c The composition degree of f for cmp (reduce the output error).
 * @param dg_i The composition degree of g for ind (reduce the input gap).
 * @param df_i The composition degree of f for ind (reduce the output error).
 * @return Ciphertext<DCRTPoly> A ciphertext representing the sorted vector.
 */
std::vector<Ciphertext<DCRTPoly>> sortWithCorrectionFG(
    const std::vector<Ciphertext<DCRTPoly>> &c,
    const size_t subVectorLength,
    uint32_t dg_c,
    uint32_t df_c,
    uint32_t dg_i,
    uint32_t df_i
);