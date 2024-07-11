#pragma once

#include "openfhe.h"


using namespace lbcrypto;


/**
 * @brief Rank the elements in a vector.
 * 
 * This function computes the ranks of elements in the input vector `vec`. In
 * the fractional mode, elements in a tie are assigned the average of the ranks
 * they span.
 * 
 * @param vec The input vector of elements for which to compute fractional
 * ranks.
 * @param fractional Whether to compute the fractioncal rank or not (default is
 * true).
 * @param epsilon The tolerance level for considering two elements as equal.
 * @return std::vector<double> A vector containing the fractional ranks of
 * elements in the input vector.
 */
std::vector<double> rank(
    const std::vector<double> &vec,
    const bool fractional = true,
    const double epsilon = 0.0
);


/**
 * @brief Computes the rank of elements in a ciphertext vector.
 * 
 * This function computes the rank of elements in a ciphertext vector `c`.
 * 
 * @param c The ciphertext vector for which to compute the rank.
 * @param vectorLength The length of the vector.
 * @param leftBoundC The left bound for comparison's approximation.
 * @param rightBoundC The right bound for comparison's approximation.
 * @param degreeC The degree of the comparison's approximation.
 * @param cmpGt Flag indicating whether to compute standard (true) or
 * fractional rank (false), default is false.
 * @return Ciphertext<DCRTPoly> The ciphertext vector containing the computed
 * ranks.
 */
Ciphertext<DCRTPoly> rank(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    double leftBoundC,
    double rightBoundC,
    uint32_t degreeC,
    bool cmpGt = false
);


/**
 * @brief Computes the rank of elements in a vector stored in multiple
 * ciphertexts.
 * 
 * This function computes the rank of elements in a ciphertext vector `c`.
 * 
 * @param c A vector containing ciphertexts, each representing a portion of the
 * input vector.
 * @param subVectorLength The length of each sub-vector stored across the
 * ciphertexts.
 * @param leftBoundC The left bound for comparison's approximation.
 * @param rightBoundC The right bound for comparison's approximation.
 * @param degreeC The degree of the comparison's approximation.
 * @param cmpGt Flag indicating whether to compute standard (true) or
 * fractional rank (false), default is false.
 * @param complOpt Flag indicating whether to use the complementary comparison
 * optimization.
 * @return std::vector<Ciphertext<DCRTPoly>> A vector containing ciphertexts
 * representing the computed ranks.
 */
std::vector<Ciphertext<DCRTPoly>> rank(
    const std::vector<Ciphertext<DCRTPoly>> &c,
    const size_t subVectorLength,
    double leftBoundC,
    double rightBoundC,
    uint32_t degreeC,
    bool cmpGt = false,
    bool complOpt = true
);


/**
 * @brief Computes the rank of elements in a ciphertext vector.
 * 
 * This function computes the rank of elements in a ciphertext vector `c`,
 * while handling ties at the same time.
 * 
 * @param c The ciphertext vector for which to compute the rank.
 * @param vectorLength The length of the vector.
 * @param leftBoundC The left bound for comparison's approximation.
 * @param rightBoundC The right bound for comparison's approximation.
 * @param degreeC The degree of the comparison's approximation.
 * @param parallel Flag indicating whether to compute the correction offset in
 * parallel or not, default is false.
 * @return Ciphertext<DCRTPoly> The ciphertext vector containing the computed
 * ranks.
 */
Ciphertext<DCRTPoly> rankWithCorrection(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    double leftBoundC,
    double rightBoundC,
    uint32_t degreeC,
    bool parallel = false
);


/**
 * @brief Computes the rank of elements in a vector stored in multiple
 * ciphertexts.
 * 
 * This function computes the rank of elements in a ciphertext vector `c`,
 * while handling ties at the same time.
 * 
 * @param c A vector containing ciphertexts, each representing a portion of the
 * input vector.
 * @param subVectorLength The length of each sub-vector stored across the
 * ciphertexts.
 * @param leftBoundC The left bound for comparison's approximation.
 * @param rightBoundC The right bound for comparison's approximation.
 * @param degreeC The degree of the comparison's approximation.
 * @return std::vector<Ciphertext<DCRTPoly>> A vector containing ciphertexts
 * representing the computed ranks.
 */
std::vector<Ciphertext<DCRTPoly>> rankWithCorrection(
    const std::vector<Ciphertext<DCRTPoly>> &c,
    const size_t subVectorLength,
    double leftBoundC,
    double rightBoundC,
    uint32_t degreeC
);
