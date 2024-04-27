#pragma once

#include "utils-basics.h"


using namespace lbcrypto;


/**
 * Compute the highest polynomial degree that can be evaluated in a circuit of
 * given depth. See https://github.com/openfheorg/openfhe-development/blob/main/src/pke/examples/FUNCTION_EVALUATION.md
 * @param depth
 * @return degree
 */
usint depth2degree(
    const usint depth
);


/**
 * @brief Compares two ciphertexts c1 and c2 as (c1 > c2).
 * 
 * This function compares two ciphertexts using the Chebyshev approximation of
 * the sign function. The result is a ciphertext, where a value of 1 indicates
 * that the corresponding components of c1 are greater than those of c2, a
 * value of 0 indicates that the corresponding components of c1 are less than
 * those of c2, and a value of 0.5 indicates that the corresponding components
 * of c1 are approximately equal to those of c2.
 * 
 * @param c1 The first ciphertext to be compared.
 * @param c2 The second ciphertext to be compared.
 * @param a The lower bound of the approximation interval.
 * @param b The upper bound of the approximation interval.
 * @param degree The degree of the Chebyshev polynomial approximation.
 * @param error The threshold for considering values close to zero (default is
 * 0.00001).
 * @return Ciphertext<DCRTPoly> A ciphertext representing the output of the
 * comparison (c1 > c2).
 */
Ciphertext<DCRTPoly> compare(
    const Ciphertext<DCRTPoly> &c1,
    const Ciphertext<DCRTPoly> &c2,
    double a,
    double b,
    uint32_t degree,
    double error = 0.00001
);



/**
 * @brief Compares two ciphertexts c1 and c2 as (c1 > c2).
 * 
 * This function compares two ciphertexts using the Chebyshev approximation of
 * the sign function. The result is a ciphertext, where a value of 1 indicates
 * that the corresponding components of c1 are greater than those of c2, and a
 * value of 0 indicates that the corresponding components of c1 are less or
 * equal than those of c2.
 * 
 * @param c1 The first ciphertext to be compared.
 * @param c2 The second ciphertext to be compared.
 * @param a The lower bound of the approximation interval.
 * @param b The upper bound of the approximation interval.
 * @param degree The degree of the Chebyshev polynomial approximation.
 * @param error The threshold for considering values close to zero (default is
 * 0.00001).
 * @return Ciphertext<DCRTPoly> A ciphertext representing the output of the
 * comparison (c1 > c2).
 */
Ciphertext<DCRTPoly> compareGt(
    const Ciphertext<DCRTPoly> &c1,
    const Ciphertext<DCRTPoly> &c2,
    double a,
    double b,
    uint32_t degree,
    double error = 0.00001
);


/**
 * @brief Evaluate the indicator function for an interval [a1, b1].
 * 
 * This function evaluate an indicator function using Chebyshev approximation.
 * If the value falls within the interval [a1, b1], the indicator function
 * returns 1; otherwise, it returns 0.
 * 
 * @param c The input ciphertext.
 * @param a1 The lower bound of the indicator interval.
 * @param b1 The upper bound of the indicator interval.
 * @param a The lower bound of the approximation interval.
 * @param b The upper bound of the approximation interval.
 * @param degree The degree of the Chebyshev polynomial approximation.
 * @return Ciphertext<DCRTPoly> A ciphertext representing the output of the
 * indicator function.
 */
Ciphertext<DCRTPoly> indicator(
    const Ciphertext<DCRTPoly> &c,
    double a1,
    double b1,
    double a,
    double b,
    uint32_t degree
);
