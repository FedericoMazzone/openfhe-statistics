#pragma once

#include "openfhe.h"


using namespace lbcrypto;


template <typename Element>
Ciphertext<Element> operator<<(const Ciphertext<Element> &a, int32_t index) {
    return a->GetCryptoContext()->EvalRotate(a, index);
}

template <typename Element>
Ciphertext<Element> operator>>(const Ciphertext<Element> &a, int32_t index) {
    return a->GetCryptoContext()->EvalRotate(a, -index);
}

template <typename Element>
Ciphertext<Element> operator+(const Ciphertext<Element> &a, const double &b) {
    return a->GetCryptoContext()->EvalAdd(a, b);
}

template <typename Element>
Ciphertext<Element> operator+(const Ciphertext<Element> &a, const std::vector<double> &b) {
    CryptoContext<Element> cc = a->GetCryptoContext();
    auto plaintext = cc->MakeCKKSPackedPlaintext(b);
    return cc->EvalAdd(a, plaintext);
}

template <typename Element>
Ciphertext<Element> operator-(const double &a, const Ciphertext<Element> &b) {
    return b->GetCryptoContext()->EvalSub(a, b);
}

template <typename Element>
Ciphertext<Element> operator*(const Ciphertext<Element> &a, const std::vector<double> &b) {
    CryptoContext<Element> cc = a->GetCryptoContext();
    auto plaintext = cc->MakeCKKSPackedPlaintext(b);
    return cc->EvalMult(a, plaintext);
}

template <typename Element>
Ciphertext<Element> operator*(const std::vector<double> &a, const Ciphertext<Element> &b) {
    return b * a;
}


/**
 * Generate a cryptocontext with given parameters.
 * @param integralPrecision
 * @param decimalPrecision
 * @param multiplicativeDepth
 * @param numSlots
 * @param enableBootstrap
 * @param ringDim
 * @param verbose
 * @return cryptocontext
 */
CryptoContext<DCRTPoly> generateCryptoContext(
    const usint integralPrecision,
    const usint decimalPrecision,
    const usint multiplicativeDepth,
    const usint numSlots,
    const bool enableBootstrap=false,
    const usint ringDim=0,
    const bool verbose=true
);


/**
 * Generate a public & private keys, relinearization key, rotation keys, and 
 * bootstrapping key.
 * @param cryptoContext
 * @param indices
 * @param numSlots
 * @param enableBootstrap
 * @param verbose
 * @return public & private keys
 */
KeyPair<DCRTPoly> keyGeneration(
    const CryptoContext<DCRTPoly> &cryptoContext,
    const std::vector<int32_t> &indices,
    const usint numSlots,
    const bool enableBootstrap=false,
    const bool verbose=true
);
