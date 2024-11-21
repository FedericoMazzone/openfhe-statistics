#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-ptxt.h"

// Implementation of the argmax functionality from NEXUS paper.
// https://github.com/zju-abclab/NEXUS/blob/main/src/argmax.cpp


Ciphertext<DCRTPoly> quickMax
(
    Ciphertext<DCRTPoly> a,
    const size_t vectorLength,
    uint32_t degreeC,
    KeyPair<DCRTPoly> keyPair
)
{
    // QuickMax
    a = a + (a >> vectorLength);
    for (size_t i = 0; i < LOG2(vectorLength); i++)
    {
        std::cout << "Step " << i + 1 << "/" << LOG2(vectorLength) << std::endl;
        auto t = a << (1 << i);
        auto c = compare(a, t, -2.0, 2.0, degreeC);
        a = a * c + t * (1 - c);  // a = max(a, t)

        // Plaintext resultP;
        // a->GetCryptoContext()->Decrypt(keyPair.secretKey, a, &resultP);
        // resultP->SetLength(vectorLength);
        // std::cout << "Result: " << resultP;

        std::cout << "Bootstraping..." << std::endl;
        a = a->GetCryptoContext()->EvalBootstrap(a);
    }
    return a;
}


Ciphertext<DCRTPoly> argmaxNEXUS
(
    Ciphertext<DCRTPoly> a,
    const size_t vectorLength,
    uint32_t degreeC,
    KeyPair<DCRTPoly> keyPair
)
{
    auto aMax = quickMax(a, vectorLength, degreeC, keyPair);
    auto b = compare(a, aMax, -2.0, 2.0, degreeC);
    b = b * 2;

    return b;
}


int main()
{

    const size_t vectorLength = 16384; // 128;

    // Generating input vector (in [0,1])
    std::vector<double> v(vectorLength);
    for (size_t i = 0; i < vectorLength - 1; i++)
        v[i] = 0.2;
    v[vectorLength - 1] = 0.7;

    // Defining approximation degree of comparison function.
    const usint compareDepth = 10;

    std::cout << std::endl << "Demo: computing argmax with NEXUS' approach."
              << std::endl << std::endl;

    std::cout << "Vector: " << std::fixed << std::setprecision(3) << v << std::endl << std::endl;
    std::cout << "Compare depth: " << compareDepth << std::endl << std::endl;


    ////////////////////////////////////////////////////////////////////////
    //                       Setting up CKKS scheme                       //
    ////////////////////////////////////////////////////////////////////////

    const usint integralPrecision       = 1;
    const usint decimalPrecision        = 59;
    const usint multiplicativeDepth     = compareDepth + 5;
    const usint numSlots                = vectorLength * 2;
    const bool enableBootstrap          = true;
    const usint ringDim                 = 1 << 16; // to match up NEXUS' parameters
    const bool verbose                  = true;

    CryptoContext<DCRTPoly> cryptoContext = generateCryptoContext(
        integralPrecision,
        decimalPrecision,
        multiplicativeDepth,
        numSlots,
        enableBootstrap,
        ringDim,
        verbose
    );

    // Generating public/private key pair, relinearization, and rotation keys
    std::vector<int32_t> indices = {};
    for (size_t i = 0; i < LOG2(vectorLength); i++)
        indices.push_back(1 << i);
    KeyPair<DCRTPoly> keyPair = keyGeneration(
        cryptoContext,
        indices,
        numSlots,
        enableBootstrap,
        verbose
    );


    ////////////////////////////////////////////////////////////////////////
    //                      Encrypting input vector                       //
    ////////////////////////////////////////////////////////////////////////

    std::cout << "Encrypting input vector...          " << std::flush;
    Ciphertext<DCRTPoly> vC = cryptoContext->Encrypt(
        keyPair.publicKey,
        cryptoContext->MakeCKKSPackedPlaintext(v)
    );
    std::cout << "COMPLETED" << std::endl << std::endl;


    ////////////////////////////////////////////////////////////////////////
    //                          Computing ArgMax                          //
    ////////////////////////////////////////////////////////////////////////

    std::cout << "Computing ArgMax..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    Ciphertext<DCRTPoly> resultC = argmaxNEXUS(
        vC,
        vectorLength,
        depth2degree(compareDepth),
        keyPair
    );
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "runtime: " << elapsed_seconds.count() << "s" << std::endl;
    Plaintext resultP;
    cryptoContext->Decrypt(keyPair.secretKey, resultC, &resultP);
    resultP->SetLength(vectorLength);
    std::vector<double> result = resultP->GetRealPackedValue();
    std::cout << "Result: " << result << std::endl << std::endl;
    size_t computedArgmax = std::distance(result.begin(), std::max_element(result.begin(), result.end()));
    size_t expectedArgmax = std::distance(v.begin(), std::max_element(v.begin(), v.end()));
    std::cout << "Computed ArgMax  : " << computedArgmax << std::endl;
    std::cout << "Expected ArgMax  : " << expectedArgmax << std::endl << std::endl;


    return 0;

}
