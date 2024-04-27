#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-matrices.h"
#include "ranking.h"
#include "minimum.h"
#include "sorting.h"

// Demo: computing ranking, minimum, and sorting of a given vector.

int main()
{

    // Defining input vector (normalized in [0,1])
    std::vector<double> v = {0.83, 0.26, 0.49, 0.97, 0.12, 0.57, 0.38, 0.74};

    // Defining approximation degree of comparison and indicator function.
    const usint compareDepth = 10;
    const usint indicatorDepth = 10;



    std::cout << std::endl <<
        "Demo: computing ranking, minimum, and sorting of a given vector."
        << std::endl << std::endl;

    std::cout << "Vector: " << v << std::endl << std::endl;
    std::cout << "Compare depth: " << compareDepth << std::endl;
    std::cout << "Indicator depth: " << indicatorDepth << std::endl
              << std::endl;


    ////////////////////////////////////////////////////////////////////////
    //                       Setting up CKKS scheme                       //
    ////////////////////////////////////////////////////////////////////////

    const size_t vectorLength = v.size();

    const usint integralPrecision       = 12;
    const usint decimalPrecision        = 48;
    const usint multiplicativeDepth     = compareDepth + indicatorDepth + 1;
    const usint numSlots                = vectorLength * vectorLength;
    const bool enableBootstrap          = false;
    const usint ringDim                 = 0;
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
    std::vector<int32_t> indices = getRotationIndices(vectorLength);
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
    //                         Computing RANKING                          //
    ////////////////////////////////////////////////////////////////////////

    std::cout << "Computing RANKING..." << std::endl;
    Ciphertext<DCRTPoly> resultC = rank(
        vC,
        vectorLength,
        -1.0, 1.0,
        depth2degree(compareDepth)
    );
    Plaintext resultP;
    cryptoContext->Decrypt(keyPair.secretKey, resultC, &resultP);
    resultP->SetLength(vectorLength);
    std::vector<double> result = resultP->GetRealPackedValue();
    std::cout << "Ranking           : " << result << std::endl;
    std::cout << "Expected ranking  : " << rank(v) << std::endl << std::endl;


    ////////////////////////////////////////////////////////////////////////
    //                         Computing MINIMUM                          //
    ////////////////////////////////////////////////////////////////////////

    std::cout << "Computing MINIMUM..." << std::endl;
    resultC = min(
        vC,
        vectorLength,
        -1.0, 1.0,
        depth2degree(compareDepth),
        depth2degree(indicatorDepth)
    );
    cryptoContext->Decrypt(keyPair.secretKey, resultC, &resultP);
    resultP->SetLength(vectorLength);
    result = resultP->GetRealPackedValue();
    double computedMin = 0.0;
    double norm = 0.0;
    for (size_t i = 0; i < v.size(); i++)
    {
        computedMin += v[i] * result[i];
        norm += result[i];
    }
    computedMin /= norm;
    std::cout << "Minimum           : " << computedMin << std::endl;
    std::cout << "Expected minimum  : " << min(v) << std::endl << std::endl;


    ////////////////////////////////////////////////////////////////////////
    //                         Computing SORTING                          //
    ////////////////////////////////////////////////////////////////////////

    std::cout << "Computing SORTING..." << std::endl;
    resultC = sort(
        vC,
        vectorLength,
        -1.0, 1.0,
        depth2degree(compareDepth),
        depth2degree(indicatorDepth)
    );
    cryptoContext->Decrypt(keyPair.secretKey, resultC, &resultP);
    resultP->SetLength(vectorLength * vectorLength);
    std::vector<double> resultMatrix = resultP->GetRealPackedValue();
    for (size_t i = 0; i < vectorLength; i++)
        result[i] = resultMatrix[i * vectorLength];
    std::cout << "Sorting           : " << result << std::endl;
    std::cout << "Expected sorting  : " << sort(v) << std::endl << std::endl;


    return 0;

}
