#include "sorting.h"
#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-matrices.h"
#include "utils-ptxt.h"

#include <omp.h>



void testSorting(
    const size_t vectorLength,
    const usint compareDepth,
    const usint indicatorDepth,
    const bool tieCorrection = false
)
{

    const usint integralPrecision       = 1;
    const usint decimalPrecision        = 59;
    const usint multiplicativeDepth     = compareDepth + indicatorDepth + 1 + (tieCorrection ? 3 : 0);
    const usint numSlots                = vectorLength * vectorLength;
    const bool enableBootstrap          = false;
    const usint ringDim                 = 0;
    const bool verbose                  = true;

    std::vector<int32_t> indices = getRotationIndices(vectorLength);

    CryptoContext<DCRTPoly> cryptoContext = generateCryptoContext(
        integralPrecision,
        decimalPrecision,
        multiplicativeDepth,
        numSlots,
        enableBootstrap,
        ringDim,
        verbose
    );

    KeyPair<DCRTPoly> keyPair = keyGeneration(
        cryptoContext,
        indices,
        numSlots,
        enableBootstrap,
        verbose
    );

    std::vector<double> v = loadPoints1D(vectorLength);

    std::cout << "Vector:           " << v << std::endl;

    std::cout << "Expected sorting: " << sort(v) << std::endl;

    Ciphertext<DCRTPoly> vC = cryptoContext->Encrypt(
        keyPair.publicKey,
        cryptoContext->MakeCKKSPackedPlaintext(v)
    );

    auto start = std::chrono::high_resolution_clock::now();

    Ciphertext<DCRTPoly> resultC;
    if (tieCorrection)
        resultC = sortWithCorrection(
            vC,
            vectorLength,
            -1.0, 1.0,
            depth2degree(compareDepth),
            depth2degree(indicatorDepth)
        );
    else
        resultC = sort(
            vC,
            vectorLength,
            -1.0, 1.0,
            depth2degree(compareDepth),
            depth2degree(indicatorDepth)
        );

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Runtime: " << elapsed_seconds.count() << "s" << std::endl;

    Plaintext resultP;
    cryptoContext->Decrypt(keyPair.secretKey, resultC, &resultP);
    resultP->SetLength(vectorLength * vectorLength);
    std::vector<double> resultMatrix = resultP->GetRealPackedValue();
    std::vector<double> result(vectorLength);
    for (size_t i = 0; i < vectorLength; i++)
        result[i] = resultMatrix[i * vectorLength];
    std::cout << "Sorting:          " << result << std::endl;

}


void testSortingMultiCtxt(
    const size_t subVectorLength,
    const size_t numCiphertext,
    const usint compareDepth,
    const usint indicatorDepth,
    const bool tieCorrection = false
)
{

    const size_t vectorLength           = subVectorLength * numCiphertext;
    const usint integralPrecision       = 1;
    const usint decimalPrecision        = 59;
    const usint multiplicativeDepth     = compareDepth + indicatorDepth + 3 + (tieCorrection ? 3 : 0);
    const usint numSlots                = subVectorLength * subVectorLength;
    const bool enableBootstrap          = false;
    const usint ringDim                 = 0;
    const bool verbose                  = true;

    std::vector<int32_t> indices = getRotationIndices(subVectorLength);

    CryptoContext<DCRTPoly> cryptoContext = generateCryptoContext(
        integralPrecision,
        decimalPrecision,
        multiplicativeDepth,
        numSlots,
        enableBootstrap,
        ringDim,
        verbose
    );

    KeyPair<DCRTPoly> keyPair = keyGeneration(
        cryptoContext,
        indices,
        numSlots,
        enableBootstrap,
        verbose
    );

    std::vector<double> v = loadPoints1D(vectorLength);
    std::vector<std::vector<double>> vTokens = splitVector(v, numCiphertext);

    std::cout << "Vector: " << vTokens << std::endl;

    std::cout << "Expected sorting: " << sort(v) << std::endl;

    std::vector<Ciphertext<DCRTPoly>> vC(numCiphertext);
    for (size_t j = 0; j < numCiphertext; j++)
        vC[j] = cryptoContext->Encrypt(
            keyPair.publicKey,
            cryptoContext->MakeCKKSPackedPlaintext(vTokens[j])
        );

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<Ciphertext<DCRTPoly>> resultC;
    if (tieCorrection)
        resultC = sortWithCorrection(
            vC,
            subVectorLength,
            -1.0, 1.0,
            depth2degree(compareDepth),
            depth2degree(indicatorDepth)
        );
    else
        resultC = sort(
            vC,
            subVectorLength,
            -1.0, 1.0,
            depth2degree(compareDepth),
            depth2degree(indicatorDepth)
        );

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Runtime: " << elapsed_seconds.count() << "s" << std::endl;

    Plaintext resultP;
    std::vector<std::vector<double>> resultTokens(numCiphertext);
    for (size_t i = 0; i < numCiphertext; i++)
    {
        cryptoContext->Decrypt(keyPair.secretKey, resultC[i], &resultP);
        resultP->SetLength(subVectorLength * subVectorLength);
        std::vector<double> resultMatrix = resultP->GetRealPackedValue();
        std::vector<double> resultToken(subVectorLength);
        for (size_t j = 0; j < subVectorLength; j++)
            resultToken[j] = resultMatrix[j * subVectorLength];
        resultTokens[i] = resultToken;
    }
    std::vector<double> result = concatVectors(resultTokens);
    std::cout << "Sorting: " << result << std::endl;

}


void testSortingAdv(
    const size_t vectorLength,
    const usint dg_c,
    const usint df_c,
    const usint dg_i,
    const usint df_i,
    const bool tieCorrection = false
)
{

    const usint integralPrecision       = 1;
    const usint decimalPrecision        = 59;
    const usint multiplicativeDepth     = 4 * (dg_c + df_c + dg_i + df_i) + 4 + (tieCorrection ? 3 : 0);
    const usint numSlots                = vectorLength * vectorLength;
    const bool enableBootstrap          = false;
    const usint ringDim                 = 0;
    const bool verbose                  = true;

    std::vector<int32_t> indices = getRotationIndices(vectorLength);

    CryptoContext<DCRTPoly> cryptoContext = generateCryptoContext(
        integralPrecision,
        decimalPrecision,
        multiplicativeDepth,
        numSlots,
        enableBootstrap,
        ringDim,
        verbose
    );

    KeyPair<DCRTPoly> keyPair = keyGeneration(
        cryptoContext,
        indices,
        numSlots,
        enableBootstrap,
        verbose
    );

    std::vector<double> v = loadPoints1D(vectorLength);

    std::cout << "Vector:           " << v << std::endl;

    std::cout << "Expected sorting: " << sort(v) << std::endl;

    Ciphertext<DCRTPoly> vC = cryptoContext->Encrypt(
        keyPair.publicKey,
        cryptoContext->MakeCKKSPackedPlaintext(v)
    );

    auto start = std::chrono::high_resolution_clock::now();

    Ciphertext<DCRTPoly> resultC;
    if (tieCorrection)
        resultC = sortWithCorrectionFG(
            vC,
            vectorLength,
            dg_c, df_c,
            dg_i, df_i
        );
    else
        resultC = sortFG(
            vC,
            vectorLength,
            dg_c, df_c,
            dg_i, df_i
        );

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Runtime: " << elapsed_seconds.count() << "s" << std::endl;

    Plaintext resultP;
    cryptoContext->Decrypt(keyPair.secretKey, resultC, &resultP);
    resultP->SetLength(vectorLength * vectorLength);
    std::vector<double> resultMatrix = resultP->GetRealPackedValue();
    std::vector<double> result(vectorLength);
    for (size_t i = 0; i < vectorLength; i++)
        result[i] = resultMatrix[i * vectorLength];
    std::cout << "Sorting:          " << result << std::endl;

}


void testSortingMultiCtxtAdv(
    const size_t subVectorLength,
    const size_t numCiphertext,
    const usint dg_c,
    const usint df_c,
    const usint dg_i,
    const usint df_i,
    const bool tieCorrection = false
)
{

    const size_t vectorLength           = subVectorLength * numCiphertext;
    const usint integralPrecision       = 1;
    const usint decimalPrecision        = 59;
    const usint multiplicativeDepth     = 4 * (dg_c + df_c + dg_i + df_i) + 6 + (tieCorrection ? 3 : 0);
    const usint numSlots                = subVectorLength * subVectorLength;
    const bool enableBootstrap          = false;
    const usint ringDim                 = 0;
    const bool verbose                  = true;

    std::vector<int32_t> indices = getRotationIndices(subVectorLength);

    CryptoContext<DCRTPoly> cryptoContext = generateCryptoContext(
        integralPrecision,
        decimalPrecision,
        multiplicativeDepth,
        numSlots,
        enableBootstrap,
        ringDim,
        verbose
    );

    KeyPair<DCRTPoly> keyPair = keyGeneration(
        cryptoContext,
        indices,
        numSlots,
        enableBootstrap,
        verbose
    );

    std::vector<double> v = loadPoints1D(vectorLength);
    std::vector<std::vector<double>> vTokens = splitVector(v, numCiphertext);

    std::cout << "Vector: " << vTokens << std::endl;

    std::cout << "Expected sorting: " << sort(v) << std::endl;

    std::vector<Ciphertext<DCRTPoly>> vC(numCiphertext);
    for (size_t j = 0; j < numCiphertext; j++)
        vC[j] = cryptoContext->Encrypt(
            keyPair.publicKey,
            cryptoContext->MakeCKKSPackedPlaintext(vTokens[j])
        );

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<Ciphertext<DCRTPoly>> resultC;
    if (tieCorrection)
        resultC = sortWithCorrectionFG(
            vC,
            subVectorLength,
            dg_c, df_c,
            dg_i, df_i
        );
    else
    {}

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Runtime: " << elapsed_seconds.count() << "s" << std::endl;

    Plaintext resultP;
    std::vector<std::vector<double>> resultTokens(numCiphertext);
    for (size_t i = 0; i < numCiphertext; i++)
    {
        cryptoContext->Decrypt(keyPair.secretKey, resultC[i], &resultP);
        resultP->SetLength(subVectorLength * subVectorLength);
        std::vector<double> resultMatrix = resultP->GetRealPackedValue();
        std::vector<double> resultToken(subVectorLength);
        for (size_t j = 0; j < subVectorLength; j++)
            resultToken[j] = resultMatrix[j * subVectorLength];
        resultTokens[i] = resultToken;
    }
    std::vector<double> result = concatVectors(resultTokens);
    std::cout << "Sorting: " << result << std::endl;

}


int main(int argc, char *argv[])
{

    const size_t vectorLength = std::stoul(argv[1]);
    const bool singleThread = (argc > 2) ? (bool) std::stoi(argv[2]) : false;

    std::cout << "Vector length         : " << vectorLength << std::endl;
    std::cout << "Single thread         : " << (singleThread ? "true" : "false") << std::endl;

    std::cout << std::fixed << std::setprecision(2);
    const int numThreads = std::thread::hardware_concurrency();
    std::cout << "Number of threads     : " << numThreads << std::endl;

    const size_t subVectorLength = 256;
    const size_t numCiphertext = std::ceil((double) vectorLength / subVectorLength);

    const size_t dg_c = 3;
    const size_t df_c = 2;
    const size_t dg_i = (log2(vectorLength) + 1) / 2;
    const size_t df_i = 2;

    std::cout << "Subvector length      : " << subVectorLength << std::endl;
    std::cout << "Number of ciphertexts : " << numCiphertext << std::endl;
    std::cout << "dg_c                  : " << dg_c << std::endl;
    std::cout << "df_c                  : " << df_c << std::endl;
    std::cout << "dg_i                  : " << dg_i << std::endl;
    std::cout << "df_i                  : " << df_i << std::endl << std::endl;
    
    if (numCiphertext == 1)
    {
        if (singleThread)
        {
            #pragma omp parallel for
            for (size_t i = 0; i < 1; i++)
                testSortingAdv(vectorLength, dg_c, df_c, dg_i, df_i);
        }
        else
        {
            omp_set_num_threads(numThreads);
            if (numCiphertext <= numThreads / 16)
                omp_set_max_active_levels(10);
            testSortingAdv(vectorLength, dg_c, df_c, dg_i, df_i);
        }
    }
    else
    {
        if (singleThread)
        {
            #pragma omp parallel for
            for (size_t i = 0; i < 1; i++)
                testSortingMultiCtxtAdv(subVectorLength, numCiphertext, dg_c, df_c, dg_i, df_i);
        }
        else
        {
            omp_set_num_threads(numThreads);
            if (numCiphertext <= numThreads / 16)
                omp_set_max_active_levels(10);
            testSortingMultiCtxtAdv(subVectorLength, numCiphertext, dg_c, df_c, dg_i, df_i);
        }
    }

    return 0;

}
