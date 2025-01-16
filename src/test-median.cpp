#include "median.h"
#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-matrices.h"
#include "utils-ptxt.h"

#include <cassert>
#include <omp.h>


double evaluateMedianValue(
    const std::vector<double> &vec,
    const std::vector<double> &computedMedianMask
)
{
    assert(vec.size() == computedMedianMask.size());

    double expectedMedian = median(vec);

    double computedMedian = 0.0;
    double norm = 0.0;
    for (size_t i = 0; i < vec.size(); i++)
    {
        computedMedian += vec[i] * computedMedianMask[i];
        norm += computedMedianMask[i];
    }
    computedMedian /= norm;

    return std::abs(expectedMedian - computedMedian) / expectedMedian;
}


double evaluateMedianMask(
    const std::vector<double> &vec,
    const std::vector<double> &computedMedianMask
)
{
    assert(vec.size() == computedMedianMask.size());

    double expectedMedian = median(vec);
    size_t size = vec.size();
    size_t posMax, posMax2;
    if (computedMedianMask[0] > computedMedianMask[1]) {posMax = 0; posMax2 = 1;}
    else                                               {posMax = 1; posMax2 = 0;}
    for (size_t i = 2; i < size; i++)
        if (computedMedianMask[i] > computedMedianMask[posMax2])
        {
            if (computedMedianMask[i] > computedMedianMask[posMax])
            {
                posMax2 = posMax;
                posMax = i;
            }
            else
            {
                posMax2 = i;
            }
        }
    double computedMedian;
    if (size % 2 == 0)
        computedMedian = (vec[posMax] + vec[posMax2]) / 2.0;
    else
        computedMedian = vec[posMax];

    return std::abs(expectedMedian - computedMedian) / expectedMedian;
}


void testMedian(
    const size_t vectorLength,
    const usint compareDepth,
    const usint indicatorDepth
)
{

    const usint integralPrecision       = 10;
    const usint decimalPrecision        = 50;
    const usint multiplicativeDepth     = compareDepth + indicatorDepth + 3;
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

    std::cout << "Vector: " << v << std::endl;

    std::cout << "Expected median: " << median(v) << std::endl;

    Ciphertext<DCRTPoly> vC = cryptoContext->Encrypt(
        keyPair.publicKey,
        cryptoContext->MakeCKKSPackedPlaintext(v)
    );

    auto start = std::chrono::high_resolution_clock::now();

    Ciphertext<DCRTPoly> resultC = median(
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
    resultP->SetLength(vectorLength);
    std::vector<double> result = resultP->GetRealPackedValue();
    std::cout << "Median: " << result << std::endl;

    double errorValue = evaluateMedianValue(v, result);
    double errorMask = evaluateMedianMask(v, result);
    std::cout << "Error value: " << errorValue << std::endl;
    std::cout << "Error mask: " << errorMask << std::endl;

}


void testMedianMultiCtxt(
    const size_t subVectorLength,
    const size_t numCiphertext,
    const usint compareDepth,
    const usint indicatorDepth
)
{

    const size_t vectorLength           = subVectorLength * numCiphertext;
    const usint integralPrecision       = 12;
    const usint decimalPrecision        = 48;
    const usint multiplicativeDepth     = compareDepth + indicatorDepth + 2 + 3;
    const usint numSlots                = subVectorLength * subVectorLength;
    const bool enableBootstrap          = false;
    const usint ringDim                 = 0;
    const bool verbose                  = true;

    std::vector<int32_t> indices = getRotationIndices(subVectorLength);
    for (size_t j = 0; j < numCiphertext; j++)
        indices.push_back(-j * subVectorLength);

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

    std::cout << "Expected median: " << median(v) << std::endl;

    std::vector<Ciphertext<DCRTPoly>> vC(numCiphertext);
    for (size_t j = 0; j < numCiphertext; j++)
        vC[j] = cryptoContext->Encrypt(
            keyPair.publicKey,
            cryptoContext->MakeCKKSPackedPlaintext(vTokens[j])
        );

    auto start = std::chrono::high_resolution_clock::now();

    Ciphertext<DCRTPoly> resultC = median(
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
    cryptoContext->Decrypt(keyPair.secretKey, resultC, &resultP);
    resultP->SetLength(vectorLength);
    std::vector<double> result = resultP->GetRealPackedValue();
    std::cout << "Median: " << result << std::endl;

    double errorValue = evaluateMedianValue(v, result);
    double errorMask = evaluateMedianMask(v, result);
    std::cout << "Error value: " << errorValue << std::endl;
    std::cout << "Error mask: " << errorMask << std::endl;

}


void testMedianAdv(
    const size_t vectorLength,
    const usint dg_c,
    const usint df_c,
    const usint dg_i,
    const usint df_i
)
{

    const usint integralPrecision       = 1;
    const usint decimalPrecision        = 59;
    const usint multiplicativeDepth     = 4 * (dg_c + df_c + dg_i + df_i) + 6;
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

    std::cout << "Vector:          " << v << std::endl;

    std::cout << "Expected median: " << median(v) << std::endl;

    Ciphertext<DCRTPoly> vC = cryptoContext->Encrypt(
        keyPair.publicKey,
        cryptoContext->MakeCKKSPackedPlaintext(v)
    );

    auto start = std::chrono::high_resolution_clock::now();

    Ciphertext<DCRTPoly> resultC = medianFG(
        vC,
        vectorLength,
        dg_c, df_c,
        dg_i, df_i
    );
    std::cout << "Remaining levels: " << resultC->GetLevel() << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Runtime: " << elapsed_seconds.count() << "s" << std::endl;

    Plaintext resultP;
    cryptoContext->Decrypt(keyPair.secretKey, resultC, &resultP);
    resultP->SetLength(vectorLength);
    std::vector<double> result = resultP->GetRealPackedValue();
    std::cout << "Median: " << result << std::endl;

    double errorValue = evaluateMedianValue(v, result);
    double errorMask = evaluateMedianMask(v, result);
    std::cout << "Error value: " << errorValue << std::endl;
    std::cout << "Error mask: " << errorMask << std::endl;

}


void testMedianMultiCtxtAdv(
    const size_t subVectorLength,
    const size_t numCiphertext,
    const usint dg_c,
    const usint df_c,
    const usint dg_i,
    const usint df_i
)
{

    const size_t vectorLength           = subVectorLength * numCiphertext;
    const usint integralPrecision       = 1;
    const usint decimalPrecision        = 59;
    const usint multiplicativeDepth     = 4 * (dg_c + df_c + dg_i + df_i) + 3 + 2 + 3;
    const usint numSlots                = subVectorLength * subVectorLength;
    const bool enableBootstrap          = false;
    const usint ringDim                 = 0;
    const bool verbose                  = true;

    std::vector<int32_t> indices = getRotationIndices(subVectorLength);
    for (size_t j = 0; j < numCiphertext; j++)
        indices.push_back(-j * subVectorLength);

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

    std::cout << "Vector:          " << vTokens << std::endl;

    std::cout << "Expected median: " << median(v) << std::endl;

    std::vector<Ciphertext<DCRTPoly>> vC(numCiphertext);
    for (size_t j = 0; j < numCiphertext; j++)
        vC[j] = cryptoContext->Encrypt(
            keyPair.publicKey,
            cryptoContext->MakeCKKSPackedPlaintext(vTokens[j])
        );

    auto start = std::chrono::high_resolution_clock::now();

    Ciphertext<DCRTPoly> resultC = medianFG(
        vC,
        subVectorLength,
        dg_c, df_c,
        dg_i, df_i
    );
    std::cout << "Remaining levels: " << resultC->GetLevel() << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Runtime: " << elapsed_seconds.count() << "s" << std::endl;

    Plaintext resultP;
    cryptoContext->Decrypt(keyPair.secretKey, resultC, &resultP);
    resultP->SetLength(vectorLength);
    std::vector<double> result = resultP->GetRealPackedValue();
    std::cout << "Median: " << result << std::endl;

    double errorValue = evaluateMedianValue(v, result);
    double errorMask = evaluateMedianMask(v, result);
    std::cout << "Error value: " << errorValue << std::endl;
    std::cout << "Error mask: " << errorMask << std::endl;

}


int main(int argc, char *argv[])
{

    const size_t vectorLength = std::stoul(argv[1]);
    const bool singleThread = (argc > 2) ? (bool) std::stoi(argv[2]) : false;

    std::cout << "Vector length         : " << vectorLength << std::endl;
    std::cout << "Single thread         : " << (singleThread ? "true" : "false") << std::endl;

    std::cout << std::fixed << std::setprecision(2);
    int numThreads = std::thread::hardware_concurrency();
    std::cout << "Number of threads     : " << numThreads << std::endl;

    const size_t subVectorLength = 256;
    const size_t numCiphertext = std::ceil((double) vectorLength / subVectorLength);

    const size_t dg_c = 3;
    const size_t df_c = 2;
    const size_t dg_i = (log2(vectorLength)) / 2;
    const size_t df_i = 1;

    std::cout << "Subvector length      : " << subVectorLength << std::endl;
    std::cout << "Number of ciphertexts : " << numCiphertext << std::endl;
    std::cout << "dg_c                  : " << dg_c << std::endl;
    std::cout << "df_c                  : " << df_c << std::endl;
    std::cout << "dg_i                  : " << dg_i << std::endl;
    std::cout << "df_i                  : " << df_i << std::endl << std::endl;
    
    if (numThreads == 64)
    {
        if (numCiphertext == 32)
            numThreads = 59;
        if (numCiphertext == 64)
            numThreads = 62;
    }

    if (numCiphertext == 1)
    {
        if (singleThread)
        {
            #pragma omp parallel for
            for (size_t i = 0; i < 1; i++)
                testMedianAdv(vectorLength, dg_c, df_c, dg_i, df_i);
        }
        else
        {
            omp_set_num_threads(numThreads);
            if (numCiphertext <= numThreads / 16)
                omp_set_max_active_levels(10);
            testMedianAdv(vectorLength, dg_c, df_c, dg_i, df_i);
        }
    }
    else
    {
        if (singleThread)
        {
            #pragma omp parallel for
            for (size_t i = 0; i < 1; i++)
                testMedianMultiCtxtAdv(subVectorLength, numCiphertext, dg_c, df_c, dg_i, df_i);
        }
        else
        {
            omp_set_num_threads(numThreads);
            if (numCiphertext <= numThreads / 16)
                omp_set_max_active_levels(10);
            testMedianMultiCtxtAdv(subVectorLength, numCiphertext, dg_c, df_c, dg_i, df_i);
        }
    }

    return 0;

}
