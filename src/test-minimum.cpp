#include "minimum.h"
#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-matrices.h"
#include "utils-ptxt.h"

#include <cassert>
#include <omp.h>


size_t argmax(
    const std::vector<double>& vec
)
{
    auto maxIter = std::max_element(vec.begin(), vec.end());
    return std::distance(vec.begin(), maxIter);
}


void testMinimum(
    const size_t vectorLength,
    const usint compareDepth,
    const usint indicatorDepth
)
{

    const usint integralPrecision       = 1;
    const usint decimalPrecision        = 59;
    const usint multiplicativeDepth     = compareDepth + indicatorDepth;
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

    std::cout << "Expected minimum: " << min(v) << std::endl;

    Ciphertext<DCRTPoly> vC = cryptoContext->Encrypt(
        keyPair.publicKey,
        cryptoContext->MakeCKKSPackedPlaintext(v)
    );

    auto start = std::chrono::high_resolution_clock::now();

    Ciphertext<DCRTPoly> resultC = min(
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
    std::cout << "Computed mask: " << result << std::endl;

    std::vector<double> minima;
    for (size_t i = 0; i < vectorLength; i++)
        if (result[i] > 0.5)
            minima.push_back(v[i]);
    std::cout << "Found minima: " << minima << std::endl;

}


void testMinimumMultiCtxt(
    const size_t subVectorLength,
    const size_t numCiphertext,
    const usint compareDepth,
    const usint indicatorDepth
)
{

    const size_t vectorLength           = subVectorLength * numCiphertext;
    const usint integralPrecision       = 1;
    const usint decimalPrecision        = 59;
    const usint multiplicativeDepth     = compareDepth + indicatorDepth + 2;
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

    std::cout << "Vector:           " << vTokens << std::endl;

    std::cout << "Expected minimum: " << min(v) << std::endl;

    std::vector<Ciphertext<DCRTPoly>> vC(numCiphertext);
    for (size_t j = 0; j < numCiphertext; j++)
        vC[j] = cryptoContext->Encrypt(
            keyPair.publicKey,
            cryptoContext->MakeCKKSPackedPlaintext(vTokens[j])
        );

    auto start = std::chrono::high_resolution_clock::now();

    Ciphertext<DCRTPoly> resultC = min(
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
    std::cout << "Computed mask: " << result << std::endl;

    std::vector<double> minima;
    for (size_t i = 0; i < vectorLength; i++)
        if (result[i] > 0.5)
            minima.push_back(v[i]);
    std::cout << "Found minima: " << minima << std::endl;

}


void testMinimumMultiCtxtAdv(
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
    const usint multiplicativeDepth     = 4 * (dg_c + df_c + dg_i + df_i) + 3 + 2;
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

    std::cout << "Vector:           " << vTokens << std::endl;

    std::cout << "Expected minimum: " << min(v) << std::endl;

    std::vector<Ciphertext<DCRTPoly>> vC(numCiphertext);
    for (size_t j = 0; j < numCiphertext; j++)
        vC[j] = cryptoContext->Encrypt(
            keyPair.publicKey,
            cryptoContext->MakeCKKSPackedPlaintext(vTokens[j])
        );

    auto start = std::chrono::high_resolution_clock::now();

    Ciphertext<DCRTPoly> resultC = minFG(
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
    std::cout << "Computed mask: " << result << std::endl;

    std::vector<double> minima;
    for (size_t i = 0; i < vectorLength; i++)
        if (result[i] > 0.5)
            minima.push_back(v[i]);
    std::cout << "Found minima: " << minima << std::endl;

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

    if (vectorLength <= 256)
    {
        const size_t subVectorLength = 128;
        const size_t numCiphertext = std::ceil((double) vectorLength / subVectorLength);

        usint compareDepth = 0;
        usint indicatorDepth = 0;
        if (vectorLength <= 32)
        {
            compareDepth = 7;
            indicatorDepth = 7;
        }
        else if (vectorLength <= 256)
        {
            compareDepth = 9;
            indicatorDepth = 7;
        }
        else if (vectorLength <= 512)
        {
            compareDepth = 10;
            indicatorDepth = 8;
        }
        else if (vectorLength <= 1024)
        {
            compareDepth = 11;
            indicatorDepth = 9;
        }

        std::cout << "Subvector length      : " << subVectorLength << std::endl;
        std::cout << "Number of ciphertexts : " << numCiphertext << std::endl;
        std::cout << "Compare depth         : " << compareDepth << std::endl << std::endl;

        if (numCiphertext == 1)
        {
            if (singleThread)
            {
                #pragma omp parallel for
                for (size_t i = 0; i < 1; i++)
                    testMinimum(vectorLength, compareDepth, indicatorDepth);
            }
            else
            {
                testMinimum(vectorLength, compareDepth, indicatorDepth);
            }
        }
        else
        {
            if (singleThread)
            {
                #pragma omp parallel for
                for (size_t i = 0; i < 1; i++)
                    testMinimumMultiCtxt(subVectorLength, numCiphertext, compareDepth, indicatorDepth);
            }
            else
            {
                omp_set_num_threads(numThreads);
                if (numCiphertext <= numThreads / 16)
                    omp_set_max_active_levels(10);
                testMinimumMultiCtxt(subVectorLength, numCiphertext, compareDepth, indicatorDepth);
            }
        }
    }
    else
    {
        const size_t subVectorLength = 256;
        const size_t numCiphertext = std::ceil((double) vectorLength / subVectorLength);

        const usint dg_c = 3;
        const usint df_c = 2;
        usint dg_i = 0;
        usint df_i = 0;
        if (vectorLength <= 2048)
        {
            dg_i = 4;
            df_i = 2;
        }
        else if (vectorLength <= 8192)
        {
            dg_i = 5;
            df_i = 2;
        }
        else if (vectorLength <= 16384)
        {
            dg_i = 6;
            df_i = 2;
        }

        std::cout << "Subvector length      : " << subVectorLength << std::endl;
        std::cout << "Number of ciphertexts : " << numCiphertext << std::endl;
        std::cout << "dg_c                  : " << dg_c << std::endl;
        std::cout << "df_c                  : " << df_c << std::endl;
        std::cout << "dg_i                  : " << dg_i << std::endl;
        std::cout << "df_i                  : " << df_i << std::endl << std::endl;

        if (numCiphertext > 1)
        {
            if (singleThread)
            {
                #pragma omp parallel for
                for (size_t i = 0; i < 1; i++)
                    testMinimumMultiCtxtAdv(subVectorLength, numCiphertext, dg_c, df_c, dg_i, df_i);
            }
            else
            {
                omp_set_num_threads(64);
                if (numCiphertext <= 4)
                    omp_set_max_active_levels(10);
                testMinimumMultiCtxtAdv(subVectorLength, numCiphertext, dg_c, df_c, dg_i, df_i);
            }
        }
    }
    
    return 0;

}
