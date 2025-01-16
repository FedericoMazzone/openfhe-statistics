#include "ranking.h"
#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-matrices.h"
#include "utils-ptxt.h"

#include <omp.h>


void testRanking(
    const size_t vectorLength,
    const usint compareDepth,
    const bool tieCorrection = false
)
{

    const usint integralPrecision       = 1;
    const usint decimalPrecision        = (vectorLength <= 16 && !tieCorrection) ? 30 : 35;
    const usint multiplicativeDepth     = compareDepth + 1 + (tieCorrection ? 3 : 0);
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

    Ciphertext<DCRTPoly> vC = cryptoContext->Encrypt(
        keyPair.publicKey,
        cryptoContext->MakeCKKSPackedPlaintext(v)
    );

    auto start = std::chrono::high_resolution_clock::now();

    Ciphertext<DCRTPoly> resultC;
    if (tieCorrection)
        resultC = rankWithCorrection(
            vC,
            vectorLength,
            -1.0, 1.0,
            depth2degree(compareDepth),
            false
        );
    else
        resultC = rank(
            vC,
            vectorLength,
            -1.0, 1.0,
            depth2degree(compareDepth)
        );

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Runtime: " << elapsed_seconds.count() << "s" << std::endl;

    Plaintext resultP;
    cryptoContext->Decrypt(keyPair.secretKey, resultC, &resultP);
    resultP->SetLength(vectorLength);
    std::vector<double> result = resultP->GetRealPackedValue();
    std::cout << "Computed ranking: " << result << std::endl;

    std::vector<double> roundedResult(result.size());
    for (size_t i = 0; i < result.size(); i++)
        if (tieCorrection)
            roundedResult[i] = std::round(result[i]);
        else
            roundedResult[i] = std::round(result[i] * 2) / 2;
    std::cout << "Rounded ranking:  " << roundedResult << std::endl;

    std::vector<double> expectedRanking = rank(v, !tieCorrection);
    std::cout << "Expected ranking: " << rank(v, !tieCorrection) << std::endl;

    int numErrors = 0;
    for (size_t i = 0; i < result.size(); i++)
        if (roundedResult[i] != expectedRanking[i])
            numErrors++;
    if (numErrors > 0)
        std::cout << "Number of errors: " << numErrors << std::endl;
    else
        std::cout << "No errors" << std::endl;

}


void testRankingMultiCtxt(
    const size_t subVectorLength,
    const size_t numCiphertext,
    const usint compareDepth,
    const bool tieCorrection = false
)
{

    const size_t vectorLength           = subVectorLength * numCiphertext;
    const usint integralPrecision       = 1;
    const usint decimalPrecision        = 35;
    const usint multiplicativeDepth     = compareDepth + 2 + (tieCorrection ? 3 : 0);
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

    std::vector<Ciphertext<DCRTPoly>> vC(numCiphertext);
    for (size_t j = 0; j < numCiphertext; j++)
        vC[j] = cryptoContext->Encrypt(
            keyPair.publicKey,
            cryptoContext->MakeCKKSPackedPlaintext(vTokens[j])
        );

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<Ciphertext<DCRTPoly>> resultC;
    if (tieCorrection)
        resultC = rankWithCorrection(
            vC,
            subVectorLength,
            -1.0, 1.0,
            depth2degree(compareDepth)
        );
    else
        resultC = rank(
            vC,
            subVectorLength,
            -1.0, 1.0,
            depth2degree(compareDepth)
        );

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Runtime: " << elapsed_seconds.count() << "s" << std::endl;

    Plaintext resultP;
    std::vector<std::vector<double>> resultTokens(numCiphertext);
    for (size_t i = 0; i < numCiphertext; i++)
    {
        cryptoContext->Decrypt(keyPair.secretKey, resultC[i], &resultP);
        resultP->SetLength(subVectorLength);
        resultTokens[i] = resultP->GetRealPackedValue();
    }
    std::vector<double> result = concatVectors(resultTokens);
    std::cout << "Computed ranking: " << result << std::endl;

    std::vector<double> roundedResult(result.size());
    for (size_t i = 0; i < result.size(); i++)
        if (tieCorrection)
            roundedResult[i] = std::round(result[i]);
        else
            roundedResult[i] = std::round(result[i] * 2) / 2;
    std::cout << "Rounded ranking:  " << roundedResult << std::endl;

    std::vector<double> expectedRanking = rank(v, !tieCorrection);
    std::cout << "Expected ranking: " << rank(v, !tieCorrection) << std::endl;

    int numErrors = 0;
    for (size_t i = 0; i < result.size(); i++)
        if (roundedResult[i] != expectedRanking[i])
            numErrors++;
    if (numErrors > 0)
        std::cout << "Number of errors: " << numErrors << std::endl;
    else
        std::cout << "No errors" << std::endl;

}


void testRankingMultiCtxtAdv(
    const size_t subVectorLength,
    const size_t numCiphertext,
    const usint dg,
    const usint df,
    const bool tieCorrection = false
)
{

    const size_t vectorLength           = subVectorLength * numCiphertext;
    const usint integralPrecision       = 1;
    const usint decimalPrecision        = 45;
    const usint multiplicativeDepth     = 4 * (dg + df) + 3 + 1 + (tieCorrection ? 3 : 0);
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

    std::cout << "Vector:           " << vTokens << std::endl;

    std::vector<Ciphertext<DCRTPoly>> vC(numCiphertext);
    for (size_t j = 0; j < numCiphertext; j++)
        vC[j] = cryptoContext->Encrypt(
            keyPair.publicKey,
            cryptoContext->MakeCKKSPackedPlaintext(vTokens[j])
        );

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<Ciphertext<DCRTPoly>> resultC;
    if (tieCorrection)
        resultC = rankWithCorrectionFG(
            vC,
            subVectorLength,
            dg, df
        );
    else
        resultC = rankFG(
            vC,
            subVectorLength,
            dg, df
        );

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Runtime: " << elapsed_seconds.count() << "s" << std::endl;

    std::cout << "Result's level: " << resultC[0]->GetLevel() << std::endl;

    Plaintext resultP;
    std::vector<std::vector<double>> resultTokens(numCiphertext);
    for (size_t i = 0; i < numCiphertext; i++)
    {
        cryptoContext->Decrypt(keyPair.secretKey, resultC[i], &resultP);
        resultP->SetLength(subVectorLength);
        resultTokens[i] = resultP->GetRealPackedValue();
    }
    std::vector<double> result = concatVectors(resultTokens);
    std::cout << "Computed ranking: " << result << std::endl;

    std::vector<double> roundedResult(result.size());
    for (size_t i = 0; i < result.size(); i++)
        if (tieCorrection)
            roundedResult[i] = std::round(result[i]);
        else
            roundedResult[i] = std::round(result[i] * 2) / 2;
    std::cout << "Rounded ranking:  " << roundedResult << std::endl;

    std::vector<double> expectedRanking = rank(v, !tieCorrection);
    std::cout << "Expected ranking: " << rank(v, !tieCorrection) << std::endl;

    int numErrors = 0;
    for (size_t i = 0; i < result.size(); i++)
        if (roundedResult[i] != expectedRanking[i])
            numErrors++;
    if (numErrors > 0)
        std::cout << "Number of errors: " << numErrors << std::endl;
    else
        std::cout << "No errors" << std::endl;

}


int main(int argc, char *argv[])
{

    const size_t vectorLength = std::stoul(argv[1]);
    const bool tieCorrection = (argc > 2) ? (bool) std::stoi(argv[2]) : false;
    const bool singleThread = (argc > 3) ? (bool) std::stoi(argv[3]) : false;

    std::cout << "Vector length         : " << vectorLength << std::endl;
    std::cout << "Tie correction        : " << (tieCorrection ? "true" : "false") << std::endl;
    std::cout << "Single thread         : " << (singleThread ? "true" : "false") << std::endl;

    std::cout << std::fixed << std::setprecision(2);
    const size_t numThreads = std::thread::hardware_concurrency();
    std::cout << "Number of threads     : " << numThreads << std::endl;

    if (!tieCorrection)
        if (vectorLength <= 256)
        {
            const size_t subVectorLength = 128;
            const size_t numCiphertext = std::ceil((double) vectorLength / subVectorLength);

            usint compareDepth = 0;
            if      (vectorLength <= 8)     compareDepth = 7;
            else if (vectorLength <= 16)    compareDepth = 8;
            else if (vectorLength <= 64)    compareDepth = 10;
            else if (vectorLength <= 256)   compareDepth = 11;

            std::cout << "Subvector length      : " << subVectorLength << std::endl;
            std::cout << "Number of ciphertexts : " << numCiphertext << std::endl;
            std::cout << "Compare depth         : " << compareDepth << std::endl << std::endl;

            if (numCiphertext == 1)
            {
                if (singleThread)
                {
                    #pragma omp parallel for
                    for (size_t i = 0; i < 1; i++)
                        testRanking(vectorLength, compareDepth, tieCorrection);
                }
                else
                {
                    testRanking(vectorLength, compareDepth, tieCorrection);
                }
            }
            else
            {
                if (singleThread)
                {
                    #pragma omp parallel for
                    for (size_t i = 0; i < 1; i++)
                        testRankingMultiCtxt(subVectorLength, numCiphertext, compareDepth, tieCorrection);
                }
                else
                {
                    omp_set_num_threads(numThreads);
                    if (numCiphertext <= numThreads / 16)
                        omp_set_max_active_levels(10);
                    testRankingMultiCtxt(subVectorLength, numCiphertext, compareDepth, tieCorrection);
                }
            }
        }
        else
        {
            const size_t subVectorLength = 128;
            const size_t numCiphertext = std::ceil((double) vectorLength / subVectorLength);

            const usint dg = 3;
            const usint df = 2;

            std::cout << "Subvector length      : " << subVectorLength << std::endl;
            std::cout << "Number of ciphertexts : " << numCiphertext << std::endl;
            std::cout << "dg                    : " << dg << std::endl;
            std::cout << "df                    : " << df << std::endl << std::endl;

            if (singleThread)
            {
                #pragma omp parallel for
                for (size_t i = 0; i < 1; i++)
                    testRankingMultiCtxtAdv(subVectorLength, numCiphertext, dg, df, tieCorrection);
            }
            else
            {
                omp_set_num_threads(numThreads);
                if (numCiphertext <= numThreads / 16)
                    omp_set_max_active_levels(10);
                testRankingMultiCtxtAdv(subVectorLength, numCiphertext, dg, df, tieCorrection);
            }
        }
    else
        if (vectorLength <= 256)
        {
            const size_t subVectorLength = 128;
            const size_t numCiphertext = std::ceil((double) vectorLength / subVectorLength);

            usint compareDepth = 0;
            if      (vectorLength <= 8)     compareDepth = 7;
            else if (vectorLength <= 16)    compareDepth = 9;
            else if (vectorLength <= 64)    compareDepth = 10;
            else if (vectorLength <= 256)   compareDepth = 12;

            std::cout << "Subvector length      : " << subVectorLength << std::endl;
            std::cout << "Number of ciphertexts : " << numCiphertext << std::endl;
            std::cout << "Compare depth         : " << compareDepth << std::endl << std::endl;

            if (numCiphertext == 1)
            {
                if (singleThread)
                {
                    #pragma omp parallel for
                    for (size_t i = 0; i < 1; i++)
                        testRanking(vectorLength, compareDepth, tieCorrection);
                }
                else
                {
                    testRanking(vectorLength, compareDepth, tieCorrection);
                }
            }
            else
            {
                if (singleThread)
                {
                    #pragma omp parallel for
                    for (size_t i = 0; i < 1; i++)
                        testRankingMultiCtxt(subVectorLength, numCiphertext, compareDepth, tieCorrection);
                }
                else
                {
                    omp_set_num_threads(numThreads);
                    if (numCiphertext <= numThreads / 16)
                        omp_set_max_active_levels(10);
                    testRankingMultiCtxt(subVectorLength, numCiphertext, compareDepth, tieCorrection);
                }
            }
        }
        else
        {
            const size_t subVectorLength = 128;
            const size_t numCiphertext = std::ceil((double) vectorLength / subVectorLength);

            const usint dg = 3;
            const usint df = 2;

            std::cout << "Subvector length      : " << subVectorLength << std::endl;
            std::cout << "Number of ciphertexts : " << numCiphertext << std::endl;
            std::cout << "dg                    : " << dg << std::endl;
            std::cout << "df                    : " << df << std::endl << std::endl;

            if (singleThread)
            {
                #pragma omp parallel for
                for (size_t i = 0; i < 1; i++)
                    testRankingMultiCtxtAdv(subVectorLength, numCiphertext, dg, df, tieCorrection);
            }
            else
            {
                omp_set_num_threads(numThreads);
                if (numCiphertext <= numThreads / 16)
                    omp_set_max_active_levels(10);
                testRankingMultiCtxtAdv(subVectorLength, numCiphertext, dg, df, tieCorrection);
            }
        }

    return 0;

}
