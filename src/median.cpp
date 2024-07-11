#include "median.h"
#include "ranking.h"
#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-matrices.h"
#include "utils-ptxt.h"
#include <cassert>
#include <omp.h>


Ciphertext<DCRTPoly> median(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    double leftBoundC,
    double rightBoundC,
    uint32_t degreeC,
    uint32_t degreeI
)
{
    c = rankWithCorrection(
        c,
        vectorLength,
        leftBoundC,
        rightBoundC,
        degreeC
    );
    c = indicator(
        c,
        (vectorLength % 2 == 0) ? 0.5 * vectorLength - 0.5 : 0.5 * vectorLength,
        (vectorLength % 2 == 0) ? 0.5 * vectorLength + 1.5 : 0.5 * vectorLength + 1.0,
        0.5, vectorLength + 0.5,
        degreeI
    );

    return c;
}


Ciphertext<DCRTPoly> median(
    const std::vector<Ciphertext<DCRTPoly>> &c,
    const size_t subVectorLength,
    double leftBoundC,
    double rightBoundC,
    uint32_t degreeC,
    uint32_t degreeI
)
{
    const size_t numCiphertext = c.size();
    const size_t vectorLength = subVectorLength * numCiphertext;

    static std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    std::chrono::duration<double> elapsed_seconds;

    std::cout << "===================================\n";
    std::cout << "Ranking\n";
    std::cout << "===================================\n";

    std::vector<Ciphertext<DCRTPoly>> ranking(numCiphertext);

    start = std::chrono::high_resolution_clock::now();

    ranking = rankWithCorrection(
        c,
        subVectorLength,
        leftBoundC,
        rightBoundC,
        degreeC
    );

    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "COMPLETED (" << std::fixed << std::setprecision(3) <<
            elapsed_seconds.count() << "s)" << std::endl << std::endl;

    std::cout << "===================================\n";
    std::cout << "Merging\n";
    std::cout << "===================================\n";

    Ciphertext<DCRTPoly> s;
    bool sinitialized = false;

    start = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for
    for (size_t j = 0; j < numCiphertext; j++)
    {
        #pragma omp critical
        {std::cout << "Merging - j = " << j << " - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

        ranking[j] = maskRow(ranking[j], subVectorLength, 0);
        ranking[j] = ranking[j] >> (j * subVectorLength);

        #pragma omp critical
        {
        if (!sinitialized) { s = ranking[j];     sinitialized = true; }
        else               { s = s + ranking[j];                      }
        }

        #pragma omp critical
        {std::cout << "Merging - j = " << j << " - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "COMPLETED (" << std::fixed << std::setprecision(3) <<
            elapsed_seconds.count() << "s)" << std::endl << std::endl;

    std::cout << "===================================\n";
    std::cout << "Indicator\n";
    std::cout << "===================================\n";

    Ciphertext<DCRTPoly> m;

    start = std::chrono::high_resolution_clock::now();

    m = indicator(
        s,
        (vectorLength % 2 == 0) ? 0.5 * vectorLength - 0.5 : 0.5 * vectorLength,
        (vectorLength % 2 == 0) ? 0.5 * vectorLength + 1.5 : 0.5 * vectorLength + 1.0,
        -0.01 * vectorLength, 1.01 * vectorLength,
        degreeI
    );

    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "COMPLETED (" << std::fixed << std::setprecision(3) <<
            elapsed_seconds.count() << "s)" << std::endl << std::endl;

    return m;

}


double median(
    std::vector<double> vec
)
{
    std::sort(vec.begin(), vec.end());
    size_t size = vec.size();
    if (size % 2 == 0) 
        return (vec[size / 2 - 1] + vec[size / 2]) / 2.0;
    else
        return vec[size / 2];
}


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


std::vector<double> testMedian(
    const size_t vectorLength = 8,
    const usint compareDepth = 7,
    const usint indicatorDepth = 11
)
{

    std::cout << "Vector length: " << vectorLength << std::endl;
    std::cout << "Compare depth: " << compareDepth << std::endl;
    std::cout << "Indicator depth: " << indicatorDepth << std::endl;

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

    std::vector<double> v(vectorLength);
    for (size_t i = 0; i < vectorLength; i++)
        v[i] = (double) rand() / RAND_MAX / 2.0 + 0.5;
        // v[i] = (double) (i + 1) / (vectorLength + 1);

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
    std::cout << elapsed_seconds.count() << "s" << std::endl;

    Plaintext resultP;
    cryptoContext->Decrypt(keyPair.secretKey, resultC, &resultP);
    resultP->SetLength(vectorLength);
    std::vector<double> result = resultP->GetRealPackedValue();
    std::cout << "Median: " << result << std::endl;

    double errorValue = evaluateMedianValue(v, result);
    double errorMask = evaluateMedianMask(v, result);
    std::cout << "Error value: " << errorValue << std::endl;
    std::cout << "Error mask: " << errorMask << std::endl;

    return {
        (double) vectorLength, (double) compareDepth, (double) indicatorDepth,
        elapsed_seconds.count(), errorValue, errorMask
    };

}


std::vector<double> testMedianMultiCtxt(
    const size_t subVectorLength = 128,
    const size_t numCiphertext = 2,
    const usint compareDepth = 7,
    const usint indicatorDepth = 11
)
{

    std::cout << "SubVector length: " << subVectorLength << std::endl;
    std::cout << "Number of ciphertexts: " << numCiphertext << std::endl;
    std::cout << "Compare depth: " << compareDepth << std::endl;
    std::cout << "Indicator depth: " << indicatorDepth << std::endl;

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

    std::vector<double> v(vectorLength);
    for (size_t i = 0; i < vectorLength; i++)
        v[i] = (double) rand() / RAND_MAX / 2.0 + 0.5;
        // v[i] = (double) (i + 1) / (vectorLength + 1);
    
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
    std::cout << elapsed_seconds.count() << "s" << std::endl;

    Plaintext resultP;
    cryptoContext->Decrypt(keyPair.secretKey, resultC, &resultP);
    resultP->SetLength(vectorLength);
    std::vector<double> result = resultP->GetRealPackedValue();
    std::cout << "Median: " << result << std::endl;

    double errorValue = evaluateMedianValue(v, result);
    double errorMask = evaluateMedianMask(v, result);
    std::cout << "Error value: " << errorValue << std::endl;
    std::cout << "Error mask: " << errorMask << std::endl;

    return {
        (double) vectorLength, (double) compareDepth, (double) indicatorDepth,
        elapsed_seconds.count(), errorValue, errorMask
    };

}


// int main()
// {

//     srand(time(NULL));

//     const size_t TEST_REPS = 10;

//     std::ofstream logFile;
//     logFile.open("log-median.txt");

//     logFile << "vector length vs. run time" << std::endl;
//     usint compareDepth = 9;
//     usint indicatorDepth = 9;
//     logFile << "compareDepth=" << compareDepth << std::endl;
//     logFile << "indicatorDepth=" << indicatorDepth << std::endl;
//     for (size_t vectorLength = 8; vectorLength <= 128; vectorLength *= 2)
//     {
//         std::vector<std::vector<double>> results;
//         for (size_t i = 0; i < TEST_REPS; i++)
//             try { results.push_back(testMedian(vectorLength, compareDepth, indicatorDepth)); }
//             catch (const std::exception& e) { std::cout << "Exception caught: " << e.what() << std::endl; }
//         if (results.size() > 0)
//         {
//             std::vector<double> average = averageVectors(results);
//             std::cout << "************************" << std::endl;
//             std::cout << average << std::endl;
//             std::cout << "************************" << std::endl;
//             logFile << average << std::endl;
//         }
//     }
//     for (size_t numCiphertext = 2; numCiphertext <= 8; numCiphertext *= 2)
//     {
//         std::vector<std::vector<double>> results;
//         for (size_t i = 0; i < TEST_REPS; i++)
//             try { results.push_back(testMedianMultiCtxt(128, numCiphertext, compareDepth, indicatorDepth)); }
//             catch (const std::exception& e) { std::cout << "Exception caught: " << e.what() << std::endl; }
//         if (results.size() > 0)
//         {
//             std::vector<double> average = averageVectors(results);
//             std::cout << "************************" << std::endl;
//             std::cout << average << std::endl;
//             std::cout << "************************" << std::endl;
//             logFile << average << std::endl;
//         }
//     }

//     // logFile << "depths vs. run time & errors" << std::endl;
//     // size_t vectorLength = 32;
//     // logFile << "vectorLength=" << vectorLength << std::endl;
//     // for (usint compareDepth = 6; compareDepth <= 14; compareDepth++)
//     //     for (usint indicatorDepth = 6; indicatorDepth <= 14; indicatorDepth++)
//     //     {
//     //         std::vector<std::vector<double>> results;
//     //         for (size_t i = 0; i < TEST_REPS; i++)
//     //             try { results.push_back(testMedian(vectorLength, compareDepth, indicatorDepth)); }
//     //             catch (const std::exception& e) { std::cout << "Exception caught: " << e.what() << std::endl; }
//     //         if (results.size() > 0)
//     //         {
//     //             std::vector<double> average = averageVectors(results);
//     //             std::cout << "************************" << std::endl;
//     //             std::cout << average << std::endl;
//     //             std::cout << "************************" << std::endl;
//     //             logFile << average << std::endl;
//     //         }
//     //     }

//     logFile.close();

//     return 0;

// }
