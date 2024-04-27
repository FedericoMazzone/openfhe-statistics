#include "ranking.h"
#include "minimum.h"
#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-matrices.h"
#include "utils-ptxt.h"
#include <cassert>
#include <omp.h>


Ciphertext<DCRTPoly> min(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    double leftBoundC,
    double rightBoundC,
    uint32_t degreeC,
    uint32_t degreeI
)
{
    c = rank(
        c,
        vectorLength,
        leftBoundC,
        rightBoundC,
        degreeC,
        true
    );
    c = indicator(
        c,
        0.5, 1.5,
        0.5, vectorLength + 0.5,
        degreeI
    );

    return c;
}


Ciphertext<DCRTPoly> min(
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

    ranking = rank(
        c,
        subVectorLength,
        leftBoundC,
        rightBoundC,
        degreeC,
        true,
        false
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
        0.5, 1.5,
        -0.01 * vectorLength, 1.01 * vectorLength,
        degreeI
    );

    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "COMPLETED (" << std::fixed << std::setprecision(3) <<
            elapsed_seconds.count() << "s)" << std::endl << std::endl;

    return m;

}


double min(
    const std::vector<double> &vec
)
{
    auto minIter = std::min_element(vec.begin(), vec.end());
    return *minIter;
}


size_t argmax(
    const std::vector<double>& vec
)
{
    auto maxIter = std::max_element(vec.begin(), vec.end());
    return std::distance(vec.begin(), maxIter);
}


double evaluateMinValue(
    const std::vector<double> &vec,
    const std::vector<double> &computedMinMask
)
{
    assert(vec.size() == computedMinMask.size());

    double expectedMin = min(vec);

    double computedMin = 0.0;
    double norm = 0.0;
    for (size_t i = 0; i < vec.size(); i++)
    {
        computedMin += vec[i] * computedMinMask[i];
        norm += computedMinMask[i];
    }
    computedMin /= norm;

    return std::abs(expectedMin - computedMin) / expectedMin;
}


double evaluateMinMask(
    const std::vector<double> &vec,
    const std::vector<double> &computedMinMask
)
{
    assert(vec.size() == computedMinMask.size());

    double expectedMin = min(vec);
    double computedMin = vec[argmax(computedMinMask)];

    return std::abs(expectedMin - computedMin) / expectedMin;
}


std::vector<double> testMinimum(
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

    std::vector<double> v(vectorLength);
    for (size_t i = 0; i < vectorLength; i++)
        v[i] = (double) rand() / RAND_MAX / 2.0 + 0.5;
        // v[i] = (double) (i + 1) / (vectorLength + 1);

    std::cout << "Vector: " << v << std::endl;

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
    std::cout << elapsed_seconds.count() << "s" << std::endl;

    Plaintext resultP;
    cryptoContext->Decrypt(keyPair.secretKey, resultC, &resultP);
    resultP->SetLength(vectorLength);
    std::vector<double> result = resultP->GetRealPackedValue();
    std::cout << "Minimum: " << result << std::endl;

    double errorValue = evaluateMinValue(v, result);
    double errorMask = evaluateMinMask(v, result);
    std::cout << "Error value: " << errorValue << std::endl;
    std::cout << "Error mask: " << errorMask << std::endl;

    return {
        (double) vectorLength, (double) compareDepth, (double) indicatorDepth,
        elapsed_seconds.count(), errorValue, errorMask
    };

}


std::vector<double> testMinimumMultiCtxt(
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

    std::vector<double> v(vectorLength);
    for (size_t i = 0; i < vectorLength; i++)
        v[i] = (double) rand() / RAND_MAX / 2.0 + 0.5;
        // v[i] = (double) (i + 1) / (vectorLength + 1);
    
    std::vector<std::vector<double>> vTokens = splitVector(v, numCiphertext);

    std::cout << "Vector: " << vTokens << std::endl;

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
    std::cout << elapsed_seconds.count() << "s" << std::endl;

    Plaintext resultP;
    cryptoContext->Decrypt(keyPair.secretKey, resultC, &resultP);
    resultP->SetLength(vectorLength);
    std::vector<double> result = resultP->GetRealPackedValue();
    std::cout << "Minimum: " << result << std::endl;

    double errorValue = evaluateMinValue(v, result);
    double errorMask = evaluateMinMask(v, result);
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
//     logFile.open("log-minimum.txt");

//     logFile << "vector length vs. run time" << std::endl;
//     usint compareDepth = 9;
//     usint indicatorDepth = 9;
//     logFile << "compareDepth=" << compareDepth << std::endl;
//     logFile << "indicatorDepth=" << indicatorDepth << std::endl;
//     for (size_t vectorLength = 8; vectorLength <= 128; vectorLength *= 2)
//     {
//         std::vector<std::vector<double>> results;
//         for (size_t i = 0; i < TEST_REPS; i++)
//             try { results.push_back(testMinimum(vectorLength, compareDepth, indicatorDepth)); }
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
//             try { results.push_back(testMinimumMultiCtxt(128, numCiphertext, compareDepth, indicatorDepth)); }
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

//     logFile << "depths vs. run time & errors" << std::endl;
//     size_t vectorLength = 32;
//     logFile << "vectorLength=" << vectorLength << std::endl;
//     for (usint compareDepth = 13; compareDepth <= 14; compareDepth++)
//         for (usint indicatorDepth = (compareDepth == 13) ? 13 : 6; indicatorDepth <= 14; indicatorDepth++)
//         {
//             std::vector<std::vector<double>> results;
//             for (size_t i = 0; i < TEST_REPS; i++)
//                 try { results.push_back(testMinimum(vectorLength, compareDepth, indicatorDepth)); }
//                 catch (const std::exception& e) { std::cout << "Exception caught: " << e.what() << std::endl; }
//             if (results.size() > 0)
//             {
//                 std::vector<double> average = averageVectors(results);
//                 std::cout << "************************" << std::endl;
//                 std::cout << average << std::endl;
//                 std::cout << "************************" << std::endl;
//                 logFile << average << std::endl;
//             }
//         }

//     logFile.close();

//     return 0;

// }
