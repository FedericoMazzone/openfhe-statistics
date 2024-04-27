#include "ranking.h"
#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-matrices.h"
#include "utils-ptxt.h"
#include <cassert>
#include <omp.h>


Ciphertext<DCRTPoly> rank(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    double leftBoundC,
    double rightBoundC,
    uint32_t degreeC,
    bool cmpGt
)
{
    if (!cmpGt)
    {
        c = compare(
            replicateRow(c, vectorLength),
            replicateColumn(transposeRow(c, vectorLength, true), vectorLength),
            leftBoundC, rightBoundC, degreeC
        );
    }
    else
    {
        c = compareGt(
            replicateRow(c, vectorLength),
            replicateColumn(transposeRow(c, vectorLength, true), vectorLength),
            leftBoundC, rightBoundC, degreeC,
            0.01
        );
    }

    c = sumRows(c, vectorLength);
    c = c + (!cmpGt ? 0.5 : 1.0);

    return c;
}


std::vector<Ciphertext<DCRTPoly>> rank(
    const std::vector<Ciphertext<DCRTPoly>> &c,
    const size_t subVectorLength,
    double leftBoundC,
    double rightBoundC,
    uint32_t degreeC,
    bool cmpGt,
    bool complOpt
)
{
    const size_t numCiphertext = c.size();

    static std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    std::chrono::duration<double> elapsed_seconds;

    std::cout << "===================================\n";
    std::cout << "Replicate\n";
    std::cout << "===================================\n";

    std::vector<Ciphertext<DCRTPoly>> replR(numCiphertext);
    std::vector<Ciphertext<DCRTPoly>> replC(numCiphertext);

    start = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for collapse(2)
    for (size_t loopID = 0; loopID < 2; loopID++)
    {
        for (size_t j = 0; j < numCiphertext; j++)
        {
            if (loopID == 0)
            {
                #pragma omp critical
                {std::cout << "ReplicateRow - j = " << j << " - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

                replR[j] = replicateRow(c[j], subVectorLength);

                #pragma omp critical
                {std::cout << "ReplicateRow - j = " << j << " - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
            }
            else
            {
                #pragma omp critical
                {std::cout << "ReplicateColumn - j = " << j << " - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

                replC[j] = replicateColumn(transposeRow(c[j], subVectorLength, true), subVectorLength);

                #pragma omp critical
                {std::cout << "ReplicateColumn - j = " << j << " - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
            }
        }
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "COMPLETED (" << std::fixed << std::setprecision(3) <<
            elapsed_seconds.count() << "s)" << std::endl << std::endl;
    
    if (!complOpt)
    {
        std::cout << "===================================\n";
        std::cout << "Compare\n";
        std::cout << "===================================\n";

        std::vector<Ciphertext<DCRTPoly>> C(numCiphertext);
        std::vector<bool> Cinitialized(numCiphertext, false);

        start = std::chrono::high_resolution_clock::now();

        #pragma omp parallel for collapse(2)
        for (size_t j = 0; j < numCiphertext; j++)
        {
            for (size_t k = 0; k < numCiphertext; k++)
            {
                #pragma omp critical
                {std::cout << "(j, k) = (" << j << ", " << k << ") - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

                Ciphertext<DCRTPoly> Cjk;
                if (!cmpGt)
                {
                    Cjk = compare(
                        replR[j],
                        replC[k],
                        leftBoundC, rightBoundC, degreeC
                    );
                }
                else
                {
                    Cjk = compareGt(
                        replR[j],
                        replC[k],
                        leftBoundC, rightBoundC, degreeC,
                        0.01
                    );
                }

                #pragma omp critical
                {
                if (!Cinitialized[j]) { C[j] = Cjk; Cinitialized[j] = true; }
                else                  { C[j] = C[j] + Cjk;                  }
                }
                
                #pragma omp critical
                {std::cout << "(j, k) = (" << j << ", " << k << ") - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
            }
        }

        end = std::chrono::high_resolution_clock::now();
        elapsed_seconds = end - start;
        std::cout << "COMPLETED (" << std::fixed << std::setprecision(3) <<
                elapsed_seconds.count() << "s)" << std::endl << std::endl;

        std::cout << "===================================\n";
        std::cout << "Sum\n";
        std::cout << "===================================\n";

        std::vector<Ciphertext<DCRTPoly>> s(numCiphertext);

        start = std::chrono::high_resolution_clock::now();

        #pragma omp parallel for
        for (size_t j = 0; j < numCiphertext; j++)
        {
            #pragma omp critical
            {std::cout << "SumRows - j = " << j << " - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

            s[j] = sumRows(C[j], subVectorLength) + (!cmpGt ? 0.5 : 1.0);
            
            #pragma omp critical
            {std::cout << "SumRows - j = " << j << " - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
        }

        end = std::chrono::high_resolution_clock::now();
        elapsed_seconds = end - start;
        std::cout << "COMPLETED (" << std::fixed << std::setprecision(3) <<
                elapsed_seconds.count() << "s)" << std::endl << std::endl;
        
        return s;
    }
    else
    {
        std::cout << "===================================\n";
        std::cout << "Compare\n";
        std::cout << "===================================\n";

        std::vector<Ciphertext<DCRTPoly>> Cv(numCiphertext);
        std::vector<Ciphertext<DCRTPoly>> Ch(numCiphertext);
        std::vector<bool> Cvinitialized(numCiphertext, false);
        std::vector<bool> Chinitialized(numCiphertext, false);

        start = std::chrono::high_resolution_clock::now();

        const size_t numReqThreads = numCiphertext * (numCiphertext + 1) / 2;
        std::cout << "Number of required threads: " << numReqThreads << std::endl;

        // for (size_t j = 0; j < numCiphertext; j++)
        // {
        //     for (size_t k = j; k < numCiphertext; k++)
        //     {
        // Collapse(2) with two nested for-loops creates issues here.
        #pragma omp parallel for
        for (size_t i = 0; i < numReqThreads; i++)
        {
            // Computing the indeces
            size_t j, k, counter = 0;
            bool loopCond = true;
            for (j = 0; j < numCiphertext && loopCond; j++)
                for (k = j; k < numCiphertext && loopCond; k++)
                    if (counter++ == i) loopCond = false;
            j--; k--;

            #pragma omp critical
            {std::cout << "(j, k) = (" << j << ", " << k << ") - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

            Ciphertext<DCRTPoly> Cjk;
            if (!cmpGt)
            {
                Cjk = compare(
                    replR[j],
                    replC[k],
                    leftBoundC, rightBoundC, degreeC
                );
            }
            else
            {
                Cjk = compareGt(
                    replR[j],
                    replC[k],
                    leftBoundC, rightBoundC, degreeC,
                    0.001
                );
            }

            #pragma omp critical
            {
            if (!Cvinitialized[j]) { Cv[j] = Cjk; Cvinitialized[j] = true; }
            else                   { Cv[j] = Cv[j] + Cjk;                  }
            }

            if (j != k)
            {
                Ciphertext<DCRTPoly> Ckj = 1.0 - Cjk;

                #pragma omp critical
                {
                if (!Chinitialized[k]) { Ch[k] = Ckj; Chinitialized[k] = true; }
                else                   { Ch[k] = Ch[k] + Ckj;                  }
                }
            }
            
            #pragma omp critical
            {std::cout << "(j, k) = (" << j << ", " << k << ") - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
        }

        end = std::chrono::high_resolution_clock::now();
        elapsed_seconds = end - start;
        std::cout << "COMPLETED (" << std::fixed << std::setprecision(3) <<
                elapsed_seconds.count() << "s)" << std::endl << std::endl;

        std::cout << "===================================\n";
        std::cout << "Sum\n";
        std::cout << "===================================\n";

        std::vector<Ciphertext<DCRTPoly>> sv(numCiphertext);
        std::vector<Ciphertext<DCRTPoly>> sh(numCiphertext);
        std::vector<Ciphertext<DCRTPoly>> s(numCiphertext);
        std::vector<bool> sinitialized(numCiphertext, false);

        start = std::chrono::high_resolution_clock::now();

        #pragma omp parallel for collapse(2)
        for (size_t loopID = 0; loopID < 2; loopID++)
        {
            for (size_t j = 0; j < numCiphertext; j++)
            {
                if (loopID == 0)
                {
                    #pragma omp critical
                    {std::cout << "SumV - j = " << j << " - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

                    sv[j] = sumRows(Cv[j], subVectorLength, true);
                    
                    #pragma omp critical
                    {
                    if (!sinitialized[j]) { s[j] = sv[j] + (!cmpGt ? 0.5 : 1.0); sinitialized[j] = true; }
                    else                  { s[j] = s[j] + sv[j];                                             }
                    }

                    #pragma omp critical
                    {std::cout << "SumV - j = " << j << " - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
                }
                else
                {
                    #pragma omp critical
                    {std::cout << "SumH - j = " << j << " - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

                    if (j > 0)
                    {
                        sh[j] = sumColumns(Ch[j], subVectorLength, true);
                        sh[j] = transposeColumn(sh[j], subVectorLength);

                        #pragma omp critical
                        {
                        if (!sinitialized[j]) { s[j] = sh[j] + (!cmpGt ? 0.5 : 1.0); sinitialized[j] = true; }
                        else                  { s[j] = s[j] + sh[j];                                             }
                        }
                    }

                    #pragma omp critical
                    {std::cout << "SumH - j = " << j << " - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
                }
            }
        }

        end = std::chrono::high_resolution_clock::now();
        elapsed_seconds = end - start;
        std::cout << "COMPLETED (" << std::fixed << std::setprecision(3) <<
                elapsed_seconds.count() << "s)" << std::endl << std::endl;
        
        return s;
    }
    
}


// fractional rank
std::vector<double> rank(
    const std::vector<double> &vec,
    double epsilon
)
{
    std::vector<double> ranking(vec.size(), 0.5);
    for (size_t i = 0; i < vec.size(); i++)
        for (size_t j = 0; j < vec.size(); j++)
        { 
            if (std::abs(vec[i] - vec[j]) <= epsilon)
                ranking[i] += 0.5;
            else if (vec[i] > vec[j])
                ranking[i] += 1.0;
        }

    return ranking;
}


double evaluateRankingAvg(
    const std::vector<double> &vec,
    const std::vector<double> &computedRanking
)
{
    std::vector<double> ranking = rank(vec);
    double error = 0.0;
    for (size_t i = 0; i < vec.size(); ++i)
        error += std::abs(computedRanking[i] - ranking[i]);

    return error / vec.size();
}


double evaluateRankingMax(
    const std::vector<double> &vec,
    const std::vector<double> &computedRanking
)
{
    std::vector<double> ranking = rank(vec);
    double error = 0.0;
    for (size_t i = 0; i < vec.size(); ++i)
        error = std::max(error, std::abs(computedRanking[i] - ranking[i]));

    return error;
}


double evaluateRankingRoundAvg(
    const std::vector<double> &vec,
    const std::vector<double> &computedRanking
)
{
    std::vector<double> ranking = rank(vec);
    double error = 0.0;
    for (size_t i = 0; i < vec.size(); ++i)
        error += std::abs(std::round(computedRanking[i] * 2) / 2 - ranking[i]);

    return error / vec.size();
}


double evaluateRankingRoundMax(
    const std::vector<double> &vec,
    const std::vector<double> &computedRanking
)
{
    std::vector<double> ranking = rank(vec);
    double error = 0.0;
    for (size_t i = 0; i < vec.size(); ++i)
        error = std::max(error, std::abs(std::round(computedRanking[i] * 2) / 2 - ranking[i]));

    return error;
}


bool evaluateRankingOrder(
    const std::vector<double> &vec,
    const std::vector<double> &computedRanking
)
{
    std::vector<double> ranking = rank(vec);
    std::vector<double> rr1 = rank(ranking);
    std::vector<double> rr2 = rank(computedRanking, 0.001);

    return std::equal(rr1.begin(), rr1.end(), rr2.begin());
}


double evaluateRankingWeighted(
    const std::vector<double> &vec,
    const std::vector<double> &computedRanking
)
{
    std::vector<double> ranking = rank(vec);

    std::vector<double> sortedVec = vec;
    std::sort(sortedVec.begin(), sortedVec.end());
    double rankingError, proximity, weightedError = 0.0;
    for (size_t i = 0; i < vec.size(); i++)
    {
        rankingError = std::abs(std::round(computedRanking[i] * 2) / 2 - (double) ranking[i]);

        proximity = 1.7e+308;
        for (size_t j = 0; j < vec.size(); j++)
            if (j != i)
                proximity = std::min(proximity, std::abs(vec[i] - vec[j]));
        
        weightedError += rankingError * proximity * proximity;
    }

    return weightedError / vec.size();
}


double evaluateRankingWeighted2(
    const std::vector<double> &vec,
    const std::vector<double> &computedRanking
)
{
    std::vector<double> ranking = rank(vec);

    std::vector<double> sortedVec = vec;
    std::sort(sortedVec.begin(), sortedVec.end());
    double rankingError, proximity, weightedError = 0.0;
    for (size_t i = 0; i < vec.size(); i++)
    {
        rankingError = std::round(computedRanking[i] * 2) / 2 - (double) ranking[i];

        proximity = 1.7e+308;
        for (size_t j = 0; j < vec.size(); j++)
            if (j != i)
                proximity = std::min(proximity, std::abs(vec[i] - vec[j]));
        
        weightedError += rankingError * rankingError * proximity * proximity;
    }

    return weightedError / vec.size();
}


std::vector<double> testRanking(
    const size_t vectorLength = 8,
    const usint compareDepth = 11
)
{

    std::cout << "Vector length: " << vectorLength << std::endl;
    std::cout << "Compare depth: " << compareDepth << std::endl;

    const usint integralPrecision       = 10;
    const usint decimalPrecision        = 50;
    const usint multiplicativeDepth     = compareDepth;
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
        v[i] = (double) rand() / RAND_MAX;
        // v[i] = (double) (i + 1) / (vectorLength + 1);
    // std::vector<double> v = {0.7, 0.2, 0.2, 0.8, 0.2, 0.5, 0.1, 0.3};

    std::cout << "Vector: " << v << std::endl;

    std::cout << "Expected ranking: " << rank(v) << std::endl;

    Ciphertext<DCRTPoly> vC = cryptoContext->Encrypt(
        keyPair.publicKey,
        cryptoContext->MakeCKKSPackedPlaintext(v)
    );

    auto start = std::chrono::high_resolution_clock::now();

    Ciphertext<DCRTPoly> resultC = rank(
        vC,
        vectorLength,
        -1.0, 1.0,
        depth2degree(compareDepth)
    );

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << elapsed_seconds.count() << "s" << std::endl;

    Plaintext resultP;
    cryptoContext->Decrypt(keyPair.secretKey, resultC, &resultP);
    resultP->SetLength(vectorLength);
    std::vector<double> result = resultP->GetRealPackedValue();
    std::cout << "Ranking: " << result << std::endl;

    double errorAvg = evaluateRankingAvg(v, result);
    double errorMax = evaluateRankingMax(v, result);
    double errorRoundAvg = evaluateRankingRoundAvg(v, result);
    double errorRoundMax = evaluateRankingRoundMax(v, result);
    double errorOrder = evaluateRankingOrder(v, result);
    double errorWeighted = evaluateRankingWeighted(v, result);
    double errorWeighted2 = evaluateRankingWeighted2(v, result);
    std::cout << "Error average: " << errorAvg << std::endl;
    std::cout << "Error max: " << errorMax << std::endl;
    std::cout << "Error round average: " << errorRoundAvg << std::endl;
    std::cout << "Error round max: " << errorRoundMax << std::endl;
    std::cout << "Error order: " << errorOrder << std::endl;
    std::cout << "Error weighted: " << std::setprecision(10) << errorWeighted << std::endl;
    std::cout << "Error weighted2: " << std::setprecision(10) << errorWeighted2 << std::endl;

    return {
        (double) vectorLength, (double) compareDepth, elapsed_seconds.count(),
        errorAvg, errorMax, errorRoundAvg, errorRoundMax, errorOrder, errorWeighted, errorWeighted2
    };

}


std::vector<double> testRankingMultiCtxt(
    const size_t subVectorLength = 128,
    const size_t numCiphertext = 2,
    const usint compareDepth = 14
)
{

    std::cout << "SubVector length: " << subVectorLength << std::endl;
    std::cout << "Number of ciphertexts: " << numCiphertext << std::endl;
    std::cout << "Compare depth: " << compareDepth << std::endl;

    const size_t vectorLength           = subVectorLength * numCiphertext;
    const usint integralPrecision       = 12;
    const usint decimalPrecision        = 48;
    const usint multiplicativeDepth     = compareDepth + 1;
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

    std::vector<double> v(vectorLength);
    for (size_t i = 0; i < vectorLength; i++)
        v[i] = (double) rand() / RAND_MAX;
        // v[i] = (double) (i + 1) / (vectorLength + 1);
    
    std::vector<std::vector<double>> vTokens = splitVector(v, numCiphertext);

    std::cout << "Vector: " << vTokens << std::endl;

    std::cout << "Expected ranking: " << rank(v) << std::endl;

    std::vector<Ciphertext<DCRTPoly>> vC(numCiphertext);
    for (size_t j = 0; j < numCiphertext; j++)
        vC[j] = cryptoContext->Encrypt(
            keyPair.publicKey,
            cryptoContext->MakeCKKSPackedPlaintext(vTokens[j])
        );

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<Ciphertext<DCRTPoly>> resultC = rank(
        vC,
        subVectorLength,
        -1.0, 1.0,
        depth2degree(compareDepth)
    );

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << elapsed_seconds.count() << "s" << std::endl;

    Plaintext resultP;
    std::vector<std::vector<double>> resultTokens(numCiphertext);
    for (size_t i = 0; i < numCiphertext; i++)
    {
        cryptoContext->Decrypt(keyPair.secretKey, resultC[i], &resultP);
        resultP->SetLength(subVectorLength);
        resultTokens[i] = resultP->GetRealPackedValue();
    }
    std::vector<double> result = concatVectors(resultTokens);
    std::cout << "Ranking: " << result << std::endl;

    double errorAvg = evaluateRankingAvg(v, result);
    double errorMax = evaluateRankingMax(v, result);
    double errorRoundAvg = evaluateRankingRoundAvg(v, result);
    double errorRoundMax = evaluateRankingRoundMax(v, result);
    double errorOrder = evaluateRankingOrder(v, result);
    double errorWeighted = evaluateRankingWeighted(v, result);
    double errorWeighted2 = evaluateRankingWeighted2(v, result);
    std::cout << "Error average: " << errorAvg << std::endl;
    std::cout << "Error max: " << errorMax << std::endl;
    std::cout << "Error round average: " << errorRoundAvg << std::endl;
    std::cout << "Error round max: " << errorRoundMax << std::endl;
    std::cout << "Error order: " << errorOrder << std::endl;
    std::cout << "Error weighted: " << std::setprecision(10) << errorWeighted << std::endl;
    std::cout << "Error weighted2: " << std::setprecision(10) << errorWeighted2 << std::endl;

    return {
        (double) vectorLength, (double) compareDepth, elapsed_seconds.count(),
        errorAvg, errorMax, errorRoundAvg, errorRoundMax, errorOrder, errorWeighted, errorWeighted2
    };

}


// int main()
// {

//     srand(time(NULL));

//     const size_t TEST_REPS = 10;

//     std::ofstream logFile;
//     logFile.open("log-ranking.txt");

//     for (size_t vectorLength = 8; vectorLength <= 128; vectorLength *= 2)
//         for (usint compareDepth = 6; compareDepth <= 14; compareDepth++)
//         {
//             std::vector<std::vector<double>> results;
//             for (size_t i = 0; i < TEST_REPS; i++)
//                 results.push_back(testRanking(vectorLength, compareDepth));
//             std::vector<double> average = averageVectors(results);
//             std::cout << "************************" << std::endl;
//             std::cout << average << std::endl;
//             std::cout << "************************" << std::endl;
//             logFile << average << std::endl;
//         }
//     for (size_t numCiphertext = 2; numCiphertext <= 8; numCiphertext *= 2)
//         for (usint compareDepth = 6; compareDepth <= 14; compareDepth++)
//         {
//             std::vector<std::vector<double>> results;
//             for (size_t i = 0; i < TEST_REPS; i++)
//                 results.push_back(testRankingMultiCtxt(128, numCiphertext, compareDepth));
//             std::vector<double> average = averageVectors(results);
//             std::cout << "************************" << std::endl;
//             std::cout << average << std::endl;
//             std::cout << "************************" << std::endl;
//             logFile << average << std::endl;
//         }

//     logFile.close();

//     return 0;

// }
