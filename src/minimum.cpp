#include "minimum.h"
#include "ranking.h"
#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-matrices.h"
#include "utils-ptxt.h"

#include <cassert>
#include <omp.h>


double min(
    const std::vector<double> &vec
)
{
    auto minIter = std::min_element(vec.begin(), vec.end());
    return *minIter;
}


Ciphertext<DCRTPoly> min(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    const double leftBoundC,
    const double rightBoundC,
    const uint32_t degreeC,
    const uint32_t degreeI
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
    const double leftBoundC,
    const double rightBoundC,
    const uint32_t degreeC,
    const uint32_t degreeI
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
    std::cout << "COMPLETED (" << elapsed_seconds.count() << "s)" << std::endl << std::endl;

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
    std::cout << "COMPLETED (" << elapsed_seconds.count() << "s)" << std::endl << std::endl;

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
    std::cout << "COMPLETED (" << elapsed_seconds.count() << "s)" << std::endl << std::endl;

    return m;

}


Ciphertext<DCRTPoly> minFG(
    const std::vector<Ciphertext<DCRTPoly>> &c,
    const size_t subVectorLength,
    const uint32_t dg_c,
    const uint32_t df_c,
    const uint32_t dg_i,
    const uint32_t df_i
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

    ranking = rankFG(
        c,
        subVectorLength,
        dg_c, df_c,
        true, false
    );

    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "COMPLETED (" << elapsed_seconds.count() << "s)" << std::endl << std::endl;

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
    std::cout << "COMPLETED (" << elapsed_seconds.count() << "s)" << std::endl << std::endl;

    std::cout << "===================================\n";
    std::cout << "Indicator\n";
    std::cout << "===================================\n";

    Ciphertext<DCRTPoly> m;

    start = std::chrono::high_resolution_clock::now();

    m = indicatorAdv(
        s - 1.0,
        vectorLength,
        dg_i, df_i
    );

    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "COMPLETED (" << elapsed_seconds.count() << "s)" << std::endl << std::endl;

    return m;

}
