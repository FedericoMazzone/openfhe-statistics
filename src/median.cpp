#include "median.h"
#include "ranking.h"
#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-matrices.h"
#include "utils-ptxt.h"

#include <cassert>
#include <omp.h>



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


Ciphertext<DCRTPoly> median(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    const double leftBoundC,
    const double rightBoundC,
    const uint32_t degreeC,
    const uint32_t degreeI
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


Ciphertext<DCRTPoly> medianFG(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    const uint32_t dg_c,
    const uint32_t df_c,
    const uint32_t dg_i,
    const uint32_t df_i
)
{
    c = rankWithCorrectionFG(
        c,
        vectorLength,
        dg_c,
        df_c
    );
    if (vectorLength % 2 == 0)
        c = indicatorAdvShifted(
            c,
            vectorLength,
            dg_i, df_i
        );
    else
        c = indicatorAdv(
            c - 0.5 * (vectorLength + 1),
            vectorLength,
            dg_i, df_i
        );

    return c;
}


Ciphertext<DCRTPoly> medianFG(
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

    ranking = rankWithCorrectionFG(
        c,
        subVectorLength,
        dg_c,
        df_c
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

    if (vectorLength % 2 == 0)
        m = indicatorAdvShifted(
            s,
            vectorLength,
            dg_i, df_i
        );
    else
        m = indicatorAdv(
            s - 0.5 * (vectorLength + 1),
            vectorLength,
            dg_i, df_i
        );

    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "COMPLETED (" << std::fixed << std::setprecision(3) <<
            elapsed_seconds.count() << "s)" << std::endl << std::endl;

    return m;

}
