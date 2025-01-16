#include "sorting.h"
#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-matrices.h"
#include "utils-ptxt.h"

#include <cassert>
#include <omp.h>



std::vector<double> sort(
    std::vector<double> vec
)
{
    std::sort(vec.begin(), vec.end());
    return vec;
}


Ciphertext<DCRTPoly> sort(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    double leftBoundC,
    double rightBoundC,
    uint32_t degreeC,
    uint32_t degreeI
)
{
    Ciphertext<DCRTPoly> VR = replicateRow(c, vectorLength);
    Ciphertext<DCRTPoly> VC = replicateColumn(transposeRow(c, vectorLength, true), vectorLength);
    Ciphertext<DCRTPoly> C = compare(
        VR, VC,
        leftBoundC, rightBoundC,
        degreeC
    );
    Ciphertext<DCRTPoly> R = sumRows(C, vectorLength);

    std::vector<double> subMask(vectorLength * vectorLength);
    for (size_t i = 0; i < vectorLength; i++)
        for (size_t j = 0; j < vectorLength; j++)
            subMask[i * vectorLength + j] = -1.0 * i - 0.5;
    Ciphertext<DCRTPoly> M = indicator(
        R + subMask,
        -0.5, 0.5,
        -1.0 * vectorLength, 1.0 * vectorLength,
        degreeI
    );

    Ciphertext<DCRTPoly> S = sumColumns(M * VR, vectorLength);

    return S;
}


Ciphertext<DCRTPoly> sortWithCorrection(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    double leftBoundC,
    double rightBoundC,
    uint32_t degreeC,
    uint32_t degreeI
)
{
    Ciphertext<DCRTPoly> VR = replicateRow(c, vectorLength);
    Ciphertext<DCRTPoly> VC = replicateColumn(transposeRow(c, vectorLength, true), vectorLength);
    Ciphertext<DCRTPoly> C = compare(
        VR, VC,
        leftBoundC, rightBoundC,
        degreeC
    );

    std::vector<double> triangularMask(vectorLength * vectorLength, 0.0);
    for (size_t i = 0; i < vectorLength; i++)
        for (size_t j = 0; j < vectorLength; j++)
            if (j >= i) triangularMask[i * vectorLength + j] = 1.0;
    Ciphertext<DCRTPoly> E = 4 * (1 - C) * C;
    Ciphertext<DCRTPoly> correctionOffset = sumRows(E * triangularMask, vectorLength) - 0.5 * sumRows(E, vectorLength);

    Ciphertext<DCRTPoly> R = sumRows(C, vectorLength) + correctionOffset;

    std::vector<double> subMask(vectorLength * vectorLength);
    for (size_t i = 0; i < vectorLength; i++)
        for (size_t j = 0; j < vectorLength; j++)
            subMask[i * vectorLength + j] = -1.0 * i - 1.0;
    Ciphertext<DCRTPoly> M = indicator(
        R + subMask,
        -0.5, 0.5,
        -1.0 * vectorLength, 1.0 * vectorLength,
        degreeI
    );

    Ciphertext<DCRTPoly> S = sumColumns(M * VR, vectorLength);

    return S;
}


std::vector<Ciphertext<DCRTPoly>> sort(
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
    std::cout << "COMPLETED (" <<
            elapsed_seconds.count() << "s)" << std::endl << std::endl;
    
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

        Ciphertext<DCRTPoly> Cjk = compare(
            replR[j],
            replC[k],
            leftBoundC, rightBoundC, degreeC
        );

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
    std::cout << "COMPLETED (" <<
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

                sv[j] = sumRows(Cv[j], subVectorLength);
                
                #pragma omp critical
                {
                if (!sinitialized[j]) { s[j] = sv[j]; sinitialized[j] = true; }
                else                  { s[j] = s[j] + sv[j];                  }
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
                    sh[j] = transposeColumn(sh[j], subVectorLength, true);
                    sh[j] = replicateRow(sh[j], subVectorLength);

                    #pragma omp critical
                    {
                    if (!sinitialized[j]) { s[j] = sh[j]; sinitialized[j] = true; }
                    else                  { s[j] = s[j] + sh[j];                  }
                    }
                }

                #pragma omp critical
                {std::cout << "SumH - j = " << j << " - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
            }
        }
    }

    // return {s[1]};

    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "COMPLETED (" <<
            elapsed_seconds.count() << "s)" << std::endl << std::endl;
    
    std::cout << "===================================\n";
    std::cout << "Sort\n";
    std::cout << "===================================\n";

    start = std::chrono::high_resolution_clock::now();

    std::vector<std::vector<double>> subMasks(numCiphertext);
    for (size_t i = 0; i < numCiphertext; i++)
    {
        std::vector<double> subMask(subVectorLength * subVectorLength);
        for (size_t j = 0; j < subVectorLength; j++)
            for (size_t k = 0; k < subVectorLength; k++)
                subMask[j * subVectorLength + k] = -1.0 * (i * subVectorLength + j) - 0.5;
        subMasks[i] = subMask;
    }

    std::vector<Ciphertext<DCRTPoly>> subSorted(numCiphertext);
    std::vector<bool> subSortedInitialized(numCiphertext, false);

    #pragma omp parallel for collapse(2)
    for (size_t j = 0; j < numCiphertext; j++)
    {
        for (size_t k = 0; k < numCiphertext; k++)
        {
            #pragma omp critical
            {std::cout << "(j, k) = (" << j << ", " << k << ") - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

            Ciphertext<DCRTPoly> ind = indicator(
                s[k] + subMasks[j],
                -0.5, 0.5,
                -1.01 * vectorLength, 1.01 * vectorLength,
                degreeI
            ) * replR[k];

            #pragma omp critical
            {
            if (!subSortedInitialized[j]) { subSorted[j] = ind; subSortedInitialized[j] = true; }
            else                          { subSorted[j] = subSorted[j] + ind;                  }
            }
            
            #pragma omp critical
            {std::cout << "(j, k) = (" << j << ", " << k << ") - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
        }
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "COMPLETED (" <<
            elapsed_seconds.count() << "s)" << std::endl << std::endl;
    
    std::cout << "===================================\n";
    std::cout << "Sum\n";
    std::cout << "===================================\n";

    start = std::chrono::high_resolution_clock::now();

    std::vector<Ciphertext<DCRTPoly>> result(numCiphertext);

    #pragma omp parallel for
    for (size_t j = 0; j < numCiphertext; j++)
    {
        #pragma omp critical
        {std::cout << "Sum - j = " << j << " - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

        result[j] = sumColumns(subSorted[j], subVectorLength);

        #pragma omp critical
        {std::cout << "Sum - j = " << j << " - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "COMPLETED (" <<
            elapsed_seconds.count() << "s)" << std::endl << std::endl;

    return result;

}


std::vector<Ciphertext<DCRTPoly>> sortWithCorrection(
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
    std::cout << "COMPLETED (" <<
            elapsed_seconds.count() << "s)" << std::endl << std::endl;
    
    std::cout << "===================================\n";
    std::cout << "Compare\n";
    std::cout << "===================================\n";

    std::vector<Ciphertext<DCRTPoly>> Cv(numCiphertext);
    std::vector<Ciphertext<DCRTPoly>> Ch(numCiphertext);
    std::vector<bool> Cvinitialized(numCiphertext, false);
    std::vector<bool> Chinitialized(numCiphertext, false);

    std::vector<Ciphertext<DCRTPoly>> Ev(numCiphertext);
    std::vector<Ciphertext<DCRTPoly>> Eh(numCiphertext);
    std::vector<bool> Evinitialized(numCiphertext, false);
    std::vector<bool> Ehinitialized(numCiphertext, false);

    std::vector<Ciphertext<DCRTPoly>> E(numCiphertext);
    std::vector<bool> Einitialized(numCiphertext, false);

    std::vector<double> triangularMask(subVectorLength * subVectorLength, 0.0);
    for (size_t i = 0; i < subVectorLength; i++)
        for (size_t j = 0; j < subVectorLength; j++)
            if (j <= i) triangularMask[i * subVectorLength + j] = 1.0;

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

        Ciphertext<DCRTPoly> Cjk = compare(
            replR[j],
            replC[k],
            leftBoundC, rightBoundC, degreeC
        );

        Ciphertext<DCRTPoly> Ejk = 4 * (1 - Cjk) * Cjk;

        #pragma omp critical
        {
        if (!Cvinitialized[j]) { Cv[j] = Cjk; Cvinitialized[j] = true; }
        else                   { Cv[j] = Cv[j] + Cjk;                  }
        }

        #pragma omp critical
        {
        if (!Evinitialized[j]) { Ev[j] = Ejk; Evinitialized[j] = true; }
        else                   { Ev[j] = Ev[j] + Ejk;                  }
        }

        if (j == k)
        {
            #pragma omp critical
            {
            if (!Einitialized[j]) { E[j] = Ejk * triangularMask; Einitialized[j] = true; }
            else                  { E[j] = E[j] + Ejk * triangularMask;                  }
            }
        }
        else
        {
            Ciphertext<DCRTPoly> Ckj = 1.0 - Cjk;

            #pragma omp critical
            {
            if (!Chinitialized[k]) { Ch[k] = Ckj; Chinitialized[k] = true; }
            else                   { Ch[k] = Ch[k] + Ckj;                  }
            }

            #pragma omp critical
            {
            if (!Ehinitialized[k]) { Eh[k] = Ejk; Ehinitialized[k] = true; }
            else                   { Eh[k] = Eh[k] + Ejk;                  }
            }

            #pragma omp critical
            {
            if (!Einitialized[k]) { E[k] = Ejk; Einitialized[k] = true; }
            else                  { E[k] = E[k] + Ejk;                  }
            }
        }
        
        #pragma omp critical
        {std::cout << "(j, k) = (" << j << ", " << k << ") - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "COMPLETED (" <<
            elapsed_seconds.count() << "s)" << std::endl << std::endl;

    std::cout << "===================================\n";
    std::cout << "Sum\n";
    std::cout << "===================================\n";

    std::vector<Ciphertext<DCRTPoly>> s(numCiphertext);
    std::vector<bool> sinitialized(numCiphertext, false);

    start = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for collapse(2)
    for (size_t loopID = 0; loopID < 5; loopID++)
    {
        for (size_t j = 0; j < numCiphertext; j++)
        {
            if (loopID == 0)
            {
                #pragma omp critical
                {std::cout << "SumV - j = " << j << " - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

                Ciphertext<DCRTPoly> svj = sumRows(Cv[j], subVectorLength);

                #pragma omp critical
                {
                if (!sinitialized[j]) { s[j] = svj; sinitialized[j] = true; }
                else                  { s[j] = s[j] + svj;                  }
                }

                #pragma omp critical
                {std::cout << "SumV - j = " << j << " - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
            }
            else if (loopID == 1)
            {
                #pragma omp critical
                {std::cout << "SumH - j = " << j << " - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

                if (j > 0)
                {
                    Ciphertext<DCRTPoly> shj = sumColumns(Ch[j], subVectorLength, true);
                    shj = transposeColumn(shj, subVectorLength, true);
                    shj = replicateRow(shj, subVectorLength);

                    #pragma omp critical
                    {
                    if (!sinitialized[j]) { s[j] = shj; sinitialized[j] = true; }
                    else                  { s[j] = s[j] + shj;                  }
                    }
                }

                #pragma omp critical
                {std::cout << "SumH - j = " << j << " - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
            }
            else if (loopID == 2)
            {
                #pragma omp critical
                {std::cout << "SumV - j = " << j << " - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

                Ciphertext<DCRTPoly> evj = sumRows(Ev[j], subVectorLength);
                
                #pragma omp critical
                {
                if (!sinitialized[j]) { s[j] = - 0.5 * evj; sinitialized[j] = true; }
                else                  { s[j] = s[j] - 0.5 * evj;                    }
                }

                #pragma omp critical
                {std::cout << "SumV - j = " << j << " - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
            }
            else if (loopID == 3)
            {
                #pragma omp critical
                {std::cout << "SumH - j = " << j << " - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

                if (j > 0)
                {
                    Ciphertext<DCRTPoly> ehj = sumColumns(Eh[j], subVectorLength, true);
                    ehj = transposeColumn(ehj, subVectorLength, true);
                    ehj = replicateRow(ehj, subVectorLength);

                    #pragma omp critical
                    {
                    if (!sinitialized[j]) { s[j] = - 0.5 * ehj; sinitialized[j] = true; }
                    else                  { s[j] = s[j] - 0.5 * ehj;                    }
                    }
                }

                #pragma omp critical
                {std::cout << "SumH - j = " << j << " - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
            }
            else if (loopID == 4)
            {
                #pragma omp critical
                {std::cout << "SumH - j = " << j << " - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}
                
                Ciphertext<DCRTPoly> ej = sumColumns(E[j], subVectorLength, true);
                ej = transposeColumn(ej, subVectorLength, true);
                ej = replicateRow(ej, subVectorLength);
                
                #pragma omp critical
                {
                if (!sinitialized[j]) { s[j] = ej; sinitialized[j] = true; }
                else                  { s[j] = s[j] + ej;                  }
                }

                #pragma omp critical
                {std::cout << "SumH - j = " << j << " - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
            }
        }
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "COMPLETED (" <<
            elapsed_seconds.count() << "s)" << std::endl << std::endl;
    
    std::cout << "===================================\n";
    std::cout << "Sort\n";
    std::cout << "===================================\n";

    start = std::chrono::high_resolution_clock::now();

    std::vector<std::vector<double>> subMasks(numCiphertext);
    for (size_t i = 0; i < numCiphertext; i++)
    {
        std::vector<double> subMask(subVectorLength * subVectorLength);
        for (size_t j = 0; j < subVectorLength; j++)
            for (size_t k = 0; k < subVectorLength; k++)
                subMask[j * subVectorLength + k] = -1.0 * (i * subVectorLength + j) - 1.0;
        subMasks[i] = subMask;
    }

    std::vector<Ciphertext<DCRTPoly>> subSorted(numCiphertext);
    std::vector<bool> subSortedInitialized(numCiphertext, false);

    #pragma omp parallel for collapse(2)
    for (size_t j = 0; j < numCiphertext; j++)
    {
        for (size_t k = 0; k < numCiphertext; k++)
        {
            #pragma omp critical
            {std::cout << "(j, k) = (" << j << ", " << k << ") - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

            Ciphertext<DCRTPoly> ind = indicator(
                s[k] + subMasks[j],
                -0.5, 0.5,
                -1.01 * vectorLength, 1.01 * vectorLength,
                degreeI
            ) * replR[k];

            #pragma omp critical
            {
            if (!subSortedInitialized[j]) { subSorted[j] = ind; subSortedInitialized[j] = true; }
            else                          { subSorted[j] = subSorted[j] + ind;                  }
            }
            
            #pragma omp critical
            {std::cout << "(j, k) = (" << j << ", " << k << ") - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
        }
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "COMPLETED (" <<
            elapsed_seconds.count() << "s)" << std::endl << std::endl;
    
    std::cout << "===================================\n";
    std::cout << "Sum\n";
    std::cout << "===================================\n";

    start = std::chrono::high_resolution_clock::now();

    std::vector<Ciphertext<DCRTPoly>> result(numCiphertext);

    #pragma omp parallel for
    for (size_t j = 0; j < numCiphertext; j++)
    {
        #pragma omp critical
        {std::cout << "Sum - j = " << j << " - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

        result[j] = sumColumns(subSorted[j], subVectorLength);

        #pragma omp critical
        {std::cout << "Sum - j = " << j << " - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "COMPLETED (" <<
            elapsed_seconds.count() << "s)" << std::endl << std::endl;

    return result;

}


Ciphertext<DCRTPoly> sortFG(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    uint32_t dg_c,
    uint32_t df_c,
    uint32_t dg_i,
    uint32_t df_i
)
{
    Ciphertext<DCRTPoly> VR = replicateRow(c, vectorLength);
    Ciphertext<DCRTPoly> VC = replicateColumn(transposeRow(c, vectorLength, true), vectorLength);
    Ciphertext<DCRTPoly> C = compareAdv(
        VR, VC,
        dg_c, df_c
    );
    std::cout << "C levels: " << C->GetLevel() << std::endl;
    Ciphertext<DCRTPoly> R = sumRows(C, vectorLength);

    std::vector<double> subMask(vectorLength * vectorLength);
    for (size_t i = 0; i < vectorLength; i++)
        for (size_t j = 0; j < vectorLength; j++)
            subMask[i * vectorLength + j] = -1.0 * i - 0.5;
    Ciphertext<DCRTPoly> M = indicatorAdv(
        R + subMask,
        vectorLength,
        dg_i, df_i
    );
    std::cout << "M levels: " << M->GetLevel() << std::endl;

    Ciphertext<DCRTPoly> S = sumColumns(M * VR, vectorLength);
    std::cout << "S levels: " << S->GetLevel() << std::endl;

    return S;
}


Ciphertext<DCRTPoly> sortWithCorrectionFG(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    uint32_t dg_c,
    uint32_t df_c,
    uint32_t dg_i,
    uint32_t df_i
)
{
    Ciphertext<DCRTPoly> VR = replicateRow(c, vectorLength);
    Ciphertext<DCRTPoly> VC = replicateColumn(transposeRow(c, vectorLength, true), vectorLength);
    Ciphertext<DCRTPoly> C = compareAdv(
        VR, VC,
        dg_c, df_c
    );
    std::cout << "C levels: " << C->GetLevel() << std::endl;

    std::vector<double> triangularMask(vectorLength * vectorLength, 0.0);
    for (size_t i = 0; i < vectorLength; i++)
        for (size_t j = 0; j < vectorLength; j++)
            if (j >= i) triangularMask[i * vectorLength + j] = 1.0;
    Ciphertext<DCRTPoly> E = 4 * (1 - C) * C;
    Ciphertext<DCRTPoly> correctionOffset = sumRows(E * triangularMask, vectorLength) - 0.5 * sumRows(E, vectorLength);

    Ciphertext<DCRTPoly> R = sumRows(C, vectorLength) + correctionOffset;

    std::vector<double> subMask(vectorLength * vectorLength);
    for (size_t i = 0; i < vectorLength; i++)
        for (size_t j = 0; j < vectorLength; j++)
            subMask[i * vectorLength + j] = -1.0 * i - 1.0;
    Ciphertext<DCRTPoly> M = indicatorAdv(
        R + subMask,
        vectorLength,
        dg_i, df_i
    );
    std::cout << "M levels: " << M->GetLevel() << std::endl;

    Ciphertext<DCRTPoly> S = sumColumns(M * VR, vectorLength);
    std::cout << "S levels: " << S->GetLevel() << std::endl;

    return S;
}


std::vector<Ciphertext<DCRTPoly>> sortWithCorrectionFG(
    const std::vector<Ciphertext<DCRTPoly>> &c,
    const size_t subVectorLength,
    uint32_t dg_c,
    uint32_t df_c,
    uint32_t dg_i,
    uint32_t df_i
)
{
    const size_t numCiphertext = c.size();
    const size_t vectorLength = subVectorLength * numCiphertext;

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
    std::cout << "replR levels: " << replR[0]->GetLevel() << std::endl;
    std::cout << "replC levels: " << replC[0]->GetLevel() << std::endl;

    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "COMPLETED (" <<
            elapsed_seconds.count() << "s)" << std::endl << std::endl;
    
    std::cout << "===================================\n";
    std::cout << "Compare\n";
    std::cout << "===================================\n";

    std::vector<Ciphertext<DCRTPoly>> Cv(numCiphertext);
    std::vector<Ciphertext<DCRTPoly>> Ch(numCiphertext);
    std::vector<bool> Cvinitialized(numCiphertext, false);
    std::vector<bool> Chinitialized(numCiphertext, false);

    std::vector<Ciphertext<DCRTPoly>> Ev(numCiphertext);
    std::vector<Ciphertext<DCRTPoly>> Eh(numCiphertext);
    std::vector<bool> Evinitialized(numCiphertext, false);
    std::vector<bool> Ehinitialized(numCiphertext, false);

    std::vector<Ciphertext<DCRTPoly>> E(numCiphertext);
    std::vector<bool> Einitialized(numCiphertext, false);

    std::vector<double> triangularMask(subVectorLength * subVectorLength, 0.0);
    for (size_t i = 0; i < subVectorLength; i++)
        for (size_t j = 0; j < subVectorLength; j++)
            if (j <= i) triangularMask[i * subVectorLength + j] = 1.0;

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

        Ciphertext<DCRTPoly> Cjk = compareAdv(
            replR[j],
            replC[k],
            dg_c, df_c
        );

        Ciphertext<DCRTPoly> Ejk = 4 * (1 - Cjk) * Cjk;

        #pragma omp critical
        {
        if (!Cvinitialized[j]) { Cv[j] = Cjk; Cvinitialized[j] = true; }
        else                   { Cv[j] = Cv[j] + Cjk;                  }
        }

        #pragma omp critical
        {
        if (!Evinitialized[j]) { Ev[j] = Ejk; Evinitialized[j] = true; }
        else                   { Ev[j] = Ev[j] + Ejk;                  }
        }

        if (j == k)
        {
            #pragma omp critical
            {
            if (!Einitialized[j]) { E[j] = Ejk * triangularMask; Einitialized[j] = true; }
            else                  { E[j] = E[j] + Ejk * triangularMask;                  }
            }
        }
        else
        {
            Ciphertext<DCRTPoly> Ckj = 1.0 - Cjk;

            #pragma omp critical
            {
            if (!Chinitialized[k]) { Ch[k] = Ckj; Chinitialized[k] = true; }
            else                   { Ch[k] = Ch[k] + Ckj;                  }
            }

            #pragma omp critical
            {
            if (!Ehinitialized[k]) { Eh[k] = Ejk; Ehinitialized[k] = true; }
            else                   { Eh[k] = Eh[k] + Ejk;                  }
            }

            #pragma omp critical
            {
            if (!Einitialized[k]) { E[k] = Ejk; Einitialized[k] = true; }
            else                  { E[k] = E[k] + Ejk;                  }
            }
        }
        
        #pragma omp critical
        {std::cout << "(j, k) = (" << j << ", " << k << ") - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
    }
    std::cout << "Cv levels: " << Cv[0]->GetLevel() << std::endl;
    std::cout << "Ch levels: " << Ch[1]->GetLevel() << std::endl;
    std::cout << "Ev levels: " << Ev[0]->GetLevel() << std::endl;
    std::cout << "Eh levels: " << Eh[1]->GetLevel() << std::endl;
    std::cout << "E levels: " << E[0]->GetLevel() << std::endl;

    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "COMPLETED (" <<
            elapsed_seconds.count() << "s)" << std::endl << std::endl;

    std::cout << "===================================\n";
    std::cout << "Sum\n";
    std::cout << "===================================\n";

    std::vector<Ciphertext<DCRTPoly>> s(numCiphertext);
    std::vector<bool> sinitialized(numCiphertext, false);

    start = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for collapse(2)
    for (size_t loopID = 0; loopID < 5; loopID++)
    {
        for (size_t j = 0; j < numCiphertext; j++)
        {
            if (loopID == 0)
            {
                #pragma omp critical
                {std::cout << "SumV - j = " << j << " - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

                Ciphertext<DCRTPoly> svj = sumRows(Cv[j], subVectorLength);

                #pragma omp critical
                {
                if (!sinitialized[j]) { s[j] = svj; sinitialized[j] = true; }
                else                  { s[j] = s[j] + svj;                  }
                }

                #pragma omp critical
                {std::cout << "SumV - j = " << j << " - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
            }
            else if (loopID == 1)
            {
                #pragma omp critical
                {std::cout << "SumH - j = " << j << " - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

                if (j > 0)
                {
                    Ciphertext<DCRTPoly> shj = sumColumns(Ch[j], subVectorLength, true);
                    shj = transposeColumn(shj, subVectorLength, true);
                    shj = replicateRow(shj, subVectorLength);

                    #pragma omp critical
                    {
                    if (!sinitialized[j]) { s[j] = shj; sinitialized[j] = true; }
                    else                  { s[j] = s[j] + shj;                  }
                    }
                }

                #pragma omp critical
                {std::cout << "SumH - j = " << j << " - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
            }
            else if (loopID == 2)
            {
                #pragma omp critical
                {std::cout << "SumV - j = " << j << " - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

                Ciphertext<DCRTPoly> evj = sumRows(Ev[j], subVectorLength);
                
                #pragma omp critical
                {
                if (!sinitialized[j]) { s[j] = - 0.5 * evj; sinitialized[j] = true; }
                else                  { s[j] = s[j] - 0.5 * evj;                    }
                }

                #pragma omp critical
                {std::cout << "SumV - j = " << j << " - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
            }
            else if (loopID == 3)
            {
                #pragma omp critical
                {std::cout << "SumH - j = " << j << " - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

                if (j > 0)
                {
                    Ciphertext<DCRTPoly> ehj = sumColumns(Eh[j], subVectorLength, true);
                    ehj = transposeColumn(ehj, subVectorLength, true);
                    ehj = replicateRow(ehj, subVectorLength);

                    #pragma omp critical
                    {
                    if (!sinitialized[j]) { s[j] = - 0.5 * ehj; sinitialized[j] = true; }
                    else                  { s[j] = s[j] - 0.5 * ehj;                    }
                    }
                }

                #pragma omp critical
                {std::cout << "SumH - j = " << j << " - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
            }
            else if (loopID == 4)
            {
                #pragma omp critical
                {std::cout << "SumH - j = " << j << " - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}
                
                Ciphertext<DCRTPoly> ej = sumColumns(E[j], subVectorLength, true);
                ej = transposeColumn(ej, subVectorLength, true);
                ej = replicateRow(ej, subVectorLength);
                
                #pragma omp critical
                {
                if (!sinitialized[j]) { s[j] = ej; sinitialized[j] = true; }
                else                  { s[j] = s[j] + ej;                  }
                }

                #pragma omp critical
                {std::cout << "SumH - j = " << j << " - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
            }
        }
    }
    std::cout << "S levels: " << s[0]->GetLevel() << std::endl;

    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "COMPLETED (" <<
            elapsed_seconds.count() << "s)" << std::endl << std::endl;
    
    std::cout << "===================================\n";
    std::cout << "Sort\n";
    std::cout << "===================================\n";

    start = std::chrono::high_resolution_clock::now();

    std::vector<std::vector<double>> subMasks(numCiphertext);
    for (size_t i = 0; i < numCiphertext; i++)
    {
        std::vector<double> subMask(subVectorLength * subVectorLength);
        for (size_t j = 0; j < subVectorLength; j++)
            for (size_t k = 0; k < subVectorLength; k++)
                subMask[j * subVectorLength + k] = -1.0 * (i * subVectorLength + j) - 1.0;
        subMasks[i] = subMask;
    }

    std::vector<Ciphertext<DCRTPoly>> subSorted(numCiphertext);
    std::vector<bool> subSortedInitialized(numCiphertext, false);

    #pragma omp parallel for collapse(2)
    for (size_t j = 0; j < numCiphertext; j++)
    {
        for (size_t k = 0; k < numCiphertext; k++)
        {
            #pragma omp critical
            {std::cout << "(j, k) = (" << j << ", " << k << ") - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

            Ciphertext<DCRTPoly> ind = indicatorAdv(
                s[k] + subMasks[j],
                vectorLength,
                dg_i, df_i
            ) * replR[k];

            #pragma omp critical
            {
            if (!subSortedInitialized[j]) { subSorted[j] = ind; subSortedInitialized[j] = true; }
            else                          { subSorted[j] = subSorted[j] + ind;                  }
            }
            
            #pragma omp critical
            {std::cout << "(j, k) = (" << j << ", " << k << ") - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
        }
    }
    std::cout << "subSorted levels: " << subSorted[0]->GetLevel() << std::endl;

    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "COMPLETED (" <<
            elapsed_seconds.count() << "s)" << std::endl << std::endl;
    
    std::cout << "===================================\n";
    std::cout << "Sum\n";
    std::cout << "===================================\n";

    start = std::chrono::high_resolution_clock::now();

    std::vector<Ciphertext<DCRTPoly>> result(numCiphertext);

    #pragma omp parallel for
    for (size_t j = 0; j < numCiphertext; j++)
    {
        #pragma omp critical
        {std::cout << "Sum - j = " << j << " - thread_id = " << omp_get_thread_num() << " - START" << std::endl;}

        result[j] = sumColumns(subSorted[j], subVectorLength);

        #pragma omp critical
        {std::cout << "Sum - j = " << j << " - thread_id = " << omp_get_thread_num() << " - END" << std::endl;}
    }

    std::cout << "result levels: " << result[0]->GetLevel() << std::endl;

    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << "COMPLETED (" <<
            elapsed_seconds.count() << "s)" << std::endl << std::endl;

    return result;

}
