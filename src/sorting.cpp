#include "sorting.h"
#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-matrices.h"
#include "utils-ptxt.h"
#include <cassert>
#include <omp.h>


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
    uint32_t degreeI,
    bool advanced
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

        Ciphertext<DCRTPoly> Cjk;
        if (advanced)
            Cjk = compareAdv(
                replR[j],
                replC[k],
                degreeC, degreeI
            );
        else
            Cjk = compare(
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

            Ciphertext<DCRTPoly> ind;
            if (advanced)
                ind = indicatorAdv(
                    s[k] + subMasks[j],
                    1.01 * vectorLength,
                    degreeC, degreeI
                ) * replR[k];
            else
                ind = indicator(
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


std::vector<double> sort(
    std::vector<double> vec
)
{
    std::sort(vec.begin(), vec.end());
    return vec;
}


std::vector<double> testSorting(
    const size_t vectorLength = 8,
    const usint compareDepth = 7,
    const usint indicatorDepth = 11,
    const bool tieCorrection = false
)
{

    std::cout << "Vector length: " << vectorLength << std::endl;
    std::cout << "Compare depth: " << compareDepth << std::endl;
    std::cout << "Indicator depth: " << indicatorDepth << std::endl;
    std::cout << "Tie correction: " << tieCorrection << std::endl;

    const usint integralPrecision       = 12;
    const usint decimalPrecision        = 48;
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
    std::cout << elapsed_seconds.count() << "s" << std::endl;

    Plaintext resultP;
    cryptoContext->Decrypt(keyPair.secretKey, resultC, &resultP);
    resultP->SetLength(vectorLength * vectorLength);
    std::vector<double> resultMatrix = resultP->GetRealPackedValue();
    std::vector<double> result(vectorLength);
    for (size_t i = 0; i < vectorLength; i++)
        result[i] = resultMatrix[i * vectorLength];
    std::cout << "Sorting:          " << result << std::endl;
    // Plaintext resultP;
    // cryptoContext->Decrypt(keyPair.secretKey, resultC, &resultP);
    // resultP->SetLength(vectorLength * vectorLength);
    // std::vector<double> resultMatrix = resultP->GetRealPackedValue();
    // std::cout << "Sorting:          " << vector2matrix(resultMatrix, vectorLength) << std::endl;

    return {
        (double) vectorLength, (double) compareDepth, (double) indicatorDepth,
        elapsed_seconds.count()
    };

}


std::vector<double> testSortingMultiCtxt(
    const size_t subVectorLength = 128,
    const size_t numCiphertext = 2,
    const usint compareDepth = 7,
    const usint indicatorDepth = 11,
    const bool tieCorrection = false
)
{

    std::cout << "SubVector length: " << subVectorLength << std::endl;
    std::cout << "Number of ciphertexts: " << numCiphertext << std::endl;
    std::cout << "Compare depth: " << compareDepth << std::endl;
    std::cout << "Indicator depth: " << indicatorDepth << std::endl;
    std::cout << "Tie correction: " << tieCorrection << std::endl;

    const size_t vectorLength           = subVectorLength * numCiphertext;
    const usint integralPrecision       = 12;
    const usint decimalPrecision        = 48;
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
    std::cout << elapsed_seconds.count() << "s" << std::endl;

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
    // Plaintext resultP;
    // cryptoContext->Decrypt(keyPair.secretKey, resultC[0], &resultP);
    // resultP->SetLength(subVectorLength * subVectorLength);
    // std::vector<double> resultMatrix = resultP->GetRealPackedValue();
    // std::cout << "Sorting:          " << vector2matrix(resultMatrix, subVectorLength) << std::endl;

    return {
        (double) vectorLength, (double) compareDepth, (double) indicatorDepth,
        elapsed_seconds.count()
    };

}


void testSortingMultiCtxtHighPrecision(
    const size_t subVectorLength = 128,
    const size_t numCiphertext = 2,
    const usint dg = 4,
    const usint df = 2
)
{

    std::cout << "SubVector length: " << subVectorLength << std::endl;
    std::cout << "Number of ciphertexts: " << numCiphertext << std::endl;
    std::cout << "dg: " << dg << std::endl;
    std::cout << "df: " << df << std::endl;

    const size_t vectorLength           = subVectorLength * numCiphertext;
    const usint integralPrecision       = 1;
    const usint decimalPrecision        = 59;
    const usint multiplicativeDepth     = (dg + df) * 10 + 6;
    const usint numSlots                = subVectorLength * subVectorLength;
    const bool enableBootstrap          = false;
    const usint ringDim                 = 1 << 15;
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
    resultC = sortWithCorrection(
        vC,
        subVectorLength,
        -1.0, 1.0,
        dg, df,
        true
    );

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << elapsed_seconds.count() << "s" << std::endl;

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
    // Plaintext resultP;
    // cryptoContext->Decrypt(keyPair.secretKey, resultC[0], &resultP);
    // resultP->SetLength(subVectorLength * subVectorLength);
    // std::vector<double> resultMatrix = resultP->GetRealPackedValue();
    // std::cout << "Sorting:          " << vector2matrix(resultMatrix, subVectorLength) << std::endl;
}


// int main()
// {

//     srand(time(NULL));
//     std::cout << std::fixed << std::setprecision(2);

//     // High precision cmp function
//     const size_t numCiphertext = 2;
//     const size_t dg = 5;
//     const size_t df = 2;
//     testSortingMultiCtxtHighPrecision(128, numCiphertext, dg, df);

//     const size_t TEST_REPS = 10;

//     std::ofstream logFile;
//     logFile.open("log-sorting2.txt");

//     logFile << "vector length vs. run time" << std::endl;
//     usint compareDepth = 9;
//     usint indicatorDepth = 10;
//     bool tieCorrection = true;
//     logFile << "compareDepth=" << compareDepth << std::endl;
//     logFile << "indicatorDepth=" << indicatorDepth << std::endl;
//     logFile << "tieCorrection=" << tieCorrection << std::endl;
//     for (size_t vectorLength = 8; vectorLength <= 128; vectorLength *= 2)
//     {
//         std::vector<std::vector<double>> results;
//         for (size_t i = 0; i < TEST_REPS; i++)
//             try { results.push_back(testSorting(vectorLength, compareDepth, indicatorDepth, tieCorrection)); }
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
//             try { results.push_back(testSortingMultiCtxt(128, numCiphertext, compareDepth, indicatorDepth, tieCorrection)); }
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

//     logFile.close();

//     return 0;

// }
