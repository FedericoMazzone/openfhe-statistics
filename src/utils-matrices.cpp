#include "utils-matrices.h"
#include "utils-ptxt.h"
#include <cassert>



Ciphertext<DCRTPoly> maskRow(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    const size_t rowIndex
)
{
    assert(rowIndex >= 0 && rowIndex < matrixSize && "Invalid row index");

    std::vector<double> mask(matrixSize * matrixSize, 0.0);
    for (size_t i = 0; i < matrixSize; i++)
        mask[matrixSize * rowIndex + i] = 1.0;
    
    return c * mask;
}


Ciphertext<DCRTPoly> maskColumn(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    const size_t columnIndex
)
{
    assert(columnIndex >= 0 && columnIndex < matrixSize && "Invalid column index");

    std::vector<double> mask(matrixSize * matrixSize, 0.0);
    for (size_t i = 0; i < matrixSize; i++)
        mask[matrixSize * i + columnIndex] = 1.0;
    
    return c * mask;
}


Ciphertext<DCRTPoly> replicateRow(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize
)
{
    for (size_t i = 0; i < LOG2(matrixSize); i++)
        c += c >> (1 << (LOG2(matrixSize) + i));

    return c;
}


Ciphertext<DCRTPoly> replicateColumn(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize
)
{
    for (size_t i = 0; i < LOG2(matrixSize); i++)
        c += c >> (1 << i);

    return c;
}


Ciphertext<DCRTPoly> sumRows(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    bool maskOutput,
    const size_t outputRow
)
{
    for (size_t i = 0; i < LOG2(matrixSize); i++)
        c += c >> (1 << (LOG2(matrixSize) + i));

    if (maskOutput)
        c = maskRow(c, matrixSize, outputRow);

    return c;
}


Ciphertext<DCRTPoly> sumColumns(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    bool maskOutput
)
{
    for (size_t i = 0; i < LOG2(matrixSize); i++)
        c += c << (1 << i);
    
    if (maskOutput)
        c = maskColumn(c, matrixSize, 0);
    
    return c;
}


Ciphertext<DCRTPoly> transposeRow(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    bool maskOutput
)
{
    for (size_t i = 1; i <= LOG2(matrixSize); i++)
        c += c >> (matrixSize * (matrixSize - 1) / (1 << i));
    
    if (maskOutput)
        c = maskColumn(c, matrixSize, 0);

    return c;
}


Ciphertext<DCRTPoly> transposeColumn(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    bool maskOutput
)
{
    for (size_t i = 1; i <= LOG2(matrixSize); i++)
        c += c << (matrixSize * (matrixSize - 1) / (1 << i));

    if (maskOutput)
        c = maskRow(c, matrixSize, 0);

    return c;
}


std::vector<int32_t> getRotationIndices
(
    const size_t matrixSize
)
{
    std::vector<int32_t> indices;

    int32_t index;
    for (size_t i = 0; i < LOG2(matrixSize); i++)
    {
        index = 1 << i;
        indices.push_back(index);   // sumColumns
        indices.push_back(-index);  // replicateColumn

        index = 1 << (LOG2(matrixSize) + i);
        indices.push_back(-index);  // replicateRow, sumRows

        index = (matrixSize * (matrixSize - 1) / (1 << (i + 1)));
        indices.push_back(index);   // transposeColumn
        indices.push_back(-index);  // transposeRow
    }

    return indices;
}


// int main()
// {

//     const usint matrixSize              = 8;

//     const usint integralPrecision       = 10;
//     const usint decimalPrecision        = 50;
//     const usint multiplicativeDepth     = 2;
//     const usint numSlots                = matrixSize * matrixSize;
//     const bool enableBootstrap          = false;
//     const usint ringDim                 = 0;
//     const bool verbose                  = true;

//     std::vector<int32_t> indices = getRotationIndices(matrixSize);

//     CryptoContext<DCRTPoly> cryptoContext = generateCryptoContext(
//         integralPrecision,
//         decimalPrecision,
//         multiplicativeDepth,
//         numSlots,
//         enableBootstrap,
//         ringDim,
//         verbose
//     );

//     KeyPair<DCRTPoly> keyPair = keyGeneration(
//         cryptoContext,
//         indices,
//         numSlots,
//         enableBootstrap,
//         verbose
//     );

//     std::vector<double> matrix(matrixSize * matrixSize);
//     for (size_t i = 0; i < matrixSize; i++)
//         for (size_t j = 0; j < matrixSize; j++)
//             matrix[i * matrixSize + j] = i * matrixSize + j;

//     std::cout << "Matrix: " << vector2matrix(matrix, matrixSize) << std::endl;

//     Ciphertext<DCRTPoly> matrixC = cryptoContext->Encrypt(
//         keyPair.publicKey,
//         cryptoContext->MakeCKKSPackedPlaintext(matrix)
//     );

//     auto start = std::chrono::high_resolution_clock::now();

//     Ciphertext<DCRTPoly> resultC = matrixC;
//     // Ciphertext<DCRTPoly> resultC = maskRow(matrixC, matrixSize, 0);
//     // resultC = replicateRow(resultC, matrixSize);
//     // resultC = transposeRow(resultC, matrixSize, true);
//     // Ciphertext<DCRTPoly> resultC = maskColumn(matrixC, matrixSize, 0);
//     // resultC = replicateColumn(resultC, matrixSize);
//     // resultC = transposeColumn(resultC, matrixSize, true);
//     // Ciphertext<DCRTPoly> resultC = sumRows(matrixC, matrixSize, index);
//     // Ciphertext<DCRTPoly> resultC = sumColumns(matrixC, matrixSize);

//     auto end = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> elapsed_seconds = end - start;
//     std::cout << elapsed_seconds.count() << "s" << std::endl;

//     Plaintext resultP;
//     cryptoContext->Decrypt(keyPair.secretKey, resultC, &resultP);
//     resultP->SetLength(matrixSize * matrixSize);

//     std::vector<std::vector<double>> result = vector2matrix(
//         resultP->GetRealPackedValue(),
//         matrixSize
//     );
//     std::cout << "Result: " << result << std::endl;

//     return 0;

// }
