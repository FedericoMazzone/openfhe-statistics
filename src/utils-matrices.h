#pragma once

#include "utils-basics.h"


using namespace lbcrypto;


// Assuming square matrices of size a power of 2.


/**
 * @brief Masks out everything except a specific row of a ciphertext matrix.
 * 
 * This function masks out everything except a specific row of a ciphertext
 * matrix, multiplying by a proper plaintext mask.
 * 
 * @param c The original ciphertext matrix to be masked.
 * @param matrixSize The size of the square matrix.
 * @param rowIndex The index of the row to be preserved (0-based index).
 * @return Ciphertext<DCRTPoly> The masked ciphertext matrix.
 * 
 * @throws std::out_of_range if the rowIndex is out of range (not within 0 to
 * matrixSize - 1).
 */
Ciphertext<DCRTPoly> maskRow(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    const size_t rowIndex
);


/**
 * @brief Masks out everything except a specific column of a ciphertext matrix.
 * 
 * This function masks out everything except a specific column of a ciphertext
 * matrix, multiplying by a proper plaintext mask.
 * 
 * @param c The original ciphertext matrix to be masked.
 * @param matrixSize The size of the square matrix.
 * @param columnIndex The index of the column to be preserved (0-based index).
 * @return Ciphertext<DCRTPoly> The masked ciphertext matrix.
 * 
 * @throws std::out_of_range if the columnIndex is out of range (not within 0
 * to matrixSize - 1).
 */
Ciphertext<DCRTPoly> maskColumn(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    const size_t columnIndex
);


/**
 * @brief Replicates a single row of a ciphertext matrix to all other rows.
 * 
 * This function takes a ciphertext matrix with only one row non-zero and
 * replicates this row to all other rows of the matrix. This function works
 * recursively, thus the number of rows in the matrix is expected to be a power
 * of 2.
 * 
 * @param c The original ciphertext matrix with a single non-zero row.
 * @param matrixSize The size of the square matrix (number of rows/columns).
 * @return Ciphertext<DCRTPoly> The ciphertext matrix with the single row
 * replicated to all other rows.
 */
Ciphertext<DCRTPoly> replicateRow(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize
);


/**
 * @brief Replicates the first column of a ciphertext matrix to all other
 * columns.
 * 
 * This function takes a ciphertext matrix with only the first column non-zero
 * and replicates this column to all other columns of the matrix. This function
 * works recursively, thus the number of columns in the matrix is expected to
 * be a power of 2.
 * 
 * @param c The original ciphertext matrix with the first column non-zero.
 * @param matrixSize The size of the square matrix (number of rows/columns).
 * @return Ciphertext<DCRTPoly> The ciphertext matrix with the first column
 * replicated to all other columns.
 */
Ciphertext<DCRTPoly> replicateColumn(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize
);


/**
 * @brief Sums all rows of a ciphertext matrix into a single row.
 * 
 * This function sums all rows of a ciphertext matrix. The result is a
 * ciphertext matrix where each row is the component-wise sum of all the rows
 * of the original matrix. This function works recursively, thus the number of
 * columns in the matrix is expected to be a power of 2.
 * 
 * @param c The original ciphertext matrix.
 * @param matrixSize The size of the square matrix (number of rows/columns).
 * @param maskOutput If true, masks out everything but the output row (default
 * is false).
 * @param outputRow The index of the row to store in the sum output (default is
 * 0).
 * @return Ciphertext<DCRTPoly> The ciphertext matrix with all rows summed into
 * a single row.
 * 
 * @throws std::out_of_range if the outputRow index is out of range (not within
 * 0 to matrixSize - 1).
 */
Ciphertext<DCRTPoly> sumRows(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    bool maskOutput=false,
    const size_t outputRow=0
);


/**
 * @brief Sums all columns of a ciphertext matrix into a single column.
 * 
 * This function sums all columns of a ciphertext matrix. The result is a
 * ciphertext matrix where each column is the component-wise sum of all the
 * columns of the original matrix. This function works recursively, thus the
 * number of rows in the matrix is expected to be a power of 2.
 * 
 * @param c The original ciphertext matrix.
 * @param matrixSize The size of the square matrix (number of rows/columns).
 * @param maskOutput If true, masks out everything but the output column
 * (default is false).
 * @return Ciphertext<DCRTPoly> The ciphertext matrix with all columns summed
 * into the first column.
 */
Ciphertext<DCRTPoly> sumColumns(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    bool maskOutput=false
);


/**
 * @brief Transposes a ciphertext row vector into a column vector.
 * 
 * This function takes a ciphertext matrix with only the first row non-zero and
 * transposes it. It assumes the input matrix to be a square matrix of size a
 * power of 2.
 * 
 * @param c The original ciphertext matrix with only the first row non-zero.
 * @param matrixSize The size of the square matrix (number of rows/columns).
 * @param maskOutput If true, masks out everything but the output column
 * (default is false).
 * @return Ciphertext<DCRTPoly> The transposed ciphertext matrix.
 */
Ciphertext<DCRTPoly> transposeRow(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    bool maskOutput=false
);


/**
 * @brief Transposes a ciphertext column vector into a row vector.
 * 
 * This function takes a ciphertext matrix with only the first column non-zero
 * and transposes it. It assumes the input matrix to be a square matrix of size
 * a power of 2.
 * 
 * @param c The original ciphertext matrix with only the first column non-zero.
 * @param matrixSize The size of the square matrix (number of rows/columns).
 * @param maskOutput If true, masks out everything but the output row (default
 * is false).
 * @return Ciphertext<DCRTPoly> The transposed ciphertext matrix.
 */
Ciphertext<DCRTPoly> transposeColumn(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    bool maskOutput=false
);


/**
 * @brief Generates rotation indices for homomorphic operations on a ciphertext
 * matrix.
 * 
 * This function generates rotation indices for homomorphic operations on a
 * ciphertext matrix. The indices are used to specify the number of shifts
 * required for various operations such as summing rows/columns, replicating
 * rows/columns, and transposing rows/columns. The rotation indices are
 * computed based on the size of the square matrix.
 * 
 * @param matrixSize The size of the square matrix (number of rows/columns).
 * @return std::vector<int32_t> A vector containing the rotation indices.
 */
std::vector<int32_t> getRotationIndices(
    const size_t matrixSize
);
