#pragma once

#include <iomanip>
#include <iostream>
#include <vector>

#define LOG2(X) (size_t) std::ceil(std::log2((X)))
#define MIN_VEC(V) *std::min_element(V.begin(), V.end())
#define MAX_VEC(V) *std::max_element(V.begin(), V.end())


/**
 * Calculates the average of vectors component-wise.
 *
 * @param vectors A vector of vectors of doubles to be averaged.
 * @return A vector representing the average of the input vectors
 * component-wise.
 */
std::vector<double> averageVectors(
    const std::vector<std::vector<double>>& vectors
);


/**
 * @brief Splits a vector into a specified number of subvectors.
 * 
 * This function splits a vector into a specified number of subvectors of equal size.
 * It divides the input vector into `numSubvectors` subvectors, each containing an equal
 * number of elements. If the input vector's size is not perfectly divisible by `numSubvectors`,
 * the function may not behave as expected.
 * 
 * @param vec The input vector to be split.
 * @param numSubvectors The number of subvectors to split the input vector into.
 * @return std::vector<std::vector<double>> A vector of vectors containing the split subvectors.
 */
std::vector<std::vector<double>> splitVector(
    const std::vector<double>& vec,
    const size_t numSubvectors
);


/**
 * @brief Flattens a vector of vectors into a single vector.
 * 
 * This function takes a vector of vectors and flattens it into a single
 * vector. It concatenates all subvectors into one contiguous vector. The order
 * of elements in the resulting vector follows the order of subvectors in the
 * input vector.
 * 
 * @param vecOfVecs The input vector of vectors to be flattened.
 * @return std::vector<double> A vector containing all elements from the input
 * subvectors concatenated.
 */
std::vector<double> concatVectors(
    const std::vector<std::vector<double>>& vecOfVecs
);


/**
 * @brief Overloaded stream insertion operator to print a matrix in a nice and
 * aligned format.
 * 
 * This operator allows printing a two-dimensional matrix in a nicely aligned
 * format when used with the standard output stream (`std::cout`). It
 * calculates the maximum width of elements in each column to determine the
 * width needed for alignment, then prints each element with the appropriate
 * width using `std::setw()`. Spacing between elements within a row is added to
 * ensure proper alignment.
 * 
 * @param os The output stream where the matrix will be printed.
 * @param matrix The matrix to be printed.
 * @return std::ostream& A reference to the output stream.
 */
std::ostream& operator<<(
    std::ostream& os,
    const std::vector<std::vector<double>>& matrix
);


/**
 * @brief Converts a one-dimensional vector into a two-dimensional matrix.
 * 
 * This function takes a one-dimensional vector and reshapes it into a square
 * matrix of the specified size. The elements of the input vector are filled
 * into the matrix row-wise.
 * 
 * @param vec Input vector to be converted into a matrix.
 * @param matrixSize Size of the square matrix to be generated.
 * @return std::vector<std::vector<double>> A two-dimensional vector
 * representing the converted matrix.
 */
std::vector<std::vector<double>> vector2matrix(
    const std::vector<double> &vec,
    size_t matrixSize
);


/**
 * @brief Converts a two-dimensional matrix into a one-dimensional vector
 * row-wise.
 * 
 * This function takes a two-dimensional matrix and flattens it into a
 * one-dimensional vector row-wise. The elements of the matrix are filled into
 * the vector sequentially by traversing each row.
 * 
 * @param matrix Input matrix to be flattened into a vector.
 * @param matrixSize Size of the square matrix (number of rows or columns).
 * @return std::vector<double> A one-dimensional vector representing the
 * flattened matrix.
 */
std::vector<double> matrix2vector(
    const std::vector<std::vector<double>>& matrix,
    const size_t matrixSize
);


/**
 * @brief Loads a one-dimensional vector from a CSV file.
 * 
 * This function reads a one-dimensional vector from a CSV file. It reads the
 * file line by line, parsing each line as a double value and storing it in the
 * output vector. The function stops reading the file when the output vector is
 * filled or when the end of the file is reached.
 * 
 * @param vectorLength Length of the output vector to be read from the file.
 * @return std::vector<double> A one-dimensional vector containing the values
 * read from the file.
 */
std::vector<double> loadPoints1D(
    const size_t vectorLength
);
