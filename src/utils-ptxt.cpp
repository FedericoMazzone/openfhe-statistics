#include "utils-ptxt.h"


std::vector<double> averageVectors(
    const std::vector<std::vector<double>>& vectors
)
{
    std::vector<double> average(vectors[0].size(), 0.0);
    for (const auto& vec : vectors)
        for (size_t i = 0; i < vec.size(); i++)
            average[i] += vec[i];
    for (auto& val : average)
        val /= vectors.size();

    return average;
}


std::vector<std::vector<double>> splitVector(
    const std::vector<double>& vec,
    const size_t numSubvectors
)
{
    size_t subSize = vec.size() / numSubvectors;
    std::vector<std::vector<double>> result;
    auto it = vec.begin();
    for (size_t i = 0; i < numSubvectors; ++i)
    {
        result.push_back(std::vector<double>(it, it + subSize));
        it += subSize;
    }

    return result;
}


std::vector<double> concatVectors(
    const std::vector<std::vector<double>>& vecOfVecs
)
{
    size_t totalSize = 0;
    for (const auto& subvec : vecOfVecs)
        totalSize += subvec.size();
    std::vector<double> result(totalSize);
    size_t index = 0;
    for (const auto& subvec : vecOfVecs)
        for (const auto& elem : subvec)
            result[index++] = elem;

    return result;
}


std::ostream& operator<<(
    std::ostream& os,
    const std::vector<std::vector<double>>& matrix
)
{
    size_t numRows = matrix.size();
    size_t numCols = matrix[0].size();

    // Find the maximum width of elements in each column
    std::vector<size_t> maxWidths(numCols, 0);
    for (size_t i = 0; i < numRows; ++i) {
        for (size_t j = 0; j < numCols; ++j) {
            size_t width = std::to_string(matrix[i][j]).length();
            if (width > maxWidths[j]) {
                maxWidths[j] = width;
            }
        }
    }

    // Print the matrix
    os << std::endl;
    for (size_t i = 0; i < numRows; ++i) {
        for (size_t j = 0; j < numCols; ++j) {
            os << std::setw(maxWidths[j]) << matrix[i][j];
            if (j != numCols - 1) {
                os << "  "; // Adjust this spacing as needed
            }
        }
        os << std::endl;
    }

    return os;
}


std::vector<std::vector<double>> vector2matrix(
    const std::vector<double> &vec,
    const size_t matrixSize
)
{
    return splitVector(vec, matrixSize);
}


std::vector<double> matrix2vector(
    const std::vector<std::vector<double>>& matrix,
    const size_t matrixSize
)
{
    std::vector<double> vec(matrixSize * matrixSize);

    for (size_t i = 0; i < matrixSize; ++i) {
        for (size_t j = 0; j < matrixSize; ++j) {
            vec[i * matrixSize + j] = matrix[i][j];
        }
    }

    return vec;
}