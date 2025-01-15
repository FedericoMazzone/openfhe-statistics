#include "utils-ptxt.h"

#include <fstream>
#include <sstream>
#include <string>


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


std::vector<double> loadPoints1D(
    const size_t vectorLength
)
{
    std::vector<double> v(vectorLength, 0.0);
    std::ifstream file("data/points1d.csv");

    if (!file.is_open())
    {
        file.open("../data/points1d.csv");
        if (!file.is_open())
            std::cerr << "Error: Could not open the file!" << std::endl;
    }

    std::string line;
    size_t count = 0;

    while (std::getline(file, line) && count < vectorLength)
    {
        try
        {
            v[count] = std::stod(line);
            count++;
        } catch (const std::exception& e)
        {
            std::cerr << "Error: Invalid number format in the file." << std::endl;
        }
    }

    file.close();

    return v;
}
