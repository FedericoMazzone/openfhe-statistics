#include "utils-eval.h"



usint depth2degree(
    const usint depth
)
{
    switch(depth)
    {
        case 3:     return 2;
        case 4:     return 5;
        case 5:     return 13;
        case 6:     return 27;
        case 7:     return 59;
        case 8:     return 119;
        case 9:     return 247;
        case 10:    return 495;
        case 11:    return 1007;
        case 12:    return 2031;

        case 13:    return 4031;
        case 14:    return 8127;
        default:    return -1;
    }
}


Ciphertext<DCRTPoly> equal(
    const Ciphertext<DCRTPoly> &c1,
    const Ciphertext<DCRTPoly> &c2,
    double a,
    double b,
    uint32_t degree,
    double error
)
{
    return c1->GetCryptoContext()->EvalChebyshevFunction(
        [error](double x) -> double { 
            if      (x > error)   return 0;
            else if (x >= -error) return 1.0;
            else                  return 0;
        },
        c1 - c2,
        a, b, degree
    );
}


Ciphertext<DCRTPoly> compare(
    const Ciphertext<DCRTPoly> &c1,
    const Ciphertext<DCRTPoly> &c2,
    double a,
    double b,
    uint32_t degree,
    double error
)
{
    return c1->GetCryptoContext()->EvalChebyshevFunction(
        [error](double x) -> double { 
            if      (x > error)   return 1;
            else if (x >= -error) return 0.5;
            else                  return 0;
        },
        c1 - c2,
        a, b, degree
    );
}


Ciphertext<DCRTPoly> signAdv(
    Ciphertext<DCRTPoly> &c,
    const size_t dg,
    const size_t df
)
{
    // f3(x) = (35 x - 35 x^3 + 21 x^5 - 5 x^7) / 16
    std::vector<double> coeffF3 = {0, 35.0 / 16.0, 0, -35.0 / 16.0, 0, 21.0 / 16.0, 0, -5.0 / 16.0};
    std::vector<double> coeffF3_final = {0.5, 35.0 / 32.0, 0, -35.0 / 32.0, 0, 21.0 / 32.0, 0, -5.0 / 32.0};
    // g3(x) = (4589 x - 16577 x^3 + 25614 x^5 - 12860 x^7) / 1024
    std::vector<double> coeffG3 = {0, 4589.0 / 1024.0, 0, -16577.0 / 1024.0, 0, 25614.0 / 1024.0, 0, -12860.0 / 1024.0};

    for (size_t d = 0; d < dg; d++)
        c = c->GetCryptoContext()->EvalPolyLinear(c, coeffG3);
    for (size_t d = 0; d < df - 1; d++)
        c = c->GetCryptoContext()->EvalPolyLinear(c, coeffF3);
    c = c->GetCryptoContext()->EvalPolyLinear(c, coeffF3_final);
    return c;
}


Ciphertext<DCRTPoly> compareAdv(
    const Ciphertext<DCRTPoly> &c1,
    const Ciphertext<DCRTPoly> &c2,
    const size_t dg,
    const size_t df
)
{
    auto c = c1 - c2;
    return signAdv(c, dg, df);
}


Ciphertext<DCRTPoly> compareGt(
    const Ciphertext<DCRTPoly> &c1,
    const Ciphertext<DCRTPoly> &c2,
    double a,
    double b,
    uint32_t degree,
    double error
)
{
    return c1->GetCryptoContext()->EvalChebyshevFunction(
        [error](double x) -> double { 
            if      (x > error)   return 1;
            else                  return 0;
        },
        c1 - c2,
        a, b, degree
    );
}


Ciphertext<DCRTPoly> indicator(
    const Ciphertext<DCRTPoly> &c,
    double a1,
    double b1,
    double a,
    double b,
    uint32_t degree
)
{
    return c->GetCryptoContext()->EvalChebyshevFunction(
        [a1,b1](double x) -> double {
            return (x < a1 || x > b1) ? 0 : 1; },
        c,
        a, b, degree
    );
}


Ciphertext<DCRTPoly> indicatorAdv(
    const Ciphertext<DCRTPoly> &c,
    const double b,
    const size_t dg,
    const size_t df
)
{
    // // auto r = c->GetCryptoContext()->EvalChebyshevFunction(
    // //     [](double x) -> double {
    // //         return std::sin((59.69)*(x*x*x*x)) / ((59.69)*(x*x*x*x)); },
    // //     c,
    // //     a, b, degree
    // // );
    // // r = r * r;

    // // std::vector<double> coeffF3 = {0, 35.0 / 16.0, 0, -35.0 / 16.0, 0, 21.0 / 16.0, 0, -5.0 / 16.0};
    // // r = r->GetCryptoContext()->EvalPolyLinear(r, coeffF3);
    // // std::cout << r->GetLevel() << std::endl;
    // // return r;

    auto c1 = (1.0 / b) * (c + 0.5);
    auto c2 = (1.0 / b) * (c - 0.5);
    c1 = signAdv(c1, dg, df);
    c2 = signAdv(c2, dg, df);
    return c1 * (1 - c2);
}


// int main()
// {

//     const usint vectorLength            = 4;
//     const usint functionDepth           = 10;

//     const usint integralPrecision       = 10;
//     const usint decimalPrecision        = 50;
//     const usint multiplicativeDepth     = functionDepth;
//     const usint numSlots                = vectorLength;
//     const bool enableBootstrap          = false;
//     const usint ringDim                 = 0;
//     const bool verbose                  = true;

//     std::vector<int32_t> indices = {};

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

//     std::vector<double> v = {-0.75, -0.1, 0.0, 0.75};
//     std::vector<double> zero = {0.0, 0.0, 0.0, 0.0};

//     std::cout << "Vector: " << v << std::endl;

//     Ciphertext<DCRTPoly> vC = cryptoContext->Encrypt(
//         keyPair.publicKey,
//         cryptoContext->MakeCKKSPackedPlaintext(v)
//     );
//     Ciphertext<DCRTPoly> zeroC = cryptoContext->Encrypt(
//         keyPair.publicKey,
//         cryptoContext->MakeCKKSPackedPlaintext(zero)
//     );

//     auto start = std::chrono::high_resolution_clock::now();

//     // Ciphertext<DCRTPoly> resultC = compare(
//     //     zeroC, vC,
//     //     -1.0, 1.0,
//     //     depth2degree(functionDepth)
//     // );
//     Ciphertext<DCRTPoly> resultC = indicator(
//         vC,
//         -0.5, 0.5,
//         -1.0, 1.0,
//         depth2degree(functionDepth)
//     );

//     auto end = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> elapsed_seconds = end - start;
//     std::cout << elapsed_seconds.count() << "s" << std::endl;

//     Plaintext resultP;
//     cryptoContext->Decrypt(keyPair.secretKey, resultC, &resultP);
//     resultP->SetLength(vectorLength);

//     std::vector<double> result = resultP->GetRealPackedValue();
//     std::cout << "Result: " << result << std::endl;

//     return 0;

// }
