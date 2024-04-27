#include "utils-eval.h"



usint depth2degree(
    const usint depth
)
{
    switch(depth)
    {
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
