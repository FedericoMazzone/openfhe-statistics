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
    const double a,
    const double b,
    const uint32_t degree,
    const double error
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
    const double a,
    const double b,
    const uint32_t degree,
    const double error
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
    const double a,
    const double b,
    const uint32_t degree,
    const double error
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
    const double a1,
    const double b1,
    const double a,
    const double b,
    const uint32_t degree
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
    auto tmp = (1.0 / b) * c;
    auto c1 = tmp + 0.5 / b;
    auto c2 = tmp - 0.5 / b;
    c1 = signAdv(c1, dg, df);
    c2 = signAdv(c2, dg, df);
    return c1 * (1 - c2);
}


Ciphertext<DCRTPoly> indicatorAdvShifted(
    const Ciphertext<DCRTPoly> &c,
    const double b,
    const size_t dg,
    const size_t df
)
{
    auto c1 = (2.0 / (b + 1)) * c + 2.0 / (b + 1) - 1.0;
    auto c2 = (-2.0 / (b + 1)) * c + 2.0 / (b + 1) + 1.0;
    c1 = signAdv(c1, dg, df);
    c2 = signAdv(c2, dg, df);
    return c1 * c2;
}
