#include "utils-basics.h"



CryptoContext<DCRTPoly> generateCryptoContext(
    const usint integralPrecision,
    const usint decimalPrecision,
    const usint multiplicativeDepth,
    const usint numSlots,
    const bool enableBootstrap,
    const usint ringDim,
    const bool verbose
)
{

    const usint plaintextPrecision = integralPrecision + decimalPrecision;
    const SecretKeyDist secretKeyDist = UNIFORM_TERNARY;

    const std::vector<usint> levelBudget = (numSlots <= 8) ?
            std::vector<usint>({3, 3}) : std::vector<usint>({4, 4});
    const usint approxBootstrapDepth = 8;
    const std::vector<usint> bsgsDim = {0, 0};
    const usint depth = enableBootstrap ?
        multiplicativeDepth + FHECKKSRNS::GetBootstrapDepth(
            approxBootstrapDepth, levelBudget, secretKeyDist) :
        multiplicativeDepth;

    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetSecretKeyDist(secretKeyDist);
    parameters.SetScalingModSize(decimalPrecision);
    parameters.SetScalingTechnique(ScalingTechnique::FLEXIBLEAUTO);
    parameters.SetFirstModSize(plaintextPrecision);
    parameters.SetMultiplicativeDepth(depth);
    parameters.SetKeySwitchTechnique(KeySwitchTechnique::HYBRID);
    parameters.SetBatchSize(numSlots);
    if (ringDim == 0)
    {
        parameters.SetSecurityLevel(HEStd_128_classic);
    }
    else
    {
        parameters.SetSecurityLevel(HEStd_NotSet);
        parameters.SetRingDim(ringDim);
    }
    
    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);
    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    cryptoContext->Enable(ADVANCEDSHE);
    if (enableBootstrap)
    {
        cryptoContext->Enable(FHE);
        cryptoContext->EvalBootstrapSetup(levelBudget, bsgsDim, numSlots);
    }

    if (verbose)
    {
        const BigInteger ciphertextModulus = cryptoContext->GetModulus();
        const usint ciphertextModulusBitsize = ciphertextModulus.GetLengthForBase(2);
        const usint ringDimension = cryptoContext->GetRingDimension();
        const usint maxNumSlots = ringDimension / 2;
        const auto elementParameters = cryptoContext->GetCryptoParameters()->GetElementParams()->GetParams();
        std::vector<NativeInteger> moduliChain(elementParameters.size());
        std::vector<usint> moduliChainBitsize(elementParameters.size());
        for (size_t i = 0; i < elementParameters.size(); i++)
        {
            moduliChain[i] = elementParameters[i]->GetModulus();
            moduliChainBitsize[i] = moduliChain[i].GetLengthForBase(2);
        }

        std::ostringstream logMessage;
        logMessage << "CKKS PARAMETERS"                                                              << std::endl;
        logMessage << "Integral Bit Precision        : " << integralPrecision                        << std::endl;
        logMessage << "Decimal Bit Precision         : " << decimalPrecision                         << std::endl;
        logMessage << "Ciphertext Modulus Precision  : " << ciphertextModulusBitsize                 << std::endl;
        logMessage << "Ring Dimension                : " << ringDimension                            << std::endl;
        logMessage << "Max Slots                     : " << maxNumSlots                              << std::endl;
        logMessage << "Slots                         : " << numSlots                                 << std::endl;
        logMessage << "Multiplicative Depth          : " << parameters.GetMultiplicativeDepth()      << std::endl;
        logMessage << "Security Level                : " << parameters.GetSecurityLevel()            << std::endl;
        logMessage << "Secret Key Distribution       : " << parameters.GetSecretKeyDist()            << std::endl;
        logMessage << "Scaling Technique             : " << parameters.GetScalingTechnique()         << std::endl;
        logMessage << "Encryption Technique          : " << parameters.GetEncryptionTechnique()      << std::endl;
        logMessage << "Multiplication Technique      : " << parameters.GetMultiplicationTechnique()  << std::endl;
        logMessage << "Moduli Chain Bitsize          : " << moduliChainBitsize                       << std::endl;
        // logMessage << "Moduli Chain                  : " << moduliChain                              << std::endl;
        logMessage << std::endl;

        std::cout << logMessage.str();
    }

    return cryptoContext;

}


KeyPair<DCRTPoly> keyGeneration(
    const CryptoContext<DCRTPoly> &cryptoContext,
    const std::vector<int32_t> &indices,
    const usint numSlots,
    const bool enableBootstrap,
    const bool verbose
)
{

    auto start = std::chrono::high_resolution_clock::now();

    if (verbose)
        std::cout << "Key generation protocol...          " << std::flush;
    
    // Public & private key generation
    auto keyPair = cryptoContext->KeyGen();

    // Relinearization key generation
    cryptoContext->EvalMultKeyGen(keyPair.secretKey);

    // Rotation keyPair generation
    cryptoContext->EvalRotateKeyGen(keyPair.secretKey, indices);

    if (enableBootstrap)
    {
        // Bootstrap key generation
        cryptoContext->EvalBootstrapKeyGen(keyPair.secretKey, numSlots);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    if (verbose)
        std::cout << "COMPLETED (" <<
            elapsed_seconds.count() << "s)" << std::endl << std::endl;
    
    return keyPair;

}
