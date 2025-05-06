#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <random>
#include <algorithm>
#include <map>
#include <fstream>
#include <limits>

// Класс QAM модулятора
class QAMModulator {
public:
    QAMModulator(int bitsPerSymbol) : bitsPerSymbol(bitsPerSymbol) {
        initializeConstellation();
    }

    std::vector<std::complex<double>> modulate(const std::vector<int>& bits) {
        std::vector<std::complex<double>> symbols;
        for (size_t i = 0; i < bits.size(); i += bitsPerSymbol) {
            int symbolBits = 0;
            for (int j = 0; j < bitsPerSymbol && (i + j) < bits.size(); ++j) {
                symbolBits = (symbolBits << 1) | bits[i + j];
            }
            if (constellation.find(symbolBits) != constellation.end()) {
                symbols.push_back(constellation.at(symbolBits));
            }
        }
        return symbols;
    }

private:
    int bitsPerSymbol;
    std::map<int, std::complex<double>> constellation;

    void initializeConstellation() {
        switch (bitsPerSymbol) {
        case 2: {  // QPSK
            constellation = {
                {0, {-1.0, -1.0}}, {1, {-1.0, 1.0}},
                {2, {1.0, -1.0}},  {3, {1.0, 1.0}}
            };
            // Нормализация энергии
            for (auto& [key, val] : constellation) {
                val /= std::sqrt(2.0);
            }
            break;
        }
        case 4: {  // 16-QAM
            const double levels[] = { -3.0, -1.0, 1.0, 3.0 };
            const double normalizationFactor = std::sqrt(10.0);
            for (int i = 0; i < 16; ++i) {
                int iBits = (i >> 2) & 0x3;
                int qBits = i & 0x3;
                double iLevel = levels[iBits] / normalizationFactor;
                double qLevel = levels[qBits] / normalizationFactor;
                constellation[i] = std::complex<double>(iLevel, qLevel);
            }
            break;
        }
        case 6: {  // 64-QAM
            const double levels[] = { -7.0, -5.0, -3.0, -1.0, 1.0, 3.0, 5.0, 7.0 };
            const double normalizationFactor = std::sqrt(42.0);
            for (int i = 0; i < 64; ++i) {
                int iBits = (i >> 3) & 0x7;
                int qBits = i & 0x7;
                double iLevel = levels[iBits] / normalizationFactor;
                double qLevel = levels[qBits] / normalizationFactor;
                constellation[i] = std::complex<double>(iLevel, qLevel);
            }
            break;
        }
        default:
            throw std::invalid_argument("Unsupported bits per symbol");
        }
    }
};

// Класс для добавления гауссовского шума
class AWGNChannel {
public:
    AWGNChannel(double noiseVariance) : noiseVariance(noiseVariance) {
        std::random_device rd;
        gen = std::mt19937(rd());
        dist = std::normal_distribution<>(0.0, 1.0);
    }

    std::vector<std::complex<double>> addNoise(const std::vector<std::complex<double>>& symbols) {
        std::vector<std::complex<double>> noisySymbols;
        for (const auto& symbol : symbols) {
            double noiseI = dist(gen) * std::sqrt(noiseVariance / 2.0);
            double noiseQ = dist(gen) * std::sqrt(noiseVariance / 2.0);
            noisySymbols.emplace_back(symbol.real() + noiseI, symbol.imag() + noiseQ);
        }
        return noisySymbols;
    }

private:
    double noiseVariance;
    std::mt19937 gen;
    std::normal_distribution<> dist;
};

// Класс QAM демодулятора
class QAMDemodulator {
public:
    QAMDemodulator(int bitsPerSymbol) : bitsPerSymbol(bitsPerSymbol) {
        initializeConstellation();
    }

    std::vector<int> demodulate(const std::vector<std::complex<double>>& symbols) {
        std::vector<int> bits;
        for (const auto& symbol : symbols) {
            int bestSymbol = findClosestSymbol(symbol);
            for (int i = bitsPerSymbol - 1; i >= 0; --i) {
                bits.push_back((bestSymbol >> i) & 1);
            }
        }
        return bits;
    }

private:
    int bitsPerSymbol;
    std::map<int, std::complex<double>> constellation;

    void initializeConstellation() {
        switch (bitsPerSymbol) {
        case 2: {  // QPSK
            constellation = {
                {0, {-1.0, -1.0}}, {1, {-1.0, 1.0}},
                {2, {1.0, -1.0}},  {3, {1.0, 1.0}}
            };
            for (auto& [key, val] : constellation) {
                val /= std::sqrt(2.0);
            }
            break;
        }
        case 4: {  // 16-QAM
            const double levels[] = { -3.0, -1.0, 1.0, 3.0 };
            const double normalizationFactor = std::sqrt(10.0);
            for (int i = 0; i < 16; ++i) {
                int iBits = (i >> 2) & 0x3;
                int qBits = i & 0x3;
                double iLevel = levels[iBits] / normalizationFactor;
                double qLevel = levels[qBits] / normalizationFactor;
                constellation[i] = std::complex<double>(iLevel, qLevel);
            }
            break;
        }
        case 6: {  // 64-QAM
            const double levels[] = { -7.0, -5.0, -3.0, -1.0, 1.0, 3.0, 5.0, 7.0 };
            const double normalizationFactor = std::sqrt(42.0);
            for (int i = 0; i < 64; ++i) {
                int iBits = (i >> 3) & 0x7;
                int qBits = i & 0x7;
                double iLevel = levels[iBits] / normalizationFactor;
                double qLevel = levels[qBits] / normalizationFactor;
                constellation[i] = std::complex<double>(iLevel, qLevel);
            }
            break;
        }
        default:
            throw std::invalid_argument("Unsupported bits per symbol");
        }
    }

    int findClosestSymbol(const std::complex<double>& point) {
        int bestSymbol = 0;
        double minDistance = std::numeric_limits<double>::max();
        for (const auto& [symbol, constellationPoint] : constellation) {
            double dx = point.real() - constellationPoint.real();
            double dy = point.imag() - constellationPoint.imag();
            double distance = dx * dx + dy * dy;  // Квадрат расстояния (без sqrt для оптимизации)
            if (distance < minDistance) {
                minDistance = distance;
                bestSymbol = symbol;
            }
        }
        return bestSymbol;
    }
};

// Функция для расчета BER
double calculateBER(const std::vector<int>& original, const std::vector<int>& received) {
    int errors = 0;
    for (size_t i = 0; i < std::min(original.size(), received.size()); ++i) {
        if (original[i] != received[i]) {
            ++errors;
        }
    }
    return static_cast<double>(errors) / original.size();
}

// Функция симуляции
void simulateQAM(int bitsPerSymbol, const std::string& filename) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> bitDist(0, 1);

    // Генерация случайных бит
    std::vector<int> bits(100000);
    for (auto& bit : bits) {
        bit = bitDist(gen);
    }

    QAMModulator modulator(bitsPerSymbol);
    QAMDemodulator demodulator(bitsPerSymbol);

    std::vector<double> noiseVariances = { 0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0 };
    std::vector<double> berResults;

    std::ofstream outFile(filename);
    outFile << "Noise Variance,BER\n";

    for (double variance : noiseVariances) {
        AWGNChannel channel(variance);

        auto symbols = modulator.modulate(bits);
        auto noisySymbols = channel.addNoise(symbols);
        auto receivedBits = demodulator.demodulate(noisySymbols);

        double ber = calculateBER(bits, receivedBits);
        berResults.push_back(ber);

        outFile << variance << "," << ber << "\n";
        std::cout << "Variance: " << variance << ", BER: " << ber << std::endl;
    }
    outFile.close();
}

int main() {
    // Симуляция для QPSK (2 бита на символ)
    std::cout << "=== QPSK Simulation ===" << std::endl;
    simulateQAM(2, "qpsk_results.csv");

    // Симуляция для 16-QAM (4 бита на символ)
    std::cout << "=== 16-QAM Simulation ===" << std::endl;
    simulateQAM(4, "16qam_results.csv");

    // Симуляция для 64-QAM (6 бит на символ)
    std::cout << "=== 64-QAM Simulation ===" << std::endl;
    simulateQAM(6, "64qam_results.csv");

    return 0;
}