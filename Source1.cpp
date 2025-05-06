#include "matplotlibcpp.h"
#include <vector>
#include <fstream>
#include <sstream>

namespace plt = matplotlibcpp;

void plotBER() {
    // Чтение данных из CSV
    std::ifstream file("qpsk_results.csv");
    std::vector<double> variances, bers;
    std::string line;

    std::getline(file, line); // Пропускаем заголовок
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        double var, ber;
        char comma;
        ss >> var >> comma >> ber;
        variances.push_back(var);
        bers.push_back(ber);
    }

    // Построение графика
    plt::semilogy(variances, bers, "o-");
    plt::title("Зависимость BER от дисперсии шума (QPSK)");
    plt::xlabel("Дисперсия шума");
    plt::ylabel("BER");
    plt::grid(true);
    plt::save("ber_plot.png");
    plt::show();
}

int main() {
    plotBER();
    return 0;
}