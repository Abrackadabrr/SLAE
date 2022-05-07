//
// Created by evgen on 20.03.2022.
//

#ifndef MY_PROJECT_DIRIHLETPIOSSON_H
#define MY_PROJECT_DIRIHLETPIOSSON_H

#include "my_project/matrix/sparse/CSRmatrix.h"
#include "set"
#include "my_project/solvers/SteepestGradiendDescent.h"

namespace Slae::Examples {
    std::pair<Slae::Matrix::CSR<double>, std::vector<double>>
    PoissonDirichletProblemSnake(unsigned size_n, double step,
                                 const std::function<double(double, double)> &border_value,
                                 const std::function<double(double, double)> &lap_value) {

        std::set<Slae::Matrix::Element<double>> elements{};
        std::vector<double> b(size_n * size_n, 0);
        for (unsigned i = 0; i < size_n; i++) {
            elements.insert({i, i, 1.});
            b.at(i) = border_value(i * step, 0);

            elements.insert({size_n * (size_n - 1) + i, size_n * (size_n - 1) + i, 1.});

            b.at(size_n * (size_n - 1) + i) = border_value((size_n - i - 1) * step, (size_n - 1) * step);
            if (size_n % 2 == 1) b.at(size_n * (size_n - 1) + i) = border_value(i * step, (size_n - 1) * step);

            elements.insert({i * size_n, i * size_n, 1.});
            b.at(i * size_n) = border_value((size_n - 1) * step, i * step);
            if (i % 2 == 0) b.at(i * size_n) = border_value(0, i * step);

            if (i != 0) {
                elements.insert({i * size_n - 1, i * size_n - 1, 1.});
                b.at(i * size_n - 1) = border_value(0, (i - 1) * step);
                if (i % 2 == 1) b.at(i * size_n - 1) = border_value((size_n - 1) * step, (i - 1) * step);
            }
        }
        for (int i = 1; i < size_n - 1; i++) {
            for (int j = 1; j < size_n - 1; j++) {
                unsigned number = i * size_n + j;
                elements.insert({number, number, -4});
                elements.insert({number, number - 1, 1});
                elements.insert({number, number + 1, 1});
                elements.insert({number, 2 * size_n * i - 1 - number, 1});
                elements.insert({number, 2 * size_n * (i + 1) - 1 - number, 1});
                b.at(number) = lap_value(j * step, i * step) * step * step;
                if (i % 2 == 1) b.at(number) = lap_value((size_n - 1 - j) * step, i * step) * step * step;
            }
        }
        auto A = Slae::Matrix::CSR<double>(size_n * size_n, size_n * size_n, elements);
        return {A, b};
    }


    std::pair<Slae::Matrix::CSR<double>, std::vector<double>>
    PoissonDirichletProblemWithBoarders(unsigned size_n, double step,
                                        const std::function<double(double, double)> &border_value,
                                        const std::function<double(double, double)> &lap_value) {

        std::set<Slae::Matrix::Element<double>> elements{};
        std::vector<double> b(size_n * size_n, 0);
        for (unsigned i = 0; i < size_n; i++) {
            elements.insert({i, i, 1.});
            b.at(i) = border_value(i * step, 0);

            elements.insert({size_n * (size_n - 1) + i, size_n * (size_n - 1) + i, 1.});
            b.at(size_n * (size_n - 1) + i) = border_value(i * step, (size_n - 1) * step);

            elements.insert({i * size_n, i * size_n, 1.});
            b.at(i * size_n) = border_value(0, i * step);

            if (i != 0) {
                elements.insert({i * size_n - 1, i * size_n - 1, 1.});
                b.at(i * size_n - 1) = border_value((size_n - 1) * step, (i - 1) * step);
            }
        }
        for (int i = 1; i < size_n - 1; i++) {
            for (int j = 1; j < size_n - 1; j++) {
                unsigned number = i * size_n + j;
                elements.insert({number, number, -4});
                elements.insert({number, number - 1, 1});
                elements.insert({number, number + 1, 1});
                elements.insert({number, number + size_n, 1});
                elements.insert({number, number - size_n, 1});
                b.at(number) = lap_value(j * step, i * step) * step * step;
            }
        }
        auto A = Slae::Matrix::CSR<double>(size_n * size_n, size_n * size_n, elements);
        return {A, b};
    }


    std::pair<Slae::Matrix::CSR<double>, std::vector<double>>
    PoissonDirichletProblem(unsigned size, double step,
                            const std::function<double(double, double)> &border_value,
                            const std::function<double(double, double)> &lap_value) {
        std::set<Slae::Matrix::Element<double>> elements{};
        const unsigned n = size - 2;
        std::vector<double> b(n * n, 0);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                unsigned number = i * n + j;
                elements.insert({number, number, -4});

                if (j != 0) {
                    elements.insert({number, number - 1, 1});
                } else {
                    b.at(number) -= border_value(0, (i + 1) * step);
                }

                if (j != n-1) {
                    elements.insert({number, number + 1, 1});
                } else {
                    b.at(number) -= border_value((n+1) * step, (i + 1) * step);
                }

                if (i != n-1) {
                    elements.insert({number, number + n, 1});
                } else {
                    b.at(number) -= border_value((j + 1) * step, (n+1) * step);
                }

                if (i != 0) {
                    elements.insert({number, number - n, 1});
                } else {
                    b.at(number) -= border_value((j + 1) * step, 0);
                }

                b.at(number) += lap_value((j + 1) * step, (i + 1) * step) * step * step;
            }
        }
        auto A = Slae::Matrix::CSR<double>(n*n, n*n, elements);
        return {A, b};
    }
}

#endif //MY_PROJECT_DIRIHLETPIOSSON_H
