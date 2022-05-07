//
// Created by evgen on 28.04.2021.
//

#ifndef MATRIX_CSRMATRIX_H
#define MATRIX_CSRMATRIX_H

#include <vector>
#include <iomanip>
#include <algorithm>
#include <my_project/vectors/Norm.hpp>
#include "my_project/matrix/sparse/CSRExeption.h"
#include "my_project/vectors/Vectors.hpp"
#include "my_project/matrix/utility/Element.hpp"
#include "set"

using std::vector;
using std::ostream;

namespace Slae::Matrix {
    template<typename value_t>
    class CSR {
    public:
        using ind_t = unsigned;

        template<typename T>
        friend ostream &operator<<(ostream &os, const CSR<T> &matrix);

        template<typename T>
        friend ostream &operator<<(ostream &os, const CSR<T> &&matrix);

        template<typename T>
        friend CSR operator*(value_t k, const CSR<T> &matrix);

//        template<typename T>
//        friend vector<T> Slae::Solvers::GaussSeidelIteration(const int &A, const std::vector<T> &b, std::vector<T> x) {}

    public:
        vector<value_t> _value;
        vector<ind_t> _col;
        vector<ind_t> _row_index;

        ind_t _height;
        ind_t _width;

    public:
        static CSR<value_t> Diag(value_t element, unsigned N) {
            std::vector<value_t> value{};
            std::vector<ind_t> row_index{};
            std::vector<ind_t> col{};
            row_index.resize(0);
            col.resize(0);
            value.resize(0);

            for (unsigned i = 0; i < N; ++i) {
                value.push_back(element);
                col.push_back(i);
                row_index.push_back(i);
            }
            row_index.push_back(N);
            return std::move(CSR<value_t>(value, col, row_index, N, N));
        }

        [[nodiscard]] unsigned int height() const {
            return _height;
        }

        [[nodiscard]] unsigned int width() const {
            return _width;
        }

        [[nodiscard]] ind_t n_values() const {return _value.size();};

        CSR(const vector<value_t> &value, const vector<ind_t> &col, const vector<ind_t> &row_index, ind_t H,
            ind_t W) : _width(W), _height(H) {
            if (row_index.size() != H + 1)
                throw Slae::Matrix::CSRException("Invalid creation CSR (row_index.size != martrix height + 1)");
            if (value.size() != col.size())
                throw Slae::Matrix::CSRException("Invalid creation CSR (udefined amount of not zero elements)");
            if (row_index[H] != col.size())
                throw Slae::Matrix::CSRException("Invalid creation CSR (incorrect data)");
            if (*std::max_element(col.begin(), col.end()) >= W)
                throw std::runtime_error("Invalid creation CSR (incorrect data)");

            _col = col;
            _row_index = row_index;
            _value = value;
        }

        CSR(const ind_t &h, const ind_t &w, const std::set<Slae::Matrix::Element<value_t>> &in) : _height(h), _width(w) {
            _value.resize(in.size());
            _col.resize(in.size());
            _row_index.resize(h + 1, 0);
            int countInRow = 0;
            int currRow = 0;
            auto it = in.begin();
            for (ind_t k = 0; k < in.size(); ++k) {
                while (currRow < it->i) {
                    _row_index[currRow + 1] = _row_index[currRow] + countInRow;
                    ++currRow;
                    countInRow = 0;
                }
                _value[k] = it->value;
                _col[k] = it->j;
                ++countInRow;
                it = std::next(it);
            }
            for (++currRow; currRow <= _height; ++currRow) {
                _row_index[currRow] = in.size();
            }
        }

        CSR() : _col{}, _value{}, _row_index{}, _height(0), _width(0) {};

        CSR(const CSR<value_t> &matrix) = default;

        ~CSR() = default;

        void swap(CSR<value_t> &_m) noexcept {
            _value.swap(_m._value);
            _col.swap(_m._col);
            _row_index.swap(_m._row_index);
            ind_t tmp = _height;
            _height = _m._height;
            _m._height = tmp;
            tmp = _width;
            _width = _m._width;
            _m._width = tmp;
        };

        CSR &operator=(const CSR<value_t> &matrix) {
            if (this != &matrix) {
                CSR<value_t>(matrix).swap(*this);
            }
            return *this;
        }

        value_t operator()(unsigned i, unsigned j) const {
            for (unsigned num_of_elem = _row_index[i]; num_of_elem < _row_index[i + 1]; num_of_elem++) {
                if (_col[num_of_elem] == j)
                    return _value[num_of_elem];
            }
            return 0;
        }

        CSR operator+(const CSR &first_matrix) const {
            std::vector<unsigned> using_cols(_width, -1);  // vector for mark used cols
            std::vector<unsigned> tmp_cols;  // temporary сol
            std::vector<ind_t> cols;  // col
            std::vector<ind_t> row_ind(_height + 1);
            std::vector<value_t> tmp_value;  // temporary vector contains values of added rows
            std::vector<value_t> val;
            row_ind[0] = 0;
            // merge column indices of non-null elements into one row array i
            for (size_t i = 0; i < _height; ++i) {
                tmp_value.resize(_width);
                for (unsigned j = _row_index[i]; j < _row_index[i + 1]; ++j) {
                    if (using_cols[_col[j]] != i) {
                        using_cols[_col[j]] = i;
                        tmp_cols.emplace_back(_col[j]);
                    }
                    tmp_value[_col[j]] = _value[j];
                }
                for (unsigned j = first_matrix._row_index[i]; j < first_matrix._row_index[i + 1]; ++j) {
                    if (using_cols[first_matrix._col[j]] != i) {
                        using_cols[first_matrix._col[j]] = i;
                        tmp_cols.emplace_back(first_matrix._col[j]);
                    }
                    tmp_value[first_matrix._col[j]] += first_matrix._value[j];
                }
                for (const auto &j: tmp_cols) {
                    if (std::abs(tmp_value[j]) > std::numeric_limits<double>::epsilon()) {
                        cols.emplace_back(j);
                        val.emplace_back(tmp_value[j]);
                    } else continue;
                }
                row_ind[i + 1] = cols.size();
                tmp_cols.clear();
                tmp_value.clear();
            }
            CSR newmat(val, cols, row_ind, _height, _width);
            return std::move(newmat);
        }

        CSR operator-(const CSR &matrix) const {
            return *this + (-1) * matrix;
        }

        CSR operator*=(value_t k) {
            for (int i = 0; i < _value.size(); i++) {
                _value[i] *= k;
            }
            return *this;
        }

        CSR operator*(value_t k) const {
            vector<value_t> n_val;
            n_val.reserve(_value.size());
            for (int i = 0; i < _value.size(); i++) {
                n_val.emplace_back(_value[i] * k);
            }
            return CSR(n_val, _col, _row_index, _height, _width);
        }

        friend CSR operator*(value_t k, const CSR &matrix) {
            return matrix * k;
        }

        vector<value_t> operator*(const vector<value_t> &vector) const {
            std::vector<value_t> multi;
            for (unsigned int str_number = 0; str_number < _height; ++str_number) {
                value_t element = 0;
                for (unsigned int num_of_elem = _row_index[str_number];
                     num_of_elem < _row_index[str_number + 1]; ++num_of_elem)
                    element += _value[num_of_elem] * vector[_col[num_of_elem]];
                multi.push_back(element);
            }
            return multi;
        }

        void printAsCSR(ostream &os) {
            os << "Value: " << _value << std::endl;
            os << "Col: " << _col << std::endl;
            os << "Row_indexes:" << _row_index << std::endl;
        }

        CSR transpose() const {
            ind_t NonZero = this->_value.size();
            std::vector<value_t> tVals(NonZero);
            std::vector<ind_t> tCols(NonZero);
            std::vector<ind_t> tRows(_width + 1);
            for (ind_t i = 0; i < NonZero; ++i)
                tRows[_col[i] + 1]++;  //Посчитали число ненулевых элементов в каждой строке транспонированной матрицы
            ind_t S = 0;
            ind_t tmp;
            for (ind_t i = 1; i <= _width; ++i) {          //запишем в индекс каждой из строк сумму из всех предыдущих
                tmp = tRows[i];
                tRows[i] = S;
                S = S + tmp;
            }
            ind_t j1, j2, Col, RIndex, IIndex;
            value_t V;
            for (ind_t i = 0; i < _height; ++i) { //обойдем по всем строкам исходной матрицы
                j1 = _row_index[i];
                j2 = _row_index[i + 1];
                Col = i;   //столбец в транспонированной = строка в исходной
                for (ind_t j = j1; j < j2; ++j) { // обойдем по всем элементам текущей строки
                    V = _value[j];      //значение
                    RIndex = _col[j];   //номер строки транспонированной соответствующий значению
                    IIndex = tRows[RIndex + 1]; // порядковый номер элемента в транспонированной матрице
                    tVals[IIndex] = V;  // вставим значение в транспонированную матрицу
                    tCols[IIndex] = Col;    // вставим индекс, соответствующий значению
                    tRows[RIndex + 1]++;    //дополним каждую строку числом вставленных в нее элементов
                }
            }
            return CSR(tVals, tCols, tRows, _width, _height);
        }
    };

    template<typename value_t>
    ostream &operator<<(ostream &os, const CSR<value_t> &matrix) { // <- working
        for (int i = 0; i < matrix.height(); i++) {
            os << "|| ";
            for (int j = 0; j < matrix.width(); j++)
                os << std::setw(3) << matrix(i, j) << " ";
            if (i != matrix.height() - 1) os << "||" << std::endl;
            else os << "||";
        }
        os << "(" << matrix.height() << "x" << matrix.width() << ")" << std::endl;
        return os;
    }

    template<typename value_t>
    ostream &operator<<(ostream &os, const CSR<value_t> &&matrix) { // <- working
        os << matrix;
        return os;
    }

    template<typename value_t>
    void printToFile
        (ostream &os, const CSR<value_t> &matrix) { // <- working
        for (int i = 0; i < matrix.height(); i++) {
            for (int j = 0; j < matrix.width(); j++) {
                if (j != matrix.width() - 1) os << matrix(i, j) << " ";
                else os << matrix(i, j);
            }
            if (i != matrix.height() - 1 ) os << '\n';
        }
    }

}  // Slae::Matrix

#endif //MATRIX_CSRMATRIX_H
