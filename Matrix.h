#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <thread>
#include <chrono>
#include <future>
template <class T>
class Matrix
{
private:
    std::vector<std::vector<T>> matrix;
    int n, m;
    int max_threads = std::thread::hardware_concurrency();

public:
    Matrix(const std::string &);
    Matrix(int, int);
    Matrix(){};
    ~Matrix(){};
    void print();
    void read();
    Matrix& operator=(Matrix const &);
    void operator+=(Matrix const &);
    void operator-=(Matrix const &);
    Matrix operator*(Matrix const &) const;
    Matrix operator*(int) const;
    Matrix operator/(int) const;
    Matrix operator+(Matrix const &) const;
    Matrix operator-(Matrix const &) const;
    bool operator==(Matrix const &) const;
    bool operator==(int) const;
    Matrix operator!();
    double get_minor(int,int);
    static Matrix create_zero_matrix(int);
    static Matrix create_identity_matrix(int);
    T operator()(int, int);
    Matrix &id(int);
    void transpose();
    double get_det(const Matrix& m = Matrix<T>::create_zero_matrix(0));
    void print_to_file();
    void read_from_file(const std::string &name);

    double parallel_get_det(int);
    Matrix parallel_multiply(int, int) const;
    Matrix parallel_sum(Matrix const &, int) const;
    Matrix parallel_multiply_matrices(Matrix const &, int) const;
    Matrix parallel_inverse_matrix(int);


    Matrix async_sum(Matrix const &, int) const;
    Matrix async_multiply(int, int) const;
    Matrix async_multiply_matrices(Matrix const &, int) const;



};
#endif