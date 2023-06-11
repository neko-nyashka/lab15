#include "Matrix.h"
template <typename T>
void Matrix<T>::print_to_file()
{
    std::ofstream out("output.txt");
    if (!out.is_open())
    {
        std::cerr << "failed to open output file" << std::endl;
        exit(1);
    }
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            out << matrix[i][j] << ' ';
        }
        out << '\n';
    }
    out.close();
}

template <typename T>
void Matrix<T>::read_from_file(const std::string &name)
{
    std::string str;
    std::ifstream in(name);
    if (!in.is_open())
    {
        std::cerr << "failed to open input file" << std::endl;
        exit(1);
    }
    getline(in, str);
    std::istringstream iss(str);
    int n, m;
    iss >> n >> m;
    this->n = n;
    this->m = m;
    matrix.resize(n, std::vector<T>(m));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            T elem;
            iss >> elem;
            matrix[i][j] = elem;
        }
    }
    in.close();
}

template <typename T>
Matrix<T> &Matrix<T>::id(int size)
{
    for (int i = 0; i < size; i++)
    {
        matrix[i][i] = 1;
    }
    return *this;
}

template <typename T>
Matrix<T>::Matrix(int rows, int cols) : n(rows), m(cols)
{
    matrix.resize(n, std::vector<T>(m, 0));
}

template <typename T>
Matrix<T>::Matrix(const std::string &str)
{
    std::istringstream iss(str);
    int n, m;
    iss >> n >> m;
    this->n = n;
    this->m = m;
    matrix.resize(n, std::vector<T>(m));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            T elem;
            iss >> elem;
            matrix[i][j] = elem;
        }
    }
}


template <typename T>
void Matrix<T>::print()
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            std::cout << matrix[i][j] << ' ';
        }
        std::cout << '\n';
    }
}

template <typename T>
Matrix<T> &Matrix<T>::operator=(Matrix const &a)
{
    if (matrix.empty())
    {
        matrix.resize(a.n, std::vector<T>(a.m));
        n = a.n;
        m = a.m;
    }
    if ((n != a.n || m != a.m))
    {
        std::cerr << "matrices of different sizes";
        exit(1);
    }
    else
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                matrix[i][j] = a.matrix[i][j];
            }
        }
        return *this;
    }
}

template <typename T>
void Matrix<T>::operator+=(Matrix const &a)
{
    if (n != a.n || m != a.m)
    {
        std::cerr << "matrices of different sizes";
        exit(1);
    }
    else
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                matrix[i][j] += a.matrix[i][j];
            }
        }
    }
}

template <typename T>
void Matrix<T>::operator-=(Matrix const &a)
{
    if (n != a.n || m != a.m)
    {
        std::cerr << "matrices of different sizes";
        exit(1);
    }
    else
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                matrix[i][j] -= a.matrix[i][j];
            }
        }
    }
}

template <typename T>
Matrix<T> Matrix<T>::operator*(Matrix const &a) const
{
    Matrix<T> result(std::min(n,a.n), std::min(m,a.m));
    if (m != a.n)
    {
        std::cerr << "matrices cannot be multiplied \n the number of columns of the first matrix must equal the number of rows of the second matrix";
        exit(1);
    }
    else
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < a.m; ++j)
            {
                T sum = 0;
                for (int k = 0; k < a.n; ++k)
                {
                    sum += matrix[i][k] * a.matrix[k][j];
                }
                result.matrix[i][j] = sum;
            }
        }
    }
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(int number) const
{
    Matrix<T> a(n, m);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {

            a.matrix[i][j] = number * matrix[i][j];
        }
    }
    return a;
}

template <typename T>
Matrix<T> Matrix<T>::operator/(int number) const
{
    Matrix<T> a(n, m);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {

            a.matrix[i][j] =  matrix[i][j] / number;
        }
    }
    return a;
}



template <typename T>
Matrix<T> Matrix<T>::operator+(Matrix<T> const &b) const
{
    Matrix<T> c(b.n, b.m);
    if (n != b.n || m != b.m)
    {
        std::cerr << "matrices of different sizes";
        exit(1);
    }
    else
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                c.matrix[i][j] = matrix[i][j] + b.matrix[i][j];
            }
        }
    }
    return c;
}



template <typename T>
Matrix<T> Matrix<T>::operator-(Matrix<T> const &b) const
{
    Matrix<T> c(b.n, b.m);
    if (n != b.n || m != b.m)
    {
        std::cerr << "matrices of different sizes";
        exit(1);
    }
    else
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                c.matrix[i][j] = matrix[i][j] - b.matrix[i][j];
            }
        }
    }
    return c;
}

template <typename T>
bool Matrix<T>::operator==(Matrix const &a) const
{
    if (n != a.n || m != a.m)
    {
        return false;
    }
    else
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < a.m; ++j)
            {
                if (matrix[i][j] != a.matrix[i][j])
                {
                    return false;
                }
            }
        }
        return true;
    }
}
template <typename T>
void Matrix<T>::transpose()
{
    if (n != m){
        throw std::invalid_argument("Only square matrices can be transposed.");
    } 
    T tmp;
    for (int i = 0; i < n; ++i)
    {
        for (int j = i; j < m; ++j)
        {   
            tmp = matrix[i][j];
            matrix[i][j] = matrix[j][i];
            matrix[j][i] = tmp;
        }
    }
}

template<typename T>
double Matrix<T>::get_det(const Matrix& m)
{
    if (m.n == 0)
    {
        return get_det(*this);
    }
    if (m.n == 1)
    {
        return m.matrix[0][0];
    }
    if (m.n == 2)
    {
        return m.matrix[0][0] * m.matrix[1][1] - m.matrix[0][1] * m.matrix[1][0];
    }
    else
    {
    double d = 0;
    for (int k = 0; k < m.n; ++k)
    {
        Matrix<T> newm(m.n - 1, m.m - 1); 
        for (int i = 1; i < m.n; ++i) 
        {
            int t = 0;
            for (int j = 0; j < m.m; ++j)
            {
                if (j == k)
                {
                    continue;
                }
                newm.matrix[i - 1][t] =  m.matrix[i][j];
                ++t;
            }
        }
        d += std::pow(-1, k + 2) * m.matrix[0][k] * get_det(newm);
    }
    return d;
    }
}

template<typename T>
double Matrix<T>::get_minor(int r, int c)
{
    Matrix<T> newm(n - 1, m - 1); 
    int t = 0;
    for (int i = 0; i < n; ++i) 
    {
        for (int j = 0; j < m; ++j)
        {
            if (i == r || j == c) 
            {
                continue;
            }
            newm.matrix[t / (m - 1)][t % (m - 1)] = matrix[i][j];
            ++t;                        
        }
    }
    return std::pow(-1, r + c) * get_det(newm);
}

template <typename T>
Matrix<T> Matrix<T>::operator!()
{
    double det = get_det();
    if (!det) throw std::runtime_error("Matrix cannot be reversed.");

    Matrix newm(n, m);
    Matrix<T>::transpose();

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            double minor = get_minor(i, j);
            newm.matrix[i][j] = minor;
        }
    }

    Matrix<T>::transpose();
    return newm / det;
}

template <typename T>
bool Matrix<T>::operator==(int number) const
{
    if (n != m)
    {
        return false;
    }
    if (number != 0 && number != 1)
    {
        std::cerr << "impossible to compare";
        exit(1);
    }
    else
    {
        if (number == 0)
        {
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < m; ++j)
                {
                    if (matrix[i][j] != 0)
                    {
                        return false;
                    }
                }
            }
            return true;
        }
        else
        {
            for (int i = 0; i < n; ++i)
            {
                if (matrix[i][i] != 1)
                {
                    return false;
                }
            }
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < m; ++j)
                {
                    if (i != j)
                    {
                        if (matrix[i][j] != 0)
                        {
                            return false;
                        }
                    }
                }
            }
            return true;
        }
    }
}
template <typename T>
Matrix<T> Matrix<T>::create_zero_matrix(int size)
{
    Matrix mtrx(size, size);
    return mtrx;
}

template <typename T>
Matrix<T> Matrix<T>::create_identity_matrix(int size)
{
    Matrix mtrx(size, size);
    return mtrx.id(size);
}


template <typename T>
Matrix<T> Matrix<T>::parallel_multiply(int h) const
{
    Matrix<T> c(n, m);
    int num_threads = std::thread::hardware_concurrency();
    num_threads = (num_threads > n? n : num_threads);
    std::vector<std::thread> threads(num_threads);
    int block_size = n / num_threads;
    int mod = n % num_threads;
    int start = 0;
    for (int i = 0; i < num_threads; i++)
    {
        int end = start + block_size + (i == 0 ? mod : 0);
        threads[i] = std::thread([this, h, &c, start, end]()
        {
            for (int i = start; i < end; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    c.matrix[i][j] += matrix[i][j] * h;
                }
            }
        });
        start = end;
    }
    for (int i = 0; i < num_threads; i++)
    {
        threads[i].join();
    }
    return c;
}

template <typename T>
Matrix<T> Matrix<T>::parallel_sum(Matrix<T> const &b) const
{
     Matrix<T> c(b.n, b.m);
    if (n != b.n || m != b.m)
    {
        std::cerr << "matrices of different sizes";
        exit(1);
    }
    else
    {
        int num_threads = std::thread::hardware_concurrency();
        num_threads = (num_threads > n? n : num_threads);
        std::vector<std::thread> threads(num_threads);
        int block_size = n / num_threads;
        int mod = n % num_threads;
        int start = 0;
        std::mutex mtx;
        for (int i = 0; i < num_threads; i++)
        {
            int end = start + block_size + (i == 0 ? mod : 0);
            threads[i] = std::thread([this, &b, &c, start, end, &mtx]()
            {
                for (int i = start; i < end; i++)
                {
                    std::lock_guard<std::mutex> lock(mtx);
                    for (int j = 0; j < m; j++)
                    {
                        c.matrix[i][j] = matrix[i][j] + b.matrix[i][j];

                    }
                }
            });
            start = end;
        }
        for (int i = 0; i < num_threads; i++)
        {
            threads[i].join();
        }
        return c;
    }



}



/// Algorithms with threads ///




template <typename T>
double Matrix<T>::parallel_get_det() 
{
    if (this -> n == 1)
    {
        return this-> matrix[0][0];
    }
    if (this -> n == 2)
    {
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    }
    else
    {
        double d = 0;
        std::vector<std::thread> threads;
        std::vector<double> determinants(this->n, 0);

        for (int k = 0; k < this->n; ++k)
        {
            threads.emplace_back([&, k]() 
            {
                Matrix newm(n - 1, m - 1);
                for (int i = 1; i < this->n; ++i)
                {
                    int t = 0;
                    for (int j = 0; j < this->m; ++j)
                    {
                        if (j == k) continue;
                        newm.matrix[i - 1][t] = matrix[i][j];
                        ++t;
                    }
                }
            determinants[k] = std::pow(-1, k + 2) * matrix[0][k] * get_det(newm);
            });
    }

    for (auto& thread : threads)
    {
        thread.join();
    }

    for (double determinant : determinants)
    {
        d += determinant;
    }

    return d;
    }
}


template <typename T>
Matrix<T> Matrix<T>::parallel_inverse_matrix()
{
    double det = parallel_get_det();
    if (!det) throw std::runtime_error("Matrix cannot be reversed.");

    Matrix newm(n, m);
    Matrix<T>::transpose();

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            double minor = get_minor(i, j);
            newm.matrix[i][j] = minor;
        }
    }

    Matrix<T>::transpose();
    return newm / det;
}


template <typename T>
Matrix<T> Matrix<T>::parallel_multiply_matrices(Matrix<T> const &a) const
{
    Matrix<T> result(std::min(a.n, n) , std::min(a.m,m));
     if (m != a.n)
    {
        std::cerr << "matrices cannot be multiplied \n the number of columns of the first matrix must equal the number of rows of the second matrix";
        exit(1);
    }
    else
    {
        int num_threads = std::thread::hardware_concurrency();
        num_threads = (num_threads > n? n : num_threads);
        std::vector<std::thread> threads(num_threads);
        int block_size = n / num_threads;
        int mod = n % num_threads;
        int start = 0;
        std::mutex mtx;
        for (int i = 0; i < num_threads; i++)
        {
            int end = start + block_size + (i == 0 ? mod : 0);
            threads[i] = std::thread([this, &a, &result, start, end]()
            {
                for (int i = 0; i < n; ++i)
                {
                    for (int j = 0; j < a.m; ++j)
                    {
                    T sum = 0;
                    for (int k = 0; k < a.n; ++k)
                    {
                        sum += matrix[i][k] * a.matrix[k][j];
                    }
                    result.matrix[i][j] = sum;
                    }
                }
                
            });
            start = end;
        }
        for (int i = 0; i < num_threads; i++)
        {
            threads[i].join();
        }
        return result;
    }



}




/// Algorithms with async ///




template <typename T>
Matrix<T> Matrix<T>::async_multiply(int h) const
{
    Matrix<T> result(n,m);
    int blocks = std::thread::hardware_concurrency();
    blocks = ( blocks > n ? n : blocks );
    int rows_in_block = n / blocks;
    std::vector<std::future<void>> futures(blocks);
    for (int i = 0; i < blocks; ++i) {
        futures[i] = std::async(std::launch::async, [this, h, rows_in_block, blocks, &result, i]() {
            int start_row = i * rows_in_block;
            int end_row = (i == blocks - 1) ? n : (i + 1) * rows_in_block;
            for(int row = start_row; row < end_row; row++)
            {
                for(int col = 0; col < m; col++)
                {
                    result.matrix[row][col] = matrix[row][col] * h;
                }
            }
        });
    }
    for(auto& f : futures)
    {
        f.wait();
    }
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::async_sum(Matrix<T> b) const
{
    Matrix<T> result(n,m);
    int blocks = std::thread::hardware_concurrency();
    blocks = ( blocks > n ? n : blocks );
    int rows_in_block = n / blocks;
    std::vector<std::future<void>> futures(blocks);
    for (int i = 0; i < blocks; ++i) {
        futures[i] = std::async(std::launch::async, [this, &b , rows_in_block, blocks, &result, i]() {
            int start_row = i * rows_in_block;
            int end_row = (i == blocks - 1) ? n : (i + 1) * rows_in_block;
            for(int row = start_row; row < end_row; row++)
            {
                for(int col = 0; col < m; col++)
                {
                    result.matrix[row][col] = matrix[row][col] + b.matrix[row][col];
                }
            }
        });
    }
    for(auto& f : futures)
    {
        f.wait();
    }
    return result;
}


int main()
{
    Matrix<double> a, b;
    // a=Matrix<int>::create_zero_matrix(5);
    // a.read_from_file("input.txt");
    // b.read_from_file("input.txt");
    // a.print_to_file();
    // b.print_to_file();
    int h=5789;
    std::vector<std::string> v = {"2.txt","2.txt","3.txt","4.txt","5.txt"};
    a.read_from_file(v[0]);
    // c = b.parallel_inverse_matrix();
    for(int file = 0; file <= 4; ++file)
    {
        a.read_from_file(v[file]);
        // Matrix<double> b = a;
        auto start_time = std::chrono::high_resolution_clock::now();
        a.parallel_inverse_matrix();
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        std::cout << duration.count()<< std::endl;
        wait(&h);

    }
    return 0;
}
