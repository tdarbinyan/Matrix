#pragma once
#include <algorithm>
#include <vector>

template <size_t N, size_t M, typename T = int64_t>
class Matrix {
  std::vector<std::vector<T>> matr_ =
      std::vector<std::vector<T>>(N, std::vector<T>(M, T()));

 public:
  // Constructors
  Matrix() = default;
  Matrix(T element);
  Matrix(Matrix<N, M, T>& matrix);
  Matrix(std::vector<std::vector<T>>& vector);

  // Operators
  template <size_t K, size_t L>
  bool operator==(Matrix<K, L, T> matrix) const;

  template <size_t K>
  Matrix<N, K, T> operator*(const Matrix<M, K, T>& matrix) const;

  Matrix operator*(T element) const;
  Matrix operator+(const Matrix<N, M, T>& matrix) const;
  Matrix operator-(const Matrix<N, M, T>& matrix) const;

  Matrix& operator=(const Matrix<N, M, T>& matrix);
  Matrix& operator+=(const Matrix<N, M, T>& matrix);
  Matrix& operator-=(const Matrix<N, M, T>& matrix);

  T& operator()(size_t first_index, size_t second_index);
  T operator()(size_t first_index, size_t second_index) const;

  // Methods
  Matrix<M, N, T> Transposed();
  T Trace() const;
};

template <size_t N, size_t M, typename T>
Matrix<N, M, T>::Matrix(T element) {
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {
      this->matr_[i][j] = element;
    }
  }
}

template <size_t N, size_t M, typename T>
Matrix<N, M, T>::Matrix(Matrix<N, M, T>& matrix) {
  this->matr_ = matrix.matr_;
}

template <size_t N, size_t M, typename T>
Matrix<N, M, T>::Matrix(std::vector<std::vector<T>>& vector) {
  this->matr_ = vector;
}

template <size_t N, size_t M, typename T>
template <size_t K, size_t L>
bool Matrix<N, M, T>::operator==(Matrix<K, L, T> matrix) const {
  return ((N == K) && (M == L) && (this->matr_ == matrix.matr_));
}

template <size_t N, size_t M, typename T>
Matrix<N, M, T> Matrix<N, M, T>::operator*(T element) const {
  Matrix<N, M, T> resulting_matrix;
  resulting_matrix += *this;
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {
      resulting_matrix(i, j) *= element;
    }
  }
  return resulting_matrix;
}

template <size_t N, size_t M, typename T>
Matrix<N, M, T>& Matrix<N, M, T>::operator=(const Matrix<N, M, T>& matrix) {
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {
      this->matr_[i][j] = matrix.matr_[i][j];
    }
  }
  return *this;
}

template <size_t N, size_t M, typename T>
Matrix<N, M, T> Matrix<N, M, T>::operator+(
    const Matrix<N, M, T>& matrix) const {
  Matrix<N, M, T> resulting_matrix;
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {
      resulting_matrix(i, j) = this->matr_[i][j] + matrix(i, j);
    }
  }
  return resulting_matrix;
}

template <size_t N, size_t M, typename T>
Matrix<N, M, T> Matrix<N, M, T>::operator-(
    const Matrix<N, M, T>& matrix) const {
  Matrix<N, M, T> resulting_matrix;
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {
      resulting_matrix(i, j) = this->matr_[i][j] - matrix(i, j);
    }
  }
  return resulting_matrix;
}

template <size_t N, size_t M, typename T>
Matrix<N, M, T>& Matrix<N, M, T>::operator+=(const Matrix<N, M, T>& matrix) {
  Matrix<N, M, T> temp_matrix = *this;
  temp_matrix = temp_matrix + matrix;
  this->matr_ = temp_matrix.matr_;
  return *this;
}

template <size_t N, size_t M, typename T>
Matrix<N, M, T>& Matrix<N, M, T>::operator-=(const Matrix<N, M, T>& matrix) {
  Matrix<N, M, T> temp_matrix = *this;
  temp_matrix = temp_matrix - matrix;
  this->matr_ = temp_matrix.matr_;
  return *this;
}

template <size_t N, size_t M, typename T>
template <size_t K>
Matrix<N, K, T> Matrix<N, M, T>::operator*(
    const Matrix<M, K, T>& matrix) const {
  Matrix<N, K, T> resulting_matrix;
  for (size_t i = 0; i < N; ++i) {
    for (size_t k = 0; k < K; ++k) {
      T sum;
      for (size_t j = 0; j < M; ++j) {
        (j == 0) ? sum = this->matr_[i][j] * matrix(j, k)
                 : sum += this->matr_[i][j] * matrix(j, k);
      }
      resulting_matrix(i, k) = sum;
    }
  }
  return resulting_matrix;
}

template <size_t N, size_t M, typename T>
T& Matrix<N, M, T>::operator()(size_t first_index, size_t second_index) {
  return this->matr_[first_index][second_index];
}

template <size_t N, size_t M, typename T>
T Matrix<N, M, T>::operator()(size_t first_index, size_t second_index) const {
  return this->matr_[first_index][second_index];
}

template <size_t N, size_t M, typename T>
Matrix<M, N, T> Matrix<N, M, T>::Transposed() {
  Matrix<M, N, T> resulting_matrix;
  for (size_t j = 0; j < M; ++j) {
    for (size_t i = 0; i < N; ++i) {
      resulting_matrix(j, i) = this->matr_[i][j];
    }
  }
  return resulting_matrix;
}

template <size_t N, size_t M, typename T>
T Matrix<N, M, T>::Trace() const {
  if (N == M) {
    T result;
    for (size_t i = 0; i < N; i++) {
      (i == 0) ? result = this->matr_[i][i] : result += this->matr_[i][i];
    }
    return result;
  }
  throw("Can't calculate trace for non square matrix");
}
