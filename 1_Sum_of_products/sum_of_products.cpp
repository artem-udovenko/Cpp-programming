#include <iostream>

bool Contain(int index, const int* indexes, int depth) {
  for (int i = 0; i < depth; ++i) {
    if (indexes[i] == index) {
      return true;
    }
  }
  return false;
}

long long ResultRecursiveCount(int depth, const int kMaxDepth, int* indexes, int** matrix, int* lens) {
  if (depth == kMaxDepth) {
    return 1;
  }
  long long sum = 0;
  for (int i = 0; i < lens[depth]; ++i) {
    if (!Contain(i, indexes, depth)) {
      indexes[depth] = i;
      sum += matrix[depth][i] *
             ResultRecursiveCount(depth + 1, kMaxDepth, indexes, matrix, lens);
    }
  }
  return sum;
}

int main(int argc, char *argv[]) {
  const int kMaxDepth = argc - 1;
  int* lens = new int[kMaxDepth];
  int** matrix = new int*[kMaxDepth];
  for (int i = 0; i < kMaxDepth; ++i) {
    lens[i] = atoi(argv[i + 1]);
    matrix[i] = new int[lens[i]];
  }
  for (int i = 0; i < kMaxDepth; ++i) {
    for (int j = 0; j < lens[i]; ++j) {
      std::cin >> matrix[i][j];
    }
  }
  int* indexes = new int[kMaxDepth];
  std::cout << ResultRecursiveCount(0, kMaxDepth, indexes, matrix, lens);
  delete[] lens;
  for (int i = 0; i < kMaxDepth; ++i) {
    delete[] matrix[i];
  }
  delete[] matrix;
  delete[] indexes;
  return 0;
}
