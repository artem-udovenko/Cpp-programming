#include <iostream>
#include <string.h>

void Input(char*& string) {
  string = new char[1];
  string[0] = 0;
  char ch;
  int size = 0;
  int capacity = 1;
  std::cin.get(ch);

  while (ch == ' ' || ch == '\n') {
    std::cin.get(ch);
  }

  while (ch != ' ' && ch != '\n') {
    string[size] = ch;
    ++size;

    if (size == capacity) {
      char* newstring = new char[capacity *= 2];
      memcpy(newstring, string, size);
      newstring[size] = 0;
      delete[] string;
      string = newstring;
    }

    string[size] = 0;
    std::cin.get(ch);
  }
}

void Push(char**& stack, int& size, int& capacity) {
  Input(stack[size]);
  ++size;

  if (size == capacity) {
    char** newstack = new char*[capacity *= 2];
    memcpy(newstack, stack, size * sizeof(char*));
    stack = newstack;
  }

  std::cout << "ok\n";
}

void Pop(char**& stack, int& size) {
  if (size == 0) {
    std::cout << "error\n";
    return;
  }

  std::cout << stack[--size] << '\n';
}

void Back(char** stack, int size) {
  if (size == 0) {
    std::cout << "error\n";
    return;
  }

  std::cout << stack[size - 1] << '\n';
}

void Size(int size) { std::cout << size << '\n'; }

void Clear(char**& stack, int& size, int& capacity) {
  for (int i = 0; i < size; ++i) {
    delete[] stack[i];
  }

  delete[] stack;
  size = 0;
  stack = new char*[1];
  capacity = 1;
  std::cout << "ok\n";
}

int main() {
  int capacity = 1;
  char** stack = new char*[capacity];
  int size = 0;
  char* command;
  Input(command);

  while (strcmp(command, "exit")) {
    if (strcmp(command, "push") == 0) {
      Push(stack, size, capacity);
    } else if (strcmp(command, "pop") == 0) {
      Pop(stack, size);
    } else if (strcmp(command, "back") == 0) {
      Back(stack, size);
    } else if (strcmp(command, "size") == 0) {
      Size(size);
    } else if (strcmp(command, "clear") == 0) {
      Clear(stack, size, capacity);
    }

    Input(command);
  }
  std::cout << "bye";
  delete[] command;

  for (int i = 0; i < size; ++i) {
    delete[] stack[i];
  }

  delete[] stack;
  return 0;
}