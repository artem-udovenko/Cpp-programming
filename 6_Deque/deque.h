#include <iostream>
#include <vector>

template <typename T>
class Deque {
  private:
    static const int SIZE = 16;
    int capacity;
    T** buffer;
    T** array;
    int HEAD;
    int TAIL;
    int head;
    int tail;

    T** allocate(size_t quantity) {
      buffer = new T*[quantity];
      for (size_t i = 0; i < quantity; ++i) {
        try {
          T* arr = reinterpret_cast<T*>(new char[SIZE * sizeof(T)]);
          buffer[i] = arr;
        } catch(...) {
          for (size_t k = 0; k < i; ++k) {
            delete[] reinterpret_cast<char*>(buffer[k]);
          }
          delete[] buffer;
          throw;
        }
      }
      return buffer;
    }

    void deallocate() {
      for (int i = 0; i < capacity; ++i) {
        delete[] reinterpret_cast<char*>(array[i]);
      }
      delete[] array;
    }

    void front_expansion() {
      size_t delta = TAIL - HEAD + 1;
      T** new_array = allocate(capacity + delta);
      for (int i = 0; i < capacity; ++i) {
        new_array[i + delta] = array[i];
      }
      delete[] array;
      array = new_array;
      capacity += delta;
      HEAD += delta;
      TAIL += delta;
    }

    void back_expansion() {
      size_t delta = TAIL - HEAD + 1;
      T** new_array = allocate(capacity + delta);
      for (int i = 0; i < capacity; ++i) {
        new_array[i] = array[i];
      }
      delete[] array;
      array = new_array;
      capacity += delta;
    }

  public:
    template<typename U = T, const int REV = 1> class Iterator;
    void show() const {
      for (size_t i = 0; i < capacity; ++i) {
        std::cout << i << ": ";
        for (size_t j = 0; j < SIZE; ++j) {
          if (!(HEAD > i || i > TAIL || (i == HEAD && j < head) || (i == TAIL && j >= tail))) {
            std::cout << array[i][j] << ' ';
          }
          else {
            std::cout << "}{ ";
          }
        }
        std::cout << '\n';
      }
      std::cout << "HEAD: " << HEAD << " head: " << head << " TAIL: " << TAIL << " tail: " << tail << '\n';
    }

    Deque() {
      T* arr = reinterpret_cast<T*>(new char[SIZE * sizeof(T)]);
      array = new T*[1];
      array[0] = arr;
      capacity = 1;
      HEAD = TAIL = 0;
      head = tail = SIZE / 2;
    }

    Deque(const Deque& other) {
      array = allocate(other.capacity);
      capacity = other.capacity;
      HEAD = other.HEAD;
      TAIL = other.TAIL;
      head = other.head;
      tail = other.tail;
      for (int i = HEAD; i <= TAIL; ++i) {
        for (int j = (i == HEAD ? head : 0); j < (i == TAIL ? tail : SIZE); ++j) {
          try {
            new(array[i] + j) T(other.array[i][j]);
          } catch(...) {
            for (size_t k = HEAD; k <= i; ++k) {
              for (int p = (k == HEAD ? head : 0); p < (k == i ? j : SIZE); ++p) {
                (array[k] + p)->~T();
              }
            }
            deallocate();
            throw;
          }
        }
      }
    }

    Deque(int size, const T& val) {
      array = allocate(capacity = size / SIZE + 1);
      HEAD = head = 0;
      TAIL = size / SIZE;
      tail = size % SIZE;
      for (int i = 0; i <= TAIL; ++i) {
        for (int j = 0; j < (i == TAIL ? tail : SIZE); ++j) {
          try {
            new(array[i] + j) T(val);
          } catch(...) {
            for (int k = 0; k <= i; ++k) {
              for (int p = 0; p < (k == i ? j : SIZE); ++p) {
                (array[k] + p)->~T();
              }
            }
            deallocate();
            throw;
          }
        }
      }
    }

    explicit Deque(int size) : Deque(size, T()) {}

    ~Deque() {
      for (int i = HEAD; i <= TAIL; ++i) {
        for (int j = (i == HEAD ? head : 0); j < (i == TAIL ? tail : SIZE); ++j) {
          (array[i] + j)->~T();
        }
      }
      deallocate();
    }

    Deque<T>& operator=(const Deque<T> other) {
      if (array == other.array) return *this;
      for (int i = HEAD; i <= TAIL; ++i) {
        for (int j = (i == HEAD ? head : 0); j < (i == TAIL ? tail : SIZE); ++j) {
          (array[i] + j)->~T();
        }
      }
      deallocate();
      array = allocate(other.capacity);
      capacity = other.capacity;
      HEAD = other.HEAD;
      TAIL = other.TAIL;
      head = other.head;
      tail = other.tail;
      for (int i = HEAD; i <= TAIL; ++i) {
        for (int j = (i == HEAD ? head : 0); j < (i == TAIL ? tail : SIZE); ++j) {
          try {
            new(array[i] + j) T(other.array[i][j]);
          } catch(...) {
            for (size_t k = HEAD; k <= i; ++k) {
              for (int p = (k == HEAD ? head : 0); p < (k == i ? j : SIZE); ++p) {
                (array[k] + p)->~T();
              }
            }
            deallocate();
            throw;
          }
        }
      }
      return *this;
    }

    size_t size() const { return SIZE * (TAIL - HEAD) + tail - head; }

    T& operator[](size_t index) { return array[HEAD + index / SIZE][head + index % SIZE]; }

    const T& operator[](size_t index) const { return array[HEAD + index / SIZE][head + index % SIZE]; }

    T& at(size_t index) {
      if (index >= size()) throw std::out_of_range{"You have submitted some cringe as index"};
      return array[HEAD + index / SIZE][head + index % SIZE];
    }

    const T& at(size_t index) const {
      if (index >= size()) throw std::out_of_range{"You have submitted some cringe as index"};
      return array[HEAD + index / SIZE][head + index % SIZE];
    }

    void push_front(T val) {
      if (HEAD == 0 && head == 0) front_expansion();
      if (head == 0) {
        --HEAD;
        head = SIZE - 1;
      } else {
        --head;
      }
      new(array[HEAD] + head) T(val);
    }

    void push_back(T val) {
      if (TAIL + 1 == capacity && tail + 1 == SIZE) back_expansion();
      new(array[TAIL] + tail) T(val);
      if (tail + 1 == SIZE) {
        ++TAIL;
        tail = 0;
      } else {
        ++tail;
      }
    }

    void pop_front() {
      if (HEAD == TAIL && head == tail) return;
      (array[HEAD] + head)->~T();
      if (head + 1 == SIZE) {
        ++HEAD;
        head = 0;
      } else {
        ++head;
      }
    }

    void pop_back() {
      if (HEAD == TAIL && head == tail) return;
      if (tail == 0) {
        --TAIL;
        tail = SIZE - 1;
      } else {
        --tail;
      }
      (array[TAIL] + tail)->~T();
    }

    template<typename U, const int REV>
    class Iterator {
        friend class Deque;
      protected:
        int POS;
        int pos;
        U** arr;

        Iterator(int POS, int pos, U** arr) : POS(POS), pos(pos), arr(arr) {}

      public:
        void show() {
          std::cout << "POS: " << POS << " pos: " << pos << '\n';
        }

        using value_type = U;
        using difference_type = int;
        using iterator_category = std::random_access_iterator_tag;
        using pointer = value_type*;
        using reference = value_type&;

        Iterator(const Iterator<U, REV>& other) : POS(other.POS), pos(other.pos), arr(other.arr) {}

        operator Iterator<const U, REV>() {
          return Iterator<const U, REV>(POS, pos, arr);
        }

        Iterator& operator=(Iterator<U, REV> other) {
          POS = other.POS;
          pos = other.pos;
          arr = other.arr;
          return *this;
        }

        Iterator& operator+=(int delta) {
          delta *= REV;
          POS += delta / SIZE;
          pos += delta % SIZE;
          if (delta > 0 && pos >= SIZE) {
            ++POS;
            pos -= SIZE;
          } else if (delta < 0 && pos < 0) {
            --POS;
            pos += SIZE;
          }
          return *this;
        }

        Iterator& operator-=(int delta) { return (*this += -delta); }

        Iterator& operator++() { return (*this += 1); }

        Iterator operator++(int) {
          Iterator temp(*this);
          ++*this;
          return temp;
        }

        Iterator& operator--() { return (*this -= 1); }

        Iterator operator--(int) {
          Iterator temp(*this);
          --*this;
          return temp;
        }

        Iterator operator+(int delta) const { return (Iterator(*this) += delta); }

        Iterator operator-(int delta) const { return (Iterator(*this) -= delta); }

        bool operator==(const Iterator& other) const { return POS == other.POS && pos == other.pos && arr == other.arr; }

        bool operator!=(const Iterator& other) const { return !(*this == other); }

        bool operator<(const Iterator& other) const {
          if (POS != other.POS) return (POS < other.POS && REV == 1) || (other.POS < POS && REV == -1);
          return (pos < other.pos && REV == 1) || (other.pos < pos && REV == -1);
        }

        bool operator>(const Iterator& other) const { return !(*this < other || *this == other); }

        bool operator<=(const Iterator& other) const { return !(*this > other); }

        bool operator>=(const Iterator& other) const { return !(*this < other); }

        int operator-(const Iterator& other) const { return (SIZE * (POS - other.POS) + pos - other.pos) * REV; }

        U& operator*() { return *(arr[POS] + pos); }

        U* operator->() { return arr[POS] + pos; }

        const U& operator*() const { return *(arr[POS] + pos); }

        const U* operator->() const { return arr[POS] + pos; }
    };

    typedef Iterator<T, 1> iterator;
    typedef Iterator<const T, 1> const_iterator;
    typedef Iterator<T, -1> reverse_iterator;
    typedef Iterator<const T, -1> const_reverse_iterator;

    Deque<T>::iterator begin() { return iterator(HEAD, head, array); }

    Deque<T>::iterator end() { return iterator(TAIL, tail, array); }

    Deque<T>::const_iterator begin() const { return const_iterator(HEAD, head, const_cast<const T**>(array)); }

    Deque<T>::const_iterator end() const { return const_iterator(TAIL, tail, const_cast<const T**>(array)); }

    Deque<T>::reverse_iterator rbegin() {
      if (tail > 0) return reverse_iterator(TAIL, tail - 1, array);
      return reverse_iterator(TAIL - 1, SIZE - 1, array);
    }

    Deque<T>::reverse_iterator rend() {
      if (head > 0) return reverse_iterator(HEAD, head - 1, array);
      return reverse_iterator(HEAD - 1, SIZE - 1, array);
    }

    Deque<T>::const_reverse_iterator rbegin() const {
      if (tail > 0) return const_reverse_iterator(TAIL, tail - 1, const_cast<const T**>(array));
      return const_reverse_iterator(TAIL - 1, SIZE - 1, const_cast<const T**>(array));
    }

    Deque<T>::const_reverse_iterator rend() const {
      if (head > 0) return const_reverse_iterator(HEAD, head - 1, const_cast<const T**>(array));
      return const_reverse_iterator(HEAD - 1, SIZE - 1, const_cast<const T**>(array));
    }

    Deque<T>::const_iterator cbegin() const { return const_iterator(HEAD, head, const_cast<const T**>(array)); }

    Deque<T>::const_iterator cend() const { return const_iterator(TAIL, tail, const_cast<const T**>(array)); }

    Deque<T>::const_reverse_iterator crbegin() const {
      if (tail > 0) return const_reverse_iterator(TAIL, tail - 1, const_cast<const T**>(array));
      return const_reverse_iterator(TAIL - 1, SIZE - 1, const_cast<const T**>(array));
    }

    Deque<T>::const_reverse_iterator crend() const {
      if (head > 0) return const_reverse_iterator(HEAD, head - 1, const_cast<const T**>(array));
      return const_reverse_iterator(HEAD - 1, SIZE - 1, const_cast<const T**>(array));
    }

    template<typename U, const int REV>
    void insert(const Iterator<U, REV>& iter0, const T& val) {
      this->push_back(val);
      T cur = val;
      T temp = val;
      for (iterator iter(iter0.POS, iter0.pos, array); iter < this->end(); ++iter) {
        temp = *iter;
        iter->~T();
        new(&(*iter)) T(cur);
        cur = temp;
      }
    }

    template<typename U, const int REV>
    void erase(const Iterator<U, REV>& iter0) {
      for (iterator iter(iter0.POS, iter0.pos, array); iter < this->end() - 1; ++iter) {
        iter->~T();
        new(&(*iter)) T(*(iter + 1));
      }
      this->pop_back();
    }
};
