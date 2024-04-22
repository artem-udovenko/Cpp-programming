#include <iostream>

template<size_t N>
class StackStorage {
  private:
    char stack_storage[N];
    char* ptr;

  public:
    StackStorage() : ptr(stack_storage) {}

    StackStorage(const StackStorage& other) = delete;

    ~StackStorage() = default;

    char* alloc(unsigned long long n) {
      char* tmp = ptr;
      ptr += n;
      return tmp;
    }

    char* get() {
      return ptr;
    }
};

template<typename T, size_t N>
class StackAllocator {
  private:
    StackStorage<N>* storage;

  public:
    using value_type = T;
    StackStorage<N>* get_storage() const { return storage; }

    explicit StackAllocator(StackStorage<N>& storage) : storage(&storage) {}

    StackAllocator() : storage(nullptr) {}

    template<typename U>
    StackAllocator(const StackAllocator<U, N>& other) : storage(other.get_storage()) {}

    ~StackAllocator() = default;

    template <typename U>
    StackAllocator& operator=(const StackAllocator<U, N>& other) {
      storage = other.get_storage();
      return *this;
    }

    T* allocate(size_t n) {
      if (reinterpret_cast<unsigned long long>(storage->get()) % alignof(T) > 0) {
        storage->alloc(alignof(T) - reinterpret_cast<unsigned long long>(storage->get()) % alignof(T));
      }
      return reinterpret_cast<T*>(storage->alloc(n * sizeof(T)));
    }

    void deallocate(T*, size_t) {}

    template<typename U>
    bool operator==(const StackAllocator<U, N>& other) { return get_storage() == other.get_storage(); }

    template<typename U>
    bool operator!=(const StackAllocator<U, N>& other) { return !(*this == other); }

    template<typename U>
    struct rebind {
      using other = StackAllocator<U, N>;
    };
};