#include <iostream>
#include <cstring>

class String {
  private:
    char* source_;
    size_t size_;
    size_t capacity_;

    String(char* source, size_t size, size_t capacity)
        : source_(source),
          size_(size),
          capacity_(capacity) {}

    size_t Find(const String& other, int flag) const {
      for (ssize_t ind = (flag == 1) ? 0 : (size_ - 1); (flag == 1) ? (ind < static_cast<ssize_t>(size_)) : (ind >= 0); ind += flag) {
        int start = (flag == 1) ? 0 : (other.size_ - 1);

        if (source_[ind] == other.source_[start]) {
          ssize_t jnd = 0;

          while (source_[ind + flag * jnd] == other.source_[start + flag * jnd] && ((flag == 1) ? (other.source_[jnd] != '\0') : (ind >= jnd && static_cast<ssize_t>(other.size_) >= 1 + jnd))) {
            ++jnd;
          }

          if (jnd == static_cast<ssize_t>(other.size_)) {
            return ind - start;
          }
        }
      }
      return length();
    }

  public:
    ~String() {
      delete[] source_;
    }

    String() : source_(new char[1]), size_(0), capacity_(1) { source_[0] = '\0'; }

    String(const char* cstring) : size_(strlen(cstring)){
      capacity_ = size_ + 1;
      source_ = new char[capacity_];
      memcpy(source_, cstring, size_ + 1);
    }

    String(size_t num, char ch) {
      size_ = num;
      capacity_ = size_ + 1;
      source_ = new char[capacity_];
      memset(source_, ch, num);
      source_[num] = '\0';
    }

    String(char ch) : String(1, ch) {}

    String(const String& other) : source_(new char[other.capacity_]), size_(other.size_), capacity_(other.capacity_) { memcpy(source_, other.source_, size_ + 1); }

    String& operator=(const String& other) {
      if (this == &other) {
        return *this;
      }

      if (other.size_ > size_) {
        delete[] source_;
        source_ = new char[capacity_ = other.capacity_];
      }

      memcpy(source_, other.source_, (size_ = other.size_) + 1);
      return *this;
    }

    bool operator==(const String& other) const { return strcmp(source_, other.source_) == 0; }

    bool operator<(const String& other) const { return strcmp(source_, other.source_) < 0; }

    bool operator>(const String& other) const { return !(*this < other) && !(*this == other); }

    bool operator<=(const String& other) const { return (*this < other) || (*this == other); }

    bool operator>=(const String& other) const { return (*this > other) || (*this == other); }

    bool operator!=(const String& other) const { return !(*this == other); }

    char& operator[](size_t index) { return source_[index]; }

    const char& operator[](size_t index) const { return source_[index]; }

    String& operator+=(const String& other) {
      if (capacity_ < size_ + other.size_ + 1) {
        char* new_source = new char[capacity_ = size_ + other.size_ + 1];
        memcpy(new_source, source_, size_);
        delete[] source_;
        source_ = new_source;
      }

      memcpy(source_ + size_, other.source_, other.size_ + 1);
      size_ += other.size_;
      return *this;
    }

    String operator+(const String& other) const { return String(*this) += other; }

    void push_back(char ch) {
      if (capacity_ <= size_ + 1) {
        char* new_source = new char[capacity_ = capacity_ * 2];
        memcpy(new_source, source_, size_);
        delete[] source_;
        source_ = new_source;
      }

      source_[size_++] = ch;
      source_[size_] = '\0';
    }

    void pop_back() { source_[(size_ > 0) ? --size_ : 0] = '\0'; }

    size_t find(const String& other) const { return Find(other, 1); }

    size_t rfind(const String& other) const { return Find(other, -1); }

    String substr(size_t start, size_t count) const {
      String str(count, '\0');

      if (start + count > size_) { count = size_ - start; }

      memcpy(str.source_, source_ + start, count);
      return str;
    }

    bool empty() const { return size_ == 0; }

    void clear() {
      source_[0] = '\0';
      size_ = 0;
    }

    void shrink_to_fit() {
      char* new_source = new char[size_ + 1];
      memcpy(new_source, source_, capacity_ = size_ + 1);
      delete[] source_;
      source_ = new_source;
    }

    const char& front() const { return source_[0]; }

    const char& back() const { return source_[size_ - 1]; }

    char& front() { return source_[0]; }

    char& back() { return source_[size_ - 1]; }

    char* data() { return source_; }

    const char* data() const { return source_; }

    size_t size() const { return size_; }

    size_t length() const { return size_; }

    size_t capacity() const { return capacity_ - 1; }
};

String operator+(char ch, const String& string) { return (String(ch) += string); }

bool operator==(const char* cstring, const String& other) { return strcmp(cstring, other.data()) == 0; }

bool operator<(const char* cstring, const String& other) { return strcmp(cstring, other.data()) < 0; }

bool operator>(const char* cstring, const String& other) { return !(cstring == other) && !(cstring < other); }

bool operator<=(const char* cstring, const String& other) { return (cstring == other) || (cstring < other); }

bool operator>=(const char* cstring, const String& other) { return (cstring == other) || (cstring > other); }

bool operator!=(const char* cstring, const String& other) { return !(cstring == other); }

std::ostream& operator<<(std::ostream& fout, const String& string) {return fout << string.data();}

std::istream& operator>>(std::istream& fin, String& string) {
  string.clear();
  char ch;
  fin.get(ch);

  while (!isspace(ch) && !fin.eof()) {
    string.push_back(ch);
    fin.get(ch);
  }

  return fin;
}
