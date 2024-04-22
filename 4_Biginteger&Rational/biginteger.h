#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

const int RAD =  1'000'000'000;
const int DEC = 10;
const int POW = 9;
const int PREC = 15;

class BigInteger;
class Rational;

class BigInteger {
  public:
    enum class Sign : int {
      Negative = -1,
      Neutral = 0,
      Positive = 1
    };

  private:
    Sign sign_ = Sign::Positive;
    std::vector<int> numbers_;

    BigInteger(Sign sign, std::vector<int> numbers) : sign_(sign), numbers_(numbers) {}

    BigInteger Abs(const BigInteger& biginteger) const { return BigInteger((biginteger.sign_ == Sign(0)) ? Sign(0) : Sign(1), biginteger.numbers_); }

    void Increase() {
      for (size_t i = 0; i < numbers_.size(); ++i) {
        if (numbers_[i] < RAD - 1) {
          ++numbers_[i];
          return;
        }
        numbers_[i] = 0;
      }
      numbers_.push_back(1);
    }

    void Decrease() {
      for (size_t i = 0; i < numbers_.size(); ++i) {
        if (numbers_[i] > 0) {
          --numbers_[i];
          if (numbers_.back() == 0 && numbers_.size() > 1) {
            numbers_.pop_back();
          }
          return;
        }
        numbers_[i] = RAD - 1;
      }
    }

    void Reduce(std::vector<int>& vector) const {
      if (vector.empty()) { return; }
      while (vector.back() == 0 && vector.size() > 1) {
        vector.pop_back();
      }
    }

    bool LessOrEqual(const std::vector<int>& vector1, const std::vector<int>& vector2) const {
      if (vector1.size() < vector2.size()) { return true; }
      if (vector1.size() > vector2.size()) { return false; }
      for (ssize_t i = vector1.size() - 1; i >= 0 ; --i) {
        if (vector1[i] < vector2[vector2.size() - i - 1]) {
          return true;
        }
        if (vector1[i] > vector2[vector2.size() - i - 1]) {
          return false;
        }
      }
      return true;
    }

    std::vector<int> Add(const std::vector<int>& vector1, const std::vector<int>& vector2) const {
      std::vector<int> result(vector1.size() + 1, 0);
      size_t index;
      int transfer = 0;
      for (index = 0; index < vector1.size(); ++index) {
        int value = vector1[index] + ((index < vector2.size()) ? vector2[index] : 0) + transfer;
        if (value < RAD) {
          result[index] = value;
          transfer = 0;
        } else {
          result[index] = value - RAD;
          transfer = 1;
        }
      }
      if (transfer == 1) {
        result[result.size() - 1] = 1;
      } else {
        result.pop_back();
      }
      return result;
    }

    std::vector<int> Subtract(const std::vector<int>& vector1, const std::vector<int>& vector2) const {
      std::vector<int> result(vector1.size(), 0);
      size_t index;
      int transfer = 0;
      for (index = 0; index < vector1.size(); ++index) {
        int value = vector1[index] - ((index < vector2.size()) ? vector2[index] : 0) + transfer;
        if (value >= 0) {
          result[index] = value;
          transfer = 0;
        } else {
          result[index] = RAD + value;
          transfer = -1;
        }
      }
      Reduce(result);
      return result;
    }

    std::vector<int> RSubtract(const std::vector<int>& vector1, const std::vector<int>& vector2) const {
      std::vector<int> result(vector1.size(), 0);
      size_t index;
      int transfer = 0;
      for (index = 0; index < vector1.size(); ++index) {
        int value = vector1[vector1.size() - index - 1] - ((index < vector2.size()) ? vector2[index] : 0) + transfer;
        if (value >= 0) {
          result[vector1.size() - index - 1] = value;
          transfer = 0;
        } else {
          result[vector1.size() - index - 1] = RAD + value;
          transfer = -1;
        }
      }
      while (result[0] == 0 && result.size() > 1) {
        result.erase(result.begin());
      }
      return result;
    }

    std::vector<int> Multiply(const std::vector<int>& vector, int number) const {
      std::vector<int> result(vector.size() + 1, 0);
      int transfer = 0;
      for (size_t i = 0; i < vector.size(); ++i) {
        result[i] = (1ll * vector[i] * number + transfer) % RAD;
        transfer = (1ll * vector[i] * number + transfer) / RAD;
      }
      result[vector.size()] = transfer;
      Reduce(result);
      return result;
    }

    std::vector<int> Divide(const std::vector<int>& vector1, const std::vector<int>& vector2) const {
      std::vector<int> buffer;
      std::vector<int> result(1, 0);
      bool flag = true;
      for (ssize_t i = vector1.size() - 1; i >= 0; --i) {
        buffer.push_back(vector1[i]);
        if (buffer[0] == 0 && buffer.size() > 1) {
          buffer.erase(buffer.begin());
        }
        int lnum = 0;
        int rnum = RAD;
        while (lnum != rnum - 1) {
          ((LessOrEqual(Multiply(vector2, (lnum + rnum) / 2), buffer)) ? lnum : rnum) = (lnum + rnum) / 2;
        }
        --rnum;
        if (rnum > 0 && flag) {
          result = std::vector<int>(0);
          flag = false;
        }
        result.push_back(rnum);
        buffer = RSubtract(buffer, Multiply(vector2, rnum));
      }
      std::vector<int> fresult(result.size());
      for (size_t i = 0; i < result.size(); ++i) {
        fresult[i] = result[result.size() - i - 1];
      }
      return fresult;
    }

    BigInteger(std::string string) {
      sign_ = Sign(1);
      if (string[0] == '-') {
        sign_ = Sign(-1);
        string.erase(string.begin());
      }
      for (ssize_t i = string.size() - 1; i >= 0;) {
        std::string temp;
        for (int j = 0; j < POW; ++j) {
          if (i >= 0) {
            temp.push_back(string[i]);
            --i;
          }
        }
        numbers_.push_back(atoi(std::string(temp.rbegin(), temp.rend()).c_str()));
      }
      if (numbers_ == std::vector<int>(1, 0)) {
        sign_ = Sign(0);
      }
    }

    explicit BigInteger(const char* cstring) : BigInteger(std::string(cstring)) {}

  public:
    ~BigInteger() = default;

    BigInteger() : sign_(Sign(0)), numbers_(std::vector<int>(1, 0)) {}

    BigInteger(int number) {
      sign_ = static_cast<Sign>(Signum(number));
      number = abs(number);
      if (number == 0) {
        numbers_ = std::vector<int>(1, 0);
      }
      while (number > 0) {
        numbers_.push_back(number % RAD);
        number /= RAD;
      }
    }

    BigInteger(const BigInteger& other) = default;

    BigInteger& operator=(const BigInteger& other) = default;

    BigInteger operator-() const {
      BigInteger new_int(*this);
      new_int.sign_ = static_cast<Sign>(static_cast<int>(new_int.sign_) * (-1));
      return new_int;
    }

    BigInteger& operator++() {
      (static_cast<int>(sign_) < 0) ? Decrease() : Increase();
      if (numbers_ == std::vector<int>(1, 0)) {
        sign_ = Sign(0);
      } else if (numbers_ == std::vector<int>(1, 1) && sign_ == Sign(0)) {
        sign_ = Sign(1);
      }
      return *this;
    }

    BigInteger operator++(int) {
      const BigInteger temp(*this);
      ++*this;
      return temp;
    }

    BigInteger& operator--() {
      *this = -(++(-*this));
      return *this;
    }

    BigInteger operator--(int) {
      BigInteger temp(*this);
      --*this;
      return temp;
    }

    bool operator==(const BigInteger& other) const {
      if (sign_ != other.sign_) { return false; }
      return numbers_ == other.numbers_;
    }

    bool operator!=(const BigInteger& other) const { return !(*this == other); }

    bool operator<(const BigInteger& other) const {
      if (sign_ != other.sign_) { return sign_ < other.sign_; }
      if (sign_ == Sign(0)) { return false; }
      const int size = numbers_.size();
      const int osize = other.numbers_.size();
      if (size != osize) { return sign_ == Sign(Signum(osize - size));}
      for (ssize_t i = size - 1; i >= 0 ; --i) {
        if (numbers_[i] < other.numbers_[i]) {
          return sign_ == Sign(1);
        }
        if (numbers_[i] > other.numbers_[i]) {
          return sign_ == Sign(-1);
        }
      }
      return false;
    }

    bool operator>(const BigInteger& other) const { return other < *this; }

    bool operator<=(const BigInteger& other) const { return !(*this > other); }

    bool operator>=(const BigInteger& other) const { return !(*this < other); }

    static int Signum(int num) {
      if (num == 0) { return 0; }
      if (num < 0) { return -1; }
      return 1;
    }

    std::string toString() const {
      std::string string;
      for (size_t i = 0; i < numbers_.size() ; ++i) {
        int temp = numbers_[i];
        for (int j = 0; j < POW; ++ j) {
          char ch = temp % DEC + '0';
          temp /= DEC;
          string.push_back(ch);
        }
      }
      while (string.back() == '0' && string.size() > 1) {
        string.pop_back();
      }
      if (sign_ == Sign(-1)) {
        string.push_back('-');
      }
      return {string.rbegin(), string.rend()};
    }

    explicit operator bool() const { return sign_ != Sign(0); }

    explicit operator int() const { return std::stoi(toString()); }

    BigInteger& operator+=(const BigInteger& other) {
      if (sign_ == Sign(0)) {
        return *this = other;
      }
      if (other.sign_ == Sign(0)) {
        return *this;
      }
      if (sign_ == other.sign_) {
        return (Abs(*this) > Abs(other)) ? *this = {sign_, Add(numbers_, other.numbers_)} : *this = {sign_, Add(other.numbers_, numbers_)};
      }
      std::vector<int> result;
      if (Abs(*this) > Abs(other)) {
        result = Subtract(numbers_, other.numbers_);
        return *this = {(result.back() == 0) ? Sign(0) : sign_, result};
      }
      result = Subtract(other.numbers_, numbers_);
      return *this = {(result.back() == 0) ? Sign(0) : other.sign_, result};
    }

    BigInteger& operator-=(const BigInteger& other) { return *this += -other; }

    BigInteger& operator*=(const BigInteger& other) {
      if (numbers_ == std::vector<int>(1, 0) || other.numbers_ == std::vector<int>(1, 0)) { return *this = {0}; }
      BigInteger result(0);
      for (size_t i = 0; i < other.numbers_.size(); ++i) {
        std::vector<int> temp1(i,0);
        std::vector<int> temp2 = Multiply(numbers_, other.numbers_[i]);
        temp1.insert(temp1.end(), temp2.begin(), temp2.end());
        Reduce(temp1);
        result += BigInteger((temp1.back() == 0) ? Sign(0) : Sign(1), temp1);
      }
      if (sign_ == Sign(1) ^ other.sign_ == Sign(1)) {
        return *this = -result;
      }
      return *this = result;
    }

    BigInteger& operator/=(const BigInteger& other) {
      assert(other.sign_ != Sign(0));
      if (sign_ == Sign(0) || Abs(*this) < Abs(other)) { return *this = {0}; }
      BigInteger result(Sign(1), Divide(numbers_, other.numbers_));
      if (sign_ == Sign(1) ^ other.sign_ == Sign(1)) {
        return *this = -result;
      }
      return *this = result;
    }

    BigInteger& operator%=(const BigInteger& other) { return *this -= ((BigInteger(*this) /= other) *= other); }

    friend std::istream& operator>>(std::istream& fin, BigInteger& biginteger);

    friend BigInteger operator"" _bi(const char* number);

    friend std::istream& operator>>(std::istream& fin, Rational& rational);

    friend BigInteger operator"" _bi(unsigned long long number);
};

std::ostream& operator<<(std::ostream& fout, const BigInteger& biginteger) { return fout << biginteger.toString(); }

std::istream& operator>>(std::istream& fin, BigInteger& biginteger) {
  std::string string;
  fin >> string;
  biginteger = BigInteger(string);
  return fin;
}

BigInteger operator"" _bi(const char* number) {
  return BigInteger(number);
}

BigInteger operator"" _bi(unsigned long long number) {
  return {std::to_string(number)};
}

BigInteger operator+(const BigInteger& one, const BigInteger& other) {
  BigInteger result = one;
  return result += other;
}

BigInteger operator-(const BigInteger& one, const BigInteger& other) { return BigInteger(one) -= other; }

BigInteger operator*(const BigInteger& one, const BigInteger& other) {
  BigInteger result(one);
  return result *= other;
}

BigInteger operator/(const BigInteger& one, const BigInteger& other) {
  BigInteger result(one);
  return result /= other;
}

BigInteger operator%(const BigInteger& one, const BigInteger& other) { return BigInteger(one) %= other; }

BigInteger GCD(BigInteger int1, BigInteger int2) {
  while(int1 != 0 && int2 != 0) {
    if (int1 > int2) {
      int1 %= int2;
      continue;
    }
    int2 %= int1;
  }
  return (int1 == 0) ? int2 : int1;
}

class Rational {
    using Sign = BigInteger::Sign;
  private:
    Sign sign_;
    BigInteger numerator_;
    BigInteger denominator_;


    static BigInteger Abs(const BigInteger& biginteger) { return (biginteger < 0) ? -biginteger : biginteger; }

    static Sign SignInt(const BigInteger& biginteger) {
      if (biginteger > 0) { return Sign(1); }
      if (biginteger < 0) { return Sign(-1); }
      return Sign(0);
    }

  public:
    ~Rational() = default;

    Rational() : sign_(Sign(0)), numerator_(BigInteger()), denominator_(BigInteger(1)) {}

    Rational(const Rational& other) = default;

    Rational& operator=(const Rational& other) = default;

    Rational(const BigInteger& biginteger) : sign_(SignInt(biginteger)), numerator_(Abs(biginteger)), denominator_(BigInteger(1)) {}

    Rational(int number) : Rational(BigInteger(number)) {}

    Rational operator-() const {
      Rational result(*this);
      result.sign_ = static_cast<Sign>(static_cast<int>(result.sign_) * (-1));
      return result;
    }

    bool operator==(const Rational& other) const { return (sign_ == other.sign_ && numerator_ == other.numerator_ && denominator_ == other.denominator_); }

    bool operator!=(const Rational& other) const { return !(*this == other); }

    bool operator<(const Rational& other) const {
      if (sign_ < other.sign_) { return true; }
      if (sign_ > other.sign_ || sign_ == Sign(0)) { return false; }
      return (numerator_ * other.denominator_ < denominator_ * other.numerator_) ^ (sign_ == Sign(-1));
    }

    bool operator>(const Rational& other) const { return other < *this; }

    bool operator<=(const Rational& other) const { return !(*this > other); }

    bool operator>=(const Rational& other) const { return !(*this < other); }

    explicit operator bool() const { return sign_ != Sign(0); }

    std::string toString() const {
      std::string string;
      if (sign_ == Sign(-1)) {
        string.push_back('-');
      }
      string += numerator_.toString();
      if (denominator_ != 1) {
        string += '/';
        string += denominator_.toString();
      }
      return string;
    }

    void Normalize(Sign& sign, BigInteger& numerator, BigInteger& denominator) {
      sign = SignInt(numerator);
      numerator = Abs(numerator);
      BigInteger gcd = GCD(numerator, denominator);
      numerator /= gcd;
      denominator /= gcd;
    }

    Rational& operator+=(const Rational& rational) {
      numerator_ = ((sign_ == Sign(-1)) ? (-1) : 1) * rational.denominator_ * numerator_ + ((rational.sign_ == Sign(-1)) ? (-1) : 1) * rational.numerator_ * denominator_;
      denominator_ *= rational.denominator_;
      Normalize(sign_, numerator_, denominator_);
      return *this;
    }

    Rational& operator-=(const Rational& rational) {
      numerator_ = ((sign_ == Sign(-1)) ? (-1) : 1) * rational.denominator_ * numerator_ + ((rational.sign_ == Sign(1)) ? (-1) : 1) * rational.numerator_ * denominator_;
      denominator_ *= rational.denominator_;
      Normalize(sign_, numerator_, denominator_);
      return *this;
    }

    Rational& operator*=(const Rational& rational) {
      numerator_ = ((sign_ == Sign(-1)) ? (-1) : 1) * numerator_ * ((rational.sign_ == Sign(-1)) ? (-1) : 1) * rational.numerator_;
      denominator_ *= rational.denominator_;
      Normalize(sign_, numerator_, denominator_);
      return *this;
    }

    Rational& operator/=(const Rational& rational) {
      numerator_ = ((sign_ == Sign(-1)) ? (-1) : 1) * numerator_ * ((rational.sign_ == Sign(-1)) ? (-1) : 1) * rational.denominator_;
      denominator_ *= rational.numerator_;
      Normalize(sign_, numerator_, denominator_);
      return *this;
    }

    std::string asDecimal(size_t precision = 0) {
      if (precision == 0) { return (numerator_ / denominator_).toString(); }
      BigInteger result_int(numerator_ / denominator_);
      BigInteger result_fract(numerator_ % denominator_);
      for (size_t i = 0; i < precision; ++i) { result_fract *= DEC; }
      result_fract /= denominator_;
      std::string zeros(precision - result_fract.toString().size(), '0');
      std::string result;
      if (sign_ == Sign(-1)) {
        result.push_back('-');
      }
      return result += result_int.toString() + '.' + zeros + result_fract.toString();
    }

    explicit operator double() { return std::stod(asDecimal(PREC)); }
};

Rational operator+(const Rational& rational1, const Rational& rational2) { return Rational(rational1) += rational2; }

Rational operator-(const Rational& rational1, const Rational& rational2) { return Rational(rational1) -= rational2; }

Rational operator*(const Rational& rational1, const Rational& rational2) { return Rational(rational1) *= rational2; }

Rational operator/(const Rational& rational1, const Rational& rational2) { return Rational(rational1) /= rational2; }

std::ostream& operator<<(std::ostream& fout, const Rational& rational) { return fout << rational.toString(); }