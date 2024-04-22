#include <cmath>
#include <iostream>
#include <vector>

constexpr double delta = 1e-6;
constexpr double DBL_MAX = 1.7976931348623158e+308;

bool approx_equal(double  x1, double  x2) {
  return std::abs(x1 - x2) < delta;
}

bool approx_less(double  x1, double  x2) {
  return x1 < x2 && !(approx_equal(x1, x2));
}

bool approx_loe(double  x1, double  x2) {
  return x1 < x2 || (approx_equal(x1, x2));
}

class Line;

struct Vector;

struct Point {
  double x;
  double y;

  Point() : x(0), y(0) {}

  Point(double  x, double y) : x(x), y(y) {}

  void rotate(const Point& center, double angle);

  void reflect(const Point& center);

  void reflect(const Line& axis);

  void scale(const Point& center, double coefficient);
};

struct Vector : Point {
  Vector(double x0, double y0) {
    x = x0;
    y = y0;
  }

  Vector(Point p1, Point p2) {
    x = (p2.x - p1.x);
    y = (p2.y - p1.y);
  }

  Vector operator+(Vector other) { return {x + other.x, y + other.y}; }

  Vector operator-(Vector other) { return {x - other.x, y - other.y}; }

  double operator*(Vector other) { return x * other.x + y * other.y; }

  double operator^(Vector other) { return x * other.y - other.x * y; }

  bool operator||(Vector other) { return approx_equal(x * other.y - y * other.x, 0); }

  double Modulus() { return sqrt(x * x + y * y); }

  Vector Orthogonal() { return {y, -x}; }
};

Vector operator*(double k, Vector vector) { return {vector.x * k, vector.y * k}; }

Point operator>>(Point point, Vector vector) { return {point.x + vector.x, point.y + vector.y}; }

bool operator==(Point p1, Point p2) { return approx_equal(p1.x, p2.x) && approx_equal(p1.y, p2.y); }

bool operator!=(Point p1, Point p2) { return !(p1 == p2); }

class Line {
  private:
    Point reference_point_;
    Vector guiding_vector_;

  public:
    Line(Point reference_point, Vector guiding_vector) : reference_point_(reference_point), guiding_vector_(guiding_vector) {}

    Line(Point p1, Point p2) : Line(p1, Vector(p1, p2)) {}

    Line(double coefficient, double shift) : reference_point_(Point(0, shift)), guiding_vector_(Vector(1, coefficient)) {}

    Line(Point reference_point, double coefficient) : reference_point_(reference_point), guiding_vector_(Vector(1, coefficient)) {}

    bool operator==(Line other) { return (Vector(reference_point_, other.reference_point_) || guiding_vector_) && (guiding_vector_ || other.guiding_vector_); }

    bool operator!=(Line other) { return !(*this == other); }

    Point getReferencePoint() const { return reference_point_; }

    Vector getGuidingVector() const { return guiding_vector_; }

    static Point Intersection(Line line1, Line line2) {
      if (line1.guiding_vector_ || line2.guiding_vector_) { return {DBL_MAX, DBL_MAX}; }
      double delta1 = line1.guiding_vector_.x * (-line2.guiding_vector_.y) + line1.guiding_vector_.y * line2.guiding_vector_.x;
      double delta2 = (line2.reference_point_.x - line1.reference_point_.x) * (-line2.guiding_vector_.y) + (line2.reference_point_.y - line1.reference_point_.y) * line2.guiding_vector_.x;
      return {line1.reference_point_.x + line1.guiding_vector_.x * delta2 / delta1, line1.reference_point_.y + line1.guiding_vector_.y * delta2 / delta1};
    }
};

void Point::rotate(const Point& center, double angle) {
  Vector rad(center, *this);
  angle *= M_PI / 180;
  rad = {rad.x * cos(angle) + -rad.y * sin(angle), rad.x * sin(angle) + rad.y * cos(angle)};
  rad = {round(rad.x / delta / delta) * delta * delta, round(rad.y / delta / delta) * delta * delta};
  *this = (center >> rad);
}

void Point::reflect(const Point& center) {
  *this = (center >> Vector(*this, center));
}

void Point::reflect(const Line& axis) {
  Vector rad(axis.getReferencePoint(), *this);
  Vector vector = axis.getGuidingVector();
  rad = rad + 2 * ((rad * vector / vector.Modulus() / vector.Modulus()) * vector - rad);
  *this = (axis.getReferencePoint() >> rad);
}

void Point::scale(const Point& center, double coefficient) {
  *this = (center >> coefficient * Vector(center, *this));
}

class Shape {
  public:
    virtual double perimeter() const = 0;

    virtual double area() const = 0;

    virtual bool operator==(const Shape& another) const = 0;

    bool operator!=(const Shape& another) const {
      return !(*this == another);
    }

    bool isCongruentTo(const Shape& another) const {
      return isSimilarTo(another) && approx_equal(area(), another.area());
    }

    virtual bool isSimilarTo(const Shape& another) const = 0;

    virtual bool containsPoint(const Point& point) const = 0;

    virtual void rotate(const Point& center, double angle) = 0;

    virtual void reflect(const Point& center) = 0;

    virtual void reflect(const Line& axis) = 0;

    virtual void scale(const Point& center, double coefficient) = 0;

    virtual ~Shape() = default;
};

class Polygon : public Shape {
  protected:
    std::vector<Point> vertices_;

    Polygon() = default;

  public:
    Polygon(std::vector<Point> vertices) : vertices_(vertices) {}

    Polygon(std::initializer_list<Point> vertices) : vertices_(vertices) {}

    template<typename... T>
    Polygon(T... points) : vertices_({points...}) {};

    size_t verticesCount() const { return vertices_.size(); }

    const std::vector<Point>& getVertices() const { return vertices_; }

    bool isConvex() const {
      size_t size = vertices_.size();
      size_t index = size - 1;
      Vector v1(vertices_[index % size], vertices_[(index + 1) % size]);
      Vector v2(vertices_[(index + 1) % size], vertices_[(index + 2) % size]);
      double flag = v1 ^ v2;
      for (index = 0; index < size - 1; ++index) {
        v1 = Vector(vertices_[index % size], vertices_[(index + 1) % size]);
        v2 = Vector(vertices_[(index + 1) % size], vertices_[(index + 2) % size]);
        if ((v1 ^ v2) * flag < 0) {
          return false;
        }
      }
      return true;
    }

    double perimeter() const override {
      double perimeter = 0;
      size_t size = vertices_.size();
      for (size_t i = 0; i < size; ++i) {
        perimeter += Vector(vertices_[i % size], vertices_[(i + 1) % size]).Modulus();
      }
      return perimeter;
    }

    double area() const override {
      double area = 0;
      size_t size = vertices_.size();
      for (size_t i = 1; i < size - 1; ++i) {
        area += Vector(vertices_[0], vertices_[i]) ^ Vector(vertices_[0], vertices_[i + 1]);
      }
      return std::abs(area) / 2;
    }

    bool operator==(const Shape& another) const override {
      auto polygon = dynamic_cast<const Polygon*>(&another);
      if (polygon == nullptr) { return false; }
      if (polygon->verticesCount() != verticesCount()) { return false; }
      const size_t number = verticesCount();
      for (size_t i = 0; i < number; ++i) {
        if (vertices_[i] == polygon->vertices_[0]) {
          if (vertices_[(i + 1) % number] == polygon->vertices_[1]) {
            for (size_t j = 0; j < number; ++j) {
              if (vertices_[(i + j) % number] != polygon->vertices_[j]) {
                return false;
              }
            }
            return true;
          }
          if (vertices_[(i + number - 1) % number] == polygon->vertices_[1]) {
            for (size_t j = 0; j < number; ++j) {
              if (vertices_[(i + number - j) % number] != polygon->vertices_[j]) {
                return false;
              }
            }
            return true;

          }
          return false;
        }
      }
      return false;
    }

    std::vector<double> angles_cos() const {
      const int number = verticesCount();
      std::vector<double> answer;
      for (int i = 0; i < number; ++i) {
        Vector v1(vertices_[i], vertices_[(i + 1) %  number]);
        Vector v2(vertices_[i], vertices_[(i - 1 + number) %  number]);
        answer.push_back(v1 * v2 / v1.Modulus() / v2.Modulus());
      }
      return answer;
    }

    std::vector<Vector> sides() const {
      const int number = verticesCount();
      std::vector<Vector> answer;
      for (int i = 0; i < number; ++i) {
        answer.emplace_back(vertices_[i], vertices_[(i + 1) %  number]);
      }
      return answer;
    }

    bool isSimilarTo(const Shape& another) const override {
      auto polygon = dynamic_cast<const Polygon*>(&another);
      if (polygon == nullptr) { return false; }
      if (polygon->verticesCount() != verticesCount()) { return false; }
      const int number = verticesCount();
      for (int i = 0; i < number; ++i) {
        double coef = polygon->sides()[i].Modulus() / sides()[0].Modulus();
        bool flag = true;
        for (int j = 0; j < number; ++j) {
          if (!approx_equal(polygon->sides()[(i + j) % number].Modulus() / sides()[j].Modulus(), coef) || !approx_equal(polygon->angles_cos()[(i + j) % number], angles_cos()[j])) {
            flag = false;
            break;
          }
        }
        if (flag) { return true; }
        coef = polygon->sides()[i].Modulus() / sides()[number - 1].Modulus();
        flag = true;
        for (int j = number; j > 0; --j) {
          if (!approx_equal(polygon->sides()[(i + number - j) % number].Modulus() / sides()[(j - 1) % number].Modulus(), coef) || !approx_equal(polygon->angles_cos()[(i + number - j) % number], angles_cos()[j % number])) {
            flag = false;
            break;
          }
        }
        if (flag) { return true; }
      }
      return false;
    }

    bool containsPoint(const Point& point) const override {
      for (size_t i = 0; i < verticesCount(); ++i) {
        Vector vector(vertices_[i], point);
        if ((vector || sides()[i]) && approx_loe(vector.Modulus(), sides()[i].Modulus()) && approx_loe(0, vector * sides()[i])) {
          return true;
        }
      }
      size_t count = 0;
      Line control(point, Vector(0, -1));
      for (size_t i = 0; i < verticesCount(); ++i) {
        Point intersection = Line::Intersection(control, {vertices_[i], sides()[i]});
        if (intersection.x == DBL_MAX) {
          continue;
        }
        Vector vector(vertices_[i], intersection);
        if (approx_loe(vector.Modulus(), sides()[i].Modulus()) && approx_loe(0, vector * sides()[i]) && approx_less(intersection.y, point.y)) {
          ++count;
        }
      }
      return count % 2 == 1;
    }

    void rotate(const Point& center, double angle) override {
      for (Point& point : vertices_) {
        point.rotate(center, angle);
      }
    }

    void reflect(const Point& center) override {
      for (Point& point : vertices_) {
        point.reflect(center);
      }
    }

    void reflect(const Line& axis) override {
      for (Point& point : vertices_) {
        point.reflect(axis);
      }
    }

    void scale(const Point& center, double coefficient) override {
      for (Point& point : vertices_) {
        point.scale(center, coefficient);
      }
    }
};

class Ellipse : public Shape {
  protected:
    Point focus1_;
    Point focus2_;
    double focus_;
    double eccentricity_;
    double distance_;

    Ellipse() : focus1_(), focus2_(), focus_(0), eccentricity_(0), distance_(0) {}

  public:
    std::pair<Point,Point> focuses() const { return {focus1_, focus2_}; }

    std::pair<Line, Line> directrices() const {
      Vector vector1(focus1_, focus2_);
      Vector vector2(0.5 * (1 / (eccentricity_ * eccentricity_) - 1) * vector1);
      vector1 = Vector(vector1.y, -vector1.x);
      Point point1(focus2_.x + vector2.x, focus2_.y + vector2.y);
      Point point2(focus1_.x - vector2.x, focus1_.y - vector2.y);
      return {Line(point1, vector1), Line(point2, vector1)};
    }

    double eccentricity() const { return eccentricity_; }

    Point center() const { return {(focus1_.x + focus2_.x) / 2, (focus1_.y + focus2_.y) / 2}; }

    Ellipse(Point focus1, Point focus2, double distance) : focus1_(focus1), focus2_(focus2), focus_(Vector(focus1, focus2).Modulus() / 2), distance_(distance) {
      eccentricity_ = 2 * focus_ / distance;
    }

    double perimeter() const override {
      long double p1 = focus_ / eccentricity_ + focus_ * sqrtl(1 / eccentricity_ / eccentricity_ - 1);
      long double p2 = (focus_ / eccentricity_ - focus_ * sqrtl(1 / eccentricity_ / eccentricity_ - 1)) / p1;
      return M_PI * p1 * (1 + 3 * p2 * p2 / (10 + sqrtl(4 - 3 * p2 * p2)));
    }

    double area() const override {
      if (eccentricity_ == 0) { return M_PI * distance_ * distance_ / 4; }
      return M_PI * focus_ * focus_ / eccentricity_ * sqrtl(1 / eccentricity_ / eccentricity_ - 1);
    }

    bool operator==(const Shape& another) const override {
      auto ellipse = dynamic_cast<const Ellipse*>(&another);
      if (ellipse == nullptr) { return false; }
      if (!approx_equal(eccentricity_, ellipse->eccentricity_)) { return false; }
      return (focus1_ == ellipse->focus1_ && focus2_ == ellipse->focus2_) || (focus1_ == ellipse->focus2_ && focus2_ == ellipse->focus1_);
    }

    bool isSimilarTo(const Shape& another) const override {
      auto ellipse = dynamic_cast<const Ellipse*>(&another);
      if (ellipse == nullptr) { return false; }
      return approx_equal(ellipse->eccentricity_, eccentricity_);
    }

    bool containsPoint(const Point& point) const override {
      return approx_loe(Vector(focus1_, point).Modulus() + Vector(focus2_, point).Modulus(), distance_);
    }

    void rotate(const Point& center, double angle) override {
      focus1_.rotate(center, angle);
      focus2_.rotate(center, angle);
    }

    void reflect(const Point& center) override {
      focus1_.reflect(center);
      focus2_.reflect(center);
    }

    void reflect(const Line& axis) override {
      focus1_.reflect(axis);
      focus2_.reflect(axis);
    }

    void scale(const Point& center, double coefficient) override {
      focus1_.scale(center, coefficient);
      focus2_.scale(center, coefficient);
      focus_ = Vector(focus1_, focus2_).Modulus() / 2;
      distance_ = std::abs(distance_ * coefficient);
    }
};

class Circle : public Ellipse {
  public:
    Circle() : Ellipse() {}

    Circle(Point center_, double radius_) : Ellipse(center_, center_, 2 * radius_) {}

    double radius() const { return distance_ / 2; }
};

class Rectangle : public Polygon {
  public:
    Point center() const { return {(vertices_[0].x + vertices_[2].x) / 2, (vertices_[0].y + vertices_[2].y) / 2}; }

    std::pair<Line, Line> diagonals() const { return {Line(vertices_[0], vertices_[2]), Line(vertices_[1], vertices_[3])}; }

    Rectangle(Point point1, Point point3, double tan) {
      if (tan < 1) { tan = 1 / tan; }
      double cos = sqrt(1 / (1 + tan * tan));
      Vector diagonal(point1, point3);
      Vector side = cos * cos * Vector(diagonal.x + diagonal.y * tan, -diagonal.x * tan + diagonal.y);
      Point point2(point1.x + side.x, point1.y + side.y);
      Point point4(point3.x - side.x, point3.y - side.y);
      vertices_ = {point1, point2, point3, point4};
    }

};

class Square : public Rectangle {
  public:
    Square(Point point1, Point point3) : Rectangle(point1, point3, 1) {}

    double side() const { return Vector(vertices_[0], vertices_[1]).Modulus(); }

    Circle circumscribedCircle() const { return {center(), side() / sqrt(2)}; }

    Circle inscribedCircle() const { return {center(), side() / 2}; }
};

class Triangle : public Polygon {
    using Polygon::Polygon;
  public:
    Triangle(std::vector<Point> vertices) : Polygon(vertices) {}

    Triangle(std::initializer_list<Point> vertices) : Polygon(vertices) {}

    Point centroid() const { return {(vertices_[0].x + vertices_[1].x + vertices_[2].x) / 3, (vertices_[0].y + vertices_[1].y + vertices_[2].y) / 3}; }

    Point orthocenter() const {
      Vector vector1(vertices_[0], vertices_[1]);
      Vector vector2(vertices_[1], vertices_[2]);
      vector1 = vector1.Orthogonal();
      vector2 = vector2.Orthogonal();
      Line line1(vertices_[2], vector1);
      Line line2(vertices_[0], vector2);
      return Line::Intersection(line1, line2);
    }

    Circle inscribedCircle() const {
      Vector v1(vertices_[0], vertices_[1]);
      Vector v2(vertices_[0], vertices_[2]);
      Vector v3(vertices_[1], vertices_[2]);
      v1 = (1 / v1.Modulus()) * v1;
      v2 = (1 / v2.Modulus()) * v2;
      v3 = (1 / v3.Modulus()) * v3;
      Vector b1 = v1 + v2;
      Vector b2 = v3 - v1;
      Point center = Line::Intersection(Line(vertices_[0], b1), Line(vertices_[1], b2));
      return {center, 2 * area() / perimeter()};
    }

    Circle circumscribedCircle() const {
      Point p0 = vertices_[0];
      Point p1 = vertices_[1];
      Point p2 = vertices_[2];
      Vector v1 = Vector(p0, p1).Orthogonal();
      Vector v2 = Vector(p0, p2).Orthogonal();
      Line l1({(p0.x + p1.x) / 2, (p0.y + p1.y) / 2}, v1);
      Line l2({(p0.x + p2.x) / 2, (p0.y + p2.y) / 2}, v2);
      Point center = Line::Intersection(l1, l2);
      return {center, Vector(center, p0).Modulus()};
    }

    Line EulerLine() const {
      return {centroid(), orthocenter()};
    }

    Circle ninePointsCircle() const {
      Point p0 = vertices_[0];
      Point p1 = vertices_[1];
      Point p2 = vertices_[2];
      Point m1((p0.x + p1.x) / 2, (p0.y + p1.y) / 2);
      Point m2((p0.x + p2.x) / 2, (p0.y + p2.y) / 2);
      Point m3((p2.x + p1.x) / 2, (p2.y + p1.y) / 2);
      return Triangle{m1, m2, m3}.circumscribedCircle();
    }
};