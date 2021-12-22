#pragma once

#include <util/math.h>
#include <iostream>

template <int rows, typename T = double>
struct Vector
{
    Vector() 
    : m_data{} {}

	Vector(std::initializer_list<T> list) {
		for(int i=0; const auto& l: list) {
            if (i == rows) {
                throw "Extra elements were provided to Vector";
            }
			m_data[i] = l;
			++i;
		}
	}

	T& operator[](int index) {
		return m_data[index];
	}

	const T& operator[](int index) const { return m_data[index]; }

	T* begin() { return m_data; }
    const T* begin() const { return m_data; }

    T* end() { return m_data + rows; }
    const T* end() const { return m_data + rows;  }

	constexpr int size() const { return rows;	}

	Vector& operator=(const Vector& other) {
		for(int i=0; i<rows; ++i) {
			m_data[i] = other[i];
		}
		return *this;
	}

	Vector operator-() const {
		auto vec{*this};
		for(int i=0; i<rows; ++i) {
			vec[i] = -m_data[i];
		}
		return vec;
	}

	Vector& operator+=(const Vector& other) {
		for(int i=0; i<rows; ++i) {
			m_data[i] += other[i];
		}
		return *this;
	}

	Vector operator+(const Vector& other) const {
		auto vec{*this};
		return vec+=other;
	}

	Vector& operator-=(const Vector& other) {
		return (*this) += -other;
	}

	Vector operator-(const Vector& other) const {
		auto vec{*this};
		return vec + (-other);
	}

	Vector operator*(const T& val) const {
		auto vec{*this};
		for(int i=0; i<rows; ++i) {
			vec[i] *= val;
		}
		return vec;
	}

	Vector& operator*=(const Vector& other) {
		for(int i=0; i<rows; ++i) {
			m_data[i] *= other[i];
		}
		return *this;
	}

	Vector operator*(const Vector& other) const {
		auto vec{*this};
		return vec*=other;
	}

    Vector& operator/=(const T& val) {
        return (*this) *= (1/val);
    }

	Vector operator/(const T& val) const {
		if(val == 0)
			std::runtime_error("divisibility by 0.");
		return (*this) * (1/val);
	}

	T dot(const Vector& other) const {
		T sum = 0;
		for(int i=0; i<size(); ++i) {
			sum += m_data[i] * other[i];
		}
		return sum;
	}	

	void setLength(double len) {
		(*this) = len * normalize();
	}

	void setDirection(const Vector<rows>& vec) {
		*this = length() * vec;
	}

	auto length() const {
		return sqrt(dot(*this));
	}

    auto lengthSq() const {
        return dot(*this);
    }

    static float getAngle(const Vector& a, const Vector& b) {
        auto dotRet = a.dot(b)/(a.length() * b.length());
        return std::acos(dotRet);
    }

	Vector<rows>& normalize() const {
		return (*this) /= length();
	}

    static Vector<rows> normalize(const Vector& vec) {
        return vec / vec.length();
    }

private:
	T m_data[rows];
};

template <int rows, typename T>
std::ostream& operator<<(std::ostream& out, const Vector<rows, T>& vec)
{
	std::cout << '|';
	for(int i=0; i<rows; ++i) {
		std::cout << vec[i];
		if(i<rows-1)
			std::cout << ',';
	}
	std::cout << '|';
	return out;
}

template <typename T>
class Vector2 : public Vector<2, T> {
    public:
        Vector2() {}

        Vector2(const T& x, const T& y)
        : Vector<2, T>{x, y} {}

        Vector2(const Vector2& other)
        : Vector<2, T>({other[0], other[1]}) {}

        Vector2& operator=(const Vector2& other) {
            (*this)[0] = other[0];
            (*this)[1] = other[1];
            return *this;
        }

        Vector2(Vector2&& other)
        : Vector<2, T>{other[0], other[1]} {
            other[0] = 0;
            other[1] = 0;
        }

    public:

        T& x = (*this)[0];
        T& y = (*this)[1];
};

template <typename T>
class Vector3 : public Vector<3, T> {
    public:
        Vector3() {}

        Vector3(const T& x, const T& y, const T& z)
        : Vector<3, T>{x, y, z} {}

        Vector3(const Vector3& other)
        : Vector<3, T>{other[0], other[1], other[2]} {}

        Vector3& operator=(const Vector3& other) {
            (*this)[0] = other[0];
            (*this)[1] = other[1];
            (*this)[2] = other[2];
            return *this;
        }

        Vector3(Vector3&& other)
        : Vector<3, T>{other[0], other[1], other[2]} {
            other[0] = 0;
            other[1] = 0;
            other[2] = 0;
        }

    public:

        T& x = (*this)[0];
        T& y = (*this)[1];
        T& z = (*this)[2];
};

