#pragma once

#include <math.h>
#include <iostream>

template <int rows, typename T = double>
struct Vector
{
    Vector() {}

	Vector(std::initializer_list<T> list) {
		for(int i=0; const auto& l: list) {
            if (i == rows) {
                throw "Extra elements were provided to Vector";
            }
			m_data[i] = l;
			++i;
		}
	}

	template <typename... R>
	void set(R&&... elems) {
		int index = 0;
		((m_data[index] = elems, ++index),...);
	}

	T& operator[](int index) {
		return m_data[index];
	}

	const T& operator[](int index) const {
		return m_data[index];
	}

	T* begin() {
		return m_data;
	}

    const T* begin() const {
        return m_data;
    }

    T* end() {
        return m_data + rows;
    }

    const T* end() const {
        return m_data + rows; 
    }

	constexpr int length() const {
		return rows;
	}

	Vector& operator=(const Vector& other) {
		for(int i=0; i<length(); ++i) {
			m_data[i] = other[i];
		}
		return *this;
	}

	Vector operator-() const {
		auto vec{*this};
		for(int i=0; i<length(); ++i) {
			vec[i] = -m_data[i];
		}
		return vec;
	}

	Vector& operator+=(const Vector& other) {
		for(int i=0; i<length(); ++i) {
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
		for(int i=0; i<length(); ++i) {
			vec[i] *= val;
		}
		return vec;
	}

	Vector& operator*=(const Vector& other) {
		for(int i=0; i<length(); ++i) {
			m_data[i] *= other[i];
		}
		return *this;
	}

	Vector operator*(const Vector& other) const {
		auto vec{*this};
		return vec*=other;
	}

	Vector operator/(const T& val) const {
		if(val == 0)
			std::runtime_error("divisibility by 0.");
		return (*this) * (1/val);
	}

	T dot(const Vector& other) const {
		T sum = 0;
		for(int i=0; i<length(); ++i) {
			sum += m_data[i] * other[i];
		}
		return sum;
	}	

	void setMagnitude(double mag) {
		(*this) = mag * getDirection();
	}

	template <typename... R>
	void setDirectionAngle(R&&... elems) {
		int index = 0;
		(((*this)[index] = getMagnitude() * cos(elems), ++index), ...);
	}

	template <typename... R>
	void setDirection(R&&... elems) {
		int index = 0;
		(((*this)[index] = getMagnitude() * elems, ++index), ...);
	}

	void setDirection(const Vector<rows>& vec) {
		*this = getMagnitude() * vec;
	}

	 auto getMagnitude() const {
		return sqrt(dot(*this));
	}

	Vector<rows> getDirection() const {
		return (*this)/getMagnitude();
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

template <typename T = double>
struct Vector2 : public Vector<2, T> {
    
    T& x;
    T& y;     

    Vector2(T _x = 0, T _y = 0)
    : x{(*this)[0]}, y{(*this)[1]} {
        x = _x;
        y = _y;
    }

    Vector2& operator=(const Vector2<T>& vec) {
        x = vec.x;
        y = vec.y;
        return (*this);
    }
};
