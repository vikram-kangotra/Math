#pragma once

#include <math.h>
#include <ostream>

template <int rows, typename T = double>
struct Vec
{
    Vec() {}

	Vec(std::initializer_list<T> list) {
		for(int i=0; const auto& l: list) {
            if (i == rows) {
                throw "Extra elements were provided to Vec";
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

	Vec& operator=(const Vec& other) {
		for(int i=0; i<length(); ++i) {
			m_data[i] = other[i];
		}
		return *this;
	}

	Vec operator-() const {
		auto vec{*this};
		for(int i=0; i<length(); ++i) {
			vec[i] = -m_data[i];
		}
		return vec;
	}

	Vec& operator+=(const Vec& other) {
		for(int i=0; i<length(); ++i) {
			m_data[i] += other[i];
		}
		return *this;
	}

	Vec operator+(const Vec& other) const {
		auto vec{*this};
		return vec+=other;
	}

	Vec& operator-=(const Vec& other) {
		return (*this) += -other;
	}

	Vec operator-(const Vec& other) const {
		auto vec{*this};
		return vec + (-other);
	}

	Vec operator*(const T& val) const {
		auto vec{*this};
		for(int i=0; i<length(); ++i) {
			vec[i] *= val;
		}
		return vec;
	}

	Vec& operator*=(const Vec& other) {
		for(int i=0; i<length(); ++i) {
			m_data[i] *= other[i];
		}
		return *this;
	}

	Vec operator*(const Vec& other) const {
		auto vec{*this};
		return vec*=other;
	}

	Vec operator/(const T& val) const {
		if(val == 0)
			std::runtime_error("divisibility by 0.");
		return (*this) * (1/val);
	}

	T dot(const Vec& other) const {
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

	void setDirection(const Vec<rows>& vec) {
		*this = getMagnitude() * vec;
	}

	 auto getMagnitude() const {
		return sqrt(dot(*this));
	}

	Vec<rows> getDirection() const {
		return (*this)/getMagnitude();
	}

private:
	T m_data[rows];
};

template <int rows, typename T>
std::ostream& operator<<(std::ostream& out, const Vec<rows, T>& vec)
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
