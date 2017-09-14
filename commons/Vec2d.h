#ifndef VEC2D_HPP
#define VEC2D_HPP

#include <iostream>

struct Vec2d {
	double x,y;
	
	// default + parameterized constructor
	Vec2d(double x=0, double y=0) : x(x), y(y) {}
	
	inline Vec2d& operator=(const Vec2d& a) {
            x=a.x;
            y=a.y;
            return *this;
	}
	
	//************* ARITHMATIC ****************** //
	//vector [OP] vector
	inline Vec2d operator+(const Vec2d& a) const {
		return Vec2d(x+a.x, y+a.y);
	}

	inline Vec2d operator-(const Vec2d& a) const {
		return Vec2d(x-a.x, y-a.y);
	}
	
	// vector [OP] scalar
        inline Vec2d operator*(const double& a) const {
		return Vec2d(x*a, y*a);
	}
	
        inline Vec2d operator/(const double& a) const {
		return Vec2d(x/a, y/a);
	}
	
	// ************** ASSINGMENTS **************** //
	//vector [OP] vector
	inline Vec2d& operator+=(const Vec2d& a) {
		x+=a.x;
		y+=a.y;
		return *this;
	}
	
	inline Vec2d& operator-=(const Vec2d& a) {
		x-=a.x;
		y-=a.y;
		return *this;
	}
	
	// vector [OP] scalar
        inline Vec2d& operator+=(const double& a) {
		x+=a;
		y+=a;
		return *this;
	}

        inline Vec2d& operator-=(const double& a) {
		x-=a;
		y-=a;
                return *this;
	}
	
        inline Vec2d& operator/=(const double &a) {
		x/=a;
		y/=a;
		return *this;
	}
	
        inline Vec2d& operator*=(const double &a) {
		x*=a;
		y*=a;
		return *this;
	}
	
	
};

//std::ostream& operator << (std::ostream& os, const Vec2d &p) {
//    os <<"[" << p.x << ", " << p.y << "]";
//    return os;
//}

#endif
