/*
 * Color.h
 *
 *  Created on: Feb 22, 2019
 *      Author: student
 */

#ifndef COLOR_H_
#define COLOR_H_
#include "easy_image.h"
#include <cmath>
inline int roundToInt(double d){
    return static_cast<int>(std::round(d));
}

class Color{
	double red;
	double green;
	double blue;
public:
	Color():red(0), green(0), blue(0){}
	Color(double r, double g, double b):red(r), green(g), blue(b){}

	double getBlue() const {
		return blue;
	}

	void setBlue(double blue) {
		this->blue = blue;
	}

	double getGreen() const {
		return green;
	}

	void setGreen(double green) {
		this->green = green;
	}

	double getRed() const {
		return red;
	}

	void setRed(double red) {
		this->red = red;
	}

	const img::Color getColor() const {
	    return  img::Color(roundToInt(red*255), roundToInt(green*255), roundToInt(blue*255));
	}
};




#endif /* COLOR_H_ */
