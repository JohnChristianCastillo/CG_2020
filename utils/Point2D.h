#ifndef POINT2D_H_
#define POINT2D_H_


class Point2D{
	double x;
	double y;
public:
    Point2D(){};
	Point2D(double ix, double iy):x(ix), y(iy){}

	double getX() const {
		return x;
	}

	void setX(double x) {
		this->x = x;
	}

	double getY() const {
		return y;
	}

	void setY(double y) {
		this->y = y;
	}
};

#endif /* POINT_H_ */
