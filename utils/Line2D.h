//
// Created by reed on 17.05.19.
//

#ifndef UTILS_LINE2D_H
#define UTILS_LINE2D_H


#include "Color.h"
#include "Point2D.h"

class Line2D {
public:
    Point2D p1;
    Point2D p2;
    Color color;
    double z1 = 0;
    double z2 = 0;
    Line2D();

    Line2D(Point2D point1, Point2D point2, Color colorInput = Color()): p1(point1), p2(point2), color(colorInput){};

    const img::Color getColor() const {
        return color.getColor();
    }

    const Point2D& getP1() const {
        return p1;
    }

    const Point2D& getP2() const {
        return p2;
    }

};


#endif //UTILS_LINE2D_H
