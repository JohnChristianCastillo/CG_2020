//
// Created by reed on 10.05.19.
//

#include "Figure.h"

const std::string &Figure::getType() const {
    return type;
}

void Figure::setType(const std::string &type) {
    Figure::type = type;
}


double Figure::getScale() const {
    return scale;
}

void Figure::setScale(double scale) {
    Figure::scale = scale;
}

const std::vector<Vector3D> &Figure::getPoints() const {
    return points;
}

void Figure::addPoint(Vector3D pointVar) {
    points.push_back(pointVar);
}

const std::vector<Face* > &Figure::getFaces() const {
    return faces;
}

void Figure::addFace(Face* face) {
    Figure::faces.push_back(face);
}

const Color &Figure::getColor() const {
    return color;
}

void Figure::setColor(Color &vColor) {
    Figure::color = vColor;
}
