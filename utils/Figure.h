//
// Created by reed on 10.05.19.
//

#ifndef UTILS_FIGURE_H
#define UTILS_FIGURE_H


#include <vector>
#include "vector/vector3d.h"
#include "Face.h"
#include "Color.h"

class Figure {
public:
    std::string type;
    double rotateX;
    double rotateY;
    double rotateZ;
    double scale = 1.0;
    Color color;

    Vector3D center;

    std::vector<Vector3D> points;
    std::vector<Face*> faces;



    const std::string &getType() const;

    void setType(const std::string &type);



    double getScale() const;

    void setScale(double scale);

    const std::vector<Vector3D> &getPoints() const;

    void addPoint(Vector3D pointVar);

    const std::vector<Face* > &getFaces() const;

    void addFace(Face* face);

    const Color &getColor() const;

    void setColor(Color &vColor);





};


#endif //UTILS_FIGURE_H
