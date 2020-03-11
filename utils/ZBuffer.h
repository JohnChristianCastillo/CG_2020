//
// Created by reed on 24.05.19.
//

#ifndef UTILS_ZBUFFER_H
#define UTILS_ZBUFFER_H

#include <iostream>
#include <vector>
#include <limits>

class ZBuffer: public std::vector<std::vector<double>> {

public:

    int imageWidth;
    int imageHeight;

    ZBuffer(const int width, const int height):imageWidth(width), imageHeight(height),
    vector (height, vector<double>(width, std::numeric_limits<double>::infinity())){};

    bool allowedByZBuffer(double z, unsigned int x, unsigned int y){
        return (double) (1 / z) < (*this)[y][x];
    }
};

#endif //UTILS_ZBUFFER_H
