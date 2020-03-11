//
// Created by reed on 10.05.19.
//

#include "Face.h"

const std::vector<int> &Face::getPointIndexes() const {
    return point_indexes;
}

void Face::setPointIndexes(const std::vector<int> &pointIndexes) {
    point_indexes = pointIndexes;
}
