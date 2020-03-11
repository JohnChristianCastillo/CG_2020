//
// Created by reed on 10.05.19.
//

#ifndef UTILS_FACE_H
#define UTILS_FACE_H


#include <vector>

class Face {  //Face == Vlak

public:
    //initializes Face using 2 pointIndex-pair parameters (started using for Cube (3D Figures))
    Face(int index1, int index2){
        point_indexes = {index1, index2};
    };
    Face(){};
    Face(std::vector<int> indexes):point_indexes(indexes){};

    //De indexen refereren naar
    //punten in de ‘points’ vector
    //van de Figure-klasse
    std::vector <int> point_indexes;

    const std::vector<int> &getPointIndexes() const;

    void setPointIndexes(const std::vector<int> &pointIndexes);
};


#endif //UTILS_FACE_H
