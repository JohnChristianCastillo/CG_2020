//
// Created by reed on 09.03.19.
//

#ifndef UTILS_EXTRAFUNCTIONS_H
#define UTILS_EXTRAFUNCTIONS_H

#include <iostream>
#include <cmath>
#include "vector/vector3d.h"
#include "Figure.h"
#include "Line2D.h"
#include "ZBuffer.h"
#include <list>
#include <algorithm>

typedef std::list<Figure*> Figures3D;





Matrix rotateX(const double angle, std::string radOrDeg = "Degrees"){
    double radianAngle = angle;
    if(radOrDeg == "Degrees"){
        radianAngle = (angle*M_PI)/180;
    }
    Matrix Mx;
    Mx(2,2) = cos(radianAngle);
    Mx(2,3) = sin(radianAngle);
    Mx(3,2) = -sin(radianAngle);
    Mx(3,3) = cos(radianAngle);
    //for(int i = 0; i < points.size(); i++){
    //    std::cout << "old point: " << points[i].x << ", " << points[i].y << ", " << points[i].z;
    //    points[i] = points[i]*Mx;
    //    std::cout << "  ------>" << points[i].x << ", " << points[i].y << ", " << points[i].z << std::endl;
    //}
    //std::cout << "Mx is: \n"<< Mx  << std::endl;
    return Mx;

}
Matrix rotateY(const double angle, std::string radOrDeg = "Degrees"){
    double radianAngle = angle;
    if(radOrDeg == "Degrees"){
        radianAngle = (angle*M_PI)/180;
    }
    Matrix My;
    My(1,1) = cos(radianAngle);
    My(1,3) = -sin(radianAngle);
    My(3,1) = sin(radianAngle);
    My(3,3) = cos(radianAngle);
    //for(int i = 0; i < points.size(); i++){
    //    std::cout << "old point: " << points[i].x << ", " << points[i].y << ", " << points[i].z;
    //    points[i] = points[i]*My;
    //    std::cout << "  ------>" << points[i].x << ", " << points[i].y << ", " << points[i].z << std::endl;
    //}
    //std::cout << "My is: \n"<< My  << std::endl;
    return My;

}
Matrix rotateZ(const double angle, std::string radOrDeg = "Degrees"){
    double radianAngle = angle;
    if(radOrDeg == "Degrees"){
        radianAngle = (angle*M_PI)/180;
    }
    Matrix Mz;
    Mz(1,1) = cos(radianAngle);
    Mz(1,2) = sin(radianAngle);
    Mz(2,1) = -sin(radianAngle);
    Mz(2,2) = cos(radianAngle);
    //for(int i = 0; i < points.size(); i++){
    //    std::cout << "old point: " << points[i].x << ", " << points[i].y << ", " << points[i].z;
    //    points[i] = points[i]*Mz;
    //    std::cout << "  ------>" << points[i].x << ", " << points[i].y << ", " << points[i].z << std::endl;
    //}

    //std::cout << "Mz is: \n"<< Mz  << std::endl;
    return Mz;
}
Matrix translate(const Vector3D &vector){
    //std::cout << "------------------------BEGIN TRANSLATION---------------------- \n";
    Matrix tMatrix;
    tMatrix(4,1) = vector.x;
    tMatrix(4,2) = vector.y;
    tMatrix(4,3) = vector.z;
    //for(int i = 0; i < points.size(); i++){
    //    std::cout << "old point: " << points[i].x << ", " << points[i].y << ", " << points[i].z;
    //    points[i] = points[i]*tMatrix;
    //    std::cout << "  ------>" << points[i].x << ", " << points[i].y << ", " << points[i].z << std::endl;
    //}
    return tMatrix;

}
Matrix scaleFigure(const double scale){
    Matrix scaleMatrix;
    scaleMatrix(1,1) = scale;
    scaleMatrix(2,2) = scale;
    scaleMatrix(3,3) = scale;

    return scaleMatrix;
}
void toPolar(const Vector3D &eyePoint, double &theta, double &phi, double &r){
    //x = rsin(Ф)*cos(θ)
    //y = rsin(Ф)*sin(θ)
    //z = rcos(Ф)
    //calculate r  // r = sqrt(x²+y²+z²)
    r = sqrt(pow(eyePoint.x, 2) + pow(eyePoint.y, 2) + pow(eyePoint.z, 2));
    //calculate θ  // arctan(θ) range is ]-pi/2, pi/2[
    theta = atan2(eyePoint.y, eyePoint.x);
    //calculate Ф
    phi = acos(eyePoint.z/r);
}


void applyTransformation(Figure &figure){
    Matrix transformatieMatrix;
    transformatieMatrix = scaleFigure(figure.scale)*rotateX(figure.rotateX)*rotateY(figure.rotateY)*rotateZ(figure.rotateZ)*translate(figure.center);
    for(unsigned int i = 0; i < figure.points.size();i++){
        figure.points[i] = figure.points[i]*transformatieMatrix;
    }
}


Matrix eyePointTrans(const Vector3D &eyepoint, const std::vector<double> &viewDirection, const bool vdExists){
    double theta;
    double phi;
    double r;
    if(vdExists){
        const Vector3D VD = Vector3D::vector(viewDirection[0],viewDirection[1],viewDirection[2]);
        toPolar(-VD, theta, phi, r);
                //camera wordt nieuwe O van coordinaatsysteem
        auto t = translate(Vector3D::vector(0,0,0)-eyepoint);
        auto u = rotateZ(-(theta+(M_PI/2)), "radians");
        auto v = rotateX(-phi, "radians");
        return translate(Vector3D::vector(0,0,0)-eyepoint)*rotateZ((-M_PI/2)-theta, "radians")*rotateX(-phi, "radians");
    }
    //else
    toPolar(eyepoint, theta, phi, r);
    return rotateZ((-M_PI/2)-theta, "radians")*rotateX(-phi, "radians")*translate(Vector3D::vector(0,0,-r));

}
typedef std::list<Figure*> Figures3D;

void applyTransformation(Figures3D &figures3D){
    for(const auto& figure : figures3D){
        applyTransformation(*figure);
    }

}

typedef std::list<Line2D> Lines2D;

Point2D doProjection(const Vector3D &point, const double d = 1){   // EVENTUEEL , const double d = 1  als parameter toevoegen
    //calculate the 2D x and y coordinates given the 3D x,y and z coordinates
    double x = d*point.x/(-point.z); //new x = d*oldX/(-oldZ)
    double y = d*point.y/(-point.z); //new y = d*oldy/(-oldZ)

    //make point out of the calculated 2D x and y coordinates
    return Point2D(x, y);


};
Lines2D doProjection(const Figures3D& figures3D, const Vector3D &eye,
                     const std::vector<double> &viewDirection = {}, bool vdExists = false){
    Matrix V = eyePointTrans(eye, viewDirection, vdExists);
    Lines2D lines2D;

    for(const auto& figure : figures3D){
        //iterate over every figure points and multiply V (eyepoint transformation to each point)
        /*for(unsigned int i = 0; i < figure->points.size(); i++){
            figure->points[i] = figure->points[i]*V;
        }*/

        //make lines out of every figure in figures3D and put it in Lines2D which will contain all lines in the end
        for(unsigned int j = 0; j<figure->faces.size(); j++){
            std::vector<int> currentPointPair = figure->faces[j]->point_indexes;
            //loop over every point and make a line out of it
            for(unsigned int i = 0; i < currentPointPair.size(); i++) {
                Point2D point1 = doProjection(figure->points[currentPointPair[i]]);
                Point2D point2 = doProjection(figure->points[currentPointPair[(i+1)%currentPointPair.size()]]);
                //index is like this so if we're at the last index then we make a line out of last and first index

                //make line out of the point pair
                Line2D tempLine = Line2D(point1, point2, figure->color);
                (tempLine.z1) = figure->points[currentPointPair[i]].z;   //Zbuff
                tempLine.z2 = figure->points[currentPointPair[(i+1)%currentPointPair.size()]].z;//Zbuff
                lines2D.push_back(tempLine);
            }
        }
    }
    return lines2D;
}
/*
Lines2D doClipping(const Figures3D& figures3D, const Vector3D &eye, const Vector3D &viewDirection, const double hfov
                    const double aspectRatio, const double dNear, const double dFar){

}*/


////////vanaf hier Genereren van 3D Figuren

void createCube(Figure* tempFig){
    tempFig->addPoint(Vector3D::point(1,-1,-1));
    tempFig->addPoint(Vector3D::point(-1,1,-1));
    tempFig->addPoint(Vector3D::point(1,1,1));
    tempFig->addPoint(Vector3D::point(-1,-1,1));
    tempFig->addPoint(Vector3D::point(1,1,-1));
    tempFig->addPoint(Vector3D::point(-1,-1,-1));
    tempFig->addPoint(Vector3D::point(1,-1,1));
    tempFig->addPoint(Vector3D::point(-1,1,1));
    //now add the faces to be matched
    std::vector<std::vector<int>> pointIndexVector;
    pointIndexVector.push_back({0,4,2,6});
    pointIndexVector.push_back({4,1,7,2});
    pointIndexVector.push_back({1,5,3,7});
    pointIndexVector.push_back({5,0,6,3});
    pointIndexVector.push_back({6,2,7,3});
    pointIndexVector.push_back({0,5,1,4});
    for(const auto & i : pointIndexVector){
        Face* tempFace = new Face(i);
        tempFig->faces.push_back(tempFace);
    }
}


void createTetrahedron(Figure* tempFig){
    tempFig->addPoint(Vector3D::point(-1,-1,1));
    tempFig->addPoint(Vector3D::point(-1,1,-1));
    tempFig->addPoint(Vector3D::point(1,1,1));
    tempFig->addPoint(Vector3D::point(1,-1,-1));
    //tempFig->addPoint(Vector3D::point(1,-1,-1));
    //tempFig->addPoint(Vector3D::point(-1,1,-1));
    //tempFig->addPoint(Vector3D::point(1,1,1));
    //tempFig->addPoint(Vector3D::point(-1,-1,1));

    //now add the points to be matched
    std::vector<std::vector<int>> pointIndexVector;
    pointIndexVector.push_back({2,1,0});
    pointIndexVector.push_back({2,3,1});
    pointIndexVector.push_back({1,3,0});
    pointIndexVector.push_back({3,2,0});
    for(const auto & i : pointIndexVector){
        Face* tempFace = new Face(i);
        tempFig->faces.push_back(tempFace);
    }
}

void createOctahedron(Figure* tempFig){
    tempFig->addPoint(Vector3D::point(1,0,0));
    tempFig->addPoint(Vector3D::point(0,1,0));
    tempFig->addPoint(Vector3D::point(-1,0,0));
    tempFig->addPoint(Vector3D::point(0,-1,0));
    tempFig->addPoint(Vector3D::point(0,0,-1));
    tempFig->addPoint(Vector3D::point(0,0,1));
    //now add the points to be matched
    std::vector<std::vector<int>> pointIndexVector;
    pointIndexVector.push_back({0,1,5});
    pointIndexVector.push_back({1,2,5});
    pointIndexVector.push_back({2,3,5});
    pointIndexVector.push_back({3,0,5});
    pointIndexVector.push_back({1,0,4});
    pointIndexVector.push_back({2,1,4});
    pointIndexVector.push_back({3,2,4});
    pointIndexVector.push_back({0,3,4});
    for(const auto & i : pointIndexVector){
        Face* tempFace = new Face(i);
        tempFig->faces.push_back(tempFace);
    }
}

void createIcosahedron(Figure* tempFig){
    tempFig->addPoint(Vector3D::point(0,0,sqrt(5)/2)); //punt1
    for(int i = 2;i<=6; i++){ //punten 2-6
        tempFig->addPoint(Vector3D::point(cos((i-2)*2*M_PI/5), sin((i-2)*2*M_PI/5), 0.5));
    }
    for(int i = 7; i <= 11;i++){ //punten 7-11
        tempFig->addPoint(Vector3D::point(cos((M_PI/5)+(i-7)*2*M_PI/5), sin((M_PI/5)+(i-7)*2*M_PI/5), -0.5));
    }
    tempFig->addPoint(Vector3D::point(0,0,-sqrt(5)/2)); //punt12

    //now add the points to be matched
    std::vector<std::vector<int>> pointIndexVector;
    pointIndexVector.push_back({0,1,2});
    pointIndexVector.push_back({0,2,3});
    pointIndexVector.push_back({0,3,4});
    pointIndexVector.push_back({0,4,5});
    pointIndexVector.push_back({0,5,1});
    pointIndexVector.push_back({1,6,2});
    pointIndexVector.push_back({2,6,7});
    pointIndexVector.push_back({2,7,3});
    pointIndexVector.push_back({3,7,8});
    pointIndexVector.push_back({3,8,4});
    pointIndexVector.push_back({4,8,9});
    pointIndexVector.push_back({4,9,5});
    pointIndexVector.push_back({5,9,10});
    pointIndexVector.push_back({5,10,1});
    pointIndexVector.push_back({1,10,6});
    pointIndexVector.push_back({11,7,6});
    pointIndexVector.push_back({11,8,7});
    pointIndexVector.push_back({11,9,8});
    pointIndexVector.push_back({11,10,9});
    pointIndexVector.push_back({11,6,10});


    for(const auto & i : pointIndexVector){
        Face* tempFace = new Face(i);
        tempFig->faces.push_back(tempFace);
    }
}

void createDodecahedron(Figure* tempFig){
    Figure* icosahedron = new Figure;
    createIcosahedron(icosahedron);

    //calculate the middlepoint of each triangle, this will be the new point of the dodecahedron
    // there are 20 triangles in a icosahedron, there will thus be 20 points for dodecahedron
    for(unsigned int i = 0; i < icosahedron->faces.size(); i++){
        //loop over faces(triangles) of the icosahedron and calculate the middlepoint
        //we know there will be 3 elements in every face element since it's a triangle
        std::vector<int> tIndex = icosahedron->faces[i]->point_indexes; //triangle indices
        double xcor = (icosahedron->points[tIndex[0]].x
                    + icosahedron->points[tIndex[1]].x
                    + icosahedron->points[tIndex[2]].x)/3;
        double ycor = (icosahedron->points[tIndex[0]].y
                       + icosahedron->points[tIndex[1]].y
                       + icosahedron->points[tIndex[2]].y)/3;
        double zcor = (icosahedron->points[tIndex[0]].z
                       + icosahedron->points[tIndex[1]].z
                       + icosahedron->points[tIndex[2]].z)/3;
        Vector3D tempPoint =  Vector3D::point(xcor,ycor,zcor);
        tempFig->addPoint(tempPoint);
    }
    //now add the points to be matched
    std::vector<std::vector<int>> pointIndexVector;
    pointIndexVector.push_back({0,1,2,3,4});
    pointIndexVector.push_back({0,5,6,7,1});
    pointIndexVector.push_back({1,7,8,9,2});
    pointIndexVector.push_back({2,9,10,11,3});
    pointIndexVector.push_back({3,11,12,13,4});
    pointIndexVector.push_back({4,13,14,5,0});
    pointIndexVector.push_back({19,18,17,16,15});
    pointIndexVector.push_back({19,14,13,12,18});
    pointIndexVector.push_back({18,12,11,10,17});
    pointIndexVector.push_back({17,10,9,8,16});
    pointIndexVector.push_back({16,8,7,6,15});
    pointIndexVector.push_back({15,6,5,14,19});
    for(const auto & i : pointIndexVector){
        Face* tempFace = new Face(i);
        tempFig->faces.push_back(tempFace);
    }
    delete icosahedron; //delete ongebruikte object
}
void createSphere(Figure* tempFig, int n = 1){
    //make copy of figure, empty it's points and faces, calculate new points and faces
    while(n != 0){
    std::vector<Vector3D> copyPoints = tempFig->points;
    tempFig->points = {};

    std::vector<Face*> tempFaces = {}; //this will become our final faces object
    //calculate the 4 triangles per face
    for(unsigned int i = 0; i < tempFig->faces.size(); i++){
        //we know there will be 3 elements in every face element since it's a triangle
        std::vector<int> tIndex = tempFig->faces[i]->point_indexes; //triangle indices
        //add original points
        tempFig->addPoint(copyPoints[tIndex[0]]);
        int origIndex1 = tempFig->points.size()-1;
        tempFig->addPoint(copyPoints[tIndex[1]]);
        int origIndex2 = tempFig->points.size()-1;
        tempFig->addPoint(copyPoints[tIndex[2]]);
        int origIndex3 = tempFig->points.size()-1;
        //calculate point 1
        double xcor1 = (copyPoints[tIndex[0]].x + copyPoints[tIndex[1]].x)/2;
        double ycor1 = (copyPoints[tIndex[0]].y + copyPoints[tIndex[1]].y)/2;
        double zcor1 = (copyPoints[tIndex[0]].z + copyPoints[tIndex[1]].z)/2;
        Vector3D tempPoint =  Vector3D::point(xcor1,ycor1,zcor1);
        tempFig->addPoint(tempPoint);
        int index1 = tempFig->points.size()-1; //will come in handy for pairing faces
        //calculate point 2
        double xcor2 = (copyPoints[tIndex[0]].x + copyPoints[tIndex[2]].x)/2;
        double ycor2 = (copyPoints[tIndex[0]].y + copyPoints[tIndex[2]].y)/2;
        double zcor2 = (copyPoints[tIndex[0]].z + copyPoints[tIndex[2]].z)/2;
        tempPoint =  Vector3D::point(xcor2,ycor2,zcor2);
        tempFig->addPoint(tempPoint);
        int index2 = tempFig->points.size()-1; //will come in handy for pairing faces

        //calculate point 3
        double xcor3 = (copyPoints[tIndex[1]].x + copyPoints[tIndex[2]].x)/2;
        double ycor3 = (copyPoints[tIndex[1]].y + copyPoints[tIndex[2]].y)/2;
        double zcor3 = (copyPoints[tIndex[1]].z + copyPoints[tIndex[2]].z)/2;
        tempPoint =  Vector3D::point(xcor3,ycor3,zcor3);
        tempFig->addPoint(tempPoint);
        int index3 = tempFig->points.size()-1; //will come in handy for pairing faces

        //start matching the points to faces

        Face* tempFace = new Face({origIndex1,index1,index2});
        tempFaces.push_back(tempFace);

        tempFace = new Face({origIndex2,index1,index3});
        tempFaces.push_back(tempFace);

        tempFace = new Face({origIndex3,index2,index3});
        tempFaces.push_back(tempFace);

        tempFace = new Face({index1,index2,index3});//inner triangle
        tempFaces.push_back(tempFace);
    }
    //by the end of this for-loop, we have every points in our private data member


    //now we need to overwrite our saved faces data member by our accumulated face object
    for (auto it : tempFig->faces){delete it;} //clear pointers
    tempFig->faces = tempFaces;
    n = n-1;
    }
    //re-scale every point
    for(auto & i : tempFig->points){
        Vector3D point = i;
        double radius = sqrt(pow(point.x, 2) + pow(point.y, 2) + pow(point.z, 2));
        i.x = i.x/radius;
        i.y = i.y/radius;
        i.z = i.z/radius;
    }
}

void createCone(const int n, const double h, Figure* tempFig){
    tempFig->addPoint(Vector3D::point(0,0,h));  //pn (top point)
    for(int i = 1; i < n+1; i++){
        tempFig->addPoint(Vector3D::point(cos(2*(i)*M_PI/(n)),sin(2*(i)*M_PI/(n)),0));
    }
    //now add the points to be matched
    std::vector<std::vector<int>> pointIndexVector;
    for(int i = 1; i < tempFig->points.size(); i++){
        if(i+1 == tempFig->points.size()){
            pointIndexVector.push_back({i,1,0});
            break;
        }
        pointIndexVector.push_back({i,i+1,0});
    }
    std::vector<int> tempv;
    for(int i = 0; i<n; i++){
        tempv.push_back(i+1);
    }
    pointIndexVector.push_back(tempv);
    for(const auto & i : pointIndexVector){
        Face* tempFace = new Face(i);
        tempFig->faces.push_back(tempFace);
    }
}

void createCylinder(const int n, double h, Figure* tempFig){
    for(int i = 1; i <= n+1; i++){
        tempFig->addPoint(Vector3D::point(cos(2*(i+1)*M_PI/(n)),sin(2*(i+1)*M_PI/(n)),0)); //grondvlak punten
        tempFig->addPoint(Vector3D::point(cos(2*(i+1)*M_PI/(n)),sin(2*(i+1)*M_PI/(n)),h)); //bovenvlak punten
    }
    //now add the points to be matched
    std::vector<std::vector<int>> pointIndexVector;
    for(int i = 0; i < tempFig->points.size()-1; i++){
        pointIndexVector.push_back({i,(i+1), (i+3)%int(tempFig->points.size()), (i+2)%int(tempFig->points.size())});
        //pointIndexVector.push_back({i,(i+2)%int(tempFig->points.size())});
        //pointIndexVector.push_back({i+1,(i+3)%int(tempFig->points.size())});
        i++; //point pairs are currently next to each other
    }

    //oppervlakken
    std::vector<int> tempv;
    std::vector<int> tempv2;
    for(int i = 0; i<n; i++){
        tempv.push_back(i*2+1);
        tempv2.push_back(i*2);
    }
    pointIndexVector.push_back(tempv2);
    pointIndexVector.push_back(tempv);
    for(const auto &i: pointIndexVector){
        Face* tempFace = new Face(i);
        tempFig->faces.push_back(tempFace);
    }
}

void createTorus(const double r, const double R,  int n,  int m, Figure* tempFig){
    //met R afstand van 0punt tot centrum van buis
    // en r straal van buis
    // n aantal cirkel
    // m aantal punten per cirkel

    std::vector<std::vector<int>> pointIndexVector;
    for(int i = 0; i < n; i++){ //loop over number of circles
        double u = 2.0*i*M_PI/n;
        for(int j = 0; j < m; j++){ //loop over number of points in circle
            double v = 2.0*j*M_PI/m;
            double xuv = (R+r*cos(v))*cos(u);  //met u,v E [0,2pi]
            double yuv = (R+r*cos(v))*sin(u);
            double zuv = r*sin(v);
            tempFig->addPoint(Vector3D::point(xuv,yuv,zuv));
        }
    }
    for(unsigned int i= 0; i<n; ++i){
        for(unsigned int j = 0; j<m;++j){
            pointIndexVector.push_back({i*m+j,
                                        (i+1)%n*m+j,
                                        (i+1)%n*m+(j+1)%m,
                                        i*m+(j+1)%m
                                       });
        }
    }


    for(const auto &i: pointIndexVector){
        Face* tempFace = new Face(i);
        tempFig->faces.push_back(tempFace);
    }
}

std::vector<Face*> triangulate(const Face* face){
    int n=1;
    std::vector<Face*> tempFaces;
    while(n<=face->getPointIndexes().size()-2){
        std::vector<int> indices = {face->point_indexes[0],face->point_indexes[n],face->point_indexes[n+1]};
        Face* tempFace = new Face(indices);
        tempFaces.push_back(tempFace);
        n+=1;
    }
    return tempFaces;
}

std::vector<double> calculateInitial_d_dx_dy(const Lines2D &lines, double size){
    //Initialize xmin xmax ymin ymax
    double currentXMax = lines.begin()->getP1().getX();
    double currentYMax = lines.begin()->getP1().getY();

    double currentXMin = lines.begin()->getP1().getX();
    double currentYMin = lines.begin()->getP1().getY();
    for(const Line2D& line : lines){
        if(line.getP1().getX() > currentXMax){ // find maxX in p1
            currentXMax = line.getP1().getX();
        }
        if(line.getP1().getX() < currentXMin){ // find minX in p1
            currentXMin = line.getP1().getX();
        }
        if(line.getP1().getY() > currentYMax){ // find maxX in p2
            currentYMax = line.getP1().getY();
        }
        if(line.getP1().getY() < currentYMin){ // find minX in p2
            currentYMin = line.getP1().getY();
        }

        if(line.getP2().getY() > currentYMax){ // find maxX in p2
            currentYMax = line.getP2().getY();
        }
        if(line.getP2().getY() < currentYMin){ // find minX in p2
            currentYMin = line.getP2().getY();
        }
        if(line.getP2().getX() > currentXMax){ // find maxX in p1
            currentXMax = line.getP2().getX();
        }
        if(line.getP2().getX() < currentXMin){ // find minX in p1
            currentXMin = line.getP2().getX();
        }
    }
    //vanaf hier xmin xmax ymin ymax bepaald
    //nu ranges berekenen

    double xRange = currentXMax - currentXMin;
    double yRange = currentYMax - currentYMin;

    double imageX = size*(xRange/(std::max(xRange,yRange)));
    double imageY = size*(yRange/(std::max(xRange,yRange)));

    double d = 0.95*(imageX/xRange);

    double DCx = (d*(currentXMin + currentXMax))/2;
    double DCy = (d*(currentYMin + currentYMax))/2;

    double dx = (imageX/2) - DCx;
    double dy = (imageY/2) - DCy;

    unsigned int width = roundToInt(imageX);
    unsigned int height = roundToInt(imageY);
    img::EasyImage image(width, height);

    return {d, dx, dy, width, height};
}

void generateFractal(Figure* fig, Figures3D& fractal, const int nr_iterations, const double scale){
    fractal.push_back(fig);
    Figures3D threeDFigures2 = {};
    for(int it = 0; it<nr_iterations; ++it){
        for(Figure* loopFig:fractal){
            /*if (loopFig->type.find("Fractal") == std::string::npos) {
                continue;
            }*/
            std::vector<Vector3D> addPoints = {};
        //make the scaled version of orig figure;
            for(int p=0; p<loopFig->points.size(); ++p){
                addPoints.push_back(loopFig->points[p]*scaleFigure(1/scale));
            }
            for(int k = 1; k<=loopFig->points.size(); k++){
                std::vector<Vector3D> pointsI = {};
                int j = (k-1)%(addPoints.size());
                for(int p2 = 0; p2<addPoints.size();++p2){
                    pointsI.push_back((addPoints[p2]-addPoints[j])+loopFig->points[k-1]);
                }
                Figure* tempFig = new Figure();
                tempFig->points = pointsI;
                tempFig->faces = loopFig->faces;
                tempFig->scale = loopFig->scale;
                tempFig->color = loopFig->color;
                tempFig->type = loopFig->type;
                tempFig->center = loopFig->center;
                tempFig->rotateX = loopFig->rotateX;
                tempFig->rotateY = loopFig->rotateY;
                tempFig->rotateZ = loopFig->rotateZ;
                threeDFigures2.push_back(tempFig);
            }
        }
        fractal=threeDFigures2;
        threeDFigures2={};
    }
}




























































void createBuckyBall(Figure* tempFig){

    tempFig->addPoint(Vector3D::point( -3.4306, 0.3484, 0.3630)); //punt12
    tempFig->addPoint(Vector3D::point( -3.1790, 1.1810,-0.7334)); //punt12
    tempFig->addPoint(Vector3D::point( -2.9160, 0.3690,-1.8427)); //punt12
    tempFig->addPoint(Vector3D::point( -3.0048,-0.9660,-1.4314)); //punt12
    tempFig->addPoint(Vector3D::point( -3.3229,-0.9791,-0.0682)); //punt12
    tempFig->addPoint(Vector3D::point( -1.9449, 0.7442,-2.7732)); //punt12
    tempFig->addPoint(Vector3D::point( -1.2332, 1.9362,-2.5954)); //punt12
    tempFig->addPoint(Vector3D::point( -1.4944, 2.7485,-1.4902)); //punt12
    tempFig->addPoint(Vector3D::point( -2.4680, 2.3700,-0.5585)); //punt12
    tempFig->addPoint(Vector3D::point( -2.0114, 2.7298, 0.7143)); //punt12
    tempFig->addPoint(Vector3D::point( -2.2640, 1.9005, 1.8090)); //punt12
    tempFig->addPoint(Vector3D::point( -2.9745, 0.7074, 1.6331)); //punt12
    tempFig->addPoint(Vector3D::point( -2.4129,-0.2585, 2.4753)); //punt12
    tempFig->addPoint(Vector3D::point( -2.3071,-1.5829, 2.0466)); //punt12
    tempFig->addPoint(Vector3D::point( -2.7627,-1.9425, 0.7730)); //punt12
    tempFig->addPoint(Vector3D::point( -1.8823,-2.8972, 0.2521)); //punt12
    tempFig->addPoint(Vector3D::point( -1.5635,-2.8861,-1.1071)); //punt12
    tempFig->addPoint(Vector3D::point( -2.1254,-1.9193,-1.9490)); //punt12
    tempFig->addPoint(Vector3D::point( -1.1533,-1.5397,-2.8812)); //punt12
    tempFig->addPoint(Vector3D::point( -1.0615,-0.2088,-3.2934)); //punt12
    tempFig->addPoint(Vector3D::point(  0.0879, 1.7218,-3.0046)); //punt12
    tempFig->addPoint(Vector3D::point( -0.4365, 3.3422,-0.7917)); //punt12
    tempFig->addPoint(Vector3D::point( -0.7555, 3.3307, 0.5716)); //punt12
    tempFig->addPoint(Vector3D::point( -1.2619, 1.6696, 2.7586)); //punt12
    tempFig->addPoint(Vector3D::point( -1.3534, 0.3346, 3.1711)); //punt12
    tempFig->addPoint(Vector3D::point( -1.1426,-2.3137, 2.3104)); //punt12
    tempFig->addPoint(Vector3D::point( -0.8795,-3.1267, 1.2011)); //punt12
    tempFig->addPoint(Vector3D::point( -0.2429,-3.1010,-1.5184)); //punt12
    tempFig->addPoint(Vector3D::point(  0.0114,-2.2689,-2.6155)); //punt12
    tempFig->addPoint(Vector3D::point(  0.1936, 0.3941,-3.4362)); //punt12
    tempFig->addPoint(Vector3D::point(  3.3266, 0.9826, 0.0695)); //punt12
    tempFig->addPoint(Vector3D::point(  3.0096, 0.9706, 1.4326)); //punt12
    tempFig->addPoint(Vector3D::point(  2.9177,-0.3640, 1.8443)); //punt12
    tempFig->addPoint(Vector3D::point(  3.1778,-1.1770, 0.7349)); //punt12
    tempFig->addPoint(Vector3D::point(  3.4306,-0.3451,-0.3621)); //punt12
    tempFig->addPoint(Vector3D::point(  1.9473,-0.7475, 2.7721)); //punt12
    tempFig->addPoint(Vector3D::point(  1.0645, 0.2062, 3.2916)); //punt12
    tempFig->addPoint(Vector3D::point(  1.1547, 1.5385, 2.8834)); //punt12
    tempFig->addPoint(Vector3D::point(  2.1280, 1.9203, 1.9528)); //punt12
    tempFig->addPoint(Vector3D::point(  1.5657, 2.8862, 1.1107)); //punt12
    tempFig->addPoint(Vector3D::point(  1.8834, 2.9000,-0.2489)); //punt12
    tempFig->addPoint(Vector3D::point(  2.7651, 1.9463,-0.7707)); //punt12
    tempFig->addPoint(Vector3D::point(  2.3092, 1.5866,-2.0437)); //punt12
    tempFig->addPoint(Vector3D::point(  2.4144, 0.2631,-2.4755)); //punt12
    tempFig->addPoint(Vector3D::point(  2.9758,-0.7031,-1.6328)); //punt12
    tempFig->addPoint(Vector3D::point(  2.2666,-1.8962,-1.8100)); //punt12
    tempFig->addPoint(Vector3D::point(  2.0133,-2.7276,-0.7172)); //punt12
    tempFig->addPoint(Vector3D::point(  2.4693,-2.3670, 0.5560)); //punt12
    tempFig->addPoint(Vector3D::point(  1.4973,-2.7481, 1.4877)); //punt12
    tempFig->addPoint(Vector3D::point(  1.2349,-1.9394, 2.5952)); //punt12
    tempFig->addPoint(Vector3D::point( -0.1917,-0.3936, 3.4362)); //punt12
    tempFig->addPoint(Vector3D::point( -0.0088, 2.2691, 2.6153)); //punt12
    tempFig->addPoint(Vector3D::point(  0.2447, 3.1026, 1.5191)); //punt12
    tempFig->addPoint(Vector3D::point(  0.8815, 3.1269,-1.2000)); //punt12
    tempFig->addPoint(Vector3D::point(  1.1442, 2.3150,-2.3100)); //punt12
    tempFig->addPoint(Vector3D::point(  1.3558,-0.3333,-3.1706)); //punt12
    tempFig->addPoint(Vector3D::point(  1.2640,-1.6687,-2.7595)); //punt12
    tempFig->addPoint(Vector3D::point(  0.7577,-3.3290,-0.5712)); //punt12
    tempFig->addPoint(Vector3D::point(  0.4382,-3.3422, 0.7920)); //punt12
    tempFig->addPoint(Vector3D::point( -0.0856,-1.7214, 3.0049)); //punt12

    //now add the points to be matched
    std::vector<std::vector<int>> pointIndexVector;
    for(int i = 0; i < 60; i++){
        pointIndexVector.push_back({i,i+1,i+2});
    }



    for(const auto & i : pointIndexVector){
        Face* tempFace = new Face(i);
        tempFig->faces.push_back(tempFace);
    }
}

#endif //UTILS_EXTRAFUNCTIONS_H

