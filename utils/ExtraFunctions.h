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
    //std::cout << "------------------------BEGIN SCALING---------------------- \n";
    Matrix scaleMatrix;
    scaleMatrix(1,1) = scale;
    scaleMatrix(2,2) = scale;
    scaleMatrix(3,3) = scale;
    //for(int i = 0; i < points.size(); i++){
    //    std::cout << "old point: " << points[i].x << ", " << points[i].y << ", " << points[i].z;
    //    points[i] = points[i]*scaleMatrix;
    //    std::cout << "  ------>" << points[i].x << ", " << points[i].y << ", " << points[i].z << std::endl;
    //}
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

Matrix eyePointTrans(const Vector3D &eyepoint){
    double theta;
    double phi;
    double r;
    toPolar(eyepoint, theta, phi, r);
    //Matrix V =  rotateZ((M_PI/2)-theta, "radians")*rotateX(phi, "radians")*translate(Vector3D::vector(0,0,r));
    Matrix V =  rotateZ((-M_PI/2)-theta, "radians")*rotateX(-phi, "radians")*translate(Vector3D::vector(0,0,-r));

    return V;
}
void applyTransformation(Figure &figure){
    Matrix transformatieMatrix;
    transformatieMatrix = scaleFigure(figure.scale)*rotateX(figure.rotateX)*rotateY(figure.rotateY)*rotateZ(figure.rotateZ)*translate(figure.center);
    for(unsigned int i = 0; i < figure.points.size();i++){
        figure.points[i] = figure.points[i]*transformatieMatrix;
    }
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
Lines2D doProjection(const Figures3D& figures3D, const Vector3D &eye){
    Matrix V = eyePointTrans(eye);
    Lines2D lines2D;
    for(const auto& figure : figures3D){
        //iterate over every figure points and multiply V (eyepoint transformation to each point)
        for(unsigned int i = 0; i < figure->points.size(); i++){
            figure->points[i] = figure->points[i]*V;
        }
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
    //now add the points to be matched
    std::vector<std::vector<int>> pointIndexVector;
    pointIndexVector.push_back({0,1,2});
    pointIndexVector.push_back({1,3,2});
    pointIndexVector.push_back({0,3,1});
    pointIndexVector.push_back({0,2,3});
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
        double u = 2*i*M_PI/n;
        for(int j = 0; j < m; j++){ //loop over number of points in circle
            double v = 2*j*M_PI/m;
            double xuv = (R+r*cos(v))*cos(u);  //met u,v E [0,2pi]
            double yuv = (R+r*cos(v))*sin(u);
            double zuv = r*sin(v);
            tempFig->addPoint(Vector3D::point(xuv,yuv,zuv));
        }
        if(i>=1){
            int index = tempFig->points.size()-1;
            pointIndexVector.push_back({index, index-m+1});
            for(int k = 0; k < m; k++){ //for every point in the circle
                if(k==0){
                    continue;
                }
                if(i == n||i==n-1){//this will be the last circle next indices will restart to 0 again
                    pointIndexVector.push_back({index,index-1, (index-m-1)%n, (index-m)%n});
                }
                pointIndexVector.push_back({index,index-1, index-m-1,index-m});
                index-=1;
            }
        }
    }
    pointIndexVector.push_back({0, m-1}); //1st point and last point of first circle;

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
#endif //UTILS_EXTRAFUNCTIONS_H
