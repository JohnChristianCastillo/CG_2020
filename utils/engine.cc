#include "easy_image.h"
#include "ini_configuration.h"
#include "Line2D.h"
#include <fstream>
#include <list>
#include <stdexcept>
#include <string>
#include "l_parser.h"
#include "Figure.h"
#include "ExtraFunctions.h"
#include <stack>
#include <string>
// student@lab33:~/Desktop/utils/intro$ ../engine Inro1_ColorRectangle.ini  from photo folers
// student@lab33:~/Desktop/utils$ make   from root folder

//scaling 3D matrix  https://www.youtube.com/watch?v=HEydT6_Re84
typedef std::list<Line2D> Lines2D;
typedef std::list<Figure*> Figures3D;



Figures3D threeDFigures;




img::EasyImage ColorRectangle(unsigned int width, unsigned int height){
    img::EasyImage image(width, height);
    for(unsigned int i = 0; i < width; i++)
    {
        for(unsigned int j = 0; j < height; j++)
        {
            image(i,j).red = roundToInt(i*255/width);             //round to int is in "ExtraFunctions.h" which is included in color
            image(i,j).green = roundToInt(j*255/height);
            image(i,j).blue = roundToInt((i*255)/width+(j*255)/height)%256;
        }
    }
    return image;
}




img::EasyImage draw2DLines(const Lines2D &lines, const int size = 50, Color color = Color(0,0,0), bool zBuff = false){
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
    for(unsigned int i = 0; i < width; i++)
    {
        for(unsigned int j = 0; j < height; j++)
        {
            image(i,j).red = roundToInt(color.getRed()*255);             //round to int is in "ExtraFunctions.h" which is included in color
            image(i,j).green = roundToInt(color.getGreen()*255);
            image(i,j).blue = roundToInt(color.getBlue()*255);
        }
    }
    if(zBuff){
        ZBuffer zBuffer = ZBuffer(width, height);
        for(const Line2D& line : lines){
            unsigned int p1X = roundToInt(line.getP1().getX()*d+dx);
            unsigned int p1Y = roundToInt(line.getP1().getY()*d+dy);
            unsigned int p2X = roundToInt(line.getP2().getX()*d+dx);
            unsigned int p2Y = roundToInt(line.getP2().getY()*d+dy);
            image.draw_zbuf_line(zBuffer, p1X, p1Y, line.z1, p2X, p2Y, line.z2, line.getColor());
        }
    }
    else if(!zBuff){
	    for(const Line2D& line : lines){
	    	image.draw_line(roundToInt(line.getP1().getX()*d+dx), roundToInt(line.getP1().getY()*d+dy),
	    			roundToInt(line.getP2().getX()*d+dx), roundToInt(line.getP2().getY()*d+dy), line.getColor());

	    }
    }

	return image;

	//use draw line vna easyimage  image.draw_line
}

Lines2D drawLSystem(const LParser::LSystem2D &l_system, Color& color){
    double angle = l_system.get_starting_angle();
    double x1 = 0;
    double x2 = 0;
    double y1 = 0;
    double y2 = 0;
    bool twoDLineExists = false;
    std::stack<double> xStack;
    std::stack<double> yStack;
    std::stack<double> hoekStack;
    Lines2D lines2 = {};
    std::string currentString = l_system.get_initiator();


    for(unsigned int i = 0; i<l_system.get_nr_iterations(); i++){ //iterate over the number of provided iterations
    //for(int i = 0; i < 1; i++){
        std::string tempString;
        for (char j : currentString) {while(angle < 0){
                angle += 360;
            }
            while(angle >= 360){ //make sure angle is 0 <= angle <360
                angle -= 360;
            }
            //auto temp = j;
            const bool charInAlphabet = l_system.get_alphabet().find(j) != l_system.get_alphabet().end();
            if(charInAlphabet){
                tempString += l_system.get_replacement(j);// replace the current char
            }
            else{
                tempString += j;
            }
        }

        currentString = tempString; //current string is now the current string in which the replacement rules has been used
        tempString = "";
    }
    /*
    std::vector<char> open;
    std::vector<char> close;
    for(char c: currentString){
        if(c == '('){
            open.push_back(c);
        }
        else if(c == ')'){
            close.push_back(c);
        }
    }

    int openSize = open.size();
    int closeSize = close.size();
    */
    for(char c: currentString){
        if(c == '-'){
            angle -= l_system.get_angle();
        }
        else if(c == '+'){
            angle += l_system.get_angle();
        }
        else if(c == '('){
            xStack.push(x1);
            yStack.push(y1);
            hoekStack.push(angle);
        }
        else if(c == ')'){
            x1 = xStack.top();
            y1 = yStack.top();
            angle = hoekStack.top();
            xStack.pop();
            yStack.pop();
            hoekStack.pop();
        }
        else if(l_system.draw(c)){
            //convert angle to radians

            double radianAngle = (angle*M_PI)/180;
            //calculate x2, y2
            x2 = x1+cos(radianAngle);
            y2 = y1+sin(radianAngle);

            twoDLineExists = true;
            }
        if(twoDLineExists){
            Line2D twoDLine(Point2D(x1,y1), Point2D(x2,y2), color);
            lines2.push_back(twoDLine);
            //now set the x1 to x2 after creating and saving the object
            x1 = x2;
            y1 = y2;
            twoDLineExists = false;
        }
    }
    return lines2;
}
void draw3DLSystem(const LParser::LSystem3D &l_system, Figure* tempFig);

img::EasyImage generate_image(const ini::Configuration &configuration)
{
    std::string type = configuration["General"]["type"].as_string_or_die();

    if(type == "IntroColorRectangle"){
        int width = configuration["ImageProperties"]["width"].as_int_or_die();
        int height = configuration["ImageProperties"]["height"].as_int_or_die();
        return ColorRectangle(width, height);
    }
    else if(type == "2dLines"){
        Lines2D lines;
    	return draw2DLines(lines , 50);
    }
    else if(type == "2DLSystem"){
        LParser::LSystem2D l_system;
        std::cout << configuration["2DLSystem"]["inputfile"].as_string_or_die() << std::endl;
        std::cout << configuration["General"]["size"].as_int_or_die();
        std::ifstream input_stream(configuration["2DLSystem"]["inputfile"].as_string_or_die());

        std::vector<double> backgroundColor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        std::vector<double> lineColor = configuration["2DLSystem"]["color"].as_double_tuple_or_die();

        Color lineColorObject = Color(lineColor[0], lineColor[1], lineColor[2]);
        Color backgroundColorObject = Color(backgroundColor[0], backgroundColor[1], backgroundColor[2]);

        input_stream >> l_system;
        input_stream.close();

        return draw2DLines(drawLSystem(l_system, lineColorObject), configuration["General"]["size"], backgroundColorObject);

    }
    else if(type == "Wireframe" || type == "ZBufferedWireframe"){
        bool zBuff = false;
        if(type == "ZBufferedWireframe"){zBuff = true;}
        std::vector<double> bgColorVect = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        std::vector<double> eyeVector = configuration["General"]["eye"].as_double_tuple_or_die();

        Vector3D eye = Vector3D::vector(eyeVector[0], eyeVector[1], eyeVector[2]);
        Color bgColor = Color(bgColorVect[0], bgColorVect[1], bgColorVect[2]);
        int nrFigures = configuration["General"]["nrFigures"].as_int_or_die(); //needs to positive whole number
        std::string figureType = configuration["General"]["type"].as_string_or_die();
        for(int i = 0; i<nrFigures; i++){
            std::string figure = "Figure"+std::to_string(i);
            Figure* tempFig = new Figure();
            tempFig->setType(configuration[figure]["type"].as_string_or_die());
            //save Center figure object
            std::vector<double> tempCenter = configuration[figure]["center"].as_double_tuple_or_die();
            tempFig->center = Vector3D::point(tempCenter[0],tempCenter[1],tempCenter[2]);
            //save Color figure object
            std::vector<double> figureColor = configuration[figure]["color"].as_double_tuple_or_die();
            Color figureColorObject = Color(figureColor[0], figureColor[1], figureColor[2]);
            tempFig->setColor(figureColorObject);
            //save Rotates
            tempFig->rotateX = configuration[figure]["rotateX"];
            tempFig->rotateY = configuration[figure]["rotateY"];
            tempFig->rotateZ = configuration[figure]["rotateZ"];
            tempFig->scale = configuration[figure]["scale"];
            //save Points to figure object
            if(tempFig->type == "LineDrawing"){
                for(int j = 0; j < configuration[figure]["nrPoints"].as_int_or_die(); j++){
                    std::string point = "point"+std::to_string(j);
                    std::vector<double> pointVec = configuration[figure][point].as_double_tuple_or_die();
                    Vector3D tempPoint = Vector3D::point(pointVec[0], pointVec[1], pointVec[2]);
                    tempFig->addPoint(tempPoint);
                }
                //save Faces to figure object
                for(int j = 0; j < configuration[figure]["nrLines"].as_int_or_die(); j++){
                    std::string line = "line"+std::to_string(j);
                    std::vector<int> lineVec = configuration[figure][line].as_int_tuple_or_die();
                    Face* tempFace = new Face;
                    tempFace->setPointIndexes(lineVec);
                    tempFig->addFace(tempFace);
                }
            }
            else if(tempFig->type == "Cube"){
                createCube(tempFig);
            }
            else if(tempFig->type == "Tetrahedron"){
                createTetrahedron(tempFig);
            }
            else if(tempFig->type == "Octahedron"){
                createOctahedron(tempFig);
            }
            else if(tempFig->type == "Icosahedron"){
                createIcosahedron(tempFig);
            }
            else if(tempFig->type == "Dodecahedron"){
                createDodecahedron(tempFig);
            }
            else if(tempFig->type == "Sphere"){
                Figure* icosahedron = new Figure;
                createIcosahedron(icosahedron);
                tempFig->points = icosahedron->points;
                tempFig->faces = icosahedron->faces;
                createSphere(tempFig, configuration[figure]["n"].as_int_or_die());
            }
            else if(tempFig->type == "Cone"){
                createCone(configuration[figure]["n"].as_int_or_die(), configuration[figure]["height"].as_double_or_die(), tempFig);
            }
            else if(tempFig->type == "Cylinder"){
                createCylinder(configuration[figure]["n"].as_int_or_die(), configuration[figure]["height"].as_double_or_die(), tempFig);
            }
            else if(tempFig->type == "Torus"){
                createTorus(configuration[figure]["r"], configuration[figure]["R"], configuration[figure]["n"], configuration[figure]["m"], tempFig);
            }
            else if(tempFig->type == "3DLSystem"){
                LParser::LSystem2D l_system;
                LParser::LSystem3D l_system3D;
                std::cout << configuration[figure]["inputfile"].as_string_or_die() << std::endl;
                std::cout << configuration["General"]["size"].as_int_or_die();
                std::ifstream input_stream(configuration[figure]["inputfile"].as_string_or_die());

                std::vector<double> backgroundColor = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
                std::vector<double> lineColor = configuration[figure]["color"].as_double_tuple_or_die();


                input_stream >> l_system3D;
                input_stream.close();


                draw3DLSystem(l_system3D, tempFig);
            }
            //add the parsed figure object to 3DFigures
            threeDFigures.push_back(tempFig);
        }
        //now that it's parsed we can actually apply the transformations and project it
        applyTransformation(threeDFigures); //calls "applyTransforation() to every 3Dfigures member
        Lines2D wireLines = doProjection(threeDFigures, eye); //applies the eyepointTransformation to every point of 3DFigures members
                                                             //and generates the lines that are to be drawn
        threeDFigures = {};
        return draw2DLines(wireLines, configuration["General"]["size"], bgColor, zBuff);
        //return draw2DLines(wireLines, 2000, bgColor);
    }

	return img::EasyImage();
}
void draw3DLSystem(const LParser::LSystem3D &l_system, Figure* tempFig){
    int point1 = 0;
    double sigmaRad = (l_system.get_angle()*M_PI)/180;
    Vector3D xyz = Vector3D::point(0,0,0);
    tempFig->addPoint(xyz);
    Vector3D H = Vector3D::vector(1,0,0); //x-as
    Vector3D L = Vector3D::vector(0,1,0); //y-as //links
    Vector3D U = Vector3D::vector(0,0,1); //z-as //opwaarts
    std::stack<Vector3D> xyzStack;
    std::stack<Vector3D> Hstack;
    std::stack<Vector3D> Lstack;
    std::stack<Vector3D> Ustack;
    std::stack<int> point1Stack;

    Lines2D lines2 = {};
    std::string currentString = l_system.get_initiator();


    for(unsigned int i = 0; i<l_system.get_nr_iterations(); i++){ //iterate over the number of provided iterations
        std::string tempString;
        for (char j : currentString) {
            const bool charInAlphabet = l_system.get_alphabet().find(j) != l_system.get_alphabet().end();
            if(charInAlphabet){
                tempString += l_system.get_replacement(j);// replace the current char
            }
            else{
                tempString += j;
            }
        }

        currentString = tempString; //current string is now the current string in which the replacement rules has been used
        tempString = "";
    }
    /*
    std::vector<char> open;
    std::vector<char> close;
    for(char c: currentString){
        if(c == '('){
            open.push_back(c);
        }
        else if(c == ')'){
            close.push_back(c);
        }
    }

    int openSize = open.size();
    int closeSize = close.size();
    */
    std::vector<std::vector<int>> pointIndexVector;

    for(char c: currentString){
        if(c == '-'){
           //H = H*rotateZ(l_system.get_angle());
           //L = L*rotateZ(l_system.get_angle());
            H = Vector3D::vector(H.x,H.y,H.z)*cos(sigmaRad) + Vector3D::vector(L.x,L.y,L.z)*sin(sigmaRad);
            L = -H*sin(sigmaRad) + L*cos(sigmaRad);
            //L = -Vector3D::vector(H.x,H.y,H.z)*sin(sigmaRad) + Vector3D::vector(L.x,L.y,L.z)*cos(sigmaRad);
        }

        else if(c == '+'){
            //auto h = Vector3D::point(H.x,H.y,H.z)*cos((angle*M_PI)/180) + Vector3D::point(L.x,L.y,L.z)*sin((angle*M_PI)/180);
             //H = H*rotateZ(-l_system.get_angle());
             //L = L*rotateZ(-l_system.get_angle());

           H = Vector3D::vector(H.x,H.y,H.z)*cos(-sigmaRad) + Vector3D::vector(L.x,L.y,L.z)*sin(-sigmaRad);
           L = -Vector3D::vector(H.x,H.y,H.z)*sin(-sigmaRad) + Vector3D::vector(L.x,L.y,L.z)*cos(-sigmaRad);

        }
        else if(c == '^'){
            //auto h = Vector3D::point(H.x,H.y,H.z)*cos((angle*M_PI)/180) + Vector3D::point(U.x,U.y,U.z)*sin((angle*M_PI)/180);
            //H = H*rotateY(l_system.get_angle());
            //U = U*rotateY(l_system.get_angle());
           H = Vector3D::vector(H.x,H.y,H.z)*cos(sigmaRad) + Vector3D::vector(U.x,U.y,U.z)*sin(sigmaRad);
           U = -Vector3D::vector(H.x,H.y,H.z)*sin(sigmaRad) + Vector3D::vector(U.x,U.y,U.z)*cos(sigmaRad);

        }
        else if(c == '&'){
          //H = H*rotateY(-l_system.get_angle());
          //U = U*rotateY(-l_system.get_angle());
            H = Vector3D::vector(H.x,H.y,H.z)*cos(-sigmaRad) + Vector3D::vector(U.x,U.y,U.z)*sin(-sigmaRad);
            U = -Vector3D::vector(H.x,H.y,H.z)*sin(-sigmaRad) + Vector3D::vector(U.x,U.y,U.z)*cos(-sigmaRad);

        }
        else if(c == '\\'){
           //U = U*rotateX(l_system.get_angle());
           //L = L*rotateX(l_system.get_angle());
           L = Vector3D::vector(L.x,L.y,L.z)*cos(sigmaRad) - Vector3D::vector(U.x,U.y,U.z)*sin(sigmaRad);
           U = Vector3D::vector(L.x,L.y,L.z)*sin(sigmaRad) + Vector3D::vector(U.x,U.y,U.z)*cos(sigmaRad);

        }
        else if(c == '/'){
          //U = U*rotateX(-l_system.get_angle());
          //L = L*rotateX(-l_system.get_angle());
           L = Vector3D::vector(L.x,L.y,L.z)*cos(-sigmaRad) - Vector3D::vector(U.x,U.y,U.z)*sin(-sigmaRad);
           U = Vector3D::vector(L.x,L.y,L.z)*sin(-sigmaRad) + Vector3D::vector(U.x,U.y,U.z)*cos(-sigmaRad);
        }
        else if(c == '|'){
            //H = H*rotateZ(180);
            //L = L*rotateZ(180);
            H = -H;
            L = -L;
        }
        else if(c == '('){
            xyzStack.push(xyz);
            Hstack.push(H);
            Ustack.push(U);
            Lstack.push(L);
            point1Stack.push(point1);
        }
        else if(c == ')'){
            xyz = xyzStack.top();
            H = Hstack.top();
            U = Ustack.top();
            L = Lstack.top();
            point1 = point1Stack.top();

            xyzStack.pop();
            Hstack.pop();
            Ustack.pop();
            Lstack.pop();
            point1Stack.pop();
        }
        else if(l_system.draw(c)){
            xyz = Vector3D::point(xyz.x+H.x, xyz.y+H.y, xyz.z+H.z);
            tempFig->addPoint(xyz);
            if(tempFig->points.size() == 3){
                //pointIndexVector.push_back({1,4});
                //pointIndexVector.push_back({0,3});
                break;
            }
            pointIndexVector.push_back({point1,int(tempFig->points.size()-1)});
            point1 = tempFig->points.size()-1;
        }
    }
    //now add the points to be matched

    for(const auto &i: pointIndexVector){
        Face* tempFace = new Face(i);
        tempFig->faces.push_back(tempFace);
    }
}

int main(int argc, char const* argv[])
{
        int retVal = 0;
        try
        {
                for(int i = 1; i < argc; ++i)
                {
                        ini::Configuration conf;
                        try
                        {
                                std::ifstream fin(argv[i]);
                                fin >> conf;
                                fin.close();
                        }
                        catch(ini::ParseException& ex)
                        {
                                std::cerr << "Error parsing file: " << argv[i] << ": " << ex.what() << std::endl;
                                retVal = 1;
                                continue;
                        }

                        img::EasyImage image = generate_image(conf);
                        if(image.get_height() > 0 && image.get_width() > 0)
                        {
                                std::string fileName(argv[i]);
                                std::string::size_type pos = fileName.rfind('.');
                                if(pos == std::string::npos)
                                {
                                        //filename does not contain a '.' --> append a '.bmp' suffix
                                        fileName += ".bmp";
                                }
                                else
                                {
                                        fileName = fileName.substr(0,pos) + ".bmp";
                                }
                                try
                                {
                                        std::ofstream f_out(fileName.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
                                        f_out << image;

                                }
                                catch(std::exception& ex)
                                {
                                        std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                                        retVal = 1;
                                }
                        }
                        else
                        {
                                std::cout << "Could not generate image for " << argv[i] << std::endl;
                        }
                }
        }
        catch(const std::bad_alloc &exception)
        {
    		//When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
    		//Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
    		//(Unless of course you are already consuming the maximum allowed amount of memory)
    		//If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
		//mark the test as failed while in reality it just needed a bit more memory
                std::cerr << "Error: insufficient memory" << std::endl;
                retVal = 100;
        }
        return retVal;
}
