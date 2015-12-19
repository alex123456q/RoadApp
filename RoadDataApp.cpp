// RoadDataApp.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "RoadData.h"
#include "RoadDataLoader.h"
#include <dlib/gui_widgets.h>
#include <dlib/image_transforms.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>

struct myPoint{
    RoadPoint p;
    int numb_clast;
    myPoint(RoadPoint* _p, int _numb_clast):
    p(*_p), numb_clast(_numb_clast)
    {}
};
/* INTERSECTION
     createjs.Rectangle.prototype.intersects = function(rect){
                return (this.x <= rect.x + rect.width &&
                        rect.x <= this.x + this.width &&
                        this.y <= rect.y + rect.height &&
                        rect.y <= this.y + this.height);
            }

ANOTHER
bool DoBoxesIntersect(Box a, Box b) {
  return (abs(a.x - b.x) * 2 < (a.width + b.width)) &&
         (abs(a.y - b.y) * 2 < (a.height + b.height));
}
function IntersectRect(r1:Rectangle, r2:Rectangle):Boolean {
    return !(r2.left > r1.right
        || r2.right < r1.left
        || r2.top > r1.bottom
        || r2.bottom < r1.top);
}
*/
/* BOUNDING BOX
auto xExtremes = std::minmax_element(v.begin(), v.end(),
                                     [](const ofPoint& lhs, const ofPoint& rhs) {
                                        return lhs.x < rhs.x;
                                     });

auto yExtremes = std::minmax_element(v.begin(), v.end(),
                                     [](const ofPoint& lhs, const ofPoint& rhs) {
                                        return lhs.y < rhs.y;
                                     });

ofPoint upperLeft(xExtremes.first->x, yExtremes.first->y);
ofPoint lowerRight(xExtremes.second->x, yExtremes.second->y);

#include <boost/algorithm/minmax_element.hpp>
bool compareX(ofPoint lhs, ofPoint rhs) { return lhs.x < rhs.x; };
bool compareY(ofPoint lhs, ofPoint rhs) { return lhs.y < rhs.y; };

// ....
    pair<vector<ofPoint>::iterator, vector<ofPoint>::iterator> xExtremes, yExtremes;
    xExtremes = boost::minmax_element(overlap_point.begin(), overlap_point.end(), compareX);
    yExtremes = boost::minmax_element(overlap_point.begin(), overlap_point.end(), compareY);
    ofPoint upperLeft(xExtremes.first->x, yExtremes.first->y);
    ofPoint lowerRight(xExtremes.second->x, yExtremes.second->y);
*/
#define R 30     //ширина поиска локальных сгущений - входной параметр алгоритма
/*Clustering_not_finish();      //все ли объекты кластеризованы
Rand_object();      //возвращает произвольный некластеризованный объект
Generate_same_object(type *object);     //возвращает массив объектов, расположенных на расстоянии <= R от текущего
Center_of_objects(type *mas_of_objects);      //возвращает центр тяжести указанных объектов
Delete_objects(type *mas_of_objects);    //удаляет указанные объекты из выборки (мы их уже кластеризовали)
 
while(Clustering_not_finish())
{
   currently_object = Rand_object();
   mas_of_same_objects = Generate_same_object(currently_object); 
   center_object = Center_of_objects(mas_of_same_objects);
 
   while (center_object != currently_object)  //пока центр тяжести не стабилизируется
   {
      currently_object = center_object;
      mas_of_same_objects = Generate_same_object(currently_object);
      center_object = Center_of_objects(mas_of_same_objects);
   } 
   Delete_objects(mas_of_same_objects);
}*/

double squaredDistance(myPoint& a, myPoint& b){
    return pow(a.p.x - b.p.x, 2) + pow(a.p.y-b.p.y, 2) + pow(a.p.z-b.p.z, 2);
};

myPoint* Center_of_objects(std::vector<myPoint>& points){
    double mindist = 1000000.0;
    double curdist;
    myPoint* minpoint;
    minpoint = NULL;
    for (int i = 0; i < points.size(); ++i){
        curdist = 0;
        for (int j = i + 1; j < points.size(); ++j)
            curdist += squaredDistance(points[i], points[j]);
        if (mindist < curdist){
            minpoint = &points[i];
        }
    }
};

void FOREL(std::vector<myPoint>& allpoints){
    std::vector<int> indexes(allpoints.size());
    int curind = rand()%indexes.size();
    for (int i = 0; i < indexes.size(); ++i)
        if (squaredDistance() < R)
             curclast.push_back();
    while (curobj != predobj){
    center_of_objects = ;
    }

};

int _tmain(int argc, _TCHAR* argv[])
{
    if (argc < 2) return -1;
    OpenData(argv[1]);

    //Количество блоков в las файле
    int nBlock = GetBlocksCount();
    int bBLockSize = GetBlockSize();
    //Получить метаточки описывающие блок
    std::vector<RoadPoint> points;
    points.resize(nBlock);
    GetMetaPoints(points.data(), nBlock);
    points.resize(bBLockSize);
    BlockInfoHeader pHeader;


    std::ofstream fout("New1.txt"); // создаём объект класса ofstream для записи и связываем его с файлом cppstudio.txt

    std::vector<dlib::perspective_window::overlay_dot> pointsPaint;
    
    double eps = 10e-5;
    double minZ = 22.6;
    double maxZ = 0;
    std::vector<myPoint> allpoints;
    std::vector<std::vector<myPoint*> > heights;
    
    for (int i = 0; i < 5/*nBlock*/; i++)
    {
        GetCloudPoints(i, points.data(), &pHeader);
        fout << "Frame #" << i << std::endl;
        //fout << " X: " << pHeader.FrameXCoord << std::setprecision(5);
        //fout << " Y: " << pHeader.FrameYCoord << std::setprecision(5);
        //fout << " Z: " << pHeader.FrameZCoord << std::setprecision(5) << std::endl;
        for (int j = 0; j < points.size(); ++j){
            if (points.data()[j].x <eps && points.data()[j].y < eps && points.data()[j].z < eps)
                continue;
            //minZ = min(minZ, points.data()[j].z);
            //maxZ = max(maxZ, points.data()[j].z);
            allpoints.push_back(myPoint(&points.data()[j], (int)floor(points.data()[j].z - minZ)));
            if (points.data()[j].z - minZ > heights.size())
                heights.push_back(std::vector<myPoint*>());
            heights[(int)floor(points.data()[j].z - minZ)].push_back(&allpoints[allpoints.size()-1]);
            //for visualization
            dlib::vector<double> val(points.data()[j].x, points.data()[j].y, points.data()[j].z);
            // Pick a color based on how far we are along the spiral
            dlib::rgb_pixel color = dlib::colormap_jet(1,0,20);
            // And add the point to the list of points we will display
            pointsPaint.push_back(dlib::perspective_window::overlay_dot(val, color));
        }

    }
    std::cout << minZ << maxZ << std::endl;

    fout.close();
    CloseData();

    //for visualization
    dlib::perspective_window win;
    win.set_title("perspective_window 3D point cloud");
    win.add_overlay(pointsPaint);
    win.wait_until_closed();


	return 0;
}

