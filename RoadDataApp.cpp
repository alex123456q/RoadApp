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
#include <map>

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

double squaredDistance(myPoint& a, myPoint& b){
    return pow(a.p.x - b.p.x, 2) + pow(a.p.y-b.p.y, 2) + pow(a.p.z-b.p.z, 2);
};

int/*myPoint* */ CenterObject(std::vector<myPoint>& points, std::vector<int>& indexes){
    double mindist = 1000000.0;
    double curdist;
    int minindex = 0;
    myPoint* minpoint;
    minpoint = NULL;
    for (int i = 0; i < indexes.size(); ++i){
        curdist = 0.0;
        for (int j = i + 1; j < indexes.size(); ++j)
            curdist += squaredDistance(points[indexes[i]], points[indexes[j]]);
        if (mindist > curdist){
            minpoint = &points[indexes[i]];
            mindist = curdist;
            minindex = indexes[i];
        }
    }
    return minindex;//minpoint;
};

void genClaster(std::vector<myPoint>& points, std::vector<int>& indexes, int center, std::vector<int>* curclast){
    for (int i = 0; i < indexes.size(); ++i){
        if (squaredDistance(points[indexes[i]], points[center]) < R*R)
            curclast->push_back(indexes[i]);
    }
}

int FOREL(std::vector<myPoint>& allpoints){
    std::vector<int> indexes(allpoints.size());
    for (int i = 0; i < indexes.size(); ++i)
        indexes[i] = i;
    std::vector<int> curclast;
    int predcenter, curcenter;
    int numClast = 0;
    while (indexes.size() > 0){
        predcenter = rand()%indexes.size();
        genClaster(allpoints, indexes, predcenter, &curclast);
        curcenter = CenterObject(allpoints, curclast);
        while (curcenter != predcenter){
            predcenter = curcenter;
            curclast.clear();
            genClaster(allpoints, indexes, predcenter, &curclast);
            curcenter = CenterObject(allpoints, curclast);
        }
        std::vector<int>::iterator it = indexes.begin();
        for (int i = 0; i <curclast.size(); ++i){
            indexes.erase(it + curclast[i] - i);
            allpoints[curclast[i]].numb_clast= numClast;                  //+height*10
        }
        ++numClast;
    }
    return numClast;
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
    //std::map<int, std::vector<myPoint*> > heights;
    std::vector<int, std::vector<myPoint*> > heights;
    for (int i = 0; i < 5/*nBlock*/; i++)
    {
        GetCloudPoints(i, points.data(), &pHeader);
        fout << "Frame #" << i << std::endl;
        fout << " X: " << pHeader.FrameXCoord << std::setprecision(10);
        fout << " Y: " << pHeader.FrameYCoord << std::setprecision(10);
        fout << " Z: " << pHeader.FrameZCoord << std::setprecision(10) << std::endl;
        for (int j = 0; j < points.size(); ++j){
            if (points.data()[j].x <eps && points.data()[j].y < eps && points.data()[j].z < eps)
                continue;
            //minZ = min(minZ, points.data()[j].z);
            //maxZ = max(maxZ, points.data()[j].z);
            allpoints.push_back(myPoint(&points.data()[j], (int)floor(points.data()[j].z - minZ)));
            //if (points.data()[j].z - minZ > heights.size())
            //    heights.push_back(std::vector<myPoint*>());
            heights[(int)floor((points.data()[j].z - minZ))].push_back(&allpoints[allpoints.size()-1]);
        }

    }

    for (std::map<int, std::vector<myPoint*> >::iterator it = heights.begin(); it != heights.end(); ++it)
        std::cout << FOREL(heights[(it->first)]);
        for  (int j = 0; j < heigths.size(); ++j){
            dlib::vector<double> val(allpoints[j].p.x, points.data()[j].y, points.data()[j].z);
            // Pick a color based on how far we are along the spiral
            dlib::rgb_pixel color = dlib::colormap_jet(,0,20);
            // And add the point to the list of points we will display
            pointsPaint.push_back(dlib::perspective_window::overlay_dot(val, color));
        }
    /*
     //for visualization
            dlib::vector<double> val(points.data()[j].x, points.data()[j].y, points.data()[j].z);
            // Pick a color based on how far we are along the spiral
            dlib::rgb_pixel color = dlib::colormap_jet(1,0,20);
            // And add the point to the list of points we will display
            pointsPaint.push_back(dlib::perspective_window::overlay_dot(val, color));*/
    //int FOREL(std::vector<myPoint>& allpoints)
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

