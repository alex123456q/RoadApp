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
#include <cmath>

struct myPoint{
    RoadPoint p;
    int numb_clast;
    myPoint(RoadPoint* _p, int _numb_clast = -1):
    numb_clast(_numb_clast)
    {
        p.x = _p->x;
        p.y = _p->y;
        p.z = _p->z;
    }
    myPoint(double x, double y, double z, int _numb_clast):
    numb_clast(_numb_clast)
    {
        p.x = x;
        p.y = y;
        p.z = z;
    }
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
//BOUNDING BOX
#define A 0.49
#define H 1
#define kap 0.5
void createBoundingBox(std::vector<myPoint*>& points){
    auto xExtremes = std::minmax_element(points.begin(), points.end(),
                                         [](const myPoint* lhs, const myPoint* rhs) {
                                             return lhs->p.x < lhs->p.x;
                                         });

    auto yExtremes = std::minmax_element(points.begin(), points.end(),
                                         [](const myPoint* lhs, const myPoint* rhs) {
                                             return lhs->p.y < lhs->p.y;
                                         });
    auto zExtremes = std::minmax_element(points.begin(), points.end(),
                                         [](const myPoint* lhs, const myPoint* rhs) {
                                             return lhs->p.z < lhs->p.z;
                                         });

    myPoint Left((*xExtremes.first)->p.x, (*yExtremes.first)->p.y, (*zExtremes.first)->p.z, points[0]->numb_clast);
    myPoint Right((*xExtremes.second)->p.x, (*yExtremes.second)->p.y, (*zExtremes.second)->p.z, points[0]->numb_clast);
}

bool checkBox(myPoint Left, myPoint Right){
    double Lx = fabs(Right.p.x - Left.p.x);
    double Ly = fabs(Right.p.y - Left.p.y);
    double Lz = fabs(Right.p.z - Left.p.z);
    return (Lx*Ly < A && Lz >= kap * H);
}

#define R 10     //ширина поиска локальных сгущений - входной параметр алгоритма

double squaredDistance(myPoint* a, myPoint* b){
    return pow(a->p.x - b->p.x, 2) + pow(a->p.y-b->p.y, 2) + pow(a->p.z-b->p.z, 2);
};

int CenterObject(std::vector<myPoint*>& points, std::vector<int>& indexes, std::vector<int>& curclast){
    double mindist = 1000000000.0;
    double curdist;
    int minindex = -1;
    for (int i = 0; i < curclast.size(); ++i){
        curdist = 0.0;
        for (int j = 0; j < curclast.size(); ++j)
            curdist += sqrt(squaredDistance(points[indexes[curclast[i]]], points[indexes[curclast[j]]]));
        if (mindist > curdist){
            mindist = curdist;
            minindex = indexes[curclast[i]];
        }
    }
    return minindex;//minpoint;
};

void genClaster(std::vector<myPoint*>& points, std::vector<int>& indexes, int center, std::vector<int>* curclast){
    for (int i = 0; i < indexes.size(); ++i){
        if (squaredDistance(points[indexes[i]], points[center]) < R*R)
            curclast->push_back(i);
    }
}

//void addToClaster((std::vector<myPoint*>& points, std::vector<int>& indexes, int center, std::vector<int>* curclast){
//    for (int i = 0; i < indexes.size(); ++i){
//        if (squaredDistance(points[indexes[i]], points[center]) < R*R)
//            curclast->push_back(i);
//    }
//}

int FOREL(std::vector<myPoint*>& allpoints){
    std::vector<int> indexes;
    for (int i = 0; i < allpoints.size(); ++i)
        indexes.push_back(i);
    std::vector<int> curclast;
    int predcenter, curcenter, j = -1;
    int numClast = 0;
    while (indexes.size() > 0){
        curclast.clear();
        std::cout << indexes.size() <<std::endl;
        predcenter = indexes[rand()%indexes.size()];
        genClaster(allpoints, indexes, predcenter, &curclast);
        curcenter = CenterObject(allpoints, indexes, curclast);
        while (squaredDistance(allpoints[curcenter], allpoints[predcenter]) > 5){
            std::cout << curcenter <<" " << predcenter << " " << curclast.size() << std::endl;
            predcenter = curcenter;
            curclast.clear();
            genClaster(allpoints, indexes, predcenter, &curclast);
            curcenter = CenterObject(allpoints, indexes, curclast);
        }
        std::vector<int>::iterator it = indexes.begin();
        for (int i = 0; i < curclast.size(); ++i){
            j = indexes[curclast[i]];
            //std::cout << j;
            indexes.erase(it + curclast[i] - i);
            allpoints[j]->numb_clast= numClast;                  //+height*10
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
    std::vector<myPoint*> allpoints;
    //allpoints.reserve(40000);
    std::map<int, std::vector<myPoint*> > heights;
    //std::vector<int, std::vector<myPoint*> > heights;
    for (int i = 0; i < 20/*nBlock*/; ++i)
    {
        GetCloudPoints(i, points.data(), &pHeader);
        //fout << "Frame #" << i << std::endl;
        //fout << " X: " << pHeader.FrameXCoord << std::setprecision(10);
        //fout << " Y: " << pHeader.FrameYCoord << std::setprecision(10);
        //fout << " Z: " << pHeader.FrameZCoord << std::setprecision(10) << std::endl;
        for (int j = 0; j < points.size(); ++j){
            if (points.data()[j].x <eps && points.data()[j].y < eps && points.data()[j].z < eps)
                continue;

            //allpoints.push_back(new myPoint(&points.data()[j], -1));

            heights[(int)floor(points.data()[j].z/H)].push_back(new myPoint(&points.data()[j], -1));
            //std::cout << allpoints[allpoints.size()-1].numb_clast <<" " << heights[(int)floor(points.data()[j].z)][heights[(int)floor(points.data()[j].z)].size()-1]->numb_clast <<std::endl;
           // for (int k = 0; k < heights[(int)floor(points.data()[j].z)].size(); ++ k)
            //    std::cout << heights[(int)floor(points.data()[j].z)][k]->numb_clast;
            /*if (heights.find(h) == heights.end())
                heights[h] = std::vector<myPoint*>(1, &allpoints[allpoints.size()-1]);*/
        }

    }
    //std::cout << "All" << allpoints.size()<<std::endl;

    std::map<int, std::vector< std::vector<myPoint*> > > clasters;

    for (std::map<int, std::vector<myPoint*> >::iterator it = heights.begin(); it != heights.end(); ++it){
        std::cout <<it->first <<" "<<FOREL(it->second) << std::endl;
        //std::cout << it->second.size() << std::endl;
        for (int k = 0; k < it->second.size(); ++k){
            dlib::vector<double> val(it->second[k]->p.x, it->second[k]->p.y, it->second[k]->p.z);
            // Pick a color based on how far we are along the spiral
            dlib::rgb_pixel color = dlib::colormap_jet(it->second[k]->numb_clast,0,5);
            // And add the point to the list of points we will display
            pointsPaint.push_back(dlib::perspective_window::overlay_dot(val, color));
        }
            //std::cout<<it->second[k]->numb_clast;
    }
   /* for (std::map<int, std::vector<myPoint*> >::iterator it = heights.begin(); it != heights.end(); ++it){
        //std::cout << FOREL(it->second);
        //std::vector<myPoint*>* tmpvec = &(it->second);
        int i =  it->first;
        for  (int j = 0; j < heights[i].size(); ++j){
            dlib::vector<double> val(heights[i][j]->p.x, heights[i][j]->p.y, heights[i][j]->p.z);
            // Pick a color based on how far we are along the spiral
            dlib::rgb_pixel color = dlib::colormap_jet(heights[i][j]->numb_clast,0,20);
            // And add the point to the list of points we will display
            pointsPaint.push_back(dlib::perspective_window::overlay_dot(val, color));
        }
    }*/
    /*
     //for visualization
            dlib::vector<double> val(points.data()[j].x, points.data()[j].y, points.data()[j].z);
            // Pick a color based on how far we are along the spiral
            dlib::rgb_pixel color = dlib::colormap_jet(1,0,20);
            // And add the point to the list of points we will display
            pointsPaint.push_back(dlib::perspective_window::overlay_dot(val, color));*/
    //int FOREL(std::vector<myPoint>& allpoints)
    //std::cout << minZ << maxZ << std::endl;

    fout.close();
    CloseData();

    //for visualization
    dlib::perspective_window win;
    win.set_title("perspective_window 3D point cloud");
    win.add_overlay(pointsPaint);
    win.wait_until_closed();


    return 0;
}