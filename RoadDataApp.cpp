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
#include <tuple>


struct myPoint{
    RoadPoint p;
    int numb_clast;
    myPoint(){
    p.x = 0;
    p.y = 0;
    p.z = 0;};
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

typedef std::pair<myPoint, myPoint> Box;
bool IntersectRect(Box a, Box b) {
    return !(a.first.p.x > b.second.p.x
        || a.second.p.x < b.first.p.x
        || a.first.p.y > b.second.p.y
        || a.second.p.y < b.first.p.y);
}

#define A 0.49
#define H 1
#define kap 0.5
#define R 6     //radius of clasterization

//BOUNDING BOX
std::pair<myPoint, myPoint> createBoundingBox(std::vector<myPoint*>& points){
    auto xExtremes = std::minmax_element(points.begin(), points.end(),
                                         [](const myPoint* lhs, const myPoint* rhs) {
                                             return lhs->p.x < rhs->p.x;
                                         });

    auto yExtremes = std::minmax_element(points.begin(), points.end(),
                                         [](const myPoint* lhs, const myPoint* rhs) {
                                             return lhs->p.y < rhs->p.y;
                                         });
    auto zExtremes = std::minmax_element(points.begin(), points.end(),
                                         [](const myPoint* lhs, const myPoint* rhs) {
                                             return lhs->p.z < rhs->p.z;
                                         });
    
    myPoint Left((*xExtremes.first)->p.x, (*yExtremes.first)->p.y, (*zExtremes.first)->p.z, points[0]->numb_clast);
    myPoint Right((*xExtremes.second)->p.x, (*yExtremes.second)->p.y, (*zExtremes.second)->p.z, points[0]->numb_clast);
    return std::make_pair(Left, Right);
}

bool checkBox(myPoint Left, myPoint Right){
    double Lx = fabs(Right.p.x - Left.p.x);
    double Ly = fabs(Right.p.y - Left.p.y);
    double Lz = fabs(Right.p.z - Left.p.z);
    //std::cout << "LLLL " << Lx << " "<<  Ly << " "<< Lz << std::endl;
    return (Lx*Ly < A && Lz >= kap * H);
}


double squaredDistance(myPoint* a, myPoint* b){
    return pow(a->p.x - b->p.x, 2) + pow(a->p.y-b->p.y, 2) + pow(a->p.z-b->p.z, 2);
};

myPoint CenterObject(std::vector<myPoint*>& points, std::vector<int>& indexes, std::vector<int>& curclast){
    myPoint centre;
    for (int i = 0; i < curclast.size(); ++i){
        centre.p.x += points[indexes[curclast[i]]]->p.x;
        centre.p.y += points[indexes[curclast[i]]]->p.y;
        centre.p.z += points[indexes[curclast[i]]]->p.z;
    }
    centre.p.x /= curclast.size();
    centre.p.y /= curclast.size();
    centre.p.z /= curclast.size();
    return centre;
};

void genClaster(std::vector<myPoint*>& points, std::vector<int>& indexes, myPoint center, std::vector<int>* curclast){
    for (int i = 0; i < indexes.size(); ++i){
        if (squaredDistance(points[indexes[i]], &center) < R*R)
            curclast->push_back(i);
    }
}

struct myClast{
    std::pair<myPoint, myPoint> box;
    std::vector<myPoint*> elems;

    myClast(){};
};
int FOREL(std::vector<myPoint*>& allpoints, std::vector< myClast >& clasters){
    std::vector<int> indexes;
    for (int i = 0; i < allpoints.size(); ++i)
        indexes.push_back(i);
    std::vector<int> curclast;
    myPoint predcenter, curcenter;
    int j = -1;
    int numClast = 0;
    while (indexes.size() > 0){
        curclast.clear();
        //std::cout << indexes.size() <<std::endl;
        predcenter = *allpoints[indexes[rand()%indexes.size()]];
        genClaster(allpoints, indexes, predcenter, &curclast);
        curcenter = CenterObject(allpoints, indexes, curclast);
        while (squaredDistance(&curcenter, &predcenter) > 0.00001){
            //std::cout << curcenter.p.x <<" " << predcenter.p.x << " " << curclast.size() << std::endl;
            predcenter = curcenter;
            curclast.clear();
            genClaster(allpoints, indexes, predcenter, &curclast);
            curcenter = CenterObject(allpoints, indexes, curclast);
        }
        std::vector<int>::iterator it = indexes.begin();
        if (clasters.size() < numClast + 1)
            clasters.resize(numClast + 1);
        for (int i = 0; i < curclast.size(); ++i){
            j = indexes[curclast[i]-i];
            indexes.erase(it + curclast[i] - i);
            allpoints[j]->numb_clast= numClast;
            clasters[numClast].elems.push_back(allpoints[j]);
        }
        std::cout << "SIZEofclaster " << clasters[numClast].elems.size()<<std::endl;
        std::pair<myPoint, myPoint> box = createBoundingBox(clasters[numClast].elems);
        clasters[numClast].box = box;
        if (!checkBox(std::get<0>(box), std::get<1>(box))){
            clasters[numClast].elems.clear();                 //comment this 2 lines to see
            --numClast;                                       //whole clasterization
        }
        ++numClast;
    }
    if (clasters[clasters.size() - 1].elems.empty())
        clasters.pop_back();
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

    std::map<int, std::vector<myPoint*> > heights;
    std::vector< myClast > clustersALL;
    //int outer = 0;
for (int outer = 0; outer < nBlock / 20; ++outer){                  // outer
    for (int i = outer*20; i < 20*(outer+1)/*nBlock*/; ++i)
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
        }

    }

    std::map<int, std::vector< myClast > > clasters;

    for (std::map<int, std::vector<myPoint*> >::iterator it = heights.begin(); it != heights.end(); ++it){
        int z = it->first;
        std::cout <<it->first <<" "<<FOREL(it->second, clasters[z]);
        std::cout << " " << clasters[z].size() << std::endl;
        if (!clasters[z].size())
            clasters.erase(z);

        /* PaintALL
        for (int k = 0; k < it->second.size(); ++k){
            dlib::vector<double> val(it->second[k]->p.x, it->second[k]->p.y, it->second[k]->p.z);
            dlib::rgb_pixel color = dlib::colormap_jet(it->second[k]->numb_clast,0,6);
            pointsPaint.push_back(dlib::perspective_window::overlay_dot(val, color));
        } */
    }
    bool f = false;
    for (std::map<int, std::vector< myClast > >::const_reverse_iterator itlevel1 = clasters.crbegin(); itlevel1 != clasters.crend(); ++itlevel1){
        int level1 = itlevel1->first;
        for (int clastl1 = 0; clastl1 < itlevel1->second.size(); ++clastl1){
                std::map<int, std::vector< myClast > >:: const_reverse_iterator itlevel2 = itlevel1;
                ++itlevel2;
            //for (std::map<int, std::vector< myClast > >:: const_reverse_iterator itlevel2 = itlevel1; itlevel2 != clasters.crend(); ++itlevel2){
                int level2 = itlevel2->first;
                f = false;
                for (int clastl2 = 0; clastl2 < itlevel2->second.size(); ++clastl2){
                    if (IntersectRect(itlevel1->second[clastl1].box, itlevel2->second[clastl2].box)){
                        std::cout << "yep" << " ";
                        f = true;
                        clasters[level2][clastl2].elems.insert(clasters[level2][clastl2].elems.end(), 
                                                                clasters[level1][clastl1].elems.begin(), clasters[level1][clastl1].elems.end());
                        clasters[level1][clastl1].elems.clear();
                        clasters[level2][clastl2].box = createBoundingBox(clasters[level2][clastl2].elems);
                        //it->second[l].elems.insert(it->second[l].elems.end(), it2->second[l2].elems.begin(), it2->second[l2].elems.end());
                    }
                }
                if (!f){ //copy down
                    std::cout << "none";
                    clasters[level2].push_back(clasters[level1][clastl1]);
                    clasters[level1][clastl1].elems.clear();
                }
        }
        clasters[level1].clear();
        clasters.erase(level1);
    }
    auto itlevel1 = clasters.begin();
    auto itlevel2  = clasters.begin();
    ++itlevel1;
    int level1 = itlevel1->first;
    for (int clastl1 = 0; clastl1 < itlevel1->second.size(); ++clastl1){
            //for (std::map<int, std::vector< myClast > >:: const_reverse_iterator itlevel2 = itlevel1; itlevel2 != clasters.crend(); ++itlevel2){
            int level2 = itlevel2->first;
            f = false;
            for (int clastl2 = 0; clastl2 < itlevel2->second.size(); ++clastl2){
                if (IntersectRect(itlevel1->second[clastl1].box, itlevel2->second[clastl2].box)){
                    //std::cout << "yep" << " ";
                    f = true;
                    clasters[level2][clastl2].elems.insert(clasters[level2][clastl2].elems.end(), 
                                                            clasters[level1][clastl1].elems.begin(), clasters[level1][clastl1].elems.end());
                    clasters[level1][clastl1].elems.clear();
                    clasters[level2][clastl2].box = createBoundingBox(clasters[level2][clastl2].elems);
                    //it->second[l].elems.insert(it->second[l].elems.end(), it2->second[l2].elems.begin(), it2->second[l2].elems.end());
                }
            }
            if (!f){ //copy down
                //std::cout << "none";
                clasters[level2].push_back(clasters[level1][clastl1]);
                clasters[level1][clastl1].elems.clear();
            }
    }
    clasters[level1].clear();
    clasters.erase(level1);


    std::cout <<"end " << clasters.size() << " " <<clasters.begin()->second.size() << std::endl;
    /*for (std::map<int, std::vector< myClast > >::const_reverse_iterator itlevel1 = clasters.crbegin(); itlevel1 != clasters.crend(); ++itlevel1){
        for (int l = 0; l < itlevel1->second.size(); ++l){
            for (int k = 0; k < itlevel1->second[l].elems.size(); ++k){
                dlib::vector<double> val(itlevel1->second[l].elems[k]->p.x, itlevel1->second[l].elems[k]->p.y, itlevel1->second[l].elems[k]->p.z);
                dlib::rgb_pixel color = dlib::colormap_jet(l+5,0, itlevel1->second.size()+5);
                pointsPaint.push_back(dlib::perspective_window::overlay_dot(val, color));
            }
        }
    }*/
    clustersALL.insert(clustersALL.end(), clasters.begin()->second.begin(),clasters.begin()->second.end());
    clasters.clear();
}  //outer
    
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