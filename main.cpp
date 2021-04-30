#include <iostream>
#include <spot.hpp>
#include <cmath>
#include "SIFT.h"
#include <chrono>
#define PI  3.14159265358979323846

using namespace std;
typedef vector<vector<double>> pic;




template <typename T>
void image_to_2DimVec(const spot::image &img, vector<vector<T>> &vec){
    vector<unsigned char> img_uchar = img.y();
    vector<double> image_gray_one_dim(img_uchar.begin(), img_uchar.end());
    for (decltype(image_gray_one_dim.size()) i = 0; i < image_gray_one_dim.size(); ++i) {
        vec[i / img.w][i % img.w] = image_gray_one_dim[i];
    }
}
template <typename T>
void vector_image_save(const string &filename,const vector<vector<T>> &img_vec){
    vector<unsigned char> DoG;
    for (const auto &i : img_vec)
        for (const auto &j : i) {DoG.push_back(j);}
    cout << spot::write_bmp(filename.c_str(),img_vec[0].size(),img_vec.size(),1,&DoG[0]) << endl;
}




int main(int argc, char* argv[]) {
    int size = 3;
    string filename = argv[1];
    spot::image img(filename);
    spot::color mark = spot::color(0.3,0.7,0.9,1);
    pic image_result(img.h, vector<double>(img.w));
    image_to_2DimVec(img, image_result);
    vector<pic> LoG;
    vector<pic> DoG;
    vector<keypoint> kpt;
    SIFT sift;
    sift.buildGaussianPyramid(image_result,LoG,5);
    sift.buildDoGPyramid(LoG, DoG);
    sift.findScaleSpaceExtrema(LoG, DoG, kpt);
    sift.piexl_mark(img,2,kpt,mark);
    img.save("SIFT_" + filename);
    return 0;
}
