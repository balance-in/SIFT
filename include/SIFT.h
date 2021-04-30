//
// Created by balance on 2021/4/29.
//

#ifndef SIFT_SIFT_H
#define SIFT_SIFT_H
#include <iostream>
#include <spot.hpp>
#include <cmath>
#include <spot.hpp>
#define PI  3.14159265358979323846
using namespace std;;

typedef  vector<vector<double>> pic;
vector<double> Gaussian(int r, double sigma);
void GaussianBlur(const pic &image, pic &base, double sig);
void resize(const pic &input_pic, pic &base, int height, int width);
void pic_subtract(const pic &src1, const pic &src2, pic &base);

struct keypoint{
    int x = 0;
    int y = 0;
    int layers = 0;
};

class SIFT{
public:
    int n_features = 0;
    int nOctaveLayers = 3;
    double contrastThreshold = 0.04;
    double sigma = 1.6;
public:
    explicit SIFT(int n_features=0, int nOctaveLayers=3, double contrastThreshold=0.04,double sigma=1.6) :
    sigma(sigma), n_features(n_features), nOctaveLayers(nOctaveLayers){};
    void buildGaussianPyramid(const pic &base, vector<pic> &pyr, int nOctave) const;
    void buildDoGPyramid(const vector<pic> &gpyr, vector<pic> &dogpyr) const;
    void findScaleSpaceExtrema(const vector<pic> &gauss_pyr, const vector<pic> &dog_pyr, vector<keypoint> & key_points) const;
    void piexl_mark(spot::image &img, int radius, vector<keypoint> & key_points, spot::color &hsl) const;
};



#endif //SIFT_SIFT_H
