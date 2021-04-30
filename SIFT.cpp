//
// Created by balance on 2021/4/29.
//
#include "SIFT.h"
vector<double> Gaussian(const int r, const double sigma){
    int size = 2 * r + 1;
    vector<double> kernel_vec(size);
    double sum = 0;
    static const double SQRT2PI = sqrt(2.0 * PI);
    for (int i = 0, j = - r; j <= r; ++i, ++j) {
        kernel_vec[i] = exp(-(double)(j * j) / (2 * sigma * sigma)) / (SQRT2PI * sigma);
        sum += kernel_vec[i];
    }
    for (auto &d :kernel_vec) d = d / sum;
    return kernel_vec;
}
void GaussianBlur(const pic &image, pic &base, double sig){
    int r = ( int ) ( sig * 3 + 0.5 );
    vector<double> kernel = Gaussian(r, sig);
    int image_row = image.size();
    int image_col = image[0].size();
    base = pic(image_row, vector<double>(image_col));
    pic mid = pic(image_row, vector<double>(image_col));
    for (int i = 0; i < image_row; ++i) {
        for (int j = 0; j < image_col; ++j) {
            double rec = 0;
            for (int k = 0, l = -r; l <= r; ++k, ++l) {
                if ((j + l >= 0) && (j + l < image_col)){
                    double mul = image[i][j + l] * kernel[k];
                    rec += mul;
                }
            }
            mid[i][j] += rec;
        }
    }
    for (int j = 0; j < image_col; ++j) {
        for (int i = 0; i < image_row; ++i) {
            double rec = 0;
            for (int k = 0, l = -r; l <= r; ++k, ++l) {
                if (((i + l) >= 0) && ((i + l) < image_row)){
                    double mul = mid[i + l][j] * kernel[k];
                    rec += mul;
                }
            }
            base[i][j] = rec;
        }
    }
}
void resize(const pic &input_pic, pic &base, const int height, const int width){
    base = pic(height,vector<double>(width));
    double h_rate_scale = double(input_pic.size()) / height;
    double w_rate_scale = double(input_pic[0].size()) / width;
    for (int i = 0; i < height; ++i) {
        double x = i * h_rate_scale;
        int ox = int(x);
        for (int j = 0; j < width; ++j) {
            double y = j * w_rate_scale;
            int oy = int(y);
            if (((ox+1) >= input_pic.size()) || ((oy+1) >= input_pic[0].size())){
                base[i][j] = input_pic[ox][oy];
            }
            else{
                double u = x - ox;
                double v = y - oy;
                base[i][j] = (u)*(v)*input_pic[ox][oy] + u*(1-v)*input_pic[ox][oy+1] +
                             (1-u)*v*input_pic[ox+1][oy] + (1-u)*(1-v)*input_pic[ox+1][oy+1];
            }
        }
    }
}
void pic_subtract(const pic &src1, const pic &src2, pic &base){
    if ((src1.size() != src2.size()) || (src1[0].size() != src2[0].size())){
        cerr<<"vector size not equal";
        exit(1);
    }
    base = pic(src1.size(), vector<double>(src1[0].size()));
    for (int i = 0; i < src1.size(); ++i) {
        for (int j = 0; j < src1[0].size(); ++j) {
            base[i][j] = src1[i][j] - src2[i][j];
        }
    }
}
void SIFT::buildGaussianPyramid(const pic &base, vector<pic> &pyr, int nOctave) const {
    vector<double> sig(nOctaveLayers + 3);
    pyr.resize(nOctave * (nOctaveLayers + 3));
    sig[0] = sigma;
    double k = pow(2., 1. / nOctaveLayers);
    for (int i = 1; i < nOctaveLayers + 3; ++i) {
        double sig_prev = pow(k, (double)(i-1))*sigma;
        double sig_total = sig_prev*k;
        sig[i] = sqrt(sig_total*sig_total - sig_prev*sig_prev);
    }
    for (int o = 0; o < nOctave; ++o) {
        for (int i = 0; i < nOctaveLayers + 3; ++i) {
            pic &dst = pyr[o*(nOctaveLayers + 3) + i];
            if (o == 0 && i == 0){
                dst = base;
            }
            else if (i == 0){
                const pic &src = pyr[(o-1)*(nOctaveLayers + 3) + nOctaveLayers];
                resize(src, dst, src.size() / 2, src[0].size() / 2);
            }
            else
            {
                const pic &src = pyr[o*(nOctaveLayers + 3) + i - 1];
                GaussianBlur(src, dst, sig[i]);
            }
        }
    }
}

void SIFT::buildDoGPyramid(const vector<pic> &gpyr, vector<pic> &dogpyr) const{
    int nOctaves = (int)gpyr.size()/(nOctaveLayers + 3);
    dogpyr.resize(nOctaves*(nOctaveLayers + 2));
    for (int o = 0; o < nOctaves; ++o) {
        for (int i = 0; i < nOctaveLayers + 2; ++i) {
            const pic &src1 = gpyr[o*(nOctaveLayers + 3) + i];
            const pic &src2 = gpyr[o*(nOctaveLayers + 3) + i + 1];
            pic &dst = dogpyr[o*(nOctaveLayers + 2) + i];
            pic_subtract(src1, src2, dst);
        }
    }
}

void SIFT::findScaleSpaceExtrema(const vector<pic> &gauss_pyr, const vector<pic> &dog_pyr,
                                 vector<keypoint> &key_points) const {
    int nOctaves = (int)gauss_pyr.size()/(nOctaveLayers + 3);
    int threshold = int(0.5 * contrastThreshold / nOctaveLayers * 255);
    key_points.clear();
    keypoint kpt;
    int sift_border = 2;
    for (int o = 0; o < nOctaves; ++o) {
        for (int i = 1; i <= nOctaveLayers; ++i) {
            int idx = o*(nOctaveLayers + 2) + i;
            const pic &img = dog_pyr[idx];
            const pic &prev = dog_pyr[idx-1];
            const pic &next = dog_pyr[idx+1];
            int rows = int(img.size()), cols = int(img[0].size());
            for (int r = 1; r < rows - 1; ++r) {
                for (int c = 1; c < cols - 1; ++c) {
                    double val = img[r][c];
                    if (abs(val) > threshold &&
                    (val > 0 && val >= img[r+1][c] && val >= img[r][c+1] && val >=img[r+1][c+1] &&
                    val >= img[r-1][c] && val >= img[r][c-1] && val >= img[r-1][c-1] && val >= img[r-1][c+1] && val >= img[r+1][c-1] &&
                    val >= prev[r+1][c] && val >= prev[r][c+1] && val >= prev[r+1][c+1] && val >= prev[r][c] &&
                    val >= prev[r-1][c] && val >= prev[r][c-1] && val >= prev[r-1][c-1] && val >= prev[r-1][c+1] && val >= prev[r+1][c-1] &&
                    val >= next[r+1][c] && val >= next[r][c+1] && val >= next[r+1][c+1] && val >= next[r][c] &&
                    val >= next[r-1][c] && val >= next[r][c-1] && val >= next[r-1][c-1] && val >= next[r-1][c+1] && val >= next[r+1][c-1]) ||
                    (val < 0 && val <= img[r+1][c] && val <= img[r][c+1] && val <= img[r+1][c+1] &&
                    val <= img[r-1][c] && val <= img[r][c-1] && val >= img[r-1][c-1] && val >= img[r-1][c+1] && val <= img[r+1][c-1] &&
                    val <= prev[r+1][c] && val <= prev[r][c+1] && val <= prev[r+1][c+1] && val <= prev[r][c] &&
                    val <= prev[r-1][c] && val <= prev[r][c-1] && val <= prev[r-1][c-1] && val <= prev[r-1][c+1] && val <= prev[r+1][c-1] &&
                    val <= next[r+1][c] && val <= next[r][c+1] && val <= next[r+1][c+1] && val <= next[r][c] &&
                    val <= next[r-1][c] && val <= next[r][c-1] && val <= next[r-1][c-1] && val <= next[r-1][c+1] && val <= next[r+1][c-1])){
                        keypoint res;
                        res.x = r;
                        res.y = c;
                        res.layers = o;
                        key_points.push_back(res);
                    }
                }
            }
        }
    }


}

void SIFT::piexl_mark(spot::image &img, int radius, vector<keypoint> &key_points, spot::color &hsl) const {
    for (const auto &k : key_points){
        if (k.layers == 0){
            for (int i = -radius; i <= radius; ++i) {
                for (int j = -radius; j < radius; ++j) {
                    if (k.x + i >= 0 && k.x + i < img.w && k.y + j >= 0 && k.y + j < img.h){
                        img.at(k.x + i,k.y + j) = hsl;
                    }
                }
            }
        }
        else{
            int ix = k.x * k.layers * 2;
            int iy = k.y * k.layers * 2;
            for (int i = -radius; i <= radius; ++i) {
                for (int j = -radius; j < radius; ++j) {
                    if ((ix + i >= 0) && (ix + i < img.w) && (iy + j >= 0) && (iy + j < img.h)){
                        img.at(ix + i,iy + j) = hsl;
                    }
                }
            }
        }
    }
}










//template <typename T>
//vector<vector<T>> Gaussian_two_dim_Kernel(const int size, const double sigma){
//    vector<vector<T>> kernel_matrix(size,vector<T>(size));
//    int center = size / 2;
//    T sum = 0;
//    for (int i = 0; i < size; ++i) {
//        for (int j = 0; j < size; ++j) {
//            kernel_matrix[i][j] = (1 / (2*PI*sigma*sigma))* exp(-((i-center)*(i-center)+(j-center)*(j-center)) / (2*sigma*sigma));
//            sum += kernel_matrix[i][j];
//        }
//    }
//    for (auto &row : kernel_matrix) {
//        for (auto &col : row) col = col / sum;
//    }
//    return kernel_matrix;
//}
//vector<double> Gaussian(const int r, const double sigma){
//    int size = 2 * r + 1;
//    vector<double> kernel_vec(size);
//    double sum = 0;
//    static const double SQRT2PI = sqrt(2.0 * PI);
//    for (int i = 0, j = - r; j <= r; ++i, ++j) {
//        kernel_vec[i] = exp(-(double)(j * j) / (2 * sigma * sigma)) / (SQRT2PI * sigma);
//        sum += kernel_vec[i];
//    }
//    for (auto &d :kernel_vec) d = d / sum;
//    return kernel_vec;
//}