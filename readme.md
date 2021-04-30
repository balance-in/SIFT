# SIFT特征点提取

SIFT特征对旋转、尺度缩放、亮度变化等保持不变性，是非常稳定的局部特征。本程序利用SIFT算法提取图像中的特征点。

## 1.运行环境

本程序使用c++语言，基于c++11编写，运行环境是win10系统，编译器为mingw64。由于需要读取各种图片，使用了一个轻量级的图形库spot，读取图片并存入相关的容器类。

## 2.实验运行结果

打开release文件，添加图片并打开cmd命令行，输入SIFT.exe 图片名，即可生成提取特征点之后的图片。

下图为实验结果：

![Lenna](E:\jetbrain\SIFT\picture\Lenna.jpg)

![SIFT_Lenna](E:\jetbrain\SIFT\picture\SIFT_Lenna.jpg)

上图为原图片，下图为提取了特征点的图片，黄色点即为提取的特征点。

## 3. 算法软件说明

SIFT特征点提取的具体流程主要是先对图片进行高斯模糊，构建出不同尺度空间的图片构成高斯尺度空间金字塔，高层图片由底层图片下采样得到，然后利用金字塔中上下两层图片相减得到高斯差分尺度空间DoG。极值点寻找通过比较检测点与领域的点以及相邻尺度中同位置的点，即3*3立方体中所有的点，即与26个点比较。从而得到SIFT特征点。

1. SIFT提取特征点时需要构建高斯金字塔，故首先需要对图片进行高斯模糊操作，卷积核的半径由sigma决定，一般来说超过3*sigma距离以外贡献就比较小。故二维卷积核的尺寸为（6 * sigma + 1）*（6 * sigma + 1），但由于二维卷积核耗费时间太长，再加上高斯函数具有线性可分：故可以用两个一维卷积核分别在水平和垂直方向进行卷积，相加后可得到同样的结果，且速度提升明显。示意图如下

   ![Gaussian](E:\jetbrain\SIFT\picture\Gaussian.png)

   函数声明如下：

   ```c++
   typedef  vector<vector<double>> pic;//二维容器存储图片
   ```

   ```c++
   vector<double> Gaussian(int r, double sigma);//返回一维卷积核
   void GaussianBlur(const pic &image, pic &base, double sig);//image为输入图片，base为卷积后图片，sig为卷积核参数
   ```

2. 构建高斯尺度空间金字塔时需要对图片进行下采样，对图片进行resize的方法有多种，考虑到运行时间，采用双线性插值的方法对图片进行resize。

   ```c++
   void resize(const pic &input_pic, pic &base, int height, int width);//input_pic为输入图片，base为resize后的图片，width和height为输出图片的大小
   ```

3. 定义一个SIFT类用于构建高斯差分尺度空间并提取特征点，keypoint结构体用于保存提取后的特征点。

   void buildGaussianPyramid()用于生成高斯尺度空间金字塔

   buildDoGPyramid()用于生成高斯差分尺度空间金字塔

   findScaleSpaceExtrema()用于寻找特征点

   piexl_mark()用于在图片上显示特征点

   类声明如下所示：

   ```c++
   struct keypoint{
       int x = 0;
       int y = 0;
       int layers = 0;
   };
   
   class SIFT{
   public:
       int n_features = 0;//选取特征点数目
       int nOctaveLayers = 3;//每层金字塔特征图数目
       double contrastThreshold = 0.04;//特征点阈值
       double sigma = 1.6;//卷积核参数
   public:
       explicit SIFT(int n_features=0, int nOctaveLayers=3, double contrastThreshold=0.04,double sigma=1.6) : sigma(sigma), n_features(n_features), nOctaveLayers(nOctaveLayers){};//构造函数
       void buildGaussianPyramid(const pic &base, vector<pic> &pyr, int nOctave) const;//base为输入图片,pyr为输出的高斯尺度空间金字塔
       void buildDoGPyramid(const vector<pic> &gpyr, vector<pic> &dogpyr) const;//gpyr为高斯尺度空间金字塔，dogpyr为高斯差分尺度空间金字塔
       void findScaleSpaceExtrema(const vector<pic> &gauss_pyr, const vector<pic> &dog_pykey_points为提取特征点r, vector<keypoint> & key_points) const;//gauss_pyr为高斯尺度空间金字塔，dog_pyr为高斯差分尺度空间金字塔,key_points为提取特征点
       void piexl_mark(spot::image &img, int radius, vector<keypoint> & key_points, spot::color &hsl) const;
   };//img为原图像,radius为标记点大小,key_points为提取特征点,hsl为标记点颜色
   ```

   