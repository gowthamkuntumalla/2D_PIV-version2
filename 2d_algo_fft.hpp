#ifndef _2D_ALGO_FFT_HPP_
#define _2D_ALGO_FFT_HPP_
#include <utility>
#include <string>

#define PI 3.14159265
using namespace std;
using namespace cv;

/**access in vector is [row][column]**/
/** <i><x><col> are equivalent and similarly <j><y><row> are equivalent**/
/** convention : (col,row), (x,y) order is followed**/
/** cartesian coordinates (c,r) i.e., the top left corner represents the address of subwindow**/

/*********************** Auxiliary Functions ***********************/
double avg(const int x2,const int y2, Mat image,const int x1,const int y1)//average of sub matrix
{
    double sum=0.0;
    for (int i=x1; i<x2; i++) //window starting at (x,y)
    {
        for (int j=y1; j<y2; j++) // <i><x><col> equivalent and similar
        {
            double a = (int)image.at<uchar> (j,i);//(row,column)
            sum+=a;
        }
    }
    sum/=((x2-x1)*(y2-y1));
    return sum;
}
double sd(const int x2,const int y2, Mat image,const int x1,const int y1)//standard deviation of sub matrix
{
    double var= 0.0;//variance
    double aver = avg(x2,y2,image,x1,y1);
    //myfile<<aver<<endl;
    for (int i=x1; i<x2; i++) //window starting at (x,y)
    {
        for (int j=y1; j<y2; j++)
        {
            double a = (int)image.at<uchar> (j,i);//(row,column)
            var+=((a-aver)*(a-aver));
        }
    }
    var=sqrt(var)/((x2-x1)*(y2-y1));
    return var;
}
void max_coef(vector< vector<double> > t,const int win_size,int& max_x,int& max_y )
{
    int a=-1;// all intensities are greater than this
    int m=max_x,n=max_y;
    for(int i=m; i<win_size+m; i++)
    {
        for(int j=n; j<win_size+n; j++)
        {
            if(t[j][i]>a)//access by t(row,column)
            {
                a=t[j][i];//(row,column)
                max_x=i;
                max_y=j;
            }
        }
    }
    return;
}
void twiddle(double &real, double &comp,int k, int N) // euler notation cos(a)+i*sin(a)
{
    real = cos(2*PI*k/N);// even function
    comp = sin ((-1)*2*PI*k/N);
    return;
}
/**bit_reversal**/
// preprocessed input to fft
//vector<int> bit_reversed_8= {0,4,2,6,1,5,3,7};
//vector<int> bit_reversed_16= {0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15};
//vector<int> bit_reversed_32= {0,16,8,24,4,20,12,28,2,18,10,26,6,22,14,30,1,17,9,25,5,21,13,29,3,19,11,27,7,23,15,31};
//vector<int> bit_reversed_64= {0,32,16,48,8,40,24,56,4,36,20,52,12,44,28,60,2,34,18,50,10,42,26,58,6,38,22,54,14,46,30,62,1,33,17,49,9,41,25,57,5,37,21,53,13,45,29,61,3,35,19,51,11,43,27,59,7,39,23,55,15,47,31,63};

/* // this function is used to produce above data
int bit_reversal(int i,const int bit_size)//i_bit_size
{
    int a;
    const int bit_siz=4;
    bitset<bit_siz> bit(i);// int to bit_array, bit_size= // 4 for 16, 5 for 32 etc
    for(size_t t=0; t!=bit_size/2; t++)
    {
        a=bit[t];
        bit[t]=bit[bit_size-t-1];
        bit[bit_size-t-1]=a;
    }
    return bit.to_ulong();
}
*/

/** 2D Fast Fourier Transform **/
//Cooley Tukey Algorithm - "Divide & Conquer"
void fft(Mat img,vector< vector <pair<double,double> > > &ft_img)
{

}

void fft2d(Mat img,vector< vector <pair<double,double> > > &ft_img,int c,int r,int win_size)//image to fourier transformed 2D vector
{
    // fft is done for a total of 2*win_size times ( win_size times for row & win_size times for column )
    /** row wise **/
    for(int i=0; i<win_size; i++) //'i'th row
    {
        //traverse row according to the bit reversed array
        for(int t=0; t<win_size; t++)
        {

        }
    }
    /** column wise **/
    for(int j=0; j<win_size; j++) //'j'th column
    {
        //traverse column according to the bit reversed array
        for(int t=0; t<win_size; t++)
        {

        }
    }
}

void inv_fft(vector< vector <pair<double,double> > > &ft_img,vector< vector <double> > &inv_ft_coef)//fourier transformed 2D vector to correlation coefficient(2D vector)
{
}

void vec_prod(int win_size,vector< vector <pair<double,double> > > vec_im1,vector< vector <pair<double,double> > > vec_im2,vector< vector <pair<double,double> > > &vec_im_prod) /**calculate i1*(i2')**/
{
    for(int i=0; i<win_size; i++)
    {
        for(int j=0; j<win_size; j++)
        {
            for (int k = 0; k < win_size; k++)
            {
                //(a+i*b)*(c-i*d) = ac+bd+i(bc-ad)
                vec_im_prod[i][j].first+= vec_im1[i][k].first * vec_im2[j][k].first + vec_im1[i][k].second * vec_im2[j][k].second; // since complex conjugate multiplication
                vec_im_prod[i][j].second+= vec_im1[i][k].second * vec_im2[j][k].first - vec_im1[i][k].first * vec_im2[j][k].second;
            }
        }
    }
return;
}

/** piv_2d_fft is called from main.cpp **/

void piv_2d_fft(Mat image1, Mat image2,vector< vector <pair<int,int> > > &max_points,int win_size, int i) // i= image number
{
    /******* Declarations *******/
    ofstream myfile; // like 'cout' it outputs to a file
    int initial_value = 0;
    double avg1=0,avg2=0;
    double sd1=0, sd2=0;
    /**This is valid for displacements of maximum A/3 where A is window_Size **/
    int intr_area=2*win_size;//interrogation area. say 32x32window in 64x64interrogation area
    int totrow1= image1.rows,totcol1=image1.cols; //opencv functions to get the rows and columns of image.
    int totrow2= image2.rows,totcol2=image2.cols;

    /** different vectors **/

    vector< vector <pair<int,int> > > max_coef_point;//  vector for storing of max coeff coordinates in image 2 corresponding to image 1
    max_coef_point.resize(totrow1,vector<pair<int,int> >(totcol1));//initializing the vector
    vector< vector<double> > cortable;// 2D array of correlation at various (x.y)

    vector< vector <pair<double,double> > > vec_im1; //Fourier transform For image1
    vector< vector <pair<double,double> > > vec_im2; //For image2
    vector< vector <pair<double,double> > > vec_im_prod; //For product (I1)*(I2')
    vec_im1.resize(win_size, vector <pair<double,double> > (win_size,make_pair(0,0)));//initialization
    vec_im2.resize(win_size, vector <pair<double,double> > (win_size,make_pair(0,0)));//initialization
    vec_im_prod.resize(win_size, vector <pair<double,double> > (win_size,make_pair(0,0)));//initialization

    myfile.open ("data_"+to_string(i)+".txt");// output to a text file

    /** Loop over entire image**/
    for(int c=0; c<totcol1-win_size; c+=win_size/2) // according to win_size.
    {
        for(int r=0; r<totrow1-win_size; r+=win_size/2) // 50% overlap of windows
        {
            myfile<<c<<","<<r<<", ";//initial point (x,y)
            int m=0,n=0;//the max coefficent point
            cortable = vector<vector<double> >(totrow1, vector<double>(totcol1,initial_value));
            /**do FFT for image1 & image2**/
            fft2d(image1,vec_im1,c,r,win_size);
            fft2d(image1,vec_im2,c,r,win_size);

            /**calculate i1*(i2')**/
            vec_prod(win_size,vec_im1,vec_im2,vec_im_prod); // vec_im_prod is computed.

            /**do inverse FFT**/
            inv_fft(vec_im_prod,cortable);

            /** find maximum correlation coefficient point**/
            max_coef(cortable,win_size,m,n);//myfile << "Writing this to a file.\n";

            /**Save the data in a "data_i.txt" file**/
            max_coef_point[r][c].first=n;//row index
            max_coef_point[r][c].second=m;//column index
            myfile<<m<<","<<n<<endl;//final point (x,y)
        }
    }
    max_points=max_coef_point;
    myfile.close();
    return;
}
#endif // _2D_ALGO_FFT_HPP_


