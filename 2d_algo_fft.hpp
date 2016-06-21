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
double magn_vect(vector< vector<pair<double,double> > > t,int j,int i)
{
    double ans;
    ans = (t[j][i].first*t[j][i].first)+(t[j][i].second*t[j][i].second);
    return sqrt(ans);
}

void max_coef(vector< vector<pair<double,double> > > t,const int win_size,int& max_x,int& max_y )
{
    double a=-1;// all intensities are greater than this
    for(int i=0; i<win_size; i++)
    {
        for(int j=0; j<win_size; j++)
        {
            double mag=magn_vect(t,j,i);
            if(mag>a)//access by t(row,column)
            {
                a=mag;//(row,column)
                max_x=i;
                max_y=j;
            }
        }
    }
    return;
}

/** twiddle values are used in fft_1d function**/
void twiddle(double &real, double &comp,int k, int N) // euler notation cos(a)+i*sin(a)
{
    real = cos(2*PI*k/N);// even function
    comp = sin ((-1)*2*PI*k/N);// odd function
    return;
}

/**bit_reversal**/
// preprocessed input to fft
vector<int> bit_reversed_8= {0,4,2,6,1,5,3,7};
vector<int> bit_reversed_16= {0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15};
vector<int> bit_reversed_32= {0,16,8,24,4,20,12,28,2,18,10,26,6,22,14,30,1,17,9,25,5,21,13,29,3,19,11,27,7,23,15,31};
vector<int> bit_reversed_64= {0,32,16,48,8,40,24,56,4,36,20,52,12,44,28,60,2,34,18,50,10,42,26,58,6,38,22,54,14,46,30,62,1,33,17,49,9,41,25,57,5,37,21,53,13,45,29,61,3,35,19,51,11,43,27,59,7,39,23,55,15,47,31,63};

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
void fft_1d(vector< vector <pair<double,double> > > &ft_img,int win_size,int stride, const char flag,int row,int col)// stride length =1 initially
{
    double real=0;// for twiddle function
    double comp=0;// for twiddle function
    double temp_real=0,temp_comp=0;// temporary variable for storing in inner for loop.
    if(flag=='r')//row vector fft
    {
        if(win_size==1)
        {
            //ft_img[row][col]=ft_img[row][col];
            return;//trivial case
        }
        else
        {
            fft_1d(ft_img,win_size/2,'r',2*stride,row,col);
            fft_1d(ft_img,win_size/2,'r',2*stride,row,col+stride);
            for(int k=0; k<win_size/2; k++)
            {
                twiddle(real,comp,k,win_size);
                temp_real=ft_img[row][k].first;// temporary variable
                temp_comp=ft_img[row][k].second;
                ft_img[row][k].first = temp_real+ft_img[row][k+win_size/2].first*real-ft_img[row][k+win_size/2].second*comp;
                ft_img[row][k].second = temp_comp+ft_img[row][k+win_size/2].first*comp+ft_img[row][k+win_size/2].second*real;
                ft_img[row][k+win_size/2].first = temp_real-ft_img[row][k+win_size/2].first*real-ft_img[row][k+win_size/2].second*comp;
                ft_img[row][k+win_size/2].second = temp_comp-ft_img[row][k+win_size/2].first*comp+ft_img[row][k+win_size/2].second*real;
            }
        }
        return;
    }
    if(flag=='c')//column vector fft
    {
        if(win_size==1)
        {
            //ft_img[row][col]=ft_img[row][col];
            return;//trivial case
        }
        else
        {
            fft_1d(ft_img,win_size/2,'c',2*stride,row,col);
            fft_1d(ft_img,win_size/2,'c',2*stride,row+stride,col);
            for(int k=0; k<win_size/2-1; k++)
            {
                twiddle(real,comp,k,win_size);
                temp_real=ft_img[k][col].first;// temporary variable
                temp_comp=ft_img[k][col].second;
                ft_img[k][col].first = temp_real+ft_img[k+win_size/2][col].first*real-ft_img[k+win_size/2][col].second*comp;
                ft_img[k][col].second = temp_comp+ft_img[k+win_size/2][col].first*comp+ft_img[k+win_size/2][col].second*real;
                ft_img[k+win_size/2][col].first = temp_real-ft_img[k+win_size/2][col].first*real-ft_img[k+win_size/2][col].second*comp;
                ft_img[k+win_size/2][col].second = temp_comp-ft_img[k+win_size/2][col].first*comp+ft_img[k+win_size/2][col].second*real;
            }
        }
        return;
    }
}

// fft_1d is called from fft_2d a total of 2*win_size times
/**2D forward FFT**/
void fft_2d(Mat img,vector< vector <pair<double,double> > > &ft_img,int c1,int r1,int win_size)//image to fourier transformed 2D vector
{
    // fft is done for a total of 2*win_size times ( win_size times for row & win_size times for column )
    double avg1=avg(c1+win_size,r1+win_size,img,c1,r1);
    for(int i=0; i<win_size; i++)//copy img values into ft_img
    {
        for(int j=0; j<win_size; j++)
        {
            ft_img[i][j].first = (int)img.at<uchar> (r1+i,c1+j)-avg1;
            ft_img[i][j].second = 0;
        }
    }
    /** row wise fft **/
    for(int i=0; i<win_size; i++) //'r1+i'th row
    {
        fft_1d(ft_img,win_size,'r',1,i,0);//traverse all rows .0 implies data is used just for initialization
    }
    /** column wise fft **/
    for(int j=0; j<win_size; j++) //'c1+j'th column
    {
        fft_1d(ft_img,win_size,'c',1,0,j);//traverse all columns. 0 implies data is used for initialization
    }
    return;
}

void inv_fft_1d(vector< vector <pair<double,double> > > &ft_img,int win_size,int stride, const char flag,int row,int col)// stride length =1 initially
{
    double real=0;// for twiddle function
    double comp=0;// for twiddle function
    double temp_real=0,temp_comp=0;// temporary variable for storing in inner for loop.
    if(flag=='r')//row vector fft
    {
        if(win_size==1)
        {
            ft_img[row][col].first=ft_img[row][col].first/win_size;
            ft_img[row][col].second=ft_img[row][col].second/win_size;
            return;//trivial case
        }
        else
        {
            inv_fft_1d(ft_img,win_size/2,'r',2*stride,row,col);
            inv_fft_1d(ft_img,win_size/2,'r',2*stride,row,col+stride);
            for(int k=0; k<win_size/2; k++)
            {
                twiddle(real,comp,k,win_size);
                temp_real=ft_img[row][k].first;// temporary variable
                temp_comp=ft_img[row][k].second;
                ft_img[row][k].first = (temp_real+ft_img[row][k+win_size/2].first*real+ft_img[row][k+win_size/2].second*comp)/win_size;
                ft_img[row][k].second = (temp_comp-ft_img[row][k+win_size/2].first*comp+ft_img[row][k+win_size/2].second*real)/win_size;
                ft_img[row][k+win_size/2].first = (temp_real-ft_img[row][k+win_size/2].first*real+ft_img[row][k+win_size/2].second*comp)/win_size;
                ft_img[row][k+win_size/2].second = (temp_comp+ft_img[row][k+win_size/2].first*comp+ft_img[row][k+win_size/2].second*real)/win_size;
            }
        }
        return;
    }
    if(flag=='c')//column vector fft
    {
        if(win_size==1)
        {
            ft_img[row][col].first=ft_img[row][col].first/win_size;
            ft_img[row][col].second=ft_img[row][col].second/win_size;
            return;//trivial case
        }
        else
        {
            inv_fft_1d(ft_img,win_size/2,'c',2*stride,row,col);
            inv_fft_1d(ft_img,win_size/2,'c',2*stride,row+stride,col);
            for(int k=0; k<win_size/2-1; k++)
            {
                twiddle(real,comp,k,win_size);
                temp_real=ft_img[k][col].first;// temporary variable
                temp_comp=ft_img[k][col].second;
                ft_img[k][col].first = (temp_real+ft_img[k+win_size/2][col].first*real+ft_img[k+win_size/2][col].second*comp)/win_size;
                ft_img[k][col].second = (temp_comp-ft_img[k+win_size/2][col].first*comp+ft_img[k+win_size/2][col].second*real)/win_size;
                ft_img[k+win_size/2][col].first = (temp_real-ft_img[k+win_size/2][col].first*real+ft_img[k+win_size/2][col].second*comp)/win_size;
                ft_img[k+win_size/2][col].second = (temp_comp+ft_img[k+win_size/2][col].first*comp+ft_img[k+win_size/2][col].second*real)/win_size;
            }
        }
        return;
    }
}

/**2D inverse FFT**/
void inv_fft_2d(vector< vector <pair<double,double> > > &ft_img,int win_size)//fourier transformed 2D vector to correlation coefficient(2D vector)
{
    // fft is done for a total of 2*win_size times ( win_size times for row & win_size times for column )

    /** row wise fft **/
    for(int i=0; i<win_size; i++) //'r1+i'th row
    {
        inv_fft_1d(ft_img,win_size,'r',1,i,0);//traverse all rows .0 implies data is used just for initialization
    }
    /** column wise fft **/
    for(int j=0; j<win_size; j++) //'c1+j'th column
    {
        inv_fft_1d(ft_img,win_size,'c',1,0,j);//traverse all columns. 0 implies data is used for initialization
    }
    return;
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
    double sd1=0, sd2=0;
    /**This is valid for displacements of maximum A/3 where A is window_Size **/
    int intr_size=2*win_size;//interrogation area. say 32x32window in 64x64interrogation area
    int totrow1= image1.rows,totcol1=image1.cols; //opencv functions to get the rows and columns of image.
    int totrow2= image2.rows,totcol2=image2.cols;
    int x=0,y=0;// point coordinates

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
    for(int c=0; c<(totcol1-win_size); c+=win_size/2) // according to win_size.
    {
        for(int r=0; r<(totrow1-win_size); r+=win_size/2) // 50% overlap of windows
        {
            myfile<<c<<","<<r<<", ";//initial point (x,y)
            int m=0,n=0;//displacement of the max coefficent point
            cortable = vector<vector<double> >(totrow1, vector<double>(totcol1,initial_value));

            if((c-win_size)<0)//left border
            {
                if((r-win_size)<0)//top border
                {
                    for(x=0,m=x; x<intr_size; x++)// interrogation area - column iteration
                    {
                        for(y=0,n=y; y<intr_size; y++)// interrogation area - row iteration
                        {
                            /**do FFT for image1 & image2**/
                            fft_2d(image1,vec_im1,c,r,win_size);
                            fft_2d(image2,vec_im2,c,r,win_size);

                            /**calculate i1*(i2')**/
                            vec_prod(win_size,vec_im1,vec_im2,vec_im_prod); // vec_im_prod is computed.

                            /**do inverse FFT**/
                            inv_fft_2d(vec_im_prod,win_size);

                            /** find maximum correlation coefficient point**/
                            max_coef(vec_im_prod,win_size,m,n);//myfile << "Writing this to a file.\n";

                            /**Save the data in a "data_i.txt" file**/
                            max_coef_point[r][c].first=n+r;//row index
                            max_coef_point[r][c].second=m+c;//column index
                            myfile<<m<<","<<n<<endl;//final point (x,y)
                        }
                    }
                }
                if(((r-win_size)>=0)&&((r+win_size)<=totrow2))//inner region
                {
                    for(x=0,m=x; x<intr_size; x++)
                    {
                        for(y=r-win_size,n=y; y<r+win_size; y++)
                        {
                            /**do FFT for image1 & image2**/
                            fft_2d(image1,vec_im1,c,r,win_size);
                            fft_2d(image2,vec_im2,c,r,win_size);

                            /**calculate i1*(i2')**/
                            vec_prod(win_size,vec_im1,vec_im2,vec_im_prod); // vec_im_prod is computed.

                            /**do inverse FFT**/
                            inv_fft_2d(vec_im_prod,win_size);

                            /** find maximum correlation coefficient point**/
                            max_coef(vec_im_prod,win_size,m,n);//myfile << "Writing this to a file.\n";

                            /**Save the data in a "data_i.txt" file**/
                            max_coef_point[r][c].first=n+r;//row index
                            max_coef_point[r][c].second=m+c;//column index
                            myfile<<m<<","<<n<<endl;//final point (x,y)
                        }
                    }
                }
                if((r+win_size)>totrow2)//bottom border
                {
                    for(x=0,m=x; x<intr_size; x++)
                    {
                        for(y=totrow2-intr_size,n=y; y<totrow2; y++)
                        {
                            /**do FFT for image1 & image2**/
                            fft_2d(image1,vec_im1,c,r,win_size);
                            fft_2d(image2,vec_im2,c,r,win_size);

                            /**calculate i1*(i2')**/
                            vec_prod(win_size,vec_im1,vec_im2,vec_im_prod); // vec_im_prod is computed.

                            /**do inverse FFT**/
                            inv_fft_2d(vec_im_prod,win_size);

                            /** find maximum correlation coefficient point**/
                            max_coef(vec_im_prod,win_size,m,n);//myfile << "Writing this to a file.\n";

                            /**Save the data in a "data_i.txt" file**/
                            max_coef_point[r][c].first=n+r;//row index
                            max_coef_point[r][c].second=m+c;//column index
                            myfile<<m<<","<<n<<endl;//final point (x,y)
                        }
                    }
                }
                /*
                 if((r-win_size)<=0&&(r+win_size)>=totrow2)//bad image!!
                 {
                     cerr<<"too small image";
                 }
                 */
            }
            if(((c-win_size)>=0)&&((c+win_size)<=totcol2))//inner region
            {
                if((r-win_size)<0)//top border
                {
                    for(x=c-win_size,m=x; x<c+win_size; x++)// interrogation area - column iteration
                    {
                        for(y=0,n=y; y<intr_size; y++)// interrogation area - row iteration
                        {
                            /**do FFT for image1 & image2**/
                            fft_2d(image1,vec_im1,c,r,win_size);
                            fft_2d(image2,vec_im2,c,r,win_size);

                            /**calculate i1*(i2')**/
                            vec_prod(win_size,vec_im1,vec_im2,vec_im_prod); // vec_im_prod is computed.

                            /**do inverse FFT**/
                            inv_fft_2d(vec_im_prod,win_size);

                            /** find maximum correlation coefficient point**/
                            max_coef(vec_im_prod,win_size,m,n);//myfile << "Writing this to a file.\n";

                            /**Save the data in a "data_i.txt" file**/
                            max_coef_point[r][c].first=n+r;//row index
                            max_coef_point[r][c].second=m+c;//column index
                            myfile<<m<<","<<n<<endl;//final point (x,y)
                        }
                    }
                }
                if(((r-win_size)>=0)&&((r+win_size)<=totrow2))//inner region
                {
                    for(x=c-win_size,m=x; x<c+win_size; x++)
                    {
                        for(y=r-win_size,n=y; y<r+win_size; y++)
                        {
                            /**do FFT for image1 & image2**/
                            fft_2d(image1,vec_im1,c,r,win_size);
                            fft_2d(image2,vec_im2,c,r,win_size);

                            /**calculate i1*(i2')**/
                            vec_prod(win_size,vec_im1,vec_im2,vec_im_prod); // vec_im_prod is computed.

                            /**do inverse FFT**/
                            inv_fft_2d(vec_im_prod,win_size);

                            /** find maximum correlation coefficient point**/
                            max_coef(vec_im_prod,win_size,m,n);//myfile << "Writing this to a file.\n";

                            /**Save the data in a "data_i.txt" file**/
                            max_coef_point[r][c].first=n+r;//row index
                            max_coef_point[r][c].second=m+c;//column index
                            myfile<<m<<","<<n<<endl;//final point (x,y)
                        }
                    }
                }
                if((r+win_size)>totrow2)//bottom border
                {
                    for(x=c-win_size,m=x; x<c+win_size; x++)
                    {
                        for(y=totrow2-intr_size,n=y; y<totrow2; y++)
                        {
                            /**do FFT for image1 & image2**/
                            fft_2d(image1,vec_im1,c,r,win_size);
                            fft_2d(image2,vec_im2,c,r,win_size);

                            /**calculate i1*(i2')**/
                            vec_prod(win_size,vec_im1,vec_im2,vec_im_prod); // vec_im_prod is computed.

                            /**do inverse FFT**/
                            inv_fft_2d(vec_im_prod,win_size);

                            /** find maximum correlation coefficient point**/
                            max_coef(vec_im_prod,win_size,m,n);//myfile << "Writing this to a file.\n";

                            /**Save the data in a "data_i.txt" file**/
                            max_coef_point[r][c].first=n+r;//row index
                            max_coef_point[r][c].second=m+c;//column index
                            myfile<<m<<","<<n<<endl;//final point (x,y)
                        }
                    }
                }
                /*
                 if((r-win_size)<=0&&(r+win_size)>=totrow2)//bad image!!
                 {
                     cerr<<"too small image";
                 }
                 */

            }
            if((c+win_size)>totcol2)//right border
            {
                if((r-win_size)<0)//top border
                {
                    for(x=totcol2-intr_size,m=x; x<totcol2; x++)// interrogation area - column iteration
                    {
                        for(y=0,n=y; y<intr_size; y++)// interrogation area - row iteration
                        {/**do FFT for image1 & image2**/
                            fft_2d(image1,vec_im1,c,r,win_size);
                            fft_2d(image2,vec_im2,c,r,win_size);

                            /**calculate i1*(i2')**/
                            vec_prod(win_size,vec_im1,vec_im2,vec_im_prod); // vec_im_prod is computed.

                            /**do inverse FFT**/
                            inv_fft_2d(vec_im_prod,win_size);

                            /** find maximum correlation coefficient point**/
                            max_coef(vec_im_prod,win_size,m,n);//myfile << "Writing this to a file.\n";

                            /**Save the data in a "data_i.txt" file**/
                            max_coef_point[r][c].first=n+r;//row index
                            max_coef_point[r][c].second=m+c;//column index
                            myfile<<m<<","<<n<<endl;//final point (x,y)
                        }
                    }
                }
                if(((r-win_size)>=0)&&((r+win_size)<=totrow2))//inner region
                {
                    for(x=totcol2-intr_size,m=x; x<totcol2; x++)
                    {
                        for(y=r-win_size,n=y; y<r+win_size; y++)
                        {
                            /**do FFT for image1 & image2**/
                            fft_2d(image1,vec_im1,c,r,win_size);
                            fft_2d(image2,vec_im2,c,r,win_size);

                            /**calculate i1*(i2')**/
                            vec_prod(win_size,vec_im1,vec_im2,vec_im_prod); // vec_im_prod is computed.

                            /**do inverse FFT**/
                            inv_fft_2d(vec_im_prod,win_size);

                            /** find maximum correlation coefficient point**/
                            max_coef(vec_im_prod,win_size,m,n);//myfile << "Writing this to a file.\n";

                            /**Save the data in a "data_i.txt" file**/
                            max_coef_point[r][c].first=n+r;//row index
                            max_coef_point[r][c].second=m+c;//column index
                            myfile<<m<<","<<n<<endl;//final point (x,y)
                        }
                    }
                }
                if((r+win_size)>totrow2)//bottom border
                {
                    for(x=totcol2-intr_size,m=x; x<totcol2; x++)
                    {
                        for(y=totrow2-intr_size,n=y; y<totrow2; y++)
                        {
                            /**do FFT for image1 & image2**/
                            fft_2d(image1,vec_im1,c,r,win_size);
                            fft_2d(image2,vec_im2,c,r,win_size);

                            /**calculate i1*(i2')**/
                            vec_prod(win_size,vec_im1,vec_im2,vec_im_prod); // vec_im_prod is computed.

                            /**do inverse FFT**/
                            inv_fft_2d(vec_im_prod,win_size);

                            /** find maximum correlation coefficient point**/
                            max_coef(vec_im_prod,win_size,m,n);//myfile << "Writing this to a file.\n";

                            /**Save the data in a "data_i.txt" file**/
                            max_coef_point[r][c].first=n+r;//row index
                            max_coef_point[r][c].second=m+c;//column index
                            myfile<<m<<","<<n<<endl;//final point (x,y)
                        }
                    }
                }
                /*
                 if((r-win_size)<=0&&(r+win_size)>=totrow2)//bad image!!
                 {
                     cerr<<"too small image";
                 }
                 */
            }
            /*
            if((c-win_size)<=0&&(c+win_size)>=totcol2)//bad image!!
            {
                cerr<<"too small image";
            }
            */
        }
    }
    max_points=max_coef_point;
    myfile.close();
    return;
}
#endif // _2D_ALGO_FFT_HPP_
