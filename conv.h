#include <cmath>
// #include <iostream>

using namespace std;

// enum{DEBUG_CONV = 1};
class Conv {
    
public:
    static float** fft(float* &in,float** &in_2d, 
        unsigned int &n, unsigned int pow_of_2_size);

    static float* ifft(float** &in, unsigned int &pow_of_2_size);

    static float** fft_priv(float** &in, const float &one_re, 
        const float &one_im, unsigned int &n, 
        const unsigned int &width);

    static float** ifft_priv(float** &in, const float &one_re, 
        const float &one_im, unsigned int &n, 
        const unsigned int &width);

    static float** fft_mult(float** l_cplx, float** r_cplx,
                   unsigned int n);


    // static float** doConvProcess(float**& dbl_ptr_IR1, unsigned int IR1_size,
    //     float**& dbl_ptr_IR2, unsigned int IR2_size);
   
    
};

