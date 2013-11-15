/* this is our Fofb class header file */
#include "Conv.h"
#include <iostream>


using namespace std;

// TOLERANCE = 5E-6f;

// const float ONE_RE = 1.f;
// const float ONE_IM = 0.f;
// const unsigned int WIDTH = 2;

float** Conv::fft_priv(float**& in, 
	const float &one_re, 
	const float &one_im, 
	unsigned int &n, 
	const unsigned int &width) {
	
	
	if (n == 1) {
		return in;
	}

    // Calculate w_n

    float theta = 2*(3.1415926)/n;    
    float w_n_re = cos(theta);
    float w_n_im = sin(theta);

    // Calculate w

    float w_re = one_re;
    float w_im = one_im;

    // Separate input into odds and evens.
	
	unsigned int n_half = n/2;

	float** in_evens = new float*[width];
	float** in_odds = new float*[width];
	
	for (unsigned int i = 0; i < width; i++) {
		in_evens[i] = new float[n_half];	
		in_odds[i] = new float[n_half];
	} 

	// Reals and Imag

	unsigned int idx = 0;
	
	for (unsigned int i = 0; i < width; i++) {
		for (unsigned int j = 0; j < n; j++) {
			if (j % 2 == 0) 
				in_evens[i][idx] = in[i][j];
			else 
				in_odds[i][idx++] = in[i][j];
		}

		idx = 0; // Reset
	}	

	float** out_evens = fft_priv(in_evens, one_re, one_im, n_half, width);
	float** out_odds = fft_priv(in_odds, one_re, one_im, n_half, width);

	// Initialize out array

	float** out = new float*[width];
	
	for (unsigned int i = 0; i < width; i++) {
		out[i] = new float[n];
	}    

	float out_re;
	float out_im;

	float re_1st;
	float im_1st;
	float im_2nd;
	float re_2nd;

	float prod_re;
	float prod_im;

	unsigned int k_plus_n_half;

	for (unsigned int k = 0; k < n_half; k++) {
         
        // Complex* prod = (yOdds->at(k))->times(w);  

        out_re = out_odds[0][k];
        out_im = out_odds[1][k];

        re_1st = out_re * w_re;
        im_1st = out_re * w_im;
        im_2nd = out_im * w_re;
        re_2nd = out_im * w_im;


		prod_re = re_1st - re_2nd;
		prod_im = im_1st + im_2nd;

		// y->at(k)  = (yEvens->at(k))->plus(prod);

		out[0][k] = out_evens[0][k] + prod_re;
		out[1][k] = out_evens[1][k] + prod_im;

		k_plus_n_half = k + n_half;  

		// y->at(k_plus_n_half) = (yEvens->at(k))->minus(prod);   

		out[0][k_plus_n_half] = out_evens[0][k] - prod_re;
		out[1][k_plus_n_half] = out_evens[1][k] - prod_im;

		// w->_times(w_n);

		re_1st = w_re * w_n_re;
		im_1st = w_re * w_n_im;
		im_2nd = w_im * w_n_re;
		re_2nd = w_im * w_n_im;

		prod_re = re_1st - re_2nd;
		prod_im = im_1st + im_2nd;

		w_re = prod_re;
		w_im = prod_im;

	}

	// cleanup:

	for (unsigned int i = 0; i < width; i++) {
		if (in_evens[i] != out_evens[i]) {
			delete[] in_evens[i];
			delete[] in_odds[i];
		}
	    
        delete[] out_evens[i];
        delete[] out_odds[i];	
	}

	if (in_evens != out_evens) {
		delete[] in_evens;
		delete[] in_odds;
	}

	delete[] out_evens;
	delete[] out_odds;

	return out;

}

float** Conv::ifft_priv(float**& in, 
	const float &one_re, 
	const float &one_im, 
	unsigned int &n, 
	const unsigned int &width) {
	
	
	if (n == 1) {
		return in;
	}

    // Calculate w_n

    float theta = 2*(3.1415926)/n;    
    float w_n_re = cos(theta);
    float w_n_im = sin(theta);

    // Calculate w_n_inv = one->dividedBy(n)

    float magn1 = 1;
	float magn2 = sqrt(w_n_re * w_n_re + w_n_im * w_n_im);
	float angl1 = 0;
	float angl2 = atan2(w_n_im, w_n_re);

	float t_magn = magn1 / magn2;
	float t_angl = angl1 - angl2; 

	float w_n_inv_re = t_magn*cos(t_angl);
	float w_n_inv_im = t_magn*sin(t_angl);




    // Calculate w

    float w_re = one_re;
    float w_im = one_im;

    // Separate input into odds and evens.
	
	unsigned int n_half = n/2;

	float** in_evens = new float*[width];
	float** in_odds = new float*[width];
	
	for (unsigned int i = 0; i < width; i++) {
		in_evens[i] = new float[n_half];	
		in_odds[i] = new float[n_half];
	} 

	// Reals and Imag

	unsigned int idx = 0;
	
	for (unsigned int i = 0; i < width; i++) {
		for (unsigned int j = 0; j < n; j++) {
			if (j % 2 == 0) 
				in_evens[i][idx] = in[i][j];
			else 
				in_odds[i][idx++] = in[i][j];
		}

		idx = 0; // Reset
	}	

	float** out_evens = ifft_priv(in_evens, one_re, one_im, n_half, width);
	float** out_odds = ifft_priv(in_odds, one_re, one_im, n_half, width);

	// Initialize out array

	float** out = new float*[width];
	
	for (unsigned int i = 0; i < width; i++) {
		out[i] = new float[n];
	}    

	float out_re;
	float out_im;

	float re_1st;
	float im_1st;
	float im_2nd;
	float re_2nd;

	float prod_re;
	float prod_im;

	unsigned int k_plus_n_half;

	for (unsigned int k = 0; k < n_half; k++) {
         
        // Complex* prod = (yOdds->at(k))->times(w);  

        out_re = out_odds[0][k];
        out_im = out_odds[1][k];

        re_1st = out_re * w_re;
        im_1st = out_re * w_im;
        im_2nd = out_im * w_re;
        re_2nd = out_im * w_im;


		prod_re = re_1st - re_2nd;
		prod_im = im_1st + im_2nd;

		// y->at(k)  = (yEvens->at(k))->plus(prod);

		out[0][k] = out_evens[0][k] + prod_re;
		out[1][k] = out_evens[1][k] + prod_im;

		k_plus_n_half = k + n_half;  

		// y->at(k_plus_n_half) = (yEvens->at(k))->minus(prod);   

		out[0][k_plus_n_half] = out_evens[0][k] - prod_re;
		out[1][k_plus_n_half] = out_evens[1][k] - prod_im;

		// w->_times(w_n_inv);

		re_1st = w_re * w_n_inv_re;
		im_1st = w_re * w_n_inv_im;
		im_2nd = w_im * w_n_inv_re;
		re_2nd = w_im * w_n_inv_im;

		prod_re = re_1st - re_2nd;
		prod_im = im_1st + im_2nd;

		w_re = prod_re;
		w_im = prod_im;

	}

	// cleanup:

	for (unsigned int i = 0; i < width; i++) {
		if (in_evens[i] != out_evens[i]) {
			delete[] in_evens[i];
			delete[] in_odds[i];
		}
	    
        delete[] out_evens[i];
        delete[] out_odds[i];	
	}

	if (in_evens != out_evens) {
		delete[] in_evens;
		delete[] in_odds;
	}

	delete[] out_evens;
	delete[] out_odds;

	return out;

}

// Takes in audio input in dynamic array float*
float** Conv::fft(float*& in, float**& in_2d, unsigned int &n, unsigned int pow_of_2_size) {
	const float ONE_RE = 1;
	const float ONE_IM = 0;
	const unsigned int WIDTH = 2;

	in_2d = new float*[WIDTH];

	for (unsigned int i = 0; i < WIDTH; i++) {
		in_2d[i] = new float[pow_of_2_size];
	}

	// Copy/Initialize

	for (unsigned int i = 0; i < pow_of_2_size; i++) {
		if (i < n) {
			in_2d[0][i] = in[i];
			in_2d[1][i] = 0;		// imaginary = 0, maybe this line isn't needed?
		}
		else {
			in_2d[0][i] = 0;
			in_2d[1][i] = 0;		// imaginary = 0, maybe this line isn't needed?
		}
	}

	return fft_priv(in_2d, ONE_RE, ONE_IM, pow_of_2_size, WIDTH);
}

//
float* Conv::ifft(float**& in, unsigned int &pow_of_2_size) {
	const float ONE_RE = 1;
	const float ONE_IM = 0;
	const unsigned int WIDTH = 2;

	float** ifft_norm = Conv::ifft_priv(in, ONE_RE, ONE_IM, pow_of_2_size, WIDTH);
	float* copy = new float[pow_of_2_size];

	for (unsigned int i = 0; i < pow_of_2_size; i++) {
		ifft_norm[0][i] /= pow_of_2_size;
		copy[i] = ifft_norm[0][i];
	}

	for (unsigned int i = 0; i < WIDTH; i++) {
		delete ifft_norm[i];
	}


	delete[] ifft_norm;

	return copy;
}

// float** Conv::doConvProcess(float**& dbl_ptr_IR1, unsigned int IR1_size, float**& dbl_ptr_IR2, unsigned int IR2_size) {
// 	unsigned int DEBUG_CONV = 1;
// 	unsigned int pow_of_2_size = 2;
// 	unsigned int bigger = (IR1_size > IR2_size ? IR1_size : IR2_size);

// 	while (pow_of_2_size < bigger) 
// 		pow_of_2_size *= 2;

// 	// Multiply again by 2 to ensure that pow_of_2_size is
// 	// at least twice the size of the bigger.

// 	pow_of_2_size *= 2;

// 	float* l_IR1 = dbl_ptr_IR1[0];
// 	float* r_IR1 = dbl_ptr_IR1[1];

// 	float* l_IR2 = dbl_ptr_IR2[0];
// 	float* r_IR2 = dbl_ptr_IR2[0];

// 	// Saving here just in case we want to iterate to inspect for
// 	// debugging purposes.

// 	float** l_IR1_2d;
// 	float** r_IR1_2d;
// 	float** l_IR2_2d;
// 	float** r_IR2_2d;

// 	// Take the fft.

// 	float** l_IR1_cplx = fft(l_IR1, l_IR1_2d, IR1_size, pow_of_2_size); 
// 	float** r_IR1_cplx = fft(r_IR1, r_IR1_2d, IR1_size, pow_of_2_size);
// 	float** l_IR2_cplx = fft(l_IR2, l_IR2_2d, IR2_size, pow_of_2_size);
// 	float** r_IR2_cplx = fft(r_IR2, r_IR2_2d, IR2_size, pow_of_2_size);

// 	// Point-wise multiplication.

// 	float** l_fft_prod = fft_mult(l_IR1_cplx, l_IR2_cplx, pow_of_2_size);
// 	float** r_fft_prod = fft_mult(r_IR1_cplx, r_IR2_cplx, pow_of_2_size);

// 	// Take the ifft:

// 	unsigned int width = 2;
// 	float** dbl_ptr_conv_sig = new float*[width];

// 	// Left & right channel

// 	dbl_ptr_conv_sig[0] = ifft(l_fft_prod, pow_of_2_size);
// 	dbl_ptr_conv_sig[1] = ifft(r_fft_prod, pow_of_2_size);

// 	if (DEBUG_CONV) {
		
// 		for (unsigned int j = 0; j < pow_of_2_size; j++)
// 			std::cout << "ifft[" << 0 << "][" << j << "\t"  << ifft[0][j] <<
// 				"ifft[" << 1 << "][" << j << "]\t"  << ifft[1][j] << std::endl;
		
// 	}


// 	// cleanup:
// 	// Don't need to delete l_IR1, r_IR1, l_IR2, r_IR2
// 	// Need to delete l_IR1_2d, etc.

// 	for (unsigned int i = 0; i < width; i++) {
// 		delete[] l_IR1_2d[i];
// 		delete[] r_IR1_2d[i];
// 		delete[] l_IR2_2d[i];
// 		delete[] r_IR2_2d[i];
// 		delete[] l_IR1_cplx[i];
// 		delete[] r_IR1_cplx[i];
// 		delete[] l_IR2_cplx[i];
// 		delete[] r_IR2_cplx[i];
// 		delete[] l_fft_prod[i];
// 		delete[] r_fft_prod[i];

// 	}
	
// 	delete[] l_IR1_2d;
// 	delete[] r_IR1_2d;
// 	delete[] l_IR2_2d;
// 	delete[] r_IR2_2d;
// 	delete[] l_IR1_cplx;
// 	delete[] r_IR1_cplx;
// 	delete[] l_IR2_cplx;
// 	delete[] r_IR2_cplx;
// 	delete[] l_fft_prod;
// 	delete[] r_fft_prod;

// 	return dbl_ptr_conv_sig;
// }

float** Conv::fft_mult(float** l_cplx, float** r_cplx,
                   unsigned int n) {

	unsigned int width = 2;

	float** fft_prod = new float*[width];

	for (unsigned int i = 0; i < width; i++) {
		fft_prod[i] = new float[n];
	}

	float l_re;
	float l_im;
	float r_re;
	float r_im;
	float re_1st;
	float im_1st;
	float im_2nd;
	float re_2nd;


	for (unsigned int i = 0; i < n; i++) {

		l_re = l_cplx[0][i];
		l_im = l_cplx[1][i];

		r_re = r_cplx[0][i];
		r_im = r_cplx[1][i];

		re_1st = l_re * r_re;
        im_1st = l_re * r_im;
        im_2nd = l_im * r_re;
        re_2nd = l_im * r_im;

		fft_prod[0][i] = re_1st - re_2nd;
		fft_prod[1][i] = im_1st + im_2nd;
	}

	// for (unsigned int i = 0; i < width; i++) {
	// 	delete[] 


	return fft_prod;
}


// string Complex::toString() {
// 	ostringstream s;
// 	if (this->_im < 0) {
// 		s << this->_re << " - " << this->_im*(-1) << "i";
// 	}
// 	else {
// 		s << this->_re << " + " << this->_im << "i"; 
// 	}
// 	return s.str();
// }
