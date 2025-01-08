#include <stdio.h>
#include <algorithm>
#include <math.h>
#include "CMU418intrin.h"
#include "logger.h"
using namespace std;

void logIntVec(__cmu418_vec_int e, __cmu418_mask maskExists) {
	int exponentArray[VECTOR_WIDTH]{0};
	_cmu418_vstore_int(exponentArray, e, maskExists);
	char exponentArrayStr[VECTOR_WIDTH*10] = {0};
	for (int i = 0; i < VECTOR_WIDTH; i++) {
		char buffer[10];
		sprintf(buffer, "%d,", exponentArray[i]);
		strcat(exponentArrayStr, buffer);
	}
	CMU418Logger.addLog(exponentArrayStr, maskExists, VECTOR_WIDTH);
}

void absSerial(float* values, float* output, int N) {
    for (int i=0; i<N; i++) {
	float x = values[i];
	if (x < 0) {
	    output[i] = -x;
	} else {
	    output[i] = x;
	}
    }
}

// implementation of absolute value using 15418 instrinsics
void absVector(float* values, float* output, int N) {
    __cmu418_vec_float x;
    __cmu418_vec_float result;
    __cmu418_vec_float zero = _cmu418_vset_float(0.f);
    __cmu418_mask maskAll, maskIsNegative, maskIsNotNegative;

    //  Note: Take a careful look at this loop indexing.  This example
    //  code is not guaranteed to work when (N % VECTOR_WIDTH) != 0.
    //  Why is that the case?
    for (int i=0; i<N; i+=VECTOR_WIDTH) {

	// All ones
	maskAll = _cmu418_init_ones();

	// All zeros
	maskIsNegative = _cmu418_init_ones(0);

	// Load vector of values from contiguous memory addresses
	_cmu418_vload_float(x, values+i, maskAll);               // x = values[i];

	// Set mask according to predicate
	_cmu418_vlt_float(maskIsNegative, x, zero, maskAll);     // if (x < 0) {

	// Execute instruction using mask ("if" clause)
	_cmu418_vsub_float(result, zero, x, maskIsNegative);      //   output[i] = -x;

	// Inverse maskIsNegative to generate "else" mask
	maskIsNotNegative = _cmu418_mask_not(maskIsNegative);     // } else {

	// Execute instruction ("else" clause)
	_cmu418_vload_float(result, values+i, maskIsNotNegative); //   output[i] = x; }

	// Write results back to memory
	_cmu418_vstore_float(output+i, result, maskAll);
    }
}

// Accepts an array of values and an array of exponents
// For each element, compute values[i]^exponents[i] and clamp value to
// 4.18.  Store result in outputs.
// Uses iterative squaring, so that total iterations is proportional
// to the log_2 of the exponent
void clampedExpSerial(float* values, int* exponents, float* output, int N) {
    for (int i=0; i<N; i++) {
	float x = values[i];
	float result = 1.f;
	int y = exponents[i];
	float xpower = x;
	while (y > 0) {
	    if (y & 0x1) {
			result *= xpower;
		}
	    xpower = xpower * xpower;
	    y >>= 1;
	}
	if (result > 4.18f) {
	    result = 4.18f;
	}
	output[i] = result;
    }
}

void clampedExpVector(float* values, int* exponents, float* output, int N) {
	__cmu418_vec_int zeros = _cmu418_vset_int(0);
	__cmu418_vec_int ones = _cmu418_vset_int(1);
	__cmu418_vec_float result;
	__cmu418_mask maskExists, maskClamped;

	int curIndex = 0;
	while (curIndex < N) {
		__cmu418_mask mask, maskMulPower;
		__cmu418_vec_float v, vpower;
		__cmu418_vec_int e;
		__cmu418_vec_int eOdd = _cmu418_vset_int(0);
		result = _cmu418_vset_float(1.f);
		maskExists = _cmu418_init_ones(0);
		maskClamped = _cmu418_init_ones(0);
		int numElements = min(N - curIndex, VECTOR_WIDTH);
		mask = _cmu418_init_ones(numElements);
		maskMulPower = _cmu418_init_ones(numElements);

		_cmu418_vload_float(v, values+curIndex, mask);
		_cmu418_vload_float(vpower, values+curIndex, mask);
		_cmu418_vload_int(e, exponents+curIndex, mask);

		_cmu418_vgt_int(maskExists, e, zeros, mask);
		while(_cmu418_cntbits(maskExists) > 0) {
			// Inside clampedExpVector function
			// logIntVec(e, maskExists);
			eOdd = _cmu418_vset_int(0);
			_cmu418_vbitand_int(eOdd, e, ones, maskExists);
			maskMulPower = _cmu418_init_ones(0);
			_cmu418_vgt_int(maskMulPower, eOdd, zeros, maskExists);
			_cmu418_vmult_float(result, result, vpower, maskMulPower);
			_cmu418_vmult_float(vpower, vpower, vpower, maskExists);
			_cmu418_vshiftright_int(e, e, ones, maskExists);
			_cmu418_vgt_int(maskExists, e, zeros, maskExists);
		}

		__cmu418_vec_float clamped = _cmu418_vset_float(4.18f);
		_cmu418_vgt_float(maskClamped, result, clamped, mask);
		_cmu418_vmove_float(result, clamped, maskClamped);
		_cmu418_vstore_float(output+curIndex, result, mask);

		curIndex += VECTOR_WIDTH;
	}
}


float arraySumSerial(float* values, int N) {
    float sum = 0;
    for (int i=0; i<N; i++) {
	sum += values[i];
    }

    return sum;
}

// Assume N % VECTOR_WIDTH == 0
// Assume VECTOR_WIDTH is a power of 2
float arraySumVector(float* values, int N) {
    // Implement your vectorized version here
    __cmu418_vec_float sum = _cmu418_vset_float(0.f);
    __cmu418_mask mask = _cmu418_init_ones();
	__cmu418_vec_float add = _cmu418_vset_float(0.f);

	for(int i = 0; i < N; i += VECTOR_WIDTH) {
		_cmu418_vload_float(add, values + i, mask);
		_cmu418_vadd_float(sum, sum, add, mask);
	}

	int width = VECTOR_WIDTH;
	while(width > 1) {
		_cmu418_hadd_float(sum, sum);
		_cmu418_interleave_float(sum, sum);
		width /= 2;
	}
	float result[VECTOR_WIDTH];
	_cmu418_vstore_float(result, sum, mask);
	return result[0];
}
