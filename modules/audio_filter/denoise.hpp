
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

/* VLC core API headers */
#include <vlc_common.h>
#include <vlc_aout.h>
#include <vlc_filter.h>
#include <vlc_plugin.h>

/* need it explicitly for the denoise filter */
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <vector>

// noiseMask + windowBuffer struct
struct filter_sys_t {
	//Filter static config
	std::vector<float> noiseDensity;
	// used for window overlapping
	float* window;
	float* secondWindow;
	// marker for the level of using window
	int windowZero = 0;
	// marker for the start of overlapping
	int overlappingStart = 0;
	// give the grad of complete
	int positions = 1;
	bool maskComplete = false;
	// noiseMAsk
	float* noise_xr;
	float* noise_xi;

};

// using by fft real and imaginary part
struct frequenzStruct {
	float* xr;
	float* xi;

};

// get the size of the potens of 2
unsigned createN(unsigned samples) {


	double param, fractpart, intpart;
	param = log(samples) / log(2);
	fractpart = modf (param , &intpart);
	if (fractpart > 0) {
		intpart++;
	}

	return (unsigned) pow(2.0, intpart);

}


float* fillBuffer(unsigned N, unsigned sampleSize, unsigned mono, float* temp_buffer, float* spl) {
	for (unsigned i = N; i > 0; i--) {
		// zero Pattern / set all not used points to 0
		if (i > sampleSize) {
			temp_buffer[i - 1] = 0.0f;
		} else {
			mono--;
			temp_buffer[i - 1] = *(spl + mono--);

		}

	}

	return temp_buffer;
}

frequenzStruct* fillStruct(int n, float* temp_buffer, frequenzStruct* fS) {
	// initial fill / imaginary part part is allways 0
	for (int i = 0; i < n; i++) {
		fS->xr[i] = temp_buffer[i];
		fS->xi[i] = 0.0f;
	}

	return fS;

}

frequenzStruct* inversion(int n, frequenzStruct* fS) {
	int i, j, k;
	float hr, hi;
	for (i = j = 1, k = n / 2; i < n; i++, j += k, k = n / 2) {
		if (j > i) {
			hr = fS->xr[i - 1];
			hi = fS->xi[i - 1];
			fS->xr[i - 1] = fS->xr[j - 1];
			fS->xi[i - 1] = fS->xi[j - 1];
			fS->xr[j - 1] = hr;
			fS->xi[j - 1] = hi;
		}
		for (; k < j; j -= k, k /= 2);
	}

	return fS;

}

frequenzStruct* overdriveCutOut(int n, frequenzStruct* fS) {
	float factor = 1;
	// important part reduce with the length of samples
	factor /= n;
	for (int i = 0; i < n; i++)
	{
		fS->xr[i] *= factor;
	}
	return fS;
}


frequenzStruct* allInOnefft(int n, frequenzStruct* fS, int typ) {

	int i, j, k, l, m;
	int le;
	float hr, hi;
	float wr, wi, whr, whi;
	for (i = n, m = 0; i > 1; m++, i /= 2);
	for (i = 0; i < m; i++ ) {
		le = 1 << (i + 1);
		whr = 1;
		whi = 0;
		// to the frequenz-area
		if (typ == 0) {
			wr = (float)cos(M_PI / le * 2);
			wi = (float) - sin(M_PI / le * 2);
		} else {
			// to the time-area
			wr = (float)cos(M_PI / le * 2);
			wi = (float)sin(M_PI / le * 2);
		}
		for (j = 0; j < le / 2; j++) {
			for (k = j; k < n; k += le) {
				l = k + le / 2;
				hr = fS->xr[l] * whr - fS->xi[l] * whi;
				hi = fS->xi[l] * whr + fS->xr[l] * whi;
				fS->xr[l] = fS->xr[k] - hr;
				fS->xi[l] = fS->xi[k] - hi;
				fS->xr[k] = fS->xr[k] + hr;
				fS->xi[k] = fS->xi[k] + hi;
			}
			hr = whr * wr - whi * wi;
			hi = whi * wr + whr * wi;
			whr = hr;
			whi = hi;

		}

	}

	return fS;

}

// set the average of noise data
bool averageAmplitude(int n, filter_t* filter) {

	int length = filter->p_sys->positions;
	float xr = 0;
	float xi = 0;
	float density = 0;

	// add all date divide with count
	for (int i = 0; i < n; i++) {

		for (int j = 0; j < length; j++) {

			density += filter->p_sys->noiseDensity[i + j * n];
			xr += filter->p_sys->noise_xr[i + j * n];
			xi += filter->p_sys->noise_xi[i + j * n];
		}

		filter->p_sys->noiseDensity[i] = density /= length;
		filter->p_sys->noise_xr[i] = xr /= length;
		filter->p_sys->noise_xi[i] = xi /= length;
		xr = 0;
		xi = 0;
		density = 0;
	}
	filter->p_sys->noiseDensity.resize(n);
	return filter->p_sys->maskComplete = true;

}


// charge the power of density for all re and im
std::vector<float> powerDensity(int n, frequenzStruct* str, int frames) {

	// vector for size of frames
	std::vector<float> temp(n * frames);
	for (int i = 0; i < n; i++) {
		temp[i] = sqrt(str->xr[i] * str->xr[i] + str->xi[i] * str->xi[i]);
	}
	return temp;
}

frequenzStruct* fillNoiseMask(unsigned n, int s, filter_t* filter, frequenzStruct* data, int frames) {
	// if finished get out
	if (filter->p_sys->maskComplete) {

		return data;
	}

	int flag = 0;
	unsigned firstMatch = 0;
	int position = filter->p_sys->positions;
	std::vector<float> temp = powerDensity(n, data, 1);

	// if noise mask allready existing but not finished fill with more data
	if (filter->p_sys->noiseDensity.size() != 0) {
		int value = 0;
		for (unsigned k = filter->p_sys->noiseDensity.size() ; k < n * position; k++) {

			filter->p_sys->noiseDensity[k] = temp[value];
			filter->p_sys->noise_xr[k] = data->xr[value];
			filter->p_sys->noise_xi[k] = data->xi[value++];
		}
		filter->p_sys->noiseDensity.resize(n * position);
		// if it's complete build average of data
		if (filter->p_sys->positions == frames) {
			bool average = averageAmplitude( n, filter);

			return data;
		}
		filter->p_sys->positions++;

		return data;
	}

	for (unsigned i = 0; i < n; i++) {
		// only for the first data initial noise struct
		if (filter->p_sys->noiseDensity.size() == 0 && temp[i] != 0.0f) {
			flag = 1;
			firstMatch = i;
			filter->p_sys->noise_xr = new float[n * frames];
			filter->p_sys->noise_xi = new float[n * frames];
			break;
		}
	}
	// set first data for the mask
	if (filter->p_sys->noiseDensity.size() == 0 && flag == 1) {
		std::vector<float> init(n * frames);
		filter->p_sys->noiseDensity = init;
		// if using data not start at the first position skip
		if (firstMatch != 0) {
			filter->p_sys->noiseDensity.resize(n * frames - 1);
			return data;
		}
		for (unsigned j = 0; j < n; j++) {
			filter->p_sys->noiseDensity[j] = temp[j];
			filter->p_sys->noise_xr[j] = data->xr[j];
			filter->p_sys->noise_xi[j] = data->xi[j];
		}

		filter->p_sys->noiseDensity.resize(n);
		// set on one frame finished
		filter->p_sys->positions++;

	}

	return data;

}

// frequence subtraction
float frequSub(float value, float divider) {

	divider = value / divider;
	
	if (value < 0.0f) {
		value += divider; 

	}

	if (value > 0.0f) {

		value -= divider; 

	}

	return value;
}



// denoise that signal
frequenzStruct* denoise(int n, filter_t* filter, frequenzStruct* song) {

	if (!filter->p_sys->maskComplete)
		return song;
	std::vector<float> songBuffer = powerDensity(n, song, 1);
	std::vector<float> noiseBuffer = filter->p_sys->noiseDensity;

	// all data they are smaller or equal that gate -> reduce
	for (int i = 0; i < n; i++) {

		if (songBuffer[i] <= noiseBuffer[i] * 10000) {


			song->xr[i] = frequSub(song->xr[i], 2.0f);
			song->xi[i] = frequSub(song->xi[i], 2.0f);

			// song->xr[i] = 0.0f;
			// song->xi[i] = 0.0f;
		}
	}

	return song;
}

// calculate value of fadein
float fadeIn(float i, float length) {

	length--;
	i = i / length;
	return i;
}

// calculate value of fadeout
float fadeOut(float i, float length) {

	length--;
	i = 1 - i / length;
	return i;
}


// prepare window one for the addition
float* editWindowOne(int fadingWindow, int startFade, float* windowOne) {

	for (int i = 0; i < fadingWindow; i++) {
		windowOne[startFade++] *= fadeOut(i, fadingWindow);

	}

	return windowOne;
}

// prepare window two for the addition
float* editWindowTwo(int fadingWindow, float* windowTwo) {

	for (int i = 0; i < fadingWindow; i++) {
		windowTwo[i] *= fadeIn(i, fadingWindow);

	}

	return windowTwo;
}

// create window three for the addition
float* fillWindowThree(int fadingWindow, int startFade, float* windowOne, float* windowTwo) {

	// Size of the window three
	int windowThreeSize = fadingWindow * 2;

	// initial value of window three
	float* windowThree = new float[windowThreeSize];

	for (int s = 0; s < fadingWindow; s++) {
		windowThree[s] = windowOne[startFade++];
		windowThree[fadingWindow + s] = windowTwo[s];
	}


	// all fading position in window three
	for (int i = 0; i < fadingWindow; i++) {
		windowThree[i] *= fadeIn(i, fadingWindow);
		windowThree[i + fadingWindow] *= fadeOut(i, fadingWindow);

	}

	return windowThree;
}

// overlapped the frames for a better sound
float* windowOverlapping(int n, filter_t* filter, frequenzStruct* song) {
	if (!filter->p_sys->maskComplete)
		return song->xr;

	// fill initial the windowbuffer
	if (filter->p_sys->overlappingStart == 0) {
		filter->p_sys->overlappingStart = 1;

		// initialize the windows
		filter->p_sys->window = new float[n];
		filter->p_sys->secondWindow = new float[n];

		// set window marker to in use
		filter->p_sys->windowZero = 1;

		for (int i = 0; i < n; i++) {

			filter->p_sys->window[i] = song->xr[i];
			song->xr[i] = 0.0f;
		}

		return song->xr;
	}

	// if there only a sconde window buffer
	// then fill the first window also
	if (filter->p_sys->windowZero == 0) {

		// set window marker to in use
		filter->p_sys->windowZero = 1;

		for (int i = 0; i < n; i++) {
			filter->p_sys->window[i] = song->xr[i];
		}

		return filter->p_sys->secondWindow;
	}

	// // only for fade_ in and out
	int fadingWindow = n / 2;

	// // Startpoint of fade
	int startFade = n - fadingWindow;

	// create the windows for overlapping
	float* three = fillWindowThree(fadingWindow, startFade, filter->p_sys->window, song->xr);
	float* one = editWindowOne(fadingWindow, startFade, filter->p_sys->window);
	float* two = editWindowTwo(fadingWindow, song->xr);

	// mix the windows
	int justAdditionTwo = fadingWindow;
	for (int x = 0; x < fadingWindow; x++) {
		one[startFade++] += three[x];
		two[x] += three[justAdditionTwo++];

	}

	// fill the second window buffer for the next frame
	for (int i = 0; i < n; i++) {
		filter->p_sys->secondWindow[i] = two[i];
	}

	// the first window buffer is no longer used, set it to 0
	filter->p_sys->windowZero = 0;

	return one;
}