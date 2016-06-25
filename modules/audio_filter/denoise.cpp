/*****************************************************************************
 * denoise.cpp : denoise mode
 *****************************************************************************
 * Copyright © 2016 Hagen Handtke
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston MA 02110-1301, USA.
 *****************************************************************************/



#include "denoise.hpp"

// declaration the methode the important methode, using by VLC
static int Open (vlc_object_t *);
// static void Close (vlc_object_t *);
static block_t* Process (filter_t *, block_t *);

// Module descriptor
vlc_module_begin ()
set_shortname (N_("Denoise"))
set_description (N_("Denoise that Song"))
// kind of audio / #define  CAT_AUDIO   2
set_category (CAT_AUDIO)
set_subcategory (SUBCAT_AUDIO_AFILTER)
// kind of module
set_capability ("audio filter", 0)
set_callbacks (Open, NULL)
vlc_module_end ()


// Start that module
static int Open (vlc_object_t* obj)
{
	filter_t* filter = (filter_t *)obj;

	if (filter->fmt_in.audio.i_channels != 2)
	{
		msg_Err (filter, "voice removal requires stereo");
		return VLC_EGENERIC;
	}

	// Allocate structure && definate frameBuffer
	filter_sys_t* p_sys = filter->p_sys = new filter_sys_t;
	if (!p_sys)
		return VLC_ENOMEM;

	filter->fmt_in.audio.i_format = VLC_CODEC_FL32;
	filter->fmt_out.audio = filter->fmt_in.audio;
	filter->pf_audio_filter = Process;
	return VLC_SUCCESS;
}

// static void Close (vlc_object_t* obj) {

// 	filter_t* filter = (filter_t *)obj;
// 	// free the your allocated memory
// 	delete[] filter->p_sys;

// }


// Main methode with the work
static block_t* Process (filter_t* filter, block_t* block) {
	// i_nb_samples get the size of sample
	const unsigned sampleSize = block->i_nb_samples;
	// create a size with a potenz of 2
	int n = (int) createN(block->i_nb_samples);
	// the p_buffer have all audio data
	float* spl = (float *)block->p_buffer;
	float* temp_buffer_ = new float[n];


	// real and imagin part/ for the fft
	frequenzStruct* fftStruct = new frequenzStruct;

	fftStruct->xr = new float[n];
	fftStruct->xi = new float[n];

	// to get only one chanel
	unsigned mono = sampleSize * 2;

	// fill the buffer with all data
	float* temp_buffer = fillBuffer(n, sampleSize, mono, temp_buffer_, spl);


	// build denoice mask
	// start with fft analyse
	// fill noise struct
	frequenzStruct* fft_Struct = fillStruct(n, temp_buffer, fftStruct);
	// sort the sampling points for the fft
	frequenzStruct* fftStr = inversion(n, fft_Struct);
	// // transform samples to the frequenz Area
	frequenzStruct* forNoiseMask = allInOnefft(n, fftStr, 0);
	// create the noise - mask
	frequenzStruct* frequenzArea = fillNoiseMask(n, sampleSize, filter, forNoiseMask, 39);

	// if noise mask created -> start with denoise
	frequenzStruct* frequenzAreaCo = denoise(n, filter, frequenzArea);


	// prepare fft
	frequenzStruct* fftStrTwo = inversion(n, frequenzAreaCo);
	// frequenzStruct* fftStrTwo = inversion(n, frequenzArea);

	// back to time - area
	frequenzStruct* timeArea = allInOnefft(n, fftStrTwo, 1);
	// reduce the output fo fft
	frequenzStruct* timeAreaI = overdriveCutOut(n, timeArea);

	// window overlapping function
	float* temp_bufferOverlapping = windowOverlapping(sampleSize, filter, timeAreaI);

	// set data back to the block
	for (unsigned i = 0; i < sampleSize; i++) {
		// float s = timeAreaI->xr[i];
	
		float s = *(temp_bufferOverlapping++);
		// channel stereo
		* (spl++) = s;
		* (spl++) = s;
	}

	(void) filter;
	// finish get the block back
	return block;
}