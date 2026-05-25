#include "daisy_seed.h"
#include "daisysp.h"
#include <string>
#include <cmath>
#define _USE_MATH_DEFINES
#define M_PI 3.14159265358979
#include <math.h>

#define WINDOWLEN 2048 // pitch shift window length
#define LOG2N 11 // log_2(WINDOWLEN)
#define HOPLEN 512 // pitch shift hop length
#define SAMPLERATE 96000 // sample rate in samples per second
#define MAXDELAY 5.0 // maximum delay time in seconds

using namespace daisy;
using namespace daisysp;
using daisy::Encoder;
using daisy::Switch;

DaisySeed hw;

GPIO LEDB1, LEDB2, LEDB3, LEDB4;
GPIO Pitch1, Pitch2;
GPIO EnableSwitch;
GPIO BeatSwitch;
GPIO ModeStandard, ModeTap, ModePitch;
GPIO Amp1, Amp2, Amp3, Dl1, Dl2, Dl3;
GPIO Seg1, Seg2, Seg3, Seg4, Seg5, Seg6, Seg7;
GPIO LED1, LED2, LED3, LED4;
GPIO onLED;
GPIO CharPin1, CharPin2;
AdcChannelConfig adc_config;

// LED Buttons:
// 1: pin 15
// 2: pin 13
// 3: pin 9
// 4: pin 11
// LED button LEDs:
// 1: pin 12
// 2: pin 10
// 3: pin 8
// 4: pin 5
// Pitch Buttons:
// 1: pin 6
// 2: pin 7
// Modes:
// Mode 0: pin 2
// Mode 1: pin 3
// Mode 2: pin 4
// Beat Mode: pin 29
// Enable Switch: pin 14
// Encoders:
// Amp A,B,?: pin 31, 32, 30
// Delay A,B,?: pin 34, 35, 33
// 7 Seg Segments:
// 1: pin 26
// 2: pin 25
// 3: pin 24
// 4: pin 23
// 5: pin 1
// 6: pin 36
// 7: pin 37
// Char 1: pin 27
// Char 2: pin 28
// onLED: pin 22

int pStartTime = 0;
int c1Time = 40;
int c2Time = 80;
int c3Time = 120;
int c4Time = 160;
int beatMode = 0; // 3 or 4 beats per delay cycle
int beat1 = 1; // Beat 1 activated
int beat2 = 0; // Beat 2 activated
int beat3 = 0; // Beat 3 activated
int beat4 = 0; // Beat 4 activated
int selectedBeat = 1; // During pitch shift mode - which is selected, only used bt I/O
int mode = 0; // 0 for standard, 1 for pitch shift, 2 for tap tempo
int amplitude = -60; // Amplitude shift for delay, in cB (centiBels)
int beat1Pitch = 0; // Pitch shift for beat 1
int beat2Pitch = 0; // Pitch shift for beat 2
int beat3Pitch = 0; // Pitch shift for beat 3
int beat4Pitch = 0; // Pitch shift for beat 4
int enabledF = 0; // 0 for bypass, 1 for enabled
size_t tempo; // tempo, in beats per minute. Controls delay time
size_t delaytime; // delay time, in samples
float DSY_SDRAM_BSS delaybuf[(size_t)(MAXDELAY * SAMPLERATE)] = {0.0}; // Buffer for delay
float DSY_SDRAM_BSS seqbuf[(size_t)(MAXDELAY * SAMPLERATE)] = {0.0}; // Buffer for sequencer
size_t delaybufindex; // Index of buffers for delay and sequencer
bool delayReset; // whether the buffers have just been reset within the last cycle

// ------------- Pitch Shift Code Begins --------------

bool approximate(double a, double b) {
    double eps = pow(10, -8);
    return abs(a - b) < eps;
}

int approx_ceil(double a) {
    if (approximate(a, std::floor(a)))
        return std::floor(a);
    return std::ceil(a);
}

int approx_floor(double a) {
    if (approximate(a, std::ceil(a)))
        return std::ceil(a);
    return std::floor(a);
}


float complex_magnitude(float real, float imaginary) {
    return std::sqrt(real*real + imaginary*imaginary);
}

float complex_phase(float real, float imaginary) {
    if (real == 0) {
        if (imaginary > 0) {
            return M_PI/2;
        }
        return -M_PI/2;
    }
    float output = atan(imaginary / real);
    if (real > 0) {
        return output;
    }
    if (imaginary > 0) {
        return M_PI + output;
    }
    return -M_PI + output;
}

float assign_polar_real(float m, float p) {
    return m * cos(p);
}

float assign_polar_imaginary(float m, float p) {
    return m * sin(p);
}

float rotate_180_real(float real, float imaginary) {
    return -real;
}
float rotate_180_imaginary(float real, float imaginary) {
    return -imaginary;
}

float rotate_90_real(float real, float imaginary){
    return -imaginary;
}
float rotate_90_imaginary(float real, float imaginary){
    return real;
}

float rotate_270_real(float real, float imaginary){
    return imaginary;
}
float rotate_270_imaginary(float real, float imaginary){
    return -real;
}

float complex_multiply_real(float a_real, float a_imaginary, float b_real, float b_imaginary) {
    return a_real * b_real - a_imaginary * b_imaginary;
}
float complex_multiply_imaginary(float a_real, float a_imaginary, float b_real, float b_imaginary) {
    return a_real * b_imaginary + a_imaginary * b_real;
}

size_t bit_reverse(size_t in, size_t nbits) {
    if (nbits == 1) {
        return in;
    }
    
    size_t lbit = in >> (nbits - 1);
    size_t remain = in - (lbit << (nbits - 1));
    return lbit + 2 * bit_reverse(remain, nbits - 1);
}

float DSY_SDRAM_BSS fourier_temp_r[WINDOWLEN];
float DSY_SDRAM_BSS fourier_temp_i[WINDOWLEN];

float DSY_SDRAM_BSS twiddle_dft_r[WINDOWLEN / 2];
float DSY_SDRAM_BSS twiddle_dft_i[WINDOWLEN / 2];

float DSY_SDRAM_BSS twiddle_idft_r[WINDOWLEN / 2];
float DSY_SDRAM_BSS twiddle_idft_i[WINDOWLEN / 2];

size_t DSY_SDRAM_BSS bit_reverses[WINDOWLEN];




void dft(float in_r[], float in_i[], float out_r[], float out_i[]) {

    float* left_r;
    float* left_i;

    float* right_r;
    float* right_i;

    float* helper_r;
    float* helper_i;

    float temp_r;
    float temp_i;

    if ((LOG2N % 2) == 0) {
        left_r = out_r;
        left_i = out_i;

        right_r = fourier_temp_r;
        right_i = fourier_temp_i;
    } else {
        left_r = fourier_temp_r;
        left_i = fourier_temp_i;

        right_r = out_r;
        right_i = out_i;
    }
    
    for (size_t i = 0; i < WINDOWLEN; i++) {
        left_r[i] = in_r[bit_reverses[i]];
        left_i[i] = in_i[bit_reverses[i]];
    }

    for (size_t i = 0; i < LOG2N; i++){
        size_t add_dist = 1 << i;
        size_t group_count = (WINDOWLEN >> (i + 1));

        for (size_t g = 0; g < group_count; g++) {
            for (size_t j = 0; j < add_dist; j++) {
                size_t top_index = g * (add_dist * 2) + j;
                size_t bottom_index = top_index + add_dist;

                temp_r = complex_multiply_real(left_r[bottom_index],
                    left_i[bottom_index],
                    twiddle_dft_r[j * group_count],
                    twiddle_dft_i[j * group_count]);
                temp_i = complex_multiply_imaginary(left_r[bottom_index],
                    left_i[bottom_index],
                    twiddle_dft_r[j * group_count],
                    twiddle_dft_i[j * group_count]);
                
                left_r[bottom_index] = temp_r;
                left_i[bottom_index] = temp_i;

                right_r[top_index] = left_r[top_index] + left_r[bottom_index];
                right_i[top_index] = left_i[top_index] + left_i[bottom_index];
                
                right_r[bottom_index] = left_r[top_index] - left_r[bottom_index];
                right_i[bottom_index] = left_i[top_index] - left_i[bottom_index];
            }
        }
        helper_r = left_r;
        helper_i = left_i;

        left_r = right_r;
        left_i = right_i;

        right_r = helper_r;
        right_i = helper_i;
    }

    return;
}

void idft(float in_r[], float in_i[], float out_r[], float out_i[]) {
    float* left_r;
    float* left_i;

    float* right_r;
    float* right_i;

    float* helper_r;
    float* helper_i;

    float temp_r;
    float temp_i;


    if ((LOG2N % 2) == 0) {
        left_r = out_r;
        left_i = out_i;

        right_r = fourier_temp_r;
        right_i = fourier_temp_i;
    } else {
        left_r = fourier_temp_r;
        left_i = fourier_temp_i;

        right_r = out_r;
        right_i = out_i;
    }

    for (size_t i = 0; i < WINDOWLEN; i++) {
        left_r[i] = in_r[i];
        left_i[i] = in_r[i];
    }
    
    for (size_t i = 0; i < LOG2N; i++){
        size_t add_dist = (WINDOWLEN >> (i + 1));
        size_t group_count = 1 << i;

        for (size_t g = 0; g < group_count; g++) {
            for (size_t j = 0; j < add_dist; j++) {
                size_t top_index = g * (add_dist * 2) + j;
                size_t bottom_index = top_index + add_dist;

                right_r[top_index] = left_r[top_index] + left_r[bottom_index];
                right_i[top_index] = left_i[top_index] + left_i[bottom_index];

                temp_r = left_r[top_index] - left_r[bottom_index];
                temp_i = left_i[top_index] - left_i[bottom_index];

                right_r[bottom_index] = complex_multiply_real(temp_r,
                    temp_i,
                    twiddle_idft_r[j * group_count],
                    twiddle_idft_i[j * group_count]);

                right_i[bottom_index] = complex_multiply_imaginary(temp_r,
                    temp_i,
                    twiddle_idft_r[j * group_count],
                    twiddle_idft_i[j * group_count]);

            }
        }
        helper_r = left_r;
        helper_i = left_i;

        left_r = right_r;
        left_i = right_i;

        right_r = helper_r;
        right_i = helper_i;
    }

    for (size_t i = 0; i < WINDOWLEN; i++) {
        fourier_temp_r[i] = out_r[i];
        fourier_temp_i[i] = out_i[i];
    }

    for (size_t i = 0; i < WINDOWLEN; i++) {
        out_r[i] = fourier_temp_r[bit_reverses[i]] / WINDOWLEN;
        out_i[i] = fourier_temp_i[bit_reverses[i]] / WINDOWLEN;
    }

    return;
}



size_t shift_index_in = 0; // Index of sample since reset
size_t shift_index_in_mod = 0; // shift_index_in % WINDOWLEN

size_t shift_index_out = 0; // shift_index_in - 2 * WINDOWLEN
size_t shift_index_out_mod = 0; // shift_index_out % (4 * WINDOWLEN)

float DSY_SDRAM_BSS factors[4]; // Pitch shifting factors

float DSY_SDRAM_BSS inbufs[4][WINDOWLEN] = {0.0}; // Buffers for input
float DSY_SDRAM_BSS hann[WINDOWLEN] = {0.0}; // Hann window

float DSY_SDRAM_BSS inframes_r[4][WINDOWLEN]; // Inputs of DFT (real part)
float DSY_SDRAM_BSS inframes_i[4][WINDOWLEN]; // Inputs of DFT (imaginary part)
float DSY_SDRAM_BSS outframes_r[4][WINDOWLEN]; // Outputs of DFT (real part)
float DSY_SDRAM_BSS outframes_i[4][WINDOWLEN]; // Outputs of DFT (imaginary part)

float DSY_SDRAM_BSS mags[4][2][WINDOWLEN]; // Magnitudes of frequency components of frames
float DSY_SDRAM_BSS pdifs[4][2][WINDOWLEN]; // Phase differences of frequency components of frames
float DSY_SDRAM_BSS plast[4][WINDOWLEN]; // Last phases for input
float DSY_SDRAM_BSS plast_out[4][WINDOWLEN]; // last phases for output
int circ = 0; // Circular buffer index for mags/pdifs/plast - used for second dimension [-][circ][-]

float ks[4] = {1.0}; // Proportion of first frame, for interpolating mags/pdfs/plast

bool pitch_shift_reset = false; // flag to reset when available
bool pitch_shift_reset_0 = true; // whether the pitch-shifting is just reset
bool pitch_shift_reset_1 = false; // one frame after reset

bool window_flag = false;

float DSY_SDRAM_BSS outbufs[4][4 * WINDOWLEN] = {0.0}; // Buffers for output
float DSY_SDRAM_BSS hannbufs[4][4 * WINDOWLEN] = {0.0}; // Buffers for summed Hann windows
float DSY_SDRAM_BSS outbuf_index[4]; // Current place in output buffer
int DSY_SDRAM_BSS outbuf_refresh[4]; // Place in output buffer to refresh

double semitones_to_ratio(int semitones) {
	return pow(2.0, (double) semitones / 12.0);
}

//Set scaling factors for each pitch
void r_t_update() {
	factors[0] = (float) semitones_to_ratio(beat1Pitch);
	factors[1] = (float) semitones_to_ratio(beat2Pitch);
	factors[2] = (float) semitones_to_ratio(beat3Pitch);
	factors[3] = (float) semitones_to_ratio(beat4Pitch);

	return;
}


// Setup for pitch-shifting
void r_t_setup() {

    twiddle_dft_r[0] = 1.0;
    twiddle_idft_r[0] = 1.0;
    twiddle_dft_i[WINDOWLEN/4] = -1.0;
    twiddle_idft_i[WINDOWLEN/4] = 1.0;

	// double angle;
	// float hann_val;

	// int j;

    // for (int i = 1; i < (WINDOWLEN/4); i++) {
	// 	j = i + WINDOWLEN/4;
    //     angle = -2.0 * M_PI * (double) i / (double) WINDOWLEN;
    //     twiddle_dft_r[i] = assign_polar_real(1.0, angle);
    //     twiddle_dft_i[i] = assign_polar_imaginary(1.0, angle);

    //     twiddle_idft_r[i] = twiddle_dft_r[i];
    //     twiddle_idft_i[i] = -twiddle_dft_i[i];

    //     twiddle_dft_r[j] = rotate_270_real(twiddle_dft_r[i], twiddle_dft_i[i]);
    //     twiddle_dft_i[j] = rotate_270_imaginary(twiddle_dft_r[i], twiddle_dft_i[i]);
        
    //     twiddle_idft_r[j] = twiddle_dft_r[j];
    //     twiddle_idft_i[j] = -twiddle_dft_i[j];
    // }

    // // Initialize the values of the Hann window
    // for (int i = 0; i < WINDOWLEN; i++) {
    //     hann_val = pow(sin((float) i * M_PI / ((float) WINDOWLEN - 1.0)), 2.0);
    //     hann[i] = hann_val;
	// 	bit_reverses[i] = bit_reverse(i, LOG2N);
    // }

	// r_t_update();
}

void r_t_reset() {
	// Reset everything!
	shift_index_in = 0;
	shift_index_in_mod = 0;
	shift_index_out = 0;
	shift_index_out_mod = 0;

	bool pitch_shift_reset_0 = true;
	bool pitch_shift_reset_1 = false;
	circ = 0;

	for (size_t i = 0; i < 4; i++) {
		outbuf_index[i] = 0.0; // Current place in output buffer
		outbuf_refresh[i] = 0; // Place in output buffer to refresh

		ks[i] = 1.0;
	}
	r_t_update();
}



// Upload a window for pitch-shifting
void load_window(float* w, size_t offset, size_t sel) {

    for (size_t i = 0; i < WINDOWLEN; i++) {
        // Set the value of the input frames
        inframes_r[sel][i] = w[(i + offset) % WINDOWLEN] * hann[i];
        inframes_i[sel][i] = 0.0;
    }

    dft(inframes_r[sel], inframes_i[sel], outframes_r[sel], outframes_i[sel]);

    float newphase;
    for (size_t i = 0; i < WINDOWLEN; i++) {
        mags[sel][circ][i] = complex_magnitude(outframes_r[sel][i], outframes_i[sel][i]);
        newphase = complex_phase(outframes_r[sel][i], outframes_i[sel][i]);
        if (pitch_shift_reset_0) {
            pdifs[sel][circ][i] = newphase;
        } else {
            pdifs[sel][circ][i] = std::fmod(newphase - plast[sel][i], 2 * M_PI);
        }
        
        plast[sel][i] = newphase;
    }

    circ = 1 - circ; // Swap between 0 and 1
    return;
}



void get_frame(float factor, size_t sel) {
	float k = ks[sel];
    float ko = 1 - k; // 1 - k
    int co = 1 - circ; // 1 - circ
    float mag; // magnitude of Complex
    float phs; // phase of Complex

    if (approximate(k, 1.0)) {
        // Interpolated frame lies on an input frame

        for (size_t i = 0; i < WINDOWLEN; i++) {
            mag = mags[sel][co][i];
            if (pitch_shift_reset_0) {
                phs = pdifs[sel][co][i];
            } else {
                phs = plast_out[sel][i] + pdifs[sel][co][i];
            }

            inframes_r[sel][i] = assign_polar_real(mag, phs);
            inframes_i[sel][i] = assign_polar_imaginary(mag, phs);

            plast_out[sel][i] = phs;
        }
    } else {
        // Interpolated frime lies between two input frames
        for (size_t i = 0; i < WINDOWLEN; i++) {
            mag = ko * mags[sel][circ][i] + k * mags[sel][co][i];
            if (pitch_shift_reset_0) { // This situation should not occur!
                phs = pdifs[sel][co][i];
            } else if (pitch_shift_reset_1) {
                phs = plast_out[sel][i] + pdifs[sel][co][i]; 
            } else {
                phs = plast_out[sel][i] + pdifs[sel][co][i] ;
            }

            inframes_r[sel][i] = assign_polar_real(mag, phs);
            inframes_i[sel][i] = assign_polar_imaginary(mag, phs);

            plast_out[sel][i] = phs;
        }
    }

    idft(inframes_r[sel], inframes_i[sel], outframes_r[sel], outframes_i[sel]);

    size_t lb = (size_t) approx_ceil(outbuf_index[sel] + 1 / factor); // Lower bound
    size_t ub = (size_t) approx_floor(outbuf_index[sel] + (WINDOWLEN - 1) / factor); // Upper bound

    int i2 = 1; // Output frame index
    float ir = outbuf_index[sel] + (float) i2 / factor;
    
    float k2;
    float k2o;

    int imod; // i % (4 * WINDOWLEN)


    for (size_t i = lb; i <= ub; i++){
        imod = i % (4 * WINDOWLEN);

        if (outbuf_refresh[sel] >= outbuf_index[sel]) {
            if ((imod >= outbuf_refresh[sel]) || (imod < outbuf_index[sel])) {
                outbufs[sel][imod] = 0.0;
                hannbufs[sel][imod] = 0.0;
            }
        } else {
            if ((imod >= outbuf_refresh[sel]) && (imod < outbuf_index[sel])) {
                outbufs[sel][imod] = 0.0;
                hannbufs[sel][imod] = 0.0;
            }
        }

        // Increment until ready
        while((ir < (int) i) && (!approximate(ir, i))) {
            i2++;
            ir = outbuf_index[sel] + (float) i2 / factor;
        }

        //std::cout << " " << i2 << "\n";
        if (approximate(ir, i)) {
            outbufs[sel][imod] += outframes_r[sel][i2];
            hannbufs[sel][imod] += hann[i2];
        } else { // ir > i
            k2o = factor * (ir - i); // 1 - k2
            k2 = 1 - k2o; // (i - il) / (ir - il) - Proportion of i between "il" [ir - 1/factor] and ir
            outbufs[sel][imod] += k2o * outframes_r[sel][i2 - 1] + k2 * outframes_r[sel][i2];
            hannbufs[sel][imod] += k2o * hann[i2 - 1] + k2 * hann[i2];
        }
    }
    outbuf_refresh[sel] = (ub + 1) % (4 * WINDOWLEN);

    if (pitch_shift_reset_1) {
        pitch_shift_reset_1 = false;
    } else if (pitch_shift_reset_0) {
        pitch_shift_reset_1 = true;
        pitch_shift_reset_0 = false;
    }

    return;
}




void window_prep() {
	for (size_t j = 0; j < 1; j++) {
		load_window(inbufs[j], (shift_index_in + 1) % WINDOWLEN, j);

		// Interpolate new frames
		while (ks[j] <= 1) {

			get_frame(factors[j], j);
			ks[j] += 1/factors[j];
			outbuf_index[j] += HOPLEN / factors[j];
			outbuf_index[j] = std::fmod(outbuf_index[j], 4.0 * WINDOWLEN);
		}
		ks[j]--;
	}

	if (pitch_shift_reset) {
		pitch_shift_reset = false;
	}
}


// ------------- Pitch Shift Code Ends --------------

float centiBelsToNumber(int cB) {
	return pow(10.0, (cB / 200.0));
}


void AudioCallback(AudioHandle::InputBuffer in, AudioHandle::OutputBuffer out, size_t size)
{
    float s; // Current input sample
	float o; // Current output sample
	float sdel1; // Sequencer beat 1 - present sample
	float sdel2; // Sequencer beat 2
	float sdel3; // Sequencer beat 3
	float sdel4; // Sequencer beat 4, if necessary
	float pdel1; // Pitch-shifted beat 1
	float pdel2; // Pitch-shifted beat 2
	float pdel3; // Pitch-shifted beat 3
	float pdel4; // Pitch-shifted beat 4
	float ddel; // Delay
	int ind; // Index calculation
	for (size_t i = 0; i < size; i++) // Loop over input block
	{
		// PART 1: Manage inputs before pitch-shifting
		s = in[0][i]; // Current input sample

		o = 0.0;
		ddel = 0.0;
		sdel1 = 0.0;
		sdel2 = 0.0;
		sdel3 = 0.0;
		sdel4 = 0.0;
		pdel1 = 0.0;
		pdel2 = 0.0;
		pdel3 = 0.0;
		pdel4 = 0.0;

		if (enabledF == 1) {

            if (beat1 == 1) {
                sdel1 = s;
            }

			if (beatMode == 0) {
				if (delayReset) {
					if ((beat4 == 1) && (delaybufindex >= std::round(3.0 * delaytime / 4))) {
						ind = delaybufindex - std::round(3.0 * delaytime / 4);
						if (ind < 0) ind += delaytime;
						sdel4 = seqbuf[ind];
					}
					if ((beat3 == 1) && (delaybufindex >= std::round(2.0 * delaytime / 4))) {
						ind = delaybufindex - std::round(2.0 * delaytime / 4);
						if (ind < 0) ind += delaytime;
						sdel3 = seqbuf[ind];
					}
					if ((beat2 == 1) && (delaybufindex >= std::round(1.0 * delaytime / 4))) {
						ind = delaybufindex - std::round(1.0 * delaytime / 4);
						if (ind < 0) ind += delaytime;
						sdel2 = seqbuf[ind];
					}
				} else {
					ddel = delaybuf[delaybufindex];
		
					if (beat2 == 1) {
						ind = delaybufindex - std::round(1.0 * delaytime / 4);
						if (ind < 0) ind += delaytime;
						sdel2 = seqbuf[ind];
					}
					if (beat3 == 1) {
						ind = delaybufindex - std::round(2.0 * delaytime / 4);
						if (ind < 0) ind += delaytime;
						sdel3 = seqbuf[ind];
					}
					if (beat4 == 1) {
						ind = delaybufindex - std::round(3.0 * delaytime / 4);
						if (ind < 0) ind += delaytime;
						sdel4 = seqbuf[ind];
					}
				}
			} else {
				if (delayReset) {
					if ((beat3 == 1) && (delaybufindex >= std::round(2.0 * delaytime / 3))) {
						ind = delaybufindex - std::round(2.0 * delaytime / 3);
						if (ind < 0) ind += delaytime;
						sdel3 = seqbuf[ind];
					}
					if ((beat2 == 1) && (delaybufindex >= std::round(1.0 * delaytime / 3))) {
						ind = delaybufindex - std::round(1.0 * delaytime / 3);
						if (ind < 0) ind += delaytime;
						sdel2 = seqbuf[ind];
					}
				} else {
					ddel = delaybuf[delaybufindex];
		
					if (beat2 == 1) {
						ind = delaybufindex - std::round(1.0 * delaytime / 3);
						if (ind < 0) ind += delaytime;
						sdel2 = seqbuf[ind];
					}
					if (beat3 == 1) {
						ind = delaybufindex - std::round(2.0 * delaytime / 3);
						if (ind < 0) ind += delaytime;
						sdel3 = seqbuf[ind];
					}
				}
			}

			// // PART 2: Put all inputs in pitch-shifter
			// inbufs[0][shift_index_in_mod] = sdel1;
			// inbufs[1][shift_index_in_mod] = sdel2;
			// inbufs[2][shift_index_in_mod] = sdel3;
			// inbufs[3][shift_index_in_mod] = sdel4;

			
			// // If a set of input frames is fully loaded
			// if (((((int) shift_index_in + 1 - (int) WINDOWLEN) % HOPLEN) == 0) && (shift_index_in + 1 >= WINDOWLEN)) {
			// 	window_flag = true; // Tells the main function that a window is ready
			// } 

			// // PART 3: Load output samples from output buffer

			// if (shift_index_in >= 2 * WINDOWLEN) {
			// 	if (hannbufs[0][shift_index_out_mod] != 0.0) {
			//         pdel1 = outbufs[0][shift_index_out_mod] / hannbufs[0][shift_index_out_mod];
			//     }

			// 	if (hannbufs[1][shift_index_out_mod] != 0.0) {
			//         pdel2 = outbufs[1][shift_index_out_mod] / hannbufs[1][shift_index_out_mod];
			//     }

			// 	if (hannbufs[2][shift_index_out_mod] != 0.0) {
			//         pdel3 = outbufs[2][shift_index_out_mod] / hannbufs[2][shift_index_out_mod];
			//     }

			// 	if (hannbufs[3][shift_index_out_mod] != 0.0) {
			//         pdel4 = outbufs[3][shift_index_out_mod] / hannbufs[3][shift_index_out_mod];
			//     }
			// }

			if (true) { // For debugging purposes, ignore pitch shift
				pdel1 = sdel1;
				pdel2 = sdel2;
				pdel3 = sdel3;
				pdel4 = sdel4;
			}

			o += ddel;
			o += pdel1;
			o += pdel2;
			o += pdel3;
			o += pdel4;

			seqbuf[delaybufindex] = s; // updates the sequencer buffer
			delaybuf[delaybufindex] = o * centiBelsToNumber(amplitude); // updates the delay buffer

		} else {
			o = s;
            seqbuf[delaybufindex] = 0.0;
			delaybuf[delaybufindex] = 0.0;
		}

        delaybufindex++; // Increment the index of the delay buffer
        if (delaybufindex >= delaytime) {
            delaybufindex = 0; 
            delayReset = false;
        }

        // if (shift_index_in >= 2 * WINDOWLEN) {
        // 	shift_index_out++;
        // 	shift_index_out_mod = shift_index_out % (4 * WINDOWLEN);
        // }
        // shift_index_in++;
        // shift_index_in_mod = shift_index_in % WINDOWLEN;
		out[0][i] = o;
		
	}
}


void updateDelay() {
	delaytime = (size_t)(SAMPLERATE * 60 / tempo);
	delayReset = true;
}

void Setup() {


	daisy::Pin L1Pin, L2Pin, L3Pin, L4Pin;
	daisy::Pin LEDPin, Char1, Char2;
	daisy::Pin SegL1, SegL2, SegL3, SegL4, SegL5, SegL6, SegL7;
	daisy::Pin mPitch, mSt, mTap, beatSwPin;
	daisy::Pin BT1,BT2,BT3,BT4,BTE,BTP1,BTP2;
	daisy::Pin Ap1,Ap2,Ap3,Dp1,Dp2,Dp3;
	L1Pin.pin = 8;
	L1Pin.port = PORTB;
	L2Pin.pin = 4;
	L2Pin.port = PORTB;
	L3Pin.pin = 10;
	L3Pin.port = PORTG;
	L4Pin.pin = 8;
	L4Pin.port = PORTC;
	LEDPin.pin = 0;
	LEDPin.port = PORTC;
	Char1.pin = 4;
	Char1.port = PORTC;
	Char2.pin = 1;
	Char2.port = PORTC;
	SegL1.pin = 6;
	SegL1.port = PORTA;
	SegL2.pin = 7;
	SegL2.port = PORTA;
	SegL3.pin = 1;
	SegL3.port = PORTB;
	SegL4.pin = 3;
	SegL4.port = PORTA;
	SegL5.pin = 12;
	SegL5.port = PORTB;
	SegL6.pin = 14;
	SegL6.port = PORTB;
	SegL7.pin = 15;
	SegL7.port = PORTB;
	mSt.pin = 11;
	mSt.port = PORTC;
	mPitch.pin = 10;
	mPitch.port = PORTC;
	mTap.pin = 9;
	mTap.port = PORTC;
	Ap1.pin = 1;
	Ap1.port = PORTA;
	Ap2.pin = 0;
	Ap2.port = PORTA;
	Ap3.pin = 4;
	Ap3.port = PORTA;
	Dp1.pin = 9;
	Dp1.port = PORTG;
	Dp2.pin = 2;
	Dp2.port = PORTA;
	Dp3.pin = 11;
	Dp3.port = PORTD;
	BT1.pin = 7;
	BT1.port = PORTB;
	BT2.pin = 9;
	BT2.port = PORTB;
	BT3.pin = 11;
	BT3.port = PORTG;
	BT4.pin = 5;
	BT4.port = PORTB;
	BTE.pin = 6;
	BTE.port = PORTB;
	BTP1.pin = 2;
	BTP1.port = PORTD;
	BTP2.pin = 12;
	BTP2.port = PORTC;
	beatSwPin.pin = 5;
	beatSwPin.port = PORTA;


	LED1.Init(L1Pin,GPIO::Mode::OUTPUT, GPIO::Pull::NOPULL, GPIO::Speed::LOW);
	LED2.Init(L2Pin,GPIO::Mode::OUTPUT, GPIO::Pull::NOPULL, GPIO::Speed::LOW);
	LED3.Init(L3Pin,GPIO::Mode::OUTPUT, GPIO::Pull::NOPULL, GPIO::Speed::LOW);
	LED4.Init(L4Pin,GPIO::Mode::OUTPUT, GPIO::Pull::NOPULL, GPIO::Speed::LOW);
	onLED.Init(LEDPin,GPIO::Mode::OUTPUT, GPIO::Pull::NOPULL, GPIO::Speed::LOW);
	CharPin1.Init(Char1,GPIO::Mode::OUTPUT, GPIO::Pull::NOPULL, GPIO::Speed::LOW);
	CharPin2.Init(Char2,GPIO::Mode::OUTPUT, GPIO::Pull::NOPULL, GPIO::Speed::LOW);
	Seg1.Init(SegL1,GPIO::Mode::OUTPUT, GPIO::Pull::NOPULL, GPIO::Speed::LOW);
	Seg2.Init(SegL2,GPIO::Mode::OUTPUT, GPIO::Pull::NOPULL, GPIO::Speed::LOW);
	Seg3.Init(SegL3,GPIO::Mode::OUTPUT, GPIO::Pull::NOPULL, GPIO::Speed::LOW);
	Seg4.Init(SegL4,GPIO::Mode::OUTPUT, GPIO::Pull::NOPULL, GPIO::Speed::LOW);
	Seg5.Init(SegL5,GPIO::Mode::OUTPUT, GPIO::Pull::NOPULL, GPIO::Speed::LOW);
	Seg6.Init(SegL6,GPIO::Mode::OUTPUT, GPIO::Pull::NOPULL, GPIO::Speed::LOW);
	Seg7.Init(SegL7,GPIO::Mode::OUTPUT, GPIO::Pull::NOPULL, GPIO::Speed::LOW);
	ModeStandard.Init(mSt,GPIO::Mode::INPUT, GPIO::Pull::PULLDOWN, GPIO::Speed::LOW);
	ModePitch.Init(mPitch,GPIO::Mode::INPUT, GPIO::Pull::PULLDOWN, GPIO::Speed::LOW);
	ModeTap.Init(mTap,GPIO::Mode::INPUT, GPIO::Pull::PULLDOWN, GPIO::Speed::LOW);
	Amp1.Init(Ap1,GPIO::Mode::INPUT, GPIO::Pull::PULLUP, GPIO::Speed::LOW);
	Amp2.Init(Ap2,GPIO::Mode::INPUT, GPIO::Pull::PULLUP, GPIO::Speed::LOW);
	Amp3.Init(Ap3,GPIO::Mode::INPUT, GPIO::Pull::PULLUP, GPIO::Speed::LOW);
	Dl1.Init(Dp1,GPIO::Mode::INPUT, GPIO::Pull::PULLUP, GPIO::Speed::LOW);
	Dl2.Init(Dp2,GPIO::Mode::INPUT, GPIO::Pull::PULLUP, GPIO::Speed::LOW);
	Dl3.Init(Dp3,GPIO::Mode::INPUT, GPIO::Pull::PULLUP, GPIO::Speed::LOW);
	LEDB1.Init(BT1,GPIO::Mode::INPUT, GPIO::Pull::PULLUP, GPIO::Speed::LOW);
	LEDB2.Init(BT2,GPIO::Mode::INPUT, GPIO::Pull::PULLUP, GPIO::Speed::LOW);
	LEDB3.Init(BT3,GPIO::Mode::INPUT, GPIO::Pull::PULLUP, GPIO::Speed::LOW);
	LEDB4.Init(BT4,GPIO::Mode::INPUT, GPIO::Pull::PULLUP, GPIO::Speed::LOW);
	EnableSwitch.Init(BTE,GPIO::Mode::INPUT, GPIO::Pull::PULLUP, GPIO::Speed::LOW);
	Pitch1.Init(BTP1,GPIO::Mode::INPUT, GPIO::Pull::PULLUP, GPIO::Speed::LOW);
	Pitch2.Init(BTP2,GPIO::Mode::INPUT, GPIO::Pull::PULLUP, GPIO::Speed::LOW);
	BeatSwitch.Init(beatSwPin,GPIO::Mode::INPUT, GPIO::Pull::PULLDOWN, GPIO::Speed::LOW);

    // // Code modified from
    // // https://electro-smith.github.io/libDaisy/md_doc_2md_2__a4___getting-_started-_a_d_cs.html
    // // Set up the ADC config with a connection to ADC input
    // adc_config.InitSingle(in[0][SIZE]);

    // // Initialize the ADC peripheral with that configuration
    // hw.adc.Init(&adc_config, 1);

    // // Start the ADC
    // hw.adc.Start();

	// return;
}

void SevenSegment(int input) {
	int c = 1;

	char c1 = '0';
	char c2 = '0';
	char c3 = '0';
	char c4 = '0';

	//if (input < 0) {
	//	std::string inString = std::to_string(input);
	//	c1 = '-';
	//	c2 = inString[0];
	//	if (inString.length() > 1) {
	//		c3 = inString[1];
	//	} else {c3 = 'n';} 
	//	if (inString.length() > 2){
	//		c4 = inString[2];
	//	} else {c4 = 'n';}

	//} else {
		std::string inString = std::to_string(input);
		c4 = inString[0];
		if (inString.length() > 1) {
			c3 = inString[1];
		} else {c3 = 'n';}
		if (inString.length() > 2){
			c2 = inString[2];
		} else {c2 = 'n';}
		if (inString.length() > 3) {
			c1 = inString[3];
		} else {c1 = 'n';}
	//}

	// split input into chars

	if (daisy::System::GetNow() > pStartTime + 20) {
		pStartTime = daisy::System::GetNow();
		c1Time = pStartTime + 5;
		c2Time = pStartTime + 10;
		c3Time = pStartTime + 15;
		c4Time = pStartTime + 20;
	}

	if ((daisy::System::GetNow() <= c1Time) && daisy::System::GetNow() > pStartTime) {
		CharPin1.Write(false);
		CharPin2.Write(false);
		if (c1 == '-'){
			Seg1.Write(false);
			Seg2.Write(false);
			Seg3.Write(false);
			Seg4.Write(false);
			Seg5.Write(false);
			Seg6.Write(false);
			Seg7.Write(true);
		} else if (c1 == '0') {
			Seg1.Write(true);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(true);
			Seg5.Write(true);
			Seg6.Write(true);
			Seg7.Write(false);
		} else if (c1 == '1') {
			Seg1.Write(false);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(false);
			Seg5.Write(false);
			Seg6.Write(false);
			Seg7.Write(false);
		} else if (c1 == '2') {
			Seg1.Write(true);
			Seg2.Write(true);
			Seg3.Write(false);
			Seg4.Write(true);
			Seg5.Write(true);
			Seg6.Write(false);
			Seg7.Write(true);
		} else if (c1 == '3') {
			Seg1.Write(true);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(true);
			Seg5.Write(false);
			Seg6.Write(false);
			Seg7.Write(true);
		} else if (c1 == '4') {
			Seg1.Write(false);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(false);
			Seg5.Write(false);
			Seg6.Write(true);
			Seg7.Write(true);
		} else if (c1 == '5') {
			Seg1.Write(true);
			Seg2.Write(false);
			Seg3.Write(true);
			Seg4.Write(true);
			Seg5.Write(false);
			Seg6.Write(true);
			Seg7.Write(true);
		} else if (c1 == '6') {
			Seg1.Write(true);
			Seg2.Write(false);
			Seg3.Write(true);
			Seg4.Write(true);
			Seg5.Write(true);
			Seg6.Write(true);
			Seg7.Write(true);
		} else if (c1 == '7') {
			Seg1.Write(true);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(false);
			Seg5.Write(false);
			Seg6.Write(false);
			Seg7.Write(false);
		} else if (c1 == '8') {
			Seg1.Write(true);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(true);
			Seg5.Write(true);
			Seg6.Write(true);
			Seg7.Write(true);
		} else if (c1 == '9') {
			Seg1.Write(true);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(false);
			Seg5.Write(false);
			Seg6.Write(true);
			Seg7.Write(true);
		} else {
			Seg1.Write(false);
			Seg2.Write(false);
			Seg3.Write(false);
			Seg4.Write(false);
			Seg5.Write(false);
			Seg6.Write(false);
			Seg7.Write(false);
		}

	} else if (daisy::System::GetNow() <= c2Time && daisy::System::GetNow() > c1Time) {

		CharPin2.Write(true);
		CharPin1.Write(false);

		if (c2 == '0') {
			Seg1.Write(true);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(true);
			Seg5.Write(true);
			Seg6.Write(true);
			Seg7.Write(false);
		} else if (c2 == '1') {
			Seg1.Write(false);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(false);
			Seg5.Write(false);
			Seg6.Write(false);
			Seg7.Write(false);
		} else if (c2 == '2') {
			Seg1.Write(true);
			Seg2.Write(true);
			Seg3.Write(false);
			Seg4.Write(true);
			Seg5.Write(true);
			Seg6.Write(false);
			Seg7.Write(true);
		} else if (c2 == '3') {
			Seg1.Write(true);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(true);
			Seg5.Write(false);
			Seg6.Write(false);
			Seg7.Write(true);
		} else if (c2 == '4') {
			Seg1.Write(false);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(false);
			Seg5.Write(false);
			Seg6.Write(true);
			Seg7.Write(true);
		} else if (c2 == '5') {
			Seg1.Write(true);
			Seg2.Write(false);
			Seg3.Write(true);
			Seg4.Write(true);
			Seg5.Write(false);
			Seg6.Write(true);
			Seg7.Write(true);
		} else if (c2 == '6') {
			Seg1.Write(true);
			Seg2.Write(false);
			Seg3.Write(true);
			Seg4.Write(true);
			Seg5.Write(true);
			Seg6.Write(true);
			Seg7.Write(true);
		} else if (c2 == '7') {
			Seg1.Write(true);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(false);
			Seg5.Write(false);
			Seg6.Write(false);
			Seg7.Write(false);
		} else if (c2 == '8') {
			Seg1.Write(true);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(true);
			Seg5.Write(true);
			Seg6.Write(true);
			Seg7.Write(true);
		} else if (c2 == '9') {
			Seg1.Write(true);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(false);
			Seg5.Write(false);
			Seg6.Write(true);
			Seg7.Write(true);
		} else {
			Seg1.Write(false);
			Seg2.Write(false);
			Seg3.Write(false);
			Seg4.Write(false);
			Seg5.Write(false);
			Seg6.Write(false);
			Seg7.Write(false);
		}

	} else if (daisy::System::GetNow() <= c3Time && daisy::System::GetNow() > c2Time) {

		CharPin2.Write(false);
		CharPin1.Write(true);

		if (c3 == '0') {
			Seg1.Write(true);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(true);
			Seg5.Write(true);
			Seg6.Write(true);
			Seg7.Write(false);
		} else if (c3 == '1') {
			Seg1.Write(false);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(false);
			Seg5.Write(false);
			Seg6.Write(false);
			Seg7.Write(false);
		} else if (c3 == '2') {
			Seg1.Write(true);
			Seg2.Write(true);
			Seg3.Write(false);
			Seg4.Write(true);
			Seg5.Write(true);
			Seg6.Write(false);
			Seg7.Write(true);
		} else if (c3 == '3') {
			Seg1.Write(true);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(true);
			Seg5.Write(false);
			Seg6.Write(false);
			Seg7.Write(true);
		} else if (c3 == '4') {
			Seg1.Write(false);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(false);
			Seg5.Write(false);
			Seg6.Write(true);
			Seg7.Write(true);
		} else if (c3 == '5') {
			Seg1.Write(true);
			Seg2.Write(false);
			Seg3.Write(true);
			Seg4.Write(true);
			Seg5.Write(false);
			Seg6.Write(true);
			Seg7.Write(true);
		} else if (c3 == '6') {
			Seg1.Write(true);
			Seg2.Write(false);
			Seg3.Write(true);
			Seg4.Write(true);
			Seg5.Write(true);
			Seg6.Write(true);
			Seg7.Write(true);
		} else if (c3 == '7') {
			Seg1.Write(true);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(false);
			Seg5.Write(false);
			Seg6.Write(false);
			Seg7.Write(false);
		} else if (c3 == '8') {
			Seg1.Write(true);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(true);
			Seg5.Write(true);
			Seg6.Write(true);
			Seg7.Write(true);
		} else if (c3 == '9') {
			Seg1.Write(true);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(false);
			Seg5.Write(false);
			Seg6.Write(true);
			Seg7.Write(true);
		} else {
			Seg1.Write(false);
			Seg2.Write(false);
			Seg3.Write(false);
			Seg4.Write(false);
			Seg5.Write(false);
			Seg6.Write(false);
			Seg7.Write(false);
		}

	} else if (daisy::System::GetNow() <= c4Time && daisy::System::GetNow() > c3Time) {

		CharPin1.Write(true);
		CharPin2.Write(true);

		if (c4 == '0') {
			Seg1.Write(true);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(true);
			Seg5.Write(true);
			Seg6.Write(true);
			Seg7.Write(false);
		} else if (c4 == '1') {
			Seg1.Write(false);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(false);
			Seg5.Write(false);
			Seg6.Write(false);
			Seg7.Write(false);
		} else if (c4 == '2') {
			Seg1.Write(true);
			Seg2.Write(true);
			Seg3.Write(false);
			Seg4.Write(true);
			Seg5.Write(true);
			Seg6.Write(false);
			Seg7.Write(true);
		} else if (c4 == '3') {
			Seg1.Write(true);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(true);
			Seg5.Write(false);
			Seg6.Write(false);
			Seg7.Write(true);
		} else if (c4 == '4') {
			Seg1.Write(false);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(false);
			Seg5.Write(false);
			Seg6.Write(true);
			Seg7.Write(true);
		} else if (c4 == '5') {
			Seg1.Write(true);
			Seg2.Write(false);
			Seg3.Write(true);
			Seg4.Write(true);
			Seg5.Write(false);
			Seg6.Write(true);
			Seg7.Write(true);
		} else if (c4 == '6') {
			Seg1.Write(true);
			Seg2.Write(false);
			Seg3.Write(true);
			Seg4.Write(true);
			Seg5.Write(true);
			Seg6.Write(true);
			Seg7.Write(true);
		} else if (c4 == '7') {
			Seg1.Write(true);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(false);
			Seg5.Write(false);
			Seg6.Write(false);
			Seg7.Write(false);
		} else if (c4 == '8') {
			Seg1.Write(true);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(true);
			Seg5.Write(true);
			Seg6.Write(true);
			Seg7.Write(true);
		} else if (c4 == '9') {
			Seg1.Write(true);
			Seg2.Write(true);
			Seg3.Write(true);
			Seg4.Write(false);
			Seg5.Write(false);
			Seg6.Write(true);
			Seg7.Write(true);
        } else if (c4 == '-') {
			Seg1.Write(false);
			Seg2.Write(false);
			Seg3.Write(false);
			Seg4.Write(false);
			Seg5.Write(false);
			Seg6.Write(false);
			Seg7.Write(true);
		} else {
			Seg1.Write(false);
			Seg2.Write(false);
			Seg3.Write(false);
			Seg4.Write(false);
			Seg5.Write(false);
			Seg6.Write(false);
			Seg7.Write(false);
		}

	}
	return;
}

int main(void) 
{
	delaybufindex = 0;
	tempo = 120;
	updateDelay();
	// r_t_setup();

	hw.Init();
	hw.SetAudioBlockSize(4); // number of samples handled per callback
	hw.SetAudioSampleRate(SaiHandle::Config::SampleRate::SAI_96KHZ);
	hw.StartAudio(AudioCallback);
	
	Setup();

	int b1last = 0;
	int b2last = 0;
	int b3last = 0;
	int b4last = 0;
	int bElast = 0;
	int bP1last = 0;
	int bP2last = 0;
	int lastTimeTap = 0;
	int activatedVar = 0; // Determines which var to display on 7seg
	bool Amplast = false;
	bool Dllast = false;
	int lastButtonPressTime = 0;

	while(1) {


		if (ModeStandard.Read() == true) { // Standard Mode
			mode = 0;
			//hw.SetLed(true);
		} if (ModePitch.Read() == true) { // Pitch Mode
			mode = 1;
			//hw.SetLed(false);
		} if (ModeTap.Read() == true) { // Tap Tempo
			mode = 2;
			//hw.SetLed(true);
		}


		if (BeatSwitch.Read() == true) {
			beatMode = 1;
			beat4 = 0;
			//hw.SetLed(true);
		} // 3
		else if (BeatSwitch.Read() == false) {beatMode = 0;} // 4

		


		bool incrementAmp = Amp1.Read();
		if (incrementAmp != Amplast){ // Amplitude Encoder Controls (Maybe prob doesnt work)
			if (incrementAmp != Amp2.Read()){
				if (amplitude < 0) {
					amplitude++;
				}
				activatedVar = 0;
			} else {
				if (amplitude > -999) {
					amplitude--;
				}
				activatedVar = 0;
			}
			Amplast = incrementAmp;
		} 
		int incrementDelay = Dl1.Read();
		if (incrementDelay != Dllast){ // Amplitude Encoder Controls (Maybe prob doesnt work)
			if (incrementDelay != Dl2.Read()){
				if (tempo < 500){
					tempo++;
					updateDelay();
					activatedVar = 1;
				}
			} else {
				if (tempo > 12){
					tempo--; 
					updateDelay();
					activatedVar = 1;
				}
			}
			Dllast = incrementDelay;
		} 
		

		if (mode != 2) { // if not in Tap Tempo, pedal switch controls enabled
			if (EnableSwitch.Read() == false) { 
				if (lastButtonPressTime + 500 < daisy::System::GetNow()){
					if (enabledF == 0) {
						enabledF = 1;
						pitch_shift_reset = true;
						// r_t_reset();
					}
					else if (enabledF == 1) {enabledF = 0;}
					bElast = 1;
					lastButtonPressTime = daisy::System::GetNow();
				}
			} else {bElast = 0;}
			if (enabledF == 1) {
				onLED.Write(true);
				//hw.SetLed(true);
			}
			if (enabledF == 0) {onLED.Write(false);}
		}

		if (mode == 0){ // Standard Mode Code
			hw.SetLed(true);
			if (LEDB1.Read() == false) { 
				if (b1last == 0){
					if (beat1 == 1) { 
						beat1 = 0;
					}
					else if (beat1 == 0) { 
						beat1 = 1;
					}
					b1last = 1;
				}
			} else if (LEDB1.Read() == true) {b1last = 0;}
			if (LEDB2.Read() == false) { 
				if (b2last == 0){
					if (beat2 == 1) { 
						beat2 = 0;
					}
					else if (beat2 == 0) { 
						beat2 = 1;
					}
					b2last = 1;
				}
			} else if (LEDB2.Read() == true) {b2last = 0;}
			if (LEDB3.Read() == false) { 
				if (b3last == 0){
					if (beat3 == 1) { 
						beat3 = 0;
					}
					else if (beat3 == 0) { 
						beat3 = 1;
					}
					b3last = 1;
				}
			} else if (LEDB3.Read() == true) {b3last = 0;}
			if (LEDB4.Read() == false) { 
				if (b4last == 0){
					if (beat4 == 1) { 
						beat4 = 0;
					}
					else if (beat4 == 0 && beatMode == 0) { 
						beat4 = 1;
					}
					b4last = 1;
				}
			} else if (LEDB4.Read() == true) {b4last = 0;}
	
			if (beat1 == 1) {
				LED1.Write(true);
			} else if (beat1 == 0) {
				LED1.Write(false);
			}
			if (beat2 == 1) {
				LED2.Write(true);
			} else if (beat2 == 0) {
				LED2.Write(false);
			}
			if (beat3 == 1) {
				LED3.Write(true);
			} else if (beat3 == 0) {
				LED3.Write(false);
			}
			if (beat4 == 1) {
				LED4.Write(true);
			} else if (beat4 == 0) {
				LED4.Write(false);
			}

			if (activatedVar == 0) {
				SevenSegment(amplitude);
			} else if (activatedVar == 1) {
				SevenSegment(tempo);
			}

		} else if (mode == 1) { // Pitch Shift mode Code
			hw.SetLed(false);

			if (LEDB1.Read() == false) { 
				if (b1last == 0){
					selectedBeat = 1;
					b1last = 1;
				}
			} else if (LEDB1.Read() == true) {b1last = 0;}
			if (LEDB2.Read() == false) { 
				if (b2last == 0){
					selectedBeat = 2;
					b2last = 1;
				}
			} else if (LEDB2.Read() == true) {b2last = 0;}
			if (LEDB3.Read() == false) { 
				if (b3last == 0){
					selectedBeat = 3;
					b3last = 1;
				}
			} else if (LEDB3.Read() == true) {b3last = 0;}
			if (LEDB4.Read() == false && beatMode == 0) { 
				if (b4last == 0){
					selectedBeat = 4;
					b4last = 1;
				}
			} else if (LEDB4.Read() == true) {b4last = 0;}
			

			if (Pitch1.Read() == false) { 
				if (lastButtonPressTime + 250 < daisy::System::GetNow()){
					if ((selectedBeat == 1) && (beat1Pitch < 12)) {
						beat1Pitch++;
						pitch_shift_reset = true;
						// r_t_reset();
					}
					if ((selectedBeat == 2) && (beat2Pitch < 12)) {
						beat2Pitch++;
						pitch_shift_reset = true;
						// r_t_reset();
					}
					if ((selectedBeat == 3) && (beat3Pitch < 12)) {
						beat3Pitch++;
						pitch_shift_reset = true;
						// r_t_reset();
					}
					if ((selectedBeat == 4) && (beat4Pitch < 12)) {
						beat4Pitch++;
						pitch_shift_reset = true;
						// r_t_reset();
					}
					bP1last = 1;
					lastButtonPressTime = daisy::System::GetNow();
				}
			} else if (Pitch1.Read() == true) {bP1last = 0;}
			if (Pitch2.Read() == false) { 
				if (lastButtonPressTime + 250 < daisy::System::GetNow()){
					if ((selectedBeat == 1) && (beat1Pitch > -12)) {
						beat1Pitch--;
						pitch_shift_reset = true;
						// r_t_reset();
					}
					if ((selectedBeat == 2) && (beat2Pitch > -12)) {
						beat2Pitch--;
						pitch_shift_reset = true;
						// r_t_reset();
					}
					if ((selectedBeat == 3) && (beat3Pitch > -12)) {
						beat3Pitch--;
						pitch_shift_reset = true;
						// r_t_reset();
					}
					if ((selectedBeat == 4) && (beat4Pitch > -12)) {
						beat4Pitch--;
						pitch_shift_reset = true;
						// r_t_reset();
					}
					bP2last = 1;
					lastButtonPressTime = daisy::System::GetNow();
				}
			} else if (Pitch2.Read() == true) {bP2last = 0;}

			if (selectedBeat == 1) {
				LED1.Write(true);
				LED2.Write(false);
				LED3.Write(false);
				LED4.Write(false);
				SevenSegment(beat1Pitch);
			} else if (selectedBeat == 2) {
				LED1.Write(false);
				LED2.Write(true);
				LED3.Write(false);
				LED4.Write(false);
				SevenSegment(beat2Pitch);
			} else if (selectedBeat == 3) {
				LED1.Write(false);
				LED2.Write(false);
				LED3.Write(true);
				LED4.Write(false);
				SevenSegment(beat3Pitch);
			} else if (selectedBeat == 4) {
				LED1.Write(false);
				LED2.Write(false);
				LED3.Write(false);
				LED4.Write(true);
				SevenSegment(beat4Pitch);
			}

		} else if (mode == 2) {
			if (EnableSwitch.Read() == false) {
				if (bElast = 0){
					int timeNow = daisy::System::GetNow();
					if (timeNow - lastTimeTap < 2000) {
						size_t dt = (delaytime + (timeNow - lastTimeTap) * SAMPLERATE / 1000) / 2;
						tempo = (SAMPLERATE * 60 / dt);
						updateDelay();
						activatedVar = 1;
					}
					lastTimeTap = timeNow;
					bElast = 1;
				} else if (EnableSwitch.Read() == true) {bElast = 0;}
			}
		}

		// WARNING: THIS MAY TAKE A LONG TIME!
		if ((enabledF == 1) && window_flag) {
			window_flag = false;
			// window_prep();
		}

	}
}