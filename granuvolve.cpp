#include "Granuvolve.h"
#include <cmath> 

using namespace std;

Granuvolve::Granuvolve (audioMasterCallback audioMaster)
: AudioEffectX (audioMaster, 0, PARAMS) { 
	
	setNumInputs(2);
	setNumOutputs(2);

	IR1RecWasJustOn = false;
    IR2RecWasJustOn = false;
    GR1RecWasJustOn = false;
    noteOnGR1WasJustOn = false;
    playConvWasJustOn = false;
    playBufIR1WasJustOn = false;
    playBufIR2WasJustOn = false;
    playBufGR1WasJustOn = false;

    pow_of_2_size = 2;      // only use this many samples
    _vecFrames = 0;                      // Only used to display # of frames.
    dur_convProc = 0;

	verysmall = 1.0e-4;

    isADSR = true;    
    isStdDev_GR1 = true;

    stdDev_GR1 = 1;
    randomProb_GR1 = 0.0;

    adsr_start1 = 1;
    adsr_start2 = 1;
    adsr_start3 = 1;

    adsr_end1 = 1;
    adsr_end2 = 1;
    adsr_end3 = 1;

    isNoteOn_GR1 = false;

    adsr_offset1 = 0;
    adsr_offset2 = 0;
    adsr_offset3 = 0;

    dur_sampAttack_GR1 = getSampleRate();
	dur_sampDecay_GR1 = getSampleRate();
	dur_sampRelease_GR1 = getSampleRate();

	dur_sampAttackMax_GR1 = 5 * getSampleRate();               // Maximum possible samples
	dur_sampDecayMax_GR1 = 5 * getSampleRate();
	dur_sampReleaseMax_GR1 = 5 * getSampleRate();

    endSampIdx_GR1 = 0;

    maxLength_IR1 = 10*getSampleRate();    
    maxLength_IR2 = maxLength_IR1;  
    maxLength_GR1 = maxLength_IR1;
    tailLenMax_IFFT = 1*getSampleRate();
    tailLen_IFFT = 0;      

    max_size_ifft = 0;

    dur_repeat_GR1 = 0.f;
    dur_ffts = 0.f;
    dur_freq_mult = 0.f;
    dur_iffts = 0.f;

    tolerance = 5e-4;
                    
	isIR1Rec = false;            // record Rec is set to 0 as default.
    isIR2Rec = false;
	isGR1Rec = false;

    isIR1DoneRec = false;
    isIR2DoneRec = false;
    isGR1DoneRec = false;

    isIR1RecNotFully = false;
    isIR2RecNotFully = false;
	isGR1RecNotFully = false;
 
	isIR1BufFull = false;                 // no impulse respse recorded as default.
    isIR2BufFull = false; 
	isGR1BufFull = false;
	
	isIFFTBufFull = false;
    isMuteInput = false;

    isPlayConv = false;
    isPlayBufIR1 = false;
    isPlayBufIR2 = false;
	isPlayBufGR1 = false;

    isRevBufIR1 = false;
    isRevBufIR2 = false;
	isRevBufGR1 = false;

    isRevIR1Full = false;
    isRevIR2Full = false;
    isRevGR1Full = false;
    
	isClearBufIR1 = true;
    isClearBufIR2 = true;
	isClearBufGR1 = true;

    hasJustClearedIR1 = false;
    hasJustClearedIR2 = false;
    isDoConvProc = false;
    isConvProcDone = false;
    isGranulate = false;
	isShiftingStartStopIdx_GR1 = false;
	isGR1Rec = false;

	isFindSil_IFFT = false;

    fracLength_IR1 = 1;
    fracLength_IR2 = 1;
    fracLength_GR1 = 1;

    writeIndex_IR1 = 0;
    writeIndex_IR2 = 0;
    writeIndex_GR1 = 0;

    readIndex_IR1 = 0;
    readIndex_IR2 = 0;
    readIndex_GR1 = 0;
    readIndex_ifft = 0;

    startIdx_IR1 = 0;
    startIdx_IR2 = 0;
    startIdx_GR1 = 0;

    stopIdx_IR1 = 0;
    stopIdx_IR2 = 0;
    stopIdx_GR1 = 0;

    startSampIdx_IR1 = 0;
    startSampIdx_IR2 = 0;
    startSampIdx_GR1 = 0;

    stopSampIdx_IR1 = 0;
    stopSampIdx_IR2 = 0;
    stopSampIdx_GR1 = 0;

    fac_attack_GR1 = 1;
    fac_decay_GR1 = 1;
    fac_release_GR1 = 1;

    startStopIdx_GR1 = 0;

    _IR1_size = 0;
    _IR2_size = 0;
    _GR1_size = 0;


    dur_attack_GR1 = 1;
    dur_decay_GR1 = 1;
    dur_release_GR1 = 1;

    startGainAttack_GR1 = 1;
    startGainDecay_GR1 = 1;       // display in dB?
    startGainRelease_GR1 = 1;

    procRepCount = 0;

    currLen_IR1 = maxLength_IR1;
    currLen_IR2 = maxLength_IR2;
    currLenRev_GR1 = 0;                // HMMMMMM....

}

Granuvolve::~Granuvolve() {
}

// VstInt32 Granuvolve::processEvents(VstEvents* ev) {
//     for (int i = 0; i < ev->numEvents; i++) {
//         if ((ev->events[i])->type != kVstMidiType) continue;
//         VstMidiEvent* event = (VstMidiEvent*)ev->events[i];
//         char* data = event->midiData;
//         VstInt32 status = data[0] & CHNMASK;  // OMNI

//         if (status == NOTEON || status == NOTEOFF) {
//             if (status == NOTEON) {

//             }
//             else {

//             }
// //         }
// //     }
// }



void Granuvolve::processReplacing(float** inputs,
 float** outputs, VstInt32 vecframes)
{
    #ifdef DEBUG_VST
    procRepCounter();

    // Display vecframes
    _vecFrames = vecframes;
    
    #endif

    float avg_samp_idx_GR1 = (startSampIdx_GR1 + endSampIdx_GR1) / 2;
    std::normal_distribution<float> normal_dist(avg_samp_idx_GR1, stdDev_GR1);
    std::uniform_int_distribution<int> uni_dist(1,1000);

    float *in1 = inputs[0],*in2 = inputs[1];
	float *out1 = outputs[0], *out2 = outputs[1];

    for (VstInt32 i = 0; i < vecframes; i++) {
        if (isMuteInput) {  // if mute and nothing is being recorded:
            if (isPlayConv)   // if play is on and the convolution processing is done,                                                
				playConv(out1, out2, i);
            else if (isPlayBufIR1) 
				playBufIR1(out1, out2, i);
			else if (isPlayBufIR2)
				playBufIR2(out1, out2, i);
			else if (isPlayBufGR1) {
				playBufGR1(out1, out2, i, normal_dist, uni_dist);
                if (isIR1Rec)
                    recBufIR1(out1, out2, i);
                if (isIR2Rec)
                    recBufIR2(out1, out2, i);
            }
            else {
                out1[i] = 0;   // Play silence. (Bypassed input)
                out2[i] = 0;
            }
        }
        else {
			//if (isGranulate) 

			

			if (isGR1Rec)
				recBufGR1(in1, in2, i);

            // Let audio input through

            out1[i] = in1[i];
            out2[i] = in2[i];
        }

        // Used to update display correctly.

        updateInfo();

    }
}

AudioEffect* createEffectInstance(audioMasterCallback audioMaster) {
	return new Granuvolve(audioMaster);
}


void Granuvolve::updateInfo() {

    if (isIR1DoneRec || isIR2DoneRec ||
        isGR1DoneRec) {
        isIR1DoneRec = false;
        isIR2DoneRec = false;
        isGR1DoneRec = false;
        updateDisplay();
    } 

	if (isIR1RecNotFully) {
		isIR1RecNotFully = false;
		updateDisplay();

	}

	if (isIR2RecNotFully) {
		isIR2RecNotFully = false;
		updateDisplay();
	}

	if (hasJustClearedIR1) {
		isClearingIR1 = false;
		hasJustClearedIR1 = false;
		updateDisplay();
	}

	if (hasJustClearedIR2) {
		isClearingIR2 = false;
		hasJustClearedIR2 = false;
		updateDisplay();
	}

    if (hasJustClearedGR1) {
        isClearingGR1 = false;
        hasJustClearedGR1 = false;
        updateDisplay();
    }

	if (isDoConvProc) {
		updateDisplay();
		isDoConvProc = false;
	}

	if (isConvProcDone) {
		isDoConvProc = false;
		isConvProcDone = false;
		isIFFTBufFull = true;
		updateDisplay();
	}

	
}

void Granuvolve::shiftStartStopIdx_GR1(float shiftIdx) {
    
    float stopIdx = shiftIdx + stopIdx_GR1-startIdx_GR1;
    float startIdx = shiftIdx;

	startSampIdx_GR1 = startIdx * _GR1_size;
	stopSampIdx_GR1 = stopIdx * _GR1_size;

	if (startSampIdx_GR1 == _GR1_size)
		startSampIdx_GR1 = _GR1_size - 1;
	else if (startSampIdx_GR1 < 0)
		startSampIdx_GR1 = 0;

	if (stopSampIdx_GR1 == _GR1_size)
		stopSampIdx_GR1 = _GR1_size - 1;
	else if (stopSampIdx_GR1 < 0)
		stopSampIdx_GR1 = 0;

	//isShiftingStartStopIdx_GR1 = true;
    //setStartIdxGR1(startIdx);
    //setStopIdxGR1(stopIdx);
	//updateDisplay();
    
}

void Granuvolve::recBufGR1(float* &in1, float* &in2, VstInt32 &i) {
	left_ch_GR1.push_back(in1[i]);
	right_ch_GR1.push_back(in2[i]);

	writeIndex_GR1++;

	// Have we reached max size?

	if (writeIndex_GR1 >= maxLength_GR1) {

		currLen_GR1 = writeIndex_GR1-1;
		stopSampIdx_GR1 = currLen_GR1;
		startIdx_GR1 = 0;
		dur_repeat_GR1 = 1;

		_GR1_size = currLen_GR1;
		writeIndex_GR1 = 0;             // Reset
		isGR1BufFull = true;            // Full buffer so...
		isClearBufGR1 = false;
		isGR1Rec = false;               // Turn off rec. of IR1

		updateDisplay();
	}
}

void Granuvolve::recBufIR1(float* &in1, float* &in2, VstInt32 &i) {
	left_ch_IR1.push_back(in1[i]);
	right_ch_IR1.push_back(in2[i]);

	writeIndex_IR1++;

	// Have we reached max size?

	if (writeIndex_IR1 >= maxLength_IR1) {

		currLen_IR1 = writeIndex_IR1-1;
		stopSampIdx_IR1 = writeIndex_IR1-1;
		_IR1_size = currLen_IR1-1;
		writeIndex_IR1 = 0;             // Reset
		isIR1BufFull = true;            // Full buffer so...
		isClearBufIR1 = false;
		isIR1Rec = false;               // Turn off rec. of IR1

		updateDisplay();
	}
}

void Granuvolve::recBufIR2(float* &in1, float* &in2, VstInt32 &i) {
	left_ch_IR2.push_back(in1[i]);
	right_ch_IR2.push_back(in2[i]);

	writeIndex_IR2++;

	// Have we reached max size?

	if (writeIndex_IR2 == maxLength_IR2) {

		currLen_IR2 = writeIndex_IR2;
		stopSampIdx_IR2 = writeIndex_IR2;
		writeIndex_IR2 = 0;             // Reset
		isIR2BufFull = true;            // Full buffer so...
		isClearBufIR2 = false;
		isIR2Rec = false;               // Turn off rec. of IR2

		updateDisplay();
	}
}

void Granuvolve::playBufGR1(float* &out1, float* &out2, VstInt32 &i, 
    std::normal_distribution<float> &normal_dist,
    std::uniform_int_distribution<int> &uni_dist) {
	if (isADSR) {
        if (isNoteOn_GR1) {  // replace with NOTEON. Connect to midi
			if (isNoteOnFirstTime) {
				calcAttackFac_GR1();
				calcDecayFac_GR1();
				isNoteOnFirstTime = false;
				sampGrainIdx = 0;

			}
            if (sampGrainIdx < dur_sampAttack_GR1) {
                gain_mult = adsr_offset1 + adsr_start1 * adsr_valrange1;
                adsr_start1 *= fac_attack_GR1;
            }

            else if (sampGrainIdx < dur_sampAttack_GR1 + dur_sampDecay_GR1) {
                gain_mult = adsr_offset2 + adsr_start2 * adsr_valrange2;
                adsr_start2 *= fac_decay_GR1;
            }

            else // sustain
                gain_mult = startGainRelease_GR1;

            sampGrainIdx++;
			isNoteOffFirstTime = true;
        }

        else {
			if (isNoteOffFirstTime) {
				calcReleaseFac_GR1();
				isNoteOffFirstTime = false;
			}
            gain_mult = adsr_offset3 + adsr_start3 * adsr_valrange3;
            adsr_start3 *= fac_release_GR1;
			isNoteOnFirstTime = true;
        }
    }
    else {
        gain_mult = 1;
    }

    
	// Adding randomness to the mix.
	int dice_roll = uni_dist(generator);
	if (dice_roll  < randomProb_GR1*1000) {

		int randIdx_GR1 = (int) normal_dist(generator);

		if (isRevBufGR1) {
			if (randIdx_GR1 >= stopSampIdx_GR1)
				randIdx_GR1 = startSampIdx_GR1;
            else if (randIdx_GR1 < 0)
                randIdx_GR1 = 0;

			out1[i] = gain_mult * l_rev_GR1.at(randIdx_GR1);
			out2[i] = gain_mult * r_rev_GR1.at(randIdx_GR1);
		}

		else {
			if (randIdx_GR1 >= stopSampIdx_GR1)
				randIdx_GR1 = startSampIdx_GR1;
            else if (randIdx_GR1 < 0)
                randIdx_GR1 = 0;

			out1[i] = gain_mult * left_ch_GR1.at(randIdx_GR1);
			out2[i] = gain_mult * right_ch_GR1.at(randIdx_GR1);
		}
	}

	else {
	
		if (isRevBufGR1) {
			if (readIndex_GR1 >= stopSampIdx_GR1)
				readIndex_GR1 = startSampIdx_GR1;

			out1[i] = gain_mult * l_rev_GR1.at(readIndex_GR1);
			out2[i] = gain_mult * r_rev_GR1.at(readIndex_GR1);
		}

		else {
			if (readIndex_GR1 >= stopSampIdx_GR1)
				readIndex_GR1 = startSampIdx_GR1;

			out1[i] = gain_mult * left_ch_GR1.at(readIndex_GR1);
			out2[i] = gain_mult * right_ch_GR1.at(readIndex_GR1);
		}

		readIndex_GR1++;
	}
}

void Granuvolve::playBufIR1(float* &out1, float* &out2, VstInt32 &i) {
    if (isRevBufIR1) {
        if (readIndex_IR1 < stopSampIdx_IR1) {
            out1[i] = l_rev_IR1.at(readIndex_IR1);
            out2[i] = r_rev_IR1.at(readIndex_IR1);
        }
        else {
            readIndex_IR1 = startSampIdx_IR1;
            isPlayBufIR1 = false;
            isMuteInput = false;
            updateDisplay();
        }
    }

    else {
        if (readIndex_IR1 < stopSampIdx_IR1) {
            out1[i] = left_ch_IR1.at(readIndex_IR1);
            out2[i] = right_ch_IR1.at(readIndex_IR1);      
        }
        else { 
            readIndex_IR1 = startSampIdx_IR1;
            isPlayBufIR1 = false;
            isMuteInput = false;
            updateDisplay();
        }
    }
    
    readIndex_IR1++;
    
}

void Granuvolve::playBufIR2(float* &out1, float* &out2, VstInt32 &i) {
	if (isRevBufIR2) {
		if (readIndex_IR2 < currLenRev_IR2) {
			out1[i] = l_rev_IR2.at(readIndex_IR2);
			out2[i] = r_rev_IR2.at(readIndex_IR2);
		}
		else {
			readIndex_IR2 = 0;
			isPlayBufIR2 = false;
			isMuteInput = false;
			updateDisplay();
		}
	}

	else {
		if (readIndex_IR2 < stopSampIdx_IR2) {
			out1[i] = left_ch_IR2.at(readIndex_IR2);
			out2[i] = right_ch_IR2.at(readIndex_IR2);
		}
		else {
			readIndex_IR2 = startSampIdx_IR2;
			isPlayBufIR2 = false;
			isMuteInput = false;
			updateDisplay();
		}
	}

	readIndex_IR2++;
}

void Granuvolve::playConv(float* &out1, float* &out2, VstInt32 &i) {
    // Copy the ifft of fft prod result into the output

    if (readIndex_ifft < max_size_ifft) {
        out1[i] = left_ch_ifft.at(readIndex_ifft);
        out2[i] = right_ch_ifft.at(readIndex_ifft);
        readIndex_ifft++;

        if (isFindSil_IFFT) {
            bool prev_samp_zero = false;
            
            if (isInToleranceLevel(left_ch_ifft.at(readIndex_ifft-1)) &&
                isInToleranceLevel(right_ch_ifft.at(readIndex_ifft-1))) 
                prev_samp_zero = true;
            else 
                prev_samp_zero = false;

            if (readIndex_ifft < max_size_ifft && 
                isInToleranceLevel(left_ch_ifft.at(readIndex_ifft)) &&
                isInToleranceLevel(right_ch_ifft.at(readIndex_ifft)) &&
                prev_samp_zero)
                tailLen_IFFT++;
            else {
                tailLen_IFFT = 0;
            }

            if (tailLen_IFFT == tailLenMax_IFFT) {
                isPlayConv = false;
                tailLen_IFFT = 0;
                updateDisplay();
            }
        }
    }

    else { 
        tailLen_IFFT = 0;
        readIndex_ifft = 0;
        isPlayConv = false;
        isMuteInput = false;
        updateDisplay();
    }
}

void Granuvolve::setIsIR1Rec(float value){
	if (value > 0.5 && isIR1Rec_able()) {
		stopIdx_IR1 = 1;
		startIdx_IR1 = 0;
		startSampIdx_IR1 = 0;
		isIR1Rec = true;
        IR1RecWasJustOn = true;
	}
    else if (value < 0.5 && IR1RecWasJustOn) {
        isIR1Rec = false;
        currLen_IR1 = writeIndex_IR1;
        stopSampIdx_IR1 = writeIndex_IR1;
        _IR1_size = left_ch_IR1.size();
        writeIndex_IR1 = 0;
        isIR1BufFull = true;
        isClearBufIR1 = false;
        isIR1RecNotFully = true;
		isIR1DoneRec = true;
        IR1RecWasJustOn = false;
    }
}
void Granuvolve::setIsIR2Rec(float value){
	if (value > 0.5 && isIR2Rec_able()) {
        stopIdx_IR2 = 1;
        startIdx_IR2 = 0;
        startSampIdx_IR2 = 0;
		isIR2Rec = true;
        IR2RecWasJustOn = true;
	}
	else if (value < 0.5 && IR2RecWasJustOn) {
		isIR2Rec = false;
		currLen_IR2 = writeIndex_IR2;
        stopSampIdx_IR2 = writeIndex_IR2;
		_IR2_size = left_ch_IR2.size();
        writeIndex_IR2 = 0;
		isIR2BufFull = true;
		isClearBufIR2 = false;
		isIR2RecNotFully = true;
		isIR2DoneRec = true;
        IR2RecWasJustOn = false;
	}
}

void Granuvolve::setIsGR1Rec(float value){
	if (value > 0.5 && isGR1Rec_able()) {
		stopIdx_GR1 = 1;
		startIdx_GR1 = 0;
		startSampIdx_GR1 = 0;
		isGR1Rec = true;
        GR1RecWasJustOn = true;
	}
	else if (value < 0.5 && GR1RecWasJustOn) {
		isGR1Rec = false;
		currLen_GR1 = writeIndex_GR1;
		stopSampIdx_GR1 = writeIndex_GR1;
		dur_repeat_GR1 = 1;
		if (stopSampIdx_GR1 >= left_ch_GR1.size())
			stopSampIdx_GR1 = left_ch_GR1.size() - 1;
		_GR1_size = left_ch_GR1.size();
		writeIndex_GR1 = 0;
		isGR1BufFull = true;
		isClearBufGR1 = false;
		isGR1RecNotFully = true;
		isGR1DoneRec = true;
        GR1RecWasJustOn = false;

		// Calculate the ADSR envelope
		
		calcAttackFac_GR1();
		calcDecayFac_GR1();
		calcReleaseFac_GR1();
	}
}


void Granuvolve::setIsDoConvProc(float value){
	if (value > 0.5 && !isDoConvProc) {
		// && isDoConvProc_able()) {
		isDoConvProc = true;
		updateDisplay();
		isDoConvProc = Granuvolve::convolute();
		isConvProcDone = true;
		max_size_ifft = left_ch_ifft.size();

	}
}
void Granuvolve::setIsMuteInput(float value){
	if (value > 0.5 && isMuteInput_able())
		isMuteInput = true;
	else if (value <= 0.5 && isMuteInput_Offable())
		isMuteInput = false;
}
void Granuvolve::setIsClearBufIR1(float value){
	if (value > 0.5 && isClearBufIR1_able()) {
		isClearingIR1 = true;
		isClearBufIR1 = clear_buffer_IR1();
		hasJustClearedIR1 = true;
	}
}
void Granuvolve::setIsClearBufIR2(float value){
	if (value > 0.5 && isClearBufIR2_able()) {
		isClearingIR2 = true;
		isClearBufIR2 = clear_buffer_IR2();
		hasJustClearedIR2 = true;
	}
}
void Granuvolve::setIsClearBufGR1(float value){
    if (value > 0.5 && isClearBufGR1_able()) {
        isClearingGR1 = true;
        isClearBufGR1 = clear_buffer_GR1();
        hasJustClearedGR1 = true;
    }
}
void Granuvolve::setIsPlayConv(float value){
	if (value > 0.5 && isPlayConv_able()) {
		isPlayBufIR1 = false;
		isPlayBufIR2 = false;
		readIndex_IR1 = 0;
		readIndex_IR2 = 0;

		if (!isIR1Rec && !isIR2Rec)
			isMuteInput = true;

		updateDisplay();

		isPlayConv = true;
        playConvWasJustOn = true;
	}
	else if (value <= 0.5 && playConvWasJustOn) {
		isPlayConv = false;
		tailLen_IFFT = 0;
		readIndex_ifft = 0;
		isMuteInput = false;
        playConvWasJustOn = false;
	}    
    // else
    //     isPlayConv = false; 
}
void Granuvolve::setIsPlayBufIR1(float value){
	if (value > 0.5 && isPlayBufIR1_able()) {
		isPlayBufGR1 = false;
		isPlayBufIR2 = false;
		isPlayConv = false;
        isNoteOn_GR1 = false;

        readIndex_IR1 = startSampIdx_IR1;
		readIndex_IR2 = startSampIdx_IR2;
		readIndex_ifft = 0;
		currLen_IR1 = left_ch_IR1.size();

		if (!isIR1Rec && !isIR2Rec) // && !isGR1Rec
			isMuteInput = true;

		updateDisplay();

		isPlayBufIR1 = true;
        playBufIR1WasJustOn = true;
	}
	else  if (value <= 0.5 && playBufIR1WasJustOn) {
		isPlayBufIR1 = false;
		readIndex_IR1 = startSampIdx_IR1;

		isMuteInput = false;
        playBufIR1WasJustOn = false;
	}
    
}
void Granuvolve::setIsPlayBufIR2(float value){
	if (value > 0.5 && isPlayBufIR2_able()) {
		isPlayBufIR1 = false;
		isPlayBufGR1 = false;
		isPlayConv = false;
        isNoteOn_GR1 = false;

		readIndex_IR1 = startSampIdx_IR1;
		readIndex_IR2 = startSampIdx_IR2;
		readIndex_GR1 = startSampIdx_GR1;
		readIndex_ifft = 0;
        currLen_IR2 = left_ch_IR1.size();

		if (!isIR1Rec && !isIR2Rec)
			isMuteInput = true;

		updateDisplay();

		isPlayBufIR2 = true;
        playBufIR2WasJustOn = true;
	}
	else if (value <= 0.5 && playBufIR2WasJustOn) {
		isPlayBufIR2 = false;
		readIndex_IR2 = 0;
		isMuteInput = false;
        playBufIR2WasJustOn = false;
	}    
}

void Granuvolve::setIsPlayBufGR1(float value){
    if (value > 0.5 && isPlayBufGR1_able()) {
        isMuteInput = true;
		sampGrainIdx = 0;

        isPlayBufIR1 = false;
		isPlayBufIR2 = false;
        isPlayConv = false;
        isNoteOn_GR1 = false;

		readIndex_GR1 = startSampIdx_GR1;
		readIndex_IR1 = startSampIdx_IR1;
		readIndex_IR2 = startSampIdx_IR2;
        readIndex_ifft = 0;

        if (!isIR1Rec && !isGR1Rec)
            isMuteInput = true;

        updateDisplay();

        isPlayBufGR1 = true;
        playBufGR1WasJustOn = true;
    }
    else if (value <= 0.5 && playBufGR1WasJustOn) {
        isPlayBufGR1 = false;
        readIndex_GR1 = 0;
        isMuteInput = false;

        playBufGR1WasJustOn = false;
		
    }    
}
void Granuvolve::setIsRevBufIR1(float value){
	if (value > 0.5 && isRevBufIR1_able()) {

        // isPlayBufIR1 = false;

        if (!isRevIR1Full) {
            l_rev_IR1.clear();
            r_rev_IR1.clear();
            isRevIR1Full = false;

            unsigned int lch_size = left_ch_IR1.size();

            for (unsigned int j = lch_size - 1; j > 0; j--) {

                l_rev_IR1.push_back(left_ch_IR1.at(j));
                r_rev_IR1.push_back(right_ch_IR1.at(j));
            }

            currLenRev_IR1 = l_rev_IR1.size(); // might fail
            isRevIR1Full = true;


        }

        isRevBufIR1 = true;
        updateDisplay();

    }
    else if (value <= 0.5) {
        isRevBufIR1 = false;
    }

}
void Granuvolve::setIsRevBufIR2(float value){
	if (value > 0.5 && isRevBufIR2_able()) {

		// isPlayBufIR2 = false;

		if (!isRevIR2Full) {
			l_rev_IR2.clear();
			r_rev_IR2.clear();
			isRevIR2Full = false;

			unsigned int lch_size = left_ch_IR2.size();

			for (unsigned int j = lch_size - 1; j > 0; j--) {

				l_rev_IR2.push_back(left_ch_IR2.at(j));
				r_rev_IR2.push_back(right_ch_IR2.at(j));
			}

			currLenRev_IR2 = l_rev_IR2.size(); // might fail
			isRevIR2Full = true;


		}

		isRevBufIR2 = true;
		updateDisplay();

	}
	else if (value <= 0.5) {
		isRevBufIR2 = false;
	}
}

void Granuvolve::setIsRevBufGR1(float value){
    if (value > 0.5 && isRevBufGR1_able()) {

        // isPlayBufGR1 = false;

        if (!isRevGR1Full) {
            l_rev_GR1.clear();
            r_rev_GR1.clear();
            isRevGR1Full = false;

            unsigned int lch_size = left_ch_GR1.size();

            for (unsigned int j = lch_size - 1; j > 0; j--) {

                l_rev_GR1.push_back(left_ch_GR1.at(j));
                r_rev_GR1.push_back(right_ch_GR1.at(j));
            }

            currLenRev_GR1 = l_rev_GR1.size(); // might fail
            isRevGR1Full = true;

            
        }

		isRevBufGR1 = true;
		updateDisplay();

    }
    else if (value <= 0.5) {
        isRevBufGR1 = false;
    }
}

void Granuvolve::setIsFindSilIFFT(float value) {
	(value > 0.5) ? isFindSil_IFFT = true : isFindSil_IFFT = false;
}
void Granuvolve::setIsRevBufIFFT(float value){
    
}
void Granuvolve::setMaxLenIR1(float value){
	if (isMaxLenIR1_editable()) {
		fracLength_IR1 = value;
		currLen_IR1 = maxLength_IR1 * fracLength_IR1;
	}
}
void Granuvolve::setMaxLenIR2(float value){
	if (isMaxLenIR2_editable()) {
		fracLength_IR2 = value;
		currLen_IR2 = maxLength_IR2 * fracLength_IR2;
	}
}
void Granuvolve::setMaxLenGR1(float value){
    if (isMaxLenGR1_editable()) {
        fracLength_GR1 = value;
        currLen_GR1 = maxLength_GR1 * fracLength_GR1;
    }
}
void Granuvolve::setStartFindProcDur(float value) {
	if (value > 0.5 && !startFindProcDur)
		startFindProcDur = true;
}
void Granuvolve::setStartIdxIR1(float value){
	if (value < stopIdx_IR1) {
		startIdx_IR1 = value;
		// Enforce boundaries.
		startSampIdx_IR1 = startIdx_IR1 * _IR1_size;
		if (startSampIdx_IR1 < 0)
			startSampIdx_IR1 = 0;
		else if (startSampIdx_IR1 >= _IR1_size)
			startSampIdx_IR1 = _IR1_size - 1;
		readIndex_IR1 = startSampIdx_IR1;
	}
}
void Granuvolve::setStopIdxIR1(float value){
	if (startIdx_IR1 < value) {
		stopIdx_IR1 = value;

		stopSampIdx_IR1 = stopIdx_IR1 * _IR1_size;
		if (stopSampIdx_IR1 < 0)
			stopSampIdx_IR1 = 0;
		else if (stopSampIdx_IR1 >= _IR1_size)
			stopSampIdx_IR1 = _IR1_size - 1;
	}

}

void Granuvolve::setStartIdxIR2(float value){
	if (value < stopIdx_IR2) {
		startIdx_IR2 = value;
		// Enforce boundaries.
		startSampIdx_IR2 = startIdx_IR2 * _IR2_size;
		if (startSampIdx_IR2 < 0)
			startSampIdx_IR2 = 0;
		else if (startSampIdx_IR2 >= _IR2_size)
			startSampIdx_IR2 = _IR2_size - 1;
		readIndex_IR2 = startSampIdx_IR2;
	}
}
void Granuvolve::setStopIdxIR2(float value){
	if (startIdx_IR2 < value) {
		stopIdx_IR2 = value;

		stopSampIdx_IR2 = stopIdx_IR2 * _IR2_size;
		if (stopSampIdx_IR2 < 0)
			stopSampIdx_IR2 = 0;
		else if (stopSampIdx_IR2 >= _IR2_size)
			stopSampIdx_IR2 = _IR2_size - 1;
	}
}

void Granuvolve::setStartIdxGR1(float value){
	
	float begSampIdx_GR1 = value * _GR1_size;
	endSampIdx_GR1 = begSampIdx_GR1 + dur_repeat_GR1 * _GR1_size;

	if (endSampIdx_GR1 < _GR1_size) {
		stopSampIdx_GR1 = endSampIdx_GR1;
		startIdx_GR1 = value;
		startSampIdx_GR1 = begSampIdx_GR1;
	}
}
void Granuvolve::setStopIdxGR1(float value){
	if (startIdx_GR1 < value) {
		stopIdx_GR1 = value;

		stopSampIdx_GR1 = stopIdx_GR1 * _GR1_size;
		if (stopSampIdx_GR1 < 0)
			stopSampIdx_GR1 = 0;
		else if (stopSampIdx_GR1 >= _GR1_size)
			stopSampIdx_GR1 = _GR1_size - 1;
	}

}


void Granuvolve::setDurRepeatGR1(float value) {
	dur_repeat_GR1 = value;
	stopSampIdx_GR1 = startSampIdx_GR1 + dur_repeat_GR1 * _GR1_size;
	if (stopSampIdx_GR1 >= _GR1_size)
		stopSampIdx_GR1 = _GR1_size - 1;
	else if (stopSampIdx_GR1 < 0)
		stopSampIdx_GR1 = 0;
}
void Granuvolve::procRepCounter() {
	if (procRepCount == 1) {
		clk = clock() - clk;
		dur_procRep = ((float) clk) / CLOCKS_PER_SEC;
		procRepCount = 0;
		startFindProcDur = false;

		updateDisplay();

	}

	if (startFindProcDur) {
		clk = clock();
		startFindProcDur = false;
		procRepCount++;
	}
}

void Granuvolve::setIsStdDev_GR1(float value) {
    isStdDev_GR1 = value > 0.5 ? true : false;
}
void Granuvolve::setStdDev_GR1(float value) {
    stdDev_GR1 = value > 0 ? value : 0.0001;
}
void Granuvolve::setRandomProb_GR1(float value) {
    randomProb_GR1 = value;
}


void Granuvolve::setParameter (VstInt32 index, float value) {
	switch(index) {
		case IS_IR1_REC:
            setIsIR1Rec(value);
            break;
        case IS_IR2_REC:
            setIsIR2Rec(value);
            break;
		case IS_GR1_REC:
			setIsGR1Rec(value);
			break;
       
        case IS_PLAY_CONV:
            setIsPlayConv(value);
            break;
        case IS_PLAY_BUF_IR1:
            setIsPlayBufIR1(value);
            break;
        case IS_PLAY_BUF_IR2:
            setIsPlayBufIR2(value);
            break;
		case IS_PLAY_BUF_GR1:
			setIsPlayBufGR1(value);
			break;
        case IS_DO_CONV_PROC:
            setIsDoConvProc(value);
            break;
        case IS_CLEAR_BUF_IR1:
            setIsClearBufIR1(value);
            break;
        case IS_CLEAR_BUF_IR2:
            setIsClearBufIR2(value);
            break;
		case IS_CLEAR_BUF_GR1:
			setIsClearBufGR1(value);
			break;
        case IS_REV_BUF_IR1:
            setIsRevBufIR1(value);
            break;
        case IS_REV_BUF_IR2:
            setIsRevBufIR2(value);
            break;
		case IS_REV_BUF_GR1:
			setIsRevBufGR1(value);
			break;
        case IS_FIND_SIL_IFFT:
            setIsFindSilIFFT(value);
            break;
		
        
        case START_IDX_IR1: 
            setStartIdxIR1(value);
            break;        
        case STOP_IDX_IR1: 
            setStopIdxIR1(value);
            break;
        case START_IDX_IR2: 
            setStartIdxIR2(value);
            break;        
        case STOP_IDX_IR2: 
            setStopIdxIR2(value);
            break;
        case START_IDX_GR1: 
            setStartIdxGR1(value);
            break;   

        case START_GAIN_ATTACK_GR1: 
            setStartGainAttack_GR1(value);
            break;       
        case START_GAIN_DECAY_GR1: 
            setStartGainDecay_GR1(value);
            break;       
        case START_GAIN_RELEASE_GR1: 
            setStartGainRelease_GR1(value);
            break;       
        case DUR_ATTACK_GR1: 
            setDurAttack_GR1(value);
            break;       
        case DUR_DECAY_GR1: 
            setDurDecay_GR1(value);
            break;       
        case DUR_RELEASE_GR1: 
            setDurRelease_GR1(value);
            break;            

        case IS_NOTE_ON_GR1:
            setIsNoteOnGR1(value);
            break;
   		case DUR_REPEAT_GR1:
            setDurRepeatGR1(value);
            break;

        
       
        case RANDOM_PROB_GR1:
            setRandomProb_GR1(value);
            break;

#ifdef DEBUG_VST

		case IS_STD_DEV_GR1:
			setIsStdDev_GR1(value);
			break;

		case MAX_LENGTH_IR1: // HOW ABOUT FOR GR1?
			setMaxLenIR1(value);
			break;
		case MAX_LENGTH_IR2:
			setMaxLenIR2(value);
			break;
		case MAX_LENGTH_GR1:
			setMaxLenGR1(value);
			break;
        case IS_MUTE_INPUT:
            setIsMuteInput(value);
            break;
		case STOP_IDX_GR1: 
		    setStopIdxGR1(value);
		    break;
		case SHIFT_START_STOP_IDX_GR1:
		    shiftStartStopIdx_GR1(value);
		    break;
        case START_FIND_PROC_DUR:
            setStartFindProcDur(value);
            break;
		case STD_DEV_GR1:
			setStdDev_GR1(value);
			break;
#endif
	}
}

float Granuvolve::getParameter (VstInt32 index) {
	switch(index) {
		case IS_IR1_REC: 
			return isIR1Rec;
			break;
        case IS_IR2_REC: 
            return isIR2Rec;
            break;
		case IS_GR1_REC:
			return isGR1Rec;
			break;
        
     
        case IS_PLAY_CONV:
            return isPlayConv;
            break;
        case IS_PLAY_BUF_IR1:
            return isPlayBufIR1;
            break;
        case IS_PLAY_BUF_IR2:
            return isPlayBufIR2;
            break;
		case IS_PLAY_BUF_GR1:
			return isPlayBufGR1;
			break;
        case IS_DO_CONV_PROC:
            return isDoConvProc;
            break;
        case IS_REV_BUF_IR1:
            return isRevBufIR1;
            break;
        case IS_REV_BUF_IR2:
            return isRevBufIR2;
            break;
		case IS_REV_BUF_GR1:
			return isRevBufGR1;
			break;
        case IS_FIND_SIL_IFFT:
            return isFindSil_IFFT;
            break;
        
        case IS_CLEAR_BUF_IR1:
            return isClearBufIR1;
            break;
        case IS_CLEAR_BUF_IR2:
            return isClearBufIR2;
            break;
		case IS_CLEAR_BUF_GR1:
			return isClearBufGR1;
			break;
      
        
        case START_IDX_IR1:
            return startIdx_IR1;
            break;
        case STOP_IDX_IR1:
            return stopIdx_IR1;
            break;
        case START_IDX_IR2:
            return startIdx_IR2;
            break;
        case STOP_IDX_IR2:
            return stopIdx_IR2;
            break;
        case START_IDX_GR1:
            return startIdx_GR1;
            break;
        case DUR_REPEAT_GR1:
            return dur_repeat_GR1;
            break;

        case START_GAIN_ATTACK_GR1:
            return startGainAttack_GR1;
            break;
        case START_GAIN_DECAY_GR1:
            return startGainDecay_GR1;
            break;
        case START_GAIN_RELEASE_GR1:
            return startGainRelease_GR1;
            break;
        case DUR_ATTACK_GR1:
            return dur_attack_GR1;
            break;
        case DUR_DECAY_GR1:
            return dur_decay_GR1;
            break;
        case DUR_RELEASE_GR1:
            return dur_release_GR1;
            break;

        case IS_NOTE_ON_GR1:
            return isNoteOn_GR1;
            break;
       
        case STD_DEV_GR1:
            return stdDev_GR1;
            break;
        case RANDOM_PROB_GR1:
            return randomProb_GR1;
            break;


#ifdef DEBUG_VST    

		case MAX_LENGTH_IR1:
			return fracLength_IR1;
			break;
		case MAX_LENGTH_IR2:
			return fracLength_IR2;
			break;
		case MAX_LENGTH_GR1:
			return fracLength_GR1;
			break;
		case IS_STD_DEV_GR1:
			return isStdDev_GR1;
			break;
        case IS_MUTE_INPUT:
            return isMuteInput;
            break;
        case POW_OF_2_SIZE:        
            return pow_of_2_size;
            break;
        case VEC_FRAMES:
            return _vecFrames;
            break;

        case START_FIND_PROC_DUR:
            return startFindProcDur;
            break;

        case IS_IR1_BUF_FULL:
            return isIR1BufFull;
            break;
        case IS_IR2_BUF_FULL:
            return isIR2BufFull;
            break;
        case IS_GR1_BUF_FULL:
            return isGR1BufFull;
            break;
        case IS_IFFT_BUF_FULL:
            return isIFFTBufFull;
            break;

        case STOP_IDX_GR1:
            return stopIdx_GR1;
            break;
		case SHIFT_START_STOP_IDX_GR1:
			return startStopIdx_GR1;
			break;
        case START_SAMP_IDX_IR1:
            return startSampIdx_IR1;
            break;
        case STOP_SAMP_IDX_IR1:
            return stopSampIdx_IR1;
            break;
        case START_SAMP_IDX_IR2:
            return startSampIdx_IR2;
            break;
        case STOP_SAMP_IDX_IR2:
            return stopSampIdx_IR2;
            break;
        case START_SAMP_IDX_GR1:
            return startSampIdx_GR1;
            break;
        case STOP_SAMP_IDX_GR1:
            return stopSampIdx_GR1;
            break;
#endif
    }

    return 0;
}

void Granuvolve::getParameterName(VstInt32 index, char* text) {
    
	switch(index) {
		case IS_IR1_REC:
			strcpy(text, "isIR1Rec");
			break;
        case IS_IR2_REC:
            strcpy(text, "isIR2Rec");
            break;
        case IS_GR1_REC:
            strcpy(text, "isGR1Rec");
            break;
        
        case IS_PLAY_CONV:
            strcpy(text, "isPlayConv");
            break;
        case IS_PLAY_BUF_IR1:
            strcpy(text, "isPlayBufIR1");
            break;
        case IS_PLAY_BUF_IR2:
            strcpy(text, "isPlayBufIR2");
            break;
        case IS_PLAY_BUF_GR1:
            strcpy(text, "isPlayBufGR1");
            break;
        case IS_DO_CONV_PROC:
            strcpy(text, "isDoConvProc");
            break;
        case IS_CLEAR_BUF_IR1:
            strcpy(text, "isClearBufIR1");
            break;
        case IS_CLEAR_BUF_IR2:
            strcpy(text, "isClearBufIR2");
            break;
        case IS_CLEAR_BUF_GR1:
            strcpy(text, "isClearBufGR1");
            break;
        case IS_REV_BUF_IR1:
            strcpy(text, "isRevBufIR1");
            break;
        case IS_REV_BUF_IR2:
            strcpy(text, "isRevBufIR2");
            break;
        case IS_REV_BUF_GR1:
            strcpy(text, "isRevBufGR1");
            break;
        case IS_FIND_SIL_IFFT:
            strcpy(text, "isFindSil_IFFT");
            break;
      
        case DUR_REPEAT_GR1:
            strcpy(text, "dur_repeat_GR1");
            break;

         case START_IDX_IR1:
            strcpy(text, "startIdx_IR1");
            break;
        case STOP_IDX_IR1:
            strcpy(text, "stopIdx_IR1");
            break;
        case START_IDX_IR2:
            strcpy(text, "startIdx_IR2");
            break;
        case STOP_IDX_IR2:
            strcpy(text, "stopIdx_IR2");
            break;
        case START_IDX_GR1:
            strcpy(text, "startIdx_GR1");
            break;

        case START_GAIN_ATTACK_GR1:
            strcpy(text, "sg_Attack_GR1");
            break;
        case START_GAIN_DECAY_GR1:
            strcpy(text, "sg_Decay_GR1");
            break;
        case START_GAIN_RELEASE_GR1:
            strcpy(text, "sg_Release_GR1");
            break;
        case DUR_ATTACK_GR1:
            strcpy(text, "dur_attack_GR1");
            break;
        case DUR_DECAY_GR1:
            strcpy(text, "dur_decay_GR1");
            break;
        case DUR_RELEASE_GR1:
            strcpy(text, "dur_release_GR1");
            break;
        case IS_NOTE_ON_GR1:
            strcpy(text, "isNoteOn_GR1");
            break;

        
        case STD_DEV_GR1:
            strcpy(text, "stdDev_GR1");
            break;
         case RANDOM_PROB_GR1:
            strcpy(text, "randomProb_GR1");
            break;

    #ifdef DEBUG_VST

		 case MAX_LENGTH_IR1:
			 strcpy(text, "maxLengthIR1");
			 break;
		 case MAX_LENGTH_IR2:
			 strcpy(text, "maxLengthIR2");
			 break;
		 case MAX_LENGTH_GR1:
			 strcpy(text, "maxLengthGR1");
			 break;
		 case IS_STD_DEV_GR1:
			 strcpy(text, "isStdDev_GR1");
			 break;

        case IS_MUTE_INPUT:
            strcpy(text, "isMuteInput");
            break;
        case IS_IR1_BUF_FULL:
            strcpy(text, "isIR1BufFull");
            break;
        case IS_IR2_BUF_FULL:
            strcpy(text, "isIR2BufFull");
            break;
        case IS_GR1_BUF_FULL:
            strcpy(text, "isGR1BufFull");
            break;
        case IS_IFFT_BUF_FULL:
            strcpy(text, "isIFFTBufFull");
            break;
        case DUR_FFTS:
            strcpy(text, "dur_ffts");
            break;
        case DUR_FREQ_MULT:
            strcpy(text, "dur_freq_mult");
            break;
        case DUR_IFFTS:
            strcpy(text, "dur_iffts");
            break;
        case POW_OF_2_SIZE:
            strcpy(text, "pow_of_2_size");
            break;
        case DUR_CONV_PROC:
            strcpy(text, "dur_convProc");
            break;
        case DUR_PROC_REP:
            strcpy(text, "dur_procRep");
            break;
        case VEC_FRAMES:
            strcpy(text, "vecFrames");
            break;
        case START_FIND_PROC_DUR:
            strcpy(text, "startFindProcDur");
            break;

       case STOP_IDX_GR1:
            strcpy(text, "stopIdx_GR1");
            break;

        case SHIFT_START_STOP_IDX_GR1:
            strcpy(text, "shiftStartStopIdx_GR1");
            break;
        case START_SAMP_IDX_IR1:
            strcpy(text, "startSampIdx_IR1");
            break;
        case STOP_SAMP_IDX_IR1:
            strcpy(text, "stopSampIdx_IR1");
            break;
        case START_SAMP_IDX_IR2:
            strcpy(text, "startSampIdx_IR2");
            break;
        case STOP_SAMP_IDX_IR2:
            strcpy(text, "stopSampIdx_IR2");
            break;
        case START_SAMP_IDX_GR1:
            strcpy(text, "startSampIdx_GR1");
            break;
        case STOP_SAMP_IDX_GR1:
            strcpy(text, "stopSampIdx_GR1");
            break;
    #endif
       
 
	}
}

void Granuvolve::getParameterLabel(VstInt32 index, char* label) {

	switch(index) {
		
      
        case DUR_REPEAT_GR1:
            strcpy(label, "sec");
            break;

        case START_IDX_IR1: 
            strcpy(label, "sec");
            break;
        case STOP_IDX_IR1: 
            strcpy(label, "sec");
            break;
        case START_IDX_IR2: 
            strcpy(label, "sec");
            break;
        case STOP_IDX_IR2: 
            strcpy(label, "sec");
            break;
        case START_IDX_GR1: 
            strcpy(label, "sec");
            break;

        case START_GAIN_ATTACK_GR1:
            strcpy(label, "gain");
            break;
        case START_GAIN_DECAY_GR1:
            strcpy(label, "gain");
            break;
        case START_GAIN_RELEASE_GR1:
            strcpy(label, "gain");
            break;
        case DUR_ATTACK_GR1:
            strcpy(label, "sec");
            break;
        case DUR_DECAY_GR1:
            strcpy(label, "sec");
            break;
        case DUR_RELEASE_GR1:
            strcpy(label, "sec");
            break;

#ifdef DEBUG_VST

		case MAX_LENGTH_IR1:
			strcpy(label, "sec");
			break;
		case MAX_LENGTH_IR2:
			strcpy(label, "sec");
			break;
		case MAX_LENGTH_GR1:
			strcpy(label, "sec");
			break;

        case VEC_FRAMES:
            strcpy(label, "frames");
            break;
        case POW_OF_2_SIZE:
            strcpy(label, "frames");
            break;
        case DUR_CONV_PROC:
            strcpy(label, "sec");
            break;
        case DUR_FFTS:
            strcpy(label, "sec");
            break;
        case DUR_FREQ_MULT:
            strcpy(label, "sec");
            break;
        case DUR_IFFTS:
            strcpy(label, "sec");
            break;
        case DUR_PROC_REP:
            strcpy(label, "sec");
            break;

        case STOP_IDX_GR1: 
            strcpy(label, "sec");
            break;

        case START_SAMP_IDX_IR1:
            strcpy(label, "samples");
            break;
        case STOP_SAMP_IDX_IR1:
            strcpy(label, "samples");
            break;
        case START_SAMP_IDX_IR2:
            strcpy(label, "samples");
            break;
        case STOP_SAMP_IDX_IR2:
            strcpy(label, "samples");
            break;
        case START_SAMP_IDX_GR1:
            strcpy(label, "samples");
            break;
        case STOP_SAMP_IDX_GR1:
            strcpy(label, "samples");
            break;
#endif
        
    
	}
}

void Granuvolve::getParameterDisplay (VstInt32 index, char* text) {

	switch(index) {
		case IS_IR1_REC:
            if (isIR1Rec)  
                strcpy(text, "on");
            else 
                strcpy(text, "off");
			break;
		case IS_IR2_REC:
            if (isIR2Rec)  
                strcpy(text, "on");
            else 
                strcpy(text, "off"); 
            break;
        case IS_GR1_REC:
            if (isGR1Rec)  
                strcpy(text, "on");
            else 
                strcpy(text, "off"); 
            break;
        case IS_DO_CONV_PROC:
            if (isDoConvProc)
                strcpy(text, "on");
            else
                strcpy(text, "off");
            break;
        case IS_PLAY_CONV:
            if (isPlayConv)  
                strcpy(text, "on");
            else 
                strcpy(text, "off"); 
            break;
        case IS_PLAY_BUF_IR1:
            if (isPlayBufIR1)
                strcpy(text, "on");
            else
                strcpy(text, "off");
            break;
        case IS_PLAY_BUF_IR2:
            if (isPlayBufIR2)
                strcpy(text, "on");
            else
                strcpy(text, "off");
            break;
        case IS_PLAY_BUF_GR1:
            if (isPlayBufGR1)
                strcpy(text, "on");
            else
                strcpy(text, "off");
            break;
        case IS_CLEAR_BUF_IR1:
            if (isClearBufIR1)
                strcpy(text, "on");
            else 
                strcpy(text, "off");
            break;
        case IS_CLEAR_BUF_IR2:
            if (isClearBufIR2)
                strcpy(text, "on");
            else 
                strcpy(text, "off");
            break;
        case IS_CLEAR_BUF_GR1:
            if (isClearBufGR1)
                strcpy(text, "on");
            else 
                strcpy(text, "off");
            break;
        case IS_REV_BUF_IR1:
            if (isRevBufIR1)
                strcpy(text, "on");
            else 
                strcpy(text, "off");
            break;
        case IS_REV_BUF_IR2:
            if (isRevBufIR2)
                strcpy(text, "on");
            else 
                strcpy(text, "off");
            break;
        case IS_REV_BUF_GR1:
            if (isRevBufGR1)
                strcpy(text, "on");
            else 
                strcpy(text, "off");
            break;
        case IS_FIND_SIL_IFFT:
            if (isFindSil_IFFT) 
                strcpy(text, "on");
            else 
                strcpy(text, "off");
            break;
        

        case START_IDX_IR1: 
            ms2string(startIdx_IR1 * left_ch_IR1.size() / 1000,
                text,
                kVstMaxParamStrLen);
            break;
        case STOP_IDX_IR1: 
            ms2string(stopIdx_IR1 * left_ch_IR1.size() / 1000,
                text,
                kVstMaxParamStrLen);
            break;
        case START_IDX_IR2: 
            ms2string(startIdx_IR2 * left_ch_IR2.size() / 1000,
                text,
                kVstMaxParamStrLen);
            break;
        case STOP_IDX_IR2: 
            ms2string(stopIdx_IR2 * left_ch_IR2.size() / 1000,
                text,
                kVstMaxParamStrLen);
            break;
        case START_IDX_GR1: 
            ms2string(startIdx_GR1 * left_ch_GR1.size() / 1000,
                text,
                kVstMaxParamStrLen);
            break;
        case DUR_REPEAT_GR1:
            ms2string(dur_repeat_GR1 * left_ch_GR1.size() / 1000,
                text,
                kVstMaxParamStrLen);
            break;

        case START_GAIN_ATTACK_GR1:
            float2string(startGainAttack_GR1,
                text,
                kVstMaxParamStrLen);
            break;
        case START_GAIN_DECAY_GR1:
            float2string(startGainDecay_GR1,
                text,
                kVstMaxParamStrLen);
            break;
         case START_GAIN_RELEASE_GR1:
            float2string(startGainRelease_GR1,
                text,
                kVstMaxParamStrLen);
            break;


        case DUR_ATTACK_GR1:
            ms2string(dur_attack_GR1 * dur_sampAttackMax_GR1 / 1000,
                text,
                kVstMaxParamStrLen);
            break;
        case DUR_DECAY_GR1:
            ms2string(dur_decay_GR1 * dur_sampDecayMax_GR1 / 1000,
                text,
                kVstMaxParamStrLen);
            break;
        case DUR_RELEASE_GR1:
            ms2string(dur_release_GR1 * dur_sampReleaseMax_GR1 / 1000,
                text,
                kVstMaxParamStrLen);
            break;

        case IS_NOTE_ON_GR1:
            if (isNoteOn_GR1)
                strcpy(text, "on");
            else
                strcpy(text, "off");
            break;
        
        case STD_DEV_GR1:
            float2string(stdDev_GR1, text, kVstMaxParamStrLen);
            break;
        case RANDOM_PROB_GR1:
            float2string(randomProb_GR1, text, kVstMaxParamStrLen);
            break;


#ifdef DEBUG_VST

		case MAX_LENGTH_IR1:
			ms2string(maxLength_IR1 * fracLength_IR1 / 1000,
				text,
				kVstMaxParamStrLen);
			break;
		case MAX_LENGTH_IR2:
			ms2string(maxLength_IR2 * fracLength_IR2 / 1000, 
				text, 
				kVstMaxParamStrLen);
			break;
		case MAX_LENGTH_GR1:
			ms2string(maxLength_GR1 * fracLength_GR1 / 1000, 
				text, 
				kVstMaxParamStrLen);
			break;
		case IS_STD_DEV_GR1:
			if (isStdDev_GR1)
				strcpy(text, "on");
			else 
				strcpy(text, "off");
			break;

        case START_FIND_PROC_DUR:
            if (startFindProcDur)
                strcpy(text, "on");
            else 
                strcpy(text, "off");
            break;
        case DUR_FFTS:
            float2string(dur_ffts, text, kVstMaxParamStrLen);
            break;
        case DUR_FREQ_MULT:
            float2string(dur_freq_mult, text, kVstMaxParamStrLen);
            break;
        case DUR_IFFTS:
            float2string(dur_iffts, text, kVstMaxParamStrLen);
            break;
        case DUR_PROC_REP:
            float2string(dur_procRep, text, kVstMaxParamStrLen);
            break;
        case POW_OF_2_SIZE:
            int2string(pow_of_2_size, text, kVstMaxParamStrLen);
            break;
        case DUR_CONV_PROC:
            float2string(dur_convProc, text, kVstMaxParamStrLen);
            break;
        case VEC_FRAMES:
            int2string(_vecFrames, text, kVstMaxParamStrLen);
            break;

        case IS_IR1_BUF_FULL:
            if (isIR1BufFull)  
                strcpy(text, "on");
            else 
                strcpy(text, "off"); 
            break;
        case IS_IR2_BUF_FULL:
            if (isIR2BufFull)  
                strcpy(text, "on");
            else 
                strcpy(text, "off"); 
            break;
        case IS_GR1_BUF_FULL:
            if (isGR1BufFull)  
                strcpy(text, "on");
            else 
                strcpy(text, "off"); 
            break;
        case IS_IFFT_BUF_FULL:
            if (isIFFTBufFull)  
                strcpy(text, "on");
            else 
                strcpy(text, "off"); 
            break;
   
        case IS_MUTE_INPUT:
            if (isMuteInput)  
                strcpy(text, "on");
            else 
                strcpy(text, "off"); 
            break;
        case STOP_IDX_GR1: 
            ms2string(stopIdx_GR1 * left_ch_GR1.size() / 1000,
                text,
                kVstMaxParamStrLen);
            break;
        case START_SAMP_IDX_IR1:
            int2string(startSampIdx_IR1, text, kVstMaxParamStrLen);
            break;
        case STOP_SAMP_IDX_IR1:
            int2string(stopSampIdx_IR1, text, kVstMaxParamStrLen);
            break;
        case START_SAMP_IDX_IR2:
            int2string(startSampIdx_IR2, text, kVstMaxParamStrLen);
            break;
        case STOP_SAMP_IDX_IR2:
            int2string(stopSampIdx_IR2, text, kVstMaxParamStrLen);
            break;
        case START_SAMP_IDX_GR1:
            int2string(startSampIdx_GR1, text, kVstMaxParamStrLen);
            break;
        case STOP_SAMP_IDX_GR1:
            int2string(stopSampIdx_GR1, text, kVstMaxParamStrLen);
            break;
		case SHIFT_START_STOP_IDX_GR1:
            ms2string(startIdx_GR1 * left_ch_GR1.size() / 1000, 
				text, 
				kVstMaxParamStrLen);
            break;      
#endif
	}
}



bool Granuvolve::convolute() {
    dur_convProc = 0;
    clock_t t0;
    t0 = clock();
    Granuvolve::convolve();
    t0 = clock() - t0;
    dur_convProc = ((float) t0)/CLOCKS_PER_SEC;

    isIFFTBufFull = true;

    return false;
}

bool Granuvolve::convolve() {
    

    dur_ffts = 0;
    dur_freq_mult = 0;
    dur_iffts = 0;

    updateDisplay();

    float** dbl_ptr_IR1 = new float*[WIDTH];
    float** dbl_ptr_IR2 = new float*[WIDTH];

    unsigned int IR1_size = stopSampIdx_IR1 - startSampIdx_IR1;
    unsigned int IR2_size = stopSampIdx_IR2 - startSampIdx_IR1;
    unsigned int IR1_rev_size_total = l_rev_IR1.size();
    unsigned int IR2_rev_size_total = l_rev_IR2.size();

    unsigned int idx;

    for (unsigned int i = 0; i < WIDTH; i++) {
        dbl_ptr_IR1[i] = new float[IR1_size];
        dbl_ptr_IR2[i] = new float[IR2_size];
    }
    
    if (isRevBufIR1) {
        for (unsigned int i = 0; i < IR1_size; i++) {
            idx = (i + startSampIdx_IR1) % IR1_rev_size_total;
            dbl_ptr_IR1[0][i] = l_rev_IR1.at(idx);
            dbl_ptr_IR1[1][i] = r_rev_IR1.at(idx);
        }
    }
    else {
        for (unsigned int i = 0; i < IR1_size; i++) {
            dbl_ptr_IR1[0][i] = left_ch_IR1.at(i + startSampIdx_IR1);
            dbl_ptr_IR1[1][i] = right_ch_IR1.at(i + startSampIdx_IR1);
        }
    }

    if (isRevBufIR2) {
        for (unsigned int i = 0; i < IR2_size; i++) {
            idx = i + startSampIdx_IR2 % IR2_rev_size_total;
            dbl_ptr_IR2[0][i] = l_rev_IR2.at(idx);
            dbl_ptr_IR2[1][i] = r_rev_IR2.at(idx);
        }
    }
    else {
        for (unsigned int i = 0; i < IR2_size; i++) {
            dbl_ptr_IR2[0][i] = left_ch_IR2.at(i + startSampIdx_IR1);
            dbl_ptr_IR2[1][i] = right_ch_IR2.at(i + startSampIdx_IR1);
       } 
    }
    
    pow_of_2_size = 2;
    unsigned int bigger = (IR1_size > IR2_size ? IR1_size : IR2_size);

    while (pow_of_2_size < bigger) 
        pow_of_2_size *= 2;

    // Multiply again by 2 to ensure that pow_of_2_size is
    // at least twice the size of the bigger.

    pow_of_2_size *= 2;

    float* l_IR1 = dbl_ptr_IR1[0];
    float* r_IR1 = dbl_ptr_IR1[1];

    float* l_IR2 = dbl_ptr_IR2[0];
    float* r_IR2 = dbl_ptr_IR2[1];

    // Sav here just in case we want to iterate to inspect for
    // debugg purposes.

    float** l_IR1_2d;
    float** r_IR1_2d;
    float** l_IR2_2d;
    float** r_IR2_2d;

    // Take the fft.

    clock_t t1;
    t1 = clock();

    float** l_IR1_cplx = Conv::fft(l_IR1, l_IR1_2d, IR1_size, pow_of_2_size); 
    float** r_IR1_cplx = Conv::fft(r_IR1, r_IR1_2d, IR1_size, pow_of_2_size);
    float** l_IR2_cplx = Conv::fft(l_IR2, l_IR2_2d, IR2_size, pow_of_2_size);
    float** r_IR2_cplx = Conv::fft(r_IR2, r_IR2_2d, IR2_size, pow_of_2_size);

    t1 = clock() - t1;
    dur_ffts = ((float)t1)/CLOCKS_PER_SEC;
  
    // Point-wise multiplication.

    clock_t t2;
    t2 = clock();

    float** l_fft_prod = Conv::fft_mult(l_IR1_cplx, l_IR2_cplx, pow_of_2_size);
    float** r_fft_prod = Conv::fft_mult(r_IR1_cplx, r_IR2_cplx, pow_of_2_size);


    t2 = clock() - t2;
    dur_freq_mult = ((float)t2)/CLOCKS_PER_SEC;
 
    // Take the ifft:

    unsigned int width = 2;
    float** dbl_ptr_conv_sig = new float*[width];

    // Left & right channel

    // float t5 = std::chrono::high_resolution_clock::now();
    clock_t t3;
    t3 = clock();

    dbl_ptr_conv_sig[0] = Conv::ifft(l_fft_prod, pow_of_2_size);
    dbl_ptr_conv_sig[1] = Conv::ifft(r_fft_prod, pow_of_2_size);

    t3 = clock() - t3;

    dur_iffts = ((float)t3)/CLOCKS_PER_SEC;
  
    left_ch_ifft.resize(pow_of_2_size, 0);
    right_ch_ifft.resize(pow_of_2_size, 0);

    // Find biggest absolute value

    float max_abs_left = 0;
    float max_abs_right = 0;

    for(unsigned int i = 0; i < pow_of_2_size; i++) {
        float val_l = abs(dbl_ptr_conv_sig[0][i]);
        if (max_abs_left < val_l) 
            max_abs_left = val_l;

        float val_r = abs(dbl_ptr_conv_sig[0][i]);
        if (max_abs_right < val_r) 
            max_abs_right = val_r;
    }

    max_abs_left+=0.001;
    max_abs_right+=0.001;
    
    // Normalize

    for(unsigned int i = 0; i < pow_of_2_size; i++) {
        left_ch_ifft.at(i) = dbl_ptr_conv_sig[0][i] / 
                              max_abs_left;
        right_ch_ifft.at(i) = dbl_ptr_conv_sig[1][i] /
                              max_abs_right;
    }

  
    // cleanup:
    // Don't need to delete l_IR1, r_IR1, l_IR2, r_IR2
    // Need to delete l_IR1_2d, etc.

    for (unsigned int i = 0; i < width; i++) {
        delete[] l_IR1_2d[i];
        delete[] r_IR1_2d[i];
        delete[] l_IR2_2d[i];
        delete[] r_IR2_2d[i];
        delete[] l_IR1_cplx[i];
        delete[] r_IR1_cplx[i];
        delete[] l_IR2_cplx[i];
        delete[] r_IR2_cplx[i];
        delete[] l_fft_prod[i];
        delete[] r_fft_prod[i];

        delete[] dbl_ptr_IR1[i];
        delete[] dbl_ptr_IR2[i];
        delete[] dbl_ptr_conv_sig[i];

    }

    for (unsigned int i = 0; i < width; i++) {
        l_IR1_2d[i] = 0;
        r_IR1_2d[i] = 0;
        l_IR2_2d[i] = 0;
        r_IR2_2d[i] = 0;
        l_IR1_cplx[i] = 0;
        r_IR1_cplx[i] = 0;
        l_IR2_cplx[i] = 0;
        r_IR2_cplx[i] = 0;
        l_fft_prod[i] = 0;
        r_fft_prod[i] = 0;

        dbl_ptr_IR1[i] = 0;
        dbl_ptr_IR2[i] = 0;
        dbl_ptr_conv_sig[i] = 0;

    }
    
    delete[] l_IR1_2d;
    delete[] r_IR1_2d;
    delete[] l_IR2_2d;
    delete[] r_IR2_2d;
    delete[] l_IR1_cplx;
    delete[] r_IR1_cplx;
    delete[] l_IR2_cplx;
    delete[] r_IR2_cplx;
    delete[] l_fft_prod;
    delete[] r_fft_prod;

    delete[] dbl_ptr_IR1;
    delete[] dbl_ptr_IR2;
    delete[] dbl_ptr_conv_sig;
        
    l_IR1_2d = 0;
    r_IR1_2d = 0;
    l_IR2_2d = 0;
    r_IR2_2d = 0;
    l_IR1_cplx = 0;
    r_IR1_cplx = 0;
    l_IR2_cplx = 0;
    r_IR2_cplx = 0;
    l_fft_prod = 0;
    r_fft_prod = 0;

    dbl_ptr_IR1 = 0;
    dbl_ptr_IR2 = 0;
    dbl_ptr_conv_sig = 0;



    return false;
}


bool Granuvolve::clear_buffer_IR1() {

    left_ch_IR1.clear();
    right_ch_IR1.clear();

    l_rev_IR1.clear();
    r_rev_IR1.clear();

    currLen_IR1 = maxLength_IR1 * fracLength_IR1;
    startIdx_IR1 = 0;
    stopIdx_IR1 = 0;
    startSampIdx_IR1 = 0;
    stopSampIdx_IR1 = 0;
    readIndex_IR1 = 0;

    isIR1BufFull = false;
    writeIndex_IR1 = 0;

    updateDisplay();

    return true;
}


bool Granuvolve::clear_buffer_IR2() {

    left_ch_IR2.clear();
    right_ch_IR2.clear();

    l_rev_IR2.clear();
    r_rev_IR2.clear();

    currLen_IR2 = maxLength_IR2 * fracLength_IR2;
    startIdx_IR2 = 0;
    stopIdx_IR2 = 0;
    startSampIdx_IR2 = 0;
    stopSampIdx_IR2 = 0;
    readIndex_IR2 = 0;

    isIR2BufFull = false;
    writeIndex_IR2 = 0;
    updateDisplay();

    return true;
}

bool Granuvolve::clear_buffer_GR1() {
  
    left_ch_GR1.clear();
    right_ch_GR1.clear();

    l_rev_GR1.clear();
    r_rev_GR1.clear();

    currLen_GR1 = maxLength_GR1 * fracLength_GR1;
    
    isGR1BufFull = false;
    writeIndex_GR1 = 0;
    
    updateDisplay();
    return true;
}



bool Granuvolve::isIR1Rec_able() {
    return (!isPlayBufIR1 && 
        isClearBufIR1 && !isIR1Rec && 
        !isDoConvProc && !isIR1BufFull);
}

// ******
bool Granuvolve::isIR2Rec_able() {
    return (!isPlayBufIR2 && 
        isClearBufIR2 && !isIR2Rec &&
        !isDoConvProc && !isIR2BufFull);
}

bool Granuvolve::isGR1Rec_able() {
	return (!isMuteInput && !isPlayBufGR1 &&
		isClearBufGR1 && !isGR1Rec &&
		!isDoConvProc && !isGR1BufFull);
}


bool Granuvolve::isDoConvProc_able() {
    return (!isClearingIR1 && !isClearingIR2 && 
        isIR1BufFull && isIR2BufFull && 
        !isPlayConv && !isIR1Rec && 
        !isIR2Rec);
}

bool Granuvolve::isMuteInput_able() {
    return (!isMuteInput && !isPlayBufIR1 && 
        !isPlayBufIR2 && !isPlayConv);
}

bool Granuvolve::isMuteInput_Offable() {
    return (isMuteInput && !isIR1Rec && 
        !isIR2Rec);
}

bool Granuvolve::isClearBufIR1_able() {
    return (!isPlayBufIR1 &&
        !isDoConvProc);
}

bool Granuvolve::isClearBufIR2_able() {
    return (!isPlayBufIR2 &&
        !isDoConvProc);
}

bool Granuvolve::isClearBufGR1_able() {
    return (!isPlayBufGR1 &&
        !isDoConvProc);
}

bool Granuvolve::isPlayConv_able() {
    return (!isIR1Rec && !isIR2Rec &&
        isIFFTBufFull && !playConvWasJustOn);
}

bool Granuvolve::isPlayBufIR1_able() {
    return (!isIR1Rec && !isIR2Rec && 
        isIR1BufFull && !playBufIR1WasJustOn);
}

bool Granuvolve::isPlayBufIR2_able() {
    return (!isIR1Rec && !isIR2Rec &&
        isIR2BufFull && !playBufIR2WasJustOn);
}

bool Granuvolve::isPlayBufGR1_able() {
    return (!isIR1Rec && !isGR1Rec &&
        isGR1BufFull && !playBufGR1WasJustOn);
}

bool Granuvolve::isRevBufIR1_able() {
    return (isIR1BufFull && !isRevBufIR1);
}

bool Granuvolve::isRevBufIR2_able() {
    return (isIR2BufFull && !isRevBufIR2);
}

bool Granuvolve::isRevBufGR1_able() {
	return (isGR1BufFull && !isRevBufGR1);
}

bool Granuvolve::isRevBufIFFT_able() {
    return (isIFFTBufFull);
}


// How about when we are doing convolution?

bool Granuvolve::isMaxLenIR1_editable() {
    return (!isIR1Rec && !isPlayBufIR1 &&
        !isIR1BufFull);
}

bool Granuvolve::isMaxLenIR2_editable() {
    return (!isIR2Rec && !isPlayBufIR2 &&
        !isIR2BufFull);
}

bool Granuvolve::isMaxLenGR1_editable() {
    return (!isGR1Rec && !isPlayBufGR1 &&
        !isGR1BufFull);
}

bool Granuvolve::isInToleranceLevel(float num) {
    return (abs(num) < tolerance);
}




void Granuvolve::setIsNoteOnGR1(float value) {
	if (value > 0.5 && !noteOnGR1WasJustOn) {
		isNoteOn_GR1 = true;
		noteOnGR1WasJustOn = true;
	}
	else if (value <= 0.5 && noteOnGR1WasJustOn){
		isNoteOn_GR1 = false;
		noteOnGR1WasJustOn = false;
	}
}

// Set ADR durations

void Granuvolve::setDurAttack_GR1(float value) {
	dur_attack_GR1 = value;
    dur_sampAttack_GR1 = value * dur_sampAttackMax_GR1;
    // calcAttackFac_GR1();
}

void Granuvolve::setDurDecay_GR1(float value) {
	dur_decay_GR1 = value;
    dur_sampDecay_GR1 = value * dur_sampDecayMax_GR1;
    // calcDecayFac_GR1();
}

void Granuvolve::setDurRelease_GR1(float value) {
	dur_release_GR1 = value;
    dur_sampRelease_GR1 = value * dur_sampReleaseMax_GR1;
    // calcReleaseFac_GR1();
}

void Granuvolve::calcAttackFac_GR1() {

    adsr_valrange1 = startGainDecay_GR1 - startGainAttack_GR1;

    if (startGainAttack_GR1 > startGainDecay_GR1) {
        adsr_start1 = 1.0;
        adsr_end1 = verysmall;
        adsr_valrange1 = -adsr_valrange1;
        adsr_offset1 = startGainDecay_GR1;

    }

    else {
        adsr_start1 = verysmall;
        adsr_end1 = 1.0;
        adsr_offset1 = startGainAttack_GR1;
    }

    fac_attack_GR1 = pow(adsr_end1/adsr_start1, (float) 1.0/dur_sampAttack_GR1);
}

void Granuvolve::calcDecayFac_GR1() {

    adsr_valrange2 = startGainRelease_GR1 - startGainDecay_GR1;

    if (startGainDecay_GR1 > startGainRelease_GR1) {
        adsr_start2 = 1.0;
        adsr_end2 = verysmall;
        adsr_valrange2 = -adsr_valrange2;
        adsr_offset2 = startGainRelease_GR1;

    }

    else {
        adsr_start2 = verysmall;
        adsr_end2 = 1.0;
        adsr_offset2 = startGainDecay_GR1;
    }

    fac_decay_GR1 = pow(adsr_end2/adsr_start2, (float) 1.0/dur_sampDecay_GR1);
}

void Granuvolve::calcReleaseFac_GR1() {

    adsr_valrange3 = startGainRelease_GR1 - verysmall;


    if (startGainDecay_GR1 > verysmall) {
        adsr_start3 = 1.0;
        adsr_end3 = verysmall;
        adsr_valrange3 = -adsr_valrange3;
        adsr_offset3 = verysmall;

    }

    else {
        adsr_start3 = verysmall;
        adsr_end3 = 1.0;
        adsr_offset3 = startGainRelease_GR1;
    }

    fac_release_GR1 = pow(adsr_end3/adsr_start3, (float) 1.0/dur_sampRelease_GR1);
}

void calcNormalDist_GR1() {

}

// Set ADR heights

void Granuvolve::setStartGainAttack_GR1(float value) {
    startGainAttack_GR1 = value;
}

void Granuvolve::setStartGainDecay_GR1(float value) {
    startGainDecay_GR1 = value;
}

void Granuvolve::setStartGainRelease_GR1(float value) {
    startGainRelease_GR1 = value;
}

