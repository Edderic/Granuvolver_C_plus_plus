/* this is the base class header file */
#include "VST3 SDK/public.sdk/source/vst2.x/audioeffectx.h"
#include "conv.h"
#include <vector>
#include <ctime>
#include <random>

// #define DEBUG_VST                // Uncommment to turn debug on.


enum {
#ifndef DEBUG_VST
    PARAMS = 30,
#endif 
    IS_IR1_REC = 0, 
    IS_IR2_REC = 1, 
    IS_GR1_REC,

    IS_DO_CONV_PROC,
    IS_PLAY_CONV,
    IS_PLAY_BUF_IR1,
    IS_PLAY_BUF_IR2,
    IS_PLAY_BUF_GR1,
     IS_NOTE_ON_GR1,

    IS_CLEAR_BUF_IR1,
    IS_CLEAR_BUF_IR2,
    IS_CLEAR_BUF_GR1,

    IS_REV_BUF_IR1,
    IS_REV_BUF_IR2,
    IS_REV_BUF_GR1,

    START_IDX_IR1,
    STOP_IDX_IR1,
    START_IDX_IR2,
    STOP_IDX_IR2,
    
    START_IDX_GR1,
    DUR_REPEAT_GR1,

    START_GAIN_ATTACK_GR1,
    START_GAIN_DECAY_GR1,
    START_GAIN_RELEASE_GR1,
    DUR_ATTACK_GR1,
    DUR_DECAY_GR1,
    DUR_RELEASE_GR1,
    
    IS_FIND_SIL_IFFT,

    RANDOM_PROB_GR1,
    
    STD_DEV_GR1,



    // Just for display/debug purposes (not settable).
#ifdef DEBUG_VST
	IS_STD_DEV_GR1,

	MAX_LENGTH_IR1,
	MAX_LENGTH_IR2,
	MAX_LENGTH_GR1,


    START_SAMP_IDX_IR1,
    STOP_SAMP_IDX_IR1,
    START_SAMP_IDX_IR2,
    STOP_SAMP_IDX_IR2,
    START_SAMP_IDX_GR1,
    STOP_SAMP_IDX_GR1,

    STOP_IDX_GR1,
    SHIFT_START_STOP_IDX_GR1,

    IS_IR1_BUF_FULL,  
    IS_IR2_BUF_FULL,  
    IS_GR1_BUF_FULL,
    IS_IFFT_BUF_FULL,

    IS_MUTE_INPUT,
    
   
    DUR_FFTS,
    DUR_FREQ_MULT,
    DUR_IFFTS,
    DUR_PROC_REP,
    DUR_CONV_PROC,

    POW_OF_2_SIZE,
    
    START_FIND_PROC_DUR,
    VEC_FRAMES,

    PARAMS = 51,

#endif

    WIDTH = 2};

class Granuvolve : public AudioEffectX
{

    std::default_random_engine generator;
   

	unsigned int pow_of_2_size;
    bool isIR1Rec;                 // Is the first impulse response (IR) being recorded?
    bool isIR2Rec;                 // Is the second IR being recorded?
    bool isGR1Rec;                 // is the signal to be granulated being recorded?

	bool IR1RecWasJustOn;

    bool IR2RecWasJustOn;
    bool GR1RecWasJustOn;
    bool noteOnGR1WasJustOn;
    bool playConvWasJustOn;
    bool playBufIR1WasJustOn;
    bool playBufIR2WasJustOn;
    bool playBufGR1WasJustOn;
    
    bool isIR1RecNotFully;         // Did we record up to max length?
    bool isIR2RecNotFully; 
    bool isGR1RecNotFully;

    bool isIR1DoneRec;
    bool isIR2DoneRec;
    bool isGR1DoneRec;

    bool isIR1BufFull;                   // Is the first IR fully recorded (BufFull)?
    bool isIR2BufFull;                   // Is the second IR BufFull?
    bool isGR1BufFull;
    bool isIFFTBufFull;
    
    bool isMuteInput;                    // Lets audio input pass through if FALSE
    bool isDoConvProc;                   // Are we doing the convolution process?
    bool isConvProcDone;
    
    
    bool isPlayConv;                // Is the "Play Convolution" triggered?
    bool isPlayBufIR1;
    bool isPlayBufIR2;
    bool isPlayBufGR1;
    
    bool isRevBufIR1;
    bool isRevBufIR2;
    bool isRevBufGR1;

    bool isRevIR1Full;
    bool isRevIR2Full;
    bool isRevGR1Full;

    bool isFindSil_IFFT;

    bool isClearBufIR1;  
    bool isClearBufIR2;
    bool isClearBufGR1;
    
    bool isClearingIR1;
    bool isClearingIR2;
    bool isClearingGR1;

    bool isGranulate;
    bool isADSR;
    bool isStdDev_GR1;
    bool isNoteOn_GR1;

    bool isShiftingStartStopIdx_GR1;
    
    bool hasJustClearedIR1;
    bool hasJustClearedIR2;
    bool hasJustClearedGR1;

    bool startFindProcDur;

    bool isNoteOffFirstTime;
    bool isNoteOnFirstTime;

    float fracLength_IR1;             // What % of max_length_IR1 will be used 
    float fracLength_IR2;              // What % of max_length_IR2 will be used
    float fracLength_GR1;

    float dur_repeat_GR1;

    float dur_ffts;                    // Benchmark duration of FFTs 
    float dur_freq_mult;               // Benchmark duration of point-wise multiply
    float dur_iffts;
    float dur_procRep;    

    float adsr_start1;
    float adsr_start2;
    float adsr_start3;

    float adsr_end1;
    float adsr_end2;
    float adsr_end3;

    float adsr_valrange1;
    float adsr_valrange2;
    float adsr_valrange3;

    float adsr_offset1;
    float adsr_offset2;
    float adsr_offset3;


    float stdDev_GR1;

	float verysmall;


                // Benchmark duration of IFFTs
    float dur_convProc;

    float tolerance;
    clock_t clk;
   

    vector<float> l_rev_IR1;
    vector<float> r_rev_IR1;
    vector<float> l_rev_IR2;
    vector<float> r_rev_IR2;
    vector<float> l_rev_GR1;
    vector<float> r_rev_GR1;

    vector<float> left_ch_IR1;
    vector<float> right_ch_IR1;
    vector<float> left_ch_IR2;
    vector<float> right_ch_IR2;
    vector<float> left_ch_GR1;
    vector<float> right_ch_GR1;


    vector<float> left_ch_ifft;
    vector<float> right_ch_ifft;

    VstInt32 currLenRev_IR1;
    VstInt32 currLenRev_IR2;
    VstInt32 currLenRev_GR1;


    VstInt32 writeIndex_IR1;            // Keeps track of index for IR1
    VstInt32 writeIndex_IR2;            // Keeps track of index for IR2.
    VstInt32 writeIndex_GR1;

    VstInt32 readIndex_IR1;
    VstInt32 readIndex_IR2;
    VstInt32 readIndex_GR1;

    VstInt32 readIndex_ifft;
    VstInt32 max_size_ifft;

    VstInt32 maxLength_IR1;            // Maximum specified length of IR1
    VstInt32 maxLength_IR2;            // Maximum specified length of IR2
    VstInt32 maxLength_GR1;

    VstInt32 currLen_IR1;              // Current specified length of IR1
    VstInt32 currLen_IR2;              // Current specified length of IR2
    VstInt32 currLen_GR1;

    VstInt32 startSampIdx_IR1;
    VstInt32 startSampIdx_IR2;
    VstInt32 startSampIdx_GR1;

    VstInt32 stopSampIdx_IR2;
    VstInt32 stopSampIdx_IR1;
    VstInt32 stopSampIdx_GR1;

    VstInt32 endSampIdx_GR1;

    VstInt32 _IR1_size;
    VstInt32 _IR2_size;
    VstInt32 _GR1_size;

    // ADSR

    VstInt32 sampAttackMax;
    VstInt32 sampDecayMax;
    VstInt32 sampSustainMax;
    VstInt32 sampReleaseMax;

    VstInt32 samp_A_total_max;
    VstInt32 samp_AD_total_max;
    // VstInt32 samp_ADS_total_
    VstInt32 sampGrainIdx;


    float slope_attack;
    float slope_decay;
    float slope_release;

    

    float gain_mult;  // ADSR multiplier

    float startIdx_IR1;
    float startIdx_IR2;
    float startIdx_GR1;

    float stopIdx_IR2;
    float stopIdx_IR1;
    float stopIdx_GR1;

    float startStopIdx_GR1;

    float dur_attack_GR1;           // dispcay_GR1 in sec.
    float dur_decay_GR1;
    float dur_release_GR1;

    float startGainAttack_GR1;
    float startGainDecay_GR1;       // display in dB?
    float startGainRelease_GR1;

    float fac_attack_GR1;
    float fac_decay_GR1;
    float fac_release_GR1;

    float randomProb_GR1;
    
    VstInt32 dur_sampAttack_GR1;
    VstInt32 dur_sampDecay_GR1;
    VstInt32 dur_sampRelease_GR1;

    VstInt32 dur_sampAttackMax_GR1;
    VstInt32 dur_sampDecayMax_GR1;
    VstInt32 dur_sampReleaseMax_GR1;

    VstInt32 tailLen_IFFT;
    VstInt32 tailLenMax_IFFT;

    VstInt32 procRepCount;
             // Total time to do FFT, multiply, do IFFT
    VstInt32 _vecFrames;                // For display purposes
    

    bool isIR1Rec_able();
    bool isIR2Rec_able();
    bool isGR1Rec_able();

    bool isDoConvProc_able();
    bool isMuteInput_able();
    bool isMuteInput_Offable();
    bool isClearBufIR1_able();
    bool isClearBufIR2_able();
    bool isClearBufGR1_able();

    bool isPlayConv_able();
    bool isPlayBufIR1_able();
    bool isPlayBufIR2_able();
    bool isPlayBufGR1_able();

    bool isRevBufIR1_able();
    bool isRevBufIR2_able();
    bool isRevBufGR1_able();

    bool isRevBufIFFT_able();
    
    bool isMaxLenIR1_editable();
    bool isMaxLenIR2_editable();
    bool isMaxLenGR1_editable();
    
    bool isInToleranceLevel(float num);
    

    void setIsIR1Rec(float value);
    void setIsIR2Rec(float value);
    void setIsGR1Rec(float value);

    void setIsDoConvProc(float value);
    void setIsMuteInput(float value);
    void setIsMuteInput_Offable();
    
    void setIsClearBufIR1(float value);
    void setIsClearBufIR2(float value);
    void setIsClearBufGR1(float value);

    void setIsPlayConv(float value);
    void setIsPlayBufIR1(float value);
    void setIsPlayBufIR2(float value);
    void setIsPlayBufGR1(float value);

    void setIsRevBufIR1(float value);
    void setIsRevBufIR2(float value);
    void setIsRevBufGR1(float value);
    void setIsRevBufIFFT(float value);


    void setMaxLenIR1(float value);
    void setMaxLenIR2(float value);
    void setMaxLenGR1(float value);

    void setStartIdxIR1(float value);
    void setStopIdxIR1(float value);
    void setStartIdxIR2(float value);
    void setStopIdxIR2(float value);
    void setStartIdxGR1(float value);
    void setStopIdxGR1(float value);
    void setDurRepeatGR1(float value);

    void setIsFindSilIFFT(float value);
    void setStartFindProcDur(float value);
    void setIsNoteOnGR1(float value);
	// Play buffers

    void playBufIR1(float* &out1, float* &out2, VstInt32 &i);
	void playBufIR2(float* &out1, float* &out2, VstInt32 &i);
    void playBufGR1(float* &out1, float* &out2, VstInt32 &i, 
        std::normal_distribution<float> &normal_dist,
        std::uniform_int_distribution<int> &uni_dist);
    void playConv(float* &out1, float* &out2, VstInt32 &i);

    // Rec buffers

    void recBufIR1(float* &in1, float* &in2, VstInt32 &i);
    void recBufIR2(float* &in1, float* &in2, VstInt32 &i);
    void recBufGR1(float* &in1, float* &in2, VstInt32 &i);
    
    void updateInfo();

	bool reinit_buf(VstInt32 index);
    bool convolute();
    bool convolve(); // Actual processing function
    void procRepCounter();

    // bool conv();
    // // Clears the buffers and refreshes the old_maxLength_IRs
    bool clear_buffer_IR1();
    bool clear_buffer_IR2();
    bool clear_buffer_GR1();

    void setStartGainAttack_GR1(float value);
    void setStartGainDecay_GR1(float value); 
    void setStartGainRelease_GR1(float value); 

    void setDurAttack_GR1(float value);
    void setDurDecay_GR1(float value);
    void setDurRelease_GR1(float value);
    void calcAttackFac_GR1();
    void calcDecayFac_GR1();
    void calcReleaseFac_GR1();

    void calcNormalDist_GR1();

    void shiftStartStopIdx_GR1(float shiftIdx);

    void setIsStdDev_GR1(float value);
    void setStdDev_GR1(float value);
    void setRandomProb_GR1(float value);

	public:

	/* initialization and terminati */
		Granuvolve (audioMasterCallback audioMaster);
		~Granuvolve ();

	/* processing */
		void processReplacing (float** inputs, float** outputs,
			VstInt32 sampleFrames);
        // VstInt32 processEvents(VstEvents* ev);

		void setParameter (VstInt32 index, float value);
		float getParameter (VstInt32 index);

		void getParameterLabel(VstInt32 index, char* label);
		void getParameterDisplay (VstInt32 index, char* text); 
		void getParameterName(VstInt32 index, char* text);

        
        

};