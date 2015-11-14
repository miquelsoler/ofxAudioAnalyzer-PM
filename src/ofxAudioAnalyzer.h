#pragma once

#include "ofMain.h"

#include <iostream>
#include "essentia/algorithmfactory.h"
#include "essentia/essentiamath.h" // for the isSilent function
#include "essentia/pool.h"

using namespace std;
using namespace essentia;
using namespace standard;


class ofxAudioAnalyzer
{
public:

    void setup(int bufferSize, int sampleRate,
            bool doMelbands, int numMelBands,
            bool doSilence, int silenceThreshold, unsigned int silenceQueueLength);
    void exit();
    void analyze(float * iBuffer, int bufferSize);

    void resetOnsets();

    // Getters/Setters

    float getRms(){return rms_f;}
    float getEnergy(){return energy_f;}
    float getPower(){return power_f;}

    float getPitchFreq(){return FilteredPitch_f;}
    float getPitchConf(){return YinConfidence_f;}
    float getSalience(){return salience_f;}

    float getTuningFreq(){return tFreq_f;}
    float getTuningCents(){return tCents_f;}

    float getInharmonicity(){return inharm_f;}
    float getHfc(){return hfc_f;}
    float getCentroid(){return centroid_f;}
    float getSpectralComplex(){return spectralComplex_f;}

    float getOnsetHfc(){return onsetHfc_f;}
    float getOnsetComplex(){return onsetComplex_f;}
    float getOnsetFlux(){return onsetFlux_f;}

    bool getIsOnset(){return isOnset;}

    // TODO: Silence
    bool getIsSilent() { return isSilent; }

    float* getSpectrum(){return &spectrum_f[0];}
    float* getMelBands(){return &melBands_f[0];}
    float* getDct(){return &dct_f[0];}
    float* getHpcp(){return &hpcp_f[0];}

    void setOnsetTreshold(float val){onsetSilenceTreshold =val;}
    void setOnsetAlpha(float val){alpha=val;}

    // TODO: Silence
    int getSilenceStartFrame() { return silenceStartFrame_i; }
    int getSilenceStopFrame() { return silenceStopFrame_i; }

    // For storing results casted to Floats ---------------------

    bool addHfc, addComplex, addFlux;
    bool doRms, doEnergy, doPower, doPitch,
    doTuning, doOnsets, doMelbands, doMfcc,
    doHpcp, doHfc, doCentroid, doSpcCmplx, doInharmon;
    // TODO: Silence
    bool doStartStopSilence;

 private:

    // Getter return values
    float rms_f, energy_f, power_f;
    float YinFrequency_f, YinConfidence_f;
    float FilteredPitch_f;
    float salience_f;
    float tFreq_f, tCents_f;
    float inharm_f;
    float hfc_f;
    float centroid_f;
    float spectralComplex_f;
    float onsetHfc_f;
    float onsetComplex_f;
    float onsetFlux_f;

    // TODO: Silence
    int silenceStartFrame_i;
    int silenceStopFrame_i;


    bool isOnset;
    // TODO: Silence
    bool isSilent;

    vector<float> spectrum_f;
    vector<float> melBands_f;
    vector<float> dct_f;
    vector<float> hpcp_f;
//    vector<float> FilteredPitch_f;

    bool onsetEvaluation (Real iDetectHfc, Real iDetectComplex, Real iDetectFlux);

    // TODO: Silence
    bool silenceEvaluation();
    queue<unsigned int> silenceQueue;
    unsigned int silenceQueueLength;

    Real onsetSilenceTreshold, alpha;

    int framesize;
    int hopsize;
    int sr;
    int zeropadding;
    Real framerate;
    int numMelBands;

    //Algorithm pointers-----------------

    Algorithm* window;
    Algorithm* spectrum;
    Algorithm* pitchDetect;
    Algorithm* pitchFilter;
    Algorithm* rms;
    Algorithm* energy;
    Algorithm* power;
    Algorithm* peaks;
    Algorithm* tuningFreq;
    Algorithm* inharmonicity;
    Algorithm* melBands;
    Algorithm* dct;
    Algorithm* harmonicPeaks;
    Algorithm* hpcp;
    Algorithm* hfc;
    Algorithm* pitchSalience;
    Algorithm* centroid;
    Algorithm* spectralComplex;
    Algorithm* dcremoval;
    Algorithm* fft;
    Algorithm* cartesian2polar;
    Algorithm* onsetHfc;
    Algorithm* onsetComplex;
    Algorithm* onsetFlux;
    // TODO: Silence
    Algorithm* startStopSilence;


    //For storing algorithms results-----------------

    Real thisPitch, thisConf;
    Real rmsValue;
    Real powerValue;
    Real energyValue;
    Real tunFreqValue;
    Real tunCentsValue;
    Real inharmValue;
    Real hfcValue;
    Real salienceValue;
    Real centroidValue;
    Real spectralComplexValue;
    Real onsetHfcValue;
    Real onsetComplexValue;
    Real onsetFluxValue;

    vector<Real> audio;
    vector<Real> frame;
    vector<Real> windowedframe;
    vector<Real> spec;
    vector<Real> audioBuffer;
    vector<Real> audioBuffer_dc;
    vector<Real> peakFreqValues;
    vector<Real> peakMagValues;
    vector<Real> melValues;
    vector<Real> dctValues;
    vector<Real> logMelBands;
    vector<Real> harmFreqValues;
    vector<Real> harmMagValues;
    vector<Real> hpcpValues;
    vector<complex<Real> > fftValues;
    vector<Real> c2pMagValues;
    vector<Real> c2pPhaseValues;
    
    // Filter
    vector<Real> pitchBuffer;
    vector<Real> confidenceBuffer;
    vector<Real> thisPitchFiltered;
    int bufferFillIdx;
    int pitchBufferWidth;

    // TODO: Silence
    int silenceStartFrame;
    int silenceStopFrame;


    // Internal

    // Onset Evaluation

    int detecBufferSize;
    vector<vector<Real> > detections;
    vector<Real> detection_sum;
    Real hfc_max, complex_max, flux_max;

    // TODO: Silence
    // Silence
    int silenceLastStopFrame;
    bool silenceEvaluated;
};


