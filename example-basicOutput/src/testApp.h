#pragma once

#include "ofMain.h"
#include "ofxAudioAnalyzer.h"

class testApp : public ofBaseApp{

	public:
		void setup();
		void update();
		void draw();
        void exit();

		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y );
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);

        void audioOut(float * input, int bufferSize, int nChannels);
    
        ofSoundStream soundStream;
        ofxAudioAnalyzer audioAnalyzer1, audioAnalyzer2;
    
        float *buffer_1;
        float *buffer_2;
    
        float *spectrum1;
        int spectrumSize;
};
