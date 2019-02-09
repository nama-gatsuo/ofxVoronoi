#pragma once

#include "ofMain.h"
#include "ofxVoronoi.h"

class ofApp : public ofBaseApp{
    private:
        vector <glm::vec3> generateRandomPoints(int count, int seed, ofRectangle bounds);
        ofxVoronoi voronoi;
        vector<glm::vec3> points;
    
        bool isBorder(ofPoint _pt);
    
	public:
		void setup();
		void update();
		void draw();

		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y );
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);
};
