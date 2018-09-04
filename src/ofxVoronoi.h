#pragma once

// openFrameworks
#include "ofMain.h"

using namespace glm;

class ofxVoronoiCell {
public:
    vector<vec3> pts;
    vec3 pt;
};

class ofxVoronoi {
private:
    ofRectangle bounds;
    vector<vec3> points;
    vector<ofxVoronoiCell> cells;
    
public:
    ofxVoronoi();
    ~ofxVoronoi();
    
    void clear();
    void generate(bool ordered=true);
    void draw();
    
    bool isBorder(const vec3& _pt);
    
    void setBounds(ofRectangle _bounds);
    void setPoints(const vector<vec3>& _points);
    void addPoint(const vec3& _point);
    void addPoints(const vector<vec3>& _points);
    
    const ofRectangle& getBounds() const;
	const vector<vec3>& getPoints() const;
	const vector<ofxVoronoiCell>& getCells() const;
	const ofxVoronoiCell& getCell(const vec3& _point, bool approximate=false);
    
    //borg
    void relax();
};