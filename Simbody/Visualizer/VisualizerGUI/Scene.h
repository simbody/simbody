#ifndef _SCENE_H_
#define _SCENE_H_

#include <vector>
#include "SimTKcommon.h"

using namespace std;
using namespace SimTK;

class RenderedLine {
public:
    RenderedLine(const fVec3& color, float thickness)
    :   color(color), thickness(thickness) {}

    vector<float>& getLines() {
        return lines;
    }
    const fVec3& getColor() const {
        return color;
    }
    float getThickness() const {
        return thickness;
    }
    void computeBoundingSphere(float& radius, fVec3& center) const {
//        computeBoundingSphereForVertices(lines, radius, center);
    }
private:
    fVec3 color;
    float thickness;
    vector<float> lines;
};

class RenderedText {
public:
    RenderedText(const fVec3& position, const fVec3& scale, const fVec3& color, 
                 const string& text, bool faceCamera = true) 
    :   position(position), scale(scale/119), text(text),
        faceCamera(faceCamera) {
        this->color[0] = color[0];
        this->color[1] = color[1];
        this->color[2] = color[2];
    }
    void draw() {
/*
        glPushMatrix();
        glTranslated(position[0], position[1], position[2]);
        fVec4 rot = X_GC.R().convertRotationToAngleAxis();
        if (faceCamera)
            glRotated(rot[0]*SimTK_RADIAN_TO_DEGREE, rot[1], rot[2], rot[3]);
        glScaled(scale[0], scale[1], scale[2]);
        glColor3fv(color);
        for (int i = 0; i < (int) text.size(); i++)
            glutStrokeCharacter(GLUT_STROKE_ROMAN, text[i]);
        glPopMatrix();
*/
    }
    void computeBoundingSphere(float& radius, fVec3& center) const {
//        center = position;
//        radius = glutStrokeLength(GLUT_STROKE_ROMAN, 
//                                  (unsigned char*)text.c_str())*scale[0];
    }
private:
    fVec3 position;
    fVec3 scale;
    float color[3];
    string text;
    bool faceCamera;
};

class ScreenText {
public:
    ScreenText(const string& txt) 
    :   text(txt)  {
        //this->color[0] = color[0];
        //this->color[1] = color[1];
        //this->color[2] = color[2];
    }

    const char* getText() {return text.c_str();}

private:
    //float color[3];
    string text;
};

class RenderedMesh {
public:
    RenderedMesh(const fTransform& transform, const fVec3& scale, const fVec4& color, short representation, unsigned short meshIndex, unsigned short resolution) :
            transform(transform), scale(scale), representation(representation), meshIndex(meshIndex), resolution(resolution) {
        this->color[0] = color[0];
        this->color[1] = color[1];
        this->color[2] = color[2];
        this->color[3] = color[3];
    }
    void draw(bool setColor = true) {
		/*
        glPushMatrix();
        glTranslated(transform.p()[0], transform.p()[1], transform.p()[2]);
        fVec4 rot = transform.R().convertRotationToAngleAxis();
        glRotated(rot[0]*SimTK_RADIAN_TO_DEGREE, rot[1], rot[2], rot[3]);
        glScaled(scale[0], scale[1], scale[2]);
        if (setColor) {
            if (representation == DecorativeGeometry::DrawSurface)
                glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color);
            else
                glColor3fv(color);
        }
        meshes[meshIndex][resolution]->draw(representation);
        glPopMatrix();
		*/
    }
    const fTransform& getTransform() const {
        return transform;
    }
    void computeBoundingSphere(float& radius, fVec3& center) const {
//        meshes[meshIndex][resolution]->getBoundingSphere(radius, center);
        center += transform.p();
        radius *= max(abs(scale[0]), max(abs(scale[1]), abs(scale[2])));
    }
private:
    fTransform transform;
    fVec3 scale;
    float color[4];
    short representation;
    unsigned short meshIndex, resolution;
};


class Scene {

public:
    Scene() : simTime(0), sceneHasBeenDrawn(false) {}

    float simTime; // simulated time associated with this frame

    vector<RenderedMesh> drawnMeshes;
    vector<RenderedMesh> solidMeshes;
    vector<RenderedMesh> transparentMeshes;
    vector<RenderedLine> lines;
    vector<RenderedText> sceneText;
    vector<ScreenText>   screenText;

    bool sceneHasBeenDrawn;
};

#endif //_SCENE_H_
