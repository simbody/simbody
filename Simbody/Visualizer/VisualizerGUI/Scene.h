#ifndef _SCENE_H_
#define _SCENE_H_

#include <vector>
#include "SimTKcommon.h"

using namespace std;
using namespace SimTK;

namespace Visualizer {

class Mesh {
public:
    Mesh(vector<float>& vertices, vector<float>& normals, vector<unsigned short>& faces) 
    :   vertices(vertices), normals(normals), numVertices((int)(vertices.size()/3)), faces(faces) {
		/*
        // Build OpenGL buffers.

        GLuint buffers[2];
        glGenBuffers(2, buffers);
        vertBuffer = buffers[0];
        normBuffer = buffers[1];
        glBindBuffer(GL_ARRAY_BUFFER, vertBuffer);
        glBufferData(GL_ARRAY_BUFFER, vertices.size()*sizeof(float), &vertices[0], GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, normBuffer);
        glBufferData(GL_ARRAY_BUFFER, normals.size()*sizeof(float), &normals[0], GL_STATIC_DRAW);

		*/
        // Create the list of edges.

        set<pair<unsigned short, unsigned short> > edgeSet;
        for (int i = 0; i < (int) faces.size(); i += 3) {
            unsigned short v1 = faces[i];
            unsigned short v2 = faces[i+1];
            unsigned short v3 = faces[i+2];
            edgeSet.insert(make_pair(min(v1, v2), max(v1, v2)));
            edgeSet.insert(make_pair(min(v2, v3), max(v2, v3)));
            edgeSet.insert(make_pair(min(v3, v1), max(v3, v1)));
        }
        for (set<pair<unsigned short, unsigned short> >::const_iterator iter = edgeSet.begin(); iter != edgeSet.end(); ++iter) {
            edges.push_back(iter->first);
            edges.push_back(iter->second);
        }

        // Compute the center and radius.

        //computeBoundingSphereForVertices(vertices, radius, center);
    }
    void draw(short representation) const {
/*
        glBindBuffer(GL_ARRAY_BUFFER, vertBuffer);
        glVertexPointer(3, GL_FLOAT, 0, 0);
        glBindBuffer(GL_ARRAY_BUFFER, normBuffer);
        glNormalPointer(GL_FLOAT, 0, 0);
        if (representation == DecorativeGeometry::DrawSurface)
            glDrawElements(GL_TRIANGLES, (GLsizei)faces.size(), GL_UNSIGNED_SHORT, &faces[0]);
        else if (representation == DecorativeGeometry::DrawPoints)
            glDrawArrays(GL_POINTS, 0, numVertices);
        else if (representation == DecorativeGeometry::DrawWireframe)
            glDrawElements(GL_LINES, (GLsizei)edges.size(), GL_UNSIGNED_SHORT, &edges[0]);
*/
    }
    void getBoundingSphere(float& radius, fVec3& center) {
        radius = this->radius;
        center = this->center;
    }

	vector<float>& getVertices() {
		return vertices;
	}
	
	vector<float>& getNormals() {
		return normals;
	}

	vector<unsigned short>& getEdges() {
		return edges;
	}
	vector<unsigned short>& getFaces() {
		return faces;
	}


private:
    int numVertices;
//    GLuint vertBuffer, normBuffer;
//    vector<GLushort> edges, faces;
    vector<float> vertices, normals;
    vector<unsigned short> edges, faces;
    fVec3 center;
    float radius;
};

//	PendingCommand
// Some commands received by the listener must be executed on the rendering
// thread. They are saved in concrete objects derived from this abstract
// class.
class PendingCommand {
public:
    virtual ~PendingCommand() {}
    virtual void execute(vector<vector<Mesh *> >& meshes) = 0;
};

//	PendingMesh
class PendingMesh : public PendingCommand {
public:
    PendingMesh() {
//        index = nextMeshIndex++;
    }
    void execute(vector<vector<Mesh*> >& meshes) {
        if ((int) meshes.size() <= index)
            meshes.resize(index+1);
        meshes[index].push_back(new Mesh(vertices, normals, faces));
    }
    vector<float> vertices;
    vector<float> normals;
//    vector<GLushort> faces;
    vector<unsigned short> faces;
    int index;
};

//	PendingStandardMesh
class PendingStandardMesh : public PendingCommand {
public:
    PendingStandardMesh(unsigned short meshIndex, unsigned short mResolution) : meshIndex(meshIndex), mResolution(mResolution) {
    }
    void execute(vector<vector<Mesh*> >& meshes) {
        if ((int) meshes[meshIndex].size() <= mResolution)
            meshes[meshIndex].resize(mResolution+1, NULL);
        if (meshes[meshIndex][mResolution] == NULL) {
            switch (meshIndex) {
                case MeshBox:
                    meshes[meshIndex][mResolution] = makeBox();
                    break;
                case MeshEllipsoid:
                    meshes[meshIndex][mResolution] = makeSphere(mResolution);
                    break;
                case MeshCylinder:
                    meshes[meshIndex][mResolution] = makeCylinder(mResolution);
                    break;
                case MeshCircle:
                    meshes[meshIndex][mResolution] = makeCircle(mResolution);
                    break;
            }
        }
    }

	void addVec(vector<float>& data, float x, float y, float z) {
		data.push_back(x);
		data.push_back(y);
		data.push_back(z);
	}

	void addVec(vector<unsigned short>& data, int x, int y, int z) {
		data.push_back((unsigned short) x);
		data.push_back((unsigned short) y);
		data.push_back((unsigned short) z);
	}
	
	Mesh * makeBox()
	{
		const float halfx = 1;
		const float halfy = 1;
		const float halfz = 1;
		vector<float> vertices;
		vector<float> normals;
		vector<unsigned short> faces;

		// lower x face
		addVec(vertices, -halfx, -halfy, -halfz);
		addVec(vertices, -halfx, -halfy, halfz);
		addVec(vertices, -halfx, halfy, halfz);
		addVec(vertices, -halfx, halfy, -halfz);
		addVec(normals, -1, 0, 0);
		addVec(normals, -1, 0, 0);
		addVec(normals, -1, 0, 0);
		addVec(normals, -1, 0, 0);
		addVec(faces, 0, 1, 2);
		addVec(faces, 2, 3, 0);
		// upper x face
		addVec(vertices, halfx, halfy, halfz);
		addVec(vertices, halfx, -halfy, halfz);
		addVec(vertices, halfx, -halfy, -halfz);
		addVec(vertices, halfx, halfy, -halfz);
		addVec(normals, 1, 0, 0);
		addVec(normals, 1, 0, 0);
		addVec(normals, 1, 0, 0);
		addVec(normals, 1, 0, 0);
		addVec(faces, 4, 5, 6);
		addVec(faces, 6, 7, 4);
		// lower y face
		addVec(vertices, -halfx, -halfy, -halfz);
		addVec(vertices, halfx, -halfy, -halfz);
		addVec(vertices, halfx, -halfy, halfz);
		addVec(vertices, -halfx, -halfy, halfz);
		addVec(normals, 0, -1, 0);
		addVec(normals, 0, -1, 0);
		addVec(normals, 0, -1, 0);
		addVec(normals, 0, -1, 0);
		addVec(faces, 8, 9, 10);
		addVec(faces, 10, 11, 8);
		// upper y face
		addVec(vertices, halfx, halfy, halfz);
		addVec(vertices, halfx, halfy, -halfz);
		addVec(vertices, -halfx, halfy, -halfz);
		addVec(vertices, -halfx, halfy, halfz);
		addVec(normals, 0, 1, 0);
		addVec(normals, 0, 1, 0);
		addVec(normals, 0, 1, 0);
		addVec(normals, 0, 1, 0);
		addVec(faces, 12, 13, 14);
		addVec(faces, 14, 15, 12);
		// lower z face
		addVec(vertices, -halfx, -halfy, -halfz);
		addVec(vertices, -halfx, halfy, -halfz);
		addVec(vertices, halfx, halfy, -halfz);
		addVec(vertices, halfx, -halfy, -halfz);
		addVec(normals, 0, 0, -1);
		addVec(normals, 0, 0, -1);
		addVec(normals, 0, 0, -1);
		addVec(normals, 0, 0, -1);
		addVec(faces, 16, 17, 18);
		addVec(faces, 18, 19, 16);
		// upper z face
		addVec(vertices, halfx, halfy, halfz);
		addVec(vertices, -halfx, halfy, halfz);
		addVec(vertices, -halfx, -halfy, halfz);
		addVec(vertices, halfx, -halfy, halfz);
		addVec(normals, 0, 0, 1);
		addVec(normals, 0, 0, 1);
		addVec(normals, 0, 0, 1);
		addVec(normals, 0, 0, 1);
		addVec(faces, 20, 21, 22);
		addVec(faces, 22, 23, 20);
		return new Mesh(vertices, normals, faces);
	}

	Mesh * makeSphere(unsigned short resolution)
	{
	    const int numLatitude = 4*resolution;
		const int numLongitude = 6*resolution;
		const float radius = 1.0f;
		vector<float> vertices;
		vector<float> normals;
		vector<unsigned short> faces;
		addVec(vertices, 0, radius, 0);
		addVec(normals, 0, 1, 0);
		for (int i = 0; i < numLatitude; i++) {
			float phi = (float) (((i+1)*SimTK_PI)/(numLatitude+1));
			float sphi = sin(phi);
			float cphi = cos(phi);
			float y = radius*cphi;
			float r = radius*sphi;
			for (int j = 0; j < numLongitude; j++) {
				float theta = (float) ((j*2*SimTK_PI)/numLongitude);
				float stheta = sin(theta);
				float ctheta = cos(theta);
				addVec(vertices, r*ctheta, y, r*stheta);
				addVec(normals, sphi*ctheta, cphi, sphi*stheta);
			}
		}
		addVec(vertices, 0, -radius, 0);
		addVec(normals, 0, -1, 0);
		for (int i = 1; i < numLongitude; i++)
			addVec(faces, 0, i+1, i);
		addVec(faces, 0, 1, numLongitude);
		for (int i = 1; i < numLatitude; i++) {
			int base = (i-1)*numLongitude+1;
			for (int j = 0; j < numLongitude; j++) {
				int v1 = base+j;
				int v2 = (j == numLongitude-1 ? base : v1+1);
				int v3 = v1+numLongitude;
				int v4 = v2+numLongitude;
				addVec(faces, v1, v4, v3);
				addVec(faces, v1, v2, v4);
			}
		}
		int first = (numLatitude-1)*numLongitude+1;
		int last = numLatitude*numLongitude+1;
		for (int i = first; i < last-1; i++)
			addVec(faces, i, i+1, last);
		addVec(faces, last-1, first, last);
		return new Mesh(vertices, normals, faces);


	}

	Mesh* makeCylinder(unsigned short resolution) {
		const int numSides = 6*resolution;
		const float halfHeight = 1;
		const float radius = 1;
		vector<float> vertices;
		vector<float> normals;
		vector<unsigned short> faces;

		// Create the top face.

		addVec(vertices, 0, halfHeight, 0);
		addVec(normals, 0, 1.0, 0);
		for (int i = 0; i < numSides; i++) {
			float theta = (float) ((i*2*SimTK_PI)/numSides);
			float stheta = sin(theta);
			float ctheta = cos(theta);
			addVec(vertices, radius*ctheta, halfHeight, radius*stheta);
			addVec(normals, 0, 1.0, 0);
		}
		for (int i = 1; i < numSides; i++)
			addVec(faces, i, 0, i+1);
		addVec(faces, numSides, 0, 1);

		// Create the bottom face.

		int bottomStart = numSides+1;
		addVec(vertices, 0, -halfHeight, 0);
		addVec(normals, 0, -1.0, 0);
		for (int i = 0; i < numSides; i++) {
			float theta = (float) ((i*2*SimTK_PI)/numSides);
			float stheta = sin(theta);
			float ctheta = cos(theta);
			addVec(vertices, radius*ctheta, -halfHeight, radius*stheta);
			addVec(normals, 0, -1.0, 0);
		}
		for (int i = 1; i < numSides; i++)
			addVec(faces, bottomStart+i, bottomStart+i+1, bottomStart);
		addVec(faces, bottomStart+numSides, bottomStart+1, bottomStart);

		// Create the sides.

		for (int i = 0; i < numSides; i++) {
			float theta = (float) ((i*2*SimTK_PI)/numSides);
			float stheta = sin(theta);
			float ctheta = cos(theta);
			float x = radius*ctheta;
			float z = radius*stheta;
			addVec(vertices, x, halfHeight, z);
			addVec(normals, ctheta, 0, stheta);
			addVec(vertices, x, -halfHeight, z);
			addVec(normals, ctheta, 0, stheta);
		}
		int sideStart = 2*numSides+2;
		for (int i = 0; i < numSides-1; i++) {
			int base = sideStart+2*i;
			addVec(faces, base, base+2, base+1);
			addVec(faces, base+1, base+2, base+3);
		}
		addVec(faces, sideStart+2*numSides-2, sideStart, sideStart+2*numSides-1);
		addVec(faces, sideStart+2*numSides-1, sideStart, sideStart+1);
		return new Mesh(vertices, normals, faces);
	}

	Mesh* makeCircle(unsigned short resolution) 
	{
		const int numSides = 6*resolution;
		const float radius = 1;
		vector<float> vertices;
		vector<float> normals;
		vector<unsigned short> faces;

		// Create the front face.

		addVec(vertices, 0, 0, 0);
		addVec(normals, 0, 0, -1.0);
		for (int i = 0; i < numSides; i++) {
			float theta = (float) ((i*2*SimTK_PI)/numSides);
			float stheta = sin(theta);
			float ctheta = cos(theta);
			addVec(vertices, radius*ctheta, radius*stheta, 0);
			addVec(normals, 0, 0, -1.0);
		}
		for (int i = 1; i < numSides; i++)
			addVec(faces, i, 0, i+1);
		addVec(faces, numSides, 0, 1);

		// Create the back face.

		addVec(vertices, 0, 0, 0);
		addVec(normals, 0, 0, 1.0);
		for (int i = 0; i < numSides; i++) {
			float theta = (float) ((i*2*SimTK_PI)/numSides);
			float stheta = sin(theta);
			float ctheta = cos(theta);
			addVec(vertices, radius*ctheta, radius*stheta, 0);
			addVec(normals, 0, 0, 1.0);
		}
		int backStart = numSides+1;
		for (int i = 1; i < numSides; i++)
			addVec(faces, backStart, backStart+i, backStart+i+1);
		addVec(faces, backStart, backStart+numSides, backStart+1);
		return new Mesh(vertices, normals, faces);
	}


    unsigned short meshIndex, mResolution;
};

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

	const std::string& getText() const { return text; }

	const fVec3& getPosition() const { return position; }

	const fVec3& getScale() const { return scale; }

	const fVec3& getColor() const { return color; }

private:
    fVec3 position;
    fVec3 scale;
    fVec3 color;
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
//        center += transform.p();
//        radius *= max(abs(scale[0]), max(abs(scale[1]), abs(scale[2])));
    }

	short getRepresentation() { return representation; }

	unsigned short getMeshIndex() { return meshIndex; }
	unsigned short getResolution() { return resolution; }

	const fVec4& getColor() { return color; }

	const fVec3& getScale() { return scale; }

private:
    fTransform transform;
    fVec3 scale;
    fVec4 color;
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

}

#endif //_SCENE_H_
