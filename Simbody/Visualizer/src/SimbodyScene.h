#ifndef _SIMBODY_SCENE_H_
#define _SIMBODY_SCENE_H_

#include <vector>

#include "SimTKcommon.h"

namespace SimTK
{

class CustomMesh {
public:
    CustomMesh(std::vector<float>& vertices, std::vector<float>& normals, std::vector<unsigned short>& faces) 
    :   vertices(vertices), normals(normals), numVertices((int)(vertices.size()/3)), faces(faces) {
        // Create the list of edges.

		std::set<std::pair<unsigned short, unsigned short> > edgeSet;
        for (int i = 0; i < (int) faces.size(); i += 3) {
            unsigned short v1 = faces[i];
            unsigned short v2 = faces[i+1];
            unsigned short v3 = faces[i+2];
            edgeSet.insert(std::make_pair(std::min(v1, v2), std::max(v1, v2)));
            edgeSet.insert(std::make_pair(std::min(v2, v3), std::max(v2, v3)));
            edgeSet.insert(std::make_pair(std::min(v3, v1), std::max(v3, v1)));
        }
        for (std::set<std::pair<unsigned short, unsigned short> >::const_iterator iter = edgeSet.begin(); iter != edgeSet.end(); ++iter) {
            edges.push_back(iter->first);
            edges.push_back(iter->second);
        }
    }
    void getBoundingSphere(float& radius, SimTK::fVec3& center) {
        radius = this->radius;
        center = this->center;
    }

	const std::vector<float>& getVertices() const {
		return vertices;
	}
	
	const std::vector<float>& getNormals() const {
		return normals;
	}

	const std::vector<unsigned short>& getEdges() const {
		return edges;
	}
	const std::vector<unsigned short>& getFaces() const {
		return faces;
	}


private:
    int numVertices;
	std::vector<float> vertices, normals;
	std::vector<unsigned short> edges, faces;
	SimTK::fVec3 center;
    float radius;
};

}
#endif //_SCENE_H_
