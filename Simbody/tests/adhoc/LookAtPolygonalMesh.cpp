/* -------------------------------------------------------------------------- *
 *             Simbody(tm) Adhoc test: Look at polygonal mesh                 *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

/*                      Simbody LookAtPolygonalMesh
Utility that exercises PolygonalMesh's load-from-file methods and then
displays the resulting mesh in the visualizer. It will also try to turn
the mesh into ContactGeometry and report any problems that occur. */

#include "Simbody.h"

#include <cassert>
#include <iostream>
using std::cout; using std::endl;

using namespace SimTK;

class ShowMesh : public DecorationGenerator {
public:
    ShowMesh() {} 

    // This is a shallow reference to the supplied mesh.
    void setMesh(const PolygonalMesh& mesh) {m_mesh=mesh;}

    void generateDecorations(const State&                state, 
                             Array_<DecorativeGeometry>& geometry) OVERRIDE_11
    {
        const Real TextScale = 0.1;
        DecorativeText info; info.setIsScreenText(true);
        info.setText("Faces/vertices: " + String(m_mesh.getNumFaces()) 
                     + "/" + String(m_mesh.getNumVertices()));
        geometry.push_back(info);
        if (!m_mesh.getNumFaces()) 
            return;

        DecorativeMesh dmesh(m_mesh);
        geometry.push_back(DecorativeMesh(dmesh)
                           .setOpacity(.8).setColor(Cyan));
        geometry.push_back(DecorativeMesh(dmesh)
                           .setRepresentation(DecorativeGeometry::DrawWireframe)
                           .setLineThickness(3)
                           .setColor(Black));
    }
private:
    PolygonalMesh m_mesh;
};

int main() {
  try {    
    // Create a system containing only Ground.   
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);

    system.setUseUniformBackground(true); // no ground plane in display
    system.realizeTopology();

    Visualizer viz(system);
    ShowMesh* sp = new ShowMesh();
    viz.addDecorationGenerator(sp);
    viz.report(system.getDefaultState()); // show default shape

    // The Visualizer caches meshes by their addresses so won't update 
    // properly if memory gets reused in a modified mesh. Since there's a 
    // human in the loop here we can just keep all the meshes around until
    // termination; that way there is no danger of reuse.
    Array_<PolygonalMesh> meshes;
    meshes.push_back(PolygonalMesh::createSphereMesh(1.,3));
    sp->setMesh(meshes.back());
    viz.report(system.getDefaultState());
    viz.zoomCameraToShowAllGeometry();


    std::string line, cwd = Pathname::getCurrentWorkingDirectory();
    std::cout << "Current working directory: " << cwd << std::endl;
    std::cout << "Change working directory (ENTER to keep): ";
    std::getline(std::cin, line);
    if (!line.empty()) cwd = line;

    while(true) {
        viz.report(system.getDefaultState());
        viz.zoomCameraToShowAllGeometry();

        meshes.push_back(); // make a new empty mesh
        PolygonalMesh& mesh = meshes.back();

        printf("mesh file name (or 'end'): ");
        std::getline(std::cin, line);
        if (line=="end")
            break;

        std::string dir,fn,ext;
        bool isAbsolutePath;
        Pathname::deconstructPathname(line,isAbsolutePath,dir,fn,ext);
        if (!isAbsolutePath) line = cwd + "/" + line;

        if (!Pathname::fileExists(line)) {
            if (!line.empty()) printf("'%s' doesn't exist\n", line.c_str());
            sp->setMesh(mesh);
            continue;
        }

        try {
            mesh.loadFile(line);
        } catch(const std::exception& e) {
            cout << "File loader error: " << e.what() << "\n";
            mesh.clear();
        }
        sp->setMesh(mesh);

        try {
            ContactGeometry::TriangleMesh tri(mesh);
            printf("*** Works as contact mesh! ***\n");
        } catch(const std::exception& e) {
            cout << "XXX Can't make contact mesh: " << e.what() << "\n";
        }
    }

  } catch (const std::exception& e) {
    cout << "EXCEPTION: " << e.what() << "\n";
  }
}
