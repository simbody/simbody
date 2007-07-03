#ifndef SimTK_SIMBODY_VTK_DECORATIVE_GEOMETRY_H_
#define SimTK_SIMBODY_VTK_DECORATIVE_GEOMETRY_H_

/* Portions copyright (c) 2005-6 Stanford University and Jack Middleton.
 * Contributors:
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/DecorativeGeometry.h"
#include "Reporter3D.h"
#include "simbody/internal/AnalyticGeometry.h"

#include <cmath>
#include <vector>

namespace SimTK {


class VTKDecorativeGeometry : public Reporter3DGeom {
public:
    VTKDecorativeGeometry() 
      : myHandle(0) 
    { 
    }
    virtual ~VTKDecorativeGeometry() { }

    vtkPolyData* getVTKPolyData() {
        assert(vtkObjects.size());
        return vtkPolyDataAlgorithm::SafeDownCast(vtkObjects.back())->GetOutput();
    }


    // Caller must be sure to call VTK's Delete() methods on these objects
    // when done with them.
    void createVTKPolyData();

    // Combine a SimTK Transform with a scale factor and return the equivalent
    // vtkTransform that will translate, rotate, and scale the same way.
    vtkTransform* createVTKTransform(const Transform&, const Real&);

    // Take some polygon data as input, copy it and transform it, returning the
    // transformed polygon data.
    vtkPolyData* transformVTKPolyData(const Transform&, const Real&, vtkPolyData*);

    // returns a pointer to the vtkPolyData for this DecorativeGeometry object
    void* getReporterPolyData(){ return (void *)getVTKPolyData(); }
      

    void setMyHandle(const DecorativeGeometry& h) {myHandle = &h;}
    void clearMyHandle() {myHandle=0;}

protected:
    void rememberVTKObject(vtkObject* o) {
        vtkObjects.push_back(o);
    }

    DecorativeGeometry const *myHandle;     // DecorativeGeomety object for this VTKDecorativeGeometry

    void deleteVTKGeometry() { // Delete in reverse order of allocation
//      std::cout << "deleteVTKGeometry() \n";
        for (int i=(int)vtkObjects.size()-1; i >= 0; --i) {
            vtkObject* obj = vtkObjects[i];
//            std::cout << "ABOUT TO DELETE ";
//            obj->Print(std::cout);
            obj->Delete();
            vtkObjects[i]=0;
        }
        vtkObjects.resize(0);
    }


    // As we build the pipeline, we accumulate VTK objects which must
    // have their Delete() methods called in the desctructor.
    std::vector<vtkObject*> vtkObjects;
};

    ///////////////////////
    // VTKDecorativeLine //
    ///////////////////////

class VTKDecorativeLine : public VTKDecorativeGeometry{
public:
    // no default constructor
    VTKDecorativeLine(const DecorativeGeometry& geom, const Vec3& p1, const Vec3& p2)  {
       setMyHandle(  geom );
       createVTKPolyData(p1,p2);
    }
    ~VTKDecorativeLine() {
        deleteVTKGeometry();
        clearMyHandle();
    }

    void createVTKPolyData(const Vec3&, const Vec3& );

};

    /////////////////////////
    // VTKDecorativeCircle //
    /////////////////////////

class VTKDecorativeCircle : public VTKDecorativeGeometry {
public:
    // no default constructor
    VTKDecorativeCircle(const DecorativeGeometry& geom, Real r) {
       setMyHandle(  geom );
       createVTKPolyData(r);
    }
    ~VTKDecorativeCircle() {
        deleteVTKGeometry();
        clearMyHandle();
    }

    void createVTKPolyData(Real );


};

    /////////////////////////
    // VTKDecorativeSphere //
    /////////////////////////

class VTKDecorativeSphere : public VTKDecorativeGeometry {
    static const int DefaultResolution = 15;
public:
    // no default constructor
    VTKDecorativeSphere(const DecorativeGeometry& geom, Real r)  {
       setMyHandle( geom );
       createVTKPolyData(r);
    }
    ~VTKDecorativeSphere() {
        deleteVTKGeometry();
        clearMyHandle();
    }

    void createVTKPolyData(Real);

};

    ////////////////////////
    // VTKDecorativeBrick //
    ////////////////////////

class VTKDecorativeBrick : public VTKDecorativeGeometry {
public:
    // no default constructor
    VTKDecorativeBrick(const DecorativeGeometry& geom, const Vec3& halfHeights) {
       setMyHandle(  geom );
       createVTKPolyData( halfHeights );
    }
    ~VTKDecorativeBrick() {
        deleteVTKGeometry();
        clearMyHandle();
    }

    // virtuals
    void createVTKPolyData( const Vec3& );

};


class VTKDecorativeCylinder : public VTKDecorativeGeometry {
    static const int DefaultResolution = 10;
public:
    // no default constructor
    VTKDecorativeCylinder(const DecorativeGeometry& geom, Real r, Real halfHeight)  {
       setMyHandle( geom );
       createVTKPolyData(r, halfHeight);
    }
    ~VTKDecorativeCylinder() {
        deleteVTKGeometry();
        clearMyHandle();
    }


    // virtuals
    void createVTKPolyData(Real, Real);

};

class VTKDecorativeFrame : public VTKDecorativeGeometry{
public:
    // no default constructor
    VTKDecorativeFrame(const DecorativeGeometry& geom, Real halfLength)  {
       setMyHandle(  geom );
       createVTKPolyData(halfLength);
    }
    ~VTKDecorativeFrame() {
        deleteVTKGeometry();
        clearMyHandle();
    }


    void createVTKPolyData( Real );

};

} // namespace SimTK

#endif // SimTK_SIMBODY_VTK_DECORATIVE_GEOMETRY_REP_H_
