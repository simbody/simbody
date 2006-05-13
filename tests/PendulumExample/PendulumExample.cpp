/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
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

/**@file
 * The simple 2d pendulum example from the user's manual.
 */

#include "Simbody.h"
#include "simbody/internal/NumericalMethods.h"

#include "windows.h"

#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>

using namespace std;
using namespace SimTK;

static const Real Pi = std::acos(-1.), RadiansPerDegree = Pi/180;
static const int  GroundBodyNum = 0; // ground is always body 0

static const Real m = 5;   // kg
static const Real g = 9.8; // meters/s^2; apply in –y direction
static const Real d = 0.5; // meters
static const Real initialTheta   = 30;             // degrees
static const Real expectedPeriod = 2*Pi*sqrt(d/g); // s

class MySimbodyPendulum : public MechanicalDAESystem {
public:
    MySimbodyPendulum() 
    {
        pendBodyNum =
            pend.addRigidBody(
                MassProperties(m,        // body mass, center of mass, inertia
                               Vec3(0,0,0), 
                               InertiaMat(Vec3(0,-d/2,0), m)+InertiaMat(1e-3,1e-3,1e-3)),
                Transform(Vec3(0,d/2,0)),// jt frame on body (aligned w/body frame)
                GroundBodyNum,           // parent body
                Transform(),             // jt frame on parent (origin in this case)              
                JointSpecification(JointSpecification::Ball, false)); // joint type; pin always aligns z axes

        int pendBodyNum2 =
            pend.addRigidBody(
                MassProperties(m,        // body mass, center of mass, inertia
                               Vec3(0,0,0), 
                               InertiaMat(Vec3(0,-d/2,0), m)+InertiaMat(1e-3,1e-3,1e-3)),
                Transform(Vec3(0,d/2,0)),// jt frame on body (aligned w/body frame)
                pendBodyNum,           // parent body
                Transform(Vec3(0,-d/2,0)),             // jt frame on parent (bottom)             
                JointSpecification(JointSpecification::Ball, false)); // joint type; pin always aligns z axes


        int pendBodyNum3 =
            pend.addRigidBody(
                MassProperties(m,        // body mass, center of mass, inertia
                               Vec3(0,0,0), 
                               InertiaMat(Vec3(0,-d/2,0), m)+InertiaMat(1e-3,1e-3,1e-3)),
                Transform(Vec3(0,d/2,0)),// jt frame on body (aligned w/body frame)
                pendBodyNum2,           // parent body
                Transform(Vec3(0,-d/2,0)),             // jt frame on parent (bottom)             
                JointSpecification(JointSpecification::Ball, false)); // joint type; pin always aligns z axes
       
    int theConstraint =
        pend.addConstantDistanceConstraint(0, Vec3(0.5,-0.2,0.1),
                                           pendBodyNum3, Vec3(0,-d/2,0),
                                           1.);




        pend.realize(s, Stage::Built);
       // pend.setUseEulerAngles(s, true);
        pend.realize(s, Stage::Modeled);
        nq = s.getQ().size();
        nu = s.getU().size();
        nz = 0; // these would be for controllers and other forces with states

        printf("System size nq=%d nu=%d nz=%d (ny=%d)\n", nq, nu, nz, size());
        y.resize(nq+nu+nz);       y.setToNaN();
        ydot.resize(nq+nu+nz);    ydot.setToNaN();
        weights.resize(nq+nu+nz); weights.setToNaN();
        solutionAccuracy = constraintAccuracy = 1e-3; // default is 0.1%

        // Set state defaults
        s.updTime() = 0.;
        setPendulumAngle(initialTheta);
        copyStateToY();

        calcWeights(solutionAccuracy);
        ydotReady = false;
    }

    const SimbodySubsystem& getSimbodySubsystem() const {return pend;}
    const State& getState() const {return s;}

    void setTime(Real t) {s.updTime() = t;}

    Real getPendulumAngle() const {
        return pend.getJointQ(s,pendBodyNum,0)/RadiansPerDegree;
    }

    void setPendulumAngle(Real angleInDegrees) {
        for (int i=0; i<3; ++i)
            pend.setJointQ(s,pendBodyNum,i,angleInDegrees*RadiansPerDegree);
        copyStateToY();
    }

    // Supply required virtual methods.

    int size() const {return nq+nu+nz;}

    void setAccuracy(const Real& solution, const Real& constraint) {
        solutionAccuracy = solution; constraintAccuracy = constraint;
        calcWeights(solutionAccuracy);
    }

    void setState(const Real& t, const Vector& newY) {
        s.updTime() = t;
        y = newY;   // TODO: shouldn't need this temporary copy
        copyYtoState();
        calcWeights(solutionAccuracy);
        ydotReady = false;
    }

    bool realize() const {
        pend.realize(s, Stage::Moving); // kinematics

        // calculate and apply forces
        pend.clearAppliedForces(s);
        //pend.applyGravity(s, Vec3(0., -g, 0.));
        pend.applyGravity(s, Vec3(0, -g, 0));

        // calculate Simbody derivatives
        pend.realize(s, Stage::Reacting);

        // calculate other derivatives if needed
        //    none here

        // group derivatives together (TODO: should handle this without copying)
        copyStateToYDot();

        ydotReady = true;
        return true;
    }

    bool realizeAndProject(bool& anyChange, bool force) {
        pend.realize(s, Stage::Configured);
        pend.enforceConfigurationConstraints(s);
        pend.realize(s, Stage::Moving);
        pend.enforceMotionConstraints(s);
        return realize();
    }

    Real          getTimescale() const {return 0.1;} // a tenth of a second
    const Real&   getT()         const {return s.getTime();}
    const Vector& getY()         const {return y;}
    const Vector& getWeights()   const {return weights;}
    const Vector& getYDot()      const {assert(ydotReady); return ydot;}

private:
    void calcWeights(Real acc) {
        // Assume units around 1 for everything, and relative tolerance OK
        for (int i=0; i<size(); ++i)
            weights[i] = 1./(acc*fabs(y[i]) + acc);
    }

    //TODO: these shouldn't be needed
    void copyStateToY() {
        y(0,nq)     = s.getQ();
        y(nq,nu)    = s.getU();
        //y(nq+nu,nz) = s.getZ();
    }

    void copyStateToYDot() const {
        ydot(0,nq)     = s.getQDot();
        ydot(nq,nu)    = s.getUDot();
        //ydot(nq+nu,nz) = s.getZDot();
    }

    void copyYtoState() {
        s.updQ() = y(0,nq);
        s.updU() = y(nq,nu);
        //s.updZ() = y(nq+nu,nz);
    }

private:
    SimbodySubsystem pend;

    Vector y, weights;
    int nq, nu, nz;
    int pendBodyNum;
    Real solutionAccuracy, constraintAccuracy;

    mutable State s;         // TODO: state shouldn't have to be mutable; needed
                             //   because we are (unnecessarily) using it to hold
                             //   forces.

    // These are "cache" items and must be mutable.
    mutable Vector ydot;
    mutable bool   ydotReady; // sanity check
};



#include "vtkCommand.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkAssembly.h"
#include "vtkCamera.h"
#include "vtkActor.h"
#include "vtkFollower.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"

#include "vtkLineSource.h"
#include "vtkTextSource.h"
#include "vtkVectorText.h"
#include "vtkConeSource.h"
#include "vtkSphereSource.h"
#include "vtkCubeSource.h"
#include "vtkPointSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkCaptionActor2D.h"
#include "vtkAppendPolyData.h"
#include "vtkInteractorObserver.h"
#include "vtkInteractorStyleTrackballCamera.h"


// Callback for interaction
class vtkMyCallback2 : public vtkCommand {
public:
    static vtkMyCallback2* New() {return new vtkMyCallback2;}

    /*virtual*/ void Execute(vtkObject *caller, unsigned long eventID, void*) {
        vtkRenderWindowInteractor* iren =
            reinterpret_cast<vtkRenderWindowInteractor*>(caller);
        printf("EVENT: %s\n", vtkCommand::GetStringFromEventId(eventID));
        switch(eventID) {
            case vtkCommand::LeftButtonPressEvent: printf("LEFT\n"); break;
            default: printf("???\n");
        }

    }
};

typedef std::pair<vtkProp3D*, Transform> BodyActor;
typedef std::vector<BodyActor>           ActorList;

class VTKReporter {
public:
    class Decoration {
    public:
        Decoration() : data(0), actor(0) { }
        ~Decoration() {clear();}
        Decoration(const Decoration& d) : data(0), actor(0) {
            copyFrom(d);
        }

        Decoration& operator=(const Decoration& d) {
            if (&d != this)
                copyFrom(d);
            return *this;
        }

        vtkActor* getActor() {return actor;} // not a copy; use ShallowCopy if you need one
        vtkPolyData* getPolyData() {return data;}

        void setScale  (Real s)        {assert(actor);actor->SetScale(s);}
        void setColor  (const Vec3& c) {assert(actor);actor->GetProperty()->SetColor(c[0],c[1],c[2]);}
        void setOpacity(Real o)        {assert(actor);actor->GetProperty()->SetOpacity(o);}

    protected:
        void setPolyData(vtkPolyData* output) {
            data = vtkPolyData::New();
            data->DeepCopy(output);
            createActor();
            setDefaultColor();
        }

    private:
        vtkPolyData* data;
        vtkActor*    actor;

        void clear() {
            if (actor) actor->Delete(); actor=0; 
            if (data)  data->Delete();  data=0;
        }

        void copyFrom(const Decoration& d) {
            clear();
            if (d.data) {
                setPolyData(d.data);
                actor->SetProperty(d.actor->GetProperty());
                actor->SetScale(d.actor->GetScale());
            }
        }

        void setDefaultColor() {
            setColor(Vec3(0,0,1)); // blue
        }

        void createActor() {
            if (actor) actor->Delete();
            actor = vtkActor::New();
            vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
            mapper->SetInput(data);
            actor->SetMapper(mapper);
            mapper->Delete();
        }



    };

    class Sphere : public Decoration {
    public:
        Sphere(Real radius) {
            vtkSphereSource *sphere = vtkSphereSource::New();
            sphere->SetRadius(radius);
            sphere->SetThetaResolution(18);
            sphere->SetPhiResolution(18);
            sphere->Update();
            setPolyData(sphere->GetOutput());
            sphere->Delete();
        }
    };


public:

    void addDecoration(int bodyNum, Decoration& d, Transform X_GD) {
        addActor(bodyNum, d.getActor(), X_GD);
    }

    void addActor(int bodyNum, vtkActor* a, Transform X_GA) {
        ActorList& actors = bodies[bodyNum];
        vtkActor* acopy = vtkActor::New();
        acopy->ShallowCopy(a);
        actors.push_back(BodyActor(acopy,X_GA));
        renderer->AddActor(actors[actors.size()-1].first);
        if (bodyNum==0) setConfiguration(0, Transform()); // just do this once
    }

    ~VTKReporter() {deletePointers();}

    VTKReporter(const MySimbodyPendulum& m) 
      : mbs(m) {
        zeroPointers();
        makeShapes();

        renWin = vtkRenderWindow::New();
        renWin->SetSize(600,600);
        
        // an interactor
        vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
        iren->SetRenderWindow(renWin); 
        vtkInteractorStyleTrackballCamera* style=vtkInteractorStyleTrackballCamera::New(); 
        iren->SetInteractorStyle(style);
        style->Delete();
        iren->Initialize(); // register interactor to pick up windows messages


        renderer = vtkRenderer::New();
        renderer->SetBackground(1,1,1); // white
        renderer->GetActiveCamera()->SetFocalPoint(0,0,0);
        renderer->GetActiveCamera()->SetPosition(0,1,5);


        renWin->AddRenderer(renderer);

        const SimbodySubsystem& sbs = mbs.getSimbodySubsystem();
        bodies.resize(sbs.getNBodies());

        for (int i=0; i<(int)bodies.size(); ++i) {
            vtkActor *aAxes = vtkActor::New();
            aAxes->SetMapper(axesMapper);
            aAxes->GetProperty()->SetColor(0,i==0?1:0,0); // green for ground, else black
            aAxes->GetProperty()->SetLineWidth(3);
            addActor(i, aAxes, Transform());
            aAxes->Delete();

            vtkActor *aSphere = vtkActor::New();
            //aSphere->SetMapper(sphereMapper);
            aSphere->SetMapper(cubeMapper);
            aSphere->GetProperty()->SetColor(0,0,1); // sphere color blue
            //aSphere->GetProperty()->SetRepresentationToWireframe();
            aSphere->GetProperty()->SetOpacity(0.3);
            aSphere->SetScale(0.1);

            addActor(i, aSphere, Transform()); // at the origin

            aSphere->Delete();
        }


        renWin->Render();
    }

    void report(const State& s) {
        if (!renWin) return;

        const SimbodySubsystem& sbs = mbs.getSimbodySubsystem();
        for (int i=1; i<sbs.getNBodies(); ++i) {
            const Transform& config = sbs.getBodyConfiguration(s, i);
            setConfiguration(i, config);
        }

        renWin->Render();

        // Process any window messages since last time
        //TODO: Win32 specific
        MSG msg;
        bool done = false;
        while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE)) {
            if (msg.message == WM_QUIT) {
                printf("quit!!\n");
                renWin->Delete();renWin=0;
                break;
            }
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }

    }
private:

    const MySimbodyPendulum& mbs;
    std::vector<ActorList> bodies;
    vtkRenderWindow* renWin;
    vtkRenderer*     renderer;

    vtkPolyDataMapper *sphereMapper;
    vtkPolyDataMapper *cubeMapper;
    vtkPolyDataMapper *lineMapper;
    vtkPolyDataMapper *cylinderMapper;
    vtkPolyDataMapper *axesMapper;

    void zeroPointers() {
        renWin=0; renderer=0; sphereMapper=0; cubeMapper=0;
        lineMapper=0; cylinderMapper=0; axesMapper=0;
    }
    void deletePointers() {
        for (int i=0; i<(int)bodies.size(); ++i) {
            ActorList& actors = bodies[i];
            for (int a=0; a<(int)actors.size(); ++a)
                actors[a].first->Delete(), actors[a].first=0;
        }
        if(axesMapper)axesMapper->Delete(); 
        if(cylinderMapper)cylinderMapper->Delete(); 
        if(lineMapper)lineMapper->Delete();
        if(cubeMapper)cubeMapper->Delete(); 
        if(sphereMapper)sphereMapper->Delete();
        if(renderer)renderer->Delete();
        if(renWin)renWin->Delete();
        zeroPointers();
    }
    void makeShapes() {
        // sphere
        vtkSphereSource *sphere = vtkSphereSource::New();
        sphere->SetRadius(1.);
        sphere->SetThetaResolution(18);
        sphere->SetPhiResolution(18);
        sphereMapper = vtkPolyDataMapper::New();
        sphereMapper->SetInput(sphere->GetOutput());
        sphere->Delete();

        // cube
        vtkCubeSource *cube = vtkCubeSource::New();
        cube->SetXLength(1.); cube->SetYLength(1.); cube->SetZLength(1.);
        cubeMapper = vtkPolyDataMapper::New();
        cubeMapper->SetInput(cube->GetOutput());
        cube->Delete();    

        // line
        vtkLineSource *line = vtkLineSource::New();
        line->SetPoint1(0,0,0);
        line->SetPoint2(1,0,0);
        lineMapper = vtkPolyDataMapper::New();
        lineMapper->SetInput(line->GetOutput());
        line->Delete();

        // axes
        vtkAppendPolyData* app = vtkAppendPolyData::New();
        line = vtkLineSource::New(); line->SetPoint1(0,0,0); line->SetPoint2(1,0,0);
        app->AddInput(line->GetOutput());
        line->Delete();
        //cube = vtkCubeSource::New(); 
        //cube->SetXLength(1); cube->SetYLength(0.01); cube->SetZLength(0.01);
        //cube->SetCenter(.5,0,0);
        //app->AddInput(cube->GetOutput());
        //cube->Delete();
        line = vtkLineSource::New(); line->SetPoint1(0,0,0); line->SetPoint2(0,1,0);
        app->AddInput(line->GetOutput());
        line->Delete();
        line = vtkLineSource::New(); line->SetPoint1(0,0,0); line->SetPoint2(0,0,1);
        app->AddInput(line->GetOutput());
        line->Delete();

        vtkVectorText* xtext = vtkVectorText::New();
        xtext->SetText("x"); // default size is around 1
        vtkTransform* t = vtkTransform::New();
        t->Translate(1,-0.025,0); t->Scale(.1,.1,1); 
        vtkTransformPolyDataFilter* trf = vtkTransformPolyDataFilter::New();
        trf->SetInput(xtext->GetOutput());
        trf->SetTransform(t); 
        t->Delete();
        app->AddInput(trf->GetOutput());
        trf->Delete();

        xtext->Delete();


        axesMapper = vtkPolyDataMapper::New();
        axesMapper->SetInput(app->GetOutput());
        app->Delete();
    }
    void setConfiguration(int bodyNum, const Transform& X_GB) {
        const ActorList& actors = bodies[bodyNum];
        for (int i=0; i < (int)actors.size(); ++i) {
            vtkProp3D*       actor = actors[i].first;
            const Transform& X_BA  = actors[i].second;
            const Transform  X_GA  = X_GB*X_BA;
            actor->SetPosition(X_GA.T()[0], X_GA.T()[1], X_GA.T()[2]);
            const Vec4 av = X_GA.R().convertToAngleAxis();
            actor->SetOrientation(0,0,0);
            actor->RotateWXYZ(av[0]/RadiansPerDegree, av[1], av[2], av[3]);
        }

    }

};

static void drawCone();
static void drawSphere();

int main(int argc, char** argv) {
    //drawSphere(); 
    //return 0;

    std::vector<State> saveEm;

    try { // If anything goes wrong, an exception will be thrown.
        Real start = initialTheta;
        if (argc > 1) sscanf(argv[1], "%lg", &start);
        printf("Pendulum starting at angle +%g degrees from vertical.\n", start);

        // Create a multibody system using Simbody.
        MySimbodyPendulum myPend;
        myPend.setPendulumAngle(start);

        // And a study using the Runge Kutta Merson integrator
        RungeKuttaMerson myStudy(myPend);
        myStudy.setAccuracy(1e-4);

        VTKReporter display(myPend);
        display.addDecoration(0,VTKReporter::Sphere(0.01),
            Transform(Vec3(0.5,-0.2,0.1)));

        // Run for 5 periods without output every dt seconds,
        // starting at theta=start degrees.

        const Real dt = 0.01; // output intervals

        printf("time  theta (deg)  (period should be %gs)\n", expectedPeriod);

        myStudy.setInitialConditions(myPend.getT(), myPend.getY());
        display.report(myPend.getState());
        for (;;) {
            printf("%5g %10.3g\n", myStudy.getT(), myPend.getPendulumAngle());
            display.report(myPend.getState());
            saveEm.push_back(myPend.getState());

           // if (myStudy.getT() >= 10*expectedPeriod)
             //   break;
    
            if (myStudy.getT() >= 10.)
                break;

            // TODO: should check for errors or have or teach RKM to throw. 
            myStudy.step(myStudy.getT() + dt);
        }

        for (int i=0; i < (int)saveEm.size(); ++i)
            display.report(saveEm[i]);
    } 
    catch (const exception& e) {
        printf("EXCEPTION THROWN: %s\n", e.what());
        exit(1);
    }
}



// Callback for interaction
class vtkMyCallback : public vtkCommand {
public:
    static vtkMyCallback* New() {return new vtkMyCallback;}

    /*virtual*/ void Execute(vtkObject *caller, unsigned long, void*) {
        vtkRenderer* ren =
            reinterpret_cast<vtkRenderer*>(caller);
        cout << ren->GetActiveCamera()->GetPosition()[0] << " ";
        cout << ren->GetActiveCamera()->GetPosition()[1] << " ";
        cout << ren->GetActiveCamera()->GetPosition()[2] << endl;
    }
};

static void drawCone() {
    vtkConeSource* cone = vtkConeSource::New();
    cone->SetHeight(0.5);
    cone->SetRadius(.15);
    cone->SetResolution(10);

    vtkPolyDataMapper* coneMapper = vtkPolyDataMapper::New();
    coneMapper->SetInput(cone->GetOutput());
    vtkActor* coneActor = vtkActor::New();
    coneActor->SetMapper(coneMapper);

    vtkRenderer* ren1 = vtkRenderer::New();
    ren1->AddActor(coneActor);
    ren1->SetBackground(.1,.2,.4);

    vtkRenderWindow* renWin = vtkRenderWindow::New();
    renWin->AddRenderer(ren1);
    renWin->SetSize(300,300);

    vtkMyCallback* mo1 = vtkMyCallback::New();
    ren1->AddObserver(vtkCommand::StartEvent,mo1);
    mo1->Delete();

    for (int i=0; i<360; ++i) {
        //render the image
        renWin->Render();
        // rotate the active camera by one degree
        ren1->GetActiveCamera()->Azimuth(1);
    }

    cone->Delete();
    coneMapper->Delete();
    coneActor->Delete();
    ren1->Delete();
    renWin->Delete();
}


static void drawSphere ()
{
    // a renderer & a window to put it in
    vtkRenderer *ren1 = vtkRenderer::New();
    vtkRenderWindow *renWin = vtkRenderWindow::New();
    renWin->AddRenderer(ren1);
    renWin->SetSize(600,600);

    // create sphere geometry
    vtkSphereSource *sphere = vtkSphereSource::New();
    sphere->SetRadius(1.0);
    sphere->SetThetaResolution(18);
    sphere->SetPhiResolution(18);

    // toss in a cone too
    vtkConeSource* cone = vtkConeSource::New();
    cone->SetHeight(3);
    cone->SetRadius(1);
    cone->SetResolution(10);

        // and a line?
    vtkLineSource* line = vtkLineSource::New();
    line->SetPoint1(0,0,0);
    line->SetPoint2(3,0,0);

    //vtkTextSource* xtext = vtkTextSource::New();
    vtkVectorText* xtext = vtkVectorText::New();

    xtext->SetText("x");
    //xtext->SetBacking(0);
    //xtext->SetForegroundColor(0,1,0);
    vtkPolyDataMapper* xtextMapper = vtkPolyDataMapper::New();
    xtextMapper->SetInput(xtext->GetOutput());

    // map to graphics library
    vtkPolyDataMapper *sphereMapper = vtkPolyDataMapper::New();
    sphereMapper->SetInput(sphere->GetOutput());

    vtkPolyDataMapper* coneMapper = vtkPolyDataMapper::New();
    coneMapper->SetInput(cone->GetOutput());

    vtkPolyDataMapper* lineMapper = vtkPolyDataMapper::New();
    lineMapper->SetInput(line->GetOutput());

    // actor coordinates geometry, properties, transformation
    vtkActor *aSphere = vtkActor::New();
    aSphere->SetMapper(sphereMapper);
    aSphere->GetProperty()->SetColor(0,0,1); // sphere color blue
    //aSphere->GetProperty()->SetRepresentationToWireframe();
    aSphere->GetProperty()->SetOpacity(0.5);
    aSphere->SetPosition(-1, 0, 0);

    vtkActor* aCone = vtkActor::New();
    aCone->SetMapper(coneMapper);
    aCone->GetProperty()->SetColor(1,0,0); // cone color red
    aCone->SetPosition(1, 0, 0);

    vtkActor* xLine = vtkActor::New();
    xLine->SetMapper(lineMapper);
    xLine->GetProperty()->SetColor(0,1,0); // green line
    xLine->SetPosition(0, 0, 0);

    //vtkActor* xtextActor = vtkActor::New();
    vtkFollower* xtextActor = vtkFollower::New();
    xtextActor->SetCamera(ren1->GetActiveCamera());
    xtextActor->SetMapper(xtextMapper);
    //xtextActor->SetOrigin( xtextActor->GetCenter() );
    //xtextActor->SetPosition(3,0,0);
    xtextActor->AddPosition(3,0,0);
    xtextActor->SetScale(.5);
    xtextActor->GetProperty()->SetColor(0,1,0);
    xtextActor->GetProperty()->SetSpecular(0.5);

    vtkAssembly* xlabeledAsm = vtkAssembly::New();
    xlabeledAsm->AddPart(xLine);
    //xlabeledAsm->AddPart(xtextActor);

    vtkActor* yLine = vtkActor::New();
    yLine->SetMapper(lineMapper);
    yLine->GetProperty()->SetColor(0,1,0); // green line
    yLine->SetPosition(0, 0, 0);
    yLine->RotateZ(90);

    vtkActor* zLine = vtkActor::New();
    zLine->SetMapper(lineMapper);
    zLine->GetProperty()->SetColor(0,1,0); // green line
    zLine->SetPosition(0, 0, 0);
    zLine->RotateY(-90);



    // an interactor
    vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(renWin);

    vtkAssembly* assembly = vtkAssembly::New();
    assembly->AddPart(aSphere);
    assembly->AddPart(aCone);

    // add the actors to the scene
    //ren1->AddActor(aSphere);
    //ren1->AddActor(aCone);
    ren1->AddActor(assembly);
    //ren1->AddActor(xLine);
    ren1->AddActor(xlabeledAsm);
    ren1->AddActor(xtextActor);
    ren1->AddActor(yLine);
    ren1->AddActor(zLine);
    ren1->SetBackground(1,1,1); // Background color white


    ren1->GetActiveCamera()->SetFocalPoint(0,0,0);
    ren1->GetActiveCamera()->SetPosition(0,3,20);

    // render an image (lights and cameras are created automatically)
    renWin->Render();

    // begin mouse interaction
    iren->Start();
}



