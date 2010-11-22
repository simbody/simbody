/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/Visualizer.h"
#include "simbody/internal/VisualizationEventListener.h"
#include "simbody/internal/DecorationGenerator.h"
#include "VisualizationGeometry.h"
#include "VisualizationProtocol.h"

#include <cstdlib>
#include <cstdio>
#include <pthread.h>
#include <string>

#include <ctime>

using namespace SimTK;
using namespace std;

class Visualizer::RubberBandLine {
public:
    RubberBandLine(MobilizedBodyIndex b1, const Vec3& station1, MobilizedBodyIndex b2, const Vec3& station2, const DecorativeLine& line) :
            b1(b1), station1(station1), b2(b2), station2(station2), line(line) {
    }
    MobilizedBodyIndex b1, b2;
    Vec3 station1, station2;
    DecorativeLine line;
};

static void* drawFramesWhenReady(void* visualizerAsVoidp);

static const Real DefaultFrameRateFPS = 30;

class Visualizer::VisualizerRep {
public:
    VisualizerRep(Visualizer* owner, MultibodySystem& system, 
                  const String& title) 
    :   handle(owner), system(system), protocol(*owner, title),
        mode(PassThrough), nextFrameDueInNs(-1), timeBetweenFramesInNs(-1),
        frameRateFPS(0), oldest(0),nframe(0),capacity(0) 
    {   
        pthread_mutex_init(&queueLock, NULL); 
        pthread_cond_init(&queueNotFull, NULL); 

        setMode(PassThrough);
    }

    void setMode(Mode newMode) {
        timeBetweenFramesInNs = secToNs(1/frameRateFPS);
        if (newMode==mode) {
            nextFrameDueInNs = realTimeInNs(); // i.e., now
            return;
        }
        // Mode is changing. If it was RealTime before we have
        // to clean up first.
        if (mode == RealTime) {
            pthread_cancel(drawThread);
            initializePool(0);  // existing frames are lost
        }


        // If the new mode is RealTime we have to allocate the
        // queue and initiate the buffered rendering thread.
        // TODO: calculate pool size
        if (newMode == RealTime) {
            initializePool(10);
            pthread_create(&drawThread, NULL, drawFramesWhenReady, this);
        }

        mode = newMode; // change mode
        setDesiredFrameRate(frameRateFPS);

        // Note that the next frame is due now.
        nextFrameDueInNs = realTimeInNs(); // i.e., now
    }

    void setDesiredFrameRate(Real framesPerSec) {
        framesPerSec = std::max(framesPerSec, Real(0));
        frameRateFPS = framesPerSec;
        if (frameRateFPS > 0)
            timeBetweenFramesInNs = secToNs(1/frameRateFPS);
        else { // use default frame rate for the current mode
            if (mode == PassThrough)
                timeBetweenFramesInNs = 0LL;
            else 
                timeBetweenFramesInNs = secToNs(1/DefaultFrameRateFPS);
        }
    }

    ~VisualizerRep() {
        for (int i = 0; i < (int) listeners.size(); i++)
            delete listeners[i];
        for (int i = 0; i < (int) generators.size(); i++)
            delete generators[i];
        if (mode==RealTime)
            pthread_cancel(drawThread);
        pthread_cond_destroy(&queueNotFull);
        pthread_mutex_destroy(&queueLock);
    }

    void drawFrameNow(const State& state);

    void initializePool(int sz) {
        pthread_mutex_lock(&queueLock);
        pool.resize(sz);oldest=nframe=0;capacity=sz;
        pthread_mutex_unlock(&queueLock);
    }
    int getQueueCapacity() const {return capacity;} // varies
    int getNFramesInQueue() const {return nframe;}

    bool queuingIsEnabled() const {return pool.size() > 0;}
    bool queueIsFull() const {return nframe==capacity;}
    bool queueIsEmpty() const {return nframe==0;}

    // Called from simulation thread.
    bool tryAddStateToQueue(const State& state) {
        bool added = false;
        pthread_mutex_lock(&queueLock);
        if (!queueIsFull()) {
            pool[(oldest+nframe)%pool.size()] = state;
            nframe++;
            added = true;
        }
        pthread_mutex_unlock(&queueLock);
        return added;
    }
    // Called from simulation thread.
    void addStateToQueueWithWait(const State& state) {
        pthread_mutex_lock(&queueLock);
        while(queueIsFull()) {
            // atomic unlock, long wait, relock
            pthread_cond_wait(&queueNotFull,&queueLock);
        }
        // queue is not full
        pool[(oldest+nframe)%pool.size()] = state;
        nframe++;
        pthread_mutex_unlock(&queueLock);
    }

    // Visualization thread uses this to remove the oldest frame
    // from active duty in the pool. It may then at its leisure
    // use the contained State to generate a screen image. When
    // done it should return the frame to the pool.
    int removeOldestFrameFromPool() {
        int frame = -1;
        pthread_mutex_lock(&queueLock);
        if (nframe) {
            frame = oldest;
            oldest = (oldest+1)%pool.size();
            --capacity; // we're borrowing this frame for now
            --nframe;
        }
        pthread_mutex_unlock(&queueLock);
        return frame;
    }

    void returnRemovedFrameToPool() {
        pthread_mutex_lock(&queueLock);
        ++capacity; 
        if (capacity-nframe >= (capacity/2)+1)
            pthread_cond_signal(&queueNotFull);
        pthread_mutex_unlock(&queueLock);
    }


    Visualizer*                 handle;
    MultibodySystem&            system;
    VisualizationProtocol       protocol;
    Array_<DecorativeGeometry>  addedGeometry;
    Array_<RubberBandLine>      lines;
    Array_<VisualizationEventListener*>     listeners;
    Array_<DecorationGenerator*>            generators;

    // Operating mode.
    Mode mode;

    // Real time after which we can render a frame.
    long long nextFrameDueInNs;
    // How many nanoseconds between frames?
    long long timeBetweenFramesInNs;

    Real frameRateFPS; // in frames/sec if > 0, else use default

    // The frame queue:
    Array_<State,int> pool; // fixed size, old to new order but circular
    int oldest, nframe, capacity; // oldest is index into pool, nframe is #valid entries
    pthread_mutex_t queueLock;
    pthread_cond_t queueNotFull; // tied to queueLock
    pthread_t drawThread;
};

void Visualizer::VisualizerRep::drawFrameNow(const State& state) {
    system.realize(state, Stage::Position);
    Array_<DecorativeGeometry> geometry;
    for (Stage stage = Stage::Topology; stage <= state.getSystemStage(); ++stage)
        system.calcDecorativeGeometryAndAppend(state, stage, geometry);
    for (int i = 0; i < (int) generators.size(); i++)
        generators[i]->generateDecorations(state, geometry);
    protocol.beginScene();
    VisualizationGeometry geometryCreator(protocol, system.getMatterSubsystem(), state);
    for (int i = 0; i < (int) geometry.size(); ++i)
        geometry[i].implementGeometry(geometryCreator);
    for (int i = 0; i < (int) addedGeometry.size(); ++i)
        addedGeometry[i].implementGeometry(geometryCreator);
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
    for (int i = 0; i < (int) lines.size(); ++i) {
        const RubberBandLine& line = lines[i];
        Vec3 end1 = matter.getMobilizedBody(line.b1).getBodyTransform(state)*line.station1;
        Vec3 end2 = matter.getMobilizedBody(line.b2).getBodyTransform(state)*line.station2;
        Real thickness = line.line.getLineThickness() == -1 ? 1 : line.line.getLineThickness();
        protocol.drawLine(end1, end2, VisualizationGeometry::getColor(line.line), thickness);
    }
    protocol.finishScene();
}

// This is the main function for the buffered rendering thread. Its job is to 
// pull the oldest frames off the queue and send them to the renderer at
// the right times.
static void* drawFramesWhenReady(void* visualizerRepAsVoidp) {
    Visualizer::VisualizerRep& vizRep = 
        *reinterpret_cast<Visualizer::VisualizerRep*>(visualizerRepAsVoidp);
    while (true) {
        // Wait for the next frame time.
        long long now;
        while ((now = realTimeInNs()) < vizRep.nextFrameDueInNs) {
            long long waitInNs = vizRep.nextFrameDueInNs - now;
            timespec delayTime;
            nsToTimespec(waitInNs, delayTime);
            nanosleep(&delayTime,0);
        }

        // A frame is due. Grab the oldest frame if available and render it. 
        // TODO: should wait on a condition variable for a frame; instead
        // we'll sleep for a fraction of the frame time and retry.
        int frame;
        while ((frame=vizRep.removeOldestFrameFromPool()) < 0) {
            timespec delayTime;
            nsToTimespec(vizRep.timeBetweenFramesInNs/2, delayTime);
            nanosleep(&delayTime,0);
            now = realTimeInNs();
        }

        // Got a frame to render.
        const State& state = vizRep.pool[frame];
        // TODO: replace real time here with simulated time
        vizRep.nextFrameDueInNs = now + vizRep.timeBetweenFramesInNs;
        vizRep.drawFrameNow(state);
        vizRep.returnRemovedFrameToPool();
    }
    return 0;
}

Visualizer::Visualizer(MultibodySystem& system) : rep(0) {
    // Create a default title from the name of this executable.
    bool isAbsolutePath;
    std::string directory, fileName, extension;
    Pathname::deconstructPathname(Pathname::getThisExecutablePath(),
        isAbsolutePath, directory, fileName, extension);
    rep = new VisualizerRep(this, system, fileName);
}

Visualizer::Visualizer(MultibodySystem& system, const String& title) : rep(0) {
    rep = new VisualizerRep(this, system, title);
}

Visualizer::~Visualizer() {
    if (rep->handle == this)
        delete rep;
}

void Visualizer::setMode(Visualizer::Mode mode) {updRep().setMode(mode);}
Visualizer::Mode Visualizer::getMode() const {return getRep().mode;}

void Visualizer::setDesiredFrameRate(Real framesPerSec) {
    updRep().setDesiredFrameRate(std::max(framesPerSec, Real(0)));
}
Real Visualizer::getDesiredFrameRate() const {return getRep().frameRateFPS;}

void Visualizer::drawFrameNow(const State& state) {
    updRep().drawFrameNow(state);
}


bool Visualizer::tryReport(const State& state) {
    if (!getRep().queuingIsEnabled()) {
        if (getRep().timeBetweenFramesInNs == 0LL) {
            drawFrameNow(state);
            return true;
        }

        long long now;
        while ((now=realTimeInNs()) < getRep().nextFrameDueInNs) {
            if (getRep().mode == Sampling)
                return true; // Sampling: frame too early; ignore.

            // PassThrough: wait until frame time, then send.
            const long long waitInNs = getRep().nextFrameDueInNs - now;
            timespec delayTime;
            nsToTimespec(waitInNs, delayTime);
            nanosleep(&delayTime,0);
        }

        // Frame time reached in PassThrough mode. This frame might 
        // be on time or late; we'll schedule the next time for one
        // frame time later to keep the maximum rate down to the 
        // specified rate (otherwise a late frame could be followed
        // by lots of fast frames).
        updRep().nextFrameDueInNs = now + getRep().timeBetweenFramesInNs;
        drawFrameNow(state);
        return true;
    }
   // RealTime mode.
   return updRep().tryAddStateToQueue(state);
}

void Visualizer::report(const State& state) {
    if (tryReport(state))
        return; 
    // Queue full.
    updRep().addStateToQueueWithWait(state);
}

void Visualizer::addEventListener(VisualizationEventListener* listener) {
    updRep().listeners.push_back(listener);
}

const Array_<VisualizationEventListener*>& Visualizer::getEventListeners() const {
    return getRep().listeners;
}
void Visualizer::addDecorationGenerator(DecorationGenerator* generator) {
    updRep().generators.push_back(generator);
}

void Visualizer::addMenu(const String& title, const Array_<pair<String, int> >& items) {
    updRep().protocol.addMenu(title, items);
}

void Visualizer::addDecoration(MobilizedBodyIndex mobodIx, const Transform& X_BD, const DecorativeGeometry& geom) {
    Array_<DecorativeGeometry>& addedGeometry = updRep().addedGeometry;
    addedGeometry.push_back(geom);
    DecorativeGeometry& geomCopy = addedGeometry.back();
    geomCopy.setBodyId((int)mobodIx);
    geomCopy.setTransform(X_BD * geomCopy.getTransform());
}

void Visualizer::addRubberBandLine(MobilizedBodyIndex b1, const Vec3& station1, MobilizedBodyIndex b2, const Vec3& station2, const DecorativeLine& line) {
    updRep().lines.push_back(RubberBandLine(b1, station1, b2, station2, line));
}

void Visualizer::setCameraTransform(const Transform& transform) {
    updRep().protocol.setCameraTransform(transform);
}

void Visualizer::zoomCameraToShowAllGeometry() {
    updRep().protocol.zoomCamera();
}

void Visualizer::pointCameraAt(const Vec3& point, const Vec3& upDirection) {
    updRep().protocol.lookAt(point, upDirection);
}

void Visualizer::setCameraFieldOfView(Real fov) {
    updRep().protocol.setFieldOfView(fov);
}

void Visualizer::setCameraClippingPlanes(Real nearPlane, Real farPlane) {
    updRep().protocol.setClippingPlanes(nearPlane, farPlane);
}

void Visualizer::setGroundPosition(const CoordinateAxis& axis, Real height) {
    updRep().protocol.setGroundPosition(axis, height);
}
