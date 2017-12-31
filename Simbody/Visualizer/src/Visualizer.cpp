/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-14 Stanford University and the Authors.        *
 * Authors: Peter Eastman, Michael Sherman                                    *
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

#include "simbody/internal/common.h"
#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/Visualizer.h"
#include "simbody/internal/Visualizer_InputListener.h"

#include "VisualizerGeometry.h"
#include "VisualizerProtocol.h"

#include <cstdlib>
#include <cstdio>
#include <string>
#include <ctime>
#include <iostream>
#include <limits>
#include <condition_variable>

using namespace SimTK;
using namespace std;

static void drawingThreadMain(Visualizer::Impl& vizImpl);

static const long long UsToNs = 1000LL;          // ns = us * UsToNs
static const long long MsToNs = 1000LL * UsToNs; // ns = ms * MsToNs

static const Real      DefaultFrameRateFPS             = 30;
static const Real      DefaultDesiredBufferLengthInSec = Real(0.15); // 150ms

// These are not currently overrideable.
static const long long DefaultAllowableFrameJitterInNs      = 5 * MsToNs; //5ms
static const Real      DefaultSlopAsFractionOfFrameInterval = Real(0.05); //5%

namespace { // local classes
//==============================================================================
//                              VISUALIZER IMPL
//==============================================================================
/* This is the private implementation object contained in a Visualizer handle.
See the "implementation notes" section of the Visualizer class documentation
for some information about how this works. */

// If we are buffering frames, this is the object that represents a frame
// in the queue. It consists of a copy of a reported State, and a desired
// draw time for the frame, in adjusted real time (AdjRT).
struct Frame {
    Frame() : desiredDrawTimeAdjRT(-1LL) {}
    Frame(const State& state, const long long& desiredDrawTimeAdjRT)
    :   state(state), desiredDrawTimeAdjRT(desiredDrawTimeAdjRT) {}
    // default copy constructor, copy assignment, destructor

    void clear() {desiredDrawTimeAdjRT = -1LL;}
    bool isValid() const {return desiredDrawTimeAdjRT >= 0LL;}

    State       state;
    long long   desiredDrawTimeAdjRT; // in adjusted real time
};

// This holds the specs for rubber band lines that are added directly
// to the Visualizer.
class RubberBandLine {
public:
    RubberBandLine(MobilizedBodyIndex b1, const Vec3& station1, 
                   MobilizedBodyIndex b2, const Vec3& station2, 
                   const DecorativeLine& line) 
    :   b1(b1), station1(station1), b2(b2), station2(station2), line(line) {}

    MobilizedBodyIndex  b1, b2;
    Vec3                station1, station2;
    DecorativeLine      line;
};
} // end local namespace

// Implementation of the Visualizer.
class Visualizer::Impl {
public:
    // Create a Visualizer and put it in PassThrough mode.
    Impl(Visualizer* owner, const MultibodySystem& system,
         const Array_<String>& searchPath) 
    :   m_system(system), m_protocol(*owner, searchPath),
        m_shutdownWhenDestructed(false), m_upDirection(YAxis), m_groundHeight(0),
        m_mode(PassThrough), m_frameRateFPS(DefaultFrameRateFPS), 
        m_simTimeUnitsPerSec(1), 
        m_desiredBufferLengthInSec(DefaultDesiredBufferLengthInSec), 
        m_timeBetweenFramesInNs(secToNs(1/DefaultFrameRateFPS)),
        m_allowableFrameJitterInNs(DefaultAllowableFrameJitterInNs),
        m_allowableFrameTimeSlopInNs(
            secToNs(DefaultSlopAsFractionOfFrameInterval/DefaultFrameRateFPS)),
        m_adjustedRealTimeBase(realTimeInNs()),
        m_prevFrameSimTime(-1), m_nextFrameDueAdjRT(-1), 
        m_oldest(0),m_nframe(0),
        m_drawThreadIsRunning(false), m_drawThreadShouldSuicide(false),
        m_refCount(0)
    {   
        setMode(PassThrough);
        clearStats();

        m_protocol.setMaxFrameRate(m_frameRateFPS);
        m_protocol.setBackgroundColor(White);
        m_protocol.setBackgroundType(system.getUseUniformBackground() 
                                        ? SolidColor : GroundAndSky);
        m_protocol.setSystemUpDirection(system.getUpDirection());
    }
    
    ~Impl() {
        if (m_mode==RealTime && m_pool.size()) {
            killDrawThreadIfNecessary();
        }
        for (unsigned i = 0; i < m_controllers.size(); i++)
            delete m_controllers[i];
        for (unsigned i = 0; i < m_listeners.size(); i++)
            delete m_listeners[i];
        for (unsigned i = 0; i < m_generators.size(); i++)
            delete m_generators[i];

        if (m_shutdownWhenDestructed) {
            try {
                // This throws an exception if the pipe is broken (e.g., if the
                // simbody-visualizer has already been shut down).
                m_protocol.shutdownGUI();
            } catch (...) {}
        }
    }

    void setShutdownWhenDestructed(bool shouldShutdown)
    {   m_shutdownWhenDestructed = shouldShutdown; }

    bool getShutdownWhenDestructed() const
    {   return m_shutdownWhenDestructed; }

    // Call from simulation thread.
    void startDrawThread() {
        SimTK_ASSERT_ALWAYS(!m_drawThreadIsRunning,
            "Tried to start the draw thread when it was already running.");
        m_drawThreadShouldSuicide = false;
        m_drawThread = std::thread(drawingThreadMain, std::ref(*this));
        m_drawThreadIsRunning = true;
    }

    // Call from simulation thread.
    void killDrawThread() {
        SimTK_ASSERT_ALWAYS(m_drawThreadIsRunning,
            "Tried to kill the draw thread when it wasn't running.");
        m_drawThreadShouldSuicide = true;
        // The draw thread might be waiting on an empty queue, in which
        // case we have to wake it up (see getOldestFrameInQueue()).
        // We'll do it twice 100ms apart to avoid a timing issue where
        // we signal just before the thread waits.
        m_queueNotEmpty.notify_one(); // wake it if necessary
        sleepInSec(0.1); 
        m_queueNotEmpty.notify_one();
        m_drawThread.join(); // wait for death
        m_drawThreadIsRunning = m_drawThreadShouldSuicide = false;
    }

    void startDrawThreadIfNecessary() 
    {   if (!m_drawThreadIsRunning) startDrawThread(); }

    void killDrawThreadIfNecessary() 
    {   if (m_drawThreadIsRunning) killDrawThread(); }

    // Whenever a "set" method is called that may change one of the 
    // interrelated time quantities, set all of them. We expect
    // the mode to have been set already.
    void resetTimeRelatedQuantities(Real framesPerSec,
        Real timeScale, Real desiredBufLengthSec) 
    {
        if (framesPerSec <= 0) framesPerSec = DefaultFrameRateFPS;
        if (timeScale <= 0)    timeScale = 1;
        if (desiredBufLengthSec < 0) 
            desiredBufLengthSec = DefaultDesiredBufferLengthInSec;

        // Frame rate.
        m_frameRateFPS               = framesPerSec;
        m_timeBetweenFramesInNs      = secToNs(1/m_frameRateFPS);
        m_allowableFrameTimeSlopInNs = 
            secToNs(DefaultSlopAsFractionOfFrameInterval/m_frameRateFPS);
        m_allowableFrameJitterInNs   = DefaultAllowableFrameJitterInNs;

        // Time scale.
        m_simTimeUnitsPerSec = timeScale;

        // Realtime buffer.
        m_desiredBufferLengthInSec = desiredBufLengthSec;

        int numFrames = 
            (int)(m_desiredBufferLengthInSec/nsToSec(m_timeBetweenFramesInNs) 
                  + 0.5);
        if (numFrames==0 && m_desiredBufferLengthInSec > 0)
            numFrames = 1;

        // If we're in RealTime mode and we have changed the number of
        // frames in the buffer, reallocate the pool and kill or start
        // the draw thread if necessary.
        if (m_mode == RealTime && numFrames != m_pool.size()) {
            if (m_pool.size()) {
                // draw thread isn't needed if we get rid of the buffer
                if (numFrames == 0)
                    killDrawThreadIfNecessary();
                initializePool(numFrames);
            } else {
                // draw thread is needed if we don't have one
                initializePool(numFrames);
                startDrawThreadIfNecessary();
            }
        }

        clearStats();

        // Note that the next frame we see is the first one and we'll need
        // to initialize adjusted real time then.
        m_nextFrameDueAdjRT = -1LL; // i.e., now
    }

    void setMode(Visualizer::Mode newMode) {
        // If we're not changing modes we just clear the stats and invalidate
        // the next expected frame time so that we'll take the first one that
        // shows up.
        if (newMode==m_mode) {
            resetTimeRelatedQuantities(m_frameRateFPS,
                                       m_simTimeUnitsPerSec,
                                       m_desiredBufferLengthInSec);
            return;
        }

        // Mode is changing. If it was buffered RealTime before we have
        // to clean up first.
        if (m_mode == RealTime && m_pool.size()) {
            killDrawThreadIfNecessary();
            initializePool(0);  // clear the buffer
        }

        m_mode = newMode; // change mode
        resetTimeRelatedQuantities(m_frameRateFPS,
                                   m_simTimeUnitsPerSec,
                                   m_desiredBufferLengthInSec);
    }

    void setDesiredFrameRate(Real framesPerSec) {
        resetTimeRelatedQuantities(framesPerSec,
                                   m_simTimeUnitsPerSec,
                                   m_desiredBufferLengthInSec);
        // Make sure the GUI doesn't try to outrace us when it generates
        // its own frames.
        m_protocol.setMaxFrameRate(framesPerSec);
    }

    void setDesiredBufferLengthInSec(Real bufferLengthInSec) {
        resetTimeRelatedQuantities(m_frameRateFPS,
                                   m_simTimeUnitsPerSec,
                                   bufferLengthInSec);                               
    }

    void setRealTimeScale(Real simTimePerRealSec)  {
        resetTimeRelatedQuantities(m_frameRateFPS,
                                   simTimePerRealSec,
                                   m_desiredBufferLengthInSec);                               
    }

    Real getDesiredBufferLengthInSec() const 
    {   return m_desiredBufferLengthInSec; }

    int getActualBufferLengthInFrames() const {return m_pool.size();}
    Real getActualBufferLengthInSec() const 
    {   return (Real)nsToSec(getActualBufferLengthInFrames()
                             *m_timeBetweenFramesInNs); }

    // Generate this frame and send it immediately to the renderer without
    // thinking too hard about it.
    void drawFrameNow(const State& state);

    // In RealTime mode we have a frame to draw and a desired draw time in
    // AdjRT. Draw it when the time comes, and adjust AdjRT if necessary.
    void drawRealtimeFrameWhenReady
       (const State& state, const long long& desiredDrawTimeAdjRT);

    // Queuing is used only in RealTime mode.

    // Called from the simulation thread when it wants to report a frame
    // and we are in RealTime mode.
    void reportRealtime(const State& state);

    // Set the maximum number of frames in the buffer.
    void initializePool(int sz) {
        std::lock_guard<std::mutex> lock(m_queueMutex);
        m_pool.resize(sz);m_oldest=m_nframe=0;
    }

    int getNFramesInQueue() const {return m_nframe;}

    // Queing is enabled if the pool was allocated.
    bool queuingIsEnabled() const {return m_pool.size() != 0;}
    bool queueIsFull() const {return m_nframe==m_pool.size();}
    bool queueIsEmpty() const {return m_nframe==0;}

    // Called from simulation thread. Blocks until there is room in
    // the queue, then inserts this state unconditionally, with the indicated
    // desired rendering time in adjusted real time. We then update the 
    // "time of next queue slot" to be one ideal frame interval later than 
    // the desired draw time.
    void addFrameToQueueWithWait(const State& state, 
                                 const long long& desiredDrawTimeAdjRT)
    {
        std::unique_lock<std::mutex> lock(m_queueMutex);
        ++numReportedFramesThatWereQueued;
        if (queueIsFull()) {
            ++numQueuedFramesThatHadToWait;
            // atomic: unlock, long wait, relock
            // Only wake up if queue is not full (ignore spurious wakeups).
            m_queueNotFull.wait(lock, [&] {return !queueIsFull();});
        }

        // There is room in the queue now. We're holding the lock.
        Frame& frame = m_pool[(m_oldest+m_nframe)%m_pool.size()];
        frame.state  = state;
        frame.desiredDrawTimeAdjRT = desiredDrawTimeAdjRT;

        // Record the frame time.
        m_prevFrameSimTime = state.getTime();

        // Set the expected next frame time (in AdjRT).
        m_nextFrameDueAdjRT = desiredDrawTimeAdjRT + m_timeBetweenFramesInNs;

        if (++m_nframe == 1) 
            // wake up rendering thread on first frame
            m_queueNotEmpty.notify_one();

        lock.unlock();
    }

    // Call from simulation thread to allow the drawing thread to flush
    // any frames currently in the queue.
    void waitUntilQueueIsEmpty() {
        if (   !queuingIsEnabled() || m_nframe==0 
            || !m_drawThreadIsRunning || m_drawThreadShouldSuicide)
            return;
        std::unique_lock<std::mutex> lock(m_queueMutex);
        m_queueIsEmpty.wait(lock, [&] {return m_nframe == 0;});
        lock.unlock();
    }

    // The drawing thread uses this to find the oldest frame in the buffer.
    // It may then at its leisure use the contained State to generate a screen
    // image. There is no danger of the simulation thread modifying this
    // frame; once it has been put in it stays there until the drawing thread
    // takes it out. When done it should return the frame to the pool.
    // Returns true if it gets a frame (which will always happen in normal
    // operation since it waits until one is available), false if the draw
    // thread should quit.
    bool getOldestFrameInQueue(const Frame** fp) {
        std::unique_lock<std::mutex> lock(m_queueMutex);
        if (m_nframe == 0 && !m_drawThreadShouldSuicide) {
            ++numTimesDrawThreadBlockedOnEmptyQueue;
            // atomic: unlock, long wait, relock; ignore spurious wakeups.
            m_queueNotEmpty.wait(lock,
                    [&] {return m_nframe || m_drawThreadShouldSuicide;});
        } else {
            sumOfQueueLengths        += double(m_nframe);
            sumSquaredOfQueueLengths += double(square(m_nframe));
        }
        // There is at least one frame available now, unless we're supposed
        // to quit. We are holding the lock.
        lock.unlock();
        if (m_drawThreadShouldSuicide) {*fp=0; return false;}
        else {*fp=&m_pool[m_oldest]; return true;} // sim thread won't change oldest
    }

    // Drawing thread uses this to note that it is done with the oldest
    // frame which may now be reused by the simulation thread. The queueNotFull
    // condition is posted if there is a reasonable amount of room in the pool 
    // now.
    void noteThatOldestFrameIsNowAvailable() {
        std::unique_lock<std::mutex> lock(m_queueMutex);
        m_oldest = (m_oldest+1)%m_pool.size(); // move to next-oldest
        --m_nframe; // there is now one fewer frame in use
        if (m_nframe == 0)
            m_queueIsEmpty.notify_one(); // in case we're flushing
        // Start the simulation again when the pool is about half empty.
        if (m_nframe <= m_pool.size()/2+1)
            m_queueNotFull.notify_one();
        lock.unlock();
    }

    // Given a time t in simulation time units, return the equivalent time r in
    // seconds of real time. That is the amount of real time that should have
    // elapsed since t=0 if this simulation were running at exactly the desired
    // real time rate.
    long long convertSimTimeToNs(const double& t)
    {   return secToNs(t / m_simTimeUnitsPerSec); }

    // same as ns; that's what AdjRT tries to be
    long long convertSimTimeToAdjRT(const double& t)
    {   return convertSimTimeToNs(t); } 

    double convertAdjRTtoSimTime(const long long& a)
    {   return nsToSec(a) * m_simTimeUnitsPerSec; }

    long long convertRTtoAdjRT(const long long& r)
    {   return r - m_adjustedRealTimeBase; }

    long long convertAdjRTtoRT(const long long& a)
    {   return a + m_adjustedRealTimeBase; }


    // Adjust the real time base by a given signed offset in nanoseconds. We're
    // seeing incorrect adjusted realtime a* = r - r0*, but we know the actual
    // value is a. Pass in the error e=(a*-a), then we want to calculate a
    // new base adjustment r0 such that a = r - r0. So:
    //      a = r - r0* - e => r0=r0*+e.
    void readjustAdjustedRealTimeBy(const long long& e) {
        m_adjustedRealTimeBase += e; 
        ++numAdjustmentsToRealTimeBase;
    }

    long long getAdjustedRealTime() 
    {   return realTimeInNs() - m_adjustedRealTimeBase; }

    // Increment the reference count and return its new value.
    int incrRefCount() const {return ++m_refCount;}

    // Decrement the reference count and return its new value.
    int decrRefCount() const {return --m_refCount;}

    // Get the current value of the reference counter.
    int getRefCount() const {return m_refCount;}

    const MultibodySystem&                  m_system;
    VisualizerProtocol                      m_protocol;
    bool                                    m_shutdownWhenDestructed;

    Array_<DecorativeGeometry>              m_addedGeometry;
    Array_<RubberBandLine>                  m_lines;
    Array_<DecorationGenerator*>            m_generators;
    Array_<Visualizer::InputListener*>      m_listeners;
    Array_<Visualizer::FrameController*>    m_controllers;

    CoordinateDirection                     m_upDirection;
    Real                                    m_groundHeight;

    // User control of Visualizer behavior.
    Visualizer::Mode    m_mode;
    Real    m_frameRateFPS;       // in frames/sec if > 0, else use default
    Real    m_simTimeUnitsPerSec; // ratio of sim time units to real seconds
    Real    m_desiredBufferLengthInSec; // RT only: how much delay (<0 => default)

    // How many nanoseconds between frames?
    long long m_timeBetweenFramesInNs;
    // How much accuracy should we require from sleep()?
    long long m_allowableFrameJitterInNs;
    // How much slop is allowed in matching the time of a simulation frame
    // to the real time at which its frame is drawn?
    long long m_allowableFrameTimeSlopInNs;

    // The offset r0 to subtract from the interval timer reading to produce 
    // the adjusted real time a that we expect to match the current simulation
    // time t in ns. That is a = realTimeInNs()-r0. This base is adjusted by
    // the drawing thread when we see what time we actually were able to
    // deliver a frame.
    long long m_adjustedRealTimeBase; // r0
    
    // This is when we would like the simulation to send us another frame.
    // It is optimistically set to one frame interval later than the desired
    // draw time of the most recent frame to be put in the queue. This is 
    // also used in non-RealTime mode where AdjRT==RT.
    long long m_nextFrameDueAdjRT;

    // In RealTime mode we remember the simulated time in the previous
    // supplied frame so that we can tell if we see an earlier frame,
    // meaning (most likely) that we are starting a new simulation or
    // seeing a playback of an old one.
    double m_prevFrameSimTime;

    // The frame buffer:
    Array_<Frame,int> m_pool; // fixed size, old to new order but circular
    int m_oldest, m_nframe;   // oldest is index into pool, nframe is #valid entries
    std::mutex              m_queueMutex;
    std::condition_variable m_queueNotFull;  // these must use m_queueMutex
    std::condition_variable m_queueNotEmpty;
    std::condition_variable m_queueIsEmpty;

    std::thread         m_drawThread;    // the rendering thread
    bool                m_drawThreadIsRunning;
    bool                m_drawThreadShouldSuicide;

    mutable int         m_refCount; // how many Visualizer handles reference
                                    //   this Impl object?

    // Statistics
    int numFramesReportedBySimulation;
    int   numReportedFramesThatWereIgnored;
    int   numReportedFramesThatHadToWait;
    int   numReportedFramesThatSkippedAhead;
    int   numReportedFramesThatArrivedTooLate;
    int   numReportedFramesThatWereQueued;
    int     numQueuedFramesThatHadToWait;

    int numFramesSentToRenderer;
    int   numFramesDelayedByRenderer;
    int numTimesDrawThreadBlockedOnEmptyQueue;
    int numAdjustmentsToRealTimeBase;

    double sumOfAllJitter;        // updated at time sent to renderer (ms)
    double sumSquaredOfAllJitter; // ms^2

    // These are updated by the drawing thread each time it looks at the
    // queue to pull off a frame.
    double sumOfQueueLengths; // for computing the average length
    double sumSquaredOfQueueLengths; // for std deviation


    void clearStats() {
        numFramesReportedBySimulation=0;
          numReportedFramesThatWereIgnored=0;
          numReportedFramesThatHadToWait=0;
          numReportedFramesThatSkippedAhead=0;
          numReportedFramesThatArrivedTooLate=0;
          numReportedFramesThatWereQueued=0;
            numQueuedFramesThatHadToWait=0;

        numFramesSentToRenderer=0;
          numFramesDelayedByRenderer=0;
        numTimesDrawThreadBlockedOnEmptyQueue=0;
        numAdjustmentsToRealTimeBase=0;

        sumOfAllJitter        = 0;
        sumSquaredOfAllJitter = 0;

        sumOfQueueLengths = 0;
        sumSquaredOfQueueLengths = 0;
    }

    void dumpStats(std::ostream& o) const {
        o << "Visualizer stats:\n";
        o << "  Mode: "; 
        switch(m_mode) {
        case PassThrough: o << "PassThrough\n"; break;
        case Sampling: o << "Sampling\n"; break;
        case RealTime: 
            o << "RealTime, TimeScale=" << m_simTimeUnitsPerSec 
              << " sim time units/real second\n"; 
            o << "  Desired/actual buffer length(s): " 
              << getDesiredBufferLengthInSec() << "/" 
              << getActualBufferLengthInSec() << " (" 
              << getActualBufferLengthInFrames() << " frames)\n";
            break;
        };
        o << "  Desired frame rate: " << m_frameRateFPS << endl;
        o << "  reported frames: " << numFramesReportedBySimulation << endl;
        o << "  |       ignored: " << numReportedFramesThatWereIgnored << endl;
        o << "  |   had to wait: " << numReportedFramesThatHadToWait << endl;
        o << "  | skipped ahead: " << numReportedFramesThatSkippedAhead << endl;
        o << "  | came too late: " << numReportedFramesThatArrivedTooLate << endl;
        o << "  | were buffered: " << numReportedFramesThatWereQueued << endl;
        o << "  | | full buffer: " << numQueuedFramesThatHadToWait << endl;
        o << "  frames sent to renderer: " << numFramesSentToRenderer << endl;
        o << "  | delayed by renderer  : " << numFramesDelayedByRenderer << endl;
        if (numReportedFramesThatWereQueued && numFramesSentToRenderer) {
            const double avg = sumOfQueueLengths/numFramesSentToRenderer;
            o << "  | average buffer length (frames): " << avg << endl;
            o << "  | std dev buffer length (frames): " 
              << std::sqrt(std::max(0.,
                              sumSquaredOfQueueLengths/numFramesSentToRenderer 
                              - square(avg))) << endl;
        }
        o << "  draw blocked for empty buffer: " 
          << numTimesDrawThreadBlockedOnEmptyQueue << endl;
        o << "  adjustments to real time base: " 
          << numAdjustmentsToRealTimeBase << endl;
        if (numFramesSentToRenderer > 0) {
            const double avg = sumOfAllJitter/numFramesSentToRenderer;
            o << "  average jitter (ms): " << avg << endl;
            o << "  jitter std dev (ms): " 
              << std::sqrt(sumSquaredOfAllJitter/numFramesSentToRenderer
                 - square(avg)) << endl;
        }
    }

};

// Generate geometry for the given state and send it to the visualizer using
// the VisualizerProtocol object. In buffered mode this is called from the
// rendering thread; otherwise, this is just the main simulation thread.
void Visualizer::Impl::drawFrameNow(const State& state) {
    m_system.realize(state, Stage::Position);

    // Collect up the geometry that constitutes this scene.
    Array_<DecorativeGeometry> geometry;
    for (Stage stage = Stage::Topology; stage <= state.getSystemStage(); ++stage)
        m_system.calcDecorativeGeometryAndAppend(state, stage, geometry);
    for (unsigned i = 0; i < m_generators.size(); i++)
        m_generators[i]->generateDecorations(state, geometry);

    // Execute frame controls (e.g. camera positioning).
    for (unsigned i = 0; i < m_controllers.size(); ++i)
        m_controllers[i]->generateControls(Visualizer(this), state, geometry);

    // Calculate the spatial pose of all the geometry and send it to the
    // renderer.
    m_protocol.beginScene(state.getTime());
    VisualizerGeometry geometryCreator
        (m_protocol, m_system.getMatterSubsystem(), state);
    for (unsigned i = 0; i < geometry.size(); ++i)
        geometry[i].implementGeometry(geometryCreator);
    for (unsigned i = 0; i < m_addedGeometry.size(); ++i)
        m_addedGeometry[i].implementGeometry(geometryCreator);
    const SimbodyMatterSubsystem& matter = m_system.getMatterSubsystem();
    for (unsigned i = 0; i < m_lines.size(); ++i) {
        const RubberBandLine& line = m_lines[i];
        const MobilizedBody& B1 = matter.getMobilizedBody(line.b1);
        const MobilizedBody& B2 = matter.getMobilizedBody(line.b2);
        const Transform&  X_GB1 = B1.getBodyTransform(state);
        const Transform&  X_GB2 = B2.getBodyTransform(state);
        const Vec3 end1 = X_GB1*line.station1;
        const Vec3 end2 = X_GB2*line.station2;
        const Real thickness = line.line.getLineThickness() == -1 
                               ? 1 : line.line.getLineThickness();
        m_protocol.drawLine(end1, end2, 
            VisualizerGeometry::getColor(line.line), thickness);
    }
    m_protocol.finishScene();

    ++numFramesSentToRenderer;
}

// This is called from the drawing thread if we're buffering, otherwise
// directly from the simulation thread.
void Visualizer::Impl::drawRealtimeFrameWhenReady
   (const State& state, const long long& desiredDrawTimeAdjRT)
{
    const long long earliestDrawTimeAdjRT = 
        desiredDrawTimeAdjRT - m_allowableFrameJitterInNs;
    const long long latestDrawTimeAdjRT = 
        desiredDrawTimeAdjRT + m_allowableFrameJitterInNs;

    // Wait for the next frame time, allowing for a little jitter 
    // since we can't expect sleep to wake us up at the exact time.
    long long now = getAdjustedRealTime();
    if (now < earliestDrawTimeAdjRT) {
        ++numFramesDelayedByRenderer;
        do {sleepInNs(desiredDrawTimeAdjRT-now);}
        while ((now=getAdjustedRealTime()) < earliestDrawTimeAdjRT);

        // Keep stats on the jitter situation.
        const long long jitterInNs = now - desiredDrawTimeAdjRT;
        const double jitterInMs = jitterInNs*1e-6;
        sumOfAllJitter        += jitterInMs;
        sumSquaredOfAllJitter += square(jitterInMs);
    }

    // timingError is signed with + meaning we sent the frame late.
    const long long timingError = now - desiredDrawTimeAdjRT;

    // If we sent this frame more than one frame time late we're going to 
    // admit we're not making real time and adjust the
    // AdjRT base to  match.
    if (timingError > m_timeBetweenFramesInNs)
        readjustAdjustedRealTimeBy(now - desiredDrawTimeAdjRT);
   
    // It is time to render the frame.
    drawFrameNow(state);   
}

// Attempt to report a frame while we're in realtime mode. 
void Visualizer::Impl::reportRealtime(const State& state) {
    // If we see a simulation time that is earlier than the last one, 
    // we are probably starting a new simulation or playback. Flush the
    // old one. Note that we're using actual simulation time; we don't
    // want to get tricked by adjustments to the real time base.
    if (state.getTime() < m_prevFrameSimTime) {
        waitUntilQueueIsEmpty();
        m_nextFrameDueAdjRT = -1; // restart time base
    }

    // scale, convert to ns (doesn't depend on real time base)
    const long long t = convertSimTimeToAdjRT(state.getTime()); 

    // If this is the first frame, or first since last setMode(), then
    // we synchronize Adjusted Real Time to match. Readjustments will occur
    // if the simulation doesn't keep up with real time.
    if (m_nextFrameDueAdjRT < 0 || t == 0) {
        m_adjustedRealTimeBase = realTimeInNs() - t;
        // now getAdjustedRealTime()==t
        m_nextFrameDueAdjRT = t; // i.e., now
    }

    // "timeSlop" is the amount we'll allow a frame's simulation time to 
    // deviate from the real time at which we draw it. That is, if we're 
    // expecting a frame at time t_f and the simulator instead delivers a
    // frame at t_s, we'll consider that a match if |t_s-t_f|<=slop.
    // The reason for this is that we prefer to issue frames at regular 
    // intervals, so if the frame time and sim time match closely enough
    // we won't reschedule the frames. Otherwise, a sim frame whose time
    // is too early (t_s<t_f-slop) gets thrown away (or used in desperation
    // if we're not keeping up with real time), and a sim frame
    // whose time is too late (t_s>t_f+slop) causes us to delay drawing
    // that frame until real time catches up with what's in it. Typically 
    // timeSlop is set to a small fraction of the frame time, like 5%.
    const long long timeSlop = m_allowableFrameTimeSlopInNs;
    const long long next     = m_nextFrameDueAdjRT;
    const long long earliest = next - timeSlop;
    const long long latest   = next + timeSlop;

    if (t < earliest) {
        ++numReportedFramesThatWereIgnored; // we don't need this one
        return;
    } 
    
    long long desiredDrawTimeAdjRT = next;
    if (t > latest) {
        ++numReportedFramesThatSkippedAhead;
        desiredDrawTimeAdjRT = t;
    }

    // If buffering is enabled, push this onto the queue. Note that we
    // might have to wait until the queue has some room.
    if (queuingIsEnabled()) {
        // This also sets expectations for the next frame.
        addFrameToQueueWithWait(state, desiredDrawTimeAdjRT);
        return;
    }

    // There is no buffer. We'll just render this frame as soon as its
    // drawing time arrives. No need to copy the state here. Note that 
    // the simulation thread is doing the drawing as well as the simulating.
    // This method will also readjust adjusted real time if the frame came
    // too late.
    drawRealtimeFrameWhenReady(state, desiredDrawTimeAdjRT);

    // Now set expectations for the next frame.
    m_nextFrameDueAdjRT = t + m_timeBetweenFramesInNs;
}



//==============================================================================
//                                VISUALIZER
//==============================================================================
Visualizer::Visualizer(Visualizer::Impl* srcImpl) : impl(srcImpl) {
    if (impl) impl->incrRefCount();
}

Visualizer::Visualizer(const MultibodySystem& system) : impl(0) {
    impl = new Impl(this, system, Array_<String>());
    impl->incrRefCount();
}

Visualizer::Visualizer(const MultibodySystem& system,
                       const Array_<String>& searchPath) : impl(0) {
    impl = new Impl(this, system, searchPath);
    impl->incrRefCount();
}

Visualizer::Visualizer(const Visualizer& source) : impl(0) {
    if (source.impl) {
        impl = source.impl;
        impl->incrRefCount();
    }
}

Visualizer& Visualizer::operator=(const Visualizer& source) {
    if (impl != source.impl) {
        if (impl&& impl->decrRefCount()==0) delete impl;
        impl = source.impl;
        impl->incrRefCount();
    }
    return *this;
}

Visualizer::~Visualizer() {
    if (impl && impl->decrRefCount()==0)
        delete impl;
}

void Visualizer::shutdown() 
{   updImpl().m_protocol.shutdownGUI(); }

Visualizer& Visualizer::setShutdownWhenDestructed(bool shouldShutdown)
{   updImpl().setShutdownWhenDestructed(shouldShutdown); return *this; }

bool Visualizer::getShutdownWhenDestructed() const
{   return getImpl().getShutdownWhenDestructed(); }

int Visualizer::getRefCount() const
{   return impl ? impl->getRefCount() : 0; }

       // Frame drawing methods

void Visualizer::drawFrameNow(const State& state) const
{   const_cast<Visualizer*>(this)->updImpl().drawFrameNow(state); }

void Visualizer::flushFrames() const
{   const_cast<Visualizer*>(this)->updImpl().waitUntilQueueIsEmpty(); }

// The simulation thread normally delivers frames here. Handling is dispatched
// according the current visualization mode.
void Visualizer::report(const State& state) const {
    Visualizer::Impl& rep = const_cast<Visualizer*>(this)->updImpl();

    ++rep.numFramesReportedBySimulation;
    if (rep.m_mode == RealTime) {
        rep.reportRealtime(state);
        return;
    }

    // We're in Sampling or PassThrough mode. AdjRT and RT are the same in
    // these modes; they are just real time as determined by realTimeInNs(),
    // the current value of the interval counter. We don't care at all what
    // time the simulation thinks it is.

    // If this is the first frame, or first since last setMode(), then
    // we set our expected next frame arrival time to now.
    if (rep.m_nextFrameDueAdjRT < 0LL)
        rep.m_nextFrameDueAdjRT = realTimeInNs(); // now

    // If someone asked for an infinite frame rate just send this along now.
    if (rep.m_timeBetweenFramesInNs == 0LL) {
        drawFrameNow(state);
        return;
    }

    const long long earliestDrawTime = rep.m_nextFrameDueAdjRT
                                       - rep.m_allowableFrameJitterInNs;
    long long now = realTimeInNs();
    if (now < earliestDrawTime) {
        // Too early to draw this frame. In Sampling mode that means we
        // just ignore it.
        if (rep.m_mode == Sampling) {
            ++rep.numReportedFramesThatWereIgnored;
            return;
        }

        // We're in PassThrough mode.
        ++rep.numReportedFramesThatHadToWait;
        do {sleepInNs(rep.m_nextFrameDueAdjRT - now);}
        while ((now=realTimeInNs()) < earliestDrawTime);

        // We're not going to wake up exactly when we wanted to; keep stats.
        const double jitterInMs = (now - rep.m_nextFrameDueAdjRT)*1e-6;
        rep.sumOfAllJitter        += jitterInMs;
        rep.sumSquaredOfAllJitter += square(jitterInMs);
    }

    // Frame time reached in Sampling or PassThrough modes. Draw the frame.
    drawFrameNow(state);

    // This frame might have been on time or late; we'll schedule the next 
    // time for one ideal frame interval later to keep the maximum rate down 
    // to the specified rate. Otherwise a late frame could be followed by lots
    // of fast frames playing catch-up.
    if (now-rep.m_nextFrameDueAdjRT <= rep.m_timeBetweenFramesInNs)
        rep.m_nextFrameDueAdjRT += rep.m_timeBetweenFramesInNs;
    else { // a late frame; delay the next one
        rep.m_nextFrameDueAdjRT = now + rep.m_timeBetweenFramesInNs;
        ++rep.numAdjustmentsToRealTimeBase;
    }
}

        // Visualizer display options

Visualizer& Visualizer::setBackgroundType(BackgroundType type) 
{   updImpl().m_protocol.setBackgroundType(type); return *this; }

const Visualizer& Visualizer::setBackgroundColor(const Vec3& color) const 
{   getImpl().m_protocol.setBackgroundColor(color); return *this; }

const Visualizer& Visualizer::setShowShadows(bool showShadows) const 
{   getImpl().m_protocol.setShowShadows(showShadows); return *this; }

const Visualizer& Visualizer::setShowFrameRate(bool showFrameRate) const 
{   getImpl().m_protocol.setShowFrameRate(showFrameRate); return *this; }

const Visualizer& Visualizer::setShowSimTime(bool showSimTime) const 
{   getImpl().m_protocol.setShowSimTime(showSimTime); return *this; }

const Visualizer& Visualizer::setShowFrameNumber(bool showFrameNumber) const 
{   getImpl().m_protocol.setShowFrameNumber(showFrameNumber); return *this; }

const Visualizer& Visualizer::setWindowTitle(const String& title) const 
{   getImpl().m_protocol.setWindowTitle(title); return *this; }

        // Visualizer options

Visualizer& Visualizer::setSystemUpDirection(const CoordinateDirection& upDir)
{   updImpl().m_upDirection = upDir; 
    updImpl().m_protocol.setSystemUpDirection(upDir); return *this; }
CoordinateDirection Visualizer::getSystemUpDirection() const
{   return getImpl().m_upDirection; }


Visualizer& Visualizer::setGroundHeight(Real height) {
    updImpl().m_groundHeight = height;
    updImpl().m_protocol.setGroundHeight(height); return *this;
}
Real Visualizer::getGroundHeight() const
{   return getImpl().m_groundHeight; }

Visualizer& Visualizer::setMode(Visualizer::Mode mode) 
{   updImpl().setMode(mode); return *this; }
Visualizer::Mode Visualizer::getMode() const {return getImpl().m_mode;}

Visualizer& Visualizer::setDesiredFrameRate(Real fps) 
{   updImpl().setDesiredFrameRate(std::max(fps, Real(0))); return *this; }
Real Visualizer::getDesiredFrameRate() const 
{   return getImpl().m_frameRateFPS; }

Visualizer& Visualizer::setRealTimeScale(Real simTimePerRealSec) 
{   updImpl().setRealTimeScale(simTimePerRealSec); return *this; }
Real Visualizer::getRealTimeScale() const 
{   return getImpl().m_simTimeUnitsPerSec; }

Visualizer& Visualizer::setDesiredBufferLengthInSec(Real bufferLengthInSec)
{   updImpl().setDesiredBufferLengthInSec(bufferLengthInSec); return *this; }
Real Visualizer::getDesiredBufferLengthInSec() const
{   return getImpl().getDesiredBufferLengthInSec(); }
int Visualizer::getActualBufferLengthInFrames() const 
{   return getImpl().getActualBufferLengthInFrames(); }
Real Visualizer::getActualBufferLengthInSec() const 
{   return getImpl().getActualBufferLengthInSec(); }


int Visualizer::addInputListener(Visualizer::InputListener* listener) {
    Impl& impl = updImpl();
    const int nxt = (int)impl.m_listeners.size();
    impl.m_listeners.push_back(listener); 
    return nxt; 
}
int Visualizer::getNumInputListeners() const 
{   return (int)getImpl().m_listeners.size(); }
const Visualizer::InputListener& Visualizer::getInputListener(int i) const 
{   return *getImpl().m_listeners[i]; }
Visualizer::InputListener& Visualizer::updInputListener(int i) 
{   return *updImpl().m_listeners[i]; }

int Visualizer::addFrameController(Visualizer::FrameController* fc) {
    Impl& impl = updImpl();
    const int nxt = (int)impl.m_controllers.size();
    impl.m_controllers.push_back(fc); 
    return nxt; 
}
int Visualizer::getNumFrameControllers() const 
{   return (int)getImpl().m_controllers.size(); }
const Visualizer::FrameController& Visualizer::getFrameController(int i) const 
{   return *getImpl().m_controllers[i]; }
Visualizer::FrameController& Visualizer::updFrameController(int i) 
{   return *updImpl().m_controllers[i]; }


        // Scene-building methods

Visualizer& Visualizer::
addMenu(const String& title, int menuId, 
        const Array_<pair<String, int> >& items) 
{
    SimTK_ERRCHK2_ALWAYS(menuId >= 0, "Visualizer::addMenu()",
        "Assigned menu ids must be nonnegative, but an attempt was made to create"
        " a menu %s with id %d.", title.c_str(), menuId);

    updImpl().m_protocol.addMenu(title, menuId, items);
    return *this;
}

Visualizer& Visualizer::
addSlider(const String& title, int sliderId, 
          Real minVal, Real maxVal, Real value) 
{
    SimTK_ERRCHK2_ALWAYS(sliderId >= 0, "Visualizer::addSlider()",
        "Assigned slider ids must be nonnegative, but an attempt was made to create"
        " a slider %s with id %d.", title.c_str(), sliderId);
    SimTK_ERRCHK4_ALWAYS(minVal <= value && value <= maxVal, "Visualizer::addSlider()", 
        "Initial slider value %g for slider %s was outside the specified range [%g,%g].",
        value, title.c_str(), minVal, maxVal);

    updImpl().m_protocol.addSlider(title, sliderId, minVal, maxVal, value);
    return *this;
}

int Visualizer::
addDecoration(MobilizedBodyIndex mobodIx, const Transform& X_BD, 
              const DecorativeGeometry& geom) 
{
    Array_<DecorativeGeometry>& addedGeometry = updImpl().m_addedGeometry;
    const int nxt = (int)addedGeometry.size();
    addedGeometry.push_back(geom);
    DecorativeGeometry& geomCopy = addedGeometry.back();
    geomCopy.setBodyId((int)mobodIx);
    geomCopy.setTransform(X_BD * geomCopy.getTransform());
    return nxt;
}
int Visualizer::getNumDecorations() const 
{   return (int)getImpl().m_addedGeometry.size(); }
const DecorativeGeometry& Visualizer::getDecoration(int i) const 
{   return getImpl().m_addedGeometry[i]; }
DecorativeGeometry& Visualizer::updDecoration(int i) const
{   return const_cast<Visualizer*>(this)->updImpl().m_addedGeometry[i]; }

int Visualizer::
addRubberBandLine(MobilizedBodyIndex b1, const Vec3& station1, 
                  MobilizedBodyIndex b2, const Vec3& station2, 
                  const DecorativeLine& line) 
{   
    Impl& impl = updImpl();
    const int nxt = (int)impl.m_lines.size();
    impl.m_lines.push_back(RubberBandLine(b1,station1, b2,station2, line)); 
    return nxt; 
}
int Visualizer::getNumRubberBandLines() const 
{   return (int)getImpl().m_lines.size(); }
const DecorativeLine& Visualizer::getRubberBandLine(int i) const 
{   return getImpl().m_lines[i].line; }
DecorativeLine& Visualizer::updRubberBandLine(int i) const
{   return const_cast<Visualizer*>(this)->updImpl().m_lines[i].line; }

int Visualizer::
addDecorationGenerator(DecorationGenerator* generator) 
{   
    Impl& impl = updImpl();
    const int nxt = (int)impl.m_generators.size();
    impl.m_generators.push_back(generator); 
    return nxt;
}
int Visualizer::
getNumDecorationGenerators() const 
{   return (int)getImpl().m_generators.size(); }
const DecorationGenerator& Visualizer::
getDecorationGenerator(int i) const 
{   return *getImpl().m_generators[i]; }
DecorationGenerator& Visualizer::
updDecorationGenerator(int i) 
{   return *updImpl().m_generators[i]; }

        // Frame control methods
const Visualizer& Visualizer::
setCameraTransform(const Transform& transform) const 
{   getImpl().m_protocol.setCameraTransform(transform); return *this; }

const Visualizer& Visualizer::zoomCameraToShowAllGeometry() const 
{   getImpl().m_protocol.zoomCamera(); return *this; }

const Visualizer& Visualizer::
pointCameraAt(const Vec3& point, const Vec3& upDirection) const 
{   getImpl().m_protocol.lookAt(point, upDirection); return *this; }

const Visualizer& Visualizer::setCameraFieldOfView(Real fov) const 
{   getImpl().m_protocol.setFieldOfView(fov); return *this; }

const Visualizer& Visualizer::
setCameraClippingPlanes(Real nearPlane, Real farPlane) const 
{   getImpl().m_protocol.setClippingPlanes(nearPlane, farPlane);
    return *this; }


const Visualizer& Visualizer::setSliderValue(int slider, Real newValue) const 
{   getImpl().m_protocol.setSliderValue(slider, newValue); return *this; }

const Visualizer& Visualizer::
setSliderRange(int slider, Real newMin, Real newMax) const 
{   getImpl().m_protocol.setSliderRange(slider, newMin, newMax); return *this; }

        // Debugging and statistics
void Visualizer::dumpStats(std::ostream& o) const {getImpl().dumpStats(o);}
void Visualizer::clearStats() {updImpl().clearStats();}

        // Internal use only
const Array_<Visualizer::InputListener*>& Visualizer::getInputListeners() const
{   return getImpl().m_listeners; }
const Array_<Visualizer::FrameController*>& Visualizer::getFrameControllers() const
{   return getImpl().m_controllers; }
const MultibodySystem& Visualizer::getSystem() const {return getImpl().m_system;}


//==============================================================================
//                             BODY FOLLOWER
//==============================================================================
Visualizer::BodyFollower::BodyFollower(
        const MobilizedBody& mobodB,
        const Vec3&          stationPinB,
        const Vec3&          offset,
        const UnitVec3&      upDirection)
    :   m_mobodB(mobodB), m_stationPinB(stationPinB), m_offset(offset),
        m_upDirection(upDirection) {}
    
void Visualizer::BodyFollower::generateControls(
        const Visualizer&             viz,
        const State&                  state,
        Array_< DecorativeGeometry >& geometry)
{
    // Offset.
    Vec3 offset(m_offset);
    if (m_offset.isNaN()) {
        // Default: offset is based on system up direction and ground height.
        offset = Vec3(1, 1, 1);
        offset[viz.getSystemUpDirection().getAxis()] += viz.getGroundHeight();
    }

    // Up direction. Default: use System up direction.
    const UnitVec3& upDirection = m_upDirection.isNaN() ?
        UnitVec3(viz.getSystemUpDirection()) : m_upDirection;

    const Vec3 P = m_mobodB.findStationLocationInGround(state, m_stationPinB);
    // Position of camera (C) from ground origin (G), expressed in ground.
    const Vec3 p_GC = P + offset;
    // Rotation of camera frame (C) in ground frame (G).
    // To get the camera to point at P, we require the camera's z direction
    // (which points "back") to be parallel to the offset. We also want the
    // camera's y direction (which points to the top of the screen) to be as
    // closely aligned with the provided up direction as is possible.
    const Rotation R_GC(UnitVec3(offset), ZAxis, upDirection, YAxis);
    viz.setCameraTransform(Transform(R_GC, p_GC));
}


//==============================================================================
//                             THE DRAWING THREAD
//==============================================================================
/* When we're in RealTime mode, we run a separate thread that actually sends
frames to the renderer. It pulls frames off the back (oldest) end of the 
frame buffer queue, while the simulation thread is pushing frames onto the 
front (newest) end of the queue. The rendering thread is created whenever 
the Visualizer enters RealTime mode, and canceled whenever it leaves 
RealTime mode or is destructed.

This is the main function for the buffered drawing thread. Its job is to 
pull the oldest frames off the queue and send them to the renderer at
the right real times. We use adjusted real time (AdjRT) which should match
the simulation time kept with the frame as its desired draw time. If the
frame is ahead of AdjRT, we'll sleep to let real time catch up. If the 
frame is substantially behind, we'll render it now and then adjust the AdjRT 
base to acknowledge that we have irretrievably slipped from real time and need 
to adjust our expectations for the future. */
static void drawingThreadMain(Visualizer::Impl& vizImpl) {

    do {
        // Grab the oldest frame in the queue.
        // This will wait if necessary until a frame is available.
        const Frame* framep;
        if (vizImpl.getOldestFrameInQueue(&framep)) {
            // Draw this frame as soon as its draw time arrives, and readjust
            // adjusted real time if necessary.
            vizImpl.drawRealtimeFrameWhenReady
               (framep->state, framep->desiredDrawTimeAdjRT); 

            // Return the now-rendered frame to circulation in the pool. This may
            // wake up the simulation thread if it was waiting for space.
            vizImpl.noteThatOldestFrameIsNowAvailable();
        }
    } while (!vizImpl.m_drawThreadShouldSuicide);

    // Attempt to wake up the simulation thread if it is waiting for
    // the draw thread since there won't be any more notices!
    vizImpl.m_queueNotFull.notify_one(); // wake up if waiting for queue space
    vizImpl.m_queueIsEmpty.notify_one(); // wake up if flushing
}
