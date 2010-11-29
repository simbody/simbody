/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
 * Authors: Peter Eastman, Michael Sherman                                    *
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
#include "simbody/internal/Visualizer_EventListener.h"
#include "simbody/internal/DecorationGenerator.h"
#include "VisualizationGeometry.h"
#include "VisualizationProtocol.h"

#include <cstdlib>
#include <cstdio>
#include <pthread.h>
#include <string>
#include <ctime>
#include <iostream>
#include <limits>

using namespace SimTK;
using namespace std;


static void* drawingThreadMain(void* visualizerAsVoidp);

static const long long UsToNs = 1000LL;          // ns = us * UsToNs
static const long long MsToNs = 1000LL * UsToNs; // ns = ms * MsToNs

static const Real      DefaultFrameRateFPS             = 30;
static const Real      DefaultDesiredBufferLengthInSec = Real(0.15); // 150ms

// These are not currently overrideable.
static const long long DefaultAllowableFrameJitterInNs      = 5 * MsToNs; //5ms
static const Real      DefaultSlopAsFractionOfFrameInterval = Real(0.05); //5%

//==============================================================================
//                              VISUALIZER REP
//==============================================================================
/* This is the private implementation object contained in a Visualizer handle.
RealTime mode is worth some discussion. There is a simulation thread that
produces frames at a variable rate, and a draw thread that consumes frames at a
variable rate (by sending them to the renderer). We want to engineer things so 
that frames are sent to the renderer at a steady rate that is synchronized with
simulation time (possibly after scaling). When a thread is running too fast, 
that is easily handled by blocking the speeding thread for a while. The "too 
slow" case takes careful handling.

In normal operation, we expect the simulation to take varying amounts of
real time to generate fixed amounts of simulation time, because we prefer
to use variable time-step integrators that control errors by taking smaller
steps in more difficult circumstances, and large steps through the easy
parts of the simulation. For real time operation, the simulation must of
course *average* real time performance; we use a frame buffer to smooth
out variable delivery times. That is, frames go into the buffer at an
irregular rate but are pulled off at a regular rate. A longer buffer can
mask wider deviations in frame time, at the expense of interactive response.
In most circumstances people cannot perceive delays below about 200ms, so
for good response the total delay should be kept around that level.

Despite the buffering, there will be occasions when the simulation can't
keep up with real time. A common cause of that is that a user has paused
either the simulation or the renderer, such as by hitting a breakpoint while
debugging. In that case we deem the proper behavior to be that when we 
resume we should immediately resume real time behavior at a new start time, 
*not* attempt to catch up to the previous real time by running at high speed. 
As much as possible, we would like the simulation to behave just as it would 
have without the interruption, but with a long pause where interrupted. We
deal with this situation by introducing a notion of "adjusted real time"
(AdjRT). That is a clock that tracks the real time interval counter, but uses
a variable base offset that is used to match it to the expected simulation 
time. When the simulation is long delayed, we modify the AdjRT base when we
resume so that AdjRT once again matches the simulation time t. Adjustments
to the AdjRT base occur at the time we deliver frames to the renderer; at that
moment we compare the AdjRT reading to the frame's simulation time t and 
correct AdjRT for future frames.

You can also run in RealTime mode without buffering. In that case frames are
sent directly from the simulation thread to the renderer, but the above logic
still applies. Simulation frames that arrive earlier than the corresponding
AdjRT are delayed; frames that arrive later are drawn immediately but cause
AdjRT to be readjusted to resynchronize. Overall performance can be better
in unbuffered RealTime mode because the States provided by the simulation do
not have to be copied before being drawn. However, intermittent slower-than-
realtime frame times cannot be smoothed over; they will cause rendering delays.

PassThrough and Sampling modes are much simpler because no synchronization
is done to the simulation times. There is only a single thread and draw
time scheduling works in real time without adjustment.
*/

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

// Implementation of the Visualizer.
class Visualizer::VisualizerRep {
public:
    // Create a Visualizer and put it in PassThrough mode.
    VisualizerRep(Visualizer* owner, MultibodySystem& system, 
                  const String& title) 
    :   handle(owner), system(system), protocol(*owner, title),
        m_mode(PassThrough), m_frameRateFPS(DefaultFrameRateFPS), 
        m_simTimeUnitsPerSec(1), 
        m_desiredBufferLengthInSec(DefaultDesiredBufferLengthInSec), 
        m_timeBetweenFramesInNs(secToNs(1/DefaultFrameRateFPS)),
        m_allowableFrameJitterInNs(DefaultAllowableFrameJitterInNs),
        m_allowableFrameTimeSlopInNs(
            secToNs(DefaultSlopAsFractionOfFrameInterval/DefaultFrameRateFPS)),
        m_adjustedRealTimeBase(realTimeInNs()),
        m_nextFrameDueAdjRT(-1), oldest(0),nframe(0)
    {   
        pthread_mutex_init(&queueLock, NULL); 
        pthread_cond_init(&queueNotFull, NULL); 
        pthread_cond_init(&queueNotEmpty, NULL); 

        setMode(PassThrough);
        clearStats();
    }
    
    ~VisualizerRep() {
        if (m_mode==RealTime && pool.size()) {
            pthread_cancel(drawThread);
        }
        for (unsigned i = 0; i < controllers.size(); i++)
            delete controllers[i];
        for (unsigned i = 0; i < listeners.size(); i++)
            delete listeners[i];
        for (unsigned i = 0; i < generators.size(); i++)
            delete generators[i];
        pthread_cond_destroy(&queueNotEmpty);
        pthread_cond_destroy(&queueNotFull);
        pthread_mutex_destroy(&queueLock);
    }

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
        if (m_mode == RealTime && numFrames != pool.size()) {
            if (pool.size()) {
                // draw thread is already running
                if (numFrames == 0) {
                    pthread_cancel(drawThread);
                    sleepInSec(getActualBufferLengthInSec());
                }
                initializePool(numFrames);
            } else {
                // draw thread wasn't running but now should be
                initializePool(numFrames);
                pthread_create(&drawThread, NULL, drawingThreadMain, this);
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
        if (m_mode == RealTime && pool.size()) {
            pthread_cancel(drawThread);
            sleepInSec(getActualBufferLengthInSec());
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
    }

    // TODO: must recalculate buffer length
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

    int getActualBufferLengthInFrames() const {return pool.size();}
    Real getActualBufferLengthInSec() const 
    {   return nsToSec(getActualBufferLengthInFrames()*m_timeBetweenFramesInNs); }

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
        pthread_mutex_lock(&queueLock);
        pool.resize(sz);oldest=nframe=0;
        pthread_mutex_unlock(&queueLock);
    }

    int getNFramesInQueue() const {return nframe;}

    // Queing is enabled if the pool was allocated.
    bool queuingIsEnabled() const {return pool.size() != 0;}
    bool queueIsFull() const {return nframe==pool.size();}
    bool queueIsEmpty() const {return nframe==0;}

    // I'm capitalizing these methods because they are VERY important!
    void LOCK_Queue()   {pthread_mutex_lock(&queueLock);}
    void UNLOCK_Queue() {pthread_mutex_unlock(&queueLock);}
    void WAIT_QueueNotFull() {pthread_cond_wait(&queueNotFull,&queueLock);}
    void POST_QueueNotFull() {pthread_cond_signal(&queueNotFull);}
    void WAIT_QueueNotEmpty() {pthread_cond_wait(&queueNotEmpty,&queueLock);}
    void POST_QueueNotEmpty() {pthread_cond_signal(&queueNotEmpty);}

    // Called from simulation thread. Blocks until there is room in
    // the queue, then inserts this state unconditionally, with the indicated
    // desired rendering time in adjusted real time. We then update the 
    // "time of next queue slot" to be one ideal frame interval later than 
    // the desired draw time.
    void addFrameToQueueWithWait(const State& state, 
                                 const long long& desiredDrawTimeAdjRT)
    {
        LOCK_Queue();
        ++numReportedFramesThatWereQueued;
        if (queueIsFull()) {
            ++numQueuedFramesThatHadToWait;
            do {WAIT_QueueNotFull();} // atomic: unlock, long wait, relock
            while (queueIsFull()); // must recheck condition
        }

        // There is room in the queue now. We're holding the lock.
        Frame& frame = pool[(oldest+nframe)%pool.size()];
        frame.state  = state;
        frame.desiredDrawTimeAdjRT = desiredDrawTimeAdjRT;

        // Set the expected next frame time.
        m_nextFrameDueAdjRT = desiredDrawTimeAdjRT + m_timeBetweenFramesInNs;

        if (++nframe == 1) 
            POST_QueueNotEmpty(); // wake up rendering thread on first frame

        UNLOCK_Queue();
    }

    // The drawing thread uses this to find the oldest frame in the buffer.
    // It may then at its leisure use the contained State to generate a screen
    // image. There is no danger of the simulation thread modifying this
    // frame; once it has been put in it stays there until the drawing thread
    // takes it out. When done it should return the frame to the pool.
    const Frame& getOldestFrameInQueue() {
        LOCK_Queue();
        if (nframe == 0) {
            ++numTimesDrawThreadBlockedOnEmptyQueue;
            do {WAIT_QueueNotEmpty();} // atomic: unlock, long wait, relock
            while (nframe==0); // must recheck condition
        } else {
            sumOfQueueLengths        += double(nframe);
            sumSquaredOfQueueLengths += double(square(nframe));
        }
        // There is at least one frame available now. We are holding the lock.
        UNLOCK_Queue();
        return pool[oldest]; // simulation thread doesn't change oldest
    }

    // Drawing thread uses this to note that it is done with the oldest
    // frame which may now be reused by the simulation thread. The queueNotFull
    // condition is posted if there is a reasonable amount of room in the pool 
    // now.
    void noteThatOldestFrameIsNowAvailable() {
        LOCK_Queue();
        oldest = (oldest+1)%pool.size(); // move to next-oldest
        --nframe; // there is now one fewer frame in use
        // Start the simulation again when the pool is about half empty.
        if (nframe <= pool.size()/2+1)
            POST_QueueNotFull();
        UNLOCK_Queue();
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

    Visualizer*                             handle;
    MultibodySystem&                        system;
    VisualizationProtocol                   protocol;

    Array_<DecorativeGeometry>              addedGeometry;
    Array_<RubberBandLine>                  lines;
    Array_<DecorationGenerator*>            generators;
    Array_<Visualizer::EventListener*>      listeners;
    Array_<Visualizer::FrameController*>    controllers;

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

    // The frame buffer:
    Array_<Frame,int> pool; // fixed size, old to new order but circular
    int oldest, nframe;     // oldest is index into pool, nframe is #valid entries
    pthread_mutex_t     queueLock;
    pthread_cond_t      queueNotFull;   // must use with queueLock
    pthread_cond_t      queueNotEmpty;  // must use with queueLock
    pthread_t           drawThread;     // the rendering thread


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
// the VisualizationProtocol object. In buffered mode this is called from the
// rendering thread; otherwise, this is just the main simulation thread.
void Visualizer::VisualizerRep::drawFrameNow(const State& state) {
    system.realize(state, Stage::Position);

    // Collect up the geometry that constitutes this scene.
    Array_<DecorativeGeometry> geometry;
    for (Stage stage = Stage::Topology; stage <= state.getSystemStage(); ++stage)
        system.calcDecorativeGeometryAndAppend(state, stage, geometry);
    for (unsigned i = 0; i < generators.size(); i++)
        generators[i]->generateDecorations(state, geometry);

    // Execute frame controls (e.g. camera positioning).
    for (unsigned i = 0; i < controllers.size(); ++i)
        controllers[i]->generateControls(*handle, state, geometry);

    // Calculate the spatial pose of all the geometry and send it to the
    // renderer.
    protocol.beginScene();
    VisualizationGeometry geometryCreator
        (protocol, system.getMatterSubsystem(), state);
    for (unsigned i = 0; i < geometry.size(); ++i)
        geometry[i].implementGeometry(geometryCreator);
    for (unsigned i = 0; i < addedGeometry.size(); ++i)
        addedGeometry[i].implementGeometry(geometryCreator);
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
    for (unsigned i = 0; i < lines.size(); ++i) {
        const RubberBandLine& line = lines[i];
        const MobilizedBody& B1 = matter.getMobilizedBody(line.b1);
        const MobilizedBody& B2 = matter.getMobilizedBody(line.b2);
        const Transform&  X_GB1 = B1.getBodyTransform(state);
        const Transform&  X_GB2 = B2.getBodyTransform(state);
        const Vec3 end1 = X_GB1*line.station1;
        const Vec3 end2 = X_GB2*line.station2;
        const Real thickness = line.line.getLineThickness() == -1 
                               ? 1 : line.line.getLineThickness();
        protocol.drawLine(end1, end2, 
            VisualizationGeometry::getColor(line.line), thickness);
    }
    protocol.finishScene();

    ++numFramesSentToRenderer;
}

// This is called from the drawing thread if we're buffering, otherwise
// directly from the simulation thread.
void Visualizer::VisualizerRep::drawRealtimeFrameWhenReady
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
void Visualizer::VisualizerRep::reportRealtime(const State& state) {
    const long long t = convertSimTimeToAdjRT(state.getTime()); // scale, convert to ns

    // If this is the first frame, or first since last setMode(), then
    // we synchronize Adjusted Real Time to match. Readjustments will occur
    // if the simulation doesn't keep up with real time.
    if (m_nextFrameDueAdjRT < 0LL) {
        m_adjustedRealTimeBase = realTimeInNs() - t;
        // now getAdjustedRealTime()==tInNs
        m_nextFrameDueAdjRT = t; // now
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
Visualizer::Mode Visualizer::getMode() const {return getRep().m_mode;}

void Visualizer::setDesiredFrameRate(Real framesPerSec) 
{   updRep().setDesiredFrameRate(std::max(framesPerSec, Real(0))); }
Real Visualizer::getDesiredFrameRate() const {return getRep().m_frameRateFPS;}

void Visualizer::setRealTimeScale(Real simTimePerRealSec) 
{   updRep().setRealTimeScale(simTimePerRealSec); }
Real Visualizer::getRealTimeScale() const 
{   return getRep().m_simTimeUnitsPerSec; }

void Visualizer::setDesiredBufferLengthInSec(Real bufferLengthInSec)
{   updRep().setDesiredBufferLengthInSec(bufferLengthInSec); }
Real Visualizer::getDesiredBufferLengthInSec() const
{   return getRep().getDesiredBufferLengthInSec(); }
int Visualizer::getActualBufferLengthInFrames() const 
{   return getRep().getActualBufferLengthInFrames(); }
Real Visualizer::getActualBufferLengthInSec() const 
{   return getRep().getActualBufferLengthInSec(); }


void Visualizer::drawFrameNow(const State& state) 
{   updRep().drawFrameNow(state); }

// The simulation thread normally delivers frames here. Handling is dispatched
// according the current visualization mode.
void Visualizer::report(const State& state) {
    Visualizer::VisualizerRep& rep = updRep();

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

void Visualizer::addEventListener(Visualizer::EventListener* listener) {
    updRep().listeners.push_back(listener);
}

void Visualizer::addFrameController(Visualizer::FrameController* controller) {
    updRep().controllers.push_back(controller);
}

void Visualizer::addDecorationGenerator(DecorationGenerator* generator) {
    updRep().generators.push_back(generator);
}

void Visualizer::setGroundPosition(const CoordinateAxis& axis, Real height) {
    updRep().protocol.setGroundPosition(axis, height);
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

void Visualizer::setCameraTransform(const Transform& transform) const {
    getRep().protocol.setCameraTransform(transform);
}

void Visualizer::zoomCameraToShowAllGeometry() const {
    getRep().protocol.zoomCamera();
}

void Visualizer::pointCameraAt(const Vec3& point, const Vec3& upDirection) const {
    getRep().protocol.lookAt(point, upDirection);
}

void Visualizer::setCameraFieldOfView(Real fov) const {
    getRep().protocol.setFieldOfView(fov);
}

void Visualizer::setCameraClippingPlanes(Real nearPlane, Real farPlane) const {
    getRep().protocol.setClippingPlanes(nearPlane, farPlane);
}

void Visualizer::dumpStats(std::ostream& o) const {getRep().dumpStats(o);}
void Visualizer::clearStats() {updRep().clearStats();}

const Array_<Visualizer::EventListener*>& Visualizer::getEventListeners() const
{   return getRep().listeners; }
const Array_<Visualizer::FrameController*>& Visualizer::getFrameControllers() const
{   return getRep().controllers; }


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
static void* drawingThreadMain(void* visualizerRepAsVoidp) {
    Visualizer::VisualizerRep& vizRep = 
        *reinterpret_cast<Visualizer::VisualizerRep*>(visualizerRepAsVoidp);

    while (true) {
        // Grab the oldest frame in the queue.
        // This will wait if necessary until a frame is available.
        const Frame& frame = vizRep.getOldestFrameInQueue();

        // Draw this frame as soon as its draw time arrives, and readjust
        // adjusted real time if necessary.
        vizRep.drawRealtimeFrameWhenReady
           (frame.state, frame.desiredDrawTimeAdjRT); 

        // Return the now-rendered frame to circulation in the pool. This may
        // wake up the simulation thread if it was waiting for space.
        vizRep.noteThatOldestFrameIsNowAvailable();
    }
    return 0;
}
