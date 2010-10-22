#include "SimTKsimbody.h"
#include "simbody/internal/Visualizer.h"
#include "simbody/internal/VisualizationEventListener.h"
#include "simbody/internal/VisualizationGeometry.h"
#include "simbody/internal/VisualizationProtocol.h"
#include <cstdlib>
#include <cstdio>
#include <pthread.h>
#include <string>

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

class Visualizer::VisualizerRep {
public:
    VisualizerRep(Visualizer* owner, MultibodySystem& system) : handle(owner), system(system), protocol(*owner) {
    }
    Visualizer* handle;
    MultibodySystem& system;
    Array_<DecorativeGeometry> addedGeometry;
    VisualizationProtocol protocol;
    Array_<VisualizationEventListener*> listeners;
    Array_<RubberBandLine> lines;
};

Visualizer::Visualizer(MultibodySystem& system) : rep(new VisualizerRep(this, system)) {
}

Visualizer::~Visualizer() {
    if (rep->handle == this)
        delete rep;
}

void Visualizer::report(const State& state) const {
    MultibodySystem& system = updRep().system;
    system.realize(state, Stage::Position);
    Array_<DecorativeGeometry> geometry;
    for (Stage stage = Stage::Topology; stage < state.getSystemStage(); ++stage)
        system.calcDecorativeGeometryAndAppend(state, stage, geometry);
    VisualizationProtocol& protocol = updRep().protocol;
    protocol.beginScene();
    VisualizationGeometry geometryCreator(protocol, system.getMatterSubsystem(), state);
    for (int i = 0; i < (int) geometry.size(); ++i)
        geometry[i].implementGeometry(geometryCreator);
    const Array_<DecorativeGeometry>& addedGeometry = getRep().addedGeometry;
    for (int i = 0; i < (int) addedGeometry.size(); ++i)
        addedGeometry[i].implementGeometry(geometryCreator);
    const Array_<RubberBandLine>& lines = getRep().lines;
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

void Visualizer::addEventListener(VisualizationEventListener* listener) {
    updRep().listeners.push_back(listener);
}

const Array_<VisualizationEventListener*>& Visualizer::getEventListeners() const {
    return getRep().listeners;
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
