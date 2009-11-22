#include "SimTKsimbody.h"        //SimTK principal Library
#include "SimTKsimbody_aux.h"    //for VTK stuff
using namespace SimTK;

int main()
{

        MultibodySystem system;
        SimbodyMatterSubsystem matter(system); //Subsystem that contains all the bodies of the system
        GeneralForceSubsystem forces(system);  //Subsystem that contains all the forces of the system


        Mat33 inertiaMatrix(1,0,0,0,1,0,0,0,1);
        Vec3 centerOfMass(0);

        //Body Creations
        Body::Rigid body(MassProperties(120, centerOfMass, Inertia(inertiaMatrix)));
        Body::Rigid arm(MassProperties(15, centerOfMass, Inertia(inertiaMatrix)));
        Body::Rigid foreArm(MassProperties(15, centerOfMass, Inertia(inertiaMatrix)));
        Body::Rigid hand(MassProperties(15, centerOfMass, Inertia(inertiaMatrix)));

        MobilizedBody::Weld attach(matter.Ground(), Transform(Vec3(0,0,0)), body, Transform(Vec3(0,0,0)));
        MobilizedBody::Ball shoulder(attach, Transform(Vec3(2,3,0)), arm, Transform(Vec3(0,0,0)));
        MobilizedBody::Ball elbow(shoulder, Transform(Vec3(2,0,0)), foreArm, Transform(Vec3(0,0,0)));
        MobilizedBody::Ball wrist(elbow, Transform(Vec3(2,0,0)), hand, Transform(Vec3(0,0,0)));
        MobilizedBody::Ball wrist2(wrist, Transform(Vec3(2,0,0)), hand, Transform(Vec3(0,0,0)));

        matter.Ground().addBodyDecoration(Vec3(2,7,0), DecorativeFrame().setColor(Green));
        elbow.addBodyDecoration(Vec3(0), DecorativeText("foreArm").setScale(.25));
        try
        {
            // Initialize the system and state.
            system.realizeTopology(); 
            State state = system.getDefaultState();

            VTKVisualizer viz(system);
            viz.report(state);
            char ch=getchar();

            std::vector<MobilizedBodyIndex> indexVector;
            std::vector< std::vector<Vec3> > originPositionVector;
            std::vector< std::vector<Vec3> > targetPositionVector;

            std::vector<Vec3> originPositionOfForearm, originPositionOfHand;
            originPositionOfForearm.push_back(Vec3(0));
            originPositionOfHand.push_back(Vec3(0));

            std::vector<Vec3> targetPositionOfForearm, targetPositionOfHand;
            targetPositionOfForearm.push_back(Vec3(2,7,0));
            targetPositionOfHand.push_back(Vec3(4,7,0));

            // UNCOMMENTING LINES BELOW AVOIDS BUG
            indexVector.push_back(elbow.getMobilizedBodyIndex());
            //indexVector.push_back(wrist2.getMobilizedBodyIndex());
            originPositionVector.push_back(originPositionOfForearm);
            //originPositionVector.push_back(originPositionOfHand);
            targetPositionVector.push_back(targetPositionOfForearm);
            //targetPositionVector.push_back(targetPositionOfHand);
            
            double rms=ObservedPointFitter::findBestFit(system, state, 
                                                        indexVector, 
                                                        originPositionVector, 
                                                        targetPositionVector, 0.001);
            system.realize(state,Stage::Position);
            std::cout<<"forearm origin now at "<<elbow.getBodyOriginLocation(state)<<std::endl;
            std::cout<<"hand origin now at "<<wrist2.getBodyOriginLocation(state)<<std::endl;
            std::cout<<"Rms is: "<<rms<<" ENTER when done\n";
            viz.report(state);
            ch=getchar();
        }
        catch (const std::exception& e) 
        {
            std::cout<<"Exception: "<<std::string(e.what())+"\n";
        }
            
        

    return 0;
}

