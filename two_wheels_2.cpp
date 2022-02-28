#include "simbody/Simbody.h"

#include <cstdio>
#include <exception>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>

using std::cout; using std::endl;

using namespace SimTK;

Array_<State> saveEm;

static const Real TimeScale = 1;
static const Real FrameRate = 30;
static const Real ReportInterval = TimeScale/FrameRate;
static const Real ForceScale = .25;
static const Real MomentScale = .5;


class ForceArrowGenerator : public DecorationGenerator {
public:
    ForceArrowGenerator(const MultibodySystem& system,
                        const CompliantContactSubsystem& complCont)
            :   m_system(system), m_compliant(complCont) {}

    virtual void generateDecorations(const State& state, Array_<DecorativeGeometry>& geometry) override {
        const Vec3 frcColors[] = {Red,Orange,Cyan};
        const Vec3 momColors[] = {Blue,Green,Purple};
        m_system.realize(state, Stage::Velocity);

        const int ncont = m_compliant.getNumContactForces(state);
        for (int i=0; i < ncont; ++i) {
            const ContactForce& force = m_compliant.getContactForce(state,i);
            const ContactId     id    = force.getContactId();
            const Vec3& frc = force.getForceOnSurface2()[1];
            const Vec3& mom = force.getForceOnSurface2()[0];
            Real  frcMag = frc.norm(), momMag=mom.norm();
            int frcThickness = 1, momThickness = 1;
            Real frcScale = ForceScale, momScale = ForceScale;
            while (frcMag > 10)
                frcThickness++, frcScale /= 10, frcMag /= 10;
            while (momMag > 10)
                momThickness++, momScale /= 10, momMag /= 10;
            DecorativeLine frcLine(force.getContactPoint(),
                                   force.getContactPoint() + frcScale*frc);
            DecorativeLine momLine(force.getContactPoint(),
                                   force.getContactPoint() + momScale*mom);
            frcLine.setColor(frcColors[id%3]);
            momLine.setColor(momColors[id%3]);
            frcLine.setLineThickness(2*frcThickness);
            momLine.setLineThickness(2*momThickness);
            geometry.push_back(frcLine);
            geometry.push_back(momLine);

            ContactPatch patch;
            const bool found = m_compliant.calcContactPatchDetailsById(state,id,patch);
            //cout << "patch for id" << id << " found=" << found << endl;
            //cout << "resultant=" << patch.getContactForce() << endl;
            //cout << "num details=" << patch.getNumDetails() << endl;
            for (int i=0; i < patch.getNumDetails(); ++i) {
                const ContactDetail& detail = patch.getContactDetail(i);
                const Real peakPressure = detail.getPeakPressure();
                // Make a black line from the element's contact point in the normal
                // direction, with length proportional to log(peak pressure)
                // on that element.
                DecorativeLine normal(detail.getContactPoint(),
                                      detail.getContactPoint()+ std::log10(peakPressure)
                                                                * detail.getContactNormal());
                normal.setColor(Black);
                geometry.push_back(normal);
                // Make a red line that extends from the contact
                // point in the direction of the slip velocity, of length 3*slipvel.
                DecorativeLine slip(detail.getContactPoint(),
                                    detail.getContactPoint()+3*detail.getSlipVelocity());
                slip.setColor(Red);
                geometry.push_back(slip);
            }
        }
    }
private:
    const MultibodySystem&              m_system;
    const CompliantContactSubsystem&    m_compliant;
};

class MyReporter : public PeriodicEventReporter {
public:
    MyReporter(const MultibodySystem& system,
               const CompliantContactSubsystem& complCont,
               Real reportInterval)
            :   PeriodicEventReporter(reportInterval), m_system(system),
                m_compliant(complCont)
    {}

    ~MyReporter() {}

    void handleEvent(const State& state) const override {
        m_system.realize(state, Stage::Dynamics);
        cout << state.getTime() << ": E = " << m_system.calcEnergy(state)
             << " Ediss=" << m_compliant.getDissipatedEnergy(state)
             << " E+Ediss=" << m_system.calcEnergy(state)
                               +m_compliant.getDissipatedEnergy(state)
             << endl;
        const int ncont = m_compliant.getNumContactForces(state);
        cout << "Num contacts: " << m_compliant.getNumContactForces(state) << endl;
        for (int i=0; i < ncont; ++i) {
            const ContactForce& force = m_compliant.getContactForce(state,i);
            //cout << force;
        }
        saveEm.push_back(state);
    }
private:
    const MultibodySystem&           m_system;
    const CompliantContactSubsystem& m_compliant;
};

// These are the item numbers for the entries on the Run menu.
static const int RunMenuId = 3, HelpMenuId = 7;
static const int GoItem = 1, ReplayItem=2, QuitItem=3;

// This is a periodic event handler that interrupts the simulation on a regular
// basis to poll the InputSilo for user input. If there has been some, process it.
// This one does nothing but look for the Run->Quit selection.
class UserInputHandler : public PeriodicEventHandler {
public:
    UserInputHandler(Visualizer::InputSilo& silo, Real interval)
            :   PeriodicEventHandler(interval), m_silo(silo) {}

    virtual void handleEvent(State& state, Real accuracy,
                             bool& shouldTerminate) const override
    {
        int menuId, item;
        if (m_silo.takeMenuPick(menuId, item) && menuId==RunMenuId && item==QuitItem)
            shouldTerminate = true;
    }
private:
    Visualizer::InputSilo& m_silo;
};


int main() {
    time_t now = time(0);
    char* now_str = std::ctime(&now);

    try
    {
        MultibodySystem         system;
        SimbodyMatterSubsystem  matter(system);
        GeneralForceSubsystem   forces(system);
        Force::Gravity          gravity(forces, matter, Vec3(0,-10,0));
        Force::MobilityDiscreteForce    torqueDiscreteController;
        Force::MobilityConstantForce    torqueConstantController;

        ContactTrackerSubsystem  tracker(system);
        CompliantContactSubsystem contactForces(system, tracker);

        contactForces.setTrackDissipatedEnergy(true);
        contactForces.setTransitionVelocity(1e-3);

        //
        PolygonalMesh cylinderMesh_1;
        cylinderMesh_1 = PolygonalMesh::createCylinderMesh(ZAxis, 1, 0.5, 4);
        ContactGeometry::TriangleMesh cylinder_1(cylinderMesh_1);
        DecorativeMesh showCylinder_1(cylinder_1.createPolygonalMesh());

        PolygonalMesh cylinderMesh_2;
        cylinderMesh_2 = PolygonalMesh::createCylinderMesh(ZAxis, 1, 0.5, 4);
        ContactGeometry::TriangleMesh cylinder_2(cylinderMesh_1);
        DecorativeMesh showCylinder_2(cylinder_2.createPolygonalMesh());

        PolygonalMesh removalMesh;
        removalMesh = PolygonalMesh::createCylinderMesh(XAxis, 0.1, 0.5, 4);
        ContactGeometry::TriangleMesh removal(removalMesh);
        DecorativeMesh showRemoval(removal.createPolygonalMesh());

        PolygonalMesh platformMesh;
        removalMesh = PolygonalMesh::createCylinderMesh(XAxis, 0.1, 0.5, 4);
        ContactGeometry::TriangleMesh platform(removalMesh);
        DecorativeMesh showPlatform(platform.createPolygonalMesh());

        PolygonalMesh pinMesh;
        removalMesh = PolygonalMesh::createCylinderMesh(XAxis, 0.1, 0.5, 4);
        ContactGeometry::TriangleMesh pin(removalMesh);
        DecorativeMesh showPin(platform.createPolygonalMesh());


        const Real fFac = 1;      // friction
        const Real fDis = .5*0.2; // dissipation
        const Real fVis = .1*.1;  // viscous friction
        const Real fK = 100*1e6;  // pascals
        const Real rad = .4;
        const Real ballMass = 200;
        const Real wheelMass = 200;

        matter.Ground().updBody().addDecoration(Transform(Rotation(-Pi/2, ZAxis),Vec3(0,-3-.01,0)),DecorativeBrick(Vec3(.01,4,8)).setColor(Gray).setOpacity(.1));
        matter.Ground().updBody().addContactSurface(Transform(Rotation(-Pi/2, ZAxis), Vec3(0,-3,0)),ContactSurface(ContactGeometry::HalfSpace(),ContactMaterial(fK*.1,fDis*.9,fFac*1,fFac*.7,fVis*10)));

        Body::Rigid cylinder_1_Body(MassProperties(wheelMass, Vec3(0),wheelMass * UnitInertia::cylinderAlongX(1, 0.5))); // (Iz, Iy) ?
        cylinder_1_Body.addDecoration(Transform(), showCylinder_1.setColor(Cyan).setOpacity(.2));
        cylinder_1_Body.addContactSurface(Transform(),ContactSurface(cylinder_1,ContactMaterial(fK*.1,fDis*.9, .1*fFac*.8,.1*fFac*.7, fVis*1),.5 /*thickness*/));

        Body::Rigid cylinder_2_Body(MassProperties(wheelMass, Vec3(0),wheelMass * UnitInertia::cylinderAlongX(1, 0.5))); // (Iz, Iy) ?
        cylinder_2_Body.addDecoration(Transform(), showCylinder_2.setColor(Green).setOpacity(.2));
        cylinder_2_Body.addContactSurface(Transform(),ContactSurface(cylinder_2,ContactMaterial(fK*.1,fDis*.9, .1*fFac*.8,.1*fFac*.7, fVis*1),.5 /*thickness*/));

        Body::Rigid removalBody(MassProperties(1e-4, Vec3(0),ballMass*UnitInertia::cylinderAlongX(1e-5, 1e-5))); // (Iz, Iy) ?
        removalBody.addDecoration(Transform(),showRemoval.setColor(Cyan).setOpacity(.2));
        removalBody.addContactSurface(Transform(),ContactSurface(removal,ContactMaterial(fK*.1, fDis*.9, .1*fFac*.8,.1*fFac*.7, fVis*1),.5));

        Body::Rigid platformBody(MassProperties(1e-4, Vec3(0),ballMass*UnitInertia::cylinderAlongX(1e-5, 1e-5))); // (Iz, Iy) ?
        platformBody.addDecoration(Transform(),showPlatform.setColor(Cyan).setOpacity(.2));
        platformBody.addContactSurface(Transform(),ContactSurface(platform,ContactMaterial(fK*.1, fDis*.9, .1*fFac*.8,.1*fFac*.7, fVis*1),.5));

        Body::Rigid pinBody(MassProperties(1e-4, Vec3(0),ballMass*UnitInertia::cylinderAlongX(1e-5, 1e-5))); // (Iz, Iy) ?
        pinBody.addDecoration(Transform(),showPin.setColor(Cyan).setOpacity(.2));
        pinBody.addContactSurface(Transform(),ContactSurface(platform,ContactMaterial(fK*.1, fDis*.9, .1*fFac*.8,.1*fFac*.7, fVis*1),.5));


        // 1 вар - еще одно тело [x]
        // 2 вар - custom MobilizedBody [x]
        // 3 вар - добавить тело Translation + Pin [v]

        // Два колеса
        MobilizedBody::Translation platformFree(matter.Ground(),Transform (Rotation(-Pi/2, XAxis), Vec3(0, 0, 0)), platformBody,Vec3(0));
//        MobilizedBody::Pin pinFree(platformFree, Transform( Rotation(Pi/2, YAxis), Vec3(0, 0, -1)), pinBody, Transform (Vec3(0)));  // ???

        MobilizedBody::Planar removalFree(platformFree,Transform( Vec3(0, 0, -1)), removalBody,Vec3(0));
        MobilizedBody::Pin cylinderFree1(removalFree, Transform(Rotation(Pi/2, YAxis), Vec3(-2,0,0)), cylinder_1_Body, Transform(Vec3(0)));
        MobilizedBody::Pin cylinderFree2(removalFree, Transform(Rotation(Pi/2, YAxis), Vec3(+2,0,0)), cylinder_2_Body, Transform(Vec3(0)));

        Visualizer viz(system);
        viz.addDecorationGenerator(new ForceArrowGenerator(system,contactForces));
        viz.setMode(Visualizer::RealTime);
        viz.setDesiredBufferLengthInSec(1);
        viz.setDesiredFrameRate(FrameRate);
        viz.setGroundHeight(-3);
        viz.setShowShadows(true);

        Visualizer::InputSilo* silo = new Visualizer::InputSilo();
        viz.addInputListener(silo);
        Array_<std::pair<String,int> > runMenuItems;
        runMenuItems.push_back(std::make_pair("Go", GoItem));
        runMenuItems.push_back(std::make_pair("Replay", ReplayItem));
        runMenuItems.push_back(std::make_pair("Quit", QuitItem));
        viz.addMenu("Run", RunMenuId, runMenuItems);

        Array_<std::pair<String,int> > helpMenuItems;
        helpMenuItems.push_back(std::make_pair("TBD - Sorry!", 1));
        viz.addMenu("Help", HelpMenuId, helpMenuItems);

        system.addEventReporter(new MyReporter(system, contactForces, ReportInterval));
        system.addEventReporter(new Visualizer::Reporter(viz, ReportInterval));
        system.addEventHandler(new UserInputHandler(*silo, .25));

        // === Силы и моменты ===
        // Force::ConstantTorque(forces, cylinderFree1, Vec3(1000, 0, 0));
        // Force::ConstantForce force1(forces, cylinderFree1, Vec3(0), Vec3(0, 0, 1000));
        // Force::ConstantForce force2(forces, cylinderFree2, Vec3(0), Vec3(100, 0, 0));
        // Force::MobilityDiscreteForce torque(forces, cylinderFree2, MobilizerUIndex(0), 0);
        // Force::MobilityConstantForce(forces, cylinderFree2, MobilizerUIndex(0), 0.01);
        // torqueController.setDisabledByDefault(true);

        // Постоянные моменты
        torqueDiscreteController = Force::MobilityDiscreteForce(forces, cylinderFree1, MobilizerUIndex(0), 100);
//        torqueConstantController = Force::MobilityConstantForce(forces, cylinderFree2,MobilizerUIndex(0), 100);

        system.realizeTopology();

        // Show ContactSurfaceIndex for each contact surface
        for (MobilizedBodyIndex mbx(0); mbx < matter.getNumBodies(); ++mbx) {
            const MobilizedBody& mobod = matter.getMobilizedBody(mbx);
            const int nsurfs = mobod.getBody().getNumContactSurfaces();
            printf("mobod %d has %d contact surfaces\n", (int)mbx, nsurfs);
            for (int i=0; i<nsurfs; ++i) {
                printf("%2d: index %d\n", i,
                       (int)tracker.getContactSurfaceIndex(mbx,i));
            }
        }

        State state = system.getDefaultState();

        viz.report(state);
        printf("Default state\n");
        cout << "t=" << state.getTime()
             << " q=" << cylinderFree1.getQAsVector(state)
             << " u=" << cylinderFree1.getUAsVector(state)
             << endl;

        cout << "\nChoose 'Go' from Run menu to simulate:\n";
        int menuId, item;
        do { silo -> waitForMenuPick(menuId, item);
            if (menuId != RunMenuId || item != GoItem)
                cout << "\aDude ... follow instructions!\n";
        } while (menuId != RunMenuId || item != GoItem);


        //ExplicitEulerIntegrator integ(system);
        // CPodesIntegrator integ(system,CPodes::BDF,CPodes::Newton);
        SemiExplicitEuler2Integrator integ(system);
        //RungeKuttaFeldbergIntegrator integ(system);
        //RungeKuttaMersonIntegrator integ(system);
        //RungeKutta3Integrator integ(system);
        //VerletIntegrator integ(system);
        //integ.setMaximumStepSize(1e-0001);
        integ.setAccuracy(1e-3); // minimum for CPodes
        integ.setAllowInterpolation(false);
        //integ.setAccuracy(.01);


        // Автоматическое интегрирование
        // TimeStepper ts(system, integ);
        // ts.initialize(state);
        // ts.stepTo(20.0);


        // Ручное интегрирование
        const Real MaxStepSize    = Real(1/30.); // 33 1/3 ms (30 Hz)

        std::ofstream file;
        file.open ("/Users/alexeyzhuchkov/CLionProjects/simbody/data/" + std::string(now_str) + ".txt", std::ios_base::app);
        file << "t;cyl1_qs;cyl1_us;cyl2_qs;cyl2_us;removal_qs;removal_us;platform_qs" << endl;
        file.close();

        integ.initialize(state);
        unsigned stepNum = 0;
        while (true)
        {

            State& localState = integ.updAdvancedState();

            std::ofstream file;
            file.open ("/Users/alexeyzhuchkov/CLionProjects/simbody/data/" + std::string(now_str) + ".txt", std::ios_base::app);
            file << localState.getTime() << ";"
                 << cylinderFree1.getQAsVector(localState) << ";"
                 << cylinderFree1.getUAsVector(localState) << ";"
                 << cylinderFree2.getQAsVector(localState) << ";"
                 << cylinderFree2.getUAsVector(localState) << ";"
                 << removalFree.getQAsVector(localState) << ";"
                 << removalFree.getUAsVector(localState) << ";"
                 << platformFree.getQAsVector(localState) << ";"
                 << endl;
            file.close();

//            Постоянный момент на каждом шаге
            torqueDiscreteController.setMobilityForce(localState, 50);
            system.realize(localState, Stage::Velocity);

            viz.report(localState);

            integ.stepBy(MaxStepSize);
            cout << stepNum <<  endl;
            ++stepNum;

        }


    } catch (const std::exception& e) {
        std::printf("EXCEPTION THROWN: %s\n", e.what());
        exit(1);

    } catch (...) {
        std::printf("UNKNOWN EXCEPTION THROWN\n");
        exit(1);
    }

    return 0;
}