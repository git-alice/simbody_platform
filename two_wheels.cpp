/* -------------------------------------------------------------------------- *
 *                   Simbody(tm) Example: Contact Playground                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-13 Stanford University and the Authors.        *
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

/* This example is for experimenting with the new Simbody contact implementation,
which was in beta test in the Simbody 2.1 release, with first official release
in Simbody 2.2. The previous contact implementation is still present and
functional but will be removed soon.

The example shows how the new system tracks contact events and how you can
extract contact forces. It also shows off a number of features of the new
Simbody Visualizer, new in release 2.2. also. Here we display the forces and
torques as colored lines which remain the same color as long as a particular
contact event continues. We also track the energy dissipated by the contacts
and use it to display an energy quantity that should be conserved throughout
the simulation (that is, the current energy plus the dissipated energy
should be a constant).

The simulation uses very expensive, detailed contact surfaces using dense
meshes and the elastic foundation model. Consequently it runs with highly
variable step sizes, and fails to keep up with real time for some short
periods. We use the Visualizer's RealTime mode to buffer up some frames and
smooth out these rough spots so the simulation appears to run at an almost
steady real time rate, displayed at 30fps (depending on how fast your
computer is). Then at the end you can watch the action replay.

You can use this example to see how the different integrators behave when
confronted with a very stiff problem; depending on material properties CPodes
can be *much* faster than the explicit integrators, and it also exhibits very
high stability after the motion damps out. However, by changing material
properties and accuracy setting you can get reasonably good performance out
of the explicit integrators here, which will scale better to large systems.
*/

#include "simbody/Simbody.h"

#include <cstdio>
#include <exception>
#include <algorithm>
#include <iostream>
#include <fstream>

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


static void makeCube(Real h, PolygonalMesh& cube);
static void makeTetrahedron(Real r, PolygonalMesh& tet);
static void makePyramid(Real baseSideLength, PolygonalMesh& pyramid);
static void makeOctahedron(Real radius, PolygonalMesh& pyramid);


int main() {
    try
    {
        // ========== Создание окружения ==========

        MultibodySystem         system;
        SimbodyMatterSubsystem  matter(system);
        GeneralForceSubsystem   forces(system);
        Force::Gravity          gravity(forces, matter, Vec3(0,-10,0));
        Force::MobilityDiscreteForce    torqueDiscreteController;
        Force::MobilityConstantForce    torqueConstantController;

        ContactTrackerSubsystem  tracker(system);
        CompliantContactSubsystem contactForces(system, tracker);
        contactForces.setTrackDissipatedEnergy(true);
        contactForces.setTransitionVelocity(1e-3); // 1e-3

        GeneralContactSubsystem OLDcontact(system);
        const ContactSetIndex OLDcontactSet = OLDcontact.createContactSet();

        // ========== Создание декораций и контактной геомертрии ==========

        // -----------Колеса----------

        PolygonalMesh cylinderMesh_1;
        cylinderMesh_1 = PolygonalMesh::createCylinderMesh(ZAxis, 1, 0.5, 4); // (axis, radius, halfLength)

        ContactGeometry::TriangleMesh cylinder_1(cylinderMesh_1);
        DecorativeMesh showCylinder_1(cylinder_1.createPolygonalMesh());


        PolygonalMesh cylinderMesh_2;
        cylinderMesh_2 = PolygonalMesh::createCylinderMesh(ZAxis, 1, 0.5, 4); // (axis, radius, halfLength)

        ContactGeometry::TriangleMesh cylinder_2(cylinderMesh_1);
        DecorativeMesh showCylinder_2(cylinder_2.createPolygonalMesh());

        // ----------------------------


        // -----------~Вынос----------

        PolygonalMesh removalMesh;
        removalMesh = PolygonalMesh::createCylinderMesh(XAxis, 0.1, 0.5, 4); // (axis, radius, halfLength)

        ContactGeometry::TriangleMesh removal(removalMesh);
        DecorativeMesh showRemoval(removal.createPolygonalMesh());

        // ----------------------------


        // ------------Сфера-----------
        PolygonalMesh sphereMesh;
        sphereMesh = PolygonalMesh::createSphereMesh(1,4);

        ContactGeometry::TriangleMesh sphere(sphereMesh);
        DecorativeMesh showSphere(sphere.createPolygonalMesh());
        Array_<DecorativeLine> normals;
        const Real NormalLength = .02;
        for (int fx=0; fx < sphere.getNumFaces(); ++fx)
            normals.push_back(
                    DecorativeLine(sphere.findCentroid(fx),
                                   sphere.findCentroid(fx)
                                   + NormalLength*sphere.getFaceNormal(fx)));

        ContactCliqueId clique1 = ContactSurface::createNewContactClique();
        ContactCliqueId clique2 = ContactSurface::createNewContactClique();
        ContactCliqueId clique3 = ContactSurface::createNewContactClique();

//        const Real fFac = 0.8; // 1; // to turn off friction
//        const Real fDis = .5*0.2; // to turn off dissipation
//        const Real fVis =  0.1; // .1*.1; // to turn off viscous friction
//        const Real fK = 100*1e6; // pascals

        const Real fFac =1; // to turn off friction
        const Real fDis = .5*0.2; // to turn off dissipation
        const Real fVis =  .1*.1; // to turn off viscous friction
        const Real fK = 100*1e6; // pascals

        // Halfspace floor
        const Rotation R_xdown(-Pi/2,ZAxis);
        matter.Ground().updBody().addDecoration(
                Transform(R_xdown, Vec3(0,-3-.01,0)),
                DecorativeBrick(Vec3(.01,4,8)).setColor(Gray).setOpacity(.1));
        matter.Ground().updBody().addContactSurface(
                Transform(R_xdown, Vec3(0,-3,0)),
                ContactSurface(ContactGeometry::HalfSpace(),
                               ContactMaterial(fK*.1,fDis*.9,fFac*1,fFac*.7,fVis*10))
                        //ContactMaterial(2e6,.01,.1,.05,.01))
                        .joinClique(clique1));


        const Real rad = .4;
        const Real ballMass = 200;
        const Real wheelMass = 200;


        // -----------Колеса----------

        Body::Rigid cylinder_1_Body(MassProperties(wheelMass, Vec3(0),
                                                   wheelMass * UnitInertia::cylinderAlongX(1, 1))); // (Iz, Iy) ?
        cylinder_1_Body.addDecoration(Transform(), showCylinder_1.setColor(Cyan).setOpacity(.2));
        cylinder_1_Body.addContactSurface(
                Transform(),
                ContactSurface(
                        cylinder_1,
                        ContactMaterial(fK*.1,fDis*.9, .1*fFac*.8,.1*fFac*.7, fVis*1),
                        .5 /*thickness*/
                )
                //ContactMaterial(2e6,.01,.1,.05,.01))
                //.joinClique(clique2)
        );

        Body::Rigid cylinder_2_Body(MassProperties(wheelMass, Vec3(0),
                                                   wheelMass * UnitInertia::cylinderAlongX(1, 1))); // (Iz, Iy) ?
        cylinder_2_Body.addDecoration(Transform(), showCylinder_2.setColor(Green).setOpacity(.2));
        cylinder_2_Body.addContactSurface(
                Transform(),
                ContactSurface(
                        cylinder_2,
                        ContactMaterial(fK*.1,fDis*.9, .1*fFac*.8,.1*fFac*.7, fVis*1),
                        .5 /*thickness*/
                )
                //ContactMaterial(2e6,.01,.1,.05,.01))
                //.joinClique(clique2)
        );
        // ----------------------------


        // -----------~Вынос----------

        Body::Rigid removalBody(MassProperties(1e-4, Vec3(0),
                                               ballMass*UnitInertia::cylinderAlongX(1e-5, 1e-5))); // (Iz, Iy) ?
        removalBody.addDecoration(Transform(),
                                  showRemoval.setColor(Cyan).setOpacity(.2));
        removalBody.addContactSurface(Transform(),
                                      ContactSurface(removal,
                                                     ContactMaterial(fK*.1, fDis*.9, .1*fFac*.8,.1*fFac*.7, fVis*1),
                                                     .5)
        );

        // ----------------------------


//        // ------------Сфера-----------
//
//        Body::Rigid ballBody(MassProperties(ballMass, Vec3(0),
//                                            ballMass*UnitInertia::sphere(1)));
//        ballBody.addDecoration(Transform(),
//                               showSphere.setColor(Cyan).setOpacity(.2));
//        ballBody.addDecoration(Transform(),
//                               showSphere.setRepresentation(DecorativeGeometry::DrawWireframe));
//        ballBody.addDecoration(Transform(), DecorativeSphere(1).setColor(Gray)
//                .setOpacity(.1).setResolution(10));
//        ballBody.addContactSurface(Transform(),
//                                   ContactSurface(sphere,
//                                                  ContactMaterial(fK*.1,fDis*.9,
//                                                                  .1*fFac*.8,.1*fFac*.7,fVis*1),
//                                                  .5 /*thickness*/)
//                //ContactMaterial(2e6,.01,.1,.05,.01))
//                //.joinClique(clique2)
//        );
//        // ----------------------------


        // ============================================
        // ========== Инициализация связей ============

//        MobilizedBody::Free ball(matter.Ground(), Transform(Vec3(-2,0,0)), ballBody, Transform(Vec3(0)));

        MobilizedBody::Planar removalFree(matter.Ground(),
                                        Transform(Rotation(-Pi/2, XAxis), Vec3(0, -2, 0)), // !!!!!
                                        removalBody,
                                        Vec3(0));

        MobilizedBody::Pin cylinderFree1(removalFree, Transform(Rotation(Pi/2, YAxis), Vec3(-2,0,0)),
                                         cylinder_1_Body, Transform(Vec3(0))); // Free | matter.Ground()

        MobilizedBody::Pin cylinderFree2(removalFree, Transform(Rotation(Pi/2, YAxis), Vec3(+2,0,0)),
                                         cylinder_2_Body, Transform(Vec3(0))); // Free | matter.Ground()
        // ----------------------------

        // ================================================================
        // ========== Инициализация визцализатора и интерфейса ============

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

        system.addEventReporter(new MyReporter(system,contactForces,ReportInterval));
        system.addEventReporter(new Visualizer::Reporter(viz, ReportInterval));

        // Check for a Run->Quit menu pick every 1/4 second.
        system.addEventHandler(new UserInputHandler(*silo, .25));

        // ========== Инициализация системы и состояния ============

          Force::ConstantTorque(forces, cylinderFree1, Vec3(100, 0, 0));
        // Force::ConstantTorque(forces, cylinderFree2, Vec3(0, 0 , -100));

        // Force::MobilityDiscreteForce torque(forces, cylinderFree2, MobilizerUIndex(0));
        // Force::MobilityConstantForce(forces, cylinderFree2, MobilizerUIndex(0), 0.01);
//        torqueDiscreteController = Force::MobilityDiscreteForce(forces, cylinderFree1, MobilizerUIndex(0), 0.00001);
        torqueConstantController = Force::MobilityConstantForce(forces, cylinderFree1,
                                                                MobilizerUIndex(0), 0);
//        torqueController.setDisabledByDefault(true);

//         Force::ConstantForce force1(forces, cylinderFree1, Vec3(0), Vec3(0, 0, 100));
//        Force::ConstantForce force2(forces, cylinderFree2, Vec3(0), Vec3(0, 0, 0));

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

        // ball.setQToFitTransform(state, Transform(Rotation(Pi/2,XAxis),Vec3(-3,-1.8,0)));

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

        // -----------Цилиндр----------

        // cylinderFree1.setQToFitTransform(state, Transform(Rotation(Pi/2,YAxis)));
        // removalFree.setQToFitTransform(state, Transform(Rotation(Pi/2,ZAxis)));


        // cylinderFree1.setOneU(state, 0, 4);
        // cylinderFree2.setOneU(state, 0, 4);
        // cylinderFree2.setOneU(state, 0, 10);

        // ----------------------------


        // Сумуляция

        // The system as parameterized is very stiff (mostly due to friction)
        // and thus runs best with CPodes which is extremely stable for
        // stiff problems. To get reasonable performance out of the explicit
        // integrators (like the RKs) you'll have to run at a very loose
        // accuracy like 0.1, or reduce the friction coefficients and
        // maybe the stiffnesses.

        // ============ Интегрирование ============

        //ExplicitEulerIntegrator integ(system);
         CPodesIntegrator integ(system,CPodes::BDF,CPodes::Newton);
        // SemiExplicitEuler2Integrator integ(system);
        //RungeKuttaFeldbergIntegrator integ(system);
        //RungeKuttaMersonIntegrator integ(system);
        //RungeKutta3Integrator integ(system);
        //VerletIntegrator integ(system);
        //integ.setMaximumStepSize(1e-0001);
        integ.setAccuracy(1e-3); // minimum for CPodes
        integ.setAllowInterpolation(false);
        //integ.setAccuracy(.01);


        // ------- Автоматическое интегрирование -------

        // TimeStepper ts(system, integ);

        // ts.initialize(state);
        double cpuStart = cpuTime();
        double realStart = realTime();

        // ts.stepTo(20.0);
        // ---------------------------------------------

        // -------------Ручное интегрирование-----------


        const Real MaxStepSize    = Real(1/30.); // 33 1/3 ms (30 Hz)
        const int  DrawEveryN     = 1;           // 33 1/3 ms frame update (30 Hz)

        integ.initialize(state);
        unsigned stepNum = 0;
        while (true)
        {
            State& localState = integ.updAdvancedState();

            cout << "t=" << localState.getTime()
                 << " q_c1=" << cylinderFree1.getQAsVector(localState)
                 << " u_c1=" << cylinderFree1.getUAsVector(localState)
                 << " q_c2=" << cylinderFree2.getQAsVector(localState)
                 << " u_c2=" << cylinderFree2.getUAsVector(localState)
                 << endl;


//            torqueConstantController.setForce(localState, -40);

//            Real torqueValue = torqueConstantController.getForce(localState);
//            cout << "torque " << torqueValue << endl;
            if (stepNum > 30) {
                // torqueController.enable(localState);
                torqueConstantController.setForce(localState, 0);
                system.realize(localState, Stage::Velocity);

            }
//            torqueController.enable(state);
//            torqueController.setForce(localState, -100);

            //if (stepNum % DrawEveryN == 0) {
            viz.report(localState);
            //}

            integ.stepBy(MaxStepSize);
            cout << stepNum <<  endl;
//            cout << "  " << torqueValue << endl;
            ++stepNum;
        }

        // ---------------------------------------------


        const double timeInSec = realTime() - realStart;
        const int evals = integ.getNumRealizations();
        // cout << "Done -- took " << integ.getNumStepsTaken() << " steps in " <<
        //     timeInSec << "s elapsed for " << ts.getTime() << "s sim (avg step="
        //     << (1000*ts.getTime())/integ.getNumStepsTaken() << "ms) "
        //     << (1000*ts.getTime())/evals << "ms/eval\n";
        cout << "  CPU time was " << cpuTime() - cpuStart << "s\n";

        printf("Using Integrator %s at accuracy %g:\n", integ.getMethodName(), integ.getAccuracyInUse());
        printf("# STEPS/ATTEMPTS = %d/%d\n", integ.getNumStepsTaken(), integ.getNumStepsAttempted());
        printf("# ERR TEST FAILS = %d\n", integ.getNumErrorTestFailures());
        printf("# REALIZE/PROJECT = %d/%d\n", integ.getNumRealizations(), integ.getNumProjections());

        viz.dumpStats(std::cout);

        // Add as slider to control playback speed.
        viz.addSlider("Speed", 1, 0, 4, 1);
        viz.setMode(Visualizer::PassThrough);

        silo->clear(); // forget earlier input
        double speed = 1; // will change if slider moves
        while(true) {
            cout << "Choose Run/Replay to see that again ...\n";

            int menuId, item;
            silo -> waitForMenuPick(menuId, item);


            if (menuId != RunMenuId) {
                cout << "\aUse the Run menu!\n";
                continue;
            }

            if (item == QuitItem)
                break;
            if (item != ReplayItem) {
                cout << "\aHuh? Try again.\n";
                continue;
            }

            for (double i=0; i < (int)saveEm.size(); i += speed ) {
                int slider; Real newValue;
                if (silo->takeSliderMove(slider,newValue)) {
                    speed = newValue;
                }
                viz.report(saveEm[(int)i]);
            }
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