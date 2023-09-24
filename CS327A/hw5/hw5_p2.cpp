#include <iostream>
#include <string>
#include <csignal>
#include <thread>
#include <chrono>
#include <math.h>

#include "Sai2Model.h"
#include "Sai2Graphics.h"
#include "Sai2Simulation.h"
#include "redis/RedisClient.h"
#include "timer/LoopTimer.h"

#include <dynamics3d.h>
#include "timer/LoopTimer.h"
#include <GLFW/glfw3.h> 

using namespace std;

const double time_slowdown_factor = 10;

constexpr const char *sim_title = "cs327a - HW 5 Problem 2";
const string world_file = "./resources/world_hw5_p2.urdf";
const string robot_file = "../../resources/puma/puma_gripper.urdf";
const string robot1_name = "Puma1";
const string robot2_name = "Puma2";
const string object_name = "CoordObject";
const string object_file = "./resources/object.urdf";
const string object_link_name = "object";
const string ee_link_name = "end-effector";
const string gripper_joint_name = "gripper";

const string camera_name = "camera_front";

RedisClient redis_client;

const std::string OBJECT_CARTESIAN_POSITION_KEY = "cs327a::hw5_p2::object::cartesian_pos";

/* ----------------------------------------------------------------------------------
	Simulation and Control Loop Setup
-------------------------------------------------------------------------------------*/
// state machine setup
enum ControlMode {
	CONTROL_GRASP_STABILIZE = 0,
	CONTROL_AUGMENTED_OBJECT
};

// simulation loop
bool fSimulationRunning = false;
void control(Sai2Model::Sai2Model* robot1, Sai2Model::Sai2Model* robot2, Sai2Model::Sai2Model* object_model, Simulation::Sai2Simulation* sim);
void simulation(Simulation::Sai2Simulation* sim);

// initialize window manager
GLFWwindow* glfwInitialize();

// callback to print glfw errors
void glfwError(int error, const char* description);

// callback when a key is pressed
void keySelect(GLFWwindow* window, int key, int scancode, int action, int mods);

/* =======================================================================================
   MAIN LOOP
========================================================================================== */
int main (int argc, char** argv) {

	cout << "Loading URDF world model file: " << world_file << endl;

	// start redis client
	redis_client = RedisClient();
	redis_client.connect();

 	// set up redis callbacks
    redis_client.createWriteCallback(0);

	// load graphics scene
	auto graphics = new Sai2Graphics::Sai2Graphics(world_file, false);

	// load robots
	auto robot1 = new Sai2Model::Sai2Model(robot_file, false);
	auto robot2 = new Sai2Model::Sai2Model(robot_file, false);

	// load object
	auto coobject = new Sai2Model::Sai2Model(object_file, false);

	// load simulation world
	auto sim = new Simulation::Sai2Simulation(world_file, false);
	// set co-efficient of restition to zero to avoid bounce
    // see issue: https://github.com/manips-sai/sai2-simulation/issues/1
    sim->setCollisionRestitution(0.0);
    sim->setCoeffFrictionStatic(0.5);
    sim->setCoeffFrictionDynamic(0.5);

    // set joint damping on grippers: 
    auto base_1 = sim->_world->getBaseNode(robot1_name);
    auto gripper_1 = base_1->getJoint(gripper_joint_name);
    gripper_1->setDamping(10.0);
    gripper_1->setJointLimits(-0.005, 0.068, 0.005);
    auto base_2 = sim->_world->getBaseNode(robot2_name);
    auto gripper_2 = base_2->getJoint(gripper_joint_name);
    gripper_2->setDamping(10.0);
    gripper_2->setJointLimits(-0.005, 0.068, 0.005);

    // set initial conditions
	robot1->_q << 90/180.0*M_PI,
				-22.5/180.0*M_PI,
				212/180.0*M_PI,
				90.0/180.0*M_PI,
				100/180.0*M_PI,
				180/180.0*M_PI,
				0.0;
	sim->setJointPositions(robot1_name, robot1->_q);
	robot1->updateModel();
	robot2->_q << 90/180.0*M_PI,
				202.5/180.0*M_PI,
				-28/180.0*M_PI,
				-90.0/180.0*M_PI,
				97/180.0*M_PI,
				180/180.0*M_PI,
				0.0;
	sim->setJointPositions(robot2_name, robot2->_q);
	robot2->updateModel();
	Eigen::Affine3d ee_trans;

	// set up error callback
    glfwSetErrorCallback(glfwError);

    // initialize GLFW
    glfwInit();

    // retrieve resolution of computer display and position window accordingly
    GLFWmonitor *primary = glfwGetPrimaryMonitor();
    const GLFWvidmode *mode = glfwGetVideoMode(primary);

    // information about computer screen and GLUT display window
    int screenW = mode->width;
    int screenH = mode->height;
    int windowW = 0.8 * screenH;
    int windowH = 0.5 * screenH;
    int windowPosY = (screenH - windowH) / 2;
    int windowPosX = windowPosY;

    // create window and make it current
    glfwWindowHint(GLFW_VISIBLE, 0);
    GLFWwindow *window = glfwCreateWindow(windowW, windowH, sim_title, NULL, NULL);
    glfwSetWindowPos(window, windowPosX, windowPosY);
    glfwShowWindow(window);
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // set callbacks
    glfwSetKeyCallback(window, keySelect);

    // cache variables
    double last_cursorx, last_cursory;

	// start the simulation thread first
    fSimulationRunning = true;
	thread sim_thread(simulation, sim);

	// next start the control thread
	thread ctrl_thread(control, robot1, robot2, coobject, sim);
	
    // while window is open:
    while (!glfwWindowShouldClose(window)) {

		// update graphics. this automatically waits for the correct amount of time
		int width, height;
		glfwGetFramebufferSize(window, &width, &height);
		graphics->updateGraphics(robot1_name, robot1);
		graphics->updateGraphics(robot2_name, robot2);
		graphics->updateGraphics(object_name, coobject);
		graphics->render(camera_name, width, height);
		glfwSwapBuffers(window);
		glFinish();

	    // poll for events
	    glfwPollEvents();
	}

	// stop simulation
	fSimulationRunning = false;
	sim_thread.join();
	ctrl_thread.join();

    // destroy context
    glfwDestroyWindow(window);

    // terminate
    glfwTerminate();

	return 0;
}

/* ----------------------------------------------------------------------------------
	Utility functions
-------------------------------------------------------------------------------------*/
// Calculate the cross product matrix
Eigen::Matrix3d getCrossProductMat(const Eigen::Vector3d& t) {
	Eigen::Matrix3d ret;
	ret <<  0, -t(2), t(1),
			t(2), 0, -t(0),
			-t(1), t(0), 0;
	return ret;
}

/* =======================================================================================
   CONTROL LOOP
========================================================================================== */
void control(Sai2Model::Sai2Model* robot1, Sai2Model::Sai2Model* robot2, Sai2Model::Sai2Model* object_model, Simulation::Sai2Simulation* sim) {
	
	// create a timer
	LoopTimer timer;
	timer.initializeTimer();
	timer.setLoopFrequency(1000); //1000Hz timer
	double last_time = timer.elapsedTime()/time_slowdown_factor; //secs

	// load robot global frame to robot base transformations: todo: move out to a function
	Eigen::Affine3d robot1_base_frame = sim->getRobotBaseTransform(robot1_name);
	Eigen::Affine3d robot2_base_frame = sim->getRobotBaseTransform(robot2_name);

	// cache variables
	bool fTimerDidSleep = true;
	bool fTorqueUseDefaults = false; // set true when torques are overriden for the first time
	Eigen::VectorXd tau1 = Eigen::VectorXd::Zero(robot1->dof());
	Eigen::VectorXd tau2 = Eigen::VectorXd::Zero(robot2->dof());

	Eigen::Affine3d object_com_frame;
	Eigen::Vector3d object_current_pos;
	Eigen::MatrixXd object_inertia(6,6);
	Eigen::MatrixXd object_j(6,6);
	Eigen::VectorXd object_p(6);

	Eigen::VectorXd robot1_g(robot1->dof());
	Eigen::VectorXd robot2_g(robot2->dof());

	// **  Other sugested variables and gains **
	// Eigen::Vector3d object_com_in_robot1_ee_frame;
	// Eigen::Vector3d object_com_in_robot2_ee_frame;
	// Eigen::MatrixXd robot1_j0_objcom(6, robot1->dof());
	// Eigen::MatrixXd robot2_j0_objcom(6, robot2->dof());
	// Eigen::MatrixXd robot1_j0_objectcom_bar(robot1->dof(), 6);
	// Eigen::MatrixXd robot2_j0_objectcom_bar(robot2->dof(), 6);
	// Eigen::MatrixXd robot1_objcom_inertia(6,6);
	// Eigen::MatrixXd robot2_objcom_inertia(6,6);

	// Eigen::MatrixXd augmented_object_inertia(6,6);
	// Eigen::VectorXd augmented_object_p(6);

	// Eigen::MatrixXd G(2*6, 2*6);
	// Eigen::MatrixXd W(6, 2*6);

	// Eigen::Vector3d obj_des_pos;
	// Eigen::Vector3d obj_ori_error;
	// Eigen::VectorXd obj_task_err(6);
	// Eigen::VectorXd force_des_vec(12);
	// Eigen::VectorXd force_ee_vec(12);

	double kp = 30;
	double kv = 10;

	// ** Control Mode **
	// 		0 = grasp stabilizing controller
	//		1 = augmented object controller
	ControlMode control_mode = CONTROL_GRASP_STABILIZE; 

	int count = 0; 

	while (fSimulationRunning) { //automatically set to false when simulation is quit
		fTimerDidSleep = timer.waitForNextLoop();

		// update time
		double curr_time = timer.elapsedTime()/time_slowdown_factor;
		double loop_dt = curr_time - last_time;

        // read joint positions, velocities
        sim->getJointPositions(robot1_name, robot1->_q);
        sim->getJointVelocities(robot1_name, robot1->_dq);
        robot1->updateModel();
        sim->getJointPositions(robot2_name, robot2->_q);
        sim->getJointVelocities(robot2_name, robot2->_dq);
        robot2->updateModel();

        // read object position
        sim->getJointPositions(object_name, object_model->_q);
        sim->getJointVelocities(object_name, object_model->_dq);
        object_model->updateModel();

        // update object dynamics
		// - find object COM frame in global frame
		object_model->transform(object_com_frame, object_link_name);
		object_current_pos << object_com_frame.translation();
		redis_client.addEigenToWriteCallback(0, OBJECT_CARTESIAN_POSITION_KEY, object_current_pos);
		// - obtain object inertia matrix in world frame
		object_model->J_0(object_j, object_link_name, Eigen::Vector3d::Zero());
		object_inertia = (object_j*object_model->_M_inv*object_j.transpose()).inverse();
		// - obtain object p
		object_p << 0, 0, 9.8*0.5, 0, 0, 0;

        // ---------------------------------------------------------------------------------
        /* ---------------------- FILL ME IN ----------------------------------------------- */

		// --------------------------------------------------------------------
		// (1) Manipulator Jacobians
		//---------------------------------------------------------------------
		// ** FILL ME IN **

		Eigen::MatrixXd robot1_j0_objcom_world(6, robot1->dof());
		Eigen::MatrixXd robot2_j0_objcom_world(6, robot2->dof());

		// find transformation from robot ee to base
		Eigen::Affine3d T_ee1_to_base;
		Eigen::Affine3d T_ee2_to_base;

		robot1->transform(T_ee1_to_base, ee_link_name);
		robot2->transform(T_ee2_to_base, ee_link_name);

		// find transformation from ee to COM
		Eigen::Affine3d T_com_to_ee1 = (T_ee1_to_base.inverse() * robot1_base_frame.inverse() * object_com_frame);
		Eigen::Affine3d T_com_to_ee2 = (T_ee2_to_base.inverse() * robot2_base_frame.inverse() * object_com_frame);

		// Eigen::Affine3d T_ee1_to_com = T_base1_to_ee * robot1_base_frame * object_com_frame.inverse();
		// Eigen::Affine3d T_ee2_to_com = T_base2_to_ee * robot2_base_frame * object_com_frame.inverse();
		
		// find positon of COM in ee frame
		Eigen::Vector3d pos_com_ee1 = T_com_to_ee1.translation();
		Eigen::Vector3d pos_com_ee2 = T_com_to_ee2.translation();


		robot1->J_0(robot1_j0_objcom_world, ee_link_name, pos_com_ee1);
		robot2->J_0(robot2_j0_objcom_world, ee_link_name, pos_com_ee2);

		// --------------------------------------------------------------------
		// (2) Augmented Object Model
		//---------------------------------------------------------------------
		// ** FILL ME IN **

		//update joint gravity vector
		VectorXd g1(robot1->dof()); //joint space gravity vector
		VectorXd g2(robot2->dof()); //joint space gravity vector

        robot1 ->gravityVector(g1);
		robot2 ->gravityVector(g2);

		//update basic jacobian of the robot
		Eigen::MatrixXd robot1_j0(6, robot1->dof());
		Eigen::MatrixXd robot2_j0(6, robot2->dof());
        // robot1->J_0(robot1_j0,ee_link_name);
		// robot2->J_0(robot2_j0,ee_link_name);


        //update the mass matrix
        MatrixXd A1 = robot1->_M;
		MatrixXd A2 = robot2->_M;

		// find kinetic energy matrix 6 x 6
        MatrixXd L_hat1 = (robot1_j0_objcom_world * A1.inverse() * robot1_j0_objcom_world.transpose()).inverse();
		MatrixXd L_hat2 = (robot2_j0_objcom_world * A2.inverse() * robot2_j0_objcom_world.transpose()).inverse();

		// find dynamically consistent pseudo inverse of j0, j0_inv 7 x 6 
        // MatrixXd robot1_j0_inv = robot1_j0.transpose() * (robot1_j0 * robot1_j0.transpose()).inverse();
        // MatrixXd robot2_j0_inv = robot2_j0.transpose() * (robot2_j0 * robot2_j0.transpose()).inverse();
		MatrixXd robot1_jbar = A1.inverse() * robot1_j0_objcom_world.transpose() * L_hat1;
		MatrixXd robot2_jbar = A2.inverse() * robot2_j0_objcom_world.transpose() * L_hat2;

		VectorXd px1 = robot1_jbar.transpose() * g1;
		VectorXd px2 = robot2_jbar.transpose() * g2;

		MatrixXd augumented_object_inertia = object_inertia + L_hat1 + L_hat2;

		VectorXd augumented_object_p = object_p + px1 + px2;


		// --------------------------------------------------------------------
		// (3) Grasp Matrix
		//---------------------------------------------------------------------
		// ** FILL ME IN **

		//use equation 9.17

		MatrixXd G(12, 12);

		MatrixXd W_f(6, 6);
		MatrixXd W_m(6, 6);
		MatrixXd E_inv(1, 6);
		MatrixXd I(5, 6);
		MatrixXd zero(5, 6);

		MatrixXd eye3(3,3); 
        eye3.setIdentity();

		//find r1 and r2 operator		
		// Eigen::Vector3d r1 = -T_ee1_to_com.linear().transpose() * T_ee1_to_com.translation();
		// Eigen::Vector3d r2 = -T_ee2_to_com.linear().transpose() * T_ee2_to_com.translation();
		Eigen::Vector3d r1 = T_com_to_ee1.inverse().translation();
		Eigen::Vector3d r2 = T_com_to_ee2.inverse().translation();

		// cout << r1 << endl;
		// cout << r2 << endl;
		
		Eigen::Matrix3d r1_op = getCrossProductMat(r1);
		Eigen::Matrix3d r2_op = getCrossProductMat(r2);

        W_f.topLeftCorner(3,3) = eye3;
		W_f.topRightCorner(3,3) = eye3;  
		W_f.bottomLeftCorner(3,3) = r1_op;
        W_f.bottomRightCorner(3,3) = r2_op;

		W_m.bottomLeftCorner(3,3) = eye3;
        W_m.bottomRightCorner(3,3) = eye3;

		E_inv << -0.5,0,0,0.5,0,0;

		I << -0.5,0,0,0.5,0,0,
			 0,1,0,0,0,0,
			 0,0,1,0,0,0,
			 0,0,0,0,1,0,
			 0,0,0,0,0,1;

		zero << 0,0,0,0,0,0,
			0,0,0,0,0,0,
			0,0,0,0,0,0,
			0,0,0,0,0,0,
			0,0,0,0,0,0;
		
		

		G.topLeftCorner(6,6) = W_f;
		G.topRightCorner(6,6) = W_m;  
        G.bottomRightCorner(5,6) = I;
		G.bottomLeftCorner(5,6) = zero;

		G.row(6) << E_inv, 0,0,0,0,0,0;

		// --------------------------------------------------------------------
		// (4) Force Control 
		//---------------------------------------------------------------------
		if (control_mode == CONTROL_AUGMENTED_OBJECT) {
		
		// ** FILL ME IN ** Compute tau1 and tau2 

 		//================= Compute F_motion===============================

		//current object position
		object_model->transform(object_com_frame, object_link_name);
		object_current_pos << object_com_frame.translation();

		//desired object position
		Eigen::Vector3d object_desired_pos;
		object_desired_pos << 0, 0.15 * sin(2*M_PI*curr_time/3), 0.4;

        //update quaterniond repr of the rotation
		Matrix3d object_rot_mat; //end effector rotation
		//object_rot_mat << object_com_frame.linear();
		object_model->rotation(object_rot_mat, object_link_name);
		Eigen::Quaterniond object_rot_lambda(object_rot_mat);

		// find desired_orientation in quaterniond
        Quaterniond object_desired_lambda;
        object_desired_lambda.w() = 1;
        object_desired_lambda.x() = 0;
        object_desired_lambda.y() = 0;
        object_desired_lambda.z() = 0;

		// solve for instantaneous angular error
        Vector3d delta_phi;
        Sai2Model::orientationError(delta_phi, object_desired_lambda, object_rot_lambda);


		double lamb0d = object_desired_lambda.w();
        double lamb1d = object_desired_lambda.x();
        double lamb2d = object_desired_lambda.y();
        double lamb3d = object_desired_lambda.z();

		double lamb0d_d = 0;
        double lamb1d_d = 0;
        double lamb2d_d = 0;
        double lamb3d_d = 0;

		double lamb0d_dd = 0;
        double lamb1d_dd = 0;
        double lamb2d_dd = 0;
        double lamb3d_dd = 0;

		// find desired position vector x_pd
        VectorXd obj_pd(3);
        obj_pd << 0, 0.15 * sin(2*M_PI*curr_time/3), 0.4;

        // find 1st derivative of desired position vector x_pd_d
        VectorXd obj_pd_d(3);
        obj_pd_d << 0, M_PI/10 * cos(2*M_PI*curr_time/3), 0;

        // find 2nd derivative of desired position vector x_pd_d
        VectorXd obj_pd_dd(3);
        obj_pd_dd << 0, -2* pow(M_PI,2)/30 * sin(2*M_PI*curr_time/3), 0;

		
		//update velocity in 6 axis (3 linear and 3 rotational)
		VectorXd v0(6); //object velocity
        object_model-> velocity6d(v0, object_link_name);

		cout << v0 << endl;

		//find current position vector x_p and linear velcity vector x_p_d
        MatrixXd obj_p = object_current_pos;
        MatrixXd obj_p_d = v0(seq(0,2));

        //determine force input of decouple end effector F*
        MatrixXd F_s = obj_pd_dd - kp * (obj_p - obj_pd) - kv * (obj_p_d - obj_pd_d);
        
        //left inverse of desired rotation x_rd, equation 2.45
        MatrixXd Er_inv(3,4);
        Er_inv(0,0) = -lamb1d, Er_inv(0,1)= lamb0d, Er_inv(0,2) = -lamb3d, Er_inv(0,3) = lamb2d;
        Er_inv(1,0) = -lamb2d, Er_inv(1,1)= lamb3d, Er_inv(1,2) = lamb0d, Er_inv(1,3) = -lamb1d;
        Er_inv(2,0) = -lamb3d, Er_inv(2,1)= -lamb2d, Er_inv(2,2) = lamb1d, Er_inv(2,3) = lamb0d;

        Er_inv *= 2;

		VectorXd x_rd_d(4);
        x_rd_d << lamb0d_d, lamb1d_d, lamb2d_d, lamb3d_d;
        
        MatrixXd omega_des = Er_inv * x_rd_d;

        VectorXd lamb_des_dd(4);
        lamb_des_dd << lamb0d_dd, lamb1d_dd, lamb2d_dd, lamb3d_dd;
        
        //find desired angular acceleration omega_des_d and current 
        //angular velocity omega , eqn 4.68 
        MatrixXd omega_des_d = Er_inv * lamb_des_dd;
        
        //MatrixXd omega1 = v0(seq(3, last));
        Vector3d omega;
        object_model->angularVelocity(omega, object_link_name);

        //determine moment input of decouple end effector M*
        MatrixXd M_s = omega_des_d - kp * delta_phi - kv * (omega - omega_des);

		VectorXd input(F_s.size() + M_s.size());
        input << F_s, M_s;

		Vector3d tension_int;
		tension_int << 0,-15,0;
		Vector3d moment_int;
		moment_int << 0,0,0;
		
		VectorXd force_result;
		force_result = augumented_object_inertia * input + object_p;


		VectorXd force_des_vec(force_result.size() + 3 + 3);
		force_des_vec << force_result, tension_int, moment_int;

		VectorXd force_ee_vec = G.inverse() * force_des_vec;

		VectorXd f1_ee = force_ee_vec(seq(0,2));
		VectorXd f2_ee = force_ee_vec(seq(3,5));
		VectorXd m1_ee = force_ee_vec(seq(6,8));
		VectorXd m2_ee = force_ee_vec(seq(9,11));

		VectorXd world_F1_ee(6);
		world_F1_ee << f1_ee, m1_ee;
		VectorXd world_F2_ee(6);
		world_F2_ee << f2_ee, m2_ee;

		// cout << "===========" << endl;

		// cout << "ee1 force" << endl;
		// cout << f1_ee.transpose() << endl;
		// cout << "ee2 force" << endl;
		// cout << f2_ee.transpose() << endl;
		// cout << object_current_pos<< endl;

		Eigen::MatrixXd robot1_j0_ee_world(6, robot1->dof());
		Eigen::MatrixXd robot2_j0_ee_world(6, robot2->dof());

		robot1->J_0WorldFrame(robot1_j0_ee_world, ee_link_name);
		robot2->J_0WorldFrame(robot2_j0_ee_world, ee_link_name);
		
		tau1 = robot1_j0_ee_world.transpose() * world_F1_ee + g1;
		tau2 = robot2_j0_ee_world.transpose() * world_F2_ee + g2;


		// count ++;
		


		} else if (control_mode == CONTROL_GRASP_STABILIZE) { // initial grasp stabilization
			
			Eigen::MatrixXd robot1_j0_ee(6, robot1->dof());
			Eigen::MatrixXd robot2_j0_ee(6, robot2->dof());
			robot1->J_0(robot1_j0_ee, ee_link_name, Eigen::Vector3d::Zero());
			robot2->J_0(robot2_j0_ee, ee_link_name, Eigen::Vector3d::Zero());

			// Joint Torques
			tau1 = robot1_j0_ee.transpose()*(object_p/2) + robot1->_M*(-10.0*robot1->_dq) + robot1_g;
			tau2 = robot2_j0_ee.transpose()*(object_p/2) + robot2->_M*(-10.0*robot2->_dq) + robot2_g;
			
			// Grasp stabilization
			static uint grasp1Counter = 0;
			static uint grasp2Counter = 0;
			if (robot1->_dq[6] < 0.1) {
				grasp1Counter += 1;
			} else {
				grasp1Counter = 0;
			}
			if (robot2->_dq[6] < 0.1) {
				grasp2Counter += 1;
			} else {
				grasp2Counter = 0;
			}
			if (grasp1Counter > 40 && grasp2Counter > 40) {
				cout << " ** Switch Control Mode to Augmented Object Model ** " << endl;
				control_mode = CONTROL_AUGMENTED_OBJECT;
			}
		}

		/* ----------------------------------------------------------------------------------
			Safety torques 
		-------------------------------------------------------------------------------------*/ 
		
		// Set constant gripper forces
		tau1[6] = 15;
		tau2[6] = 15;

        // Default values if torques are exceeded:
        bool fTorqueOverride = false; // to avoid robot blow-ups
        const double tau1_max = 200;
        const double tau2_max = 200;
        if (!fTorqueUseDefaults) {
        	if (tau1.cwiseAbs().maxCoeff() > tau1_max || tau2.cwiseAbs().maxCoeff() > tau2_max) {
	        	fTorqueOverride = true;
	        	cerr << "Torque overriden. User asked torques beyond safe limit: \n";
	        	cerr << "Robot 1: " << tau1.transpose() << "\n";
	        	cerr << "Robot 2: " << tau2.transpose() << "\n";
	        	fTorqueUseDefaults = true;
	        }
	        // Also default values if object is dropped
	        const double object_thickness = 0.05;
	        bool fRobot1DroppedObject = robot1->_q[6] > object_thickness/2;
	        bool fRobot2DroppedObject = robot2->_q[6] > object_thickness/2;
	        if (fRobot1DroppedObject || fRobot2DroppedObject) {
	        	cerr << "Torque overriden. Robot 1 or 2 dropped object. \n";
	        	fTorqueUseDefaults = true;
	        }
        }
        else {
        	robot1->gravityVector(tau1);
			tau1 = tau1 + robot1->_M*(-10.0*robot1->_dq);
			tau1 = (tau1.array() >= tau1_max).select(tau1_max*Eigen::VectorXd::Ones(robot1->dof()), tau1);
			tau1 = (tau1.array() <= -tau1_max).select(-tau1_max*Eigen::VectorXd::Ones(robot1->dof()), tau1);
			robot2->gravityVector(tau2);
			tau2 = tau2 + robot2->_M*(-10.0*robot2->_dq);
			tau2 = (tau2.array() >= tau2_max).select(tau2_max*Eigen::VectorXd::Ones(robot2->dof()), tau2);
			tau2 = (tau2.array() <= -tau2_max).select(-tau2_max*Eigen::VectorXd::Ones(robot2->dof()), tau2);
        }

        /* ----------------------------------------------------------------------------------
			Send robot torques from controller
		-------------------------------------------------------------------------------------*/ 
		sim->setJointTorques(robot1_name, tau1);
		sim->setJointTorques(robot2_name, tau2);
		
		// write object location key out to redis
		redis_client.executeWriteCallback(0);

		// update last time
		last_time = curr_time;

		// if (count > 5) { break;}
	}
}

/* =======================================================================================
   SIMULATION SETUP
   -----------------------
   * Simulation loop
   * Select trajectory
   * Window initialization
   * Window error
   * Mouse click commands
========================================================================================== */
void simulation(Simulation::Sai2Simulation* sim) {
	fSimulationRunning = true;

	// create a timer
	LoopTimer timer;
	timer.initializeTimer();
	timer.setLoopFrequency(10000); //10000Hz timer

	double last_time = timer.elapsedTime()/time_slowdown_factor; //secs
	bool fTimerDidSleep = true;
	while (fSimulationRunning) {
		fTimerDidSleep = timer.waitForNextLoop();
		// integrate forward
		double curr_time = timer.elapsedTime()/time_slowdown_factor;
		double loop_dt = curr_time - last_time; 
		// sim->integrate(loop_dt);
		sim->integrate(loop_dt);
		// update last time
		last_time = curr_time;
	}
}

//------------------------------------------------------------------------------

void glfwError(int error, const char* description) {
	cerr << "GLFW Error: " << description << endl;
	exit(1);
}

//------------------------------------------------------------------------------

void keySelect(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    // option ESC: exit
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    {
        // exit application
         glfwSetWindowShouldClose(window, 1);
    }
}
