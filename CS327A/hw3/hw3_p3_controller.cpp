#include <Sai2Model.h>
#include "redis/RedisClient.h"
#include "timer/LoopTimer.h"

#include <iostream>
#include <string>

using namespace std;
using namespace Eigen;

#include <signal.h>
bool runloop = false;
void sighandler(int){runloop = false;}


// Convenience functions for converting degrees to radians and vice-versa
#define DEG2RAD(deg) ((double)(deg) * M_PI / 180.0)
#define RAD2DEG(rad) ((double)(rad) * 180.0 / M_PI)

// Location of URDF files specifying world and robot information
const string world_file = "./resources/world_hw3_p3.urdf";
const string robot_file = "../../resources/kuka_iiwa/kuka_iiwa.urdf";
const string robot_name = "Kuka-IIWA"; 
const string ee_link_name = "link6";

// Redis is just a key value store, publish/subscribe is also possible
// The visualizer and simulator will have keys like "cs225a::robot::{ROBOTNAME}::sensors::q"
// You can hardcode the robot name in like below or read them in from cli
// redis keys:
// - read (from simulation/robot):
const std::string JOINT_ANGLES_KEY  = "cs327a::hw3::robot::Kuka-IIWA::sensors::q";
const std::string JOINT_VELOCITIES_KEY = "cs327a::hw3::robot::Kuka-IIWA::sensors::dq";

// - write:
const std::string JOINT_TORQUES_COMMANDED_KEY  = "cs327a::hw3::robot::Kuka-IIWA::actuators::fgc";

int main(int argc, char** argv) {
	// Make sure redis-server is running at localhost with default port 6379
	// start redis client
	RedisClient redis_client = RedisClient();
	redis_client.connect();

    // set up redis read + write callbacks
    redis_client.createReadCallback(0);
    redis_client.createWriteCallback(0);

	// set up signal handler
	signal(SIGABRT, &sighandler);
	signal(SIGTERM, &sighandler);
	signal(SIGINT, &sighandler);

    // load robots, read current state and update the model
	auto robot = new Sai2Model::Sai2Model(robot_file, false);
	robot->_q = redis_client.getEigenMatrixJSON(JOINT_ANGLES_KEY);
	robot->_dq = redis_client.getEigenMatrixJSON(JOINT_VELOCITIES_KEY);
	robot->updateModel();

    // end effector tip position in link frame
    Eigen::Vector3d ee_pos_local;
    ee_pos_local << 0.0, 0.0, 0.0;

    // final control torques to send to the robot
    VectorXd command_torques = VectorXd::Zero(robot->dof());

    // Suggested controller variables
    Eigen::VectorXd g(robot->dof()); //joint space gravity vector
    Eigen::MatrixXd J0(6, robot->dof()); //end effector basic Jacobian
    Eigen::MatrixXd L0(6, 6); //Lambda_0 at end effector
    Eigen::VectorXd p(6); //gravity vector at end-effector
    const Eigen::MatrixXd In = Eigen::MatrixXd::Identity(robot->dof(), robot->dof()); // n x n identity matrix
    Eigen::MatrixXd J_bar(robot->dof(), 6);  //A-weighted generalized inverse of J0
    Eigen::MatrixXd N_bar(robot->dof(), robot->dof()); //I - Jbar*J0, null space projection matrix for Jbar
    Eigen::Vector3d ee_pos; //end effector position
    Eigen::Matrix3d ee_rot_mat; //end effector rotation
    robot->rotation(ee_rot_mat, ee_link_name); // initialize
    Eigen::Quaterniond ee_rot_lambda(ee_rot_mat); // end effector rotation in quaternion form
    Eigen::VectorXd ee_error(6); //end effector operational space instantaneous error
    Eigen::VectorXd v0(6); //end effector velocity
    Eigen::VectorXd v0d(6); //end effector desired velocity
    Eigen::VectorXd dv0d(6); //end effector desired acceleration
    Eigen::MatrixXd select_motion(6,6); // selection matrix for motion, Omega
    Eigen::MatrixXd select_forces(6,6); // selection matrix for forces, Omega_bar
    Eigen::VectorXd F_d(6); // desired end effector force
    Eigen::MatrixXd N_proj(robot->dof(), robot->dof()); //I - J_pseudo*J0, null space projection matrix

    // suggested starting gains
    double op_task_kp = 50; //operational space task proportional gain
    double op_task_kv = 20; //operational space task velocity gain
    double joint_damping_kv = 20; // joint damping velocity gain

    // initialize redis keys

    // prepare redis read callback
    redis_client.addEigenToReadCallback(0, JOINT_ANGLES_KEY, robot->_q);
    redis_client.addEigenToReadCallback(0, JOINT_VELOCITIES_KEY, robot->_dq);

    // prepare redis write callback
    redis_client.addEigenToWriteCallback(0, JOINT_TORQUES_COMMANDED_KEY, command_torques);

	// create a loop timer
    double control_freq = 200;
    LoopTimer timer;
    timer.setLoopFrequency(control_freq);
    double last_time = timer.elapsedTime(); //secs
    bool fTimerDidSleep = true;
    timer.initializeTimer(2000);
    double start_time = timer.elapsedTime();

	unsigned long long counter = 0;

    cout << "Starting control loop" << endl;

	runloop = true;
	while (runloop) 
	{ 
		fTimerDidSleep = timer.waitForNextLoop();
        double curr_time = timer.elapsedTime() - start_time;
        double loop_dt = curr_time - last_time;

        // read controller parameters from redis
        redis_client.executeReadCallback(0);

        // update robot model
        robot->updateModel();
        robot->gravityVector(g);

        // reset command torques
        command_torques.setZero();

		/* ------------------------------------------------------------------------------------
            FILL ME IN: set joint torques
        -------------------------------------------------------------------------------------*/

        // ---------------------------------------------------------------------------------
        // (1) Update current robot configuration parameters and desired trajectory values 
        //----------------------------------------------------------------------------------
        //update ee position of the ee
        robot->position(ee_pos, ee_link_name);

        //update quaterniond repr of the rotation
        robot->rotation(ee_rot_mat, ee_link_name);

        Quaterniond ee_rot_lambda(ee_rot_mat);

        //update velocity in 6 axis (3 linear and 3 rotational)
        robot-> velocity6d(v0, ee_link_name);

        // std::cout << "velocity v6 vector" <<std::endl;
        // std::cout << v0 <<std::endl;

        double lamb0 = ee_rot_lambda.w();
        double lamb1 = ee_rot_lambda.x();
        double lamb2 = ee_rot_lambda.y();
        double lamb3 = ee_rot_lambda.z();

        //update joint gravity vector
        robot ->gravityVector(g);

        //update basic jacobian of the robot
        robot->J_0(J0,ee_link_name);

        //update the mass matrix
        MatrixXd A = robot->_M; // 7 x 7 matrix  

        //define J bar
        J_bar = A.inverse() * J0.transpose() * (J0 * A.inverse() * J0.transpose()).inverse();     

        
        // --------------------------------------------------------------------
        // (2) Compute desired operational space trajectory values and errors 
        //---------------------------------------------------------------------

        // find desired_orientation in quaterniond
        Quaterniond desired_lambda;
        desired_lambda.w() = 1/sqrt(2) *sin(M_PI/4 * cos(2*M_PI*curr_time/5));
        desired_lambda.x() = 1/sqrt(2) *cos(M_PI/4 * cos(2*M_PI*curr_time/5));
        desired_lambda.y() = 1/sqrt(2) *sin(M_PI/4 * cos(2*M_PI*curr_time/5));
        desired_lambda.z() = 1/sqrt(2) *cos(M_PI/4 * cos(2*M_PI*curr_time/5));

        double lamb0d = desired_lambda.w();
        double lamb1d = desired_lambda.x();
        double lamb2d = desired_lambda.y();
        double lamb3d = desired_lambda.z();

        double lamb0d_d = -M_PI * M_PI/(10*sqrt(2)) * cos(M_PI/4 * cos(2*M_PI*curr_time/5))* sin(2*M_PI*curr_time/5);
        double lamb1d_d = M_PI * M_PI/(10*sqrt(2)) * sin(M_PI/4 * cos(2*M_PI*curr_time/5))* sin(2*M_PI*curr_time/5);       
        double lamb2d_d = -M_PI * M_PI/(10*sqrt(2)) * cos(M_PI/4 * cos(2*M_PI*curr_time/5))* sin(2*M_PI*curr_time/5);
        double lamb3d_d = M_PI * M_PI/(10*sqrt(2)) * sin(M_PI/4 * cos(2*M_PI*curr_time/5))* sin(2*M_PI*curr_time/5); 

        double lamb0d_dd = -M_PI * M_PI/(10*sqrt(2)) * (2*M_PI/5*cos(2*M_PI/5*curr_time)* cos(M_PI/4*cos(2*M_PI/5*curr_time)) + M_PI*M_PI/10 * pow(sin(2*M_PI/5*curr_time),2) * sin(M_PI/4*cos(2*M_PI/5*curr_time)));    
        double lamb1d_dd = M_PI * M_PI/(10*sqrt(2)) * (2*M_PI/5*cos(2*M_PI/5*curr_time)* sin(M_PI/4*cos(2*M_PI/5*curr_time)) - M_PI*M_PI/10 * pow(sin(2*M_PI/5*curr_time),2) * cos(M_PI/4*cos(2*M_PI/5*curr_time))); 
        double lamb2d_dd = -M_PI * M_PI/(10*sqrt(2)) * (2*M_PI/5*cos(2*M_PI/5*curr_time)* cos(M_PI/4*cos(2*M_PI/5*curr_time)) + M_PI*M_PI/10 * pow(sin(2*M_PI/5*curr_time),2) * sin(M_PI/4*cos(2*M_PI/5*curr_time)));    
        double lamb3d_dd = M_PI * M_PI/(10*sqrt(2)) * (2*M_PI/5*cos(2*M_PI/5*curr_time)* sin(M_PI/4*cos(2*M_PI/5*curr_time)) - M_PI*M_PI/10 * pow(sin(2*M_PI/5*curr_time),2) * cos(M_PI/4*cos(2*M_PI/5*curr_time))); 

        
        
        //sigma and sigma_bar

        Eigen::MatrixXd sigma_f(3,3);
        sigma_f << 0,0,0,
                   0,1,0,
                   0,0,1;  

        cout << sigma_f << endl;
        
        Eigen::MatrixXd sigma_m(3,3);
        sigma_m << 1,0,0,
                   0,0,0,
                   0,0,0; 
        
        cout << sigma_m << endl;

        MatrixXd eye3(3,3); 
        eye3.setIdentity();

        select_motion.topLeftCorner(3,3) = sigma_f ;
        select_motion.bottomRightCorner(3,3) = sigma_m;

        select_forces.topLeftCorner(3,3) = eye3 - sigma_f;
        select_forces.bottomRightCorner(3,3) = eye3 - sigma_m;

        
        
        // solve for instantaneous angular error
        Vector3d delta_phi;
        Sai2Model::orientationError(delta_phi, desired_lambda, ee_rot_lambda);
    
        // use equations in course reader page80 to solve for input of the 
        //decoupled ee F* and M* for trajectory tracking tasks
        
        // find desired position vector x_pd
        VectorXd x_pd(3);
        x_pd << 0.05, 0.5 + 0.1 * cos(2*M_PI*curr_time/5), 0.65 - 0.05 * cos(4*M_PI*curr_time/5);

        // find 1st derivative of desired position vector x_pd_d
        VectorXd x_pd_d(3);
        x_pd_d << 0, -M_PI/25 * sin(2*M_PI*curr_time/5), M_PI/25*sin(4*M_PI*curr_time/5);

        // find 2nd derivative of desired position vector x_pd_d
        VectorXd x_pd_dd(3);
        x_pd_dd << 0, -2*pow(M_PI, 2)/125 * cos(2*M_PI*curr_time/5), 4*pow(M_PI, 2)/125 * cos(4*M_PI*curr_time/5);

        //find current position vector x_p and linear velcity vector x_p_d
        MatrixXd x_p = ee_pos;
        MatrixXd x_p_d = v0(seq(0,2));

        //std::cout << x_p_d <<std::endl;
       
        //determine force input of decouple end effector F*
        MatrixXd F_s = x_pd_dd - op_task_kp * (x_p - x_pd) - op_task_kv * (x_p_d - x_pd_d);
        

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
        robot->angularVelocity(omega, ee_link_name);

       
        //determine moment input of decouple end effector M*
        MatrixXd M_s = omega_des_d - op_task_kp * delta_phi - op_task_kv * (omega - omega_des);
        //MatrixXd M_s = 2*(desired_lambda_ddot*ee_rot_lambda.conjugate()).vec() - op_task_kp * delta_phi - op_task_kv * (omega - 2*(desired_lambda_dot*ee_rot_lambda.conjugate()).vec());


        // find pseudo inverse of j0, j0_inv 7 x 6 
        MatrixXd J0_inv = J0.transpose() * (J0 * J0.transpose()).inverse(); 


        // find kinetic energy matrix 6 x 6
        MatrixXd kinetic_mat = J_bar.transpose() * A * J_bar;
        // MatrixXd kinetic_mat(6,6);


        //define 6 x 6 identity matrix
        MatrixXd eye7(7,7); 
        eye7.setIdentity();

        //find nullspace with j_bar
        N_proj = eye7 - J0.transpose() * J_bar.transpose();

        cout << "==================" << endl;
        cout << N_proj << endl;


        //find joint damping
        MatrixXd joint_damping = N_proj * (-1 * robot->_M * robot->_dq * joint_damping_kv + g);

        //calculate vector of gravity forces (operational space gravity vector) page 63
        VectorXd px = J_bar.transpose() * g;
     
        VectorXd input(F_s.size() + M_s.size());
        input << F_s, M_s;

        VectorXd f_motion = kinetic_mat * select_motion * input + px;

        F_d << 10,0,0,0,0,0;

        VectorXd f_force = /*kinetic_mat * */ select_forces * F_d;


        //operational space damping term
        VectorXd operation_damping = kinetic_mat * select_forces * v0 * -op_task_kv;

        // cout << kinetic_mat.size() << endl;

        VectorXd F0 = f_motion + f_force + operation_damping;

        // ---------------------------------------------------------------------------------
        // (3) Compute joint torques
        //----------------------------------------------------------------------------------

       command_torques = J0.transpose() * F0 + joint_damping;
		/* ------------------------------------------------------------------------------------
            END OF FILL ME IN
        -------------------------------------------------------------------------------------*/

        // send torques to redis
        redis_client.executeWriteCallback(0);

		counter++;
		last_time = curr_time;
	}

    command_torques.setZero();
    redis_client.setEigenMatrixJSON(JOINT_TORQUES_COMMANDED_KEY, command_torques);

	double end_time = timer.elapsedTime();
    cout << "\n";
    cout << "Control Loop run time  : " << end_time << " seconds\n";
    cout << "Control Loop updates   : " << timer.elapsedCycles() << "\n";
    cout << "Control Loop frequency : " << timer.elapsedCycles()/end_time << "Hz\n";

    return 0;
}


