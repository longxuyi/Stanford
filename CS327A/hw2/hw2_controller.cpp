#include <Sai2Model.h>
#include "redis/RedisClient.h"
#include "timer/LoopTimer.h"

#include <iostream>
#include <string>
#include <cmath>

using namespace std;
using namespace Eigen;

#include <signal.h>
bool runloop = false;
void sighandler(int){runloop = false;}


// Convenience functions for converting degrees to radians and vice-versa
#define DEG2RAD(deg) ((double)(deg) * M_PI / 180.0)
#define RAD2DEG(rad) ((double)(rad) * 180.0 / M_PI)

// Location of URDF files specifying world and robot information
const string world_file = "./resources/world_hw2.urdf";
const string robot_file = "../../resources/kuka_iiwa/kuka_iiwa.urdf";
const string robot_name = "Kuka-IIWA"; 
const string ee_link_name = "link6";

// Redis is just a key value store, publish/subscribe is also possible
// The visualizer and simulator will have keys like "cs225a::robot::{ROBOTNAME}::sensors::q"
// You can hardcode the robot name in like below or read them in from cli
// redis keys:
// - read (from simulation/robot):
const std::string JOINT_ANGLES_KEY  = "cs327a::robot::Kuka-IIWA::sensors::q";
const std::string JOINT_VELOCITIES_KEY = "cs327a::robot::Kuka-IIWA::sensors::dq";
// - read (from user defined task space control):
const std::string TASK_KP_GAINS_KEY = "cs327a::robot::Kuka-IIWA::inputs::op_task_kp";
const std::string TASK_KV_GAINS_KEY = "cs327a::robot::Kuka-IIWA::inputs::op_task_kv";

// - write:
const std::string JOINT_TORQUES_COMMANDED_KEY = "cs327a::robot::Kuka-IIWA::actuators::fgc";
const std::string PLOT_CURRENT_EE_POSITION_KEY = "cs327a::plot::Kuka-IIWA::sensors::ee_pos";
const std::string PLOT_DESIRED_EE_POSITION_KEY = "cs327a::plot::Kuka-IIWA::inputs::ee_pos_des";

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

    // Suggested controller variables
    int dof = robot->dof();
    VectorXd command_torques = VectorXd::Zero(dof);
    VectorXd g(robot->dof()); //joint space gravity vector
    MatrixXd J0(6, robot->dof()); //end effector Basic Jacobian,
    MatrixXd L_hat(6, 6); //Kinetic Energy Matrix at end effector (op. space)
    VectorXd p_hat(6); //gravity vector at end-effector (op. space)
    const Eigen::MatrixXd In = Eigen::MatrixXd::Identity(robot->dof(), robot->dof()); // n x n identity matrix
    Eigen::MatrixXd J_pseudo(robot->dof(), 6);  //Right Psuedoinverse of J0
    Eigen::MatrixXd N_proj(robot->dof(), robot->dof()); //I - J_pseudo*J0, null space projection matrix
    Vector3d ee_pos; //end effector position
    Matrix3d ee_rot_mat; //end effector rotation, in matrix form
    robot->rotation(ee_rot_mat, ee_link_name); // initialize the end effector rotation
    Quaterniond ee_rot_lambda(ee_rot_mat); // end effector rotation in quaternion form

    VectorXd v0(6); //end effector velocity
    Vector3d ee_pos_des; // end effector desired position
    VectorXd v0_des(6); //end effector desired velocity
    VectorXd dv0_des(6); //end effector desired acceleration
    VectorXd ee_error(6); //end effector operational space instantaneous error (position and orientation)

    // suggested starting gains
    double op_task_kp = 50; //operational space task proportional gain
    double op_task_kv = 20; //operational space task velocity gain
    double damping_kv = 20; // joint damping velocity gain

    // initialize redis keys

    // op task
    // redis_client.setEigenMatrixJSON(DESIRED_TASK_POSITION_KEY, desired_position);
    redis_client.set(TASK_KP_GAINS_KEY, std::to_string(op_task_kp));
    redis_client.set(TASK_KV_GAINS_KEY, std::to_string(op_task_kv));

    // prepare redis read callback
    redis_client.addEigenToReadCallback(0, JOINT_ANGLES_KEY, robot->_q);
    redis_client.addEigenToReadCallback(0, JOINT_VELOCITIES_KEY, robot->_dq);
    redis_client.addDoubleToReadCallback(0, TASK_KP_GAINS_KEY, op_task_kp);
    redis_client.addDoubleToReadCallback(0, TASK_KV_GAINS_KEY, op_task_kv);

    // prepare redis write callback
    redis_client.addEigenToWriteCallback(0, JOINT_TORQUES_COMMANDED_KEY, command_torques);
    redis_client.addEigenToWriteCallback(0, PLOT_CURRENT_EE_POSITION_KEY, ee_pos);
    redis_client.addEigenToWriteCallback(0, PLOT_DESIRED_EE_POSITION_KEY, ee_pos_des);

	// create a loop timer
	double control_freq = 1000;
	LoopTimer timer;
	timer.setLoopFrequency(control_freq);   // 1 KHz
	double last_time = timer.elapsedTime(); //secs
	bool fTimerDidSleep = true;
	timer.initializeTimer(1000000); // 1 ms pause before starting loop
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
        
        //curr_time = 1;


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
        

        // //wesley's code
        // double R = 0.1; //in meters
        // double T = 5.0; //in seconds

        // double lamb0d_dd = 1/sqrt(2)*(-sin(M_PI/4.0*cos(2.0*M_PI*curr_time/T)) * pow(-pow(M_PI, 2)/(2.0*T)*sin(2.0*M_PI*curr_time/T), 2) - cos(M_PI/4.0*cos(2.0*M_PI*curr_time/T)) * pow(M_PI, 3)/(T*T)*cos(2.0*M_PI*curr_time/T));    

        // double lamb1d_dd = 1/sqrt(2)*(-cos(M_PI/4.0*cos(2.0*M_PI*curr_time/T)) * pow(-pow(M_PI, 2)/(2.0*T)*sin(2.0*M_PI*curr_time/T), 2) + sin(M_PI/4.0*cos(2.0*M_PI*curr_time/T)) * pow(M_PI, 3)/(T*T)*cos(2.0*M_PI*curr_time/T));

        // double lamb2d_dd = 1/sqrt(2)*(-sin(M_PI/4.0*cos(2.0*M_PI*curr_time/T)) * pow(-pow(M_PI, 2)/(2.0*T)*sin(2.0*M_PI*curr_time/T), 2) - cos(M_PI/4.0*cos(2.0*M_PI*curr_time/T)) * pow(M_PI, 3)/(T*T)*cos(2.0*M_PI*curr_time/T));    

        // double lamb3d_dd = 1/sqrt(2)*(-cos(M_PI/4.0*cos(2.0*M_PI*curr_time/T)) * pow(-pow(M_PI, 2)/(2.0*T)*sin(2.0*M_PI*curr_time/T), 2) + sin(M_PI/4.0*cos(2.0*M_PI*curr_time/T)) * pow(M_PI, 3)/(T*T)*cos(2.0*M_PI*curr_time/T));

        
        // solve for instantaneous angular error
        Vector3d delta_phi;
        Sai2Model::orientationError(delta_phi, desired_lambda, ee_rot_lambda);
    
        // use equations in course reader page80 to solve for input of the 
        //decoupled ee F* and M* for trajectory tracking tasks
        
        // find desired position vector x_pd
        VectorXd x_pd(3);
        x_pd << 0, 0.5 + 0.1 * cos(2*M_PI*curr_time/5), 0.65 - 0.05 * cos(4*M_PI*curr_time/5);

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

        //std::cout << Er_inv << std::endl;        

        VectorXd x_rd_d(4);
        x_rd_d << lamb0d_d, lamb1d_d, lamb2d_d, lamb3d_d;

        //std::cout << x_rd_d << std::endl;


        
        MatrixXd omega_des = Er_inv * x_rd_d;

        VectorXd lamb_des_dd(4);
        lamb_des_dd << lamb0d_dd, lamb1d_dd, lamb2d_dd, lamb3d_dd;

        //std::cout << lamb_des_dd << std::endl;


        
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
        // MatrixXd J0_inv(7,6);
        // robot->dynConsistentInverseJacobian(J0_inv, J0);


        // find kinetic energy matrix 6 x 6
        MatrixXd kinetic_mat = J0_inv.transpose() * A * J0_inv;
        // MatrixXd kinetic_mat(6,6);
        // robot->taskInertiaMatrixWithPseudoInv(kinetic_mat, J0);

        //find joint damping
        robot->nullspaceMatrix(N_proj, J0);
        MatrixXd joint_damping = N_proj.transpose() * (-1 * robot->_M * robot->_dq * damping_kv + g);

        //calculate vector of gravity forces (operational space gravity vector) page 63
        VectorXd px = J0_inv.transpose() * g;
     
        VectorXd input(F_s.size() + M_s.size());
        input << F_s, M_s;


        // std::cout << "x_p_d - x_pd_d" <<std::endl;
        // std::cout << x_p_d - x_pd_d << std::endl;
        // std::cout << "x_pd_d" <<std::endl;
        // std::cout << x_pd_d << std::endl;
        std::cout << "x_p_d" <<std::endl;
        std::cout << x_p_d << std::endl;
        // std::cout << "x_pd_d" <<std::endl;
        // std::cout << x_pd_d << std::endl;
        // std::cout << "x_pd_dd" <<std::endl;
        // std::cout << x_pd_dd << std::endl;

        // std::cout << "delta_phi" <<std::endl;
        // std::cout << delta_phi << std::endl;
        // std::cout << "omega_des_d" <<std::endl;
        // std::cout << omega_des_d << std::endl;
        // std::cout << "omega" <<std::endl;
        // std::cout << omega << std::endl;
        // std::cout << "omega_des" <<std::endl;
        // std::cout << omega_des << std::endl;




        //define 7 x 7 identity matrix
        MatrixXd eye6(6,6); 
        eye6.setIdentity();

        VectorXd F0 = eye6 * input + px;
        //eye6 * input 

        // std::cout << input << std::endl;

        // std::cout << eye6 <<std::endl;

        // std::cout << F0 <<std::endl;

        // ---------------------------------------------------------------------------------
        // (3) Compute joint torques
        //----------------------------------------------------------------------------------

        //command_torques = J0.transpose() * F0;
         command_torques = J0.transpose() * F0 + joint_damping;

        //J0.transpose()

        /* ------------------------------------------------------------------------------------
            END OF FILL ME IN
        -------------------------------------------------------------------------------------*/

		// send torques to redis
        redis_client.executeWriteCallback(0);

		counter++;
		last_time = curr_time;

        //break;

        // if (counter > 50){
        //     break;
        // } 

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
