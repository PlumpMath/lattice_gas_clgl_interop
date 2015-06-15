//include the opencl header
#ifdef __APPLE__
#include "OpenCL/opencl.h"
#else
//using AMD APP SDK, but opencl 1.1 headers
//because theres no support for opencl 1.2 for my GT630 card
#include "cl.h"
#include "cl_gl.h"
#include <GL/glew.h>
#include "GL/glut.h"
#endif

//SFML for window
#include <SFML/System.hpp>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <GL/glew.h>

//now its only for windows actually
#include<windows.h>

//SHARING EXTENSION Defint
#if defined (__APPLE__) || defined(MACOSX)
static const char* CL_GL_SHARING_EXT = "cl_APPLE_gl_sharing";
#else
static const char* CL_GL_SHARING_EXT = "cl_khr_gl_sharing";
#endif


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <initializer_list>
#include <stdlib.h>
#include <stdio.h>
#include <random>
#include<thread>
#include <string.h>
#include<math.h>



/*/
Description of states
- Every cell can have up to 4 particles, each of them corresponding to a direction(up, left, down, right)
- A state of a pixel is 4 length binary number, which represents the number of particles in the pixel : (up, left, down, right), eg. : a pixel with one partcle going to the left is(0, 1, 0, 0)
- I will represent the pixel with the integer corresponding to the binary series.eg. : (0, 1, 0, 0) is 4
*/

//lookup table to transform int to state
	//state[0] particle going up
	//state[1] particle going right
	//state[2] particle going down
	//state[3] particle going left

int int_2_state[16][4] =
{ { 0, 0, 0, 0 },
{ 0, 0, 0, 1 },
{ 0, 0, 1, 0 },
{ 0, 0, 1, 1 },
{ 0, 1, 0, 0 },
{ 0, 1, 0, 1 },
{ 0, 1, 1, 0 },
{ 0, 1, 1, 1 },
{ 1, 0, 0, 0 },
{ 1, 0, 0, 1 },
{ 1, 0, 1, 0 },
{ 1, 0, 1, 1 },
{ 1, 1, 0, 0 },
{ 1, 1, 0, 1 },
{ 1, 1, 1, 0 },
{ 1, 1, 1, 1 } };

/*
	Lattice gas class for cpu sumulation
*/
class Lattice_gas{

public:
	Lattice_gas(int N) : states_vec(2, std::vector<int>((N)*(N))), temp(N*N) {
		this->N = N;
		this->N_small = N - 2; //outer layer is just for convienience
		iterno = 0;
	}
	~Lattice_gas(){}

private:
	int N; //size of the table
	int N_small; // interesting part size of the table

	//states in cells, represented by integer
	//I do synchronized update so 2 copies needed
	std::vector< std::vector<int> > states_vec; 
	//"image" used to write into texture
	std::vector<uint8_t> temp; 
	int iterno; //(using alternating vectors) 1.iter: vec1 is input-> vec2 is output
	//2.iter: vec2 is input->vec1 is output

	//helper function transfrom state to int
	int state_2_int(int* state) {
		return 8 * state[0] + 4 * state[1] + 2 * state[2] + 1 * state[3];
	}

	/*
		Detect collisions
		//state[0] particle going up
		//state[1] particle going right
		//state[2] particle going down
		//state[3] particle going left
	*/
	int collision(int state_repr){
		if (state_repr == 10){
			return 5;
		}
		else if (state_repr == 5){
			return 10;
		}
		else {
			return state_repr;
		}
	}


	/*
		detect collision for all cells
	*/
	void collision_all(){
		for (int i = 0; i < N*N; i++){
			states_vec[iterno % 2][i] = collision(states_vec[iterno % 2][i]);
		}
	}


	/*
		reflect particles at edges
	*/
	int edges_reflection(int cell_pos, int state){
		//top
		if (cell_pos/N == 1  && ((state & 8) == 8)){
			state = state | 2;
			state = state & 7;
		}
		//bottom
		if (cell_pos/N == N-2 && ((state & 2) == 2)){
			state = state | 8;
			state = state & 13;
		}
		//left
		if (cell_pos%N == 1 && ((state & 1) == 1)){
			state = state | 4;
			state = state & 14;
		}
		//right
		if ( cell_pos%N == N-2 && ((state & 4) == 4)){
			state = state | 1;
			state = state & 11;
		}
		return state;
	}

	/*
		reflect for all cells
			//should be just boundary?
			//how to balance it
				//multiple kernels?
	*/
	void edges_reflection_all(){
		for (int i = 0; i < N*N; i++){
			states_vec[iterno % 2][i] = edges_reflection(i, states_vec[iterno % 2][i]);
		}
	}


	/*
		Transport particles
		states[0]=top
		states[1]= left
		states[2]= right
		states[3]= bottom

		rememeber:
		//state[0] particle going up
		//state[1] particle going right
		//state[2] particle going down
		//state[3] particle going left
	*/

	int transport(int* states){
		int snew = 0;
		if ((states[0] & 2) == 2){ // from top going down
			snew = snew | 2;
		}
		if ((states[1] & 4) == 4){ // from left going right
			snew = snew | 4;
		}
		if ((states[2] & 1) == 1){ // from right going left
			snew = snew | 1;
		}
		if ((states[3] & 8) == 8){ // from bottom going up
			snew = snew | 8;
		}
		return snew;
	}


	/*
	transport for all cells
	*/
	void transport_all(){
		//no edges
		for (int i = 1; i < N-1; i++){
			for (int j = 1; j < N-1 ; j++){
				int pos = i*N + j;
				int states[4] = { states_vec[iterno % 2][pos - N], states_vec[iterno % 2][pos - 1],
					states_vec[iterno % 2][pos + 1], states_vec[iterno % 2][pos + N] };
				states_vec[(iterno + 1) % 2][pos] = transport(states);
			}
		}
		iterno++;
	}

	/*
	count the number of paritcles in a position
	*/
	int num_of_parts(int state_repr){
		int* state = int_2_state[state_repr];
		return state[0] + state[1] + state[2] + state[3];
	}



public:

	/*
		initialize a block + plus background noise
	*/
	void init_block(){
		//set background noise
		srand(42); //42 for testing
		double streng = 1.5;
		for (int i = 1; i < N - 1; i++){
			for (int j = 1; j < N - 1; j++){
				int state[4] = { (int)(streng * rand() / RAND_MAX),
					(int)(streng * rand() / RAND_MAX),
					(int)(streng * rand() / RAND_MAX),
					(int)(streng * rand() / RAND_MAX) };
				states_vec[0][i*N + j] = state_2_int(state);
			}
		}
		//set block going every direction
		for (int i = N / 2; i < N/2 + N/3; i++){
			for (int j = N / 2; j < N/2+ N/3; j++){
				states_vec[0][i*N + j] = 15;
			}
		}
	}

	/*
	update the state
	*/
	void update(){
		collision_all();
		edges_reflection_all();
		transport_all();
	}

	/*
		Print the state
	*/
	void print(){
		for (int i = 0; i < N; i++){
			for (int j = 0; j < N; j++){
				std::cout << num_of_parts(states_vec[iterno%2][i*N+j] )<< " ";
			}
			std::cout << "\n";
		}
		std::cout << std::endl;
	}

	void write(std::string filename){
		/* C++ style very slow
		std::ofstream file;
		file.open(filename);
		for (int i = 0; i < N; i++){
			for (int j = 0; j < N; j++){
				file << num_of_parts(states_vec[iterno % 2][i*N + j]) << " ";
			}
			file << "\n";
		}
		file << std::endl;
		file.close();
		*/

		// do it C
		FILE *f = fopen(filename.c_str(), "w");
		if (f == NULL)
		{
			printf("Error opening file!\n");
			exit(1);
		}
		for (int i = 0; i < N; i++){
			for (int j = 0; j < N; j++){
				fprintf(f, "%d ", num_of_parts(states_vec[iterno % 2][i*N + j]));
			}
			fprintf(f, "\n");
		}
		fclose(f);
	}

	/*
		pointer to data
	*/
	uint8_t* get_states_vec_pointer(){

		for (int i = 0; i < states_vec[iterno % 2].size(); i++){
			temp[i] = (uint8_t)(255-32*num_of_parts(states_vec[iterno % 2][i]));
		}
		return temp.data();
	}
};







//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
/*
	Lattice gas class for OpenCL simulation
*/
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//blocks size on GPU
//has to be set here, and in kernel too!!!
	//here for worksize calculation
#define CL_BLOCK_SIZE 10


//////////////////////////////////////////////////////////////////////////
/*
	Class for easier work with opencl
*/
//////////////////////////////////////////////////////////////////////////

class Lattice_gas_gpu{

	/*
		Constr, destr
			no copy, this object should not be copied
	*/
public:

	Lattice_gas_gpu(int N) : states_vec(N*N) {
		this->N = N;
		this->N_small = N - 2; //outer layer is just for convienience
		iterno = 0;
	}
	~Lattice_gas_gpu(){}

	/*
		data
	*/

private:
	int N; //size of the table
	int N_small; // interesting part size of the table

	//states in cells, represented by integer
	//I do synchronized update so 2 copies needed
	std::vector<int> states_vec;
	int iterno; //(using alternating vectors) 1.iter: vec1 is input-> vec2 is output
	//2.iter: vec2 is input->vec1 is output


private:
	//opencl related objects, and data
	cl_context context;
	cl_program program;
	cl_command_queue commandQueue;
	cl_kernel my_kernel;
	//buffers on the GPU (other device)
		//2 is needed for synchronous update
	cl_mem dev_buf[2],mem;

	/*
		Methods
	*/

	/*
		Init function, to deal with most opencl stuff
			- choose opencl platfrom, and device
			- create context, commandque
			- load kernel files, and build them

		Its quite long and verbose
			- might be changed in the future
	*/
	int opencl_initialize(size_t platform_choice, size_t device_choice)
	{


		cl_int	status = 0;

		//Getting OpenCL platforms and choose an available one.
		cl_uint numPlatforms;				//the NO. of platforms
		cl_platform_id* platforms = NULL; 	//id of available platforms
		cl_platform_id 	platform = NULL;	//id of the chosen platform

		//getting NO. of platforms
		status = clGetPlatformIDs(0, NULL, &numPlatforms);
		if (status != CL_SUCCESS)
		{
			std::cerr << "Error: Getting platforms!" << std::endl;
			std::cerr << "Error number= " << status << std::endl;
			return status;
		}

		//Choosing platform
		if (numPlatforms > 0)
		{
			//getting platform ids
			platforms = new cl_platform_id[numPlatforms];
			status = clGetPlatformIDs(numPlatforms, platforms, NULL);

			//printing platform names
			std::cout << "\nPlatform info:" << std::endl;
			for (unsigned int i = 0; i<numPlatforms; i++)
			{
				//get platform name size
				size_t platform_name_size;
				status = clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, 0, NULL, &platform_name_size);

				//get platform name
				char* platform_name = new char[platform_name_size];
				status = clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, platform_name_size, platform_name, NULL);


				//get platform version size
				size_t platform_version_size;
				status = clGetPlatformInfo(platforms[i], CL_PLATFORM_VERSION, 0, NULL, &platform_version_size);

				//get platform version
				char* platform_version = new char[platform_version_size];
				status = clGetPlatformInfo(platforms[i], CL_PLATFORM_VERSION, platform_version_size, platform_version, NULL);


				//print info
				std::cout << i << ". platform:\t" << platform_name << "\n";
				std::cout << i << ". version:\t " << platform_version << "\n";
				std::cout << std::endl;

				delete[] platform_name;
				delete[] platform_version;
			}


			//choosing platform
			std::cout << "\nChoose platform: (0)" << std::endl;
			/*int platform_choice = 0;
			std::string temp_line;
			getline(std::cin, temp_line);
			std::stringstream temp_sstr;
			temp_sstr << temp_line;
			temp_sstr >> platform_choice;
			std::cout << "platform choice" << platform_choice << std::endl;
			*/
			platform = platforms[platform_choice];

			delete[] platforms;
		}

		//Query the platform and choose the  device
		cl_uint		numDevices = 0; 	//NO. of devices
		cl_device_id	*devices;		// device ids
		cl_device_id	device;			//id of chosen device

		//getting number of devices
		status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 0, NULL, &numDevices);
		devices = new cl_device_id[numDevices];

		//getting device ids
		status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, numDevices, devices, NULL);

		//printing device info
		std::cout << "\nDevice info:" << std::endl;
		for (unsigned int i = 0; i<numDevices; i++)
		{
			//get device vendor size
			size_t device_vendor_size;
			status = clGetDeviceInfo(devices[i], CL_DEVICE_VENDOR, 0, NULL, &device_vendor_size);

			//get device vendor
			char* device_vendor = new char[device_vendor_size];
			status = clGetDeviceInfo(devices[i], CL_DEVICE_VENDOR, device_vendor_size, device_vendor, NULL);


			//get device name size
			size_t device_name_size;
			status = clGetDeviceInfo(devices[i], CL_DEVICE_NAME, 0, NULL, &device_name_size);

			//get device name
			char* device_name = new char[device_name_size];
			status = clGetDeviceInfo(devices[i], CL_DEVICE_NAME, device_name_size, device_name, NULL);


			//get devicetype 
			cl_device_type device_type;
			status = clGetDeviceInfo(devices[i], CL_DEVICE_TYPE, sizeof(cl_device_type), &device_type, NULL);


			//get device version size
			size_t device_version_size;
			status = clGetDeviceInfo(devices[i], CL_DEVICE_VERSION, 0, NULL, &device_version_size);

			//get device version
			char* device_version = new char[device_version_size];
			status = clGetDeviceInfo(devices[i], CL_DEVICE_VERSION, device_version_size, device_version, NULL);


			//print info
			std::cout << i << ". device vendor:\t" << device_vendor << std::endl;
			std::cout << "            name:\t" << device_name << std::endl;

			//device type 
			if (device_type == CL_DEVICE_TYPE_CPU)
				std::cout << "            type:\tCPU" << std::endl;
			if (device_type == CL_DEVICE_TYPE_GPU)
				std::cout << "            type:\tGPU" << std::endl;
			if (device_type == CL_DEVICE_TYPE_ACCELERATOR)
				std::cout << "            type:\tACCELERATOR" << std::endl;
			if (device_type == CL_DEVICE_TYPE_DEFAULT)
				std::cout << "            type:\tDEFAULT" << std::endl;


			std::cout << "            version:\t" << device_version << std::endl;

			delete[] device_vendor;
			delete[] device_name;
			delete[] device_version;
		}

		//choosing device
		std::cout << "\nChoose device: (0)" << std::endl;
		/*int device_choice = 0;
		std::string temp_line1;
		getline(std::cin, temp_line1);
		std::stringstream temp_sstr1;
		temp_sstr1 << temp_line1;
		temp_sstr1 >> device_choice;
		*/
		device = devices[device_choice];



		// Create CL context properties, add WGL context & handle to DC
		cl_context_properties properties[] = {
			CL_GL_CONTEXT_KHR, (cl_context_properties)wglGetCurrentContext(), // WGL Context
			CL_WGL_HDC_KHR, (cl_context_properties)wglGetCurrentDC(), // WGL HDC
			CL_CONTEXT_PLATFORM, (cl_context_properties)platform, // OpenCL platform
			0
		};
	
		
		//Create context
		context = clCreateContext(properties, 1, &devices[device_choice], NULL, NULL, &status);
		if (status != 0)
		{
			std::cerr << "ERROR creating context: " << status << std::endl;
			return status;
		}
		

		//Creating command queue associated with the context
		commandQueue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status);
		if (status != 0)
		{
			std::cerr << "ERROR creating commandqueue: " << status << std::endl;
			return status;
		}

		//open kernel file and convert it to char array
		//const char *filename = kernel_filename.c_str();
		std::string sourceStr = convertToString("lattice_gas.cl");
		const char *source = sourceStr.c_str();
		size_t sourceSize[] = { strlen(source) };

		//Create program object
		program = clCreateProgramWithSource(context, 1, &source, sourceSize, &status);
		if (status != 0)
		{
			std::cout << "ERROR creating program: " << status << std::endl;
			return status;
		}

		//Building program 
		//only prints log if there was an error
		//if there are only warnings it is not printed
		status = clBuildProgram(program, numDevices, devices, NULL, NULL, NULL);
		if (status != 0)
		{
			//print ERROR but do not quit, there may be just warnings
			std::cerr << "ERROR building program: " << status << std::endl;

			//Getting build log size
			size_t logsize = 0;
			clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &logsize);
			std::cout << logsize << std::endl;

			//Getting build log
			char* log = new char[logsize];
			clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, logsize, log, NULL);

			//print log info
			std::cout << "log:\n " << log << std::endl;
			delete[] log;

			return status;
		}

		//print log
		std::cerr << "log from building program: " << status << std::endl;

		//Getting build log size
		size_t logsize = 0;
		clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &logsize);
		std::cout << logsize << std::endl;

		//Getting build log
		char* log = new char[logsize];
		clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, logsize, log, NULL);

		//print log info
		std::cout << "log:\n " << log << std::endl;
		delete[] log;

		return status;
	}

	//create kernel
	int opencl_create_kernel(std::string kernel_name)
	{
		//error variable
		cl_int status = 0;

		// Create kernel object
		my_kernel = clCreateKernel(program, kernel_name.c_str(), &status);
		if (status != 0)
		{
			std::cerr << "ERROR creating kernel: " << status << std::endl;
			return status;
		}

		return status;
	}

private:
	//load the kernel file into a null terminated string 
	std::string convertToString(std::string infilename)
	{
		std::string str;

		// load the kernel file into a null terminated string 
		//open file in binary i/o 
		std::ifstream infile(infilename.c_str(), std::ios::binary | std::ios::ate);
		//check file 
		if (!(infile))
		{
			std::cout << "\nERROR CAN'T OPEN KERNEL FILE: " << infilename << "\n" << std::endl;
			return NULL;
		}

		//get the size of the file 
		std::ifstream::pos_type size;
		size = infile.tellg();
		//go to the begginging of file								 
		infile.seekg(0, std::ios::beg);

		//read file
		str.resize(size);
		infile.read((char*)(str.c_str()), size);
		//append "\0"
		str += '\0';

		return str;
	}

	//allocates buffers on devices
	int opencl_copy_mem(GLuint texture)
	{
		//error variable
		cl_int status = 0;
		//create from opengl texture
		mem = clCreateFromGLTexture2D(context, CL_MEM_WRITE_ONLY, GL_TEXTURE_2D, 0, texture, &status);
		if (status != 0)
		{
			std::cerr << "ERROR creating buffers: " << status << std::endl;
			return status;
		}
		
		

		// Allocate memory on device

		//data
		dev_buf[0] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, N * N * sizeof(cl_int), states_vec.data(), &status);
		//error check
		if (status != 0)
		{
			std::cerr << "ERROR creating buffers: " << status << std::endl;
			return status;
		}
		//data
		dev_buf[1] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, N * N * sizeof(cl_int), states_vec.data(), &status);
		//error check
		if (status != 0)
		{
			std::cerr << "ERROR creating buffers: " << status << std::endl;
			return status;
		}
		return status;

	}

	//Setting kernel arguments
	//with the buffers on the GPU (device)
	int set_kern_arg()
	{
		//error variable
		cl_int status = 0;

		//there is no size_t in opencl kernel
		unsigned int Nui = (unsigned int)N;

		//kernel
		status = clSetKernelArg(my_kernel, 0, sizeof(cl_mem), &dev_buf[iterno%2]);
		status |= clSetKernelArg(my_kernel, 1, sizeof(cl_mem), &dev_buf[(iterno+1) % 2]);
		status |= clSetKernelArg(my_kernel, 2, sizeof(cl_mem), &mem);
		status |= clSetKernelArg(my_kernel, 3, sizeof(unsigned int), &Nui);

		//error check
		if (status != 0)
		{
			std::cerr << "ERROR setting kernel arguments: " << status << std::endl;
			return status;
		}

		return status;
	}

	//this function launches the kernel
	int call_kernel(size_t worksize)
	{
		//error variable
		cl_int status = 0;

		//get  control
		glFinish();
		clEnqueueAcquireGLObjects(commandQueue, 1, &mem, 0, 0, NULL);

		//size_t worksize = d*d;
		// Running the kernel.
		size_t global_work_size[1] = { worksize };
		status = clEnqueueNDRangeKernel(commandQueue, my_kernel, 1, NULL, global_work_size, NULL, 0, NULL, NULL);
		if (status != 0)
		{
			std::cerr << "ERROR running kernel1: " << status << std::endl;
			return status;
		}

		//release control
		clFinish(commandQueue);
		clEnqueueReleaseGLObjects(commandQueue, 1, &mem, 0, 0, NULL);

		iterno++;

		return status;
	}

	//reads back result from device to host
	int read_from_dev()
	{
		cl_int status = 0;
		//Read the result back to host memory.
		status = clEnqueueReadBuffer(commandQueue, dev_buf[iterno%2], CL_TRUE, 0, N*N * sizeof(cl_int), this->states_vec.data(), 0, NULL, NULL);
		if (status != 0)
		{
			std::cout << "ERROR reading buffer: " << status << std::endl;
		}
		return status;
	}

	//helper function transfrom state to int
	int state_2_int(int* state) {
		return 8 * state[0] + 4 * state[1] + 2 * state[2] + 1 * state[3];
	}

	/*
	count the number of paritcles in a position
	*/
	int num_of_parts(int state_repr){
		int* state = int_2_state[state_repr];
		return state[0] + state[1] + state[2] + state[3];
	}

//////////////////////////////////////////////////////////////////////////
public:

	/*
	initialize a block + plus background noise
	*/
	void init_block(int platform_choice, int device_choice,std::string kernel_name,GLuint texture){
	
		/* 
			data on cpu
		*/
		//set background noise
		srand(42); //42 for testing
		double streng = 1.5;
		for (int i = 1; i < N - 1; i++){
			for (int j = 1; j < N - 1; j++){
				int state[4] = { (int)(streng * rand() / RAND_MAX),
					(int)(streng * rand() / RAND_MAX),
					(int)(streng * rand() / RAND_MAX),
					(int)(streng * rand() / RAND_MAX) };
				states_vec[i*N + j] = state_2_int(state);
			}
		}
		//set block going every direction
		for (int i = N / 2; i < N / 2 + N / 3; i++){
			for (int j = N / 2; j < N / 2 + N / 3; j++){
				states_vec[i*N + j] = 15;
			}
		}

		/*
			gpu
		*/
		//init
		opencl_initialize(platform_choice, device_choice);
		//compile kernel, copy data
		opencl_copy_mem(texture);
		//create kernel
		opencl_create_kernel(kernel_name);
		//set args
		set_kern_arg();
	
	}

	/*
		update state
	*/
	void update(){

		//worksize calculation
		size_t worksize = (N-2) / (CL_BLOCK_SIZE-2) *(N-2) / (CL_BLOCK_SIZE-2);

		//reset kern arg 
			//to change input, and output
		set_kern_arg();

		//call kernel
		call_kernel(worksize);
	}

	/*
		write state
	*/
	void write(std::string filename){

		//read back memory from device
		read_from_dev();

		// write it
		FILE *f = fopen(filename.c_str(), "w");
		if (f == NULL)
		{
			printf("Error opening file!\n");
			exit(1);
		}
		for (int i = 0; i < N; i++){
			for (int j = 0; j < N; j++){
				fprintf(f, "%d ", num_of_parts(states_vec[i*N + j]));
			}
			fprintf(f, "\n");
		}
		fclose(f);
	}

};




/////////////////////////////////////////////////////////////////////////
//visualisation with opengl
/////////////////////////////////////////////////////////////////////////

//vertex shader
const GLchar* vertexSource =
"#version 150 core\n"
"in vec2 position;"
"in vec3 color;"
"in vec2 texcoord;"
"out vec3 Color;"
"out vec2 Texcoord;"
"void main() {"
"   Color = color;"
"   Texcoord = texcoord;"
"   gl_Position = vec4(position, 0.0, 1.0);"
"}";

//texture shader
const GLchar* fragmentSource =
"#version 150 core\n"
"in vec3 Color;"
"in vec2 Texcoord;"
"out vec4 outColor;"
"uniform sampler2D tex;"
"void main() {"
"   outColor = texture(tex, Texcoord) * vec4(1.0,1.0,1.0, 1.0);"
"}";


int simul_gpu(){ 

	///////////////////////////////////////////////////////////////
	//opengl stuff
	///////////////////////////////////////////////////////////////

	sf::ContextSettings settings;
	settings.depthBits = 24;
	settings.stencilBits = 8;
	sf::Window window(sf::VideoMode(800, 800, 32), "OpenGL", sf::Style::Titlebar | sf::Style::Close, settings);


	// Initialize GLEW
	glewExperimental = GL_TRUE;
	glewInit();

	// Create Vertex Array Object
	GLuint vao;
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	// Create a Vertex Buffer Object and copy the vertex data to it
	GLuint vbo;
	glGenBuffers(1, &vbo);

	GLfloat vertices[] = {
		//  Position   Color             Texcoords
		-0.5f, 0.5f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, // Top-left
		0.5f, 0.5f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, // Top-right
		0.5f, -0.5f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f, // Bottom-right
		-0.5f, -0.5f, 1.0f, 1.0f, 1.0f, 0.0f, 1.0f  // Bottom-left
	};

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

	// Create an element array
	GLuint ebo;
	glGenBuffers(1, &ebo);

	GLuint elements[] = {
		0, 1, 2,
		2, 3, 0
	};

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(elements), elements, GL_STATIC_DRAW);

	// Create and compile the vertex shader
	GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vertexShader, 1, &vertexSource, NULL);
	glCompileShader(vertexShader);

	// Create and compile the fragment shader
	GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fragmentShader, 1, &fragmentSource, NULL);
	glCompileShader(fragmentShader);

	// Link the vertex and fragment shader into a shader program
	GLuint shaderProgram = glCreateProgram();
	glAttachShader(shaderProgram, vertexShader);
	glAttachShader(shaderProgram, fragmentShader);
	//this was in tutorial, but gives some error
	//glBindFragDataLocation(shaderProgram, 0, "outColor");
	glLinkProgram(shaderProgram);
	glUseProgram(shaderProgram);

	// Specify the layout of the vertex data
	GLint posAttrib = glGetAttribLocation(shaderProgram, "position");
	glEnableVertexAttribArray(posAttrib);
	glVertexAttribPointer(posAttrib, 2, GL_FLOAT, GL_FALSE, 7 * sizeof(GLfloat), 0);

	GLint colAttrib = glGetAttribLocation(shaderProgram, "color");
	glEnableVertexAttribArray(colAttrib);
	glVertexAttribPointer(colAttrib, 3, GL_FLOAT, GL_FALSE, 7 * sizeof(GLfloat), (void*)(2 * sizeof(GLfloat)));

	GLint texAttrib = glGetAttribLocation(shaderProgram, "texcoord");
	glEnableVertexAttribArray(texAttrib);
	glVertexAttribPointer(texAttrib, 2, GL_FLOAT, GL_FALSE, 7 * sizeof(GLfloat), (void*)(5 * sizeof(GLfloat)));



	//////////////////////////////////////////////////////////////////
	// simul
	// CPU
	//////////////////////////////////////////////////////////////////
	int N = 100 * 8 + 2; //table size
	int N_iter = 10; //number of iteration
	int N_write = 10; //write result at every N_write iteration
	int i_file = 0; //filenane counter

	//lattice gas class
	Lattice_gas_gpu my_lat_gas_gpu(N);

	// Create texture
	GLuint tex;
	glGenTextures(1, &tex);
	glBindTexture(GL_TEXTURE_2D, tex);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, N, N, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, 0);
	//glBindTexture(GL_TEXTURE_2D, 0);


	//set initial stage
	my_lat_gas_gpu.init_block(3, 0, "lattice_gas_0", tex);
	//load initial into texture?

	int i = 0;
	while (window.isOpen())
	{

		//report to stdout
		//std::cout << i << "\n";
		//i++;

		// Clear the screen to black
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);

		// Draw a rectangle from the 2 triangles using 6 indices
		glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

		// Swap buffers
		window.display();

		//update the state
		my_lat_gas_gpu.update();

		sf::Event windowEvent;
		while (window.pollEvent(windowEvent))
		{
			switch (windowEvent.type)
			{
			case sf::Event::Closed:
				window.close();
				break;
			}
		}
	}

	//release objects
	glDeleteTextures(1, &tex);
	glDeleteProgram(shaderProgram);
	glDeleteShader(fragmentShader);
	glDeleteShader(vertexShader);

	glDeleteBuffers(1, &ebo);
	glDeleteBuffers(1, &vbo);

	glDeleteVertexArrays(1, &vao);

	return 0;
}


int simul_cpu(){
	///////////////////////////////////////////////////////////////
	//opengl stuff
	///////////////////////////////////////////////////////////////

	sf::ContextSettings settings;
	settings.depthBits = 24;
	settings.stencilBits = 8;
	sf::Window window(sf::VideoMode(800, 800, 32), "OpenGL", sf::Style::Titlebar | sf::Style::Close, settings);


	// Initialize GLEW
	glewExperimental = GL_TRUE;
	glewInit();

	// Create Vertex Array Object
	GLuint vao;
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	// Create a Vertex Buffer Object and copy the vertex data to it
	GLuint vbo;
	glGenBuffers(1, &vbo);

	GLfloat vertices[] = {
		//  Position   Color             Texcoords
		-0.5f, 0.5f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, // Top-left
		0.5f, 0.5f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, // Top-right
		0.5f, -0.5f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f, // Bottom-right
		-0.5f, -0.5f, 1.0f, 1.0f, 1.0f, 0.0f, 1.0f  // Bottom-left
	};

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

	// Create an element array
	GLuint ebo;
	glGenBuffers(1, &ebo);

	GLuint elements[] = {
		0, 1, 2,
		2, 3, 0
	};

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(elements), elements, GL_STATIC_DRAW);

	// Create and compile the vertex shader
	GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vertexShader, 1, &vertexSource, NULL);
	glCompileShader(vertexShader);

	// Create and compile the fragment shader
	GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fragmentShader, 1, &fragmentSource, NULL);
	glCompileShader(fragmentShader);

	// Link the vertex and fragment shader into a shader program
	GLuint shaderProgram = glCreateProgram();
	glAttachShader(shaderProgram, vertexShader);
	glAttachShader(shaderProgram, fragmentShader);
	//this was in tutorial, but gives some error
	//glBindFragDataLocation(shaderProgram, 0, "outColor");
	glLinkProgram(shaderProgram);
	glUseProgram(shaderProgram);

	// Specify the layout of the vertex data
	GLint posAttrib = glGetAttribLocation(shaderProgram, "position");
	glEnableVertexAttribArray(posAttrib);
	glVertexAttribPointer(posAttrib, 2, GL_FLOAT, GL_FALSE, 7 * sizeof(GLfloat), 0);

	GLint colAttrib = glGetAttribLocation(shaderProgram, "color");
	glEnableVertexAttribArray(colAttrib);
	glVertexAttribPointer(colAttrib, 3, GL_FLOAT, GL_FALSE, 7 * sizeof(GLfloat), (void*)(2 * sizeof(GLfloat)));

	GLint texAttrib = glGetAttribLocation(shaderProgram, "texcoord");
	glEnableVertexAttribArray(texAttrib);
	glVertexAttribPointer(texAttrib, 2, GL_FLOAT, GL_FALSE, 7 * sizeof(GLfloat), (void*)(5 * sizeof(GLfloat)));

	// Create texture
	GLuint tex;
	glGenTextures(1, &tex);


	//////////////////////////////////////////////////////////////////
	// simul
	// CPU
	//////////////////////////////////////////////////////////////////
	int N = 100 * 8;// +2; //table size
	int N_iter = 10; //number of iteration
	int N_write = 10; //write result at every N_write iteration
	int i_file = 0; //filenane counter

	//lattice gas class
	Lattice_gas my_lat_gas(N);

	//set initial stage
	my_lat_gas.init_block();


	int i = 0;
	while (window.isOpen())
	{

		//load state into texture
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, N, N, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, my_lat_gas.get_states_vec_pointer());
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

		//report to stdout
		//std::cout << i << "\n";
		//i++;

		// Clear the screen to black
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);

		// Draw a rectangle from the 2 triangles using 6 indices
		glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);


		// Swap buffers
		window.display();

		//update the state
		my_lat_gas.update();

		sf::Event windowEvent;
		while (window.pollEvent(windowEvent))
		{
			switch (windowEvent.type)
			{
			case sf::Event::Closed:
				window.close();
				break;
			}
		}
	}

	//release objects
	glDeleteTextures(1, &tex);
	glDeleteProgram(shaderProgram);
	glDeleteShader(fragmentShader);
	glDeleteShader(vertexShader);

	glDeleteBuffers(1, &ebo);
	glDeleteBuffers(1, &vbo);

	glDeleteVertexArrays(1, &vao);

	return 0;
}
