#include "Node.cuh"
#include <iostream>
#include <algorithm>
#include <stack>
#include <vector>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// ------- $TD : 'WRITE MORE ORGANIZED WITH INITIAL FUNCS';
Node::Node(int size, int s, int s0)
{
	this->size = size;


	this->S = s;
	this->S0 = s0;
	this->points = 0;
	

	optimalX = new int[this->size];
	optimalY = new int[this->size];
}


void Node::setInitials(int iSize){
	this->P = new int[iSize];
	for (int i = 0; i < iSize; i++)
		this->P[i] = -1;
}

// --- $ID : 'MAYBE ALLOC MEMORY ON ONE WAY';
void Node::setInitialMatrix(int *sourceMatrix) {
	this->baseSize = this->size;
	this->setInitials(this->baseSize);
	//this->M = new int*[size];
	this->M = (int *)malloc(this->size * this->size * sizeof(int));
	this->translateX = new int[size];
	this->translateY = new int[size];

	for (int i = 0; i < size; i++) {
		//this->M[i] = new int[size];
		this->translateX[i] = i;
		this->translateY[i] = i;
		for (int j = 0; j < size; j++) {
			//this->M[i][j] = sourceMatrix[i][j];
			int offset = i * this->size + j;
			this->M[offset] = sourceMatrix[offset];
		}
	}

}

void Node::setMatrix(int *m) {
	//this->M = new int*[size];
	this->M = (int*)malloc(this->size * this->size * sizeof(int));
	/*for (int i = 0; i < size; i++)
		this->M[i] = new int[size];*/

	for (int i = 0; i < size; i++) {
		//this->M[i] = new int[size];
		for (int j = 0; j < size; j++) {

			int offset = i * this->size + j;
			this->M[offset] = m[offset];
		}
	}
}

// ---- $TD: 'REWRITE WITH MORE PRODUCITY BY BINARY COPYING'
void Node::setMatrixWithRemoveExclude(int *source, int row, int col) {

	this->M = (int *)malloc(this->size * this->size * sizeof(int));
	int originalSize = this->size + 1;

	int ni = 0;

	for (int i = 0; i < originalSize; i++){
		if (i != row){
			for (int j = 0; j < originalSize; j++){
				if (j != col){
					this->M[ni] = source[i * originalSize + j];
					ni++;
				}
			}
		}
	}

}

Node::~Node()
{

}

void Node::printMatrix() {
	std::cout << std::endl;
	std::cout << std::endl;

	for (int i = 0; i < this->size; i++) {
		for (int j = 0; j < this->size; j++)
			std::cout << this->M[i * this->size + j] << " ";
		std::cout << std::endl;
	}

	std::cout << std::endl;
	std::cout << std::endl;
}

// $TD: 'REPLACE BY STANDART STD OR BOOST';
int Node::getArrayMinValue(int restrictVal, /*int *row*/ int row) {
	int min = InfityMaxValue;
	for (int i = 0; i < size; i++) {
		int offset = row * this->size + i;
		if (this->M[offset] > restrictVal && this->M[offset] < min)
			min = this->M[offset];
	}
	return min;
}

// $TD : 'REWRITE BY INIT WAY : CHANGES 4 CYCLES BY 2 AND
//								SPEED UP  getPathForRemove ---> BY INDEXING OPERATION
//								SPEED UP  subMRows ---> return this.minRowsEls[row];
//								SPPED UP  subMCols ---> return this.minRowsEls[col] - this.minRowsEls[row] - checkin;
void Node::subMinRowsAndCorrect() {

	//this->printMatrix();

	//int matrixSize = this->size * this->size;
	//int allocatedSize = matrixSize * sizeof(int);

	//int *d_m;
	//int *d_t;

	//int s = 100;
	//int *d_s;
	//int *d_a;

	//cudaMalloc((void**)&d_s, this->size * sizeof(int));
	//cudaMalloc((void**)&d_a, this->size * this->size * sizeof(int));
	/*cudaMalloc((void**)&d_m, allocatedSize);
	cudaMalloc((void**)&d_t, this->size * sizeof(int));

	cudaMemcpy(d_m, this->M, allocatedSize, cudaMemcpyHostToDevice);
	cudaMemcpy(d_s, &s, sizeof(int), cudaMemcpyHostToDevice);

	matrixRowCorrect << <this->size, this->size >> >(d_m, d_t, d_s);

	int *response = (int*)malloc(allocatedSize);
	int *temp = (int*)malloc(this->size * sizeof(int));
	cudaMemcpy(response, d_m, allocatedSize, cudaMemcpyDeviceToHost);
	cudaMemcpy(temp, d_t, this->size * sizeof(int), cudaMemcpyDeviceToHost);
	*/
	//this->testMatrixAdduction(this->M);

	/* PREV STABLE VERSION*/

	//int *mins = new int[this->size];
	//int localSDelta = 0;
	//int min = InfityMaxValue;
	//for (int i = 0; i < size; i++) {
	//	min = getArrayMinValue(-1, i);
	//	localSDelta += min < InfityMaxValue ? min : 0;
	//	mins[i] = min;
	//	for (int j = 0; j < size; j++) {
	//		int offset = i * this->size + j;
	//		if (this->M[offset] > -1)
	//			this->M[offset] -= min;
	//		/*if (this->M[i][j] > -1)
	//			this->M[i][j] -= min;*/
	//	}
	//	//r = subMinRowsAndCorrect(min, this->M[i], this->size);
	//}
	//S += localSDelta;


	int *mins = new int[this->size];
	int localSDelta = 0;
	int min = InfityMaxValue;
	for (int i = 0; i < size; i++) {
		min = getArrayMinValue(-1, i);
		localSDelta += min < InfityMaxValue ? min : 0;
		mins[i] = min;
		if (min > 0){
			for (int j = 0; j < size; j++) {
				int offset = i * this->size + j;
				if (this->M[offset] > -1)
					this->M[offset] -= min;
				/*if (this->M[i][j] > -1)
				this->M[i][j] -= min;*/
			}
		}
		//r = subMinRowsAndCorrect(min, this->M[i], this->size);
	}
	S += localSDelta;

	/*cudaMemcpy(d_s, mins, this->size * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_a, this->M, this->size * this->size * sizeof(int), cudaMemcpyHostToDevice);
	arrayReduce << <this->size, this->size >> >(d_a, d_s);
	cudaMemcpy(this->M, d_a, this->size * this->size *  sizeof(int), cudaMemcpyDeviceToHost);*/

	//	this->printMatrix();

	/*cudaFree(d_a);
	cudaFree(d_s);*/

	//this->printMatrix();
}

// REWRITE BY subMinRowsAndCorrect
void Node::subMinColsAndCorrect() {
	int correlation = 0;
	int localMin = InfityMaxValue;
	/*int *colMins = new int[size];*/
	int offsetj = 0;

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			 offsetj = i + this->size * j;
			if (this->M[offsetj] < localMin && this->M[offsetj] > -1)
				localMin = this->M[offsetj];
		}
		//colMins[i] = localMin;

		if (localMin != 0){
			for (int j = 0; j < size; j++)
				if (this->M[offsetj] > -1)
					this->M[offsetj] -= localMin;
		}


		correlation += localMin == InfityMaxValue ? 0 : localMin;
		localMin = InfityMaxValue;
	}

	/*for (int i = 0; i < size; i++)
	if (colMins[i] != 0){
		for (int j = 0; j < size; j++)
			if (this->M[i + this->size * j] > -1)
				this->M[i + this->size * j] -= colMins[i];
			}*/


	S += correlation;
}

// REWRITE WITH PREVIOS COMMENT LOGIC : SPEED UP
void Node::getPathForRemove(int &rowE, int &colE) {
	int max = -1;
	int rowMin = InfityMaxValue;
	int colMin = InfityMaxValue;
	int ioffset, loffset;

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
				ioffset = i*this->size;
			if (this->M[ioffset + j] == 0) {

				rowMin = InfityMaxValue;
				colMin = InfityMaxValue;

				for (int r = 0; r < size; r++) {
					if (this->M[ioffset + r] != -1 && j != r && this->M[ioffset + r] < rowMin) {
						rowMin = this->M[ioffset + r];
					}

						loffset = r * this->size;
					if (this->M[loffset + j] != -1 && r != i && this->M[loffset + j] < colMin) {
						colMin = this->M[loffset + j];
					}

				}
			
				if ((colMin + rowMin) > max)
				{
					max = colMin + rowMin;
					rowE = i;
					colE = j;
				}
			}
		}
	}
}

void Node::invokeAdduction() {
	this->subMinRowsAndCorrect();
	this->subMinColsAndCorrect();
}

void Node::copySessionDescription(int *rowD, int *colD) {
	this->translateX = new int[this->size];
	this->translateY = new int[this->size];



	int ri = 0;
	int rj = 0;
	for (int i = 0; i < this->size + 1; i++)
	{
		if (rowD[i] != -1){
			this->translateX[ri] = rowD[i];
			ri++;
		}

		if (colD[i] != -1){
			this->translateY[rj] = colD[i];
			rj++;
		}
	}
}


Node* Node::leftBranching(int row, int col) {
	Node *node = new Node(this->size, this->S, this->S);
	node->baseSize = this->baseSize;

	node->cudaCleanCopy(this->M);
	//node->setMatrix(this->M);

	node->copySessionDescription(this->translateX, this->translateY);
	node->M[row * this->size + col] = -1;
	node->invokeAdduction();

	//node->printMatrix();

	node->setInitials(this->baseSize);
	//for (int i = 0; i < this->baseSize; i++) {
	//	if (this->P[i] != -1)
	//		node->P[i] = this->P[i];
	//	else
	//		node->P[i] = -1;
	//}

	memcpy(node->P, this->P, this->baseSize * sizeof(int));

	return node;
}

int Node::getRealElement(int *dataDescription, int ind) {

	for (int i = 0; i < this->size; i++)
	{
		if (dataDescription[i] == ind)
			return i;
	}
	return -1;
}

void Node::printArray(int size, int *arr){
	std::cout << std::endl;
	std::cout << std::endl;
	for (int i = 0; i < size; i++)
		std::cout << arr[i] << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
}

// REMOVE SOME FUNCTIONS
// $DIRTY CODE
Node* Node::rightBranching(int row, int col) {
	Node *node = new Node(this->size - 1, this->S, this->S);
	node->baseSize = this->baseSize;
	node->setInitials(this->baseSize);

	int realRow = this->translateX[row];
	int realCol = this->translateY[col];

	this->P[realRow] = realCol;

	memcpy(node->P, this->P, this->baseSize * sizeof(int));

	//printArray(this->baseSize, node->P);

	
	this->points++;
	node->points = this->points;


	///	printArray(this->baseSize, node->P);


	if (realRow > -1 && realCol > -1) {
		if (points > 1)
			node->handlePodcycles(realRow, realCol);

		/*printArray(this->size, translateX);*/

		//auto tmp = *std::max_element(this->translateX, this->translateX + sizeof(this->translateX) / sizeof(int));

		if (row < this->size && col < this->size){
			this->translateX[row] = -1;
			this->translateY[col] = -1;
		}

		/*this->printArray(this->size, this->translateX);
		this->printArray(this->size, this->translateY);*/

		int returnRow = getRealElement(this->translateX, realCol);
		int returnCol = getRealElement(this->translateY, realRow);

		if (returnRow > -1 && returnCol > -1)
			this->M[returnRow * this->size + returnCol] = -1;

		node->setMatrixWithRemoveExclude(this->M, row, col);

		//node->printMatrix();

		node->invokeAdduction();

		//node->printMatrix();

		node->copySessionDescription(this->translateX, this->translateY);

	}

	return node;
}

void Node::handleStraightforwardMatrix() {
	for (int i = 0; i < this->size; i++) {
		if (P[i] == -1) {
			int offsetX = i*this->size;
			for (int j = 0; j < this->size; j++)
			if (this->M[offsetX + j] == 0)
				this->P[i] = j;
		}
	}
}

/* --> NO USE IN PRODUCTION !!! <-------------------- HANDLE LOCAL CYCLES (REWRITE) $DIRTY CODE ----------------*/
int Node::getHead(int tail) {
	for (int i = 0; i < this->baseSize; i++) {
		if (this->P[i] == tail)
			return i;
	}
	return -1;
}

int Node::getTail(int head) {
	return this->P[head];
}

void Node::handlePodcycles(int &a, int &b) {
	int source = a;
	int destiny = b;


	std::vector<int> localCycle;
	int head = 0;
	int tail = 0;
	bool finish = false;

	localCycle.push_back(a);
	localCycle.push_back(b);
	int countIterations = 0;
	while (!finish)
	{
		countIterations++;
		//if (localCycle.size() > this->baseSize + 10)
		//{
		//	std::cout << "CYCLE ERROR !!!!!!" << std::endl;
		//	throw "ss";
		//}

		if (head != -1) {
			head = getHead(a);
			if (b == head) // whaaaat ?
				break;
			if (head != -1) {
				localCycle.insert(localCycle.begin(), head);
				a = head;
			}
		}

		if (tail != -1) {
			tail = getTail(b);
			if (tail == head) // whaaaat ?
				break;
			if (tail != -1) {
				localCycle.push_back(tail);
				b = tail;
			}
		}

		if (head == -1 && tail == -1)
			finish = true;
	}
	a = localCycle.front();
	b = localCycle.back();
	localCycle.clear();
}





__global__ void modifyArrayKernel(int *val, int *arr){
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < 6 && arr[i] > -1)
		arr[i] = arr[i] - *val;
}

#define N 2

__global__ void MatAdd(int A[][N], int B[][N], int C[][N]){
	int i = threadIdx.x;
	int j = threadIdx.y;

	C[i][j] = A[i][j] + B[i][j];
}

__global__ void testMatrix(int **M, int **R){
	int i = threadIdx.x;
	int j = threadIdx.y;

	R[i][j] = M[i][j] - 10;
}

__global__ void testKernel(int *s, const int *re){

	__shared__ int temp[1];

	int i = threadIdx.x;
	if (re[i] > -1 && re[i] < temp[0])
		temp[0] = re[i];

	__syncthreads();

	*s = temp[0];
}

void Node::testMatrixAdduction(int *M){


	int A[N][N] = { { 1, 2 }, { 3, 4 } };
	int B[N][N] = { { 5, 6 }, { 7, 8 } };
	int C[N][N] = { { 0, 0 }, { 0, 0 } };

	int(*pA)[N], (*pB)[N], (*pC)[N];

	cudaMalloc((void**)&pA, (N*N)*sizeof(int));
	cudaMalloc((void**)&pB, (N*N)*sizeof(int));
	cudaMalloc((void**)&pC, (N*N)*sizeof(int));

	cudaMemcpy(pA, A, (N*N)*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(pB, B, (N*N)*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(pC, C, (N*N)*sizeof(int), cudaMemcpyHostToDevice);

	int numBlocks = 1;
	dim3 threadsPerBlock(N, N);
	MatAdd << <numBlocks, threadsPerBlock >> >(pA, pB, pC);

	cudaMemcpy(C, pC, (N*N)*sizeof(int), cudaMemcpyDeviceToHost);




	/*int **R = new int*[this->size];
	for (int i = 0; i < this->size; i++)
	R[i] = new int[this->size];

	int **dev_M = (int **)malloc(this->size * sizeof(int*));
	int **dev_R = (int **)malloc(this->size * sizeof(int*));

	for (int i = 0; i < this->size; i++){
	dev_M[i] = (int *)malloc(this->size * sizeof(int));
	dev_R[i] = (int *)malloc(this->size * sizeof(int));
	}

	int size = (this->size * this->size) * sizeof(int);

	cudaMalloc((void**)dev_M, size);
	cudaMalloc((void**)dev_R, size);

	testMatrix << <1, this->size >> >(dev_M, dev_R);


	cudaMemcpy(R, dev_R, size, cudaMemcpyDeviceToHost);*/

}

int* Node::subMinRowsAndCorrect(int s, const int *row, const int size){


	int *dev_s = 0;
	int *dev_re = 0;
	int *arr = (int *)malloc(size * sizeof(int));
	cudaError_t cudaStatus;

	cudaStatus = cudaSetDevice(0);
	if (cudaStatus == cudaSuccess){
		cudaMalloc((void**)&dev_s, sizeof(int));
		cudaMalloc((void**)&dev_re, size * sizeof(int));

		cudaMemcpy(dev_s, &s, sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(dev_re, row, size * sizeof(int), cudaMemcpyHostToDevice);

		//testKernel <<<1, size >>>(dev_s, dev_re);

		modifyArrayKernel << <1, size >> >(dev_s, dev_re);

		cudaDeviceSynchronize();

		int *c = (int *)malloc(sizeof(int));


		cudaMemcpy(c, dev_s, sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(arr, dev_re, size * sizeof(int), cudaMemcpyDeviceToHost);

		cudaFree(dev_s);
		cudaFree(dev_re);
	}

	return arr;
}

__global__ void minValue(int *source, int *val){
	__shared__ int temp[1];

	int currentValue = source[threadIdx.x];
	if (currentValue > -1 && currentValue < *val){
		temp[0] = currentValue;
	}

	__syncthreads();

	*val = temp[0];
}

__device__ int minVal = 100;

__device__ int blockChange = 0;

__global__ void matrixRowCorrect(int *arr, int *tmp, int *s){

	int id = threadIdx.x + blockIdx.x * blockDim.x;

	if (blockIdx.x != blockChange){
		blockChange = blockIdx.x;
		minVal = 100;
	}

	int currentValue = arr[id];

	if (currentValue < minVal){
		tmp[blockIdx.x] = currentValue;
		minVal = currentValue;
	}

}

__global__ void arrayReduce(int *m, int *ms){
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	if (m[id] > -1)
		m[id] = m[id] - ms[blockIdx.x];
}

__global__ void cleanCopy(int *S, int *D){
	D[threadIdx.x] = S[threadIdx.x];
}

int * Node::cudaCleanCopy(int *source){

	int matrixSize = this->size * this->size;
	int allocatedSize = matrixSize * sizeof(int);

	//int *d_source = (int*)malloc(allocatedSize);
	//int *d_destiny = (int*)malloc(allocatedSize);

	////memcpy(d_source, source, allocatedSize);

	//cudaMalloc((void**)&d_source, allocatedSize);
	//cudaMalloc((void**)&d_destiny, allocatedSize);

	//cudaMemcpy(d_source, source, allocatedSize, cudaMemcpyHostToDevice);

	//cleanCopy << <1, matrixSize >> >(d_source, d_destiny);

	this->M = (int*)malloc(allocatedSize);

	//cudaMemcpy(this->M, d_destiny, allocatedSize, cudaMemcpyDeviceToHost);

	//cudaFree(d_destiny);
	//cudaFree(d_source);
	//free(source);

	memcpy(this->M, source, allocatedSize);

	return this->M;
}

