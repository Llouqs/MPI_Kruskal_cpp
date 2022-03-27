#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <inttypes.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#define VERBOSE 0

using namespace std; 

typedef struct dsNodeType {
	struct dsNodeType* parent;
	int rank;
} dsNode;

// Set that contains all vertices as dsNodes.
dsNode *dsSet;

// MakeSet operation.
inline void dsMakeSet(int nVerts) {
	// Callocates (init to 0) room for all vertices.
	// parent = NULL and rank = 0.
	dsSet = new dsNode[nVerts];
}

// Find operation.
inline dsNode* dsFind(dsNode *n) {
	if (n->parent == NULL) return n;
	n->parent = dsFind(n->parent);
	return n->parent;
}

// Union operation.
// Assuming that n and m belong to different sets.
inline void dsUnion(dsNode* n, dsNode* m) {
	if (n->rank < m->rank) {
		n->parent = m;
	}
	else if (n->rank > m->rank) {
		m->parent = n;
	}
	else {
		n->parent = m;
		n->rank += 1;
	}
}
struct edge {
	int u;
	int v;
	int weight;

	edge(){}

	edge(int _u, int _v, int _weight):u(_u), v(_v), weight(_weight){}
};

int compare(const void* a, const void* b)
{
	const edge* _a = reinterpret_cast<const edge*>(a);
	const edge* _b = reinterpret_cast<const edge*>(b);
	return _a->weight < _b->weight ? -1 : _a->weight > _b->weight ? 1 : 0;
}

MPI_Datatype mpiEdge;
int mpiNp;
int mpiRank;
int nVerts; 
int nEdges;
edge *edges; // Contains the edges of a processor.
edge *edges_sort_buffer; // Used for parallel quicksorting of 'edges' for performance.
int edgeCount;

ifstream digits;

edge *mst; // Contains the MST edges of a processor. Also used as send/recv MPI buffer.
int mstEdgeCount;
int mstLength;

// Timing variables.
double commTime;
double procTime;
double parseTime;


// Checks if (l  <=  x  <  r)
inline int in(int x, int l, int r) {
	return ((x >= l) && (x < r));
}

// Finalizes MPI and frees buffers.
inline void finalize() {
	MPI_Type_free(&mpiEdge);
	MPI_Finalize();
	free(edges);
	free(mst);
	free(dsSet);
}

// Force kill due to initialization error.
inline void die(const char *msg) {
	printf("%d: %s\n", mpiRank, msg);
	MPI_Type_free(&mpiEdge);
	MPI_Finalize();
	exit(0);
}

// Initializes MPI and mpiEdge datatype.
inline void initMPI(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpiNp);
	MPI_Type_contiguous(3, MPI_INT, &mpiEdge);
	MPI_Type_commit(&mpiEdge);
}

bool read_num(ifstream &stream, int &number)
{
	char c = 0;
	std::string buf;
	while (true) {
		stream.read(&c, 1);
		if (c == ' ' || stream.eof() || c == '\n' || c == '\t') {
			if (!buf.empty()) {
				number = atoi(buf.c_str());
				return true;
			}
		}
		else {
			buf += c;
		}
	}

	return false;
}

// Parses the edges that correspond to the current processor.
inline void parseEdges() {
	digits.open("correct4.txt", ios_base::in);
	int from = -1, tmp = 0;
	while (!digits.eof())
	{
		if (read_num(digits, tmp))
			if (from != tmp) {
				from = tmp;
				nVerts++;
			}
		read_num(digits, tmp);
		read_num(digits, tmp);
		nEdges++;
	}
	digits.close();

	int vPerProc = nVerts / mpiNp; // Vertices per processor.
	int firstVert = mpiRank * vPerProc; // First vertex index.
	int lastVert = (mpiRank + 1) * vPerProc; // Last vertex index.

	// Last processor case.
	if (mpiRank == mpiNp - 1) {
		lastVert += nVerts % mpiNp;
		vPerProc += nVerts % mpiNp;
	}

	// Mem allocation.
	edges = new edge[vPerProc * nVerts];
	edges_sort_buffer = new edge[vPerProc * nVerts];;

	// Allocating x2 because of merging. The first part of mst buffer contains the
	// MST edges and the second part is used as a receive MPI buffer.
	mst = new edge[2 * (nVerts - 1)];

	// Edge parsing.
	edgeCount = 0;
	edge e;
	digits.open("correct4.txt", ios_base::in);
	if (!digits.is_open())
		return;
	int to, weight;
	while (!digits.eof())
	{
		if (read_num(digits, from))
			if (read_num(digits, to))
				if (read_num(digits, weight)) {
					edge* e = new edge(from, to, weight);
					/*e.u = from;
					e.v = to;
					e.weight = weight;
					edges[edgeCount].u = from;
					edges[edgeCount].v = to;
					edges[edgeCount].weight = weight;
					edgeCount++;
					//cout << e.u << " - " << e.v << " : " << e.weight << endl;
					*/
					edges[edgeCount] = *e;
					//cout << edges[edgeCount].u << " - " << edges[edgeCount].v << " : " << edges[edgeCount].weight << endl;
				}
	}
	digits.close();

	//for (i = 0; i < nEdges; i++) {
	//	fread(&e, sizeof(edge), 1, input);
	//	if (in(e.u, firstVert, lastVert) || in(e.v, firstVert, lastVert)) {
	//		edges[edgeCount++] = e;
	//		cout << "u: " << e.u << " v: " << e.v << " w: " << e.weight << endl;
	//	}
	//}
	//fclose(input);
	parseTime = MPI_Wtime() - parseTime;
	// Waiting for all processors to complete parsing.
	MPI_Barrier(MPI_COMM_WORLD);
}

// Calculates the MST.
inline void calculateMst() {
	// Parallel quicksort.
	qsort(edges, 6, sizeof(edge), compare);

	free(dsSet);
	dsMakeSet(nVerts);

	// Looping trhough all edges in increasing weight order.
	mstEdgeCount = 0;
	mstLength = 0;
	int i = 0;
	for (; i < edgeCount; i++) {
		// Find parent sets of the nodes.
		dsNode *vParent = dsFind(&dsSet[edges[i].v]);
		dsNode *uParent = dsFind(&dsSet[edges[i].u]);
		// If they are from two different sets, merge them.
		if (vParent != uParent) {
			mst[mstEdgeCount++] = edges[i];
			mstLength += edges[i].weight;
			dsUnion(vParent, uParent);
		}
	}
}

int main(int argc, char **argv) {
	int rank, size;
	// Initialization and input file parsing.
	initMPI(argc, argv);
	parseTime = MPI_Wtime();
	parseEdges();
	// Processing.
	double tmpTime;
	procTime = MPI_Wtime();

	//cout << sizeof(edges) / sizeof(edge) << endl;
	// In case of 1 processor, act as serial application.
	if (mpiNp == 1) {
		/*for (int i = 0; i < 15; i++)
		{
			cout << "u: " << edges[i].u << " v: " << edges[i].v << " w: " << edges[i].weight << endl;
		}*/
		calculateMst();
	}
	else {
		int processors = mpiNp;
		int pow2 = 1; // Used to find out where to send/recv.
		MPI_Status mpiStatus;
		int recvEdgeCount;
		while (processors > 1) {
			// Calculate the local MST using 'edges' buffer.
			calculateMst();

			// Communication part.
			if ((mpiRank / pow2) % 2 != 0) {
				// Send from the first half of 'mst'.
				MPI_Send(mst, mstEdgeCount, mpiEdge, (mpiRank - pow2), 0, MPI_COMM_WORLD);
				break; // Processor did his job and can now exit.
			}
			else {
				// Receive into the second half of 'mst'.
				tmpTime = MPI_Wtime();
				MPI_Recv(mst + mstEdgeCount, nVerts - 1, mpiEdge, (mpiRank + pow2), 0, MPI_COMM_WORLD, &mpiStatus);
				MPI_Get_count(&mpiStatus, mpiEdge, &recvEdgeCount);
				commTime += MPI_Wtime() - tmpTime; // Communication time is equal to the waiting-to-recv time.
			}
			// Pointer swap between 'edges' and 'mst'.
			// 'edges' will be reused to calculate a new local MST.
			// 'mst' will contain that new local MST.
			edge *tmp = edges;
			edges = mst;
			mst = tmp;
			edgeCount = mstEdgeCount + recvEdgeCount; // New edge count after merging.

			processors /= 2;
			pow2 *= 2;
		}

		// Only the root will execute this last calculation.
		if (mpiRank == 0) calculateMst();
	}
	procTime = (MPI_Wtime() - procTime) - commTime;


	if (VERBOSE)
		printf("%d: Parse time: %.3fs\n%d: Comm time: %.3fs\n%d: Proc time: %.3fs\n"
			, mpiRank, parseTime, mpiRank, commTime, mpiRank, procTime);
	if (mpiRank == 0) {
		printf("MST length: %d\n", mstLength);
		printf("Total Time %.3fs\n", parseTime + commTime + procTime);
		printf("Without I/O %.3fs\n", commTime + procTime);
	}
	finalize();
	return 0;
}