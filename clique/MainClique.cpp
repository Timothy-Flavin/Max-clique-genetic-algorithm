#include<iostream>
#include<string>
#include<cstdlib>
#include<ctime>
#include<vector>
//using namespace std;

struct connection {
	int v1;
	int v2;
};

struct vertex {
	int value;
	std::vector<connection> connections;
	int indexInGraph;
};

struct genome {
	int genomeLength;
	vertex*vertices;
	genome() {

	}
	genome(int gl, vertex*genePool, int genePoolLength) {
		genomeLength = gl;
		int*usedNums = new int[gl];
		int usedLength = 0;
		bool used = false;
		vertices = new vertex[genomeLength];
		for (int i = 0; i < genomeLength; i++) {
			do {
				used = false;
				int vNum= rand() % genePoolLength;
				//cout << vNum << "vnum\n";
				vertices[i] = genePool[vNum];
				usedNums[i] = vNum;
				for (int j = 0; j < usedLength; j++) {
					if (usedNums[j] == vNum) {
						used = true;
					}
				}
			} while (used == true);
			usedLength++;
		}
	}
};

bool containsConnection(connection, connection*, int);
void printVertex(vertex, int);
void setConnections(int, int, connection*);
void makeGraph(int, int, connection*, vertex*);
double getFitness(genome*, int);
double* makeFitness(genome**, short, int);
double fitnessTotal(double*, int);
genome* selectClique(double*, double, int, genome**);
int getMostFit(double*, int);
void crossover(genome*, genome*, int);
void crossSwap(vertex*, vertex*, int, int);
genome* copyOf(genome*);
void mutate(genome*, int, float, int, vertex*);

int main() {
	srand(time(NULL));
	int numVertices = 200, maxConnections=0, popSize = 300, genomeLength = 3, maxGenerations = 200;
	float density = 0.3, mutationRate = 0.05;
	genome**population;
	population = new genome*[popSize];
	vertex* graph = new vertex[numVertices];

	for (int i = numVertices-1; i > 0; i--) {
		maxConnections += i;
	}

	int numConnections = static_cast<int>(density * maxConnections);
	std::cout << "Max connections" << maxConnections << std::endl;
	//cout << "Num connections" << numConnections << endl;
	//POPULATION PARAMETERS
	
	connection * connections = new connection[numConnections];
	setConnections(numConnections, numVertices, connections);
	makeGraph(numVertices, numConnections, connections, graph);

	int biggestClique = 2;
	int totalGens = 0;
	genome* maxClique = new genome(genomeLength, graph, numVertices);
	do {
		for (int i = 0; i < popSize; i++) {
			population[i] = new genome(genomeLength, graph, numVertices);
			//cout <<"before: " << i << ", " << getFitness(population[i], genomeLength) << endl;
		}

		double* fitness = makeFitness(population, popSize, genomeLength);
		double fitTotal = fitnessTotal(fitness, popSize);
		int winner;
		genome* best = copyOf(population[getMostFit(fitness, popSize)]);

		while(getFitness(best, genomeLength) < 1 && totalGens<maxGenerations) {
			genome ** gen2 = new genome*[popSize];
			for (int i = 0; i < popSize; i += 2) {
				gen2[i] = selectClique(fitness, fitTotal, popSize, population);
				gen2[i + 1] = selectClique(fitness, fitTotal, popSize, population);
			}

			for (int i = 0; i < popSize; i += 2) {
				crossover(gen2[i], gen2[i + 1], genomeLength);
				mutate(gen2[i], genomeLength, mutationRate, numVertices, graph);
				mutate(gen2[i + 1], genomeLength, mutationRate, numVertices, graph);
			}
			delete[] fitness;
			for (int i = 0; i < popSize; i++) {
				delete[] population[i]->vertices;
				delete population[i];

			}
			population = gen2;
			fitness = makeFitness(population, popSize, genomeLength);
			fitTotal = fitnessTotal(fitness, popSize);
			winner = getMostFit(fitness, popSize);
			std::cout << "winner: " << winner;
			genome*bestOfGen = copyOf(population[winner]);
			std::cout <<", with fitness: " << getFitness(bestOfGen, genomeLength) <<" at genome length: " << genomeLength<< std::endl;
			for(int q = 0; q < genomeLength; q++)
				std::cout << bestOfGen->vertices[q].indexInGraph << ", ";
			std::cout << std::endl;
			if (getFitness(best, best->genomeLength) < getFitness(bestOfGen, bestOfGen->genomeLength)) {
				best = bestOfGen;
			}
			else {
				delete[] bestOfGen->vertices;
			}
			totalGens++;
			//delete gen2;
		}
		std::cout <<"END OF CLIQUE------------------- "<< getFitness(best, genomeLength) << "----------------\n";
		std::cout << "total gens: " << totalGens <<std::endl;
		if (getFitness(best, best->genomeLength) == 1) {
			biggestClique = genomeLength;
			maxClique = copyOf(best);
			genomeLength++;
		}
		for (int i = 0; i < popSize; i++) {
			delete[] population[i]->vertices;
			delete population[i];

		}
		//delete population;
	} while (totalGens < maxGenerations);

	std::cout << "Max clique found: " << maxClique->genomeLength;
	/*
	genome testG(5, graph, numVertices);
	for (int i = 0; i < testG.genomeLength; i++) {
		cout << "Genome length: " << testG.genomeLength << ", vertex " << i << ": " << testG.vertices[i].value << endl;
		printVertex(testG.vertices[i], i);
		cout << "Compare stuff\n";
		printVertex(graph[testG.vertices[i].indexInGraph], testG.vertices[i].indexInGraph);
	}
	*/
	std::cin.get();
	std::cin.get();
	return 0;
}

void makeGraph(int numVertices, int numConnections, connection* connections, vertex* graph) {
	for (int i = 0; i < numVertices; i++) {
		graph[i].value = 0;
		graph[i].indexInGraph = i;
		for (int j = 0; j < numConnections; j++) {
			if (connections[j].v1 == i || connections[j].v2 == i) {
				graph[i].value++;
				graph[i].connections.push_back(connections[j]);
			}
		}
		//printVertex(graph[i],i);
	}
}
void setConnections(int numConnections, int numVertices, connection* connections) {
	for (int i = 0; i < numConnections; i++) {
		do {
			connections[i].v1 = rand() % numVertices;
			do {
				connections[i].v2 = rand() % numVertices;
			} while (connections[i].v2 == connections[i].v1);
		} while (containsConnection(connections[i], connections, i));
		//cout << "Connection: " << i << ", v1: " << connections[i].v1 << ", v2: " << connections[i].v2<<endl;
	}
}
double getFitness(genome*clique, int genomeLength) {
	int max = genomeLength * (genomeLength - 1);
	//cout << "max fit num: "<< max << endl;
	double fitNum = 0;
	for (int i = 0; i <genomeLength; i++) {
		for (int j = 0; j < clique->vertices[i].value;j++) {
			for (int k = 0; k < genomeLength; k++) {
				if (k != i) {
					if (clique->vertices[i].connections[j].v1 == clique->vertices[k].indexInGraph 
						|| clique->vertices[i].connections[j].v2 == clique->vertices[k].indexInGraph)
						fitNum++;
				}
			}
		}
	}
	fitNum = (static_cast<double>(fitNum) / max)*(static_cast<double>(fitNum) / max);
	//fitNum = 1.0 / (1 << static_cast<int>((max - static_cast<int>(fitNum)) / 2.0 ));
	return fitNum;
}
bool containsConnection(connection cn, connection* cns, int cnsLength) {
	for (int i = 0; i < cnsLength; i++) {
		if ((cn.v1 == cns[i].v1 && cn.v2 == cns[i].v2) || (cn.v2 == cns[i].v1 && cn.v1 == cns[i].v2))
			return true;
	}
	return false;
}
void printVertex(vertex v, int vnum) {
	std::cout << "Vertex " << vnum << " connections:\n";
	for(int i = 0; i < v.value; i++) {
		std::cout << "Connection " << i << " v1 " << v.connections[i].v1 <<", v2 "<<v.connections[i].v2<<std::endl;
	}
}
double* makeFitness(genome** population, short popSize, int genomeLength) {
	double* fitnessArray = new  double[popSize];
	for (int i = 0; i < popSize; i++) {
		fitnessArray[i] = getFitness(population[i], genomeLength);
		//fitnessArray[i] = 1.0 / (pow(pow(abs(population[i].x - goal.x),2) + pow(abs(population[i].y - goal.y),2),0.5));
		//will never be 1/0 because this funcion will not be called if a mouse made it
	}
	return fitnessArray;
}
double fitnessTotal(double* fitnessArray, int popSize) {
	double total = 0;
	for (int i = 0; i < popSize; i++)
		total += fitnessArray[i];
	return total;
}
genome* selectClique(double* fitnessArray, double fitnessTotal, int popSize, genome** population) {
	double randSelect = static_cast <double> (rand()) / (static_cast <double> ((RAND_MAX) / fitnessTotal)); //makes a random number from 0 to total fitness
	double errorCheck = randSelect;
	while (randSelect >= fitnessTotal || randSelect <=0) {
		randSelect = static_cast <double> (rand()) / (static_cast <double> ((RAND_MAX) / fitnessTotal));
		//cout << "shit boi";
	}
	int i = 0;
	while (randSelect >= fitnessArray[i]) { ////////////////////////////////////
		randSelect -= fitnessArray[i];
		i++;
	}
	bool error = false;
#ifdef DEBUG
	bool error = false;

	if (i == popSize) {
		error = true;
		cout << "i = popSize";
		i--;
		cin.get();
	}
	if (i < 0) {
		error = true;
		cout << "i < 0";
	}
	if (population[i] == NULL) {
		error = true;
		cout << "population i = null\n";
	}
	else if (population[i]->genomeLength == NULL) {
		error = true;
		cout << "population genome length = null\n";
	}
	if (error) {
		cout << "i: " << i << endl << "popSize: " << popSize << endl;
		cout << "fitnessTotal: " << fitnessTotal<<endl;
		cout << "og number" << errorCheck<<endl;
	}
#endif
	if (population[i] == NULL) {
		error = true;
		std::cout << "population i = null\n";
	}
	else if (population[i]->genomeLength == NULL) {
		std::cout << "hwat te fuk" << std::endl;
	}
	if (error) {
		std::cout << "i: " << i << std::endl << "popSize: " << popSize << std::endl;
		std::cout << "fitnessTotal: " << fitnessTotal << std::endl;
		std::cout << "og number" << errorCheck << std::endl;
	}
	return copyOf(population[i]);
}
genome* copyOf(genome* og) {
	genome* copy = new genome();
	copy->genomeLength = og->genomeLength;
	copy->vertices = new vertex[copy->genomeLength];
	for (int i = 0; i < copy->genomeLength; i++) {
		copy->vertices[i] = og->vertices[i];
	}
	return copy;
}
int getMostFit(double*fitArray, int popSize) {
	int mostFit = 0;
	double current = 0;
	for (int i = 0; i < popSize; i++) {
		//cout << "in fit tester, fit [" << i << "] = "<< fitArray[i]<<endl;
		if (fitArray[i] > current) {
			mostFit = i;
			current = fitArray[i];
		}
	}
	return mostFit;
}
void crossover(genome* parent1, genome* parent2, int genomeLength) {
	//vertex* kid1 = new vertex[parent1->genomeLength];
	//vertex* kid2 = new vertex[parent1->genomeLength];
	int crossPoint = rand() % genomeLength;

	for (int i = crossPoint; i < genomeLength; i++) {
		crossSwap(parent1->vertices, parent2->vertices, i, genomeLength);
	}

}
void crossSwap(vertex* list1, vertex* list2, int swapP, int genomeLength) {
#ifdef DEBUG
	for (int i = 0; i < genomeLength; i++) {
		cout << "starting array1: " << list1[i].indexInGraph << ",";
}
	cout << endl;
	for (int i = 0; i < genomeLength; i++) {
		cout << "starting array2: " << list2[i].indexInGraph << ",";
	}
	cout << endl;
#endif // DEBUG
	vertex v1 = list1[swapP]; // the point in list 1 at the swap point
	vertex v2 = list2[swapP]; // the point in list 2 at the swap point
							  // the veretx that has the same value as v2 in list 1
	
	int v2inl1 = -1;          // index of v2 in l1
	int v1inl2 = -1;          // index of v1 in l2

	for (int i = 0; i < genomeLength; i++) {
		if (list1[i].indexInGraph == v2.indexInGraph) {
			v2inl1 = i;
		}
		if (list2[i].indexInGraph == v1.indexInGraph) {
			v1inl2 = i;
		}
	}
	if (v2inl1 == -1) {
		list1[swapP] = v2;
	}
	else {
		//cout << "v2in1" << v2inl1;
		list1[v2inl1] = v1;
		list1[swapP] = v2;
	}
	if (v1inl2 == -1) {
		list2[swapP] = v1;
	}
	else {
		list2[swapP] = v1;
		list2[v1inl2] = v2;
	}

#ifdef DEBUG
	for (int i = 0; i < genomeLength; i++) {
		cout << "starting array1: " << list1[i].indexInGraph << ",";
	}
	cout << endl;
	for (int i = 0; i < genomeLength; i++) {
		cout << "starting array2: " << list2[i].indexInGraph << ",";
	}
#endif
}
void mutate(genome* gen, int genomeLength, float mutationRate, int graphSize, vertex* graph) {
	//cout << "new one ---------------------------------------------------------" << endl;
	for (int i = 0; i < genomeLength; i++) {
		//cout << gen->vertices[i].indexInGraph<<endl;
		if (rand() % 100 < static_cast<int>(mutationRate * 100)) {
			bool used;
			int vNum;
			do {
				used = false;
				vNum = rand() % graphSize;
				//cout << vNum << "vnum\n";
				for (int j = 0; j < genomeLength; j++) {
					if (gen->vertices[j].indexInGraph == vNum) {
						used = true;
					}
				}
			} while (used == true);
			gen->vertices[i] = graph[vNum];
			//cout << "changed to: " << gen->vertices[i].indexInGraph << endl;
		}
	}
}