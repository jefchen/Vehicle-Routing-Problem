/*
 * main.cpp
 *
 *  Created on: Jun 12, 2012
 *      Author: jeffreychen_55
 */

using namespace std;

#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <math.h>


const int FINAL_EVAL = 1285;
const int NUMSITES = 20;

// simulated annealing
const float INITIAL_TEMP = 100;
const int ITERATION = 200;
const float FINAL_TEMP = 0.01;
const float DECREMENT_RULE=0.95;

// genetic algorithm
const float CROSSOVER_PROB = 0.7;
const float MUTATION_PROB = 0.2;
const int POPULATION = 500;
const int GENERATIONS = 100;

struct Roulette_piece {
	float end;
	int index;
};

/*returns a command line parameter if it exist*/
char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

/*load the initial permutation p into its array*/
void load_initial_permutation(int arr[NUMSITES], int line_num)
{
    ifstream ifile("perm.txt");
    string line;
    int row = 0;
    if (ifile.is_open()){
        while(!ifile.eof())
        {
            getline(ifile,line);
            if (line_num == row)
            {

                int line_index = 0;
                int column = 0;
                for (;line_index<line.length();){
                    int pos = line.find(",", line_index);
                    if (pos != string::npos)
                    {
                        arr[column] = atoi((line.substr(line_index, pos-line_index)).c_str());
                        line_index = pos + 1;
                    }
                    else
                    {
                        arr[column] = atoi((line.substr(line_index, line.length() - line_index - 1)).c_str());
                        line_index = line.length();
                    }
                    column++;
                }
                return;
            }
            row++;

        }
    }
    cout<<"Fatal: no such perm entry" << endl;
    exit(1);
}

/*load the initial data matrices (flow and distance) from file into passed in array*/
void populate_array(int arr[NUMSITES][NUMSITES], string file){

    ifstream ifile(file.c_str());
    string line;
    int row = 0;
    if (ifile.is_open()){
        while(!ifile.eof()){
            getline(ifile,line);
            int line_index = 0;
            int column = 0;
            for (;line_index<line.length();){
                int pos = line.find(" ", line_index);
                if (pos != string::npos)
                {
                    arr[row][column] = atoi((line.substr(line_index, pos-line_index)).c_str());
                    line_index = pos + 1;
                }
                else
                {
                    arr[row][column] = atoi((line.substr(line_index, line.length() - line_index - 1)).c_str());
                    line_index = line.length();
                }
                column++;
            }
            row++;
        }
    }

    for (int i = 0;i<NUMSITES;i++){
        for (int j = 0; j< NUMSITES;j++)
            ;//arr[i][j] = 0;;
    }
}

void print(int p[NUMSITES]) {
	cout << "[";
	for (int i=0; i<NUMSITES; i++) {
		cout << p[i] << ",";
	}
	cout << "]";
}

int compute_result(int p[NUMSITES], int flow[NUMSITES][NUMSITES], int dist[NUMSITES][NUMSITES]){
    int sum = 0;
    for (int i = 0;i<NUMSITES;i++){
        for (int j = i; j< NUMSITES;j++)
        {
            sum += flow[p[i]-1][p[j]-1] * dist[i][j];
        }
    }
    return sum;
}

//swap two items in the permutation array
void swap(int p[NUMSITES], int i, int j)
{
    int temp = p[i];
    p[i] = p[j];
    p[j] = temp;
}

// generate random number between 0 and 1
float randomZeroAndOne()
{
	float scale=RAND_MAX+1.;
	float base=rand()/scale;
	float fine=rand()/scale;
	return base+fine/scale;
}

void simulated_annealing(int p[NUMSITES], int flows[NUMSITES][NUMSITES], int dists[NUMSITES][NUMSITES]) {

	// compute initial state valuation
	int previous_eval = compute_result(p,flows,dists);
	int current_eval = 10000;
	int iteration = 0;
	float current_temp = INITIAL_TEMP;
	int total_iteration = 0;

	while (current_temp > FINAL_TEMP) {
	//while (current_eval > FINAL_EVAL + 100) {
		while (iteration < ITERATION) {
			// get solution (randomly get a solution)
			int i = rand() % NUMSITES;
			int j = 0;
			do {
				j = rand() % NUMSITES;
			} while (j == i);

			swap(p,i,j);

			// calculate result
			current_eval = compute_result(p,flows,dists);
			if (current_eval - previous_eval >= 0) {
				// calculate probabilities
				float prob = randomZeroAndOne();
				float temp_prob = exp(-((float)(current_eval-previous_eval))/current_temp);
				if (prob > temp_prob) {
					// revert move
					swap(p,i,j);
				} else {
					// keep move and update evaluation
					previous_eval = current_eval;
				}
			} else {
				// keep move and update evaluation
				previous_eval = current_eval;
			}

			total_iteration++;
			iteration++;
			cout << "evaluation = " << current_eval << endl;
		}

		// decrease temperature
		current_temp *= DECREMENT_RULE;

		// reset iteration count
		iteration = 0;

	}

	cout << "Iteration Count = " << total_iteration << endl;
    cout << "Final Evaluation = " << current_eval << endl;
}

void crossover(int individual1[], int individual2[]) {
	// randomly select 2 points
	int a = rand() % NUMSITES;
	int b = rand() % NUMSITES;
	if (a > b) {
		// swap a and b
		int temp = a;
		a = b;
		b = temp;
	}

	// copy the stuff in between
	int child1[NUMSITES];
	int child2[NUMSITES];
	memset(&child1, 0, sizeof(child1));
	memset(&child2, 0, sizeof(child2));
	for (int j=0; j<NUMSITES; j++) {
		if (j >= a && j < b) {
			child1[j] = individual1[j];
			child2[j] = individual2[j];
		}
	}

	// copy the rest in order
	int combined1[NUMSITES*2];
	int combined2[NUMSITES*2];
	for (int i=0; i<NUMSITES; i++) {
		combined1[i] = individual1[i];
		combined1[NUMSITES+i] = individual1[i];
		combined2[i] = individual2[i];
		combined2[NUMSITES+i] = individual2[i];
	}
	// child1
	bool flag = false;
	int k = b;
	for (int i=b; i<b+NUMSITES; i++) {
		for (int j=a; j<b; j++) {
			if (combined2[i] == combined1[j]) {
				flag = true;
				break;
			}
		}
		if (!flag) {
			child1[k] = combined2[i];
			k++;
			if (k >= NUMSITES) {
				k = 0;
			}
		}
		flag = false;
	}
	// child2
	flag = false;
	k = b;
	for (int i=b; i<b+NUMSITES; i++) {
		for (int j=a; j<b; j++) {
			if (combined1[i] == combined2[j]) {
				flag = true;
				break;
			}
		}
		if (!flag) {
			child2[k] = combined1[i];
			k++;
			if (k >= NUMSITES) {
				k = 0;
			}
		}
		flag = false;
	}

	// copy array into individual 1 and individual 2
	for (int i=0; i<NUMSITES; i++) {
		individual1[i] = child1[i];
		individual2[i] = child2[i];
	}
}

void mutate(int individual[]) {
	// perform insert mutation
	int a = rand() % NUMSITES;
	int b = 0;
	do {
		b = rand() % NUMSITES;
	} while (a == b);
	if (a > b) {
		int temp = a;
		a = b;
		b = temp;
	}

	int i = b;
	int temp = individual[b];
	for (; i>a+1; i--) {
		individual[i] = individual[i-1];
	}
	individual[i] = temp;
}

void copy(int to[], int from[], int size) {
	for (int i=0; i<size; i++) {
		to[i] = from[i];
	}
}

void copy2D(int to[POPULATION][NUMSITES], int from[POPULATION][NUMSITES]) {
	for (int i=0; i<POPULATION; i++) {
		for (int j=0; j<NUMSITES; j++) {
			to[i][j] = from[i][j];
		}
	}
}

void shuffle(int p[NUMSITES]) {
	// shuffle
	for (int j=NUMSITES-1; j>0; j--) {
		int temp = rand() % (j+1);
		swap(p,temp,j);
	}
}

void genetic_algorithm(int p[NUMSITES], int flows[NUMSITES][NUMSITES], int dists[NUMSITES][NUMSITES]) {

	cout << "genetic algorithm started" << endl;

	// setup initial solutions (randomly pick solutions)
	int population[POPULATION][NUMSITES];
	for (int i=0; i<POPULATION; i++) {
		// copy
		for (int j = 0; j<NUMSITES; j++) {
			population[i][j] = p[j];
		}
		// shuffle
		shuffle(population[i]);
	}

	cout << "setup initial population" << endl;

	for (int genCount=0; genCount<GENERATIONS; genCount++) {
		// calculate fitness
		int fitness[POPULATION];
		float totalFitness = 0;
		for (int i=0; i<POPULATION; i++) {
			fitness[i] = compute_result(population[i], flows, dists);
			totalFitness += 2000/fitness[i];
		}

		// print best solution
		int best = INT_MAX;
		int bestIndividual = 0;
		for (int i=0; i<POPULATION; i++) {
			if (fitness[i] < best) {
				best = fitness[i];
				bestIndividual = i;
			}
		}
		cout << "best individual: fitness=" << best << ", solution=";
		print(population[bestIndividual]);
		cout << endl;

		// setup roulette
		Roulette_piece roulette[POPULATION];
		float count = 0;
		for (int i=0; i<POPULATION; i++) {
			count += (float)2000/(float)fitness[i]/totalFitness;
			roulette[i].end = count;
			roulette[i].index = i;
		}

		// perform roulette parent selection
		int nextGen[POPULATION][NUMSITES];
		float temp = 0;
		for (int j=0; j<POPULATION; j++) {
			temp = randomZeroAndOne();
			for (int i=0; i<POPULATION; i++) {
				if (temp <= roulette[i].end) {
					copy(nextGen[j],population[roulette[i].index],NUMSITES);
				}
			}
		}

		// cross over and mutate
		for (int i=0; i<POPULATION; i+=2) {
			// perform oder-1 crossover
			if (randomZeroAndOne() <= CROSSOVER_PROB) {
				crossover(nextGen[i], nextGen[i+1]);
			}

			// perform mutation
			if (randomZeroAndOne() <= MUTATION_PROB) {
				mutate(nextGen[i]);
				mutate(nextGen[i+1]);
			}
		}

		// replace population with next generation
		bool eliteCopied = false;
		for (int i=0; i<POPULATION; i++) {
			if (i == bestIndividual && !eliteCopied) {
				eliteCopied = true;
			} else {
				for (int j=0; j<NUMSITES; j++) {
					population[i][j] = nextGen[i][j];
				}
			}
		}
	}

	cout << "genetic algorithm finished" << endl;
}

int main(int argc, char* argv[]) {
    int p[NUMSITES];
    int flows[NUMSITES][NUMSITES];
    int dists[NUMSITES][NUMSITES];

    memset( &flows, 0, sizeof(flows));
    memset( &dists, 0, sizeof(dists));
    memset( &p, 0, sizeof(p));

	// initialize random seed
	srand(time(NULL));

    // get initial permutation line number
    int perm = atoi(getCmdOption(argv,argv+argc,"-perm"));

	// setup problem
    load_initial_permutation(p, perm);
    populate_array(flows, "flow.txt");
    populate_array(dists, "distance.txt");

    // SA
    simulated_annealing(p, flows, dists);
    //genetic_algorithm(p, flows, dists);

    /*
    int p2[NUMSITES];
    copy(p2,p,NUMSITES);
    shuffle(p2);
    cout << "orignal" << endl;
    print(p);
    cout << endl;
    print(p2);
    cout << endl;
    crossover(p,p2);
    cout << "after" << endl;
    print(p);
    cout << endl;
    print(p2);
    cout << endl;
    */

	return 0;
}
