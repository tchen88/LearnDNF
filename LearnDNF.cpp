/* A brief overview of the upcoming code:
   
   Let n be the number of variables. Then each term in the DNF can be represented by
   an n-bit string of 0's and 1's, where a 1 in the ith position corresponds to the
   presence of the variable x_i in the term. For example, if n=7 then x_1x_2x_5
   translates to 1100100 (after being translated to x_0x_1x_4 first since list indices
   start from 0 instead of 1). An entire monotone DNF expression is represented by
   multiple n-bit strings delimited by a whitespace character. We can thus read a
   monotone DNF expression and easily find the corresponding truth table.
   
   The hypothesis is represented by a vector of structs as follows: For each term
   F-hat * chi_bitstring, we create an entry with key "bitstring" and value F-hat.
   This is also the current format of the output.
*/

#include <cstdio>
#include <iostream>
#include <sstream>
#include <cmath>
#include <cstring>
#include <string>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <vector>
#include <climits>

using namespace std;

int counter = 0; // Global debug variable

// Constants
const int MAX_EXAMPLES = (int)pow(2.0, 30.0);
double TARGET_ERROR = 0.05;
int NUM_TERMS;
int NUM_VARS;
int VARS_PER_TERM;
int SEED_NUM = 0;
int TIME_LIMIT = 1500;
const int INVALID = -1;
const int UNIDIST = 0;
const int POSCORR = 1;
const int NEGCORR = 2;
const int SUPERPOSCORR = 3;
const int SUPERNEGCORR = 4;
const int DIFFERENCE = 0;
const int RATIO = 1;
const int CONSTANT = 0;
const int PROBRATIO = 1;
const int FLENGTH = 2;
vector<int> SEEDS;

struct funStruct {
	string index;
	int value;
};

struct doubleFunStruct {
	string index;
	double value;
};

struct indexStruct {
	int index;
	int value;
};

bool compareIndexStruct(indexStruct a, indexStruct b) { return (a.value < b.value); }
void findAllBitstrings(int length, vector<string>& bitstrings);
double findFourierCoef(vector<funStruct>& targetFunction, string index, vector<string>& allBitstrings);
int parityFunction(string index, string example);
double error(vector<funStruct>& targetFunction, vector<funStruct>& hypothesis);
int hypothesisOutput(vector<funStruct>& hypothesis, string example);
void printInfo(int loopIter, double currentError, double minError, vector<funStruct>& currentHypothesis, bool isNew);
void findUpNeighbors(vector<string>& bitstrings, vector<string>& neighbors);
void findDownNeighbors(string bitstring, vector<string>& neighbors);
string upgraded(string bitstring, int index) { return bitstring.substr(0,index) + "1" + bitstring.substr(index+1,bitstring.size()-index-1); }
string downgraded(string bitstring, int index) { return bitstring.substr(0,index) + "0" + bitstring.substr(index+1,bitstring.size()-index-1); }
int parity(string bitstring);
int valueAt(vector<funStruct>& s, string idx);
double valueAt(vector<doubleFunStruct>& s, string idx);
bool valueAt(vector<string>& s, string idx);
double difference(vector<funStruct>& targetFunction, vector<funStruct>& hypothesis, vector<doubleFunStruct>& fourierCoef, vector<string>& allBitstrings, int diffMethod);
double expectedValue(vector<funStruct>& hypothesis, vector<string>& allBitstrings);

int main(int argc, char** argv) {
	string filename = "";
	double alpha;
	int diffMethod;
	int initHyp;

	// Build the random seeds
	ifstream seedFile("Seeds.txt");
	string seedLine;
	while (getline(seedFile, seedLine)) {
		istringstream iss(seedLine);
		do {
			string sub;
			iss >> sub;
			if (sub != "") {
				SEEDS.push_back(atoi(sub.c_str()));
			}
		} while (iss);
	}

	// Read the options
	for (int i = 1; i < argc-1; ++i) {
		if (strstr(argv[i], "-f") != NULL) { // Text file containing DNF expression, if it exists. Overrides -a
			filename = argv[++i];
		}
		else if (strstr(argv[i], "-a") != NULL) { // Correlation between terms in the DNF expression. n turns this option off.
			++i;
			if (strchr(argv[i], 'n') != NULL || strchr(argv[i], 'N') != NULL) {
				alpha = -1.0;
			}
			else {
				alpha = strtod(argv[i], NULL);
			}
		}
		else if (strstr(argv[i], "-m") != NULL) { // Difference method
			++i;
			if (strchr(argv[i], 'd') != NULL || strchr(argv[i], 'D') != NULL) {
				diffMethod = DIFFERENCE;
			}
			else if (strchr(argv[i], 'r') != NULL || strchr(argv[i], 'R') != NULL) {
				diffMethod = RATIO;
			}
			else {
				diffMethod = INVALID;
			}
		}
		else if (strstr(argv[i], "-h") != NULL) { // Initial hypothesis
			++i;
			if (strchr(argv[i], '1') != NULL) {
				initHyp = CONSTANT;
			}
			else if (strchr(argv[i], 'p') != NULL || strchr(argv[i], 'P') != NULL) {
				initHyp = PROBRATIO;
			}
			else if (strchr(argv[i], 't') != NULL || strchr(argv[i], 'T') != NULL) {
				initHyp = FLENGTH;
			}
			else {
				initHyp = INVALID;
			}
		}
		else if (strstr(argv[i], "-e") != NULL) { // Target error
			TARGET_ERROR = strtod(argv[++i], NULL);
		}
		else if (strstr(argv[i], "-s") != NULL) { // Seed index
			SEED_NUM = strtol(argv[++i], NULL, 10);
		}
		else if (strstr(argv[i], "-t") != NULL) { // Number of terms
			NUM_TERMS = strtol(argv[++i], NULL, 10);
		}
		else if (strstr(argv[i], "-v") != NULL) { // Number of variables
			NUM_VARS = strtol(argv[++i], NULL, 10);
		}
		else if (strstr(argv[i], "-z") != NULL) { // Number of variables per term
			VARS_PER_TERM = strtol(argv[++i], NULL, 10);
		}
		else if (strstr(argv[i], "-o") != NULL) { // Number of variables per term
			TIME_LIMIT = strtol(argv[++i], NULL, 10);
		}
	}
	

	int dnf[NUM_TERMS][VARS_PER_TERM];
	if (filename != "") {
		ifstream inputFile(filename.c_str());
		string line;
		int lineNumber = 0;
		while (getline(inputFile, line)) {
			int varNumber = 0;
			for (int i = 0; i < line.size(); ++i) {
				if (line.at(i) != '0') {
					dnf[lineNumber][varNumber] = i;
					++varNumber;
				}
			}
			++lineNumber;
		}
	}

	else {
		// Create the DNF with the given random seed
		srand(SEEDS[SEED_NUM]);
		for (int i = 0; i < NUM_TERMS; ++i) {
			// Case 1: Correlation option is turned off. If the correlation option is on, we still need to go here for the first term
			if (i == 0 || alpha < 0 || alpha > 1) {
				// Determine the variables in the term by giving each variable a random number and picking the top s variables
				indexStruct randValues[NUM_VARS];
				for (int j = 0; j < NUM_VARS; ++j) {
					randValues[j].index = j;
					randValues[j].value = rand();
				}
				sort(randValues, randValues + NUM_VARS, compareIndexStruct);
				for (int j = 0; j < VARS_PER_TERM; ++j) {
					dnf[i][j] = randValues[j].index;
				}
			}

			// Case 2: Correlation option is turned on
			else {
				for (int j = 0; j < VARS_PER_TERM; ++j) {
					// Each variable can only be use once per term. Also, we mark the variable from the previous term, as well as whether it can be used
					vector<int> validChoices;
					bool validAlpha = false;
					for (int k = 0; k < NUM_VARS; ++k) {
						bool isAlphaVar = (dnf[i-1][j] == k);
						bool stillValid = true;
						for (int l = 0; l < j; ++l) {
							if (dnf[i][l] == k) {
								stillValid = false;
								break;
							}
						}
						if (stillValid) {
							if (isAlphaVar) {
								validAlpha = true;
							}
							else {
								validChoices.push_back(k);
							}
						}
					}

					double r = rand() * 1.0 / RAND_MAX;
					if (validAlpha) {
						if (r <= alpha) {
							dnf[i][j] = dnf[i-1][j];
						}
						else {
							dnf[i][j] = validChoices[(int)(validChoices.size() * (r-alpha) / (1-alpha))];
						}
					}
					else {
						dnf[i][j] = validChoices[(int)(validChoices.size() * r)];
					}
				}
			}
		}
	}

	// Create the truth table
	vector<string> allBitstrings;
	findAllBitstrings(NUM_VARS, allBitstrings);
	vector<funStruct> targetFunction;
	for (int i=0; i < allBitstrings.size(); ++i) {
		targetFunction.push_back(funStruct());
		targetFunction[i].index = allBitstrings[i];
		targetFunction[i].value = -1; // f(x) evaluates to False by default. Once we find a term that evaluates to True, f(x) flips to True
		for (int k = 0; k < NUM_TERMS; ++k) {
			bool termValue = true; // Each term evaluates to True by default. Once we find a false variable, the term flips to False
			for (int j = 0; j < VARS_PER_TERM; ++j) {
				if (allBitstrings[i].at(dnf[k][j]) == '0') {
					termValue = false;
					break;
				}
			}
			if (termValue) {
				targetFunction[i].value = 1;
				break;
			}
		}
	}

	// Start the clock
	clock_t t = clock();

	// Create initial hypothesis values
	vector<doubleFunStruct> fourierCoef;
	fourierCoef.push_back(doubleFunStruct());
	fourierCoef[0].index = string(NUM_VARS, '0'); // String of n 0's
	fourierCoef[0].value = findFourierCoef(targetFunction, string(NUM_VARS, '0'), allBitstrings);
	vector<string> ourSet; // S
	ourSet.push_back(string(NUM_VARS, '0'));
	vector<funStruct> currentHypothesis;
	currentHypothesis.push_back(funStruct());
	currentHypothesis[0].index = string(NUM_VARS, '0');
	currentHypothesis[0].value = 1;
	double currentError = error(targetFunction, currentHypothesis);
	if (initHyp == FLENGTH) {
		currentHypothesis[0].value = 2*NUM_TERMS - 1;
	}
	else if (initHyp == PROBRATIO) {
		currentHypothesis[0].value = ((int) ((1-currentError)/currentError) / 2) * 2 + 1;
	}

	// Start looping to minimize the error
	double minError = currentError;
	int loopIter = 0;
	printInfo(loopIter, currentError, minError, currentHypothesis, true);
	while (currentError > TARGET_ERROR) {
		vector<string> neighbors;
		findUpNeighbors(ourSet, neighbors);
		
		double minDiffValue = pow(2.0, 30.0); // Keep running track of the minimum difference value and it's index
		string indexOfMDV;

		// Case 1: hypotheses with new index
		for (int i = 1; i < ourSet.size(); ++i) { // i starts at 1 because we don't want 0^n as a valid hypothesis
			if (parity(ourSet[i]) == 1) {
				currentHypothesis[i].value -= 2;
			}
			else {
				currentHypothesis[i].value += 2;
			}
			// Make sure our new hypothesis is legal
			vector<string> downNeighbors;
			findDownNeighbors(ourSet[i], downNeighbors);
			bool isBadHypothesis = false;
			for (int j = 0; j < downNeighbors.size(); ++j) {
				if (downNeighbors[j].find("1") != string::npos && abs(valueAt(currentHypothesis, downNeighbors[j])) < abs(currentHypothesis[i].value)) {
					isBadHypothesis = true;
				}
			}
			if (isBadHypothesis) {
				// Restore original hypothesis and move on
				if (parity(ourSet[i]) == 1) {
					currentHypothesis[i].value += 2;
				}
				else {
					currentHypothesis[i].value -= 2;
				}
				continue;
			}
			// If the hypothesis is legal, find out if the difference value is a minimum
			double diffValue = difference(targetFunction, currentHypothesis, fourierCoef, allBitstrings, diffMethod);
			if (diffValue <= minDiffValue) {
				indexOfMDV = ourSet[i];
				minDiffValue = diffValue;
			}
			// Restore original hypothesis
			if (parity(ourSet[i]) == 1) {
				currentHypothesis[i].value += 2;
			}
			else {
				currentHypothesis[i].value -= 2;
			}
		}

		// Case 2: hypotheses with no new index
		for (int i = 0; i < neighbors.size(); ++i) {
			// Add the Fourier coefficient to our list if it hasn't been added already
			if (valueAt(fourierCoef, neighbors[i]) == (double) INT_MAX) {
				fourierCoef.push_back(doubleFunStruct());
				fourierCoef[fourierCoef.size()-1].index = neighbors[i];
				fourierCoef[fourierCoef.size()-1].value = findFourierCoef(targetFunction, neighbors[i], allBitstrings);
			}

			// Find new hypothesis
			currentHypothesis.push_back(funStruct());
			currentHypothesis[currentHypothesis.size()-1].index = neighbors[i];
			if (parity(neighbors[i]) == 1) {
				currentHypothesis[currentHypothesis.size()-1].value = -2;
			}
			else {
				currentHypothesis[currentHypothesis.size()-1].value = 2;
			}
			double diffValue = difference(targetFunction, currentHypothesis, fourierCoef, allBitstrings, diffMethod);
			if (diffValue <= minDiffValue) {
				indexOfMDV = neighbors[i];
				minDiffValue = diffValue;
			}
			// Restore original hypothesis
			currentHypothesis.pop_back();
		}

		// Update our hypothesis
		bool isNew = false;
		
		// New coefficient
		if (!valueAt(ourSet, indexOfMDV)) {
			ourSet.push_back(indexOfMDV);
			isNew = true;
			currentHypothesis.push_back(funStruct());
			currentHypothesis[currentHypothesis.size()-1].index = indexOfMDV;
			if (parity(indexOfMDV) == 1) {
				currentHypothesis[currentHypothesis.size()-1].value = -2;
			}
			else {
				currentHypothesis[currentHypothesis.size()-1].value = 2;
			}
		}

		// Existing coefficient
		else {
			for (int i = 0; i < currentHypothesis.size(); ++i) {
				if (currentHypothesis[i].index == indexOfMDV) {
					if (parity(indexOfMDV) == 1) {
						currentHypothesis[i].value -= 2;
					}
					else {
						currentHypothesis[i].value += 2;
					}
					break;
				}
			}
		}

		// Update loop variables
		++loopIter;
		currentError = error(targetFunction, currentHypothesis);
		if (currentError < minError) {
			minError = currentError;
		}
		printInfo(loopIter, currentError, minError, currentHypothesis, isNew);
		if (clock() - t > TIME_LIMIT * CLOCKS_PER_SEC) {
			cout << "Timed out: " << ((double)TIME_LIMIT) / 3600 << " hours" << endl;
			break;
		}
	}

	// Print out the DNF
	cout << "DNF: ";
	for (int i = 0; i < NUM_TERMS; ++i) {
		if (i != 0) {
			cout << " + ";
		}
		for (int j = 0; j < VARS_PER_TERM; ++j) {
			cout << 'x' << dnf[i][j];
		}
	}
	cout << endl;

	// Print some ending stats
	t = clock() - t;
	cout << "Running time: " << ((double)t) / CLOCKS_PER_SEC << "s" << endl;
	return 0;
}

void findAllBitstrings(int length, vector<string>& bitstrings) {
/* Returns a list of all bitstrings of the specified length */
	for (int i = 0; i < (int)(pow(2.0, length)); ++i) {
		int r = i;
		string result = "";
		int numOnes = 0;
		for (int j=0; j < length; ++j) {
			result += (r % 2 + '0');
			numOnes += r % 2;
			r = r / 2;
		}

		int multiplicity = 1;
		for (int j = 0; j < multiplicity; ++j) {
			bitstrings.push_back(result);
		}
	}
}

double findFourierCoef(vector<funStruct>& targetFunction, string index, vector<string>& allBitstrings) {
/* Returns the F-hat of the given bitstring */
	int result = 0;
	for (int i = 0; i < allBitstrings.size(); ++i) {
		result += targetFunction[i].value * parityFunction(index, allBitstrings[i]);
	}
	return result * 1.0 / allBitstrings.size();
}

int parityFunction(string index, string example) {
/* Returns -1 if odd number of 1's in example indexed by index; 1 if even number of 1's */
	int numOnes = 0;
	for (int i = 0; i < min(index.size(), example.size()); ++i) {
		if (index.at(i) == '1' && example.at(i) == '1') {
			++numOnes;
		}
	}
	return (numOnes % 2 == 1 ? -1 : 1);
}

double error(vector<funStruct>& targetFunction, vector<funStruct>& hypothesis) {
/* Probability that our hypothesis does not match the target function */
	int numExamples = targetFunction.size(); 
	int numHalfErrors = 0; // Keeping track of twice the number of errors

	// Probability = number of errors / number of total examples
	for (int i = 0; i < numExamples; ++i) {
		int output = hypothesisOutput(hypothesis, targetFunction[i].index);
		if (targetFunction[i].value * output < 0) { // If our hypothesis has a different sign than the target value
			numHalfErrors += 2;
		}
		else if (targetFunction[i].value * output == 0) {
			numHalfErrors += 1;
		}
	}
	return numHalfErrors * 0.5 / numExamples;
}

int hypothesisOutput(vector<funStruct>& hypothesis, string example) {
/* Returns F(x) for an example x */
	int sum = 0;
	for (int i = 0; i < hypothesis.size(); ++i) {
		sum += hypothesis[i].value * parityFunction(hypothesis[i].index, example);
	}
	return sum;
}

void printInfo(int loopIter, double currentError, double minError, vector<funStruct>& currentHypothesis, bool isNew) {
	cout << "Loop iteration: " << loopIter << endl;
	cout << "Current error: " << currentError << endl;
	cout << "Minimum error: " << minError << endl;
	cout << "Current hypothesis: {";
	int i;
	for (i = 0; i < currentHypothesis.size()-1; ++i) {
		cout << currentHypothesis[i].index << ": " << currentHypothesis[i].value << ", ";
	}
	cout << currentHypothesis[i].index << ": " << currentHypothesis[i].value << "}" << endl;
	if (isNew) {
		cout << "New coefficient!" << endl;
	}
	cout << endl;
}

void findUpNeighbors(vector<string>& bitstrings, vector<string>& neighbors) {
/* Similar to find_down_neighbors, with the folloing additions:
       1. This lakes a list of bitstrings as input
       2. Each of the results' down neighbors must be in this list.
       3. The result itself must not be in this list. */
	for (int i = 0; i < bitstrings.size(); ++i) {
		for (int j = 0; j < bitstrings[i].size(); ++j) { // Find all of the up neighbors of this string
			if (bitstrings[i].at(j) == '0') {
				string neighborCandidate = upgraded(bitstrings[i], j);

				// Check if candidate is already in the current set or the neighbor set
				bool alreadyExists = false;
				for (int k = 0; k < bitstrings.size(); ++k) {
					if (neighborCandidate == bitstrings[k]) {
						alreadyExists = true;
					}
				}
				for (int k = 0; k < neighbors.size(); ++k) {
					if (neighborCandidate == neighbors[k]) {
						alreadyExists = true;
					}
				}

				// If not, check if candidate is a valid neighbor, and add it to our neighbor list if true
				if (!alreadyExists) {
					bool isNeighbor = true;
					vector<string> downNeighbors;
					findDownNeighbors(neighborCandidate, downNeighbors);
					for (int k = 0; k < downNeighbors.size(); ++k) {
						bool downNeighborExists = false;
						for (int m = 0; m < bitstrings.size(); ++m) {
							if (downNeighbors[k] == bitstrings[m]) {
								downNeighborExists = true;
							}
						}
						if (!downNeighborExists) {
							isNeighbor = false;
						}
					}
					if (isNeighbor) {
						neighbors.push_back(neighborCandidate);
					}
				}
			}
		}
	}
}

void findDownNeighbors(string bitstring, vector<string>& neighbors) {
/* Finds all neighbors of the given bitstring with one fewer 1. */
	for (int i = 0; i < bitstring.size(); ++i) {
		if (bitstring.at(i) == '1') {
			neighbors.push_back(downgraded(bitstring, i));
		}
	}
}

int parity(string bitstring) {
// Returns 1 if odd number of 1's; 0 otherwise
	int numOnes = 0;
	for (int i = 0; i < bitstring.size(); ++i) {
		if (bitstring.at(i) == '1') {
			++numOnes;
		}
	}
	return numOnes % 2;
}

int valueAt(vector<funStruct>& s, string idx) {
	for (int i = 0; i < s.size(); ++i) {
		if (s[i].index == idx) {
			return s[i].value;
		}
	}
	return INT_MAX;
}

double valueAt(vector<doubleFunStruct>& s, string idx) {
	for (int i = 0; i < s.size(); ++i) {
		if (s[i].index == idx) {
			return s[i].value;
		}
	}
	return (double) INT_MAX;
}

bool valueAt(vector<string>& s, string idx) {
	for (int i = 0; i < s.size(); ++i) {
		if (s[i] == idx) {
			return true;
		}
	}
	return false;
}

double difference(vector<funStruct>& targetFunction, vector<funStruct>& hypothesis, vector<doubleFunStruct>& fourierCoef, vector<string>& allBitstrings, int diffMethod) {
// Returns the difference (or some other relevant metric) between the expected value of our hypothesis and the expected value of the "ideal" hypothesis.
	double ev = expectedValue(hypothesis, allBitstrings);
	double idealEV = 0.0;
	for (int i = 0; i < hypothesis.size(); ++i) {
		idealEV += valueAt(fourierCoef, hypothesis[i].index) * hypothesis[i].value;
	}
	if (diffMethod == DIFFERENCE) {
		return ev - idealEV;
	}
	if (diffMethod == RATIO) {
		return (ev - idealEV) / ev;
	}
	return 0.0;
}

double expectedValue(vector<funStruct>& hypothesis, vector<string>& allBitstrings) {
// Returns the expected magnitude of F over all possible x
	int sumMagnitude = 0;
	for (int i = 0; i < allBitstrings.size(); ++i) {
		sumMagnitude += abs(hypothesisOutput(hypothesis, allBitstrings[i]));
	}
	return sumMagnitude * 1.0 / allBitstrings.size();
}