#ifndef SIMSTRUCT_H
#define SIMSTRUCT_H

#define INT_MAX 32767

#include <vector>
#include <algorithm>
#include <queue>
#include <functional>
#include <iostream>
#include <thread>
#include <string>
#include <sstream>
#include <fstream>
#include "Graph.h"
#include "Random.h"
#include "alias.h"
#include <unordered_map>
#include <unordered_set>
#include "Timer.h"
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

/*
Make a directory
s: path of the directory
*/
int mkpath(string s, mode_t mode = 0755) {
    size_t pre = 0, pos;
    string dir;
    int mdret;
    if (s[s.size() - 1] != '/') {
        s += '/';
    }
    while ((pos = s.find_first_of('/', pre)) != string::npos) {
        dir = s.substr(0, pos++);
        pre = pos;
        if (dir.size() == 0) continue;
        if ((mdret = ::mkdir(dir.c_str(), mode)) && errno != EEXIST) {
            return mdret;
        }
    }
    return mdret;
}

// sort int in increasing order
bool comp1(const int &a, const int &b) {
	return a < b;
}

// for min heap inplemented by priority_queue
struct comp2 {
    bool operator() (const int& a, const int& b)
    {
        return a > b; // min heap
    }
};

// sort pair<int, double> in decreasing order
bool maxScoreCmp(const pair<int, double>& a, const pair<int, double>& b){
    return a.second > b.second;
}

// sort pair<int, double> in increasing order
bool minScoreCmp(const pair<int, double>& a, const pair<int, double>& b) {
    return a.second < b.second;
}

class SimStruct{
  public:
    Graph g; // Class Graph
    Random R; // Class Random
    int vert; // the number of vertices
    string filelabel; 
    double sqrtC;                                                      
    double C_value;
    double epsilon;
    double rmax;
    double nr;
    double back_nr;
    double forward_nr;
    double forward_rmax; // forward_push part's parameter 
    double backward_rmax; // backward_search's parameter

    double avg_pref_time, avg_mab_time;
    double avg_chk_nbr_time;

    double *H[2];//hash map,the same as U,C,UC
    int* U[2];
    int* C[1];
    int* UC[1];
    int* Parent; // 
    //double *reserve;
    double *residue;
    double *newReserve;
    double *newResidue;
    bool* isInArray;

    // for estimation of D
    unordered_map<int, long> hitMap;
    unordered_map<int, long> totalMap;

	/* for prefiltering */
	vector<pair<int, double> > estMean;
	vector<pair<int, double> > prefCB;
	double *answer;
	double *ans2;
	int *ansIdx;
	int curIdx;
    //double* tempAnswer;
    //double* tempDeterAnswer;
    //int* tempAnsIdx;
    //int* tempDeterAnsIdx;
    //int tempCurIdx, tempDeterCurIdx;
    unordered_map<int, double> largePartAnswer, smallPartAnswer;

	double eps_p;
	int totalSample; 
	// time
	double t_pref, t_single;

    //unordered_map<int, double> dValue;

    unordered_map<int, vector<int> > sameInNbrMap;

    /* data structures for MAB */
    //            query/cand   ->    level              -> w  -> \pi_l(s,w), r_l(s,w)
    unordered_map<int, unordered_map<int, unordered_map<int, double> > > adaptReserveMap, adaptResidueMap;
    unordered_map<int, double> part1Map, part2Map, part3Map;

    vector<int> candList; // ansIdx + curIdx
    // unordered_map<int, double> mean; // answer
    // unordered_map<int, double> sum_Xi2; // ans2
    unordered_map<int, long> cntMABSample;
    long totalMABSample;
    unordered_map<int, int> budgetMap;
    unordered_map<int, double> rsumMap;
    // unordered_map<int, Alias> aliasMap;
    vector<Alias> aliasMap;
    unordered_map<int, int> nodeToAliasPos;
    double adaptPushTime, adaptDeterQueryTime, adaptRandQueryTime;

    SimStruct(string name, string file_lable, double eps, double c) {
        filelabel = file_lable;
        g.inputGraph(name);
        R = Random();
        vert = g.n;
        C_value = c;
        sqrtC = sqrt(C_value);
        epsilon = eps;
		rmax = (1-sqrtC)*(1-sqrtC)/25*epsilon;
		back_nr = (int)(20*log(vert)/epsilon/epsilon);
		nr = (int)(0.1*back_nr);

        forward_rmax = rmax * 2;
        backward_rmax = rmax;
        avg_pref_time = 0;
        avg_mab_time = 0;
        avg_chk_nbr_time = 0;
        H[0] = new double[vert];
        H[1] = new double[vert];
        U[0] = new int[vert];
        U[1] = new int[vert];
        C[0] = new int[vert];
        UC[0] = new int[vert];
        Parent = new int[vert]; //
        //reserve = new double[vert];
        residue = new double[vert];
        newReserve = new double[vert];
        newResidue = new double[vert];
        isInArray = new bool[vert];
        for(int i = 0; i < vert; i++){
            isInArray[i] = false;
            H[0][i] = 0;
            H[1][i] = 0;
            U[0][i] = 0;
            U[1][i] = 0;
            C[0][i] = 0;
            UC[0][i] = -1;
            Parent[i] = -1;
        }
        srand(unsigned(time(0)));
        
		// initialization
		answer = new double[vert];
		ans2 = new double[vert];
        //tempAnswer = new double[vert];
        //tempDeterAnswer = new double[vert];
        ansIdx = new int[vert];
        //tempAnsIdx = new int[vert];
        //tempDeterAnsIdx = new int[vert];
        for(int i = 0; i < vert; i++){
            answer[i] = 0;
			ans2[i] = 0;
            isInArray[i] = false;

            //tempAnswer[i] = 0;
            //tempDeterAnswer[i] = 0;
        }		
		curIdx = 0;
        //tempCurIdx = 0;
        //tempDeterCurIdx = 0;
        cout << "====init done!====" << endl;
    }
    ~SimStruct() {
        delete[] H[0];
        delete[] H[1];
        delete[] U[0];
        delete[] U[1];
        delete[] C[0];
        delete[] UC[0];
        delete[] Parent; //
        //delete[] reserve;
        delete[] residue;
        delete[] newReserve;
        delete[] newResidue;
        delete[] isInArray;
		// for prefiltering
		delete[] answer;
		delete[] ans2;
		delete[] ansIdx;
        //delete[] tempAnswer;
        //delete[] tempAnsIdx;
        //delete[] tempDeterAnswer;
        //delete[] tempDeterAnsIdx;
    }

    //Estimate pi(u,w) using method FORA    
    unordered_map<int, vector<pair<int, double> > > foraMap(int u) {
        unordered_map<int, unordered_map<int, double> > answer;
        unordered_map<int, vector<pair<int, double> > > realAnswer;

        residue[u] = 1;
        vector<int> candidate_set;
        candidate_set.push_back(u);
        isInArray[u] = true;
        vector<int> new_candidate_set;
        int tempLevel = 0;

        Timer fora_timer;
        fora_timer.start();
        vector<pair<pair<int, int>, double> > aliasP;
        double rsum = 0;
        while (candidate_set.size() > 0) {
            for (int j = 0; j < candidate_set.size(); j++) {
                isInArray[candidate_set[j]] = false;
            }
            int residue_count = 0;
            for (int j = 0; j < candidate_set.size(); j++) {
                int tempNode = candidate_set[j];
                double tempR = residue[tempNode];
                newReserve[tempNode] += (1 - sqrtC) * tempR;
                int inSize = g.getInSize(tempNode);
                for (int k = 0; k < inSize; k++) {
                    int newNode = g.getInVert(tempNode, k);
                    newResidue[newNode] += residue[tempNode] * sqrtC / (double)inSize;
                    if (U[1][newNode] == 0) {
                        U[0][residue_count++] = newNode;
                        U[1][newNode] = 1;
                    }
                    if (!isInArray[newNode] && newResidue[newNode] > forward_rmax) {
                        isInArray[newNode] = true;
                        new_candidate_set.push_back(newNode);
                    }
                }
                residue[tempNode] = 0;
            }
            for (int j = 0; j < candidate_set.size(); j++) {
                int tempNode = candidate_set[j];
                if (newReserve[tempNode] > 0) {
                    answer[tempLevel][tempNode] = newReserve[tempNode];
                }
                newReserve[tempNode] = 0;
            }
            for (int j = 0; j < residue_count; j++) {
                if (newResidue[U[0][j]] <= forward_rmax) {
                    rsum += newResidue[U[0][j]];
                    aliasP.push_back(pair<pair<int, int>, double>(pair<int, int>(tempLevel + 1, U[0][j]), newResidue[U[0][j]]));
                }
                else {
                    residue[U[0][j]] = newResidue[U[0][j]];
                }
                newResidue[U[0][j]] = 0;
                U[1][U[0][j]] = 0;
            }
            candidate_set = new_candidate_set;
            new_candidate_set.clear();
            tempLevel++;
        }
        Alias alias = Alias(aliasP);
        fora_timer.end();
        cout << "forward push : " << fora_timer.timeCost << "s" << endl;
        if (rsum > 0) {
            int fora_walk_num = fora_timer.timeCost * 400000;	// 
            cout << "fora_walk_num = " << fora_walk_num << endl;
            fora_timer.start();
            double increment = rsum * (1 - sqrtC) / (double)fora_walk_num;
            // double sqrtSqrtC = sqrt(sqrtC);
            for (int j = 0; j < fora_walk_num; j++) {
                pair<int, int> tempPair = alias.generateRandom(R);
                int tempLevel = tempPair.first;
                int tempNode = tempPair.second;
                int tempCount = 0; // do not use 1 - \sqrt(1 - alpha) walk

                if (answer.find(tempLevel) != answer.end() && answer[tempLevel].find(tempNode) != answer[tempLevel].end())
                    answer[tempLevel][tempNode] += increment;
                else
                    answer[tempLevel][tempNode] = increment;
                while (R.drand() < sqrtC) { // 
                    int length = g.getInSize(tempNode);
                    if (length > 0) {
                        int r = R.generateRandom() % length;
                        tempNode = g.getInVert(tempNode, r);
                        tempCount++;
                    }
                    else {
                        break;
                    }
                    if (answer.find(tempLevel + tempCount) != answer.end() && answer[tempLevel + tempCount].find(tempNode) != answer[tempLevel + tempCount].end())
                        answer[tempLevel + tempCount][tempNode] += increment;
                    else
                        answer[tempLevel + tempCount][tempNode] = increment;
                }
            }
            fora_timer.end();
            cout << "MC : " << fora_timer.timeCost << "s" << endl;
        }
        int numFORAItems = 0, numFORAItemsPruned = 0;
        for (auto& c : answer) {
            for (auto& d : c.second) {
                realAnswer[c.first].push_back(pair<int, double>(d.first, d.second));
                numFORAItems++;
                if (d.second > forward_rmax)
                    numFORAItemsPruned++;
            }
        }
        cout << "back_nr = " << back_nr << " , numFORAItems = " << numFORAItems << " , numFORAItemsPruned = " << numFORAItemsPruned << endl;
        return realAnswer;
    }

    void foraMap(int u, unordered_map<int, unordered_map<int, double> > &foraRealAnswer, unordered_map<int, unordered_map<int, double> >& foraPrunedAnswer) {
        //unordered_map<int, unordered_map<int, double> > answer;
        //unordered_map<int, vector<pair<int, double> > > realAnswer;
        foraRealAnswer.clear();
        foraPrunedAnswer.clear();

        residue[u] = 1;
        vector<int> candidate_set;
        candidate_set.push_back(u);
        isInArray[u] = true;
        vector<int> new_candidate_set;
        int tempLevel = 0;

        Timer fora_timer;
        fora_timer.start();
        vector<pair<pair<int, int>, double> > aliasP;
        double rsum = 0;
        while (candidate_set.size() > 0) {
            for (int j = 0; j < candidate_set.size(); j++) {
                isInArray[candidate_set[j]] = false;
            }
            int residue_count = 0;
            for (int j = 0; j < candidate_set.size(); j++) {
                int tempNode = candidate_set[j];
                double tempR = residue[tempNode];
                newReserve[tempNode] += (1 - sqrtC) * tempR;
                int inSize = g.getInSize(tempNode);
                for (int k = 0; k < inSize; k++) {
                    int newNode = g.getInVert(tempNode, k);
                    newResidue[newNode] += residue[tempNode] * sqrtC / (double)inSize;
                    if (U[1][newNode] == 0) {
                        U[0][residue_count++] = newNode;
                        U[1][newNode] = 1;
                    }
                    if (!isInArray[newNode] && newResidue[newNode] > forward_rmax) {
                        isInArray[newNode] = true;
                        new_candidate_set.push_back(newNode);
                    }
                }
                residue[tempNode] = 0;
            }
            for (int j = 0; j < candidate_set.size(); j++) {
                int tempNode = candidate_set[j];
                if (newReserve[tempNode] > 0) {
                    foraRealAnswer[tempLevel][tempNode] = newReserve[tempNode];
                }
                newReserve[tempNode] = 0;
            }
            for (int j = 0; j < residue_count; j++) {
                if (newResidue[U[0][j]] <= forward_rmax) {
                    rsum += newResidue[U[0][j]];
                    aliasP.push_back(pair<pair<int, int>, double>(pair<int, int>(tempLevel + 1, U[0][j]), newResidue[U[0][j]]));
                }
                else {
                    residue[U[0][j]] = newResidue[U[0][j]];
                }
                newResidue[U[0][j]] = 0;
                U[1][U[0][j]] = 0;
            }
            candidate_set = new_candidate_set;
            new_candidate_set.clear();
            tempLevel++;
        }
        Alias alias = Alias(aliasP);
        fora_timer.end();
        cout << "forward push : " << fora_timer.timeCost << "s" << endl;
        if (rsum > 0) {
            int fora_walk_num = fora_timer.timeCost * 400000;	// original: 400000
            cout << "fora_walk_num = " << fora_walk_num << endl;
            fora_timer.start();
            double increment = rsum * (1 - sqrtC) / (double)fora_walk_num;
            // double sqrtSqrtC = sqrt(sqrtC);
            for (int j = 0; j < fora_walk_num; j++) {
                pair<int, int> tempPair = alias.generateRandom(R);
                int tempLevel = tempPair.first;
                int tempNode = tempPair.second;
                int tempCount = 0; // do not use 1 - \sqrt(1 - alpha) walk

                if (foraRealAnswer.find(tempLevel) != foraRealAnswer.end() && foraRealAnswer[tempLevel].find(tempNode) != foraRealAnswer[tempLevel].end())
                    foraRealAnswer[tempLevel][tempNode] += increment;
                else
                    foraRealAnswer[tempLevel][tempNode] = increment;
                while (R.drand() < sqrtC) { // 
                    int length = g.getInSize(tempNode);
                    if (length > 0) {
                        int r = R.generateRandom() % length;
                        tempNode = g.getInVert(tempNode, r);
                        tempCount++;
                    }
                    else {
                        break;
                    }
                    if (foraRealAnswer.find(tempLevel + tempCount) != foraRealAnswer.end() && foraRealAnswer[tempLevel + tempCount].find(tempNode) != foraRealAnswer[tempLevel + tempCount].end())
                        foraRealAnswer[tempLevel + tempCount][tempNode] += increment;
                    else
                        foraRealAnswer[tempLevel + tempCount][tempNode] = increment;
                }
            }
            fora_timer.end();
            cout << "MC : " << fora_timer.timeCost << "s" << endl;
        }
        int numFORAItems = 0, numFORAItemsPruned = 0;
        for (auto& c : foraRealAnswer) {
            for (auto& d : c.second) {
                //realAnswer[c.first].push_back(pair<int, double>(d.first, d.second));
                numFORAItems++;
                if (d.second > forward_rmax) { // forward_rmax
                    numFORAItemsPruned++;
                    foraPrunedAnswer[c.first][d.first] = d.second;
                }
            }
        }
        cout << "back_nr = " << back_nr << " , numFORAItems = " << numFORAItems << " , numFORAItemsPruned = " << numFORAItemsPruned << endl;
    }

    // estimate d(w)
    int sampleD(int nodeId, int walk_num){
        int tempCount = 0;
        for(int i = 0; i < walk_num; i++){
            int u_newNode = nodeId, v_newNode = nodeId, u_nextNode, v_nextNode;
            double meet_level = 0;
            while(R.drand() < C_value){
                int length = g.getInSize(u_newNode);
                if(length == 0)
                    break;
                int r = R.generateRandom() % length;
                u_nextNode = g.getInVert(u_newNode, r);
                length = g.getInSize(v_newNode);
                if(length == 0)
                    break;
				r = R.generateRandom() % length;
                v_nextNode = g.getInVert(v_newNode, r);
                meet_level++;
                if(u_nextNode == v_nextNode){
                    if(meet_level > 1)
                        tempCount += 1;
                    break;
                }
                u_newNode = u_nextNode;
                v_newNode = v_nextNode;
            }
        }
        return tempCount;
    }

    /*
    int sampleD2(int nodeId, int walk_num) {
        int tempCount = 0;
        for (int i = 0; i < walk_num; i++) {
            int u_newNode = nodeId, v_newNode = nodeId, u_nextNode, v_nextNode;
            while (R.drand() < C_value) {
                int length = g.getInSize(u_newNode);
                if (length == 0)
                    break;
                int r = R.generateRandom() % length;
                u_nextNode = g.getInVert(u_newNode, r);
                length = g.getInSize(v_newNode);
                if (length == 0)
                    break;
                r = R.generateRandom() % length;
                v_nextNode = g.getInVert(v_newNode, r);
                if (u_nextNode == v_nextNode) {
                    tempCount++;
                    break;
                }
                u_newNode = u_nextNode;
                v_newNode = v_nextNode;
            }
        }
        return tempCount;
    }
    */
    
    // Calculating pi(v,w)
    void randomProbe(int w, int targetLevel, double dw, double* answer) {
        //cout << "w = " << w << " , targetLevel = " << targetLevel << endl;
        double sqrtSqrtC = sqrt(sqrtC);
        double increment = pow(sqrtSqrtC, targetLevel) * (1 - sqrtC);
        unordered_map<int, double> walkMap;

        int ind = 0;
        H[ind][w] = 1;
        int Ucount = 1;
        int Ucount1 = 0;
        int UCcount = 0;
        U[0][0] = w;
        for (int i = 0; i < targetLevel; i++) {
            //cout << "level (i) = " << i << endl;
            for (int j = 0; j < Ucount; j++) {
                int tempNode = U[ind][j];
                int outCount = g.getOutSize(tempNode);
                for (int k = 0; k < outCount; k++) {
                    int newNode = g.getOutVert(tempNode, k);
                    if (C[0][newNode] == 0) {
                        C[0][newNode] = 1;
                        UC[0][UCcount] = newNode;
                        UCcount++;
                    }
                    else {
                        C[0][newNode]++;
                    }
                }
            }

            for (int j = 0; j < UCcount; j++) {
                int tempNode = UC[0][j];
                if (R.drand() < C[0][tempNode] * sqrtSqrtC / (double)g.getInSize(tempNode)) {
                    H[1 - ind][tempNode] = 1;
                    U[1 - ind][Ucount1] = tempNode;
                    Ucount1++;
                    //cout << tempNode << endl;
                }
                C[0][UC[0][j]] = 0;
                UC[0][j] = -1;
            }
            for (int j = 0; j < Ucount; j++) {
                H[ind][U[ind][j]] = 0;
                U[ind][j] = -1;
            }
            Ucount = Ucount1;
            Ucount1 = 0;
            UCcount = 0;
            ind = 1 - ind;
            if (Ucount == 0)
                break;
        }
        
        for (int i = 0; i < Ucount; i++) {
            int tempNode = U[ind][i];
            if (walkMap.find(tempNode) == walkMap.end()) {
                walkMap[tempNode] = H[ind][tempNode] * increment;
            }
            else {
                walkMap[tempNode] += H[ind][tempNode] * increment;
            }
            U[ind][i] = 0;
            H[ind][tempNode] = 0;
        }
        Ucount = 0;

        for (auto& key : walkMap) {
            int tempNode = key.first;
            double tempSim = key.second;
            if (answer[tempNode] == 0) {
                ansIdx[curIdx] = tempNode;
                curIdx++;
            }
            double score = dw * tempSim / (1 - sqrtC) / (1 - sqrtC);
            answer[tempNode] += score;
            ans2[tempNode] += score * score;
        }
    }

    void randomProbe(int w, int targetLevel, double dw, unordered_map<int, double>& partialAnswer) {
        //cout << "w = " << w << " , targetLevel = " << targetLevel << endl;
        partialAnswer.clear();
        double sqrtSqrtC = sqrt(sqrtC);
        double increment = pow(sqrtSqrtC, targetLevel) * (1 - sqrtC);
        unordered_map<int, double> walkMap;

        int ind = 0;
        H[ind][w] = 1;
        int Ucount = 1;
        int Ucount1 = 0;
        int UCcount = 0;
        U[0][0] = w;
        for (int i = 0; i < targetLevel; i++) {
            //cout << "level (i) = " << i << endl;
            for (int j = 0; j < Ucount; j++) {
                int tempNode = U[ind][j];
                int outCount = g.getOutSize(tempNode);
                for (int k = 0; k < outCount; k++) {
                    int newNode = g.getOutVert(tempNode, k);
                    if (C[0][newNode] == 0) {
                        C[0][newNode] = 1;
                        UC[0][UCcount] = newNode;
                        UCcount++;
                    }
                    else {
                        C[0][newNode]++;
                    }
                }
            }

            for (int j = 0; j < UCcount; j++) {
                int tempNode = UC[0][j];
                if (R.drand() < C[0][tempNode] * sqrtSqrtC / (double)g.getInSize(tempNode)) {
                    H[1 - ind][tempNode] = 1;
                    U[1 - ind][Ucount1] = tempNode;
                    Ucount1++;
                    //cout << tempNode << endl;
                }
                C[0][UC[0][j]] = 0;
                UC[0][j] = -1;
            }
            for (int j = 0; j < Ucount; j++) {
                H[ind][U[ind][j]] = 0;
                U[ind][j] = -1;
            }
            Ucount = Ucount1;
            Ucount1 = 0;
            UCcount = 0;
            ind = 1 - ind;
            if (Ucount == 0)
                break;
        }

        for (int i = 0; i < Ucount; i++) {
            int tempNode = U[ind][i];
            if (walkMap.find(tempNode) == walkMap.end()) {
                walkMap[tempNode] = H[ind][tempNode] * increment;
            }
            else {
                walkMap[tempNode] += H[ind][tempNode] * increment;
            }
            U[ind][i] = 0;
            H[ind][tempNode] = 0;
        }
        Ucount = 0;

        for (auto& key : walkMap) {
            int tempNode = key.first;
            double tempSim = key.second;
            double score = dw * tempSim / (1 - sqrtC) / (1 - sqrtC);
            partialAnswer[tempNode] = score;
        }
    }

    //Calculating pi(v,w)
    void randomProbeOneHopDeter(int w, int targetLevel, double dw, double* answer) {
        double sqrtSqrtC = sqrt(sqrtC);
        //double increment = pow(sqrtSqrtC, targetLevel) * (1 - sqrtC);
        unordered_map<int, double> walkMap;
        //cout << "w = " << w << " , targetLevel = " << targetLevel << endl;
        int ind = 0;
        H[ind][w] = (1 - sqrtC);
        int Ucount = 1;
        int Ucount1 = 0;
        int UCcount = 0;
        U[0][0] = w;
        for (int i = 0; i < targetLevel; i++) {
            //cout << "level (i) = " << i << endl;
            if (i < 1) {
                //cout << "deter" << endl;
                for (int j = 0; j < Ucount; j++) {
                    int tempNode = U[ind][j];
                    int outCount = g.getOutSize(tempNode);
                    for (int k = 0; k < outCount; k++) {
                        int newNode = g.getOutVert(tempNode, k);
                        //cout << "newNode = " << newNode << endl;
                        if (H[1 - ind][newNode] == 0) {
                            U[1 - ind][Ucount1] = newNode;
                            Ucount1++;
                        }
                        H[1 - ind][newNode] += sqrtC * H[ind][tempNode] / g.getInSize(newNode);
                    }
                }
                for (int j = 0; j < Ucount; j++) {
                    H[ind][U[ind][j]] = 0;
                    U[ind][j] = -1;
                }
                Ucount = Ucount1;
                Ucount1 = 0;
                UCcount = 0;
                ind = 1 - ind;
                if (Ucount == 0)
                    break;
            }
            else {
                //cout << "rand" << endl;
                for (int j = 0; j < Ucount; j++) {
                    int tempNode = U[ind][j];
                    int outCount = g.getOutSize(tempNode);
                    for (int k = 0; k < outCount; k++) {
                        int newNode = g.getOutVert(tempNode, k);
                        if (C[0][newNode] == 0) {
                            C[0][newNode] = 1;
                            Parent[newNode] = tempNode;
                            UC[0][UCcount] = newNode;
                            UCcount++;
                        }
                        else {
                            C[0][newNode]++;
                            if(R.drand() < 1.0 / C[0][newNode])
                                Parent[newNode] = tempNode;
                        }
                    }
                }

                for (int j = 0; j < UCcount; j++) {
                    int tempNode = UC[0][j];
                    if (R.drand() < C[0][tempNode] * sqrtSqrtC / (double)g.getInSize(tempNode)) {
                        //test
                        if (Parent[tempNode] == -1) {
                            cout << "error: tempNode = " << tempNode << " , Parent[tempNode] = " << Parent[tempNode] << endl;
                            exit(0);
                        } // 
                        //cout << "tempNode = " << tempNode << " , Parent[tempNode] = " << Parent[tempNode] << endl;
                        H[1 - ind][tempNode] = H[ind][Parent[tempNode]] * sqrtSqrtC;
                        U[1 - ind][Ucount1] = tempNode;
                        Ucount1++;
                    }
                    Parent[tempNode] = -1;
                    C[0][UC[0][j]] = 0;
                    UC[0][j] = -1;
                }
                for (int j = 0; j < Ucount; j++) {
                    H[ind][U[ind][j]] = 0;
                    U[ind][j] = -1;
                }
                Ucount = Ucount1;
                Ucount1 = 0;
                UCcount = 0;
                ind = 1 - ind;
                if (Ucount == 0)
                    break;
            } // end else           
        }

        for (int i = 0; i < Ucount; i++) {
            int tempNode = U[ind][i];
            if (walkMap.find(tempNode) == walkMap.end()) {
                walkMap[tempNode] = H[ind][tempNode];
            }
            else {
                walkMap[tempNode] += H[ind][tempNode];
            }
            U[ind][i] = 0;
            H[ind][tempNode] = 0;
        }
        Ucount = 0;

        for (auto& key : walkMap) {
            int tempNode = key.first;
            double tempSim = key.second;
            if (answer[tempNode] == 0) {
                ansIdx[curIdx] = tempNode;
                curIdx++;
            }
            double score = dw * tempSim / (1 - sqrtC) / (1 - sqrtC);
            answer[tempNode] += score;
            ans2[tempNode] += score * score;
        }
    }
   
    //Calculating pi(v,w) using sequential probe method
    void getRandomBackList(int w, int targetLevel, int backWalkNum, double forwardSim, double dw, double* answer) {
        unordered_map<int, double> walkMap;        
        double increment = 1 / (double)backWalkNum * pow(sqrtC, targetLevel) * (1 - sqrtC);
        for (int j = 0; j < backWalkNum; j++) {
            int ind = 0;
            H[ind][w] = 1;
            int Ucount = 1;
            int Ucount1 = 0;
            int UCcount = 0;
            U[0][0] = w;
            for (int i = 0; i < targetLevel; i++) {
                for (int j = 0; j < Ucount; j++) {
                    int tempNode = U[ind][j];
                    int outCount = g.getOutSize(tempNode);
                    double tempMaxInSize = 1 / R.drand();
                    for (int k = 0; k < outCount; k++) {
                        int newNode = g.getOutVert(tempNode, k);
                        if (g.getInSize(newNode) > tempMaxInSize) {
                            break;
                        }
                        if (H[1 - ind][newNode] == 0) {
                            U[1 - ind][Ucount1++] = newNode;
                        }
                        H[1 - ind][newNode] += H[ind][tempNode];
                    }
                }
                for (int j = 0; j < Ucount; j++) {
                    H[ind][U[ind][j]] = 0;
                    U[ind][j] = 0;
                }
                Ucount = Ucount1;
                Ucount1 = 0;
                ind = 1 - ind;
                if (Ucount == 0)
                    break;
            }
            for (int i = 0; i < Ucount; i++) {
                int tempNode = U[ind][i];
                if (walkMap.find(tempNode) == walkMap.end()) {
                    walkMap[tempNode] = H[ind][tempNode] * increment;
                }
                else {
                    walkMap[tempNode] += H[ind][tempNode] * increment;
                }
                U[ind][i] = 0;
                H[ind][tempNode] = 0;
            }
            Ucount = 0;

        }
        for (auto& key : walkMap) {
            int tempNode = key.first;
            double tempSim = key.second;
            if (answer[tempNode] == 0) {
                ansIdx[curIdx] = tempNode;
                curIdx++;
            }
            double score = dw * tempSim / (1 - sqrtC) / (1 - sqrtC);
            answer[tempNode] += score;
            ans2[tempNode] += score * score;
        }
    }

    void getRandomBackList3(int w, int targetLevel, int backWalkNum, double forwardSim, double dw, int back_nr) {
        unordered_map<int, double> walkMap;
        double increment = 1 / (double)backWalkNum * pow(sqrtC, targetLevel) * (1 - sqrtC);
        for (int j = 0; j < backWalkNum; j++) {
            int ind = 0;
            H[ind][w] = 1;
            int Ucount = 1;
            int Ucount1 = 0;
            int UCcount = 0;
            U[0][0] = w;
            for (int i = 0; i < targetLevel; i++) {
                for (int j = 0; j < Ucount; j++) {
                    int tempNode = U[ind][j];
                    int outCount = g.getOutSize(tempNode);
                    double tempMaxInSize = 1 / R.drand();
                    for (int k = 0; k < outCount; k++) {
                        int newNode = g.getOutVert(tempNode, k);
                        if (g.getInSize(newNode) > tempMaxInSize) {
                            break;
                        }
                        if (H[1 - ind][newNode] == 0) {
                            U[1 - ind][Ucount1++] = newNode;
                        }
                        H[1 - ind][newNode] += H[ind][tempNode];
                    }
                }
                for (int j = 0; j < Ucount; j++) {
                    H[ind][U[ind][j]] = 0;
                    U[ind][j] = 0;
                }
                Ucount = Ucount1;
                Ucount1 = 0;
                ind = 1 - ind;
                if (Ucount == 0)
                    break;
            }
            for (int i = 0; i < Ucount; i++) {
                int tempNode = U[ind][i];
                if (walkMap.find(tempNode) == walkMap.end()) {
                    walkMap[tempNode] = H[ind][tempNode] * increment;
                }
                else {
                    walkMap[tempNode] += H[ind][tempNode] * increment;
                }
                U[ind][i] = 0;
                H[ind][tempNode] = 0;
            }
            Ucount = 0;

        }
        for (auto& key : walkMap) {
            int tempNode = key.first;
            double tempSim = key.second;
            if (answer[tempNode] == 0) {
                ansIdx[curIdx] = tempNode;
                curIdx++;
            }
            double score = dw * tempSim * forwardSim / (1 - sqrtC) / (1 - sqrtC);
            answer[tempNode] += back_nr * score;
            ans2[tempNode] += back_nr * score * score;
        }
    }
       
    void getRandomBackList_PRSim(int w, int targetLevel, int backWalkNum, double forwardSim, double dw, unordered_map<int, double>& largePartAnswer) {
        unordered_map<int, double> walkMap;
        double increment = 1 / (double)backWalkNum * pow(sqrtC, targetLevel) * (1 - sqrtC);
        for (int j = 0; j < backWalkNum; j++) {
            int ind = 0;
            H[ind][w] = 1;
            int Ucount = 1;
            int Ucount1 = 0;
            int UCcount = 0;
            U[0][0] = w;
            for (int i = 0; i < targetLevel; i++) {
                for (int j = 0; j < Ucount; j++) {
                    int tempNode = U[ind][j];
                    int outCount = g.getOutSize(tempNode);
                    double tempMaxInSize = 1 / R.drand();
                    for (int k = 0; k < outCount; k++) {
                        int newNode = g.getOutVert(tempNode, k);
                        if (g.getInSize(newNode) > tempMaxInSize) {
                            break;
                        }
                        if (H[1 - ind][newNode] == 0) {
                            U[1 - ind][Ucount1++] = newNode;
                        }
                        H[1 - ind][newNode] += H[ind][tempNode];
                    }
                }
                for (int j = 0; j < Ucount; j++) {
                    H[ind][U[ind][j]] = 0;
                    U[ind][j] = 0;
                }
                Ucount = Ucount1;
                Ucount1 = 0;
                ind = 1 - ind;
                if (Ucount == 0)
                    break;
            }
            for (int i = 0; i < Ucount; i++) {
                int tempNode = U[ind][i];
                if (walkMap.find(tempNode) == walkMap.end()) {
                    walkMap[tempNode] = H[ind][tempNode] * increment;
                }
                else {
                    walkMap[tempNode] += H[ind][tempNode] * increment;
                }
                U[ind][i] = 0;
                H[ind][tempNode] = 0;
            }
            Ucount = 0;
        }
        for (auto& key : walkMap) {
            int tempNode = key.first;
            double tempSim = key.second;
            double score = dw * forwardSim * tempSim / (1 - sqrtC) / (1 - sqrtC);
            if(largePartAnswer.find(tempNode) == largePartAnswer.end() || largePartAnswer[tempNode] == 0)
                largePartAnswer[tempNode] = score;
            else
                largePartAnswer[tempNode] += score;
        }
    }

    void getRandomBackList_PRSim(int w, int targetLevel, int backWalkNum, double forwardSim, double dw, double* answer) {
        unordered_map<int, double> walkMap;
        double increment = 1 / (double)backWalkNum * pow(sqrtC, targetLevel) * (1 - sqrtC);
        for (int j = 0; j < backWalkNum; j++) {
            int ind = 0;
            H[ind][w] = 1;
            int Ucount = 1;
            int Ucount1 = 0;
            int UCcount = 0;
            U[0][0] = w;
            for (int i = 0; i < targetLevel; i++) {
                for (int j = 0; j < Ucount; j++) {
                    int tempNode = U[ind][j];
                    int outCount = g.getOutSize(tempNode);
                    double tempMaxInSize = 1 / R.drand();
                    for (int k = 0; k < outCount; k++) {
                        int newNode = g.getOutVert(tempNode, k);
                        if (g.getInSize(newNode) > tempMaxInSize) {
                            break;
                        }
                        if (H[1 - ind][newNode] == 0) {
                            U[1 - ind][Ucount1++] = newNode;
                        }
                        H[1 - ind][newNode] += H[ind][tempNode];
                    }
                }
                for (int j = 0; j < Ucount; j++) {
                    H[ind][U[ind][j]] = 0;
                    U[ind][j] = 0;
                }
                Ucount = Ucount1;
                Ucount1 = 0;
                ind = 1 - ind;
                if (Ucount == 0)
                    break;
            }
            for (int i = 0; i < Ucount; i++) {
                int tempNode = U[ind][i];
                if (walkMap.find(tempNode) == walkMap.end()) {
                    walkMap[tempNode] = H[ind][tempNode] * increment;
                }
                else {
                    walkMap[tempNode] += H[ind][tempNode] * increment;
                }
                U[ind][i] = 0;
                H[ind][tempNode] = 0;
            }
            Ucount = 0;
        }
        for (auto& key : walkMap) {
            int tempNode = key.first;
            double tempSim = key.second;
            double score = dw * forwardSim * tempSim / (1 - sqrtC) / (1 - sqrtC);
            if (answer[tempNode] == 0) {
                ansIdx[curIdx] = tempNode;
                curIdx++;
            }
            answer[tempNode] += score;
        }
    }

    /* sample-all-arms implementation */
    void singleSourceQuery3(int u) {
        Timer ss3_timer;
        ss3_timer.start();
        cout << "singleSourceQuery3 : u = " << u << endl;
        double single_walk_num = 10 * nr; // for estimating D
        for (int i = 0; i < curIdx; i++) { // we do not reuse the estimation of the previous round(s)
            int nodeId = ansIdx[i];
            answer[nodeId] = 0;
            ans2[nodeId] = 0;
        }
        curIdx = 0;
        totalSample = 0;

        // first time PRSim estimation
        unordered_map<int, vector<pair<int, double> > > foraAnswer;
        foraAnswer = foraMap(u);
        vector<pair<pair<int, int>, double> > aliasP;
        double rsum = 0;
        for (auto& c : foraAnswer) {
            int tempLevel = c.first;
            for (auto& d : c.second) {
                int w = d.first;
                double forwardSim = d.second;
                rsum += forwardSim;
                aliasP.push_back(pair<pair<int, int>, double>(pair<int, int>(tempLevel, w), forwardSim));
            }
        }
        Alias alias = Alias(aliasP);
        cout << "rsum = " << rsum << endl;

        Timer sample_d_timer;
        sample_d_timer.start();
        for (int i = 0; i < (int)single_walk_num; i++) {
            pair<int, int> tempPair = alias.generateRandom(R);
            int tempLevel = tempPair.first;
            int tempNode = tempPair.second;
            int hitCount = sampleD(tempNode, 1);
            if (hitMap.find(tempNode) == hitMap.end())
                hitMap[tempNode] = 0;
            if (totalMap.find(tempNode) == totalMap.end())
                totalMap[tempNode] = 0;
            hitMap[tempNode] += hitCount;
            totalMap[tempNode] += 1;
        }
        sample_d_timer.end();
        t_single = sample_d_timer.timeCost / (int)single_walk_num;
        cout << "sample_d_timer.timeCost = " << sample_d_timer.timeCost << endl;
        cout << "t_single = " << t_single << endl;
        cout << "cal D done!" << endl;

        Timer backwalk_timer;
        backwalk_timer.start();
        for (auto& c : foraAnswer) { 
            int tempLevel = c.first;
            for (auto& d : c.second) {
                int w = d.first;
                double forwardSim = d.second;
                double dw = 0;
                if (g.getInSize(w) == 0)
                    dw = 1;
                else if (hitMap.find(w) != hitMap.end())
                    dw = 1 - C_value / (double)g.getInSize(w) - (double)hitMap[w] / totalMap[w];
                else
                    dw = 1 - C_value / (double)g.getInSize(w);
                getRandomBackList3(w, tempLevel, back_nr * forwardSim / rsum + 1, forwardSim, dw, back_nr);
            }
        }
        backwalk_timer.end();
        cout << "backwalk_timer.timeCost = " << backwalk_timer.timeCost << endl;
        totalSample += back_nr;

        // second time PRSim estimation
        foraAnswer = foraMap(u);
        rsum = 0;
        aliasP.clear();
        for (auto& c : foraAnswer) {
            int tempLevel = c.first;
            for (auto& d : c.second) {
                int w = d.first;
                double forwardSim = d.second;
                rsum += forwardSim;
                aliasP.push_back(pair<pair<int, int>, double>(pair<int, int>(tempLevel, w), forwardSim));
            }
        }
        alias = Alias(aliasP);
        cout << "rsum = " << rsum << endl;

        sample_d_timer.start();
        for (int i = 0; i < (int)single_walk_num; i++) {
            pair<int, int> tempPair = alias.generateRandom(R);
            int tempLevel = tempPair.first;
            int tempNode = tempPair.second;
            int hitCount = sampleD(tempNode, 1);
            if (hitMap.find(tempNode) == hitMap.end())
                hitMap[tempNode] = 0;
            if (totalMap.find(tempNode) == totalMap.end())
                totalMap[tempNode] = 0;
            hitMap[tempNode] += hitCount;
            totalMap[tempNode] += 1;
        }
        sample_d_timer.end();
        t_single = sample_d_timer.timeCost / (int)single_walk_num;
        cout << "sample_d_timer.timeCost = " << sample_d_timer.timeCost << endl;
        cout << "t_single = " << t_single << endl;
        cout << "cal D done!" << endl;

        backwalk_timer.start();
        for (auto& c : foraAnswer) { 
            int tempLevel = c.first;
            for (auto& d : c.second) {
                int w = d.first;
                double forwardSim = d.second;
                double dw = 0;
                if (g.getInSize(w) == 0)
                    dw = 1;
                else if (hitMap.find(w) != hitMap.end())
                    dw = 1 - C_value / (double)g.getInSize(w) - (double)hitMap[w] / totalMap[w];
                else
                    dw = 1 - C_value / (double)g.getInSize(w);
                getRandomBackList3(w, tempLevel, back_nr * forwardSim / rsum + 1, forwardSim, dw, back_nr);
            }
        }
        backwalk_timer.end();
        cout << "backwalk_timer.timeCost = " << backwalk_timer.timeCost << endl;
        totalSample += back_nr; 
        
        cout << "totalSample = " << totalSample << " , back_nr = " << back_nr << endl;
        ss3_timer.end();
        cout << "ss3_timer.timeCost = " << ss3_timer.timeCost << endl;
        t_pref = ss3_timer.timeCost / totalSample;
        answer[u] = 1.0;
    }

    void singleSourceQuery2(int u) {
        Timer ss2_timer;
        ss2_timer.start();
        cout << "singleSourceQuery2 : u = " << u << endl;
        double single_walk_num = 10 * nr;
        for (int i = 0; i < curIdx; i++) {
            int nodeId = ansIdx[i];
            answer[nodeId] = 0;
            ans2[nodeId] = 0;
        }
        curIdx = 0;

        double forward_time = 0, alias_time = 0, backward_time = 0, calD_time = 0, calAns_time = 0;
        unordered_map<int, unordered_map<int, double> > foraRealAnswer;
        unordered_map<int, unordered_map<int, double> > foraPrunedAnswer;
        foraMap(u, foraRealAnswer, foraPrunedAnswer); // depend on forward_rmax
        vector<pair<pair<int, int>, double> > aliasP;
        double rsum = 0;
        for (auto& c : foraRealAnswer) {
            int tempLevel = c.first;
            for (auto& d : c.second) {
                int w = d.first;
                double forwardSim = d.second;
                rsum += forwardSim;
                aliasP.push_back(pair<pair<int, int>, double>(pair<int, int>(tempLevel, w), forwardSim));
            }
        }
        Alias alias = Alias(aliasP);
        cout << "singleSourceQuery2: rsum = " << rsum << endl;
        // estimate d(w)
        //unordered_map<int, int> hitMap;
        //unordered_map<int, int> totalMap;
        Timer sample_d_timer;
        sample_d_timer.start();
        for (int i = 0; i < (int)single_walk_num; i++) {
            pair<int, int> tempPair = alias.generateRandom(R);
            int tempLevel = tempPair.first;
            int tempNode = tempPair.second;
            int hitCount = sampleD(tempNode, 1);
            if (hitMap.find(tempNode) == hitMap.end())
                hitMap[tempNode] = 0;
            if (totalMap.find(tempNode) == totalMap.end())
                totalMap[tempNode] = 0;
            hitMap[tempNode] += hitCount;
            totalMap[tempNode] += 1;
        }
        sample_d_timer.end();
        t_single = sample_d_timer.timeCost / (int)single_walk_num;
        cout << "t_single = " << t_single << endl;
        cout << "cal D done!" << endl;

        Timer step1_timer;
        step1_timer.start();
        // step 1. 
        cout << "Step 1" << endl;
        largePartAnswer.clear();
        for (auto& c : foraPrunedAnswer) { //foraPrunedAnswer
            int tempLevel = c.first;
            for (auto& d : c.second) {
                int w = d.first;
                double forwardSim = d.second;
                double dw = 0;
                if (g.getInSize(w) == 0)
                    dw = 1;
                else if (hitMap.find(w) != hitMap.end())
                    dw = 1 - C_value / (double)g.getInSize(w) - (double)hitMap[w] / totalMap[w];
                else
                    dw = 1 - C_value / (double)g.getInSize(w);
                getRandomBackList_PRSim(w, tempLevel, back_nr * forwardSim / rsum + 1, forwardSim, dw, largePartAnswer);
            }
        }
        step1_timer.end();
        cout << "step1_timer.timeCost = " << step1_timer.timeCost << endl;

        /*
        stringstream ss_test;
        ss_test << "/home/liuyu/SimTab-release-baseline/candidates/dblp-author/k=50/" << "sorted.txt";
        vector<pair<int, double> > test_vec;
        for (auto& y : largePartAnswer) {
            int nodeId = y.first;
            double largePart = y.second;
            if (answer[nodeId] == 0) {
                ansIdx[curIdx] = nodeId;
                curIdx++;
            }
            answer[nodeId] += largePart;
            ans2[nodeId] += largePart * largePart;
            test_vec.push_back(pair<int, double>(nodeId, largePart));
        }
        sort(test_vec.begin(), test_vec.end(), maxScoreCmp);
        ofstream test_of(ss_test.str());
        for (int x = 0; x < test_vec.size(); x++)
            test_of << test_vec[x].first << "\t" << test_vec[x].second << endl;
        test_of.close();
        */
        Timer step2_timer;
        step2_timer.start();
        int count = 0;        
        // step 2
        cout << "Step 2" << endl;        
        for (int i = 0; i < back_nr; i++) { // back_nr
            // first sample a walk from u
            int tempNode = u;
            int tempLevel = 0;
            while (R.drand() < sqrtC) {
                int length = g.getInSize(tempNode);
                if (length == 0)
                    break;
                int r = R.generateRandom() % length;
                tempNode = g.getInVert(tempNode, r);
                tempLevel++;
            }
            double dw = 0;
            if (g.getInSize(tempNode) == 0)
                dw = 1;
            else if (hitMap.find(tempNode) != hitMap.end())
                dw = 1 - C_value / (double)g.getInSize(tempNode) - (double)hitMap[tempNode] / totalMap[tempNode];
            else
                dw = 1 - C_value / (double)g.getInSize(tempNode);
            
            if (foraPrunedAnswer.find(tempLevel) == foraPrunedAnswer.end() || foraPrunedAnswer[tempLevel].find(tempNode) == foraPrunedAnswer[tempLevel].end()) {
                count++;
                //cout << "tempLevel = " << tempLevel << " , tempNode = " << tempNode << " , " << (foraPrunedAnswer.find(tempLevel) == foraPrunedAnswer.end())
                    //<< " , " << (foraPrunedAnswer[tempLevel].find(tempNode) == foraPrunedAnswer[tempLevel].end()) << endl; 
                    //<< (foraPrunedAnswer[tempLevel][tempNode] == 0) << " , " << foraPrunedAnswer[tempLevel][tempNode] << endl;
                randomProbe(tempNode, tempLevel, dw, smallPartAnswer);
                // merge largePartAnswer and smallPartAnswer into answer
                
                for (auto& c : largePartAnswer) {
                    int nodeId = c.first;
                    double largePart = c.second;
                    double smallPart = 0;
                    if (smallPartAnswer.find(nodeId) != smallPartAnswer.end())
                        smallPart = smallPartAnswer[nodeId];
                    if (answer[nodeId] == 0) {
                        ansIdx[curIdx] = nodeId;
                        curIdx++;
                    }
                    answer[nodeId] += largePart + smallPart;
                    ans2[nodeId] += (largePart + smallPart) * (largePart + smallPart);
                }
                
                for(auto& c : smallPartAnswer) {
                    int nodeId = c.first;
                    double smallPart = c.second;
                    if (largePartAnswer.find(nodeId) == largePartAnswer.end() || largePartAnswer[nodeId] == 0) {
                        if (answer[nodeId] == 0) {
                            ansIdx[curIdx] = nodeId;
                            curIdx++;
                        }
                        answer[nodeId] += smallPart;
                        ans2[nodeId] += smallPart * smallPart;
                    }
                }
            }
            else {
                for (auto& c : largePartAnswer) {
                    int nodeId = c.first;
                    double largePart = c.second;
                    if (answer[nodeId] == 0) {
                        ansIdx[curIdx] = nodeId;
                        curIdx++;
                    }
                    answer[nodeId] += largePart;
                    ans2[nodeId] += largePart * largePart;
                }
            }
        }
        step2_timer.end();
        cout << "step2_timer.timeCost = " << step2_timer.timeCost << endl;
        cout << "back_nr = " << back_nr << " , count = " << count << endl;
        
        ss2_timer.end();
        cout << "ss2_timer.timeCost = " << ss2_timer.timeCost << endl;
        totalSample = back_nr;
        t_pref = ss2_timer.timeCost / totalSample;
        answer[u] = 1.0;
    }

    void singleSourceQuery(int u){
        cout << "singleSourceQuery : u = " << u << endl;
		double single_walk_num = 10 * nr;
		for(int i = 0; i < curIdx; i++){
			int nodeId = ansIdx[i];
			answer[nodeId] = 0;
			ans2[nodeId] = 0;
		}
		curIdx = 0;

        double forward_time = 0, alias_time=0, backward_time = 0, calD_time = 0, calAns_time = 0;
        unordered_map<int, vector<pair<int, double> > > forward_map = foraMap(u); // depend on forward_rmax
        vector<pair<pair<int, int>, double> > aliasP;
        double rsum = 0;
        for(auto& c: forward_map){
            int tempLevel = c.first;
            for(auto& d: c.second){
                int w = d.first;
                double forwardSim = d.second;
                rsum += forwardSim;
                aliasP.push_back(pair<pair<int, int>, double>(pair<int, int>(tempLevel,w), forwardSim));
            }
        }
        Alias alias = Alias(aliasP);
        // estimate d(w)
		//unordered_map<int, int> hitMap;
        //unordered_map<int, int> totalMap;
		Timer sample_d_timer;
		sample_d_timer.start();
		for(int i = 0; i < (int)single_walk_num; i++){
            pair<int, int> tempPair = alias.generateRandom(R);
			int tempLevel = tempPair.first;
            int tempNode = tempPair.second;
            int hitCount = sampleD(tempNode,1); 
			if(hitMap.find(tempNode) == hitMap.end())
                hitMap[tempNode] = 0;
            if(totalMap.find(tempNode) == totalMap.end())
                totalMap[tempNode] = 0;
            hitMap[tempNode] += hitCount;
			totalMap[tempNode] += 1;
        }
		sample_d_timer.end();
		t_single = sample_d_timer.timeCost / (int)single_walk_num;
		cout << "t_single = " << t_single << endl;

        // unordered_map<int, double> dValue;
        //dValue.clear();
        //int cnt_w = 0;
		//for(auto& w : hitMap){
            //cnt_w += 1;
		//	if(g.getInSize(w.first) == 0)
        //       dValue[w.first] = 1;
        //    else 
        //       dValue[w.first] = 1 - C_value / (double) g.getInSize(w.first) - w.second / (double)totalMap[w.first];
        //}
        cout << "cal D done!" << endl;
        
		for(int i = 0; i < back_nr; i++){
			// first sample a walk from u
			int tempNode = u;
			int tempLevel = 0;
			while(R.drand() < sqrtC){
				int length = g.getInSize(tempNode);
                if(length == 0)
                    break;
                int r = R.generateRandom() % length;
                tempNode = g.getInVert(tempNode, r);
				tempLevel ++;
			}
//			pair<int, int> tempPair = alias.generateRandom(R);
//			int tempLevel = tempPair.first;
//          int tempNode = tempPair.second;
			double dw = 0;
            if(g.getInSize(tempNode) == 0)
                dw = 1;
            else if(hitMap.find(tempNode) != hitMap.end())
                dw = 1 - C_value / (double)g.getInSize(tempNode) - (double)hitMap[tempNode] / totalMap[tempNode];
            else
                dw = 1 - C_value / (double) g.getInSize(tempNode);
			
			Timer pref_timer;
			pref_timer.start();
			//cout << "tempNode = " << tempNode << " , dw = " << dw << endl;
            randomProbe(tempNode, tempLevel, dw, answer);
            //randomProbeOneHopDeter(tempNode, tempLevel, dw, answer);
            //getRandomBackList(tempNode, tempLevel, 1, 1, dw, answer);
            /*
            randomProbe(tempNode, tempLevel, dw, smallPartAnswer);
            for (auto& c : smallPartAnswer) {
                int nodeId = c.first;
                double smallPart = c.second;
                //if (largePartAnswer.find(nodeId) == largePartAnswer.end() || largePartAnswer[nodeId] == 0) {
                if (answer[nodeId] == 0) {
                    ansIdx[curIdx] = nodeId;
                    curIdx++;
                }
                answer[nodeId] += smallPart;
                ans2[nodeId] += smallPart * smallPart;
                //}
            }
            */
			pref_timer.end();
            t_pref += pref_timer.timeCost;
		}
		totalSample = back_nr;

		t_pref /= totalSample;
		cout << "t_pref = " << t_pref << endl;

     	answer[u] = 1.0;
    }

    double computeCB(int nodeId, double *ans2, double estMean, int totalSample, int vert){
		// invariant:
		// ans2 = \sum_{non-zero answer}{answer * answer}
		// estMean = average over all (including zero) estimations
		double failProb = 0.001; 
		double std_dev = sqrt( max((ans2[nodeId] - totalSample * estMean * estMean), 0.0) / totalSample );
		return std_dev * sqrt(2.0 * log(3.0 / failProb) / totalSample) + 3 * C_value * log(3.0 / failProb) / totalSample;
	}

    /*
    The prefiltering phase
    u: query node
    */
	void prefilter(int u, int k){
		cout << "prefilter : qv = " << u << endl;
        Timer prefilter_timer;
        prefilter_timer.start();
		eps_p = 1.0e-5;
		epsilon = 0.5;
		estMean.clear();
		prefCB.clear();
		
        /* iteratively apply sample-all-arms, i.e., singleSourceQuery(u) */
		do{
			//cout << "epsilon = " << epsilon << endl;
			rmax = (1 - sqrtC) * (1 - sqrtC) / 25 * epsilon;
			back_nr = (int)(20 * log(vert) / epsilon / epsilon);
			nr = (int)(0.1 * back_nr); 
			forward_rmax = rmax * 2;
			backward_rmax = rmax;
			t_pref = 0;
			t_single = 0;
			totalSample = 0;
			// invoke singleQourceQuery, should update t_pref, and t_single, totalSample, 
			clock_t local_ts = clock();
            cout << "forward_rmax = " << forward_rmax << endl;
			//singleSourceQuery(u);  // baseline
            //singleSourceQuery2(u);
            singleSourceQuery3(u);
			clock_t local_te = clock();
			cout << "epsilon = " << epsilon << " , time (not accurate): " << (local_te - local_ts) / (double)CLOCKS_PER_SEC << endl;
			// now answer and ans2 has values
            cout << "t_pref = " << t_pref << endl;
            cout << "t_single = " << t_single << endl;
            if (!(t_pref > 0 && t_single > 0)) {
                cout << "time calculator error" << endl;
                exit(0);
            }
		    // assume the k-th ground truth is not 0
            cout << "curIdx = " << curIdx << endl;
            if (curIdx >= k) {
                vector<pair<int, double> > lbList;
                unordered_map<int, double> ubMap;
                unordered_map<int, double> estMap;
                for (int i = 0; i < curIdx; i++) {
                    int nodeId = ansIdx[i];
                    if (nodeId != u) {
                        double score = answer[nodeId] / totalSample;
                        //if (nodeId == 69161)
                        //    cout << "in prefilter phase, score of 69161 = " << score << endl;
                        double cb = computeCB(nodeId, ans2, score, totalSample, vert);
                        lbList.push_back(pair<int, double>(nodeId, score - cb));
                        ubMap[nodeId] = score + cb;
                        estMap[nodeId] = score;
                    }
                }
                sort(lbList.begin(), lbList.end(), maxScoreCmp); // in descending order, O(nlogn)
                //for (int i = 0; i < lbList.size(); i++)
                //    cout << lbList[i].first << " , " << lbList[i].second << " , " << estMap[lbList[i].first] << endl;
                int kthNode = lbList[k - 1].first;
                double kthLB = lbList[k - 1].second;
                double cb_0 = 3 * C_value * log(1.0 / 0.001) / totalSample;
                cout << "kthNode = " << kthNode << " , kthLB = " << kthLB << " , cb_0 = " << cb_0 << endl;
                if (kthLB + eps_p > cb_0) { // can prune 0
                    vector<int> candidates; // store candidate nodes
                    for (int i = 0; i < k; i++)
                        candidates.push_back(lbList[i].first);
                    for (int i = k; i < lbList.size(); i++) {
                        int cand = lbList[i].first;
                        if (kthLB + eps_p < ubMap[cand]) {
                            candidates.push_back(cand);
                        }
                    }
                    cout << "numCandidate = " << candidates.size() << endl;
                    if (t_pref * max(5 * k, 100) > t_single * candidates.size()) {
                        for (int i = 0; i < candidates.size(); i++) {
                            int cand = candidates[i];
                            estMean.push_back(pair<int, double>(cand, answer[cand] / totalSample));
                            //if (cand == 69161)
                            //    cout << "cand = " << cand << " , " << answer[cand] / totalSample << endl;
                            double ithCB = computeCB(cand, ans2, answer[cand] / totalSample, totalSample, vert);
                            prefCB.push_back(pair<int, double>(cand, ithCB));
                        }
                        break;
                    }
                }
            } // otherwise keep prefiltering
			epsilon /= 2; 
		}while(true);
        prefilter_timer.end();
        avg_pref_time += prefilter_timer.timeCost;
        cout << "u = " << u << " , time(sec) : " << prefilter_timer.timeCost << endl;
	}

	/* output candidate set */
	void outputCandidate(int u, int k){
		stringstream ssout;
		ssout << "candidates/" << filelabel << "/k=" << k << "/" << u << ".txt";
		ofstream fout(ssout.str());
		fout.precision(15);
		for(int i = 0; i < estMean.size(); i++){
            int nodeId = estMean[i].first;
            fout << estMean[i].first << "\t" << estMean[i].second << "\t" << totalSample << "\t" << ans2[nodeId] << "\t" << estMean.size() << endl;
		}
	}

    /*
    sample-one-arm implementation
    */
    void sampleOneArm_Naive(int queryNode, int candNode, int batchSample) {
        int tempCount = 0;
        for (int i = 0; i < batchSample; i++) {
            int u_newNode = queryNode, v_newNode = candNode, u_nextNode, v_nextNode;
            while (R.drand() < C_value) {
                int length = g.getInSize(u_newNode);
                if (length == 0)
                    break;
                int r = R.generateRandom() % length;
                u_nextNode = g.getInVert(u_newNode, r);
                length = g.getInSize(v_newNode);
                if (length == 0)
                    break;
                r = R.generateRandom() % length;
                v_nextNode = g.getInVert(v_newNode, r);

                if (u_nextNode == v_nextNode) {
                    tempCount++;
                    break;
                }
                u_newNode = u_nextNode;
                v_newNode = v_nextNode;
            }
        }
        answer[candNode] += tempCount;
        ans2[candNode] += tempCount;

        cntMABSample[queryNode] += batchSample;
        cntMABSample[candNode] += batchSample;
        totalMABSample += batchSample;
    }
    
    /* The top-k MAB phase */
    vector <pair<int, double> > UGapE(int queryNode, int k) {
        Timer ugape_timer;
        ugape_timer.start();

        double eps_min = 1e-4;
        vector<pair<int, double> > B;
        unordered_map<int, double> UB, LB;
        vector<pair<int, double> > vecUB, vecLB;
        rsumMap.clear();
        aliasMap.clear();
        nodeToAliasPos.clear();
        // clock_t ts_ugape = clock();
        double sortTime = 0, sampleTime = 0;
        adaptPushTime = 0; adaptDeterQueryTime = 0; adaptRandQueryTime = 0;
        Timer sort_timer, sample_timer;
        // mean, cntSample, and sum_Xi2 already have values
        // clock_t tsSample = clock();
        int round = 0;
        int batch = 10;
        int batchSample = 10000; // max(10000, (int)(candList.size()));
        cout << "batchSample = " << batchSample << endl;
        // test
        //for (int i = 0; i < candList.size(); i++) {
        //    int armId = candList[i];
        //    sampleOneArm_Naive(queryNode, armId, 100);
        //} //
        while (true) {
            round++;
            sort_timer.start();
            UB.clear();
            LB.clear();
            vecUB.clear();
            vecLB.clear();
            B.clear();
            //遍历arm统一用candList
            for (int i = 0; i < candList.size(); i++) {
                int armId = candList[i];
                double ugapCB = ugapeCB(candList.size(), answer[armId], ans2[armId], cntMABSample[armId], 1.0e-4, round);
                UB[armId] = answer[armId] / cntMABSample[armId] + ugapCB;
                LB[armId] = answer[armId] / cntMABSample[armId] - ugapCB;
                vecUB.push_back(pair<int, double>(armId, answer[armId] / cntMABSample[armId] + ugapCB));
            }
            sort(vecUB.begin(), vecUB.end(), maxScoreCmp); // descreasing order
            double kthUB = vecUB[k - 1].second; // k-th largest UB
            double k1thUB = vecUB[k].second;    // k+1-th largest UB
            
            for (int i = 0; i < vecUB.size(); i++) {
                int armId = vecUB[i].first;
                if (i < k)     //(i < k)
                    B.push_back(pair<int, double>(armId, k1thUB - LB[armId]));
                else
                    B.push_back(pair<int, double>(armId, kthUB - LB[armId]));
            }
            vecUB.clear();
            sort(B.begin(), B.end(), minScoreCmp); // increasing order

            double largestB;
            for (int i = 0; i < B.size(); i++) {
                int armId = B[i].first;
                if (i < k) {
                    vecLB.push_back(pair<int, double>(armId, LB[armId]));
                    if (i == k-1)
                        largestB = B[i].second;
                }
                else
                    vecUB.push_back(pair<int, double>(armId, UB[armId]));
            }
            sort(vecLB.begin(), vecLB.end(), minScoreCmp);
            sort(vecUB.begin(), vecUB.end(), maxScoreCmp);
            sort_timer.end();
            sortTime += sort_timer.timeCost;
            
            if (largestB <= eps_min) { // B[k - 1].second
                //vector<int> upageTopk;
                vector <pair<int, double> > ugapeAnswer;
                for (int i = 0; i < vecLB.size(); i++) {
                    int armId = vecLB[i].first;
                    ugapeAnswer.push_back(pair<int, double>(armId, answer[armId] / cntMABSample[armId]));
                    cout << "armId = " << armId << " , cntMABSample[armId] = " << cntMABSample[armId] << " , est = " << answer[armId] / cntMABSample[armId] << endl;
                    //if (sameInNbr.find(armId) != sameInNbr.end()) {
                    //    vector<int>& sameInNbrList = sameInNbr[armId];
                    //    for (int j = 0; j < sameInNbrList.size(); j++)
                    //        upageTopk.push_back(sameInNbrList[j]);
                    //}
                }
                sort(ugapeAnswer.begin(), ugapeAnswer.end(), maxScoreCmp);
                //for (int i = 0; i < candList.size(); i++)
                   //delete aliasMap[i];
                
                //clock_t te_ugape = clock();
                //sampleTime = (te_ugape - tsSample) / (double)CLOCKS_PER_SEC;
                //cout << "UGapE: " << (te_ugape - ts_ugape) / (double)CLOCKS_PER_SEC << endl;
                //cout << "initTime = " << initTime << endl;
                //cout << "sampleTime = " << sampleTime << endl;
                //cout << "totalSample = " << totalSample << endl;
                //cout << "sortTime = " << sortTime << endl;
                ugape_timer.end();
                avg_mab_time += ugape_timer.timeCost;
                stringstream ssout;
                ssout << "ugape/" << filelabel << "/k=" << k << "/";
                mkpath(ssout.str());
                ssout << queryNode << ".txt";
                ofstream fout(ssout.str());
                fout.precision(15);
                for (int i = 0; i < ugapeAnswer.size(); i++) {
                    fout << ugapeAnswer[i].first << "\t" << ugapeAnswer[i].second << endl;
                }
                cout << "ugape_timer.timeCost = " << ugape_timer.timeCost << " , sortTime = " << sortTime << " , sampleTime = " << sampleTime << " , adaptPushTime = " << adaptPushTime  << endl;
                cout << "adaptDeterQueryTime = " << adaptDeterQueryTime << " , adaptRandQueryTime = " << adaptRandQueryTime << endl;
                return ugapeAnswer;
            }
            // sample high
            sample_timer.start();
            for (int i = 0; i < batch; i++) {
                int armId = vecLB[i].first;
                sampleOneArm_Naive(queryNode, armId, batchSample);
                //sampleOneArm(queryNode, armId, batchSample);
            }
            // sample low
            for (int i = 0; i < batch; i++) {
                int armId = vecUB[i].first;
                sampleOneArm_Naive(queryNode, armId, batchSample);
                //sampleOneArm(queryNode, armId, batchSample);
            }
            sample_timer.end();
            sampleTime += sample_timer.timeCost;
        }
    }

    /* confidence bound */
    double ugapeCB(int numOfArms, double answer, double sum_x2, int numSample, double failProb, int round) {
        double est_mean;
        if (numSample == 0)
            est_mean = 0;
        else
            est_mean = answer / numSample;
        double estStd = sqrt((max(sum_x2 - numSample * est_mean * est_mean, 0.0)) / (numSample - 1));
        //cout << "estStd = " << estStd << endl;
        double part1 = estStd * sqrt(log(numOfArms * pow(round, 3.0) / failProb) / numSample);
        double part2 = 7.0 / 6 * log(numOfArms * pow(round, 3.0) / failProb) / (numSample - 1);
        //cout << "part1 = " << part1 << endl;
        //cout << "part2 = " << part2 << endl;
        return part1 + part2;
    }

    bool checkSameInNbr(int cnode) {
        bool found = false;
        for (auto& a : sameInNbrMap) {
            int exist_cnode = a.first;
            vector<int>& v = a.second;
            if (hasSameInNbr(exist_cnode, cnode)) {
                v.push_back(cnode);
                found = true;
                //cout << "(" << cnode << " , " << exist_cnode << ")" << endl;
                break;
            }
        }
        if (!found) {
            vector<int> v;
            v.push_back(cnode);
            sameInNbrMap[cnode] = v;
        }
        return found;
    }

    bool hasSameInNbr(int exist_cnode, int cnode) {
        if (g.getInSize(exist_cnode) != g.getInSize(cnode))
            return false;
        unordered_set<int> exist_cnode_nbr;
        for (int i = 0; i < g.getInSize(exist_cnode); i++)
            exist_cnode_nbr.insert(g.getInVert(exist_cnode, i));
        bool isSame = true;
        for (int i = 0; i < g.getInSize(cnode); i++) {
            int tempNode = g.getInVert(cnode, i);
            if (exist_cnode_nbr.find(tempNode) == exist_cnode_nbr.end()) {
                isSame = false;
                break;
            }
        }
        return isSame;
    }

    void readCandidatesFromPrefilter(string candpath, int qnode) {
        candList.clear();
        totalMABSample = 0;
        cntMABSample.clear();
        if (curIdx != 0) {
            for (int i = 0; i < curIdx; i++) {
                int nodeId = ansIdx[i];
                answer[nodeId] = 0;
                ans2[nodeId] = 0;
            }
            curIdx = 0;
        }
        // file format:
        // nodeId estScore totalAllArmsSample ans2 candidateSize
        stringstream ss;
        ss << candpath << "/" << qnode << ".txt";
        ifstream ifcand(ss.str());
        while (ifcand.good()) {
            int cnode;
            double estScore;
            int cntAllArmsSample;
            double ans2val;
            int candsz;
            ifcand >> cnode >> estScore >> cntAllArmsSample >> ans2val >> candsz;
            if (find(candList.begin(), candList.end(), cnode) == candList.end()) {
                candList.push_back(cnode); // equivalent to ansIdx? no must update ansIdx
                
                answer[cnode] = estScore * cntAllArmsSample;
                cntMABSample[cnode] = cntAllArmsSample;
                totalMABSample += cntAllArmsSample;
                ans2[cnode] = ans2val;

                budgetMap[cnode] = 1;
                adaptResidueMap[cnode][0][cnode] = 1;
            }
        }
        ifcand.close();
        cout << "qnode = " << qnode << " , candList.size() = " << candList.size() << endl;
        budgetMap[qnode] = 1;
        adaptResidueMap[qnode][0][qnode] = 1;

        nodeToAliasPos.clear();
        for (int i = 0; i < candList.size(); i++) {
            int nodeId = candList[i];
            ansIdx[curIdx] = nodeId;
            curIdx++;
        }
    }

    void readCandidatesFromPrefilter2(string candpath, int qnode) {
        // candList.clear(); // do not use candList anymore
        totalMABSample = 0;
        cntMABSample.clear();
        sameInNbrMap.clear();
        if (curIdx != 0) {
            cout << "curIdx = " << curIdx << " , clear the results of the prefiltering phase..." << endl;
            for (int i = 0; i < curIdx; i++) {
                int nodeId = ansIdx[i];
                answer[nodeId] = 0;
                ans2[nodeId] = 0;
            }
            curIdx = 0;
        }
        // file format:
        // nodeId estScore totalAllArmsSample ans2 candidateSize
        stringstream ss;
        ss << candpath << "/" << qnode << ".txt";
        ifstream ifcand(ss.str());
        double chk_nbr_time = 0;
        Timer chkNbrTimer;
        double prev_cnode = -1;
        while (ifcand.good()) {
            int cnode;
            double estScore;
            int cntAllArmsSample;
            double ans2val;
            int candsz;
            ifcand >> cnode >> estScore >> cntAllArmsSample >> ans2val >> candsz;
            if (cnode != prev_cnode) {
                chkNbrTimer.start();
                bool found = checkSameInNbr(cnode);
                chkNbrTimer.end();
                chk_nbr_time += chkNbrTimer.timeCost;
                if (!found) {
                    // candList.push_back(cnode); // equivalent to ansIdx? no must update ansIdx
                    answer[cnode] = estScore * cntAllArmsSample;
                    cntMABSample[cnode] = cntAllArmsSample;
                    totalMABSample += cntAllArmsSample;
                    ans2[cnode] = ans2val;
                    ansIdx[curIdx] = cnode;
                    curIdx++;

                    //budgetMap[cnode] = 1;
                    //adaptResidueMap[cnode][0][cnode] = 1;
                }
                prev_cnode = cnode;
            }
            else
                cout << "end of candidates : " << cnode << " , " << estScore << " , " << cntAllArmsSample << " , " << ans2val << " , " << candsz << endl;
        }
        ifcand.close();
        cntMABSample[qnode] = totalSample;
        cout << "qnode = " << qnode << " , candidate size = " << curIdx << " , sameInNbrMap.size() = " << sameInNbrMap.size() << endl;
        cout << "chk_nbr_time = " << chk_nbr_time << endl;
        avg_chk_nbr_time += chk_nbr_time;
        //budgetMap[qnode] = 1;
        //adaptResidueMap[qnode][0][qnode] = 1;

        //nodeToAliasPos.clear();
        //for (int i = 0; i < candList.size(); i++) {
        //    int nodeId = candList[i];
        //    ansIdx[curIdx] = nodeId;
        //    curIdx++;
        //}
    }

    void MC(int queryNode, int k) {
        Timer mcTimer;
        mcTimer.start();
        double eps_mc = 1.0e-4;
        int numSample = log(curIdx) / eps_mc / eps_mc;
        cout << "MC : numSample = " << numSample << endl;
        for (int i = 0; i < curIdx; i++) {
            int cnode = ansIdx[i];
            sampleOneArm_Naive(queryNode, cnode, numSample);
        }
        mcTimer.end();
        vector <pair<int, double> > ugapeAnswer;
        for (auto& item : sameInNbrMap) {
            int armId = item.first;
            vector<int>& sameNbrVec = sameInNbrMap[armId];
            for (int j = 0; j < sameNbrVec.size(); j++) {
                int tempNode = sameNbrVec[j];
                ugapeAnswer.push_back(pair<int, double>(tempNode, answer[armId] / cntMABSample[armId]));
                cout << "armId = " << tempNode << " , cntMABSample[armId] = " << cntMABSample[armId] << " , est = " << answer[armId] / cntMABSample[armId] << endl;
            }
        }
        sort(ugapeAnswer.begin(), ugapeAnswer.end(), maxScoreCmp);
        stringstream ssout;
        ssout << "MC/" << filelabel << "/k=" << k << "/";
        mkpath(ssout.str());
        ssout << queryNode << ".txt";
        ofstream fout(ssout.str());
        fout.precision(15);
        for (int i = 0; i < ugapeAnswer.size(); i++) {
            fout << ugapeAnswer[i].first << "\t" << ugapeAnswer[i].second << endl;
        }
        cout << "mcTimer.timeCost = " << mcTimer.timeCost << " , # total sample = " << numSample * curIdx << endl;
    }

    vector <pair<int, double> > UGapE2(int queryNode, int k) {
        Timer ugape_timer;
        ugape_timer.start();
        double sortTime = 0, sampleTime = 0;
        adaptPushTime = 0; adaptDeterQueryTime = 0; adaptRandQueryTime = 0;
        Timer sort_timer, sample_timer;

        // first should check if there exist exact k candidates
        int numCand = 0;
        for (auto& item : sameInNbrMap)
            numCand += item.second.size();
        cout << "numCand = " << numCand << endl;
        if (numCand <= k) {
            cout << "less than k candidates, skip sample-one-arm phase..." << endl;
            vector <pair<int, double> > ugapeAnswer;
            for (auto& item : sameInNbrMap) {
                int armId = item.first;
                vector<int>& sameNbrVec = sameInNbrMap[armId];
                for (int j = 0; j < sameNbrVec.size(); j++) {
                    int tempNode = sameNbrVec[j];
                    ugapeAnswer.push_back(pair<int, double>(tempNode, answer[armId] / cntMABSample[armId]));
                    cout << "armId = " << tempNode << " , cntMABSample[armId] = " << cntMABSample[armId] << " , est = " << answer[armId] / cntMABSample[armId] << endl;
                }
            }
            sort(ugapeAnswer.begin(), ugapeAnswer.end(), maxScoreCmp);
            ugape_timer.end();
            avg_mab_time += ugape_timer.timeCost;
            stringstream ssout;
            ssout << "ugape/" << filelabel << "/k=" << k << "/";
            mkpath(ssout.str());
            ssout << queryNode << ".txt";
            ofstream fout(ssout.str());
            fout.precision(15);
            for (int i = 0; i < ugapeAnswer.size(); i++) {
                fout << ugapeAnswer[i].first << "\t" << ugapeAnswer[i].second << endl;
            }
            cout << "ugape_timer.timeCost = " << ugape_timer.timeCost << " , sortTime = " << sortTime << " , sampleTime = " << sampleTime << " , adaptPushTime = " << adaptPushTime << endl;
            cout << "adaptDeterQueryTime = " << adaptDeterQueryTime << " , adaptRandQueryTime = " << adaptRandQueryTime << endl;
            return ugapeAnswer;
        }

        // second
        double eps_min = 1e-3;
        vector<pair<int, double> > B;
        unordered_map<int, double> UB, LB;
        vector<pair<int, double> > vecUB, vecLB;
        unordered_map<int, int> rank; 
        //rsumMap.clear();
        //aliasMap.clear();
        //nodeToAliasPos.clear();       
        
        // mean, cntSample, and sum_Xi2 already have values
        int round = 0;
        int batch = 10;
        int batchSample = 10000; // max(10000, (int)(candList.size()));
        cout << "batchSample = " << batchSample << endl;

        while (true) {
            round++;
            sort_timer.start();
            UB.clear();
            LB.clear();
            vecUB.clear();
            vecLB.clear();
            B.clear();
            for (int i = 0; i < curIdx; i++) {
                int armId = ansIdx[i];
                double ugapCB = ugapeCB(curIdx, answer[armId], ans2[armId], cntMABSample[armId], 1.0e-4, round);
                UB[armId] = answer[armId] / cntMABSample[armId] + ugapCB;
                LB[armId] = answer[armId] / cntMABSample[armId] - ugapCB;
                vecUB.push_back(pair<int, double>(armId, answer[armId] / cntMABSample[armId] + ugapCB));
                //cout << "armId = " << armId << " , ugapCB = " << ugapCB << " , UB[armId] = " << UB[armId] << " , LB[armId] = " << LB[armId] << endl;
            }
            sort(vecUB.begin(), vecUB.end(), maxScoreCmp); // descreasing order
            double kthUB; // k-th largest UB
            double k1thUB;    // k+1-th largest UB
            int accum_rank = 1; // ranking starting from 1
            for (int i = 0; i < vecUB.size(); i++) {
                int armId = vecUB[i].first;
                rank[armId] = accum_rank;
                accum_rank += sameInNbrMap[armId].size();
                //cout << "armId = " << armId << " , rank[armId] = " << rank[armId] << endl;
                if (rank[armId] <= k && accum_rank > k) {
                    kthUB = vecUB[i].second;
                    if(i == vecUB.size() - 1)
                        k1thUB = vecUB[i].second; // strange
                    else
                        k1thUB = vecUB[i + 1].second;
                }
            }
            //cout << "kthUB = " << kthUB << endl;
            //cout << "k1thUB = " << k1thUB << endl;

            for (int i = 0; i < vecUB.size(); i++) {
                int armId = vecUB[i].first;
                if (rank[armId] <= k)     //(i < k)
                    B.push_back(pair<int, double>(armId, k1thUB - LB[armId]));
                else
                    B.push_back(pair<int, double>(armId, kthUB - LB[armId]));
            }
            vecUB.clear();
            sort(B.begin(), B.end(), minScoreCmp); // increasing order
            accum_rank = 1; // ranking starting from 1
            rank.clear();
            double kthB;
            for (int i = 0; i < B.size(); i++) {
                int armId = B[i].first;
                rank[armId] = accum_rank;
                //cout << "armId = " << armId << " , rank[armId] = " << rank[armId] << endl;
                accum_rank += sameInNbrMap[armId].size();
                if (rank[armId] <= k && accum_rank > k)
                    kthB = B[i].second;
            }
            //cout << "kthB = " << kthB << endl;
            
            for (int i = 0; i < B.size(); i++) {
                int armId = B[i].first;
                if (rank[armId] <= k)
                    vecLB.push_back(pair<int, double>(armId, LB[armId]));
                else
                    vecUB.push_back(pair<int, double>(armId, UB[armId]));
            }
            sort(vecLB.begin(), vecLB.end(), minScoreCmp);
            sort(vecUB.begin(), vecUB.end(), maxScoreCmp);
            //cout << "vecLB.size() = " << vecLB.size() << " , vecUB.size() = " << vecUB.size() << endl;
            //cout << "vecLB contains:" << endl;
            //for (int i = 0; i < vecLB.size(); i++)
            //    cout << vecLB[i].first << " , " << endl;
            //cout << "vecUB contains:" << endl;
            //for (int i = 0; i < vecUB.size(); i++)
            //    cout << vecUB[i].first << " , " << endl;

            sort_timer.end();
            sortTime += sort_timer.timeCost;

            if (kthB <= eps_min) { 
                vector <pair<int, double> > ugapeAnswer;
                for (int i = 0; i < vecLB.size(); i++) {
                    int armId = vecLB[i].first;
                    vector<int>& sameNbrVec = sameInNbrMap[armId];
                    for (int j = 0; j < sameNbrVec.size(); j++) {
                        int tempNode = sameNbrVec[j];
                        ugapeAnswer.push_back(pair<int, double>(tempNode, answer[armId] / cntMABSample[armId]));
                        cout << "armId = " << tempNode << " , cntMABSample[armId] = " << cntMABSample[armId] << " , est = " << answer[armId] / cntMABSample[armId] << endl;
                    }
                }
                sort(ugapeAnswer.begin(), ugapeAnswer.end(), maxScoreCmp);
                //for (int i = 0; i < candList.size(); i++)
                   //delete aliasMap[i];

                //clock_t te_ugape = clock();
                //sampleTime = (te_ugape - tsSample) / (double)CLOCKS_PER_SEC;
                //cout << "UGapE: " << (te_ugape - ts_ugape) / (double)CLOCKS_PER_SEC << endl;
                //cout << "initTime = " << initTime << endl;
                //cout << "sampleTime = " << sampleTime << endl;
                //cout << "totalSample = " << totalSample << endl;
                //cout << "sortTime = " << sortTime << endl;
                ugape_timer.end();
                avg_mab_time += ugape_timer.timeCost;
                stringstream ssout;
                ssout << "ugape/" << filelabel << "/k=" << k << "/";
                mkpath(ssout.str());
                ssout << queryNode << ".txt";
                ofstream fout(ssout.str());
                fout.precision(15);
                for (int i = 0; i < ugapeAnswer.size(); i++) {
                    fout << ugapeAnswer[i].first << "\t" << ugapeAnswer[i].second << endl;
                }
                cout << "ugape_timer.timeCost = " << ugape_timer.timeCost << " , sortTime = " << sortTime << " , sampleTime = " << sampleTime << " , adaptPushTime = " << adaptPushTime << endl;
                cout << "adaptDeterQueryTime = " << adaptDeterQueryTime << " , adaptRandQueryTime = " << adaptRandQueryTime << endl;
                return ugapeAnswer;
            }
            // sample high
            sample_timer.start();
            for (int i = 0; i < min(batch, (int)(vecLB.size())); i++) {
                int armId = vecLB[i].first;
                sampleOneArm_Naive(queryNode, armId, batchSample);
                //sampleOneArm(queryNode, armId, batchSample);
            }
            // sample low
            for (int i = 0; i < min(batch, (int)(vecUB.size())); i++) {
                int armId = vecUB[i].first;
                sampleOneArm_Naive(queryNode, armId, batchSample);
                //sampleOneArm(queryNode, armId, batchSample);
            }
            sample_timer.end();
            sampleTime += sample_timer.timeCost;
        }
    }

    void adaptivePush(int node, int budget, int qNode) {
        cout << "adaptivePush : node = " << node << " , prev_budget = " << budgetMap[node] << " , budget = " << budget << " , qNode = " << qNode << endl;
        double cur_rmax = 1.0 / budget * 1; // tunable
        //cout << "cur_rmax = " << cur_rmax << endl;
        //if (adaptResidueMap.find(node) == adaptResidueMap.end())
            //cout << "impossible!" << endl;
        unordered_set<int> lvlSet;
        priority_queue<int, vector<int>, comp2> minQ;
        for (auto& item : adaptResidueMap[node]) {
            int level = item.first;
            if (lvlSet.find(level) == lvlSet.end()) {
                lvlSet.insert(level);
                minQ.push(level);
            }
        }
        //cout << "lvlSet.size() = " << lvlSet.size() << endl;
        //vector<int> removeList;
        while (!minQ.empty()) {
            int lvl = minQ.top();
            minQ.pop();
            unordered_map<int, double>& submap = adaptResidueMap[node][lvl];
            //cout << "lvl = " << lvl << " , submap.size() = " << submap.size() << endl;
            for (auto& item : submap) {
                int w = item.first;
                double residue_w = item.second;
                //cout << "w = " << w << " , residue_w = " << residue_w << endl;
                if ((g.getInSize(w) > 0) && (residue_w / g.getInSize(w) >= cur_rmax)) { //push
                    //cout << "push" << endl;
                    double reserve_w = residue_w * (1 - sqrtC);
                    if (adaptReserveMap[node][lvl][w] == 0)
                        adaptReserveMap[node][lvl][w] = reserve_w;
                    else
                        adaptReserveMap[node][lvl][w] += reserve_w;
                    double res_inc = residue_w * sqrtC / g.getInSize(w);
                    for (int j = 0; j < g.getInSize(w); j++) {
                        int nbr = g.getInVert(w, j);
                        if (adaptResidueMap[node][lvl + 1][nbr] == 0)
                            adaptResidueMap[node][lvl + 1][nbr] = res_inc;
                        else
                            adaptResidueMap[node][lvl + 1][nbr] += res_inc;

                        if ((g.getInSize(nbr) > 0) && (adaptResidueMap[node][lvl + 1][nbr] / g.getInSize(nbr) >= cur_rmax) && (lvlSet.find(lvl + 1) == lvlSet.end())) {
                            lvlSet.insert(lvl + 1);
                            minQ.push(lvl + 1);
                        }
                    }
                    adaptResidueMap[node][lvl][w] = 0;
                }
            }
        }
        // update budget
        budgetMap[node] = budget;
        double rsum_node = 0;
        unordered_map<int, unordered_map<int, double>> copySubmap; // remove items with residue = 0
        vector<pair<pair<int, int>, double> > aliasP;
        for (auto& item : adaptResidueMap[node]) {
            int lvl = item.first;
            unordered_map<int, double>& subsubmap = item.second;
            for (auto& subitem : subsubmap) {
                int w = subitem.first;
                double res_w = subitem.second;

                if (res_w > 0) {
                    //cout << "lvl = " << lvl << " , w = " << w << " , res_w = " << res_w << endl;
                    rsum_node += res_w;
                    aliasP.push_back(pair<pair<int, int>, double>(pair<int, int>(lvl, w), res_w));
                    copySubmap[lvl][w] = res_w;
                }
            }
        }
        adaptResidueMap[node] = copySubmap;
        //test 
        /*
        for (auto& item : adaptResidueMap[node]) {
            int lvl = item.first;
            unordered_map<int, double>& subsubmap = item.second;
            for (auto& subitem : subsubmap) {
                int w = subitem.first;
                double res_w = subitem.second;

                if (res_w == 0) {
                    cout << "error, still has 0" << endl;
                    exit(0);
                }
            }
        }
        */


        if (nodeToAliasPos.find(node) == nodeToAliasPos.end()) {
            int pos = aliasMap.size();
            nodeToAliasPos[node] = pos;
            aliasMap.push_back(Alias(aliasP)); // at position pos
        }
        else
            aliasMap[nodeToAliasPos[node]] = Alias(aliasP);

        rsumMap[node] = rsum_node;
        cout << "adaptivePush done" << endl;
        cout << "rsumMap[node] = " << rsumMap[node] << endl;

        // pre-compute some parts of the SimRank value
        // must guarantee that adaptReserveMap.find(queryNode) != adaptReserveMap.end() && adaptReserveMap.find(candNode) != adaptReserveMap.end()
        if (node != qNode) {
            precompute_part1(node, qNode);
            precompute_part2(node, qNode, part2Map, node);
            precompute_part2(qNode, node, part3Map, node);
        }
        else { // query node
            // for all candidates, precompute
            for (int i = 0; i < curIdx; i++) {
                int cnode = ansIdx[i];
                precompute_part1(cnode, qNode);
                precompute_part2(cnode, qNode, part2Map, cnode);
                precompute_part2(qNode, cnode, part3Map, cnode);
            }
        }
    }

    /* compute \sum{pi_l(cnode, w) * pi_l(qnode, w)} */
    void precompute_part1(int cnode, int qnode) {
        part1Map[cnode] = 0; // clear previous value
        double simPart1 = 0;
        unordered_map<int, unordered_map<int, double> >& qNodeMap = adaptReserveMap[qnode];
        unordered_map<int, unordered_map<int, double> >& cNodeMap = adaptReserveMap[cnode];
        for (auto& item : cNodeMap) {
            int level = item.first;
            unordered_map<int, double>& cNodeSubmap = item.second;
            for (auto& subitem : cNodeSubmap) {
                int w = subitem.first;
                double pi_l_w = subitem.second;
                if (qNodeMap[level][w] != 0) {
                    // sample dw
                    double dw;
                    int hitCount = sampleD(w, 1);
                    if (hitMap.find(w) == hitMap.end())
                        hitMap[w] = 0;
                    if (totalMap.find(w) == totalMap.end())
                        totalMap[w] = 0;
                    hitMap[w] += hitCount;
                    totalMap[w]++;
                    if (hitMap[w] < 0 || totalMap[w] < 0) {
                        cout << "sample d : long overflow" << endl;
                        exit(0);
                    }

                    if (g.getInSize(w) == 0)
                        dw = 1;
                    else if (hitMap.find(w) != hitMap.end())
                        dw = 1 - C_value / (double)g.getInSize(w) - (double)hitMap[w] / totalMap[w];
                    else
                        dw = 1 - C_value / (double)g.getInSize(w);

                    simPart1 += pi_l_w * qNodeMap[level][w] * dw;
                }
            }
        }
        simPart1 /= (1 - sqrtC) * (1 - sqrtC);
        part1Map[cnode] = simPart1;
    }

    /* compute \sum{pi_l(cnode, w) * delta_l(qnode, w)} and \sum{delta_l(cnode, w) * pi_l(qnode, w)} */
    void precompute_part2(int cnode, int qnode, unordered_map<int, double> &partMap, int key) {
        /* compute part 2/3 of SimRank value */
        partMap[key] = 0; // clear previous value
        double simPart2 = 0;
        int numSample = 10000; // argument
        //            level           -> w  -> value
        //unordered_map<int, unordered_map<int, double> > deltaCNodeMap; 
        unordered_map<int, unordered_map<int, double> > rwMap;
        for (auto& item : adaptResidueMap[cnode]) {
            int level = item.first;
            for (auto& subItem : item.second) {
                int x = subItem.first; 
                double rx = subItem.second;

                int numRW = 0;
                rwMap.clear(); // local
                if (rsumMap[cnode] > 0) {
                    numRW = (int)(numSample * rx / rsumMap[cnode]) + 1;
                    for (int i = 0; i < numRW; i++) {
                        int tempNode = x;
                        int walk_level = 0;
                        while (R.drand() < sqrtC && walk_level < log(1.0e-6) / log(C_value)) {
                            int length = g.getInSize(tempNode);
                            if (length == 0)
                                break;
                            int r = R.generateRandom() % length;
                            x = g.getInVert(x, r);
                            walk_level++;
                        }
                        if (rwMap[walk_level][tempNode] == 0)
                            rwMap[walk_level][tempNode] = 1 / (double)numRW;
                        else
                            rwMap[walk_level][tempNode] += 1 / (double)numRW;
                    }
                }
                for (auto& item : rwMap) {
                    int walk_level = item.first;
                    for (auto& subItem : item.second) {
                        int w = subItem.first;
                        double pi_x_w = subItem.second;
                        if (adaptReserveMap[qnode][level + walk_level][w] != 0) {
                            double dw;

                            int hitCount = sampleD(w, 1);
                            if (hitMap.find(w) == hitMap.end())
                                hitMap[w] = 0;
                            if (totalMap.find(w) == totalMap.end())
                                totalMap[w] = 0;
                            hitMap[w] += hitCount;
                            totalMap[w]++;

                            if (hitMap[w] < 0 || totalMap[w] < 0) {
                                cout << "sample d : long overflow" << endl;
                                exit(0);
                            }

                            if (g.getInSize(w) == 0)
                                dw = 1;
                            else if (hitMap.find(w) != hitMap.end())
                                dw = 1 - C_value / (double)g.getInSize(w) - (double)hitMap[w] / totalMap[w];
                            else
                                dw = 1 - C_value / (double)g.getInSize(w);

                            simPart2 += rx * pi_x_w * dw * adaptReserveMap[qnode][level + walk_level][w] / (1 - sqrtC) / (1 - sqrtC);
                        }
                    }
                }
            }
        }
        partMap[key] = simPart2;
    }

    void sampleOneArm(int queryNode, int candNode, int batchSample) {
        if (adaptReserveMap.find(queryNode) == adaptReserveMap.end() || adaptReserveMap.find(candNode) == adaptReserveMap.end()) {
            cout << "adaptReserveMap does not contain queryNode / candNode" << endl;
            exit(0);
        }
        cout << "invoke sampleOneArm : queryNode = " << queryNode << " , candNode = " << candNode << endl;
        double simPart1 = 0, simPart2 = 0, simPart3 = 0, simPart4 = 0;
        if (part1Map.find(candNode) != part1Map.end())
            simPart1 = part1Map[candNode];
        if (part2Map.find(candNode) != part2Map.end())
            simPart2 = part2Map[candNode];
        if (part3Map.find(candNode) != part3Map.end())
            simPart3 = part3Map[candNode];
        cout << "simPart1 = " << simPart1 << " , simPart2 = " << simPart2 << " , simPart3 = " << simPart3 << endl;

        Alias& aliasQNode = aliasMap[nodeToAliasPos[queryNode]];
        Alias& aliasCNode = aliasMap[nodeToAliasPos[candNode]];
        for (int sample = 0; sample < batchSample; sample++) {            
            pair<int, int> qNodePair = aliasQNode.generateRandom(R);
            int level_x = qNodePair.first;
            int node_x = qNodePair.second;
            while (R.drand() < sqrtC && level_x < 100) { // 
                int length = g.getInSize(node_x);
                if (length > 0) {
                    int r = R.generateRandom() % length;
                    node_x = g.getInVert(node_x, r);
                    level_x++;
                }
                else {
                    break;
                }
            }
            
            pair<int, int> cNodePair = aliasCNode.generateRandom(R);
            int level_y = cNodePair.first;
            int node_y = cNodePair.second;
            while (R.drand() < sqrtC && level_y < 100) { // 
                int length = g.getInSize(node_y);
                if (length > 0) {
                    int r = R.generateRandom() % length;
                    node_y = g.getInVert(node_y, r);
                    level_y++;
                }
                else {
                    break;
                }
            }
            if (level_x == level_y && node_x == node_y) {
                double dx;

                int hitCount = sampleD(node_x, 1);
                if (hitMap.find(node_x) == hitMap.end())
                    hitMap[node_x] = 0;
                if (totalMap.find(node_x) == totalMap.end())
                    totalMap[node_x] = 0;
                hitMap[node_x] += hitCount;
                totalMap[node_x]++;

                if (hitMap[node_x] < 0 || totalMap[node_x] < 0) {
                    cout << "sample d : long overflow" << endl;
                    exit(0);
                }

                if (g.getInSize(node_x) == 0)
                    dx = 1;
                else if (hitMap.find(node_x) != hitMap.end())
                    dx = 1 - C_value / (double)g.getInSize(node_x) - (double)hitMap[node_x] / totalMap[node_x];
                else
                    dx = 1 - C_value / (double)g.getInSize(node_x);
                simPart4 = rsumMap[queryNode] * rsumMap[candNode] * dx / (1 - sqrtC) / (1 - sqrtC);
            }
            double simValue = simPart1 + simPart2 + simPart3 + simPart4;
            answer[candNode] += simValue;
            ans2[candNode] += simValue * simValue;
        } // end for
        cntMABSample[queryNode] += batchSample;
        cntMABSample[candNode] += batchSample;
        totalMABSample += batchSample;

        Timer adapt_push_timer;
        adapt_push_timer.start();
        if (cntMABSample[queryNode] >= 2 * budgetMap[queryNode]) {
            cout << "queryNode = " << queryNode << " , cntMABSample[queryNode] = " << cntMABSample[queryNode] << " , budgetMap[queryNode] = " << budgetMap[queryNode] << endl;
            adaptivePush(queryNode, cntMABSample[queryNode], queryNode);
        }
        if (cntMABSample[candNode] >= 2 * budgetMap[candNode]) {
            cout << "candNode = " << candNode << " , cntMABSample[candNode] = " << cntMABSample[candNode] << " , budgetMap[candNode] = " << budgetMap[candNode] << endl;
            adaptivePush(candNode, cntMABSample[candNode], queryNode);
        }
        adapt_push_timer.end();
        adaptPushTime += adapt_push_timer.timeCost;
    }

    vector <pair<int, double> > UGapE_with_adptive_push(int queryNode, int k) {
        Timer ugape_timer;
        ugape_timer.start();
        double sortTime = 0, sampleTime = 0;
        adaptPushTime = 0; adaptDeterQueryTime = 0; adaptRandQueryTime = 0;
        Timer sort_timer, sample_timer;

        // first should check if there exist exact k candidates
        int numCand = 0;
        for (auto& item : sameInNbrMap)
            numCand += item.second.size();
        cout << "numCand = " << numCand << endl;
        if (numCand <= k) {
            cout << "less than k candidates, skip sample-one-arm phase..." << endl;
            vector <pair<int, double> > ugapeAnswer;
            for (auto& item : sameInNbrMap) {
                int armId = item.first;
                vector<int>& sameNbrVec = sameInNbrMap[armId];
                for (int j = 0; j < sameNbrVec.size(); j++) {
                    int tempNode = sameNbrVec[j];
                    ugapeAnswer.push_back(pair<int, double>(tempNode, answer[armId] / cntMABSample[armId]));
                    cout << "armId = " << tempNode << " , cntMABSample[armId] = " << cntMABSample[armId] << " , est = " << answer[armId] / cntMABSample[armId] << endl;
                }
            }
            sort(ugapeAnswer.begin(), ugapeAnswer.end(), maxScoreCmp);
            ugape_timer.end();
            avg_mab_time += ugape_timer.timeCost;
            stringstream ssout;
            ssout << "ugape/" << filelabel << "/k=" << k << "/";
            mkpath(ssout.str());
            ssout << queryNode << ".txt";
            ofstream fout(ssout.str());
            fout.precision(15);
            for (int i = 0; i < ugapeAnswer.size(); i++) {
                fout << ugapeAnswer[i].first << "\t" << ugapeAnswer[i].second << endl;
            }
            cout << "ugape_timer.timeCost = " << ugape_timer.timeCost << " , sortTime = " << sortTime << " , sampleTime = " << sampleTime << " , adaptPushTime = " << adaptPushTime << endl;
            cout << "adaptDeterQueryTime = " << adaptDeterQueryTime << " , adaptRandQueryTime = " << adaptRandQueryTime << endl;
            return ugapeAnswer;
        }

        // second
        double eps_min = 1e-3;
        vector<pair<int, double> > B;
        unordered_map<int, double> UB, LB;
        vector<pair<int, double> > vecUB, vecLB;
        unordered_map<int, int> rank;
        rsumMap.clear();
        aliasMap.clear();
        nodeToAliasPos.clear();       
        adaptReserveMap.clear();
        adaptResidueMap.clear();
        budgetMap.clear();
        aliasMap.clear();

        adaptResidueMap[queryNode][0][queryNode] = 1; // initialize
        for (int i = 0; i < curIdx; i++) {
            int armId = ansIdx[i];
            budgetMap[armId] = 0;
            adaptResidueMap[armId][0][armId] = 1; // initialize
            adaptivePush(armId, cntMABSample[armId], queryNode);
        }

        int round = 0;
        int batch = 10;
        int batchSample = 10000; // max(10000, (int)(candList.size()));
        cout << "batchSample = " << batchSample << endl;

        while (true) {
            round++;
            sort_timer.start();
            UB.clear();
            LB.clear();
            vecUB.clear();
            vecLB.clear();
            B.clear();
            for (int i = 0; i < curIdx; i++) {
                int armId = ansIdx[i];
                double ugapCB = ugapeCB(curIdx, answer[armId], ans2[armId], cntMABSample[armId], 1.0e-4, round);
                UB[armId] = answer[armId] / cntMABSample[armId] + ugapCB;
                LB[armId] = answer[armId] / cntMABSample[armId] - ugapCB;
                vecUB.push_back(pair<int, double>(armId, answer[armId] / cntMABSample[armId] + ugapCB));
                //cout << "armId = " << armId << " , ugapCB = " << ugapCB << " , UB[armId] = " << UB[armId] << " , LB[armId] = " << LB[armId] << endl;
            }
            sort(vecUB.begin(), vecUB.end(), maxScoreCmp); // descreasing order
            double kthUB; // k-th largest UB
            double k1thUB;    // k+1-th largest UB
            int accum_rank = 1; // ranking starting from 1
            for (int i = 0; i < vecUB.size(); i++) {
                int armId = vecUB[i].first;
                rank[armId] = accum_rank;
                accum_rank += sameInNbrMap[armId].size();
                //cout << "armId = " << armId << " , rank[armId] = " << rank[armId] << endl;
                if (rank[armId] <= k && accum_rank > k) {
                    kthUB = vecUB[i].second;
                    if (i == vecUB.size() - 1)
                        k1thUB = vecUB[i].second; // strange
                    else
                        k1thUB = vecUB[i + 1].second;
                }
            }
            //cout << "kthUB = " << kthUB << endl;
            //cout << "k1thUB = " << k1thUB << endl;

            for (int i = 0; i < vecUB.size(); i++) {
                int armId = vecUB[i].first;
                if (rank[armId] <= k)     //(i < k)
                    B.push_back(pair<int, double>(armId, k1thUB - LB[armId]));
                else
                    B.push_back(pair<int, double>(armId, kthUB - LB[armId]));
            }
            vecUB.clear();
            sort(B.begin(), B.end(), minScoreCmp); // increasing order
            accum_rank = 1; // ranking starting from 1
            rank.clear();
            double kthB;
            for (int i = 0; i < B.size(); i++) {
                int armId = B[i].first;
                rank[armId] = accum_rank;
                //cout << "armId = " << armId << " , rank[armId] = " << rank[armId] << endl;
                accum_rank += sameInNbrMap[armId].size();
                if (rank[armId] <= k && accum_rank > k)
                    kthB = B[i].second;
            }
            //cout << "kthB = " << kthB << endl;

            for (int i = 0; i < B.size(); i++) {
                int armId = B[i].first;
                if (rank[armId] <= k)
                    vecLB.push_back(pair<int, double>(armId, LB[armId]));
                else
                    vecUB.push_back(pair<int, double>(armId, UB[armId]));
            }
            sort(vecLB.begin(), vecLB.end(), minScoreCmp);
            sort(vecUB.begin(), vecUB.end(), maxScoreCmp);
            //cout << "vecLB.size() = " << vecLB.size() << " , vecUB.size() = " << vecUB.size() << endl;
            //cout << "vecLB contains:" << endl;
            //for (int i = 0; i < vecLB.size(); i++)
            //    cout << vecLB[i].first << " , " << endl;
            //cout << "vecUB contains:" << endl;
            //for (int i = 0; i < vecUB.size(); i++)
            //    cout << vecUB[i].first << " , " << endl;

            sort_timer.end();
            sortTime += sort_timer.timeCost;

            if (kthB <= eps_min) {
                vector <pair<int, double> > ugapeAnswer;
                for (int i = 0; i < vecLB.size(); i++) {
                    int armId = vecLB[i].first;
                    vector<int>& sameNbrVec = sameInNbrMap[armId];
                    for (int j = 0; j < sameNbrVec.size(); j++) {
                        int tempNode = sameNbrVec[j];
                        ugapeAnswer.push_back(pair<int, double>(tempNode, answer[armId] / cntMABSample[armId]));
                        cout << "armId = " << tempNode << " , cntMABSample[armId] = " << cntMABSample[armId] << " , est = " << answer[armId] / cntMABSample[armId] << endl;
                    }
                }
                sort(ugapeAnswer.begin(), ugapeAnswer.end(), maxScoreCmp);
                //for (int i = 0; i < candList.size(); i++)
                   //delete aliasMap[i];

                //clock_t te_ugape = clock();
                //sampleTime = (te_ugape - tsSample) / (double)CLOCKS_PER_SEC;
                //cout << "UGapE: " << (te_ugape - ts_ugape) / (double)CLOCKS_PER_SEC << endl;
                //cout << "initTime = " << initTime << endl;
                //cout << "sampleTime = " << sampleTime << endl;
                //cout << "totalSample = " << totalSample << endl;
                //cout << "sortTime = " << sortTime << endl;
                ugape_timer.end();
                avg_mab_time += ugape_timer.timeCost;
                stringstream ssout;
                ssout << "ugape/" << filelabel << "/k=" << k << "/";
                mkpath(ssout.str());
                ssout << queryNode << ".txt";
                ofstream fout(ssout.str());
                fout.precision(15);
                for (int i = 0; i < ugapeAnswer.size(); i++) {
                    fout << ugapeAnswer[i].first << "\t" << ugapeAnswer[i].second << endl;
                }
                cout << "ugape_timer.timeCost = " << ugape_timer.timeCost << " , sortTime = " << sortTime << " , sampleTime = " << sampleTime << " , adaptPushTime = " << adaptPushTime << endl;
                cout << "adaptDeterQueryTime = " << adaptDeterQueryTime << " , adaptRandQueryTime = " << adaptRandQueryTime << endl;
                return ugapeAnswer;
            }
            // sample high
            sample_timer.start();
            for (int i = 0; i < min(batch, (int)(vecLB.size())); i++) {
                int armId = vecLB[i].first;
                //sampleOneArm_Naive(queryNode, armId, batchSample);
                sampleOneArm(queryNode, armId, batchSample);
            }
            // sample low
            for (int i = 0; i < min(batch, (int)(vecUB.size())); i++) {
                int armId = vecUB[i].first;
                //sampleOneArm_Naive(queryNode, armId, batchSample);
                sampleOneArm(queryNode, armId, batchSample);
            }
            sample_timer.end();
            sampleTime += sample_timer.timeCost;
        }
    }
};



#endif
