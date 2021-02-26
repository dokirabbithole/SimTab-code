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
    // int* Parent; // 
    double *residue;
    double *newReserve;
    double *newResidue;
    bool* isInArray;

    // for estimation of D
    unordered_map<int, double> hitMap;
    unordered_map<int, double> totalMap;

	/* for prefiltering */
	vector<pair<int, double> > estMean;
	vector<pair<int, double> > prefCB;
	double *answer;
	double *ans2;
	int *ansIdx;
	int curIdx;

	double eps_p;
	int totalSample; 
	// time
	double t_pref, t_single;

    unordered_map<int, vector<int> > sameInNbrMap;
    //            indegree,          armId,     
    unordered_map<int, unordered_map<int, vector<int> > > sameInNbrMap2;
    bool* exist_cnode_nbr;

    /* data structures for MAB */
    //            query/cand   ->    level              -> w  -> \pi_l(s,w), r_l(s,w)
    unordered_map<int, unordered_map<int, unordered_map<int, double> > > adaptReserveMap, adaptResidueMap;
    unordered_map<int, unordered_map<int, double> > ppr;
    unordered_map<int, double> part1Map, part2Map, part3Map, part4Map, cntMap;

    vector<int> candList; // ansIdx + curIdx
    unordered_map<int, long> cntMABSample;
    long totalMABSample;
    unordered_map<int, int> budgetMap;
    unordered_map<int, double> rsumMap;
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
        //Parent = new int[vert]; //
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
            //Parent[i] = -1;
        }
        srand(unsigned(time(0)));
        
		// initialization
		answer = new double[vert];
		ans2 = new double[vert];
        ansIdx = new int[vert];
        for(int i = 0; i < vert; i++){
            answer[i] = 0;
			ans2[i] = 0;
            isInArray[i] = false;
        }		
		curIdx = 0;
        cout << "====init done!====" << endl;
        exist_cnode_nbr = new bool[vert];
        for (int i = 0; i < vert; i++)
            exist_cnode_nbr[i] = false;
    }
    ~SimStruct() {
        delete[] H[0];
        delete[] H[1];
        delete[] U[0];
        delete[] U[1];
        delete[] C[0];
        delete[] UC[0];
        //delete[] Parent; 
        delete[] residue;
        delete[] newReserve;
        delete[] newResidue;
        delete[] isInArray;
		// for prefiltering
		delete[] answer;
		delete[] ans2;
		delete[] ansIdx;

        delete[] exist_cnode_nbr;
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
                //realAnswer[c.first].push_back(pair<int, double>(d.first, d.second));
                numFORAItems++;
                if (d.second > 1.0e-6) {
                    numFORAItemsPruned++;
                    realAnswer[c.first].push_back(pair<int, double>(d.first, d.second));
                }
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
		eps_p = 1.0e-6;
		epsilon = 0.5;
		estMean.clear();
		prefCB.clear();
		
        /* iteratively apply sample-all-arms, i.e., singleSourceQuery(u) */
		do{
			//cout << "epsilon = " << epsilon << endl;
			rmax = (1 - sqrtC) * (1 - sqrtC) / 25 * epsilon;
			back_nr = (int)(20 * log(vert) / epsilon / epsilon);
			nr = (int)(0.1 * back_nr); 
			forward_rmax = rmax; // rmax * 2
			backward_rmax = rmax; // not important; not used
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
                for (int i = 0; i < curIdx; i++) {
                    int nodeId = ansIdx[i];
                    if (nodeId != u) {
                        double score = answer[nodeId] / totalSample;
                        double cb = computeCB(nodeId, ans2, score, totalSample, vert);
                        lbList.push_back(pair<int, double>(nodeId, score - cb));
                        ubMap[nodeId] = score + cb;
                    }
                }
                sort(lbList.begin(), lbList.end(), maxScoreCmp); // in descending order, O(nlogn)
                //for (int i = 0; i < lbList.size(); i++)
                //    cout << lbList[i].first << " , " << lbList[i].second << " , " << estMap[lbList[i].first] << endl;
                int kthNode = lbList[k - 1].first;
                double kthLB = lbList[k - 1].second;
                double cb_0 = 3 * C_value * log(1.0 / 0.0001) / totalSample;
                cout << "kthNode = " << kthNode << " , kthLB = " << kthLB << " , cb_0 = " << cb_0 << endl;
                if (kthLB + eps_p > cb_0) { // can prune 0
                    vector<int> candidates; // store candidate nodes
                    int count_cand = k;
                    for (int i = 0; i < k; i++)
                        candidates.push_back(lbList[i].first);
                    for (int i = k; i < lbList.size(); i++) {
                        int cand = lbList[i].first;
                        if (kthLB + eps_p < ubMap[cand]) {
                            candidates.push_back(cand);
                            count_cand++;
                        }
                        else if(i < 5000)
                            candidates.push_back(cand); // not counted
                    }
                    cout << "numCandidate = " << candidates.size() << endl;
                    if (t_pref * max(2 * k, 500) > t_single * count_cand) {
                        for (int i = 0; i < candidates.size(); i++) { // candidates
                            int cand = candidates[i];
                            estMean.push_back(pair<int, double>(cand, answer[cand] / totalSample));
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
		ssout << "candidates/" << filelabel << "/k=" << k << "/";
		mkpath(ssout.str());
		ssout << u << ".txt";
		ofstream fout(ssout.str());
		fout.precision(15);
		for(int i = 0; i < estMean.size(); i++){
            int nodeId = estMean[i].first;
            double ithCB = computeCB(nodeId, ans2, estMean[i].second, totalSample, vert);
            fout << estMean[i].first << "\t" << estMean[i].second << "\t" << totalSample << "\t" << ithCB << "\t" << estMean.size() << endl;
		}
	}
	
	// -----------------------------------------------------------------------------------------------

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
    
    

    /* confidence bound */
    double ugapeCB(int numOfArms, double answer, double sum_x2, long numSample, double failProb, double round) {
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
        return (part1 + part2);
    }

    bool checkSameInNbr(int cnode) {
        bool found = false;
        int ind = g.getInSize(cnode);
        if (sameInNbrMap2.find(ind) == sameInNbrMap2.end()) {
            unordered_map<int, vector<int> > subMap;
            vector<int> tempVec;
            tempVec.push_back(cnode);
            subMap[cnode] = tempVec;
            sameInNbrMap2[ind] = subMap;
            return found;
        }
        else {
            unordered_map<int, vector<int> >& subMap = sameInNbrMap2[ind];
            for (auto& a : subMap) {
                int exist_node = a.first;
                vector<int>& v = a.second;
                if (hasSameInNbr(exist_node, cnode)) {
                    v.push_back(cnode);
                    found = true;
                    break;
                }
            }
            if (!found) {
                vector<int> v;
                v.push_back(cnode);
                subMap[cnode] = v;
            }
        }

        /*
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
        */
        
        return found;
    }

    bool hasSameInNbr(int exist_cnode, int cnode) {
        //if (g.getInSize(exist_cnode) != g.getInSize(cnode))
        //    return false;
        //unordered_set<int> exist_cnode_nbr;
        /* exist_cnode and cnode has same indegree */
        for (int i = 0; i < g.getInSize(exist_cnode); i++)
            exist_cnode_nbr[g.getInVert(exist_cnode, i)] = true;
        bool isSame = true;
        for (int i = 0; i < g.getInSize(cnode); i++) {
            int tempNode = g.getInVert(cnode, i);
            if (exist_cnode_nbr[tempNode] == false) {
                isSame = false;
                break;
            }
        }
        for (int i = 0; i < g.getInSize(exist_cnode); i++)
            exist_cnode_nbr[g.getInVert(exist_cnode, i)] = false;
        return isSame;
    }

    

	void readCandidatesFromPrefilter3(int qnode, int k) {
        stringstream ss_cand_path;
        ss_cand_path << "candidates/" << filelabel << "/k=" << k << "/";
        string candpath = ss_cand_path.str();
		// candList.clear(); // do not use candList anymore
		totalMABSample = 0;
		cntMABSample.clear();
		sameInNbrMap.clear();
        sameInNbrMap2.clear(); //
		nodeToAliasPos.clear(); // newly added
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
        cntMABSample[qnode] = 0; // totalSample;
        int countRead = 0;
		while (ifcand.good() && countRead < max(5 * k, 500)) { //   && countRead < max(5 * k, 500) // for testing MC
            countRead++;
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
					answer[cnode] = 0; // estScore * cntAllArmsSample;
					cntMABSample[cnode] = 0; // cntAllArmsSample;
					budgetMap[cnode] = cntAllArmsSample;
					totalMABSample = 0; // += cntAllArmsSample;
					ans2[cnode] = 0; // ans2val;
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
		cout << "qnode = " << qnode << " , candidate size = " << curIdx << " , sameInNbrMap.size() = " << sameInNbrMap.size() << endl;
		cout << "chk_nbr_time = " << chk_nbr_time << endl;
		avg_chk_nbr_time += chk_nbr_time;
		//budgetMap[qnode] = 1;
		//adaptResidueMap[qnode][0][qnode] = 1;

		//
		//for (int i = 0; i < candList.size(); i++) {
		//    int nodeId = candList[i];
		//    ansIdx[curIdx] = nodeId;
		//    curIdx++;
		//}
	}

    void MC(int queryNode, int k) {
        Timer mcTimer;
        mcTimer.start();
        double eps_mc = 1.0e-2;
        int numSample = (int)(log(curIdx) / eps_mc / eps_mc);
        cout << "MC : numSample = " << numSample << endl;
        vector<pair<int, double> > mcAnswer;
        for (auto& item : sameInNbrMap2) {
            unordered_map<int, vector<int>> subMap = item.second;
            for (auto& subItem : subMap) {
                int cnode = subItem.first;
                sampleOneArm_Naive(queryNode, cnode, numSample);
                vector<int>& sameNbrVec = subItem.second;
                for (int j = 0; j < sameNbrVec.size(); j++) {
                    int tempNode = sameNbrVec[j];
                    mcAnswer.push_back(pair<int, double>(tempNode, answer[cnode] / cntMABSample[cnode]));
                    cout << "cnode = " << tempNode << " , cntMABSample[tempNode] = " << cntMABSample[cnode] << " , est = " << answer[cnode] / cntMABSample[cnode] << endl;
                }
            }
        }
        mcTimer.end();
        sort(mcAnswer.begin(), mcAnswer.end(), maxScoreCmp);
        stringstream ssout;
        ssout << "MC/" << filelabel << "/k=" << k << "/";
        mkpath(ssout.str());
        ssout << queryNode << ".txt";
        ofstream fout(ssout.str());
        fout.precision(15);
        for (int i = 0; i < mcAnswer.size(); i++) {
            fout << mcAnswer[i].first << "\t" << mcAnswer[i].second << endl;
        }
        fout.close();
        cout << "mcTimer.timeCost = " << mcTimer.timeCost << " , # total sample = " << numSample * curIdx << endl;
    }

    vector <pair<int, double> > UGapE2(int queryNode, int k) {
        // TEST
        // read ground truth
        unordered_map<int, double> gtMap;
        stringstream ss_gt_in;
        ss_gt_in << "/home/liuyu/SimRace-liuyu/GroundTruth_ExactSim/" << filelabel << "/" << queryNode << "+0.0000001000.txt";
        ifstream if_gt(ss_gt_in.str());
        while (if_gt.good()) {
            int node;
            double gtScore;
            if_gt >> node >> gtScore;
            gtMap[node] = gtScore;
        }
        if_gt.close();
        cout << "reading ground truth done" << endl;

        Timer ugape_timer;
        ugape_timer.start();
        double sortTime = 0, sampleTime = 0;
        Timer sort_timer, sample_timer;

        // first should check if there exist exact k candidates
        int numCand = 0;
        for (auto& item : sameInNbrMap2) {
            for(auto& subItem : item.second)
                numCand += subItem.second.size();
        }
        cout << "numCand = " << numCand << endl;
        if (numCand <= k) {
            cout << "less than k candidates, skip sample-one-arm phase..." << endl;
            vector <pair<int, double> > ugapeAnswer;
            for (auto& item : sameInNbrMap2) {
                unordered_map<int, vector<int>>& subMap = item.second;
                for (auto& subItem : subMap) {
                    int armId = subItem.first;
                    vector<int>& sameNbrVec = subItem.second;
                    for (int j = 0; j < sameNbrVec.size(); j++) {
                        int tempNode = sameNbrVec[j];
                        ugapeAnswer.push_back(pair<int, double>(tempNode, answer[armId] / cntMABSample[armId]));
                        cout << "armId = " << tempNode << " , cntMABSample[armId] = " << cntMABSample[armId] << " , est = " << answer[armId] / cntMABSample[armId] << endl;
                    }
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
            cout << "ugape_timer.timeCost = " << ugape_timer.timeCost << " , sortTime = " << sortTime << " , sampleTime = " << sampleTime << endl;
            return ugapeAnswer;
        }

        // second
        double eps_min = 1e-5;
        vector<pair<int, double> > B;
        unordered_map<int, double> UB, LB;
        vector<pair<int, double> > vecUB, vecLB;
        unordered_map<int, int> rank;        
        
        // mean, cntSample, and sum_Xi2 already have values
        double round = 0;
        int batch = 10;
        int batchSample = 10000; 
        cout << "batchSample = " << batchSample << endl;

        for (int i = 0; i < curIdx; i++) {
            int armId = ansIdx[i];
            sampleOneArm_Naive(queryNode, armId, 1000);
        }

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
                accum_rank += sameInNbrMap2[g.getInSize(armId)][armId].size();
                //cout << "armId = " << armId << " , rank[armId] = " << rank[armId] << endl;
                if (rank[armId] <= k && accum_rank > k) {
                    kthUB = vecUB[i].second;
                    if(i == vecUB.size() - 1)
                        k1thUB = vecUB[i].second; // strange
                    else
                        k1thUB = vecUB[i + 1].second;
                }
            }

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
                accum_rank += sameInNbrMap2[g.getInSize(armId)][armId].size();
                if (rank[armId] <= k && accum_rank > k)
                    kthB = B[i].second;
            }
            
            for (int i = 0; i < B.size(); i++) {
                int armId = B[i].first;
                if (rank[armId] <= k)
                    vecLB.push_back(pair<int, double>(armId, LB[armId]));
                else
                    vecUB.push_back(pair<int, double>(armId, UB[armId]));
            }
            sort(vecLB.begin(), vecLB.end(), minScoreCmp);
            sort(vecUB.begin(), vecUB.end(), maxScoreCmp);
            sort_timer.end();
            sortTime += sort_timer.timeCost;

            if (kthB <= eps_min) { 
                vector <pair<int, double> > ugapeAnswer;
                for (int i = 0; i < vecLB.size(); i++) {
                    int armId = vecLB[i].first;
                    vector<int>& sameNbrVec = sameInNbrMap2[g.getInSize(armId)][armId];
                    for (int j = 0; j < sameNbrVec.size(); j++) {
                        int tempNode = sameNbrVec[j];
                        ugapeAnswer.push_back(pair<int, double>(tempNode, answer[armId] / cntMABSample[armId]));
                        cout << "armId = " << tempNode << " , cntMABSample[armId] = " << cntMABSample[armId] << " , est = " << answer[armId] / cntMABSample[armId] << endl;
                    }
                }
                // TEST
                for (int i = 0; i < vecUB.size(); i++) {
                    int armId = vecUB[i].first;
                    vector<int>& sameNbrVec = sameInNbrMap2[g.getInSize(armId)][armId];
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
                    fout << ugapeAnswer[i].first << "\t" << ugapeAnswer[i].second << "\t" 
                         << ugapeCB(curIdx, answer[ugapeAnswer[i].first], ans2[ugapeAnswer[i].first], cntMABSample[ugapeAnswer[i].first], 1.0e-4, round) 
                         << "\t" << ugapeAnswer[i].second  - gtMap[ugapeAnswer[i].first] << "\t" << cntMABSample[ugapeAnswer[i].first] << endl;
                }
                cout << "ugape_timer.timeCost = " << ugape_timer.timeCost << " , sortTime = " << sortTime << " , sampleTime = " << sampleTime << endl;
                return ugapeAnswer;
            }
            // sample high
            sample_timer.start();
            for (int i = 0; i < min(batch, (int)(vecLB.size())); i++) {
                int armId = vecLB[i].first;
                sampleOneArm_Naive(queryNode, armId, batchSample);
            }
            // sample low
            for (int i = 0; i < min(batch, (int)(vecUB.size())); i++) {
                int armId = vecUB[i].first;
                sampleOneArm_Naive(queryNode, armId, batchSample);
            }
            sample_timer.end();
            sampleTime += sample_timer.timeCost;
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
                    int total_sample_d = budgetMap[cnode];
                    //cout << "totalMap[w] = " << totalMap[w] << " , total_sample_d * ppr[cnode][w] * ppr[qnode][w] = " << total_sample_d * ppr[cnode][w] * ppr[qnode][w] << endl;
                    if (totalMap[w] < max(total_sample_d * ppr[cnode][w] * ppr[qnode][w], 1.0)) {
                        for (int count = 0; count < max(total_sample_d * ppr[cnode][w] * ppr[qnode][w] - totalMap[w], 1.0); count++) {
                            int hitCount = sampleD(w, 1);
                            if (hitMap.find(w) == hitMap.end())
                                hitMap[w] = 0;
                            if (totalMap.find(w) == totalMap.end())
                                totalMap[w] = 0;
                            hitMap[w] += hitCount;
                            totalMap[w]++;
                        }
                    }

                    double dw;
                    if (g.getInSize(w) == 0)
                        dw = 1;
                    else if (hitMap.find(w) != hitMap.end())
                        dw = 1 - C_value / (double)g.getInSize(w) - hitMap[w] / totalMap[w];
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
		double simPart2 = 0;
        int numSample = max(min(budgetMap[cnode], budgetMap[qnode]) / 5000, 1000); // max(min(budgetMap[cnode], budgetMap[qnode]) / 500, 1000); // argument
		unordered_map<int, unordered_map<int, double> > rwMap;
		for (auto& item : adaptResidueMap[cnode]) {
			int level = item.first;
			for (auto& subItem : item.second) {
				int x = subItem.first;
				double rx = subItem.second;
				rwMap.clear(); // local
				int numRW = (int)(numSample * rx / rsumMap[cnode]) + 1;
				for (int i = 0; i < numRW; i++) {
					int tempNode = x;
					int walk_level = 0;
					while (R.drand() < sqrtC && walk_level < log(1.0e-6) / log(C_value)) {
						int length = g.getInSize(tempNode);
						if (length == 0)
							break;
						int r = R.generateRandom() % length;
						tempNode = g.getInVert(tempNode, r);
						walk_level++;
					}
					if (rwMap[walk_level][tempNode] == 0) {
						rwMap[walk_level][tempNode] = 1 / (double)numRW;
					}
					else
						rwMap[walk_level][tempNode] += 1 / (double)numRW;
				}
                for (auto& item : rwMap) {
                    int l = item.first;
                    for (auto& subItem : item.second) {
                        int w = subItem.first;
                        double pi_x_w = subItem.second;

                        if (adaptReserveMap[qnode][level + l][w] != 0) {
                            int total_sample_d = 100; // max(min(budgetMap[cnode], budgetMap[qnode]) / 500, 10000);
                            if (totalMap[w] < max(total_sample_d * ppr[cnode][w] * ppr[qnode][w], 1.0)) {
                                for (int count = 0; count < max(total_sample_d * ppr[cnode][w] * ppr[qnode][w] - totalMap[w], 1.0); count++) {
                                    int hitCount = sampleD(w, 1);
                                    if (hitMap.find(w) == hitMap.end())
                                        hitMap[w] = 0;
                                    if (totalMap.find(w) == totalMap.end())
                                        totalMap[w] = 0;
                                    hitMap[w] += hitCount;
                                    totalMap[w]++;
                                }
                            }
                            double dw;
                            if (g.getInSize(w) == 0)
                                dw = 1;
                            else if (hitMap.find(w) != hitMap.end())
                                dw = 1 - C_value / (double)g.getInSize(w) - (double)hitMap[w] / totalMap[w];
                            else
                                dw = 1 - C_value / (double)g.getInSize(w);

                            simPart2 += pi_x_w * rx * dw * adaptReserveMap[qnode][level + l][w] / (1 - sqrtC) / (1 - sqrtC);
                        }
                    }
                }
			}
		}
		partMap[key] = simPart2;

		
		/*
		if (nodeToAliasPos.find(cnode) != nodeToAliasPos.end() && adaptReserveMap.find(qnode) != adaptReserveMap.end()) {
			double simPart2 = 0;
			int numSample = 10000;
			Alias& aliasNode = aliasMap[nodeToAliasPos[cnode]];
			for (int sample = 0; sample < numSample; sample++) {
				pair<int, int> nodePair = aliasNode.generateRandom(R);
				int level_x = nodePair.first;
				int node_x = nodePair.second;
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
				//TEST
				if (globalMap2[level_x][node_x] == 0)
					globalMap2[level_x][node_x] = rsumMap[cnode] / numSample;
				else
					globalMap2[level_x][node_x] += rsumMap[cnode] / numSample;

				if (adaptReserveMap[qnode][level_x][node_x] != 0) {
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
					simPart2 += rsumMap[cnode] * adaptReserveMap[qnode][level_x][node_x] * dx / (1 - sqrtC) / (1 - sqrtC);

				}
			}
			simPart2 /= numSample;
			partMap[key] = simPart2;
		}
		*/
		//TEST
		/*
		cout << "=====================================TEST=====================================" << endl;
		for (auto& item : globalMap2) {
			int lvl = item.first;
			cout << "lvl = " << lvl << " , (2)" << globalMap2[lvl].size() << " , (1)" << globalMap1[lvl].size() << endl;
			//for (auto& subItem : item.second) {
			//	int w = subItem.first;
			//	double p = subItem.second;
			//	cout << "\t w = " << w << " , p2 = " << p << " , p1 = " << globalMap1[lvl][w] << endl;
			//}
		}
		cout << "==============================================================================" << endl;
		*/
    }

    void sampleOneArm(int queryNode, int candNode, int batchSample) {
        if (adaptReserveMap.find(queryNode) == adaptReserveMap.end() || adaptReserveMap.find(candNode) == adaptReserveMap.end()) {
            cout << "adaptReserveMap does not contain queryNode / candNode" << endl;
            exit(0);
        }
        //cout << "invoke sampleOneArm : queryNode = " << queryNode << " , candNode = " << candNode << endl;
        double simPart1 = 0, simPart2 = 0, simPart3 = 0, simPart4 = 0;
        if (part1Map.find(candNode) != part1Map.end())
            simPart1 = part1Map[candNode];
        if (part2Map.find(candNode) != part2Map.end())
            simPart2 = part2Map[candNode];
        if (part3Map.find(candNode) != part3Map.end())
            simPart3 = part3Map[candNode];
        //cout << "simPart1 = " << simPart1 << " , simPart2 = " << simPart2 << " , simPart3 = " << simPart3 << endl;

        Alias& aliasQNode = aliasMap[nodeToAliasPos[queryNode]];
        Alias& aliasCNode = aliasMap[nodeToAliasPos[candNode]];
        /*
        for (int sample = 0; sample < batchSample / 10; sample++) {
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
            if (adaptReserveMap[candNode][level_x][node_x] != 0) {
                int total_sample_d = 1000 * 1000 * 100;
                if (totalMap[node_x] < max(total_sample_d * ppr[candNode][node_x] * ppr[queryNode][node_x], 1.0)) {
                    for (int count = 0; count < max(total_sample_d * ppr[candNode][node_x] * ppr[queryNode][node_x] - totalMap[node_x], 1.0); count++) {
                        int hitCount = sampleD(node_x, 1);
                        if (hitMap.find(node_x) == hitMap.end())
                            hitMap[node_x] = 0;
                        if (totalMap.find(node_x) == totalMap.end())
                            totalMap[node_x] = 0;
                        hitMap[node_x] += hitCount;
                        totalMap[node_x]++;
                    }
                }
                double dx;
                if (g.getInSize(node_x) == 0)
                    dx = 1;
                else if (hitMap.find(node_x) != hitMap.end())
                    dx = 1 - C_value / (double)g.getInSize(node_x) - (double)hitMap[node_x] / totalMap[node_x];
                else
                    dx = 1 - C_value / (double)g.getInSize(node_x);
                simPart2 = rsumMap[queryNode] * adaptReserveMap[candNode][level_x][node_x] * dx / (1 - sqrtC) / (1 - sqrtC);
            }
            if (adaptReserveMap[queryNode][level_y][node_y] != 0) {
                int total_sample_d = 1000 * 1000 * 100;
                if (totalMap[node_y] < max(total_sample_d * ppr[candNode][node_y] * ppr[queryNode][node_y], 1.0)) {
                    for (int count = 0; count < max(total_sample_d * ppr[candNode][node_y] * ppr[queryNode][node_y] - totalMap[node_y], 1.0); count++) {
                        int hitCount = sampleD(node_y, 1);
                        if (hitMap.find(node_y) == hitMap.end())
                            hitMap[node_y] = 0;
                        if (totalMap.find(node_y) == totalMap.end())
                            totalMap[node_y] = 0;
                        hitMap[node_y] += hitCount;
                        totalMap[node_y]++;
                    }
                }
                double dy;
                if (g.getInSize(node_y) == 0)
                    dy = 1;
                else if (hitMap.find(node_y) != hitMap.end())
                    dy = 1 - C_value / (double)g.getInSize(node_y) - (double)hitMap[node_y] / totalMap[node_y];
                else
                    dy = 1 - C_value / (double)g.getInSize(node_y);
                simPart3 = rsumMap[candNode] * adaptReserveMap[queryNode][level_y][node_y] * dy / (1 - sqrtC) / (1 - sqrtC);
            }
            part2Map[candNode] += simPart2;
            part3Map[candNode] += simPart3;
            cntMap[candNode]++;
        }
        */
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
            /*
            if (adaptReserveMap[candNode][level_x][node_x] != 0) {
                int total_sample_d = 1000 * 100;
                if (totalMap[node_x] < max(total_sample_d * ppr[candNode][node_x] * ppr[queryNode][node_x], 1.0)) {
                    for (int count = 0; count < max(total_sample_d * ppr[candNode][node_x] * ppr[queryNode][node_x] - totalMap[node_x], 1.0); count++) {
                        int hitCount = sampleD(node_x, 1);
                        if (hitMap.find(node_x) == hitMap.end())
                            hitMap[node_x] = 0;
                        if (totalMap.find(node_x) == totalMap.end())
                            totalMap[node_x] = 0;
                        hitMap[node_x] += hitCount;
                        totalMap[node_x]++;
                    }
                }
                double dx;
                if (g.getInSize(node_x) == 0)
                    dx = 1;
                else if (hitMap.find(node_x) != hitMap.end())
                    dx = 1 - C_value / (double)g.getInSize(node_x) - (double)hitMap[node_x] / totalMap[node_x];
                else
                    dx = 1 - C_value / (double)g.getInSize(node_x);
                simPart2 = rsumMap[queryNode] * adaptReserveMap[candNode][level_x][node_x] * dx / (1 - sqrtC) / (1 - sqrtC);
            }
            if (adaptReserveMap[queryNode][level_y][node_y] != 0) {
                int total_sample_d = 1000 * 100;
                if (totalMap[node_y] < max(total_sample_d * ppr[candNode][node_y] * ppr[queryNode][node_y], 1.0)) {
                    for (int count = 0; count < max(total_sample_d * ppr[candNode][node_y] * ppr[queryNode][node_y] - totalMap[node_y], 1.0); count++) {
                        int hitCount = sampleD(node_y, 1);
                        if (hitMap.find(node_y) == hitMap.end())
                            hitMap[node_y] = 0;
                        if (totalMap.find(node_y) == totalMap.end())
                            totalMap[node_y] = 0;
                        hitMap[node_y] += hitCount;
                        totalMap[node_y]++;
                    }
                }
                double dy;
                if (g.getInSize(node_y) == 0)
                    dy = 1;
                else if (hitMap.find(node_y) != hitMap.end())
                    dy = 1 - C_value / (double)g.getInSize(node_y) - (double)hitMap[node_y] / totalMap[node_y];
                else
                    dy = 1 - C_value / (double)g.getInSize(node_y);
                simPart3 = rsumMap[candNode] * adaptReserveMap[queryNode][level_y][node_y] * dy / (1 - sqrtC) / (1 - sqrtC);
            }
            */
            if (level_x == level_y && node_x == node_y) {
                
                int total_sample_d = 1000 * 10;
                if (totalMap[node_x] < 1) {
                    for (int count = 0; count < max(total_sample_d * ppr[candNode][node_x] * ppr[queryNode][node_x] - totalMap[node_x], 1.0); count++) {
                        int hitCount = sampleD(node_x, 1);
                        if (hitMap.find(node_x) == hitMap.end())
                            hitMap[node_x] = 0;
                        if (totalMap.find(node_x) == totalMap.end())
                            totalMap[node_x] = 0;
                        hitMap[node_x] += hitCount;
                        totalMap[node_x]++;
                    }
                }
                double dx;
                if (g.getInSize(node_x) == 0)
                    dx = 1;
                else if (hitMap.find(node_x) != hitMap.end())
                    dx = 1 - C_value / (double)g.getInSize(node_x) - (double)hitMap[node_x] / totalMap[node_x];
                else
                    dx = 1 - C_value / (double)g.getInSize(node_x);
                simPart4 = rsumMap[queryNode] * rsumMap[candNode] * dx / (1 - sqrtC) / (1 - sqrtC);
            }
			double simValue = simPart1 + simPart2 + simPart3 + simPart4;
            //double simValue = simPart1 + part2Map[candNode] / cntMap[candNode] + part3Map[candNode] / cntMap[candNode] + simPart4;
			part4Map[candNode] += simPart4;
            //simPart2 = 0; simPart3 = 0; 
            simPart4 = 0; // clear the value            
            answer[candNode] += simValue;
            ans2[candNode] += simValue * simValue;
        } // end for
        cntMABSample[queryNode] += batchSample;
        cntMABSample[candNode] += batchSample;
        totalMABSample += batchSample;

        Timer adapt_push_timer;
        adapt_push_timer.start();
        
        //if (cntMABSample[queryNode] >= budgetMap[queryNode] && budgetMap[queryNode] < 1e6) {
        //    cout << "queryNode = " << queryNode << " , cntMABSample[queryNode] = " << cntMABSample[queryNode] << " , budgetMap[queryNode] = " << budgetMap[queryNode] << endl;
        //    forwardPush(queryNode, budgetMap[queryNode] * 10, queryNode);
        //}
        //if (cntMABSample[candNode] >= max(budgetMap[candNode], 100000) && budgetMap[candNode] < 1e6) {
        //    cout << "candNode = " << candNode << " , cntMABSample[candNode] = " << cntMABSample[candNode] << " , budgetMap[candNode] = " << budgetMap[candNode] << endl;
		//	forwardPush(candNode, max(budgetMap[candNode], 100000) * 10, queryNode);
        //}
        
        adapt_push_timer.end();
        adaptPushTime += adapt_push_timer.timeCost;
    }

    vector <pair<int, double> > UGapE_with_adptive_push(int queryNode, int k) {
        // TEST
        // read ground truth
        unordered_map<int, double> gtMap;
        stringstream ss_gt_in;
        ss_gt_in << "/home/liuyu/SimRace-liuyu/GroundTruth_ExactSim/" << filelabel << "/" << queryNode << "+0.0000001000.txt";
        ifstream if_gt(ss_gt_in.str());
        while (if_gt.good()) {
            int node;
            double gtScore;
            if_gt >> node >> gtScore;
            gtMap[node] = gtScore;
        }
        if_gt.close();
        cout << "reading ground truth done" << endl;

        Timer ugape_timer;
        ugape_timer.start();
        double sortTime = 0, sampleTime = 0;
        adaptPushTime = 0; adaptDeterQueryTime = 0; adaptRandQueryTime = 0;
        Timer sort_timer, sample_timer;

        // first should check if there exist exact k candidates
        int numCand = 0;
        for (auto& item : sameInNbrMap2) {
            for (auto& subItem : item.second)
                numCand += subItem.second.size();
        }
        cout << "numCand = " << numCand << endl;
        if (numCand <= k) {
            cout << "less than k candidates, skip sample-one-arm phase..." << endl;
            vector <pair<int, double> > ugapeAnswer;
            for (auto& item : sameInNbrMap2) {
                unordered_map<int, vector<int>>& subMap = item.second;
                for (auto& subItem : subMap) {
                    int armId = subItem.first;
                    vector<int>& sameNbrVec = subItem.second;
                    for (int j = 0; j < sameNbrVec.size(); j++) {
                        int tempNode = sameNbrVec[j];
                        ugapeAnswer.push_back(pair<int, double>(tempNode, answer[armId] / cntMABSample[armId]));
                        cout << "armId = " << tempNode << " , cntMABSample[armId] = " << cntMABSample[armId] << " , est = " << answer[armId] / cntMABSample[armId] << endl;
                    }
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
				fout << ugapeAnswer[i].first << "\t" << ugapeAnswer[i].second << "\t" << part1Map[ugapeAnswer[i].first]
					<< "\t" << part2Map[ugapeAnswer[i].first] // cntMABSample[ugapeAnswer[i].first]
					<< "\t" << part3Map[ugapeAnswer[i].first] // cntMABSample[ugapeAnswer[i].first]
					<< "\t" << part4Map[ugapeAnswer[i].first] / cntMABSample[ugapeAnswer[i].first] << endl;
            }
            cout << "ugape_timer.timeCost = " << ugape_timer.timeCost << " , sortTime = " << sortTime << " , sampleTime = " << sampleTime << " , adaptPushTime = " << adaptPushTime << endl;
            //cout << "adaptDeterQueryTime = " << adaptDeterQueryTime << " , adaptRandQueryTime = " << adaptRandQueryTime << endl;
            return ugapeAnswer;
        }

        // second
        double eps_min = 1e-4;
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
		part1Map.clear(); part2Map.clear(); part3Map.clear(); part4Map.clear(); cntMap.clear();

		Timer adapt_push_timer;
		adapt_push_timer.start();

        adaptResidueMap[queryNode][0][queryNode] = 1; // initialize
        budgetMap[queryNode] = 1;
		forwardPush(queryNode, max(budgetMap[queryNode], 10000), queryNode);
        for (int i = 0; i < curIdx; i++) {
            int armId = ansIdx[i];
            budgetMap[armId] = 1;
            adaptResidueMap[armId][0][armId] = 1; // initialize
			forwardPush(armId, max(budgetMap[armId], 10000), queryNode);
        }
		adapt_push_timer.end();
		adaptPushTime += adapt_push_timer.timeCost;

        int round = 0;
        int batch = 10;
        int batchSample = 10000; // max(10000, (int)(candList.size()));
        cout << "batchSample = " << batchSample << endl;

		for (int i = 0; i < curIdx; i++) {
			int armId = ansIdx[i];
			sampleOneArm(queryNode, armId, 1000);
		}

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
                accum_rank += sameInNbrMap2[g.getInSize(armId)][armId].size();
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
                accum_rank += sameInNbrMap2[g.getInSize(armId)][armId].size();
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
            sort_timer.end();
            sortTime += sort_timer.timeCost;

            if (kthB <= eps_min) {
                vector <pair<int, double> > ugapeAnswer;
                for (int i = 0; i < vecLB.size(); i++) {
                    int armId = vecLB[i].first;
                    vector<int>& sameNbrVec = sameInNbrMap2[g.getInSize(armId)][armId];
                    for (int j = 0; j < sameNbrVec.size(); j++) {
                        int tempNode = sameNbrVec[j];
                        ugapeAnswer.push_back(pair<int, double>(tempNode, answer[armId] / cntMABSample[armId]));
                        cout << "armId = " << tempNode << " , cntMABSample[armId] = " << cntMABSample[armId] 
							 << " , est = " << answer[armId] / cntMABSample[armId] 
							 << " , cb = " << ugapeCB(curIdx, answer[armId], ans2[armId], cntMABSample[armId], 1.0e-4, round) << endl;
                    }
                }
                // TEST
                for (int i = 0; i < vecUB.size(); i++) {
                    int armId = vecUB[i].first;
                    vector<int>& sameNbrVec = sameInNbrMap2[g.getInSize(armId)][armId];
                    for (int j = 0; j < sameNbrVec.size(); j++) {
                        int tempNode = sameNbrVec[j];
                        ugapeAnswer.push_back(pair<int, double>(tempNode, answer[armId] / cntMABSample[armId]));
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
                    fout << ugapeAnswer[i].first << "\t" << ugapeAnswer[i].second 
                       // << "\t" << part1Map[ugapeAnswer[i].first] 
						//<< "\t" << part2Map[ugapeAnswer[i].first] / cntMap[ugapeAnswer[i].first]
						//<< "\t" << part3Map[ugapeAnswer[i].first] / cntMap[ugapeAnswer[i].first]
						//<< "\t" << part4Map[ugapeAnswer[i].first] / cntMABSample[ugapeAnswer[i].first] 
                        << "\t" << ugapeCB(curIdx, answer[ugapeAnswer[i].first], ans2[ugapeAnswer[i].first], cntMABSample[ugapeAnswer[i].first], 1.0e-4, round)
                        << "\t" << ugapeAnswer[i].second - gtMap[ugapeAnswer[i].first] << "\t" << cntMABSample[ugapeAnswer[i].first] << endl;
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

	void forwardPush(int u, int budget, int qNode) {
		cout << "forwardPush : node = " << u << " , prev_budget = " << budgetMap[u] << " , budget = " << budget << " , qNode = " << qNode << endl;
		double cur_rmax = 1.0 / budget; // tunable
		if (u == qNode)
			cur_rmax = 1.0e-6;

		adaptReserveMap[u].clear();
		adaptResidueMap[u].clear(); 
        ppr[u] = unordered_map<int, double>(); // clear
		int numResidueItem = 0;
		unordered_map<int, unordered_map<int, double> > answer;

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
					if (!isInArray[newNode] && newResidue[newNode] > cur_rmax * g.getInSize(newNode)) {
						isInArray[newNode] = true;
						new_candidate_set.push_back(newNode);
					}
				}
				residue[tempNode] = 0;
			}
			for (int j = 0; j < candidate_set.size(); j++) {
				int tempNode = candidate_set[j];
				if (newReserve[tempNode] > 0) {
					if (adaptReserveMap[u][tempLevel][tempNode] == 0)
						adaptReserveMap[u][tempLevel][tempNode] = newReserve[tempNode];
					else
						adaptReserveMap[u][tempLevel][tempNode] += newReserve[tempNode];

                    if(ppr[u].find(tempNode) == ppr[u].end())
                        ppr[u][tempNode] = newReserve[tempNode];
                    else
                        ppr[u][tempNode] += newReserve[tempNode];
				}
				newReserve[tempNode] = 0;
			}
			for (int j = 0; j < residue_count; j++) {
				if (g.getInSize(U[0][j]) > 0 && newResidue[U[0][j]] / g.getInSize(U[0][j]) <= cur_rmax) {
					rsum += newResidue[U[0][j]];
					aliasP.push_back(pair<pair<int, int>, double>(pair<int, int>(tempLevel + 1, U[0][j]), newResidue[U[0][j]]));
					adaptResidueMap[u][tempLevel + 1][U[0][j]] = newResidue[U[0][j]]; //
					numResidueItem++;
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
		fora_timer.end();
		cout << "forward push : " << fora_timer.timeCost << "s" << endl;
		cout << "rsum = " << rsum << endl;

		// update budget
		budgetMap[u] = budget;
		if (nodeToAliasPos.find(u) == nodeToAliasPos.end()) {
			int pos = aliasMap.size();
			nodeToAliasPos[u] = pos;
			aliasMap.push_back(Alias(aliasP)); // at position pos
		}
		else
			aliasMap[nodeToAliasPos[u]] = Alias(aliasP);

		rsumMap[u] = rsum;
		cout << "rsumMap[node] = " << rsumMap[u] << " , cur_rmax = " << cur_rmax << endl;
		cout << "numResidueItem = " << numResidueItem << endl;

		// pre-compute some parts of the SimRank value
		if (u != qNode) {
            Timer deter_timer, rand_timer;
            deter_timer.start();
			precompute_part1(u, qNode);
            deter_timer.end();
            adaptDeterQueryTime += deter_timer.timeCost;
            rand_timer.start();
			precompute_part2(u, qNode, part2Map, u);
			precompute_part2(qNode, u, part3Map, u);
            rand_timer.end();
            adaptRandQueryTime += rand_timer.timeCost;
		}
		//else { // query node
			// for all candidates, precompute
		//	for (int i = 0; i < curIdx; i++) {
		//		int cnode = ansIdx[i];
		//		precompute_part1(cnode, qNode);
		//		precompute_part2(cnode, qNode, part2Map, cnode);
		//		precompute_part2(qNode, cnode, part3Map, cnode);
		//	}
		//}
	}

    
};



#endif
