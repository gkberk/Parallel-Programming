/**
 * Student Name: Sabri Gökberk Yılmaz
 * Student Number: 2017400144
 * Compile Status: Compiling
 * Program Status: Working
 *
 * Note1: The input file is assumed to be correct, and no precaution is taken
 * for a wrong input file since it has nothing to do with the aim of this project
 *
 * Note2: In the Manhattan Distance and diff functions, I used std::abs(). It seems
 * some systems do not allow it to return double values. However
 * my Ubuntu 20.04 with g++ 9.3.0 compiler don't seem to have a trouble with it
 * and my outputs are all correct.
 * Just in case my code has a few output differences on your system
 * like some other people who used std::abs(), I kindly ask for a re-evaluation
 * considering this fact.
 */
#include <iostream>
#include <mpi.h>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <sstream>
#include <string>
#include <cstdlib>
#include <bits/stdc++.h>
#include <set>
#include <iterator>
#include <unistd.h>
using namespace std;

/**
 *
 * @param data1 one of the operand lines
 * @param data2 the other one of the operand lines
 * @return the sum of the absolute differences of two lines' all attributes.
 */

double manhattanDistance(vector<double> &data1, vector<double> &data2){
    double distance=0;
    for(int i=0; i<data1.size()-1; i++){
        distance+=abs(data1[i]-data2[i]);
    }
    return distance;
}

/**
 * implementation of diff function given in the description
 *
 * @param atr current attribute the function will work on
 * @param i1 the operand line's index
 * @param i2 the closest hit/miss' line index
 * @param data input lines of that processor
 * @param weights W vector of the processor
 * @param rank rank of the processor
 * @return abs(operand line's Ath attribute - closest hit/miss' Ath attribute)
 * divided by max value of the attribute - min value of the attribute in the processor
 */
double diff(int atr, int i1, int i2, vector<vector<double>> &data, vector<double> &weights, int rank){

    double result;
    double max=0;
    double min=9999;
    for(int t=0; t<data.size(); t++){//loop to find max and min of atr
        if(data[t][atr]<min){
            min=data[t][atr];
        }
        if(max<data[t][atr]){
            max=data[t][atr];
        }
    }
    result=abs(data[i1][atr]-data[i2][atr])/(max-min);
    return result;
}

/**
 * implementation of the relief algorithm
 * @param data input lines of the processor
 * @param weights W vector of the processor
 * @param iterations no of iterations, M
 * @param a no of attributes
 * @param rank rank of the processor
 */
void relief(vector<vector<double>> &data, vector<double> &weights, int iterations, int a, int rank){

    vector<double> distancetemp;
    double closestHitDistance;
    double closestMissDistance;
    int closestHitIndex;
    int closestMissIndex;
    for(int itno=0; itno<iterations; itno++){ //iterations of the relief algorithm

        closestMissIndex=-3;
        closestMissIndex=-3;
        closestHitDistance=9999;
        closestMissDistance=9999;
        distancetemp.clear();

        for(int j=0; j<data.size(); j++){//calculating distances
            if(itno==j)
                distancetemp.push_back(9999);
            else {
                distancetemp.push_back(manhattanDistance(data[itno], data[j]));//pushing distances
            }
        }


        for(int k=0; k<data.size(); k++){//finding min distances
            if((int)data[k][a]!=data[itno][a]){//it's a miss
                if(closestMissDistance>distancetemp[k]){
                    closestMissDistance=distancetemp[k];
                    closestMissIndex=k;
                }
            }
            else{//it's a hit
                if(closestHitDistance>distancetemp[k]){
                    closestHitDistance=distancetemp[k];
                    closestHitIndex=k;
                }
            }
        }

        for(int atr=0; atr<a; atr++){ //modifying W vector
            weights[atr]=weights[atr] - diff(atr,itno,closestHitIndex,data,weights, rank)/iterations
                    +diff(atr,itno,closestMissIndex,data,weights, rank)/iterations;
        }
    }
}

/**
 * the main function
 * @param argc no of argumnents
 * @param argv argument array
 * @return 0, if exits gracefully
 */
int main(int argc, char* argv[]){

    int rank;
    int size;

    ///here starts the parallel programming.
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(rank==0){//master
        fstream file;
        string fname=argv[1];
        file.open(fname);
        int p, n, a, m, t;
        string s;
        if(file.is_open()) //reading the file's first 2 lines
            file>>p>>n>>a>>m>>t;
        getline(file, s);

        vector<vector<double>> lines;//the inputs are in this vector
        vector<double> temp;
        double d;

        for(int i=0; i<n; i++){//reading contents of the line
            getline(file, s);
            stringstream scan(s);

            for(int i=0; i<a+1; i++){
                scan>>d;
                temp.push_back(d);
            }
            lines.push_back(temp);
            temp.clear();
        }

        int info[]={p, n, a, m, t};
        for(int i=1; i<p; i++){
            ///here the info array is sent to all slaves which contains p,n,a,m,t values given
            /// at the first two lines of the input file with tag
            MPI_Send(&info, 5, MPI_INT, i, 0, MPI_COMM_WORLD); //data, size, type, receiver, tag
            for(int j=0; j<n/(p-1); j++){
                ///here the data is split between slave processors with tag as the sent line's index according to receiver
                MPI_Send(&lines[(i-1)*n/(p-1)+j][0], a+1, MPI_DOUBLE, i, j+1, MPI_COMM_WORLD);
            }
        }


        set<int> final;

        for(int r=1; r<p; r++){
            int temp[t];
            ///here the outputs of the slaves are received by an array then inserted into a set
            MPI_Recv(&temp[0], t, MPI_INT, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for(int e=0; e<t; e++){
                final.insert(temp[e]);
            }
        }

        vector<int> final2(final.size());
        copy(final.begin(), final.end(), final2.begin());//converting the set to vector

        int* final3=&final2[0];//converting the vector to array
        sort(final3, final3+final2.size());
        usleep(500);//0.5 second of sleep so the master prints as the last.
        cout<<"Master P0 : ";
        cout<<final3[0];
        for(int i=1; i<final2.size(); i++){//printing the master's output
            cout<<" "<<final3[i];
        }
        cout<<endl;
    }

    else{/// rank>0. Hence, slave
        int info[5];
        ///here the p,n,a,m,t numbers are received by an array
        MPI_Recv(&info, 5, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        int p,n,a,m,t;
        p=info[0];
        n=info[1];
        a=info[2];
        m=info[3];
        t=info[4];

        vector<vector<double>> lines;//the lines belonging to the slave
        vector<double> line;
        for(int j=0; j<n/(p-1); j++) {
            line.resize(a+1);
            ///here the slave receives its input and stores them
            MPI_Recv(&line[0], a + 1, MPI_DOUBLE, 0, j+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            lines.push_back(line);
            for(int i=0; i<line.size(); i++){
            }
            line.clear();
        }
        vector<double> weights;///W vector
        weights.resize(a);
        relief(lines, weights, m, a, rank);///using relief algorithm
        int output[t];
        int curmax=0;
        int maxind;

        vector<double> w;
        for(int q=0; q<weights.size(); q++){//copy weights into w
            w.push_back(weights[q]);
        }

        double max=0;
        int maxin;
        int result[t];
        for(int z=0; z<t; z++){//detect the results
            for(int v=0; v<a; v++){
                if(max<w[v]){
                    maxin=v;
                    max=w[v];
                }
            }
            result[z]=maxin;
            max=0;
            w[maxin]=-1;
        }
        int size=sizeof(result)/sizeof(result[0]);
        sort(result, result+size);
        //print the results
        cout<<"Slave P"<<rank<<" : ";
        cout<<result[0];
        for(int y=1; y<t; y++){
            cout<<" "<<result[y];
        }
        cout<<endl;
        ///here the slave sends its result array to the master
        MPI_Send(&result, t, MPI_INT, 0, 0, MPI_COMM_WORLD); //data, size, type, receiver, tag
    }
    ///here the parallel programming ends
    MPI_Finalize();
    exit(0);
}