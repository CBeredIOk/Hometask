#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
using namespace std;

double recurs( int n, vector <double> X) {
    if (n == 0) {
        return 0;
    }
    else {
        return 2 * X[n-1] - recurs(n-1, X);
    }
}
double function(int n, double h, double x, double Vx, double Vy, vector<double>& X) {
    double g = 9.81;
    double f = h + pow(-1, n) * (x - recurs(n, X)) * Vy/Vx - g / (2 * pow(Vx,2)) * pow(x - recurs(n, X), 2);
    return f;
}

int main(int argc, char** argv){
    if(argc == 2) {
        double a;
        int nn;
        vector<double> v;
        fstream F;
        F.open(argv[1]);
        if (F) {
            while (!F.eof()) {
                F >> a;
                v.push_back(a);
                nn++;
            }
            F.close();
        }
        else{
            cout << "Fail doesn't exist" << endl;
        }

        vector<double> Y_b;
        vector<double> X_b;
        for (int i = 3; i < nn; i += 2) {
            X_b.push_back(v[i]);
            Y_b.push_back(v[i + 1]);
        }

        double h = v[0];
        double Vx = v[1];
        double Vy = v[2];
        int k = 0, n = 0, r = 0;

        //string Way = " from Right ";
        int Way = 0;

        vector<double> Collision_X;
        double y = function(n, h, X_b[0], Vx, Vy, X_b);
        if (y < Y_b[0]) {
            cout << 0;
            return 0;
        }
        int High = 0;

        for (int i = 0; i < X_b.size(); i++) {
            y = function(n, h, X_b[0], Vx, Vy, X_b);
            if (y > Y_b[i]) {
                High += 1;
            }
        }

        if (High == X_b.size()) {
            cout << X_b.size();
            return 0;
        }

        for (int i = 0; i < X_b.size(); i++) {
            y = function(n, h, X_b[0], Vx, Vy, X_b);
            if (y <= Y_b[i]) {
                //Way = " from Left ";
                Way = 1;
                k = i;
                n = 1;

                Collision_X.push_back(X_b[i]);
                break;
            }
        }

        while (true) {
            //if (Way == " from Left ") {
            if (Way == 1) {
                for (int i = k - 1; i >= 0; i--) {
                    y = function(n, h, X_b[0], Vx, Vy, Collision_X);
                    if (y <= Y_b[i]) {
                        //Way = " from Right ";
                        Way = 0;
                        k = i;
                        n++;
                        Collision_X.push_back(X_b[i]);
                        break;
                    }
                    r++;
                }
                double yb = y;
                if (yb < 0) break;

                if (r == k) {
                    break;
                }
            }
            //if (Way == " from Right ") {
            if (Way == 0) {
                for (int i = k + 1; i < X_b.size(); i++) {
                    y = function(n, h, X_b[0], Vx, Vy, Collision_X);
                    if (y <= Y_b[i]) {
                        //Way = " from Left ";
                        Way = 1;
                        k = i;
                        n++;
                        Collision_X.push_back(X_b[i]);
                        break;
                    }
                }
                double yb = y;
                if (yb < 0) break;
            }
        }

        //if (Way == " from Left ") {
        if (Way == 1) {
            if (r == k) {
                cout << "<" << 0 << ">";
            } else {
                cout << "<" << k << ">";
            }
        }

        //if (Way == " from Right ") {
        if (Way == 0) {
            cout << "<" << k+1 << ">";
        }
        return 0;
    }
    else{
        cerr << " Heh " << endl;
        return 1;
    }
}
