#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
using namespace std;

double bias(int n, vector <double> X) {
    double x;
    x = 0;
    for(int i = 1; i <= n; i++){
        x = 2*X[i-1] - x;
    }
    return x;
}

double function(int n, double h, double x, double Vx, double Vy, vector<double>& X) {
    double g = 9.81;
    double y = h + pow(-1, n) * (x - bias(n, X)) * Vy/Vx - g / (2 * pow(Vx,2)) * pow(x - bias(n, X), 2);
    return y;
}

int main(int argc, char** argv){
    if(argc == 2) {
        string input;
        ifstream file("in.txt");
        double h;
        file >> h;

        double Vx, Vy;
        file >> Vx >> Vy;

        double xi, yi;
        vector<double> X_b;
        vector<double> Y_b;

        while (file >> xi >> yi) {
            X_b.push_back(xi);
            Y_b.push_back(yi);
        }

        if (X_b.size() == 0) {
            cout << 0;
        }

        // k - number of the last barrier the ball hit, n - number of collision
        // r - Checking the flight after a collision flying to the right
        int k = 0, n = 0, r = 0;

        // Way = 0 = " from Right "
        // Way = 1 = " from Left "
        int Way = 0;

        vector<double> Collision_X;

        // if the collision occurred against the first barrier
        double y = function(n, h, X_b[0], Vx, Vy, X_b);

        if (y < Y_b[0]) {
            cout << "\n" << 0;
            return 0;
        }

        int Above = 0;

        for (int i = 0; i < X_b.size(); i++) {
            y = function(n, h, X_b[i], Vx, Vy, X_b);
            if (y > Y_b[i]) {
                Above += 1;
            } else {
                Way = 1;
                k = i;
                n = 1;

                Collision_X.push_back(X_b[i]);
                break;
            }
        }

        if (Above == X_b.size()) {
            cout << X_b.size();
            return 0;
        }

        if(y > 0){
            while (true) {
                if (Way == 1) {
                    for (int i = k - 1; i >= 0; i--) {
                        y = function(n, h, X_b[0], Vx, Vy, Collision_X);
                        if (y <= Y_b[i]) {
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
                if (Way == 0) {
                    for (int i = k + 1; i < X_b.size(); i++) {
                        y = function(n, h, X_b[0], Vx, Vy, Collision_X);
                        if (y <= Y_b[i]) {
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
        }

        if (Way == 1) {
            if (r == k) {
                cout << "\n" << 0;
            } else {
                cout << "\n" << k;
            }
        }

        if (Way == 0) {
            if (r == k) {
                cout << "\n" << 0;
            } else {
                cout << "\n" << k;
            }
        }
        return 0;
    }
    else{
        cerr << " Heh " << endl;
        return 1;
    }
}

//#include <iostream>
//#include <vector>
//#include <cmath>
//#include <string>
//#include <fstream>
//using namespace std;
//
//double bias(int n, vector <double> X) {
//    double x;
//    x = 0;
//    for(int i = 1; i <= n; i++){
//        x = 2*X[i-1] - x;
//    }
//    return x;
//}
//
////double recurs(int n, vector<double> X) {
////    if (n == 0) {
////        return 0;
////    }
////    else {
////        return 2 * X[n-1] - recurs(n - 1, X);
////    }
////}
//
//double function(int n, double h, double x, double Vx, double Vy, vector<double>& X) {
//    double g = 9.81;
//    double y = h + pow(-1, n) * (x - bias(n, X)) * Vy/Vx - g / (2 * pow(Vx,2)) * pow(x - bias(n, X), 2);
//    return y;
//}
//
//int main() {
//    string input;
//    ifstream file("in.txt");
//    double h;
//    file >> h;
//
//    double Vx, Vy;
//    file >> Vx >> Vy;
//
//    double xi, yi;
//    vector<double> X_b;
//    vector<double> Y_b;
//
//    while (file >> xi >> yi) {
//        X_b.push_back(xi);
//        Y_b.push_back(yi);
//    }
//
//    if (X_b.size() == 0) {
//        cout << 0;
//    }
//
//    // k - number of the last barrier the ball hit, n - number of collision
//    // r - Checking the flight after a collision flying to the right
//    int k = 0, n = 0, r = 0;
//
//    // Way = 0 = " from Right "
//    // Way = 1 = " from Left "
//    int Way = 0;
//
//    vector<double> Collision_X;
//
//    // if the collision occurred against the first barrier
//    double y = function(n, h, X_b[0], Vx, Vy, X_b);
//
//    if (y < Y_b[0]) {
//        cout << "\n" << 0;
//        return 0;
//    }
//
//    int Above = 0;
//
//    for (int i = 0; i < X_b.size(); i++) {
//        y = function(n, h, X_b[i], Vx, Vy, X_b);
//        if (y > Y_b[i]) {
//            Above += 1;
//        } else {
//            Way = 1;
//            k = i;
//            n = 1;
//
//            Collision_X.push_back(X_b[i]);
//            break;
//        }
//    }
//
//    if (Above == X_b.size()) {
//        cout << X_b.size();
//        return 0;
//    }
//
//    if(y > 0){
//        while (true) {
//            if (Way == 1) {
//                for (int i = k - 1; i >= 0; i--) {
//                    y = function(n, h, X_b[0], Vx, Vy, Collision_X);
//                    if (y <= Y_b[i]) {
//                        Way = 0;
//                        k = i;
//                        n++;
//
//                        Collision_X.push_back(X_b[i]);
//                        break;
//                    }
//                    r++;
//                }
//                double yb = y;
//                if (yb < 0) break;
//
//                if (r == k) {
//                    break;
//                }
//            }
//            if (Way == 0) {
//                for (int i = k + 1; i < X_b.size(); i++) {
//                    y = function(n, h, X_b[0], Vx, Vy, Collision_X);
//                    if (y <= Y_b[i]) {
//                        Way = 1;
//                        k = i;
//                        n++;
//
//                        Collision_X.push_back(X_b[i]);
//                        break;
//                    }
//                }
//                double yb = y;
//                if (yb < 0) break;
//            }
//        }
//    }
//
//    if (Way == 1) {
//        if (r == k) {
//            cout << "\n" << 0;
//        } else {
//            cout << "\n" << k;
//        }
//    }
//
//    if (Way == 0) {
//        if (r == k) {
//            cout << "\n" << 0;
//        } else {
//            cout << "\n" << k;
//        }
//    }
//    return 0;
//}
