#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
using namespace std;

double bias( int n, vector <double> X) {
    double x;
    x = 0;
    for(int i = 1; i <= n; i++){
        x = 2*X[i-1] - x;
    }
    return x;
}

double function(int n, double h, double x, double Vx, double Vy, vector<double>& X) {
    double g = 9.81;
    double f = h + pow(-1, n) * (x - bias(n, X)) * Vy/Vx - g / (2 * pow(Vx,2)) * pow(x - bias(n, X), 2);
    return f;
}

int main(int argc, char** argv){
    if(argc == 2) {
        string input;
        ifstream file(argv[1]);
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

        if (X_b.size() == 0){
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
            y = function(n, h, X_b[0], Vx, Vy, X_b);
            if (y > Y_b[i]) {
                Above += 1;
            }
            else {
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

        if (Way == 1) {
            if (r == k) {
                cout << "\n" << 0 ;
            } else {
                cout << "\n" << k;
            }
        }

        if (Way == 0) {
            cout << "\n" << k+1;
        }
        return 0;
    }
    else{
        cerr << " Heh " << endl;
        return 1;
    }
}

//double bias( int n, vector <double> X) {
//    double x;
//    x = 0;
//    for(int i = 1; i <= n; i++){
//        x = 2*X[i-1] - x;
//    }
//    return x;
//}
//
//double function(int n, double h, double x, double Vx, double Vy, vector<double>& X) {
//    double g = 9.81;
//    double f = h + pow(-1, n) * (x - bias(n, X)) * Vy/Vx - g / (2 * pow(Vx,2)) * pow(x - bias(n, X), 2);
//    return f;
//}
//
//int main(){
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
//    if (X_b.size() == 0){
//        cout << 0;
//    }
//
//    cout << "\n" << X_b[0];
//    cout << "\n" << X_b[1];
//    cout << "\n" << X_b[2];
//    cout << "\n" << X_b[3];
//    cout << "\n" << " ";
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
//    if (y < Y_b[0]) {
//        cout << "\n" << 0;
//        return 0;
//    }
//
//    int Above = 0;
//
//    for (int i = 0; i < X_b.size(); i++) {
//        y = function(n, h, X_b[0], Vx, Vy, X_b);
//        if (y > Y_b[i]) {
//            Above += 1;
//        }
//        else {
//            Way = 1;
//            k = i;
//            n = 1;
//            cout << "\n" << "y: " << y;
//            cout << "\n" << "barrier_X: " << X_b[i];
//            cout << "\n" << "barrier_Y: " << Y_b[i];
//            cout << "\n" << k;
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
//    while (true) {
//        if (Way == 1) {
//            for (int i = k - 1; i >= 0; i--) {
//                y = function(n, h, X_b[0], Vx, Vy, Collision_X);
//                if (y <= Y_b[i]) {
//                    Way = 0;
//                    k = i;
//                    cout << "\n" << "y: " << y;
//                    cout << "\n" << "barrier_X: " << X_b[i];
//                    cout << "\n" << "barrier_Y: " << Y_b[i];
//                    cout << "\n" << k;
//                    n++;
//                    Collision_X.push_back(X_b[i]);
//                    break;
//                }
//                r++;
//            }
//            double yb = y;
//            if (yb < 0) break;
//
//            if (r == k) {
//                break;
//            }
//        }
//        if (Way == 0) {
//            for (int i = k + 1; i < X_b.size(); i++) {
//                y = function(n, h, X_b[0], Vx, Vy, Collision_X);
//                if (y <= Y_b[i]) {
//                    Way = 1;
//                    k = i;
//                    cout << "\n" << "y: " << y;
//                    cout << "\n" << "barrier_X: " << X_b[i];
//                    cout << "\n" << "barrier_Y: " << Y_b[i];
//                    cout << "\n" << k;
//                    n++;
//                    Collision_X.push_back(X_b[i]);
//                    break;
//                }
//            }
//            double yb = y;
//            if (yb < 0) break;
//        }
//    }
//
//    if (Way == 1) {
//        if (r == k) {
//            cout << "\n" << 0 ;
//        } else {
//            cout << "\n" << k;
//        }
//    }
//
//    if (Way == 0) {
//        cout << "\n" << k+1;
//    }
//    return 0;
//}


//#include <iostream>
//#include <fstream>
//#include <vector>
//#include <cmath>
//#include <string>
//
//using namespace std;
//
//template <typename T>
//int sgn(T val) {
//    return (T(0) < val) - (val < T(0));
//}
//
//const double g = 9.81;
//
//
//
//int main()
//{
//    string input;
//
//    ifstream cin("in.txt");
//
//    double h0;
//    cin >> h0;
//    double vx0, vy0;
//    cin >> vx0 >> vy0;
//    double xi, hi;
//    vector<double> xis;
//    vector<double> his;
//
//    double y = h0, t = 0, x = 0.0, vx = vx0, vy = vy0;
//    int ans = 0;
//
//
//    xis.push_back(0.0);
//    his.push_back(-1.);
//
//
//    while (cin >> xi >> hi) {
//        xis.push_back(xi);
//        his.push_back(hi);
//    }
//
//    xis.push_back(xis.back()*2+10.0);
//    his.push_back(-1.);
//
//    while (true) {
//        int dir = sgn(vx);
//        int next = (dir > 0) ? ans + dir : ans;
//
//        if (his[next] < 0 ) {
//            std::cout << ans;
//            return 0;
//        }
//
//
//        if (ans == 0 && dir < 0) {
//            std::cout << 0;
//            return 0;
//        }
//
//        double dt = (xis[next] - x) / vx;
//        y = y + vy * dt - g * dt * dt / 2;
//
//        if (y <= 0.0) {
//            std::cout << ans;
//            return 0;
//        }
//        else if (y > his[next]) {
//
//            x = xis[next];
//
//            ans += dir;
//        }
//        else if (y <= his[next]) {
//
//            vx = -vx;
//            x = xis[next];
//
//        }
//        vy = vy - g * dt;
//        t += dt;
//
//    }
//    return 0;
//}























