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
//int main(){
//    string input;
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

int main(){
    int nn;
    string input;
    ifstream file("in.txt");

    vector<double> v;

    double vi;

    while (file >> vi) {
        v.push_back(vi);
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
        cout << X_b.size()-1;
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
            cout << 0 ;
        } else {
            cout << k;
        }
    }

    //if (Way == " from Right ") {
    if (Way == 0) {
        cout << k+1;
    }
    return 0;
}
























