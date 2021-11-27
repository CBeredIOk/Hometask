#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <stdlib.h>
using namespace std;

int main() {
    string line;
    ifstream file("in.txt");

    double x0;
    double y0;
    file >> x0 >> y0;

    double A = x0, B = y0;

    double x_left = 0, y_left = 0, x_right = 0, y_right = 0;
    double x_left_m = 0, y_left_m = 0, x_right_m = 0, y_right_m = 0;

    double h_right_0 = fabs((B * x_right - A * y_right) / sqrt(A * A + B * B));
    double h_left_0 = fabs((B * x_left - A * y_left) / sqrt(A * A + B * B));

    double x;
    double y;
    double h_left_1;
    double h_right_1;
    double y_right_m_0;
    double x_right_m_0;

    while (file >> x >> y) {
        if ((A >= 0 and B > 0) or (A > 0 and B <= 0)) {
            if (B * x - A * y < 0) {
                x_left = x;
                y_left = y;

//                cout << "\nx_l = " << x;
//                cout << "\ny_l = " << y;

                h_left_1 = fabs((B * x_left - A * y_left) / sqrt(A * A + B * B));
                h_left_1 = round(h_left_1 * 10000000000.0) / 10000000000.0;
                if (h_left_1 >= h_left_0) {
                    h_left_0 = h_left_1;
                    x_left_m = x;
                    y_left_m = y;
                }
            } else {
                x_right = x;
                y_right = y;

//                cout << "\nx_r = " << x;
//                cout << "\ny_r = " << y;

                h_right_1 = fabs((B * x_right - A * y_right) / sqrt(A * A + B * B));
                h_right_1 = round(h_right_1 * 10000000000.0) / 10000000000.0;
                if (h_right_1 >= h_right_0) {
                    h_right_0 = h_right_1;
                    x_right_m = x;
                    y_right_m = y;
                }
                if (h_right_1 == 0){
                    x_right_m_0 = x;
                    y_right_m_0 = y;
                }
            }
        }
        if ((A <= 0 and B < 0) or (A < 0 and B >= 0)) {
            if (B * x - A * y < 0) {
                x_left = x;
                y_left = y;

//                cout << "\nx_l = " << x;
//                cout << "\ny_l = " << y;

                h_left_1 = fabs((B * x_left - A * y_left) / sqrt(A * A + B * B));
                h_left_1 = round(h_left_1 * 10000000000.0) / 10000000000.0;
                if (h_left_1 >= h_left_0) {
                    h_left_0 = h_left_1;
                    x_left_m = x;
                    y_left_m = y;
                }

            } else {
                x_right = x;
                y_right = y;

//                cout << "\nx_r = " << x;
//                cout << "\ny_r = " << y;

                h_right_1 = fabs((B * x_right - A * y_right) / sqrt(A * A + B * B));
                h_right_1 = round(h_right_1 * 10000000000.0) / 10000000000.0;
                if (h_right_1 >= h_right_0) {
                    h_right_0 = h_right_1;
                    x_right_m = x;
                    y_right_m = y;
                }
                if (h_right_1 == 0){
                    x_right_m_0 = x;
                    y_right_m_0 = y;
                }
            }
        }
    }
    file.close();
    if (h_right_1 == 0){
        cout << "\nLeftmost: " << x_left_m << " " << y_left_m;
        cout << "\nRightmost: " << x_right_m_0 << " " << y_right_m_0;
    }
    else {
        cout << "\nLeftmost: " << x_left_m << " " << y_left_m;
        cout << "\nRightmost: " << x_right_m << " " << y_right_m;
    }
}

//int main() {
//    string line;
//    ifstream file("in.txt");
//
//    int i = 0;
//    double x0;
//    double y0;
//    file >> x0 >> y0;
//
//    double A = x0, B = y0;
//
//    double x_left = 0, y_left = 0, x_right = 0, y_right = 0;
//    double x_left_m = 0, y_left_m = 0, x_right_m = 0, y_right_m = 0;
//
//    double h_right_0 = fabs((B * x_right - A * y_right) / sqrt(A * A + B * B));
//    double h_left_0 = fabs((B * x_left - A * y_left) / sqrt(A * A + B * B));
//
//    double x;
//    double y;
//    while (file >> x >> y) {
//        if ((A >= 0 and B > 0) or (A > 0 and B <= 0)) {
//            if (B * x - A * y < 0) {
//                x_left = x;
//                y_left = y;
//                //cout << "\nx_l = " << x;
//                //cout << "\ny_l = " << y;
//            } else {
//                x_right = x;
//                y_right = y;
//                //cout << "\nx_r = " << x;
//                //cout << "\ny_r = " << y;
//            }
//        }
//        if ((A <= 0 and B < 0) or (A < 0 and B >= 0)) {
//            if (B * x - A * y < 0) {
//                x_left = x;
//                y_left = y;
//                //cout << "\nx_l = " << x;
//                //cout << "\ny_l = " << y;
//            } else {
//                x_right = x;
//                y_right = y;
//                //cout << "\nx_r = " << x;
//                //cout << "\ny_r = " << y;
//            }
//        }
//        double h_right_1;
//        double h_left_1;
//        h_right_1 = fabs((B * x_right - A * y_right) / sqrt(A * A + B * B));
//        h_right_1 = round(h_right_1 * 10000000000.0) / 10000000000.0;
//        h_left_1 = fabs((B * x_left - A * y_left) / sqrt(A * A + B * B));
//        h_left_1 = round(h_left_1 * 10000000000.0) / 10000000000.0;
//
//        if (h_right_1 > h_right_0) {
//            h_right_0 = h_right_1;
//            x_right_m = x;
//            y_right_m = y;
//        }
//        if (h_left_1 > h_left_0) {
//            h_left_0 = h_left_1;
//            x_left_m = x;
//            y_left_m = y;
//        }
//        if (h_right_1 == 0 and i == 0){
//            x_right_m = x;
//            y_right_m = y;
//            i++;
//        }
//        //cout << "\ni = " << i;
//    }
//
//    file.close();
//    cout << "\nLeftmost: " << x_left_m << " " << y_left_m;
//    cout << "\nRightmost: " << x_right_m << " " << y_right_m;
//
//}