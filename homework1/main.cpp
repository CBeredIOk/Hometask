#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

std::vector<double> bubbleSort(std::vector<double>& h) {
    double temp;
    for (int i = 0; i < h.size(); i++) {
        for (int j = 0; j < h.size()-1; j++) {
            if (h[j] > h[j + 1]) {
                temp = h[j];
                h[j] = h[j+1];
                h[j + 1] = temp;
            }
        }
    }
    return h;
}

std::vector<double> distance(std::vector<double> &v_x, std::vector<double> &v_y, int n, double A, double B, double C){
    double H;
    std::vector<double> h;
    for (int i = 0; i < v_x.size(); i++) {
        H = fabs((B * v_x[i] + A * v_y[i] + C) / sqrt(A * A + B * B));
        h.push_back(H);
    }
    return h;
}

std::vector<double> divide(std::vector<double>& v, std::string type) {
    int n = v.size();
    std::vector<double> v_x;
    std::vector<double> v_y;
    for (int i = 0; i < n; i++) {
        if (i % 2 == 0) {
            v_x.push_back(v[i]);
        }
        else {
            v_y.push_back(v[i]);
        }
    }
    return (type == "X") ? v_x : v_y;
}

std::vector<double> right_or(std::vector<double>& v_x, std::vector<double>& v_y, double A, double B, std::string type) {
    std::vector<double> v_x_r;
    std::vector<double> v_y_r;
    for(int i = 1; i < v_x.size(); i++){
        if((A and B >= 0) or (A > 0 and B <= 0 )) {
            if (B * v_x[i] - A * v_y[i] <= 0) {
                v_y_r.push_back(v_y[i]);
                v_x_r.push_back(v_x[i]);
            }
        }
        if ((A and B <= 0) or (A <= 0 and B >= 0 )){
            if (B * v_x[i] - A * v_y[i] >= 0) {
                v_y_r.push_back(v_y[i]);
                v_x_r.push_back(v_x[i]);
            }
        }
    }
    return (type == "X") ? v_x_r : v_y_r;
}
std::vector<double> left_or(std::vector<double>& v_x, std::vector<double>& v_y, double A, double B, std::string type) {
    std::vector<double> v_x_l;
    std::vector<double> v_y_l;
    for(int i = 1; i < v_x.size(); i++){
        if((A and B >= 0) or (A >= 0 and B < 0 )) {
            if (B * v_x[i] - A * v_y[i] > 0) {
                v_y_l.push_back(v_y[i]);
                v_x_l.push_back(v_x[i]);
            }
        }
        if ((A and B < 0) or (A <= 0 and B > 0 )){
            if (B * v_x[i] - A * v_y[i] < 0) {
                v_y_l.push_back(v_y[i]);
                v_x_l.push_back(v_x[i]);
            }
        }
    }
    return (type == "X") ? v_x_l : v_y_l;
}

int main() {
    int n = 0;
    float a;
    std::vector<double> v;
    std::vector<double> v_r;
    std::vector<double> v_l;
    std::vector<double> h_l;
    std::vector<double> h_r;

    std::vector<double> h_0;

    std::fstream F;
    F.open("in.txt");
    if (F) {
        while (!F.eof()) {
            F >> a;
            v.push_back(a);
            n++;
        }
        F.close();
    } else std::cout << "Fail doesn't exist" << std::endl;

    double A, B, C;
    double x0 = v[0];
    double y0 = v[1];

    A = x0;
    B = y0;
    C = 0;

    std::vector<double> v_x = divide(v, "X");
    std::vector<double> v_y = divide(v, "Y");
    std::vector<double> v_x_r = right_or(v_x, v_y, A, B, "X");
    std::vector<double> v_y_r = right_or(v_x, v_y, A, B, "Y");
    std::vector<double> v_x_l = left_or(v_x, v_y, A, B, "X");
    std::vector<double> v_y_l = left_or(v_x, v_y, A, B, "Y");

    if(v_x_r.size() == 0){
        v_y_r.push_back(0);
        v_x_r.push_back(0);
    }
    if(v_x_l.size() == 0){
        v_y_l.push_back(0);
        v_x_l.push_back(0);
    }

    h_l = distance(v_x_l, v_y_l, n, A, B, C);
    std::vector<double> h_ll = h_l;
    h_r = distance(v_x_r, v_y_r, n, A, B, C);
    std::vector<double> h_rr = h_r;
    std::vector<double> h_l_b = bubbleSort(h_l);
    std::vector<double> h_r_b = bubbleSort(h_r);

    double max_r;
    int index_r;
    max_r = h_r_b[h_r_b.size()-1];
    for(int i =0; i < h_r.size(); i++) {
        if (max_r == h_rr[i]) {
            index_r = i;
        }
    }

    double max_l;
    int index_l;
    max_l = h_l_b[h_l_b.size()-1];
    for (int i = 0; i < h_l.size(); i++) {
        if (h_ll[i] == max_l) {
            index_l = i;
        }
    }

    std::cout << "\nLeftmost: " <<v_x_l[index_l] << " " << v_y_l[index_l];
    std::cout << "\nRightmost: " <<v_x_r[index_r] << " " << v_y_r[index_r];
    return 0;
}