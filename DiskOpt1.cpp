//Tasks1:
//1 -- Realize method of finite elements (FEM)  |   All right!
//2 -- Realize optimization method              |   Realize the 1D Newton method on each node -> optimum profile, epsilon about delta in mass

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>

//Constants of task:

//Al: E = 70, rho = 2700, nu = 0.34

double boundary_1 = 0.001, boundary_2 = 0.03;      //0.03

double mu = 0.33, nu = mu/*25.5 * 1e9*/ , E = 72 * 1e9, w = 1000, rho = 2850, pi = 3.141592;//w = 1000, 5000, 8000

double sigma_max = 25 * 1e9;

double h_max(double r){ return 1; };
                                        
double h_min(double r){
    if(r <= 0.1)
        return 0;
    if((0.1 <= r) and (r <= 0.3))
        return 0;
    if(r >= 0.225)
        return 0;

    return 0;
}

//Constants for numerical method:

int nodes = 101;  //101 -- is gooood

double begin = 0.1, end = 0.3;
double step  = (end - begin) * 1.0 / (nodes - 1);

std::vector<double> grid = {};   

//Variables for optimization

std::vector<double> h = {}, Dh = {}, u = {}, Du = {}, Q = {}, Sigma = {}, m = {};

//Functions for optimization

std::pair<double, int> maximum(){
    double maxx = h[0];
    int index;
    for(int i = 1; i < grid.size(); i++)
        if(maxx < h[i]){
            maxx = h[i];
            index  = i;
        }
    return {maxx, index};
}

void Init_h(){                              //Initialization profile h
    for(int i = 0; i < grid.size(); i++)
        h.push_back(0.1);
}

//Functions for vectors

std::ostream& operator<<(std::ostream& out, const std::vector<double>& a){
    out << "[ ";
    for(int i = 0; i < a.size(); i++)
        out << a[i] << " ";
    out << "]";
    return out;
}

std::vector<double> operator+(std::vector<double> a, std::vector<double> b){
    std::vector<double> res;
    res.resize(a.size());

    for(int i = 0; i < a.size(); i++)
        res[i] = a[i] + b[i];
    
    return res;
}

std::vector<double> operator-(std::vector<double> a, std::vector<double> b){
    std::vector<double> res;
    res.resize(a.size());

    for(int i = 0; i < a.size(); i++)
        res[i] = a[i] - b[i];

    return res;
}

std::vector<double> operator*(double a, std::vector<double> b){
    std::vector<double> res;
    for(int i = 0; i < b.size(); i++)
        res.push_back(b[i] * a);
    
    return res;
}

std::vector<double> Incremention(int index, double coef){
    std::vector<double> res;
    for(int i = 0; i < h.size(); i++)
        if(i == index)
            res.push_back(coef);
        else
            res.push_back(0);
    return res;
}

//Common functions:

double min(int a, int b){
    if(a < b)
        return a;
    else
        return b;
}

double max(int a, int b){
    if(a > b)
        return a;
    else
        return b;
}

//Functions for FEM:

double B(double r, double k){                           //Basis functions
    if((grid[k] - 1e-9 < r) and (r < grid[k] + 1e-9))
        return 1;
    
    if((grid[k - 1] < r) and (r < grid[k]))
        return (r - grid[k - 1]) / step;

    if((grid[k + 1] > r) and (r > grid[k]))
        return (grid[k + 1] - r) / step;

    return 0;
}

double u_func(double r){
    double sum = 0;

    for(int j = 0; j < nodes; j++)
            sum += u[j] * B(r, j);

    return sum;
}

std::vector<double> Get_Du(){
    std::vector<double> du = {};
    du.push_back((u[0] + u[1]) / step);

    for(int i = 0; i < u.size() - 1; i++)
        du.push_back(0.5 * (u[i + 1] - u[i - 1]) / step);

    du.push_back((u[u.size() - 1] - u[u.size() - 2]) / step);
    return du;
}

double h_func(double r){
    double sum = 0;

    for(int j = 0; j < nodes; j++)
            sum += h[j] * B(r, j);

    return sum;
}

double Dh_func(double r){
    if((grid[grid.size() - 1] - 1e-8 < r) and (r < grid[grid.size() - 1] + 1e-8))
        return (h_func(r) - h_func(r - 1e-8)) / 1e-8;
    if((grid[0] - 1e-8 < r) and (r < grid[0] + 1e-8))
        return (h_func(r + 1e-8) - h_func(r)) / 1e-8;

    return (h_func(r + 1e-8) - h_func(r - 1e-8)) / (2 * 1e-8);
}

double a1(double r){ return r * r * h_func(r); };

double Da1(double r){ return (a1(r + 1e-8) - a1(r - 1e-8)) / (2 * 1e-8); };

double a2(double r){ return r * (h_func(r) + r * Dh_func(r)); };

double a3(double r){ return -h_func(r) + Dh_func(r) *  mu * r; };

double  b(double r){ return -r * r * r * rho * w * w * h_func(r) * (1 - mu * mu) / E; };

double Integral_Of_a1_On_Element(int elem_i, int elem_j, int func_i, int func_j){   //on [elem_i, elem_j] integrate DB(r, func_i) * DB(r, func_j) * a1(r)
    double sum = 0;
    int num    = 100;

    double local_step = (grid[elem_j] - grid[elem_i]) / num;

    for(int i = 0; i < num; i++)
        sum += 0.5 * (a1(grid[min(elem_i, elem_j)] + i * local_step) + a1(grid[min(elem_i, elem_j)] + (i + 1) * local_step)) * local_step;

    if(func_i != func_j)
        sum *= -1;

    return sum / (step * step);
}

double Integral_Of_a1_Correct_On_Element(int elem_i, int elem_j, int func_i, int func_j){  
    double sum = 0;
    int num    = 100;

    double local_step = (grid[elem_j] - grid[elem_i]) / num;

    for(int i = 0; i < num; i++)
        sum += 0.5 * (B(grid[elem_i] +       i * local_step, func_i) * Da1(grid[elem_i] +       i * local_step)\
                    + B(grid[elem_i] + (i + 1) * local_step, func_i) * Da1(grid[elem_i] + (i + 1) * local_step))\
                    * local_step;

    if(func_j == min(elem_j, elem_i))
      sum *= -1;

    return sum / step;
}

double Integral_Of_a2_On_Element(int elem_i, int elem_j, int func_i, int func_j){   //on [elem_i, elem_j] integrate B(r, func_i) * B(r, func_j) * a3(r)
    double sum = 0;
    int num    = 100;

    double local_step = (grid[elem_j] - grid[elem_i]) / num;

    for(int i = 0; i < num; i++)
        sum += 0.5 * (a2(grid[elem_i] +       i * local_step) * B(grid[elem_i] +       i * local_step, func_i) \
                    + a2(grid[elem_i] + (i + 1) * local_step) * B(grid[elem_i] + (i + 1) * local_step, func_i))\
                    * local_step;
    
    if(func_j == min(elem_j, elem_i))
        sum *= -1;

    return sum / step;
}

double Integral_Of_a3_On_Element(int elem_i, int elem_j, int func_i, int func_j){   //on [elem_i, elem_j] integrate B(r, func_i) * B(r, func_j) * a3(r)
    double sum = 0;
    int num    = 100;

    double local_step = (grid[elem_j] - grid[elem_i]) / num;

    for(int i = 0; i < num; i++)
        sum += 0.5 * (a3(grid[elem_i] +       i * local_step) * B(grid[elem_i] +       i * local_step, func_i) * B(grid[elem_i] +       i * local_step, func_j)  \
                    + a3(grid[elem_i] + (i + 1) * local_step) * B(grid[elem_i] + (i + 1) * local_step, func_i) * B(grid[elem_i] + (i + 1) * local_step, func_j)) \
                    * local_step;
    return sum;
}

double Integral_Of_b_On_Element(int elem_i, int elem_j, int func_i){   //on [elem_i, elem_j] integrate rigth part * B(r, func_i)
    double sum = 0;
    int num    = 100;

    double local_step = (grid[elem_j] - grid[elem_i]) / num;

    for(int i = 0; i < num; i++)
        sum += 0.5 * (b(grid[elem_i] +       i * local_step) * B(grid[elem_i] +       i * local_step, func_i)  \
                    + b(grid[elem_i] + (i + 1) * local_step) * B(grid[elem_i] + (i + 1) * local_step, func_i)) \
                    * local_step;
    return sum;
}

std::vector<double> Solving_SLAE(double** a, double* b){
    double div;

    for(int i = 0; i < nodes - 1; i++){
        div = a[i + 1][i] / a[i][i];
        a[i + 1][i]     -= a[i][i] * div;
        a[i + 1][i + 1] -= a[i][i + 1] * div;
        b[i + 1]        -= b[i] * div;
    }

    std::vector<double> res;
    res.resize(nodes);

    res[nodes - 1] = b[nodes - 1] / a[nodes - 1][nodes - 1];
    for(int i = nodes - 2; 0 <= i; i--)
        res[i] = (b[i] - a[i][i + 1] * res[i + 1]) / a[i][i];

    return res;
}   

std::vector<double> FEM2(){

    double** K = new double*[nodes];
    for(int i = 0; i < nodes; i++)
        K[i] = new double[nodes];

    double* f = new double[nodes];

    for(int i = 0; i < nodes; i++){
        for(int j = 0; j < nodes; j++)
            K[i][j] = 0;
        f[i] = 0;
    }

    for(int i = 0; i < nodes - 1; i++){
        K[i][i]         += Integral_Of_a1_On_Element(i, i + 1, i, i)         + Integral_Of_a1_Correct_On_Element(i, i + 1, i, i)\
                         - Integral_Of_a2_On_Element(i, i + 1, i, i)         - Integral_Of_a3_On_Element(i, i + 1, i, i);

        K[i][i + 1]     += Integral_Of_a1_On_Element(i, i + 1, i, i + 1)     + Integral_Of_a1_Correct_On_Element(i, i + 1, i, i + 1) \
                         - Integral_Of_a2_On_Element(i, i + 1, i, i + 1)     - Integral_Of_a3_On_Element(i, i + 1, i, i + 1);
        
        K[i + 1][i]     += Integral_Of_a1_On_Element(i, i + 1, i + 1, i)     + Integral_Of_a1_Correct_On_Element(i, i + 1, i + 1, i) \
                         - Integral_Of_a2_On_Element(i, i + 1, i + 1, i)     - Integral_Of_a3_On_Element(i, i + 1, i + 1, i);
        
        K[i + 1][i + 1] += Integral_Of_a1_On_Element(i, i + 1, i + 1, i + 1) + Integral_Of_a1_Correct_On_Element(i, i + 1, i + 1, i + 1)\
                         - Integral_Of_a2_On_Element(i, i + 1, i + 1, i + 1) - Integral_Of_a3_On_Element(i, i + 1, i + 1, i + 1);
    }
       
    for(int i = 0; i < nodes - 1; i++){
        f[i]     -= Integral_Of_b_On_Element(i, i + 1, i);
        f[i + 1] -= Integral_Of_b_On_Element(i, i + 1, i + 1);
    }

    f[0]         = boundary_1;      //including Dirichlet boundary conditions
    f[nodes - 1] = boundary_2;
    K[0][0] = 1;
    K[nodes - 1][nodes - 1] = 1;
    K[0][1] = 0;
    K[nodes - 1][nodes - 2] = 0;

    f[1]         -= K[1][0] * boundary_1; 
    f[nodes - 2] -= K[nodes - 2][nodes - 1] * boundary_2;

    K[1][0]                 = 0;
    K[nodes - 2][nodes - 1] = 0;

    return Solving_SLAE(K, f);
}

std::vector<double> Q_Func(std::vector<double> U_Func){   
    std::vector<double> q = {};
    
    for(int i = 0; i < grid.size(); i++)
        q.push_back(E * h[i] * (grid[i] * Du[i] + mu * u[i]) / (1 - mu * mu));

    return q;
}

void Sigma_On_Grid(){       
    Sigma = {};
    double a, b;
    for(int i = 0; i < grid.size(); i++){
        a = Q[i] / (h[i] * grid[i]);
        b = nu * Q[i] / (h[i] * grid[i]) + u[i] * E / grid[i];
        Sigma.push_back(sqrt(pow(a, 2) + pow(b, 2) - a * b));
    }
}

double f(int i, double x){ return -pow(x * grid[i] * sigma_max, 2) + pow(Q[i], 2) + pow(nu * Q[i] + u[i] * x * E, 2) - Q[i] * (nu * Q[i] + u[i] * x * E); };

double hi_next(int i){ return h[i] - f(i, h[i]) / ((f(i, h[i] + 1e-6) - f(i, h[i] - 1e-6)) / (2 * 1e-6) ); };

void Disk_Optimization(){
    u = FEM2();
    //std::cout << u << std::endl;

    Du = Get_Du();
    //std::cout << Du < std::endl;
    
    Q = Q_Func(u);

    Sigma_On_Grid();

    for(int i = 0; i < h.size(); i++)
        //if(hi_next(i) < h_min(grid[i]))
        //    h[i] = h_min(grid[i]);
        //else if(hi_next(i) > h_max(grid[i]))
        //    h[i] = h_max(grid[i]);
        //else
            h[i] = hi_next(i);
            //if(Sigma[i] > sigma_max)
            //    h[i] = hi_next(i);
           

    //std::cout << h << std::endl;

    std::ofstream fout("/home/rus/Desktop/CourseWork/CWSemestr7/u.txt");   //Output in file
    fout << "{\n";
    for(int i = 0; i < u.size(); i++)
        if(i != u.size() - 1)
            fout << "{" << grid[i] << ", " << u[i] << "},\n";
        else
            fout << "{" << grid[i] << ", " << u[i] << "}\n}";

    
    std::ofstream fout4("/home/rus/Desktop/CourseWork/CWSemestr7/Du.txt");   //Output in file
    fout4 << "{\n";
    for(int i = 0; i < grid.size(); i++)
        if(i != grid.size() - 1)
            fout4 << "{" << grid[i] << ", " << Du[i] << "},\n";
        else
            fout4 << "{" << grid[i] << ", " << Du[i] << "}\n}";

    std::ofstream fout1("/home/rus/Desktop/CourseWork/CWSemestr7/Profile.txt");
    fout1 << "{\n";
    for(int i = 4; i < h.size(); i++)
        if(i != h.size() - 1)
            fout1 << "{" << grid[i] << ", " << h[i] << "},\n";
        else
            fout1 << "{" << grid[i] << ", " << h[i] << "}\n}";
    
    std::ofstream fout2("/home/rus/Desktop/CourseWork/CWSemestr7/Q.txt");
    fout2 << "{\n";

    for(int i = 0; i < Q.size(); i++)
        if(i != Q.size() - 1)
            fout2 << "{" << grid[i] << ", " << Q[i] << "},\n";
        else
            fout2 << "{" << grid[i] << ", " << Q[i] << "}\n}";

    Sigma_On_Grid();
   
    std::ofstream fout3("/home/rus/Desktop/CourseWork/CWSemestr7/Sigma.txt");
    fout3 << "{\n";

    for(int i = 0; i < Sigma.size(); i++)
        if(i != Sigma.size() - 1)
            fout3 << "{" << grid[i] << ", " << Sigma[i] * 1e-9 << "},\n";
        else
            fout3 << "{" << grid[i] << ", " << Sigma[i] * 1e-9 << "}\n}";
}

double mass(){
    double sum = 0;

    for(int i = 4; i < h.size() - 1; i++)
        sum += 0.5 * (h[i] * grid[i] + h[i + 1] * grid[i + 1]) * step;

    return sum * rho * 2 * pi;
}

int main(){ 
    for(int i = 0; i < nodes; i++)
        grid.push_back(begin + step * i);
   
    Init_h();

    //std::cout << maximum().first << "\t" << begin + maximum().second * step << std::endl;

    int i = 0;
    m.push_back(mass());
    do {
        Disk_Optimization();
        m.push_back(mass());
        i++;
    } while (m[i - 1] - m[i] > 10);

    std::cout << m << std::endl;

    std::ofstream foutM("/home/rus/Desktop/CourseWork/CWSemestr7/mass.txt");
    foutM << "{\n";

    for(int i = 0; i < m.size(); i++)
        if(i != m.size() - 1)
            foutM << "{" << i << ", " << m[i] << "},\n";
        else
            foutM << "{" << i << ", " << m[i] << "}\n}";
    
    //std::cout << maximum().first << "\t" << begin + maximum().second * step << std::endl;

    //Секция (1) для определения начальных данных проекта
    //1 (sigma ~ 85 -- 1, other: 25-12 GPa, m = ~73 kg)
/*    u = FEM2();             
    Du = Get_Du();
    Q = Q_Func(u);

    std::ofstream fout("/home/rus/Desktop/CourseWork/CWSemestr7/u.txt");   //Output in file
    fout << "{\n";
    for(int i = 0; i < u.size(); i++)
        if(i != u.size() - 1)
            fout << "{" << grid[i] << ", " << u[i] << "},\n";
        else
            fout << "{" << grid[i] << ", " << u[i] << "}\n}";

    Sigma_On_Grid();

    std::ofstream fout3("/home/rus/Desktop/CourseWork/CWSemestr7/Sigma.txt");
    fout3 << "{\n";

    for(int i = 0; i < Sigma.size(); i++)
        if(i != Sigma.size() - 1)
            fout3 << "{" << grid[i] << ", " << Sigma[i] * 1e-9 << "},\n";//1e-9
        else
            fout3 << "{" << grid[i] << ", " << Sigma[i] * 1e-9 << "}\n}";

    std::cout << mass() << std::endl;

    std::ofstream fout1("/home/rus/Desktop/CourseWork/CWSemestr7/Profile.txt");
    fout1 << "{\n";
    for(int i = 0; i < h.size(); i++)
        if(i != h.size() - 1)
            fout1 << "{" << grid[i] << ", " << h[i] << "},\n";
        else
            fout1 << "{" << grid[i] << ", " << h[i] << "}\n}";

    //1********/

    return 0;
}
