#include <iostream>
#include <Eigen/Dense>
#include <chrono>



// Kernel function
double KERNEL(double x, double y) {
    return 1 / fabs(x - y);
}

//**************Generate linspace nodes in [a, b]***********//
void get_linspace_node(double a, double b, Eigen::VectorXd& node) {
    int n = node.size();
    for (int i = 0; i < n; i++) {
        node(i) = a + i * (b - a) / (n - 1);               //node -----> either source_node or target_node
    }
}

//**********Generate Chebyshev nodes in [a, b]**************//
void get_cheb_node(double a, double b, Eigen::VectorXd& cheb_node){
    int n = cheb_node.size();
    for(int i=0; i<n; i++){
        cheb_node(i) = (0.5 * (a + b)) + (0.5 * (b-a) * cos(M_PI - (((2 * i + 1) * M_PI) / (2 * n))));
    }


}

//************ Shifting function form [a, b] to [-1, 1]*************//
void shift(double a, double b, Eigen::VectorXd& node, Eigen::VectorXd& shiftrd_node){
    int n = node.size();
    for(int i=0; i<n; i++){
        shiftrd_node(i) = -1 + 2 * (node(i) - a) / (b - a); 
    }
}




//************ Generating Chebyshev Polynomial of order k ****************//
double T(int k, double x){              // T_k(x)
    if(k==0){                           // x ----> either cheb_node or node for both source and target
        return 1;
    }else if(k==1){
        return x;
    }else{
        return 2 * x * T(k-1, x) - T(k-2, x);
    }
} 

//******************* Generating Chebyshev multiplier S_p(x,y) *******************//
double S(int p, double node, double cheb_node){             // S_p(node, cheb_node)
    double s = 0;
    for(int k=0; k<p-1; k++){                              //node ----> either target_node or source_node
        s = s + T(k+1, node) * T(k+1, cheb_node);         //cheb_nodeC----> eithet target_cheb_node or source_cheb_node
    }
    return (1 + 2 * s)/p;
}

 //**********Generating U_matrix and V_matrix*************//
 void get_matrix(Eigen::MatrixXd& matrix, Eigen::VectorXd& cheb_node, Eigen::VectorXd& node){
    int n = node.size(); 
    int p = cheb_node.size();                               // node ----> either source_node or target_node
    for(int i=0; i<n; i++){                                 // cheb_node ----> either source_cheb_node or target_cheb_node
        for(int j=0; j<p; j++){
            matrix(i, j) = S(p, node(i), cheb_node(j));
        }
    }
 }

  //**********Genegarte Interection matrix***********//
 void get_interection_matrix(Eigen::MatrixXd& matrix, Eigen::VectorXd& source_node, Eigen::VectorXd& target_node){
    for(int i=0; i<target_node.size(); i++){
        for(int j=0; j<source_node.size(); j++){
            matrix(i,j) = KERNEL(target_node(i) , source_node(j));
        }
    }
 }


int main() {
    int M = 2000, N=2500, p = 15;
    double a = -3, b = -1;   //source interval [-3,-1]
    double c = 1, d = 3;     //target interval [1,3]
    

    Eigen::VectorXd source_node(N), target_node(M), source_cheb_node(p), target_cheb_node(p);
    Eigen::VectorXd cheb_node(p), shifted_source_node(N), shifted_target_node(M);
    Eigen::MatrixXd U_matrix(M, p), K_matrix(p, p), V_matrix(N, p), K_actual_matrix(M, N), K_approx_martix(M, N);



    get_linspace_node(a, b, source_node);
    get_linspace_node(c, d, target_node);
    
    get_cheb_node(a, b, source_cheb_node);
    get_cheb_node(c, d, target_cheb_node);

    shift(a, b, source_node, shifted_source_node);   //source nodes are shifted from [a, b] to [-1, 1]
    shift(c,d, target_node, shifted_target_node);    //target nodes are shifted from [c, d] to [-1, 1]       
    get_cheb_node(-1, 1, cheb_node);                 //cheb_node  gives the chebyshev nodes in [-1, 1]


    // Start the timer
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    get_matrix(U_matrix, cheb_node, shifted_target_node);
    get_matrix(V_matrix, cheb_node, shifted_source_node);

    get_interection_matrix(K_matrix, source_cheb_node, target_cheb_node);
    get_interection_matrix(K_actual_matrix,source_node,target_node);

    K_approx_martix = U_matrix * K_matrix * V_matrix.transpose();

    std::cout << "\nThe relative error in approximating the interection matrix is " << (K_actual_matrix-K_approx_martix).norm()/K_actual_matrix.norm() << std::endl;

    // End the timer
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    // Calculate the duration
    std::chrono::duration<double> duration = end - start;

    // Output the runtime
    std::cout << "\nTime taken to exicute the programe is " << duration.count() << " seconds\n" << std::endl;


    return 0;
}