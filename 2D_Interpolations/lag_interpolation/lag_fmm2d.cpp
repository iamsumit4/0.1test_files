#include<iostream>
#include<Eigen/Dense>

//--------Kernel function------------//
double KERNEL(Eigen::Vector2d Xnode, Eigen::Vector2d Ynode) {
    return 1 / (Xnode - Ynode).norm();
}

//-------Generating linspace nodes in [a, b]-------//
void get_linspace_nodes(double a, double b, Eigen::VectorXd& node){
    int n = node.size();
    for(int i=0; i<n; i++){
        node(i) = a + i * (b - a) / (n-1);
    }
}

//-------Generating matrix of linspace nodes in [a, b] x [c, d]-------//
void get_2d_linspace_nodes(double a, double b, double c, double d, Eigen::MatrixX<Eigen::Vector2d>& node_matrix){
    int row = node_matrix.rows();
    int col = node_matrix.cols();

    Eigen::VectorXd X_nodes(col), Y_nodes(row);
    get_linspace_nodes(a, b, X_nodes);
    get_linspace_nodes(c, d, Y_nodes);
    for(int i=0; i<row; i++){
        for(int j=0; j<col; j++){
            node_matrix(i,j) << X_nodes(i), Y_nodes(j);
        }
    }
}

//**********Generate Chebyshev nodes in [a, b]**************//
 void get_cheb_node(double a, double b, Eigen::VectorXd& cheb_node){
    int n = cheb_node.size();
    for(int i=0; i<n; i++){
        cheb_node(i) = (0.5 * (a + b)) + (0.5 * (b-a) * cos(M_PI - (((2 * i + 1) * M_PI) / (2 * n))));
    }
}

//-------Generating matrix of Chebyshev nodes in [a, b] x [c, d]-------//
void get_2d_cheb_nodes(double a, double b, double c, double d, Eigen::MatrixX<Eigen::Vector2d>& cheb_node_matrix){
    int row = cheb_node_matrix.rows();
    int col = cheb_node_matrix.cols();

    Eigen::VectorXd X_nodes(col), Y_nodes(row);
    get_linspace_nodes(a, b, X_nodes);
    get_linspace_nodes(c, d, Y_nodes);
    for(int i=0; i<row; i++){
        for(int j=0; j<col; j++){
            cheb_node_matrix(i,j) << X_nodes(i), Y_nodes(j);
        }
    }
}



int main(){
    double a=-3, b=-1;    //source region [-3,-1] x [-3,-1]
    double c=1, d=3;      //target region [1,3] x [1,3]

    int N = 100, M = 100;
    int p = 4;

    Eigen::MatrixX<Eigen::Vector2d> source_matrix(N,N), target_matrix(M,M), source_cheb_matrix(p,p), target_cheb_matrix(p,p);

    get_2d_linspace_nodes(a, b, a, b, source_matrix);
    get_2d_linspace_nodes(c, d, c, d, target_matrix);

    get_2d_cheb_nodes(a, b, a, b, source_cheb_matrix);
    get_2d_cheb_nodes(c, d, c, d, target_cheb_matrix);

    std::cout << source_matrix(0,2) << std::endl;
    std::cout << target_matrix(0,2) << std::endl;






    return 0;
}


   // Eigen::MatrixX<Eigen::Vector2d> matrix1(2,2), matrix2(2,2); 
    // Eigen::MatrixXd kermat(2,2); 

    // for(int i=0; i<matrix1.rows(); i++){
    //     for(int j=0; j<matrix1.cols(); j++){
    //         matrix1(i,j) << 3 * (i+j+0), 4 * (i+j+1);
    //     }
    // }

    // for(int i=0; i<matrix2.rows(); i++){
    //     for(int j=0; j<matrix2.cols(); j++){
    //         matrix2(i,j) << 2 * (i+j+4), 3 * (i+j+5);
    //     }
    // }

    // std::cout << matrix1(0,1)(0) <<std::endl;
    // //std::cout << matrix2(1,1) <<std::endl;
    // //std::cout << KERNEL(matrix1(0,0),matrix2(0,0)) <<std::endl;

    // for(int i=0; i<kermat.rows(); i++){
    //     for(int j=0; j<kermat.cols(); j++){
    //         kermat(i,j) = KERNEL(matrix1(i,j),matrix2(i,j));
    //     }
    // }

    // std::cout << "The kermat is\n" << kermat << std::endl;
    // std::cout << "\nsize " << sqrt(kermat.size());