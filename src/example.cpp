#include "headers/utils.h"
#include "headers/lin_alg.h"
#include "headers/lie_algebra.h"
#include <random>

void example(){
    // Here is an example, feel free to delete it.
    g::symbol x("a"); // This creates a symbol called `a`, stored in a variable x.
    g::symbol y("y"); // This creates a symbol called `y`, stored in a variable y.
    g::symbol z("a"); // This creates a second *distinct* symbol called `a`, stored in a variable z. 
    // Note: x - x = 0, but x - z != 0, even though they both refer to expression called a! 

    g::matrix m(2,2); // This creates an empty 2x2 matrix. Matrices are indexed from 0.
    m(0,0) = 1 + x; // This sets the (0,0) entry to 1 + a.
    m(0,1) = 2; // This sets the (0,1) entry to 2.
    m(1,1) = -1 - x; // This sets the (1,1) entry to -1 - a. Note that this makes sure m is in sl(2), which most algorithms assume.
    g::matrix n = {{-1+y,2},{x,1-y}}; // This makes a 2x2 matrix with first row {1,2} and second row {a,1-y}

    lie_algebra* L1 = new lie_algebra({m,n}); // This creates a lie algebra with generators m and n
    lie_algebra* L2 = new lie_algebra({m}); // This creates a lie algebra with generator m.
    
    std::cout << "Printing Basis of L1" << std::endl; // This will write to the console, "Printing Basis of L1"
    utils::print_matrices(L1->get_basis()); // This will print the basis of L1 to the console.
    std::cout << "Dimension of the algebra L1: " << L1->get_dim() << std::endl;// This shows that the Lie algebra L1 is all of sl(2) since dim L1 = 3 = dim sl(2)
    
    //We can also make new lie algebras from old
    lie_algebra* L2_normalizer = L2->compute_normalizer(); 
    std::cout << "Printing Basis of N(L2)" << std::endl; // This will write to the console, "Printing Basis of L2_normalizer"
    utils::print_matrices(L2_normalizer->get_basis()); // This will print the basis of L2_normalizer to the console.
    std::cout << "Dimension of the algebra N(L2): " << L2_normalizer->get_dim() << std::endl;// This shows that the Lie algebra L2_normalizer is its own normalizer
    
    // You can find all the things you can do with a lie algebra in `header/lie_algebra.h` along with some simple documentation.
}

int main(){
    // Write your code in here!

    // Remove this line when you want to run your code. This is just to run the example code above.
    example();

    // Don't get rid of this line! It tells the program that it has finished correctly.
    return 0;
}