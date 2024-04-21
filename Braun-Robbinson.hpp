#include <algorithm>
#include <boost/format.hpp>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
#include <boost/range/algorithm/count.hpp>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>


class NumberGenerator {
    int current;
public:
    NumberGenerator() : current(0) {}
    int next(int n) {
        int value = current;
        current = (current + 1) % n;
        return value;
    }
};

using namespace boost::numeric::ublas;

// using namespace std;



int braun_robinson(matrix <double> H_arr, int N) {

NumberGenerator generator;


    double eps = 0.1;

    matrix<double> H_arr_obr(H_arr.size1(),H_arr.size1());
    



    permutation_matrix<std::size_t> pm(H_arr.size1());
    

    lu_factorize(H_arr, pm);
    lu_substitute(H_arr, pm, H_arr_obr);
    H_arr_obr.assign(identity_matrix<double>(N+1));

    // std::cout<<H_arr<<std::endl;
    // std::cout<<H_arr_obr<<std::endl;
    // std::cout<<pm<<std::endl;
    



vector<double> u(N+1);
    
   for(int i =0; i<N+1; i++ ){
        u(i)=1;
   }

    
    vector<double> _X(N+1);

        // std::cout<<pm<<std::endl;
        _X = prod(trans(u), H_arr_obr) / inner_prod(prod(trans(u), H_arr_obr), u);
        // std::cout<<"Analitic_X "<<_X<<std::endl;
        vector<double> _Y(N+1);
        _Y = prod(trans(H_arr_obr), trans(u)) / inner_prod(prod(trans(u), H_arr_obr), u);
        // std::cout<<"Analitic_Y "<<_Y<<std::endl;
        double _V = 1 / inner_prod(prod(u, H_arr_obr), trans(u));
        // std::cout<<"Analitic_V "<<_V<<std::endl;
// Corrected calculation of _Y

    std::vector<int> strA;
    std::vector<int> strB;
    std::vector<double> high_cost;
    std::vector<double> low_cost;
    std::vector<double> eps_fact;

    
    boost::random::random_device rng;
    boost::random::uniform_int_distribution<> dist(0, N);

    int num_strA = dist(rng);
    int num_strB = dist(rng);
    strA.push_back(num_strA);
    strB.push_back(num_strB);

    std::vector<double> winA;
    std::vector<double> loseB;


for(int i=0;i<N;i++){
    winA.push_back(H_arr(i, num_strB));
    loseB.push_back(H_arr(num_strA, i));
    // std::cout<<winA[i]<<"  winA"<<i<<std::endl;
    // std::cout<<loseB[i]<<"  loseB"<<i<<std::endl;

    };

    // std::cout<<winA[1]<<"\n"<<winA[2]<<"\n"<<winA[3]<<"\n"<<winA[4]<<"\n"<<winA[5]<<"\n"<<winA[6]<<"\n"<<winA[7]<<"\n"<<winA[8]<<"\n"<<winA[9]<<"\n"<<winA[10]<<std::endl;
    
    high_cost.push_back(*std::max_element(winA.begin(), winA.end()));

    low_cost.push_back(*std::min_element(loseB.begin(), loseB.end()));
    
    eps_fact.push_back(*std::min_element(high_cost.begin(), high_cost.end()) - *std::max_element(low_cost.begin(), low_cost.end()));

    
    
    int k = 1;
    int iteration=0;
   
// std::cout<<eps_fact[k-1]<<std::endl;
    
while ( eps_fact[k-1] > eps) {
        iteration = generator.next(N);
        num_strA = std::max_element(winA.begin(), winA.end()) - winA.begin();
        num_strB = std::min_element(loseB.begin(), loseB.end()) - loseB.begin();
        
        strA.push_back(num_strA);
        strB.push_back(num_strB);
        for(int i=0;i<N;++i){
            winA.push_back(winA[k - 2] + H_arr(i, iteration));
            loseB.push_back(loseB[k - 2] + H_arr(iteration,i));  
        }

        high_cost.push_back(1.0 / k * *std::max_element(winA.begin(), winA.end()));
        low_cost.push_back(1.0 / k * *std::min_element(loseB.begin(), loseB.end()));
        eps_fact.push_back(*std::min_element(high_cost.begin(), high_cost.end()) - *std::max_element(low_cost.begin(), low_cost.end()));
        k++;
    }



    // /* Вывод результатов */
    // std::cout << "\tAnalytical Solution:\n";
    // std::cout << "\t x = " << _X << "\n";
    // std::cout << "\t y = " << _Y << "\n";
    // std::cout << "\t v = " << _V << "\n";

    std::cout << "\t Numerical Solution using the Brown-Robinson algorithm:\n";
    std::cout << "\t x = (" << winA[k-1]<<")\n"<<std::endl;
    std::cout << "\t y = (" << loseB[k-1]<<")\n"<<std::endl;
    std::cout << "\t v = " << *std::min_element(high_cost.begin(), high_cost.end()) - *std::max_element(low_cost.begin(), low_cost.end()) << "\n"<<std::endl;
        return 0; 
}

