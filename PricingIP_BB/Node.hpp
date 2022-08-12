//
//  Node.hpp
//  PricingIP
//
//  Created by StephanP on 5/17/22.
//

#ifndef Node_hpp
#define Node_hpp

#include <stdio.h>

class Node{
   
public:
   Node(const double, const int, const long int, const long int);
   Node(const long int, const long int);
   Node();
   Node(const double, const int, const long int, const long int, const int nexti, const int nextj);
   
   double bound;
   int is_leaf;
   long int row;
   long int numinrow;
   int nexti = -1;
   int nextj = -1;
};
#endif /* Node_hpp */
