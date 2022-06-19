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
   
   double bound;
   int is_leaf;
   long int row;
   long int numinrow;
};
#endif /* Node_hpp */
