//
//  Node.cpp
//  PricingIP
//
//  Created by StephanP on 5/17/22.
//

#include "Node.hpp"

Node::Node(const double b, const int i, const long int r, const long int nir): bound(b), is_leaf(i), row(r), numinrow(nir)
{
   
}

Node::Node(): bound(0), is_leaf(1), row(-1), numinrow(-1)
{
   
}

Node::Node(const long int r, const long int nir): bound(0), is_leaf(1), row(r), numinrow(nir)
{
   
}
