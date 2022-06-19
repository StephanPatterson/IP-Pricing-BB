//
//  SuppPt.hpp
//  
//
//  Created by Stephan Patterson on 1/31/18.
//

#ifndef SuppPt_hpp
#define SuppPt_hpp

#include <iostream>

class SuppPt
{
   friend std::ostream &operator<<(std::ostream &, const SuppPt&);
   friend std::ofstream &operator<<(std::ofstream &, const SuppPt&);
//   friend void swap(SuppPt &, SuppPt &);
   
public:
   SuppPt(const double, const double, double, unsigned long int);
   SuppPt(const double, const double, double);
   SuppPt();
   void setMass(double);
   void print() const;
   SuppPt &operator-=(double);
//   bool operator<( const SuppPt &) const;
   bool lessMass(const SuppPt, const SuppPt);
   const SuppPt &operator=(const SuppPt &);
   double dist(const SuppPt );
   double operator*(const SuppPt);
   
   double loc1;
   double loc2;
   double mass;
   unsigned long int combos;

};


#endif /* SuppPt_hpp */
