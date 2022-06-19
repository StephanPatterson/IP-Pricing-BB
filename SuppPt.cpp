//
//  SuppPt.cpp
//  
//
//  Created by Stephan Patterson on 1/31/18.
//

#include "SuppPt.hpp"
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

SuppPt::SuppPt(double l1, double l2, double m, unsigned long int co): loc1(l1), loc2(l2), mass(m), combos(co)
{
   
}
SuppPt::SuppPt(double l1, double l2, double m): loc1(l1), loc2(l2), mass(m), combos(1)
{
}

SuppPt::SuppPt(): loc1(0), loc2(0), mass(0), combos(0)
{
   
}

void SuppPt::setMass(double m)
{
   if (m < 0 || m > 1)
   {
      throw std::invalid_argument("Mass out of desired range");
   }
   else
   {
      mass = m;
   }
   
}

void SuppPt::print() const
{
   std::cout << "(" << loc1 << ", " << loc2 << ")" << ", " << mass << std::endl;
}

SuppPt &SuppPt::operator-=(double m)
{
   if (m > this->mass)
   {
      std::cout << "Mass to remove will create invalid point." << std::endl;
   }
   this->mass -= m;
   return *this;
}

double SuppPt::operator*(SuppPt y)
{
   double sum = 0;
   sum += this->loc1*y.loc1;
   sum += this->loc2*y.loc2;
   return sum;
}

//bool SuppPt::operator<(const SuppPt &s) const
//{
//      if (this->loc1 < s.loc1)
//      {
//         return true;
//      }
//      else if (this->loc1 > s.loc1)
//      {
//         return false;
//      }
//      else
//      {
//         if (this->loc2 < s.loc2)
//         {
//            return true;
//         }
//         else if (this->loc2 > s.loc2)
//         {
//            return false;
//         }
//         else
//         {
//            return false;
//         }
//      }
//}


const SuppPt &SuppPt::operator=(const SuppPt &right)
{
   loc1 = right.loc1;
   loc2 = right.loc2;
   mass = right.mass;
   combos = right.combos;
   return *this;
}

std::ostream &operator<<(std::ostream &o, const SuppPt& s)
{
   o << "(" << s.loc1 << ", " << s.loc2 << ")" << ", " << s.mass << " " << s.combos;
   return o;
}

std::ofstream &operator<<(std::ofstream &o, const SuppPt& s)
{
   o << s.loc1 << " " << s.loc2 << " " << s.mass;
   return o;
}

double SuppPt::dist(const SuppPt s)
{
   double d = std::sqrt((s.loc1-this->loc1)*(s.loc1-this->loc1)+(s.loc2-this->loc2)*(s.loc2-this->loc2));
   return d;
}
/*
void swap(SuppPt &lhs, SuppPt &rhs)
{
   std::swap(lhs.loc1, rhs.loc1);
   std::swap(lhs.loc2, rhs.loc2);
   std::swap(lhs.mass, rhs.mass);
}*/
