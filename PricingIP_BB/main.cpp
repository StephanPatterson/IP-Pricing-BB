//
//  main.cpp
//  PricingIP_BB
//
//  Created by StephanP on 5/16/22.
//  This version is for implementing Branch & Bound by hand, for experimenting with choices like node selection

#include "SuppPt.hpp"
#include "Node.hpp"
#include <iostream>
#include <vector>
#include "/Library/gurobi951/macos_universal2/include/gurobi_c++.h"
#include <sstream>
#include <fstream>
#include <cmath>

double objval(const double *, const double *, std::vector< std::vector< SuppPt> >&, std::vector<std::vector<double> >&, const int, const int *);
double objval2(const double *, const double *, std::vector< std::vector< SuppPt> >&, std::vector<std::vector<double> >&, const int, const int *);
double objval3(const double *, const double *, std::vector< std::vector< SuppPt> >&, std::vector<std::vector<GRBVar> >&, std::vector<std::vector<std::vector<std::vector<GRBVar> > > >&, const int, const int *);
bool check_integral(std::vector<std::vector< GRBVar> > &, const int, const int *, const double &);

struct ij
{
   int iindex;
   int jindex;
   int setval;
   
   ij(): iindex(-1), jindex(-1), setval(-1)
   {
      
   }
   
   ij(int i,int j, int val): iindex(i), jindex(j), setval(val)
   {
      
   }
};

int main(int argc, const char * argv[]) {
   
   int Pnum =4;//If NumMoMeas = 1, Number of months to process. Will go consecutively from start month, skipping any empty months.
   std::string startmonth = "Jul";
   int startyear = 15;
   int NumMoMeas = 1;
   int Psnum[Pnum];
   std::vector< std::vector<SuppPt> > PsuppT(Pnum);
   int totalsupp = 0;
   int indices[Pnum];//The use of indices is changing (from previous code) to match the two-indexed vector representation of the support sets
   double lambda[Pnum];

   std::string temp;

   std::ifstream indata;
   std::ostringstream fname;
   int year;
   std::string month;
   fname << "/Users/spatterson/ForXcode/Barycenter/DenverCrime.csv";
   indata.open(fname.str());
   if (!indata.good())
   {
      std::cout << "File not found." << std::endl;
      return 1;
   }
   getline(indata, temp, '\n');// Advance past headers
   indata >> year;
   if (startyear < year)
   {
      std::cout << "Invalid Year Selection." << std::endl;
      return 1;
   }
   while (year < startyear)
   {
      getline(indata, temp, '\n');// Advance to next entry
      indata >> year;
   }
   indata >> month;
   while (month != startmonth)
   {
      getline(indata, temp, '\n');
      indata >> year;
      indata >> month;
      if (year != startyear)
      {
         std::cout << "Selected Month/Year not in data." << std::endl;
         return 1;
      }
      if (indata.eof())
      {
         indata.close();
         std::cout << "Invalid Input. End of file reached." << std::endl;
         return 1;
      }
   }
   
   double loc1;
   double loc2;
   indata >> loc1;
   indata >> loc2;
   std::string currentmonth = startmonth;
         
   int tempind = 0;
   for (int i = 0; i < Pnum; ++i)
   {
      Psnum[i] = 0;
      indices[i] = 0;
      for (int j = 0; j < NumMoMeas; ++j)
      {
         while (month == currentmonth )
         {
            //Adding a shift so that all coordinates are positive
            PsuppT[i].push_back(SuppPt(loc1+115, loc2, 1.0, totalsupp));
            ++tempind;
            ++Psnum[i];
            ++totalsupp;
                  
            indata >> year;
            indata >> month;
            indata >> loc1;
            indata >> loc2;
            if (indata.eof())
            {
               indata.close();
               std::cout << "Invalid Input. End of file reached." << std::endl;
               return 1;
            }
         }
         currentmonth = month;
      }
      std::cout << "Month measure: " << i+1 << " Size: " << Psnum[i] << std::endl;
      double totalmass = 0;
      for (int j = 0; j < Psnum[i]-1; ++j)
      {
         PsuppT[i][j].mass /= Psnum[i];
         totalmass += PsuppT[i][j].mass;
      }
      PsuppT[i][Psnum[i]-1].mass = 1-totalmass;
      currentmonth = month;
   }
   indata.close();
   double sum = 0;
   for (int i = 0; i < Pnum-1; ++i)
   {
      lambda[i] = (double)Psnum[i]/totalsupp;
      sum += lambda[i];
   }
   lambda[Pnum-1] = 1-sum;
   //Omitting validation code from previous scripts

   //Calculate number of possible supp pts S0
   long int S0 = 1;
   for (int i = 0; i < Pnum; ++i)
   {
      S0 *= Psnum[i];
   }
   std::cout << "Size of S0: " << S0 << std::endl;

  //Now calculate costs.
   std::vector<double> c(S0,0);
   int index = 0;
   for (unsigned long int j = 0; j < S0; ++j)
   {
      double sum1 = 0;
      double sum2 = 0;
      for (int i = 0; i < Pnum; ++i)
      {
         index = indices[i];
         sum1 += lambda[i]*PsuppT[i][index].loc1;
         sum2 += lambda[i]*PsuppT[i][index].loc2;
      }
         
      SuppPt Pbar0 = SuppPt(sum1, sum2, 0.0);
      for (int i = 0; i < Pnum; ++i)
      {
         index = indices[i];
         c[j] += lambda[i]*((Pbar0.loc1-PsuppT[i][index].loc1)*(Pbar0.loc1-PsuppT[i][index].loc1) +(Pbar0.loc2-PsuppT[i][index].loc2)*(Pbar0.loc2-PsuppT[i][index].loc2));
      }
         
      int k = Pnum-1;
      if (indices[k] < Psnum[k]-1)
      {
         ++indices[k];
      }
      else
      {
         int temp = k-1;
         while (temp >= 0)
         {
            if (indices[temp] == Psnum[temp]-1)
            {
               temp -= 1;
            }
            else
            {
               ++indices[temp];
               break;
            }
         }
         for (int l = k; l > temp; --l)
         {
            indices[l] = 0;
         }
      }
   }
   
   for (int i = 0; i < Pnum; ++i)
   {
      indices[i] = 0;
   }
   
   //Set up Master Problem. Need a y from master problem for setting up objective function for IP pricing
   GRBEnv* env = new GRBEnv();
   GRBModel* model = new GRBModel(*env);
   model->set(GRB_IntParam_Method, 0);
   model->set(GRB_IntParam_Presolve, 0);
   model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
//   model->set(GRB_IntParam_OutputFlag, 0); //Turning off display on screen. Disable to see iterations, objective value, etc.
   
   std::vector< GRBVar > w;
   GRBLinExpr  exp[totalsupp];

   std::vector< std::vector<SuppPt> >::iterator Psuppit = PsuppT.begin();
   std::vector<std::vector<SuppPt> > Psupp2(Pnum);
   std::copy(Psuppit, Psuppit+Pnum, Psupp2.begin() );
   while (Psupp2[0][Psnum[0]-1].mass >1e-15)
   {
      double minmass = 1;
      for (int i = 0; i < Pnum; ++i)
      {
         if (Psupp2[i][indices[i]].mass < minmass)
         {
            minmass = Psupp2[i][indices[i]].mass;
         }
      }
      std::cout << "Minimum mass is " << minmass <<std::endl;
      
      long int loocur = S0;
      int unki = 0;
      int index = 0;
      for (int i = 0; i < Pnum; ++i)
      {
         unsigned long int index2 = Psupp2[i][indices[i]].combos-index;
//         std::cout << Psupp2[i][indices[i]].combos << " " << index << std::endl;
//         std::cout << "Index in Pi: " << Psupp2[i][indices[i]].combos - index <<std::endl;
         loocur /= Psnum[i];
         unki += loocur*index2;
         index += Psnum[i];
         Psupp2[i][indices[i]].mass -= minmass;
         if (Psupp2[i][indices[i]].mass <= 1e-15)
         {
            ++indices[i];
         }
      }
/*      for (int i = 0; i < Pnum; ++i)
      {
         std:: cout << Psupp2[i][indices[i]].mass << " ";
      }
      std::cout << std::endl;
      
      std::cout << "current index is " << unki << " out of " << S0-1 <<std::endl;*/
      w.push_back(model->addVar(0.0, GRB_INFINITY, c[unki], GRB_CONTINUOUS));
      loocur = S0/Psnum[0];
      int startindex = Psnum[0];
      int jindex = floor(unki/loocur);
      exp[jindex] += *(w.end()-1);
      for (int l = 1; l < Pnum-1; ++l)
      {
         loocur /= Psnum[l];
         jindex = startindex+floor( (unki % (loocur*Psnum[l]))/loocur);
         exp[jindex] += *(w.end()-1);
         startindex += Psnum[l];
      }
      jindex = startindex+ unki % Psnum[Pnum-1];
      exp[jindex] += *(w.end()-1);
   }
   index = 0;
   for (int i = 0; i < Pnum; ++i)
   {
      for (int j = 0; j < Psnum[i]; ++j)
      {
         model->addConstr(exp[index] == PsuppT[i][j].mass);
         ++index;
      }
   }
   model->optimize();
   double * yhat = model->get(GRB_DoubleAttr_Pi, model->getConstrs(), totalsupp);

   std::vector<Node> Tree_to_Process;
   
   GRBModel* basemodel = new GRBModel(*env);
   basemodel->set(GRB_IntParam_Method, 0);
   basemodel->set(GRB_IntAttr_ModelSense, -1);//set to maximize

   std::vector<std::vector<double> > Best_Found_Solution(Pnum);
   double tol = 1e-4;
   
   //This will generate an initial solution from the first point in each measure
   // giving an initial lower bound
   for (int i = 0; i < Pnum; ++i)
   {
      Best_Found_Solution[i].push_back(1);
      for (int j = 1; j < Psnum[i]; ++j)
      {
         Best_Found_Solution[i].push_back(0);
      }
   }
   //Calling objective function to compute corresponding value, which is a lower bound on optimal
   double bestbound = objval(yhat, lambda, PsuppT, Best_Found_Solution, Pnum, Psnum );
   std::cout << "Initial best objective value is " << bestbound <<std::endl;
   index = 0;
   
   int Psmax = *std::max_element(Psnum, Psnum+Pnum);
   std::vector< std::vector<GRBVar > > z2(Pnum, std::vector<GRBVar>((Psmax)));
   
   std::vector< std::vector< std::vector< std::vector< GRBVar >> >>z4(Pnum,std::vector< std::vector< std::vector< GRBVar >>>(Psmax, std::vector< std::vector< GRBVar >>(Pnum, std::vector< GRBVar >(Psmax))));
   
   int currentconst = 0;
   int totalconstraints = 0;
   for (int i = 0; i < Pnum; ++i)
   {
      GRBLinExpr expp;
      for (int j = 0; j < Psnum[i]; ++j)
      {
         std::stringstream name;
         name << "z" << i+1 << "_" << j+1;
         std::string name2;
         name >> name2;
         z2[i][j] = basemodel->addVar(0.0,1.0,1,GRB_CONTINUOUS,name2);
         expp += z2[i][j];
         ++currentconst;
      }
      basemodel->addConstr(expp == 1);
      ++totalconstraints;
   }
   
   GRBLinExpr expp2[totalsupp*totalsupp*2];
   int constrcol = 0;
   for (int i = 0; i < Pnum-1; ++i)
   {
      for (int j = i+1; j < Pnum; ++j)
      {
         for (int k = 0; k < Psnum[i]; ++k)
         {
            for (int m = 0; m < Psnum[j]; ++m)
            {
               std::stringstream name;
               name << "z" << i+1 <<"_" <<k+1 << "_" <<j+1 << "_" <<m+1;
               std::string name2;
               name >> name2;
               double prod = 2*lambda[i]*lambda[j]*(PsuppT[i][k]*PsuppT[j][m]);
               z4[i][k][j][m] = basemodel->addVar(0.0,GRB_INFINITY,prod,GRB_CONTINUOUS,name2);
               expp2[constrcol]+= z4[i][k][j][m]-z2[i][k];
               ++constrcol;
               expp2[constrcol]+= z4[i][k][j][m]-z2[j][m];
               ++constrcol;
            }
         }
      }
   }

   //Add remaining constraints to model
   for (int k = 0; k < constrcol; ++k)
   {
      basemodel->addConstr(expp2[k] <= 0);
      ++totalconstraints;
   }
   // add/update the coefficient on z[i][j]
   std::vector< std::vector< double > > zikcost(Pnum, std::vector<double> ((Psmax),0));
   std::vector< std::vector< double > > zjmcost(Pnum, std::vector<double> ((Psmax),0));
   for (int i = 0; i < Pnum-1; ++i)
   {
      for (int j = i+1; j < Pnum; ++j)
      {
         for (int k = 0; k < Psnum[i]; ++k)
         {
            double prod = PsuppT[i][k]*PsuppT[i][k];
            zikcost[i][k] += lambda[i]*lambda[j]*prod;
         }
      }
   }
   for (int i = 0; i< Pnum-1; ++i)
   {
      for (int j = i+1; j < Pnum; ++j)
      {
         for (int m = 0; m < Psnum[j]; ++m)
         {
            double prod = PsuppT[j][m]*PsuppT[j][m];
            zjmcost[j][m] += lambda[i]*lambda[j]*prod;
         }
      }
   }
   currentconst = 0;
   for (int i = 0; i < Pnum; ++i)
   {
      for (int j = 0; j < Psnum[i]; ++j)
      {
         z2[i][j].set(GRB_DoubleAttr_Obj,yhat[currentconst]-zikcost[i][j]-zjmcost[i][j]);
         ++currentconst;
      }
   }
   
   //This code should be removed from final version. Using for testing/debugging.
/*   std::vector< std::vector<double > > ztest(Pnum);
   for (int i = 0; i < Pnum; ++i)
   {
      for (int j = 0; j < Psnum[i]; ++j)
      {
         std::cout << "Z_"<<i+1 << j+1 << " " << z2[i][j].get(GRB_DoubleAttr_X) <<std::endl;
         ztest[i].push_back(z2[i][j].get(GRB_DoubleAttr_X));
      }
   }*/
/*   for (int i = 0; i < Pnum-1; ++i)
   {
      for (int j = i+1; j < Pnum; ++j)
      {
         for (int k = 0; k < Psnum[i]; ++k)
         {
            for (int m = 0; m < Psnum[j]; ++m)
            {
               std::cout <<z4[i][k][j][m].get(GRB_DoubleAttr_X) << " " << z2[i][k].get(GRB_DoubleAttr_X) << " " << z2[j][m].get(GRB_DoubleAttr_X) << std::endl;
            }
         }
      }
   }*/
   
//   std::cout << objval(yhat,lambda,PsuppT,ztest,Pnum,Psnum) << " " << objval2(yhat,lambda,PsuppT,ztest,Pnum,Psnum) << " " << objval3(yhat, lambda, PsuppT, z2, z4, Pnum, Psnum) << " " << basemodel->get(GRB_DoubleAttr_ObjVal) << std::endl;
   
   //Solve full linear program (essentially initial node) and check for integrality; if integral, solution found
   basemodel->optimize();
   
   //Can use this for part of paper talking about how fractional the vertices are
   if (basemodel->get(GRB_IntAttr_Status) == 2)
   {
      if (check_integral(z2,Pnum,Psnum,tol))
      {
         std::cout << "Initial Linear Program produces integral solution. No Branch & Bound required." <<std::endl;
         return 0;
      }
   }
   
   //Branch Index will store the measure number, index within the measure, and the variable's assigned value for each step of the tree. These will likely not match the row number and column number within the tree, which is what the two indices of Branch Index are.
   //So if you initially decide to branch on the variable that goes with the first support point of measure 3, BranchIndex[0][0] = (2,0), and BranchIndex[1][0] = (2,0,0) and BranchIndex[1][1] = (2,0,1)
   std::vector< std::vector< ij > > BranchIndex(2);//Pnum is the minimum number of rows to reach a feasible solution, but starting with 2 and allocating additional space as needed
   BranchIndex[0].resize(1);
   BranchIndex[0][0]=ij(0,0,-1); //BranchIndex[0] is the initial decision. BranchIndex[1] will contain the two decisions from the first row, etc.
   BranchIndex[1].resize(2);
   
   //Begin setup of 2 scenarios for computing LPs for both child nodes
   basemodel->set(GRB_IntAttr_NumScenarios, 2);
   
   basemodel->set(GRB_IntParam_ScenarioNumber, 0);
   Node node1 = Node(1,0);
   BranchIndex[1][0] = ij(0,0,0);
   z2[0][0].set(GRB_DoubleAttr_ScenNUB, 0);
   
   basemodel->set(GRB_IntParam_ScenarioNumber, 1);
   Node node2 = Node(1,1);
   BranchIndex[1][1] = ij(0,0,1);
   z2[0][0].set(GRB_DoubleAttr_ScenNLB, 1);
   
   basemodel->optimize();
   
   int scenstatus = basemodel->get(GRB_IntAttr_Status);
   if (scenstatus == 2) //Optimal solution to LP found
   {
      //Store value as best potential of branches of the node
      node2.bound = basemodel->get(GRB_DoubleAttr_ScenNObjVal);
      std::cout << node2.bound <<std::endl;
      
      //If the optimal value of the LP isn't better than the best found solution, the node's branches do not need to be explored further
      if (node2.bound > bestbound)
      {
         //Then check for integrality
         // If integral, new best solution found
         if (check_integral(z2,Pnum,Psnum,tol))
         {
            bestbound = node2.bound;
            for (int i = 0; i < Pnum; ++i)
            {
               for (int j = 0; j < Psnum[i]; ++j)
               {
                  Best_Found_Solution[i][j] = z2[i][j].get(GRB_DoubleAttr_X);
               }
            }
         }
         else // If non-integral, add to tree for branching
         {
            node2.is_leaf = 1;
            Tree_to_Process.push_back(node2);
         }
      }
   }
   else if (scenstatus == 3 or scenstatus == 4)
   {
      node2.bound = 0;
      node2.is_leaf = 0;
   }
   else
   {
      std::cout << "Unexpected Scenario Status Reached: Status " << scenstatus << std::endl;
      return -1;
   }
   
   basemodel->set(GRB_IntParam_ScenarioNumber, 0);
   if (basemodel->get(GRB_IntAttr_Status) == 2)
   {
      node1.bound = basemodel->get(GRB_DoubleAttr_ScenNObjVal);
      std::cout << node1.bound <<std::endl;
      
      if (node1.bound > bestbound)
      {
         // If integral, new best solution found
         if (check_integral(z2,Pnum,Psnum,tol))
         {
            bestbound = node1.bound;
            for (int i = 0; i < Pnum; ++i)
            {
               for (int j = 0; j < Psnum[i]; ++j)
               {
                  Best_Found_Solution[i][j] = z2[i][j].get(GRB_DoubleAttr_X);
               }
            }
         }
         else // If non-integral, add to tree for branching
         {
            node1.is_leaf = 1;
            Tree_to_Process.push_back(node1);
         }
      }
   }
   
   //Now choose new node
   //while Tree vector isn't empty
   while (!Tree_to_Process.empty())
   {
      Node processingnode;
      bool gettingnode = true;
      while (gettingnode)
      {
         gettingnode = false;
         processingnode = Tree_to_Process.back();
         //This guarantees we process a node with the potential to improve, and that has leaves. No nodes without leaves should be added to the tree, so the check should be redundant
         if (processingnode.bound <= bestbound or processingnode.is_leaf == 0)
         {
            Tree_to_Process.pop_back();
            gettingnode = true;
         }
      }

      Node node1 = Node(processingnode.row+1,2*processingnode.numinrow);
      Node node2 = Node(processingnode.row+1,2*processingnode.numinrow+1);
         
      //If working with a new row of the tree, allocate space
      long int rowsize = pow(2,node1.row);
      if (BranchIndex.size() <= node1.row)
      {
         std::vector< ij > newvec(rowsize);
         BranchIndex.push_back(newvec);
      }
      /*
      for (int i = 0; i < rowsize; ++i)
         std::cout << BranchIndex[BranchIndex.size()-1][i].setval << std::endl;*/
         
      ij newij;
      
      //Now decide a new [i][j] to branch on
      //This choice takes us in index order: (1,1), then (1,2), then (1,3), etc., eventually then going to (2,1), making the same choice across a particular row.
      if (BranchIndex[processingnode.row][processingnode.numinrow].jindex+1 < Psnum[BranchIndex[processingnode.row][processingnode.numinrow].iindex])
      {
         newij.jindex = BranchIndex[processingnode.row][processingnode.numinrow].jindex+1;
         newij.iindex = BranchIndex[processingnode.row][processingnode.numinrow].iindex;
      }
      else
      {
         newij.jindex = 0;
         newij.iindex = BranchIndex[processingnode.row][processingnode.numinrow].iindex+1;
      }
      newij.setval = 0;
      BranchIndex[node1.row][node1.numinrow] = newij;
      newij.setval = 1;
      BranchIndex[node2.row][node2.numinrow] = newij;
      
      //Now setup each scenario for the two new nodes
      basemodel->set(GRB_IntParam_ScenarioNumber, 0);
      //Removing all previous bound changes from the model.
      for (int i = 0; i < Pnum; ++i)
      {
         for (int j = 0; j < Psnum[i]; ++j)
         {
            z2[i][j].set(GRB_DoubleAttr_ScenNUB, GRB_UNDEFINED);
            z2[i][j].set(GRB_DoubleAttr_ScenNLB, GRB_UNDEFINED);
         }
      }
      
      //now re-enter values for prior Branches. The fixed variable value is stored with the ij in BranchIndex, so we compute the previous ij
      long int bi_i = node1.row;
      long int bi_j = node1.numinrow;
      while (bi_i >= 1)
      {
         --bi_i;
         if (bi_j %2 == 1)
         {
            --bi_j;
         }
         bi_j /= 2;
         if (BranchIndex[bi_i][bi_j].setval == 1)
         {
            z2[BranchIndex[bi_i][bi_j].iindex][BranchIndex[bi_i][bi_j].jindex].set(GRB_DoubleAttr_ScenNLB, 1);
         }
         else
         {
            z2[BranchIndex[bi_i][bi_j].iindex][BranchIndex[bi_i][bi_j].jindex].set(GRB_DoubleAttr_ScenNUB, 0);
         }
      }
         
      //Repeat all for second scenario
      basemodel->set(GRB_IntParam_ScenarioNumber, 1);
      for (int i = 0; i < Pnum; ++i)
      {
         for (int j = 0; j < Psnum[i]; ++j)
         {
            z2[i][j].set(GRB_DoubleAttr_ScenNUB, GRB_UNDEFINED);
            z2[i][j].set(GRB_DoubleAttr_ScenNLB, GRB_UNDEFINED);
         }
      }
      bi_i = node2.row;
      bi_j = node2.numinrow;
      while (bi_i >= 1)
      {
         --bi_i;
         if (bi_j %2 == 1)
         {
            --bi_j;
         }
         bi_j /= 2;
         if (BranchIndex[bi_i][bi_j].setval == 1)
         {
            z2[BranchIndex[bi_i][bi_j].iindex][BranchIndex[bi_i][bi_j].jindex].set(GRB_DoubleAttr_ScenNLB, 1);
         }
         else
         {
            z2[BranchIndex[bi_i][bi_j].iindex][BranchIndex[bi_i][bi_j].jindex].set(GRB_DoubleAttr_ScenNUB, 0);
         }
      }
      //Removing the processing node from the tree can be done any time after processing is complete and before adding new nodes
      Tree_to_Process.pop_back();

      //Solve the scenarios together
      basemodel->optimize();

      int scenstatus = basemodel->get(GRB_IntAttr_Status);
      if (scenstatus == 2) //Optimal Solution to LP Reached
      {
         node2.bound = basemodel->get(GRB_DoubleAttr_ScenNObjVal);
         std::cout << node2.bound <<std::endl;
            
         if (node2.bound > bestbound)
         {
            //Then check for integrality
            // If integral, new best solution found
            if (check_integral(z2,Pnum,Psnum,tol))
            {
               bestbound = node2.bound;
               for (int i = 0; i < Pnum; ++i)
               {
                  for (int j = 0; j < Psnum[i]; ++j)
                  {
                     Best_Found_Solution[i][j] = z2[i][j].get(GRB_DoubleAttr_X);
                  }
               }
            }
            else // If non-integral, add to tree for branching
            {
               node2.is_leaf = 1;
               Tree_to_Process.push_back(node2);
            }
         }
      }
      else if (scenstatus == 3 or scenstatus == 4) //LP Infeasible
      {
         node2.bound = 0;
         node2.is_leaf = 0;
      }
      else
      {
         std::cout << "Unexpected Scenario Status Reached: Status " << scenstatus << std::endl;
         return -1;
      }

         
      basemodel->set(GRB_IntParam_ScenarioNumber, 0);
      if (basemodel->get(GRB_IntAttr_Status) == 2)
      {
         node1.bound = basemodel->get(GRB_DoubleAttr_ScenNObjVal);
         std::cout << node1.bound <<std::endl;
         
         if (node1.bound > bestbound)
         {
            // If integral, new best solution found
            if (check_integral(z2,Pnum,Psnum,tol))
            {
               bestbound = node1.bound;
               for (int i = 0; i < Pnum; ++i)
               {
                  for (int j = 0; j < Psnum[i]; ++j)
                  {
                     Best_Found_Solution[i][j] = z2[i][j].get(GRB_DoubleAttr_X);
                  }
               }
            }
            else // If non-integral, add to tree for branching
            {
               node1.is_leaf = 1;
               Tree_to_Process.push_back(node1);
            }
         }
      }
      else if (scenstatus == 3 or scenstatus == 4) //LP Infeasible
      {
         node1.bound = 0;
         node1.is_leaf = 0;
      }
      else
      {
         std::cout << "Unexpected Scenario Status Reached: Status " << scenstatus << std::endl;
         return -1;
      }
   }
   
   std::cout << "Best found solution has value: " << bestbound << std::endl;
   for (int i = 0; i < Pnum; ++i)
   {
      for (int j = 0; j < Psnum[i]; ++j)
      {
         std:: cout << "Measure: "<< i << " Point: " << j << " Value: " << Best_Found_Solution[i][j] << std::endl;
      }
   }

   
   //To write model to a file
/*      std::ostringstream filename2;
   filename2 << "/Users/spatterson/ForXcode/PricingLP/Pricing_TestBB1.lp";
   basemodel->write(filename2.str());*/

   
   
//   std::cout << basemodel->get(GRB_IntAttr_Status) <<std::endl;
   
   basemodel->terminate();
   delete basemodel;
   delete env;
   return 0;
}

double objval(const double * y, const double * lambda, std::vector< std::vector< SuppPt> >& x, std::vector<std::vector<double> >& z, const int Pnum, const int * Psnum)
{
   double val = 0;
   int yindex = 0;

   for (int i = 0; i < Pnum-1; ++i)
   {
      for (int j = i+1; j < Pnum; ++j)
      {
         for (int k = 0; k < Psnum[i]; ++k)
         {
            for (int m = 0; m < Psnum[j]; ++m)
            {
               if (z[i][k] != 0 && z[j][m] != 0)
               {
                  val -= lambda[i]*lambda[j]*(x[i][k]*x[i][k] - 2 *( x[i][k]*x[j][m]) + x[j][m]*x[j][m])*z[i][k]*z[j][m];
               }
            }
         }
      }
   }
   for (int i = 0; i < Pnum; ++i)
   {
      for (int k = 0; k < Psnum[i]; ++k)
      {
         val += y[yindex]*z[i][k];
         ++yindex;
      }
   }
   
   return val;
}
double objval2(const double * y, const double * lambda, std::vector< std::vector< SuppPt> >& x, std::vector<std::vector<double> >& z, const int Pnum, const int * Psnum)
{
   double val = 0;
   int yindex = 0;

   for (int i = 0; i < Pnum-1; ++i)
   {
      for (int j = i+1; j < Pnum; ++j)
      {
         for (int k = 0; k < Psnum[i]; ++k)
         {
            for (int m = 0; m < Psnum[j]; ++m)
            {
               if (z[i][k] != 0 && z[j][m] != 0)
               {
                  val += lambda[i]*lambda[j]* 2 *( x[i][k]*x[j][m])*z[i][k]*z[j][m];
               }
            }
         }
      }
   }
   for (int i = 0; i < Pnum-1; ++i)
   {
      for (int j = i+1; j < Pnum; ++j)
      {
         for (int k = 0; k < Psnum[i]; ++k)
         {
            val -= lambda[i]*lambda[j]*(x[i][k]*x[i][k])*z[i][k];
         }
      }
   }
   for (int i = 0; i < Pnum-1; ++i)
   {
      for (int j = i+1; j < Pnum; ++j)
      {
         for (int m = 0; m < Psnum[j]; ++m)
         {
            val -= lambda[i]*lambda[j]*(x[j][m]*x[j][m])*z[j][m];
         }
      }
   }

   for (int i = 0; i < Pnum; ++i)
   {
      for (int k = 0; k < Psnum[i]; ++k)
      {
         val += y[yindex]*z[i][k];
         ++yindex;
      }
   }
   
   return val;
}

double objval3(const double * y, const double * lambda, std::vector< std::vector< SuppPt> >& x, std::vector<std::vector<GRBVar> >& z2, std::vector<std::vector<std::vector<std::vector<GRBVar> > > > & z4, const int Pnum, const int * Psnum)
{
   double val = 0;
   int yindex = 0;

   for (int i = 0; i < Pnum-1; ++i)
   {
      for (int j = i+1; j < Pnum; ++j)
      {
         for (int k = 0; k < Psnum[i]; ++k)
         {
            for (int m = 0; m < Psnum[j]; ++m)
            {
               if (z4[i][k][j][m].get(GRB_DoubleAttr_X) != 0)
               {
                  val += lambda[i]*lambda[j]* 2 *( x[i][k]*x[j][m])*z4[i][k][j][m].get(GRB_DoubleAttr_X);
               }
            }
         }
      }
   }
   for (int i = 0; i < Pnum-1; ++i)
   {
      for (int j = i+1; j < Pnum; ++j)
      {
         for (int k = 0; k < Psnum[i]; ++k)
         {
            val -= lambda[i]*lambda[j]*(x[i][k]*x[i][k])*z2[i][k].get(GRB_DoubleAttr_X);
         }
      }
   }
   for (int i = 0; i < Pnum-1; ++i)
   {
      for (int j = i+1; j < Pnum; ++j)
      {
         for (int m = 0; m < Psnum[j]; ++m)
         {
            val -= lambda[i]*lambda[j]*(x[j][m]*x[j][m])*z2[j][m].get(GRB_DoubleAttr_X);
         }
      }
   }

   for (int i = 0; i < Pnum; ++i)
   {
      for (int k = 0; k < Psnum[i]; ++k)
      {
         val += y[yindex]*z2[i][k].get(GRB_DoubleAttr_X);
         ++yindex;
      }
   }
   
   return val;
}

//Objval and objval2 should always produce the same value. For an integer solution, objval3 also produces the same solution. However, when z2 is not restricted to integer values, neither is z4, and z2[i][k]*z2[j][m] != z4[i][k][j][m]. The values are orders of magnitude different.

bool check_integral(std::vector<std::vector< GRBVar> > &z2, const int Pnum, const int *Psnum, const double &tol)
{
   bool integral = true;
   for (int i = 0; i < Pnum; ++i)
   {
      for (int j = 0; j < Psnum[i]; ++j)
      {
         if (z2[i][j].get(GRB_DoubleAttr_X) <= 1-tol && z2[i][j].get(GRB_DoubleAttr_X) >= tol)
         {
            integral = false;
            return integral;
         }
      }
   }
   return integral;
}


