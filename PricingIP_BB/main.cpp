//
//  main.cpp
//  PricingIP_BB
//
//  Created by StephanP on 5/16/22.
//  This version is for implementing Branch & Bound by hand, for experimenting with choices like node selection

#include "SuppPt.hpp"
#include "Node.hpp"
#include "/Library/gurobi951/macos_universal2/include/gurobi_c++.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>


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

double objval(const double *, const double *, std::vector< std::vector< SuppPt> >&, std::vector<std::vector<double> >&, const int, const int *);
double objval2(const double *, const double *, std::vector< std::vector< SuppPt> >&, std::vector<std::vector<double> >&, const int, const int *);
double objval3(const double *, const double *, std::vector< std::vector< SuppPt> >&, std::vector<std::vector<GRBVar> >&, std::vector<std::vector<std::vector<std::vector<GRBVar> > > >&, const int, const int *);
bool check_integral(std::vector<std::vector< GRBVar> > &, const int, const int *, const double &);
int clear_scenario(std::vector<std::vector< GRBVar> > &z2, const int Pnum, const int *Psnum);
int input_tree_decisions(std::vector<std::vector< GRBVar> > &, long int , long int , std::vector< std::vector<ij> > &);
int process_LPstatus(const int , Node &, double &, std::vector<std::vector< GRBVar> > &, const int , const int *, const double , std::vector< std::vector<double> > &, std::vector<Node> &, GRBModel* );

int main(int argc, const char * argv[]) {
   
   int Pnum =5;//If NumMoMeas = 1, Number of months to process. Will go consecutively from start month, skipping any empty months.
   std::string startmonth = "Jul";//Three-digit code for month
   int startyear = 15;
   int NumMoMeas = 1;
   int Psnum[Pnum];
   std::vector< std::vector<SuppPt> > PsuppT(Pnum);
   int totalsupp = 0;
   int indices[Pnum];//The use of indices is changing (from previous code) to match the two-indexed vector representation of the support sets
   //Now used primarily for tracking current index in each measure while running all combinations to produce the cost vector c
   double lambda[Pnum];

   std::string temp;

   std::ifstream indata;
   std::ostringstream fname;
   int year;
   std::string month;
   if (argc > 1)
   {
      fname << argv[1];
      indata.open(fname.str());
   }
   else
   {
      std::cout << "Please enter filepath for DenverCrime file." << std::endl;
      std::string filename;
      getline( std::cin, filename);
      indata.open(filename);
   }
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
         std::cout << "Invalid Input. Selected Month/Year not in data." << std::endl;
         return 1;
      }
      if (indata.eof())
      {
         indata.close();
         std::cout << "Invalid Input. End of file reached while searching for starting Month/Year." << std::endl;
         return 1;
      }
   }
   
   double loc1;
   double loc2;
   indata >> loc1 >> loc2;
   std::string currentmonth = startmonth;
         
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
            ++Psnum[i];
            ++totalsupp;
                  
            indata >> year >> month >> loc1 >> loc2;
            if (indata.eof())
            {
               indata.close();
               std::cout << "Invalid Input. End of file reached while reading in months to measures." << std::endl;
               return 1;
            }
         }
         currentmonth = month;
      }
      std::cout << "Month measure: " << i+1 << " Size: " << Psnum[i] << std::endl;
      
      //Scale the masses for each month to sum to 1
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
   
   //Compute the weights for each month
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

   //The ultimate goal is to not compute and store all of the cost vector c, as it is exponential in size
   //In this file, all c's are not necessary
   //Computation is done later as part of generating greedy solution, only place needed
   
   //Reset to the first support point of each measure
   //Next used in generating an initial feasible solution (greedy)
   for (int i = 0; i < Pnum; ++i)
   {
      indices[i] = 0;
   }
   
   //Set up Initial Master Problem. Contains Greedy solution for initial columns
   //As of 6/20/22, the pricing problem is not adding back into this problem for complete column generation
   //Purpose: Need a y from master problem for setting up objective function for IP pricing
   GRBEnv* env = new GRBEnv();
   GRBModel* model = new GRBModel(*env);
   model->set(GRB_IntParam_Method, 0);
   model->set(GRB_IntParam_Presolve, 0);
   model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
//   model->set(GRB_IntParam_OutputFlag, 0); //Turning off display on screen. Disable to see iterations, objective value, etc.
   
   std::vector< GRBVar > w;
   GRBLinExpr  exp[totalsupp];

   //Create a temporary copy of the set of support points, will be modified
   std::vector< std::vector<SuppPt> >::iterator Psuppit = PsuppT.begin();
   std::vector<std::vector<SuppPt> > Psupp2(Pnum);
   std::copy(Psuppit, Psuppit+Pnum, Psupp2.begin() );
   //Using index as a dummy variable for readability
   int index = 0;
   //While the last support point of (any, using first) measure still has positive mass
   while (Psupp2[0][Psnum[0]-1].mass >1e-15)
   {
      //Find smallest remaining mass among support points of current combination
      double minmass = 1;
      for (int i = 0; i < Pnum; ++i)
      {
         if (Psupp2[i][indices[i]].mass < minmass)
         {
            minmass = Psupp2[i][indices[i]].mass;
         }
      }
      std::cout << "Minimum mass is " << minmass <<std::endl;
      
      //Index math based on Algorithm 1 in A Column Generation Approach to the Discrete Barycenter Problem
      long int loocur = S0;
      int unki = 0;
      int index = 0;
      
      double sum1 = 0;
      double sum2 = 0;
      for (int i = 0; i < Pnum; ++i)
      {
         unsigned long int index2 = Psupp2[i][indices[i]].combos-index;
//         std::cout << Psupp2[i][indices[i]].combos << " " << index << std::endl;
//         std::cout << "Index in Pi: " << Psupp2[i][indices[i]].combos - index <<std::endl;

         sum1 += lambda[i]*PsuppT[i][indices[i]].loc1;
         sum2 += lambda[i]*PsuppT[i][indices[i]].loc2;

         loocur /= Psnum[i];
         unki += loocur*index2;
         index += Psnum[i];
         Psupp2[i][indices[i]].mass -= minmass;
      }
/*      for (int i = 0; i < Pnum; ++i)
      {
         std:: cout << Psupp2[i][indices[i]].mass << " ";
      }
      std::cout << std::endl;
      
      std::cout << "current index is " << unki << " out of " << S0-1 <<std::endl;*/
      
      //Create a variable for this combination
      //This is the only place the c's are being used
      double newc = 0;
      SuppPt Pbar0 = SuppPt(sum1, sum2, 0.0);
      for (int i = 0; i < Pnum; ++i)
      {
         newc += lambda[i]*((Pbar0.loc1-PsuppT[i][indices[i]].loc1)*(Pbar0.loc1-PsuppT[i][indices[i]].loc1) +(Pbar0.loc2-PsuppT[i][indices[i]].loc2)*(Pbar0.loc2-PsuppT[i][indices[i]].loc2));
         //Updating indices needed to be done later than previous code projects to allow for cost to be computed here
         if (Psupp2[i][indices[i]].mass <= 1e-15)
         {
            ++indices[i];
         }
      }
//      std::cout << c[unki] << " " << newc <<std::endl;
//      w.push_back(model->addVar(0.0, GRB_INFINITY, c[unki], GRB_CONTINUOUS));
      w.push_back(model->addVar(0.0, GRB_INFINITY, newc, GRB_CONTINUOUS));
      
      //Add variable to corresponding constraints
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
   //This reduced-cost vector yhat is the needed information for the Integer Pricing Program
   double * yhat = model->get(GRB_DoubleAttr_Pi, model->getConstrs(), totalsupp);

   //Begin branch & bound for IP pricing
   //For storing unprocessed nodes during the branch & bound
   //Will sort these as part of experimenting with order of node processing
   std::vector<Node> Tree_to_Process;
   
   GRBModel* basemodel = new GRBModel(*env);
   basemodel->set(GRB_IntParam_Method, 0);
   basemodel->set(GRB_IntAttr_ModelSense, -1);//set to maximize

   //For storing the best integer solution
   std::vector<std::vector<double> > Best_Found_Solution(Pnum);
   double tol = 1e-4;
   
   //This will generate an initial solution from the first point in each measure
   //giving an initial lower bound
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
   
   //Create a variable for each support point
   //z2 are the 2-indexed variables in the integer program
   int Psmax = *std::max_element(Psnum, Psnum+Pnum);
   std::vector< std::vector<GRBVar > > z2(Pnum, std::vector<GRBVar>((Psmax)));
   
   //Create a variable for each product of the two-indexed variables
   //z4 are the four-indexed variables in the integer program
   std::vector< std::vector< std::vector< std::vector< GRBVar >> >>z4(Pnum,std::vector< std::vector< std::vector< GRBVar >>>(Psmax, std::vector< std::vector< GRBVar >>(Pnum, std::vector< GRBVar >(Psmax))));
   
   //Creating base model, which will then be split into two scenarios for each node
   int currentconst = 0;
   int totalconstraints = 0;
   //Add 2-indexed variables to model and the choose-one-per-measure constraints
   for (int i = 0; i < Pnum; ++i)
   {
      GRBLinExpr expp;
      for (int j = 0; j < Psnum[i]; ++j)
      {
         //Names can be removed, also for readability
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
   
   //Now repeat with 4-indexed variables and their inequality constraints
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
               //Names also optional
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
   // add/update the objective function coefficient on z2[i][j]
   std::vector< std::vector< double > > zikcost(Pnum, std::vector<double> ((Psmax),0));
   std::vector< std::vector< double > > zjmcost(Pnum, std::vector<double> ((Psmax),0));
   //Compute the coefficients for z_ik first
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
   //Then compute for z_jm
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
   //Add both to objective function
   currentconst = 0;
   for (int i = 0; i < Pnum; ++i)
   {
      for (int j = 0; j < Psnum[i]; ++j)
      {
         z2[i][j].set(GRB_DoubleAttr_Obj,yhat[currentconst]-zikcost[i][j]-zjmcost[i][j]);
         ++currentconst;
      }
   }
   
//   std::cout << objval(yhat,lambda,PsuppT,ztest,Pnum,Psnum) << " " << objval2(yhat,lambda,PsuppT,ztest,Pnum,Psnum) << " " << objval3(yhat, lambda, PsuppT, z2, z4, Pnum, Psnum) << " " << basemodel->get(GRB_DoubleAttr_ObjVal) << std::endl;
   
   //Solve full linear program (essentially: initial node) and check for integrality; if integral, solution found
   basemodel->optimize();
   
   //Can use the solution to this intial solve for the part of paper talking about how fractional the vertices are
   if (basemodel->get(GRB_IntAttr_Status) == 2)
   {
      if (check_integral(z2,Pnum,Psnum,tol))
      {
         std::cout << "Initial Linear Program produces integral solution. No Branch & Bound required." <<std::endl;
         return 0;
      }
   }
   
   //Branch Index will store the measure number, index within the measure, and the variable's assigned value for each step of the tree.
   //Measure/Index numbers will likely not match the row number and column number within the tree, that is; the two indices of Branch Index.
   //So if you initially decide to branch on the variable that goes with the first support point of measure 3, BranchIndex[0][0] = (2,0), and BranchIndex[1][0] = (2,0,0) and BranchIndex[1][1] = (2,0,1)
   std::vector< std::vector< ij > > BranchIndex(2);//Pnum is the minimum number of rows to reach a feasible solution, but starting with 2 and allocating additional space as needed
   //Experiment: preallocate plenty of space vs. Minimizing space utilized
   //Potential Memory savings: don't allocate entire row of space? Would need to change
   BranchIndex[0].resize(1);
   BranchIndex[0][0]=ij(0,0,-1); //BranchIndex[0] is the initial decision. BranchIndex[1] will contain the two decisions from the first row, etc.
   BranchIndex[1].resize(2);
   
   //Begin setup of 2 scenarios for computing LPs for both child nodes
   //Set Upper bound = 0 in one scenario and Lower bound = 1 in the other
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
   process_LPstatus(basemodel->get(GRB_IntAttr_Status), node2, bestbound, z2, Pnum, Psnum, tol, Best_Found_Solution, Tree_to_Process, basemodel);

   basemodel->set(GRB_IntParam_ScenarioNumber, 0);
   process_LPstatus(basemodel->get(GRB_IntAttr_Status), node1, bestbound, z2, Pnum, Psnum, tol, Best_Found_Solution, Tree_to_Process, basemodel);

   //Now choose new node
   //while Tree vector isn't empty
   while (!Tree_to_Process.empty())
   {
      Node processingnode;
      bool noimprovingnode= false;
      bool finding_next_node = true;
      while (finding_next_node)
      {
         if (Tree_to_Process.empty())
         {
            noimprovingnode = true;
            break;
         }
         processingnode = Tree_to_Process.back();
         //This guarantees we process a node with the potential to improve, and that has leaves.
         //No nodes without leaves should be added to the tree, so the check should be redundant
         if (processingnode.bound <= bestbound or processingnode.is_leaf == 0)
         {
            Tree_to_Process.pop_back();
         }
         else
         {
            finding_next_node = false;
         }
      }
      if (noimprovingnode)
      {
         break;
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
      //This choice takes us in index order: (0,0), then (0,1), then (0,2), etc., eventually then going to (1,0), making the same choice across a particular row.
      //Modify here to choose nodes in a different order, i.e., closest to integer
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
      clear_scenario(z2,Pnum,Psnum);
      
      //now re-enter values for prior Branches. The fixed variable value is stored with the ij in BranchIndex, so we compute the previous ij
      input_tree_decisions(z2,node1.row,node1.numinrow,BranchIndex);
         
      //Repeat all for second scenario
      basemodel->set(GRB_IntParam_ScenarioNumber, 1);
      clear_scenario(z2,Pnum,Psnum);
      input_tree_decisions(z2,node2.row,node2.numinrow,BranchIndex);
      //Removing the processing node from the tree can be done any time after processing is complete and before adding new nodes
      Tree_to_Process.pop_back();

      //Solve the scenarios together
      basemodel->optimize();
      //Process results of current scenario = 1
      process_LPstatus(basemodel->get(GRB_IntAttr_Status), node2, bestbound, z2, Pnum, Psnum, tol, Best_Found_Solution, Tree_to_Process, basemodel);
      //Switch to other scenario and process
      basemodel->set(GRB_IntParam_ScenarioNumber, 0);
      process_LPstatus(basemodel->get(GRB_IntAttr_Status), node1, bestbound, z2, Pnum, Psnum, tol, Best_Found_Solution, Tree_to_Process, basemodel);
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
      
   basemodel->terminate();
   delete basemodel;
   delete env;
   return 0;
}

// Functions for each formulation of computing the objective value from new paper
//Objval and objval2 should always produce the same value.
//For an integer solution, objval3 also produces the same solution. However, when z2 is not restricted to integer values, neither is z4, and z2[i][k]*z2[j][m] != z4[i][k][j][m]. The values are orders of magnitude different.

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

//Function for checking if an LP solution is integral
//function check_integral():
//z2 Vector of Vector of GRB Variables from current scenario
//Pnum number of measures
//Psnum array of ints containing number of points per measure
//Tolerance is set in main; could overload to use default 1e-4
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

//For removing all previous modifications of a specific scenario. Surprisingly, could not find this functionality in Gurobi
//Only need to reset UB/LB since our method does not modify the scenarios in other ways
//z2 Vector of Vector of GRB Variables from current scenario
//Pnum number of measures
//Psnum array of ints containing number of points per measure
int clear_scenario(std::vector<std::vector< GRBVar> > &z2, const int Pnum, const int *Psnum)
{
   for (int i = 0; i < Pnum; ++i)
   {
      for (int j = 0; j < Psnum[i]; ++j)
      {
         z2[i][j].set(GRB_DoubleAttr_ScenNUB, GRB_UNDEFINED);
         z2[i][j].set(GRB_DoubleAttr_ScenNLB, GRB_UNDEFINED);
      }
   }
   return 0;
}

//For inputting the prior decisions in the tree
//z2 Vector of Vector of GRB Variables from current scenario
//Pnum number of measures
//Psnum array of ints containing number of points per measure
//BranchIndex array of ij (struct) containing index, index j of variable z2 and the 0 or 1 decision for z2[i][j] at a particular row/index of tree
int input_tree_decisions(std::vector<std::vector< GRBVar> > &z2, long int bi_i, long int bi_j, std::vector< std::vector<ij> > &BranchIndex)
{
   if (bi_i <= 0)
   {
      std::cout << "No prior rows in tree to enter." <<std::endl;
      return 2;
   }
   if (bi_j <= 0)
   {
      std::cout << "Invalid row index." <<std::endl;
   }
   
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

   return 0;
}

//scenstatus: status code from Gurobi model. 2 is GRB_OPTIMAL, 3 and 4 will be returned for infeasible models, others may be due to time-outs
//node: current node being tested in scenario
//bestbound: transport cost associated with the best found integer solution
//z2 Vector of Vector of GRB Variables from current scenario
//Pnum number of measures
//Psnum array of ints containing number of points per measure
//Tol is small double value used in check_integral
//Best_Found_Solution stores individual 0/1 decision for each variable of the best answer
//Tree_to_Process contains unprocessed nodes
//basemodel is Gurobi model of the IP
int process_LPstatus(const int scenstatus, Node &node, double &bestbound, std::vector<std::vector< GRBVar> > &z2, const int Pnum, const int *Psnum, const double tol, std::vector< std::vector<double> > &Best_Found_Solution, std::vector<Node> &Tree_to_Process, GRBModel* basemodel)
{
   node.bound = basemodel->get(GRB_DoubleAttr_ScenNObjVal);
//   std::cout << node.bound <<std::endl;
         
   if (node.bound > bestbound)
   {
      //Then check for integrality
      // If integral, new best solution found
      if (check_integral(z2,Pnum,Psnum,tol))
      {
         bestbound = node.bound;
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
         node.is_leaf = 1;
         Tree_to_Process.push_back(node);
      }
   }
   else if (scenstatus == 3 or scenstatus == 4) //LP Infeasible or Infeasible/Unbounded, should be infeasible
   {
      node.bound = 0;
      node.is_leaf = 0;
   }
   else
   {
      std::cout << "Unexpected Scenario Status Reached: Status " << scenstatus << std::endl;
      return -1;
   }
   return 0;
}
