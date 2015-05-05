#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include <opengm/opengm.hxx>
#include <opengm/datastructures/marray/marray.hxx>
#include <opengm/functions/potts.hxx>
#include <opengm/graphicalmodel/space/discretespace.hxx>
#include <opengm/graphicalmodel/graphicalmodel.hxx>
#ifdef WITH_HDF5
#include <opengm/graphicalmodel/graphicalmodel_hdf5.hxx>
#endif
#include <opengm/operations/adder.hxx>
#include <opengm/operations/maximizer.hxx>
#include <opengm/operations/minimizer.hxx>
//#include <opengm/operations/argmin.hxx>
#include <opengm/inference/astar.hxx>
//#include <opengm/inference/messagepassing.hxx>
#include <opengm/inference/messagepassing/messagepassing.hxx>
#include <opengm/inference/messagepassing/messagepassing_bp.hxx>

using namespace std; // 'using' is used only in example code
//using namespace opengm; 

template<class T>
void createAndPrintData(size_t nrOfWorkers, size_t nrOfItems, marray::Marray<T>& data, vector<int> &actual_labels) {
   size_t shape[]={nrOfWorkers, nrOfItems};
   data.resize(shape, shape+2);
   cout << "pariwise costs:" << endl;
   srand(0);
   double workerProbability[2000]; // probability of honesty for each worker, sampled from uniform distribution   
   cout << "Worker Probabilities = " << endl; 
   for(size_t v = 0; v < nrOfWorkers; v++) { 
      workerProbability[v] = 1; // static_cast<float>(rand() % 100) * 0.01;
      cout << "workerProbabilities[" << v << "] = "<< setw(6) << setprecision(2) << workerProbability[v] << endl;  
   } 

   size_t taskLabels[2000]; // task labels: binary 0 or 1  
   cout << "Task labels = " << endl; 
   for(size_t s = 0; s < nrOfItems; s++) { 
      taskLabels[s] = rand() % 2;
      cout << "taskLabels["  << s << "] = " << setw(6) << setprecision(2) <<  taskLabels[s] << endl;
      actual_labels[s] = taskLabels[s];  
   } 	

   for(size_t v=0; v<nrOfWorkers; v++) {
      for(size_t s=0; s<nrOfItems; s++) {
         double tempProbability = static_cast<float>(rand() % 100) * 0.01;
         if( workerProbability[v] > tempProbability ) 
            data(v, s) = taskLabels[s]; 
         else 	
            data(v, s) = 1-taskLabels[s]; 
         cout << left << setw(6) << setprecision(2) << data(v, s);
      }
      cout << endl;
   }
}


void printSolution(const vector<size_t>& solution, size_t discretizationSize, size_t nrOfItems, size_t nrOfWorkers, vector<int> actual_labels) {
   set<size_t> unique;
   cout << endl << "Solution Labels :" << endl;
   for(size_t v=0;v<solution.size();++v) {
      cout << left << setw(2) << v << "  ->   " << solution[v] << endl;
   }

   // worker confidences : 
   for( int i = 0; i < nrOfWorkers; i++) 
      cout << "Worker Confidence = " << solution[i] * 1.0 / discretizationSize << endl; 

   int correctNum = 0; 
   for( int i = 0; i < nrOfItems; i++ ) { 
      if( solution[i+nrOfWorkers] == actual_labels[i] )
         correctNum++; 
      cout << "Item label = " <<  solution[i+nrOfWorkers] << endl; 
   }
   cout << " >>>>>>>>>>>>>>>>>>> Accuracy on labels = " << setprecision(4) << correctNum * 100.0 / nrOfItems << " percent " << endl; 
}

int main() {
   // model parameters
   const int discretizationSize = 10; 
   const size_t nrOfWorkers = 5;
   const size_t nrOfItems = 5; 
   //const size_t nrOfLabels = discretizationSize*2;
   cout << endl << "Matching with one to one correspondences:" << endl
        << nrOfWorkers << " variables with " 
        << nrOfItems <<" items" << endl << endl;
	
   // pairwise costs
   marray::Marray<double> data;
   vector<int> actual_labels(5);  
   cout << ">>>>>>>>>>>>>>>>>>> About to clear data" << endl; 
   createAndPrintData(nrOfWorkers, nrOfItems, data, actual_labels);
	
   // build the model with
   // - nrOfworkers variables 
   // - Each potential function: f(workerLabel, workerHonesty) = workerHonesty I { workerLabel == data(worker,item) } + (1 - workerHonesty) I {workerLabel != data(worker,item)} 
   typedef opengm::ExplicitFunction<double> ExplicitFunction;
   typedef opengm::GraphicalModel<
      double, 
      opengm::Adder,
      opengm::ExplicitFunction<float>,  
      opengm::DiscreteSpace<>
      > Model;  

   typedef Model::FunctionIdentifier FunctionIdentifier;
	
   opengm::DiscreteSpace<> space; 
   // worker variables 
   for( int workerIter = 0; workerIter < nrOfWorkers; workerIter++ ) 
      space.addVariable(discretizationSize); 
   // item variables 
   for( int itemIter = 0; itemIter < nrOfItems; itemIter++ ) 
      space.addVariable(2); 
	
   cout << " number of variables = " << space.numberOfVariables() << endl; 

   // define the graphical model 	
   Model gm(space);

   // add 1st order functions and factors
   for(size_t workerIter=0; workerIter < nrOfWorkers; workerIter++) { // loop through the workers 
      for(size_t itemIter = 0; itemIter < nrOfItems; itemIter++) {
         // define an explicit function 
         size_t shape[] = {2, discretizationSize};
         opengm::ExplicitFunction<float> f(shape, shape + 2, 1.0); 
         for( size_t itemLabel = 0; itemLabel < 2; itemLabel++ ) { 
            for(  size_t workerConfidence = 0; workerConfidence < discretizationSize; workerConfidence++ ) { 
               if( data(workerIter, itemIter) == itemLabel ) 
                  f(itemLabel, workerConfidence) = log(1-workerConfidence * 1.0 / discretizationSize + 1e-8); 	
               else if( data(workerIter, itemIter) != itemLabel ) 
                  f(itemLabel, workerConfidence) =   log(workerConfidence * 1.0 / discretizationSize + 1e-8); 
            }
         }
         Model::FunctionIdentifier fid = gm.addSharedFunction(f); 
         size_t vi[]={ workerIter, nrOfWorkers + itemIter }; 
         gm.addFactor(fid,vi,vi+2); 	
      }
   }

#ifdef WITH_HDF5
     cout << ">> Save model to model.h5: " << endl;
     opengm::hdf5::save(gm,"model.h5","gm");
#endif
   
   // set up the optimizer (A-star search)
   typedef opengm::AStar<Model, opengm::Minimizer> AstarType;
   AstarType astar(gm);
   // obtain and print the argmin
   AstarType::VerboseVisitorType verboseVisitor;
   cout << "\nA-star search:\n";
   astar.infer(verboseVisitor);
   vector<size_t> solution;
   astar.arg(solution);

/*     
   // set up the optimizer (loopy belief propagation) 
   typedef opengm::BeliefPropagationUpdateRules<Model, opengm::Minimizer> UpdateRules;
   typedef opengm::MessagePassing<Model, opengm::Minimizer, UpdateRules> BeliefPropagation; 
   const size_t maxNumberOfIterations = 100; 
   const double convergenceBound = 1.0/1000000; 
   const double damping = 0.0;
   BeliefPropagation::Parameter parameter(maxNumberOfIterations, convergenceBound, damping); 
   BeliefPropagation bp(gm, parameter);
   // optimize (approximately) 11
   BeliefPropagation::VerboseVisitorType visitor; 
   bp.infer(visitor); 
   // obtain the (approximate) ARGMax 14
   vector<size_t> solution(space.numberOfVariables()); 
   bp.arg(solution); 
*/ 

  cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>> Printing the solution : " << endl; 
  printSolution(solution, discretizationSize, nrOfItems, nrOfWorkers, actual_labels); 

}
