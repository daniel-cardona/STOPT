#include "RobOptTraj/directTranscription.h"

namespace OCSolver
{

namespace NLP {

namespace derivatives {

namespace analytical {


void obtainAnalyticalGradientInformation(Problem &problem, Alg &algorithm){



}

void computeAnalyticalJacobian(Eigen::VectorXd &x0, Problem &problem, Alg& algorithm,Eigen::VectorXd &nzValues){


    //Obtain the gradient of the time step

    Eigen::VectorXd hkGrad(2);

    hkGrad(0)=-1.0/problem.nSegments;
    hkGrad(1)=1.0/problem.nSegments;

    double t0,tf;

    int nStates=problem.nStates;
    int nControls=problem.nControls;
    int nPath= problem.nPath;

    int eventOffset=problem.nNodes*problem.nStates;
    int pathOffset=(problem.nNodes*problem.nStates)+problem.nEvents;

    utils::getDecVarMatrix(problem,algorithm,x0,problem.decVarMatrix,t0,tf);

    double hk=(tf-t0)/problem.nSegments;

    for(int k=0;k<problem.nNodes;k++){

        //Obtain the states, controls and time in the k node and evaluate the dae

       utils::getVariables(problem,problem.decVarMatrix,problem.kStates,problem.kControls,k);

       double tk=utils::convert_to_original_time(problem.snodes(k),t0,tf);

       problem.dae(problem.kStates,problem.kControls,tk,problem.dx,problem.kPath,problem);

       fGradient(problem.kStates,problem.kControls,tk,problem.dxGradient,problem);

       problem.derivativeMatrix.col(0).segment(k*problem.nStates,problem.nStates)=hkGrad(0)*problem.dx;        //Time grad wrt t0
       problem.derivativeMatrix.col(1).segment(k*problem.nStates,problem.nStates)=hkGrad(1)*problem.dx;        //Time grad wrt tf

       problem.derivativeMatrix.block(k*(nStates),2+(k*(nStates+nControls)),nStates,(nStates+nControls))=hk*problem.dxGradient;

       pathGradient(problem.kStates,problem.kControls,tk,problem.pathGradient,problem);

       problem.derivativeMatrix.block(pathOffset+(k*(nPath)),2+(k*(nStates+nControls)),nPath,(nStates+nControls))=(problem.scale.path.col(k).asDiagonal()*problem.pathGradient);


    }


    //------------------ GRADIENT OF THE EVENT CONSTRAINTS -----------

    Eigen::VectorXd initial_states(problem.nStates);
    Eigen::VectorXd final_states(problem.nStates);

    utils::getVariables(problem,problem.decVarMatrix,initial_states,problem.kControls,0);

    utils::getVariables(problem,problem.decVarMatrix,final_states,problem.kControls,problem.nNodes-1);

    eventGradient(initial_states,final_states,t0,tf,problem.eventGradient_t0,problem.eventGradient_tF,problem);

    problem.derivativeMatrix.block(eventOffset,0,problem.nEvents,1)=problem.scale.events.asDiagonal()*problem.eventGradient_t0.col(0);
    problem.derivativeMatrix.block(eventOffset,1,problem.nEvents,1)=problem.scale.events.asDiagonal()*problem.eventGradient_tF.col(0);

    problem.derivativeMatrix.block(eventOffset,2,problem.nEvents,problem.nStates)=problem.scale.events.asDiagonal()*problem.eventGradient_t0.rightCols(problem.nStates);
    problem.derivativeMatrix.block(eventOffset,problem.nDecVar-(problem.nStates+problem.nControls),problem.nEvents,problem.nStates)=problem.scale.events.asDiagonal()*problem.eventGradient_tF.rightCols(problem.nStates);

    problem.derivativeMatrix(problem.derivativeMatrix.rows()-1,0)=1*problem.scale.time*problem.scale.timeCns;
    problem.derivativeMatrix(problem.derivativeMatrix.rows()-1,1)=-1*problem.scale.time*problem.scale.timeCns;


    //--------------------Computation of the analytical sparse jacobian of the constraints------


    int iRow=0;
    int iCol=0;
    double value=0;

    typedef Eigen::Triplet<double,int> T;
    vector<T> data;
    data.reserve(problem.nnzD);


    for(int i=0;i<problem.nnzD;i++){

        iRow=problem.nzGD(i,0);
        iCol=problem.nzGD(i,1);

        value=problem.derivativeMatrix(iRow,iCol);

        data.push_back(T(iRow,iCol,value));

    }

    problem.J.setFromTriplets(data.begin(),data.end());

    problem.D=problem.sparsity.A+ (problem.sparsity.B2*(problem.J*problem.scale.invScaleMatrix));

    int externalIterator=0;

    nzValues.setZero(problem.nnz);

    for(int k=0;k<problem.D.outerSize();++k)
        for(Eigen::SparseMatrix<double>::InnerIterator it(problem.D,k); it; ++it){

            nzValues(externalIterator)=it.value();
            externalIterator++;
        }

    Eigen::MatrixXd DDense(problem.J.rows(),problem.J.cols());
    DDense=Eigen::MatrixXd(problem.J);




}



}//END analytical namespace

}//END derivatives namespace

}//END NLP namespace

}//END OCSOlver namespace
