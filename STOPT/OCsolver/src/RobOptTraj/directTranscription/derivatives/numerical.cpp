#include "RobOptTraj/directTranscription.h"

namespace OCSolver
{

namespace NLP
{

namespace derivatives
{

namespace numerical {


void scalarGradient(double(*fun)(Eigen::VectorXd &x,Problem &problem,Alg &algorithm),Eigen::VectorXd &x,Problem &problem,Alg &algorithm,Eigen::VectorXd &grad){

    int nDecVar=problem.nDecVar;

    Eigen::VectorXd xlb(nDecVar);
    Eigen::VectorXd xub(nDecVar);

    xlb=problem.xlb;
    xub=problem.xub;


    double sqreps;
    double delj;
    double xs=0.0;

    double F1=0.0;
    double F2=0.0;
    double F3=0.0;

    double dfdx=0.0;

    grad.setZero(nDecVar);


    sqreps=sqrt(MC_EPSILON);


    bool c1,c2;


    c1=(x.array()>=(xlb.array()+sqreps)).any();
    c2=(x.array()<=(xub.array()-sqreps)).any();


        if(c1||c2){

            F3=fun(x,problem,algorithm);
        }

        for(int j=0;j<nDecVar;j++){

            delj=sqreps*(1.0+fabs(x(j)));

            xs=x(j);

            if(xs< xub(j)-delj || (xub(j)==xlb(j))){

                x(j)+=delj;
                F1=fun(x,problem,algorithm);;
            }

            if(xs> xlb(j)+delj || (xub(j)==xlb(j))){

                x(j)=xs-delj;
                F2=fun(x,problem,algorithm);
            }



            if (( (xs< (xub(j)-delj)) && (xs> (xlb(j)+delj)) ) || (xub(j)==xlb(j)) ) {
              // Use central difference formula
              dfdx = ( F1 - F2 )/(2*delj);
            }
            else if (xs>= (xub(j)-delj)) {
              // Variable at upper bound, use backward difference formula
              dfdx = ( F2 - F3 )/(-delj);
            }
            else if (xs<= (xlb(j)+delj)) {
              // Variable at lower bound, use forward difference formula
              dfdx = ( F1 - F3 )/(delj);
            }
            x(j)=xs;

            grad(j)=dfdx;

        }





    }

void jacobianColumn(void(*fun)(Eigen::VectorXd &x,Eigen::VectorXd &g,Problem &problem,Alg& algorithm),Eigen::VectorXd &x, Problem &problem, Alg &algorithm,int jCol,Eigen::VectorXd &jacColumn){

    //Computes only one column of the Jacobian Matrix

    double delj;
    double sqreps;
    double xs;
    int j;


    Eigen::VectorXd xlb(problem.nDecVar);
    Eigen::VectorXd xub(problem.nDecVar);

    xlb=problem.xlb;
    xub=problem.xub;


    Eigen::VectorXd F1 (problem.nCns);
    Eigen::VectorXd F2 (problem.nCns);
    Eigen::VectorXd F3 (problem.nCns);

    Eigen::VectorXd dfdx_j(problem.nCns);

    F1.setZero();
    F2.setZero();
    F3.setZero();

    jacColumn.setZero(problem.nCns);


    sqreps=sqrt(MC_EPSILON);

    bool c1,c2;


    c1=(x.array()>(xlb.array()+sqreps)).any();
    c2=(x.array()<(xub.array()-sqreps)).any();

    if(c1||c2){

        fun(x,F3,problem,algorithm);
    }

    j=jCol;

    delj=sqreps*(1+fabs(x(j)));


    xs=x(j);

    //Central difference formula
    if((xs<xub(j)-delj && xs>xlb(j)+delj) || (xub(j)==xlb(j))){
            x(j)+=delj;
            fun(x,F1,problem,algorithm);
            x(j)=xs-delj;
            fun(x,F2,problem,algorithm);
            dfdx_j=(F1-F2)/(2*delj);

            x(j)=xs;
    }

    //Backward difference formula (variable at upper bound)
    else if(xs>=xub(j)-delj)
    {
        x(j)=xs-delj;
        fun(x,F1,problem,algorithm);
        x(j)=xs;

        dfdx_j=(F1-F3)/(-delj);
    }

    //Forward difference formula (variable at lower bound)

    else if(xs<=xlb(j)+delj){
        x(j)=xs+delj;
        fun(x,F1,problem,algorithm);
        x(j)=xs;

        dfdx_j=(F1-F3)/delj;

    }

    jacColumn=dfdx_j;

}

void computeNumericalJacobian(void(*fun)(Eigen::VectorXd &x,Eigen::VectorXd &g,Problem &problem,Alg& algorithm),Eigen::VectorXd &x,Problem &problem,Alg& algorithm, Eigen::VectorXd &nzValues){

    /* This function uses the method of Curtis, Powell and Reid (1974) to
     * evaluate efficiently the sparse Jacobian by perturbing simultaneously groups of variables.
     * Reference:
     * A. R. Curtis, M.J.D. Powell and J.K. Reid
     * "On the estimation of Sparse Jacobian Matrices"
     * J Inst Maths Applics (1974) 13, 117-119
     *
     */

    double delj;
    double sqreps;

    int nDecVar=problem.nDecVar;
    int nCns=problem.nStates*(problem.nNodes)+problem.nEvents+(problem.nPath*problem.nNodes)+1;

    int iRow=0;
    int iCol=0;

    int element=0;

    Eigen::SparseMatrix<double, Eigen::RowMajor> J(nCns,nDecVar);
    typedef Eigen::Triplet<double,int> T;
    vector<T> data;

    //---------------

    Eigen::VectorXd v(problem.nnzD);

    //---------------

    //Variables of groups

    int nGroups=problem.sizeGroupsD.size();

    Eigen::VectorXd F1(nCns);
    Eigen::VectorXd F2(nCns);

    Eigen::VectorXd xp(nDecVar);

    double valueNZ=0;

    nzValues.setZero(problem.nnz);

    sqreps=sqrt(MC_EPSILON);

    delj=sqreps;

    //delj=1e-8;

    //Iterate along the groups

    for(int i=0;i<nGroups;i++){

            xp=x;

            //Iterate along the non-zero elements that doesn't generate zeros in the jacobian

            xp+=problem.sparsity.positivePerturbationMatrix.col(i);

            fun(xp,F1,problem,algorithm);

            //Perturbate again the decisión vector

            xp-=problem.sparsity.negativePerturbationMatrix.col(i);

            fun(xp,F2,problem,algorithm);

            //Iterate to set the jacobian just in the non-zero elements

            for(int j=0;j<problem.sizeGroupsD(i);j++){


                    for(int k=0;k<problem.nnzD;k++){

                            iRow=problem.nzGD(k,0); //Row of the k th non zero element
                            iCol=problem.nzGD(k,1); //Column of the k th non zero element


                            if(iCol==problem.idxGroupsD(i,j)){

                                    valueNZ=(F1(iRow)-F2(iRow))/(2*delj);

                                    //-----------
                                    v(k)=valueNZ;
                                    //-----------

                                    data.push_back(T(iRow,iCol,valueNZ));

                            }
                    }

            }

    }

    J.setFromTriplets(data.begin(), data.end());

    Eigen::SparseMatrix<double> D(problem.nCns,nDecVar);

    D=problem.sparsity.A+ (problem.sparsity.B2*J);

    int externalIterator=0;

    for(int k=0;k<D.outerSize();++k)
        for(Eigen::SparseMatrix<double>::InnerIterator it(D,k); it; ++it){

            nzValues(externalIterator)=it.value();
            externalIterator++;
        }




    Eigen::MatrixXd DDense(J.rows(),J.cols());
    DDense=Eigen::MatrixXd(J);

    //--------------------------------------------------------
    std::ofstream myfile;
    myfile.open ("numerical.txt");
    myfile  << std::scientific << std::setprecision(20) << DDense;
    //--------------------------------------------------------


}

//%%%%%%%%%%%%%%%%%%%%%%%%%% DEPRECATED FUNCTIONS V2

void computeNumericalJacobianv2(void(*fun)(Eigen::VectorXd &x,Eigen::VectorXd &g,Problem &problem,Alg& algorithm),Eigen::VectorXd &x,Problem &problem,Alg& algorithm, Eigen::VectorXd &nzValues){


    /* This function uses the method of Curtis, Powell and Reid (1974) to
     * evaluate efficiently the sparse Jacobian by perturbing simultaneously groups of variables.
     * Reference:
     * A. R. Curtis, M.J.D. Powell and J.K. Reid
     * "On the estimation of Sparse Jacobian Matrices"
     * J Inst Maths Applics (1974) 13, 117-119
     *
     */

    double delj;
    double sqreps;

    int nDecVar=problem.nDecVar;
    int nCns=problem.nCns;

    int iRow=0;
    int iCol=0;

    int element=0;

    //Variables of groups

    int nGroups=problem.sizeGroups.size();

    Eigen::VectorXd F1(nCns);
    Eigen::VectorXd F2(nCns);

    Eigen::VectorXd xp(nDecVar);

    nzValues.setZero(problem.nnz);

    sqreps=sqrt(MC_EPSILON);

    delj=sqreps;

    //Iterate along the groups

    for(int i=0;i<nGroups;i++){

            xp=x;

            //Iterate along the non-zero elements that doesn't generate zeros in the jacobian

            for(int j=0;j<problem.sizeGroups(i);j++){

                    delj=sqreps;

                    element=problem.idxGroups(i,j);

                    xp(element)+=delj;

            }

            fun(xp,F1,problem,algorithm);

            //Perturbate again the decisión vector

            for(int j=0;j<problem.sizeGroups(i);j++){

                    delj=sqreps;

                    element=problem.idxGroups(i,j);

                    xp(element)-=2*delj;

            }

            fun(xp,F2,problem,algorithm);

            //Iterate to set the jacobian just in the non-zero elements

            for(int j=0;j<problem.sizeGroups(i);j++){


                    for(int k=0;k<problem.nnz;k++){

                            iRow=problem.nzG(k,0); //Row of the k th non zero element
                            iCol=problem.nzG(k,1); //Column of the k th non zero element


                            if(iCol==problem.idxGroups(i,j)){

                                    nzValues(k)=(F1(iRow)-F2(iRow))/(2*delj);


                            }
                    }

            }

    }







}


//%%%%%%%%%%%%%%%%%%%%%%%%%% VERSION 3 NUMERICAL DERIVATIVES FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%

void dxJacColumn(Eigen::VectorXd &x,int jCol,Eigen::VectorXd &jacColumn,Problem &problem){

    double delj;
    double sqreps;
    double xs;

    int nCns=problem.nStates;

    Eigen::VectorXd F1(nCns);
    Eigen::VectorXd F2(nCns);
    Eigen::VectorXd F3(nCns);

    Eigen::VectorXd dfdx_j(nCns);

    Eigen::VectorXd states(problem.nStates);
    Eigen::VectorXd controls(problem.nControls);

    Eigen::VectorXd dummy_path(problem.nPath);

    states=x.head(problem.nStates);
    controls=x.tail(problem.nControls);

    Eigen::VectorXd xlb(problem.nStates+problem.nControls);
    Eigen::VectorXd xub(problem.nStates+problem.nControls);

    xlb<<problem.bounds.states.lower, problem.bounds.controls.lower;
    xub<<problem.bounds.states.upper, problem.bounds.controls.upper;

    double t0=0.0;
    double tf=1.0;

    F1.setZero();
    F2.setZero();
    F3.setZero();

    jacColumn.setZero(nCns);

    sqreps=sqrt(MC_EPSILON);

    bool c1,c2;

    c1=(x.array()>(xlb.array()+sqreps)).any();
    c2=(x.array()<(xub.array()-sqreps)).any();

    if(c1||c2){

        problem.dae(states,controls,t0,F3,dummy_path,problem);
    }

    int j=jCol;

    delj=sqreps*(1+fabs(x(j)));


    xs=x(j);

    //Central difference formula
    if((xs<xub(j)-delj && xs>xlb(j)+delj) || (xub(j)==xlb(j))){

            x(j)+=delj;
            states=x.head(problem.nStates);
            controls=x.tail(problem.nControls);

            problem.dae(states,controls,t0,F1,dummy_path,problem);

            x(j)=xs-delj;
            states=x.head(problem.nStates);
            controls=x.tail(problem.nControls);

            problem.dae(states,controls,t0,F2,dummy_path,problem);

            dfdx_j=(F1-F2)/(2*delj);

            x(j)=xs;
    }

    //Backward difference formula (variable at upper bound)
    else if(xs>=xub(j)-delj)
    {
        x(j)=xs-delj;
        states=x.head(problem.nStates);
        controls=x.tail(problem.nControls);

        problem.dae(states,controls,t0,F1,dummy_path,problem);

        x(j)=xs;

        dfdx_j=(F1-F3)/(-delj);
    }

    //Forward difference formula (variable at lower bound)

    else if(xs<=xlb(j)+delj){

        x(j)=xs+delj;
        states=x.head(problem.nStates);
        controls=x.tail(problem.nControls);

        problem.dae(states,controls,t0,F1,dummy_path,problem);
        x(j)=xs;

        dfdx_j=(F1-F3)/delj;

    }

    jacColumn=dfdx_j;
}

void pathJacColumn(Eigen::VectorXd &x,int jCol,Eigen::VectorXd &jacColumn,Problem &problem){


    double delj;
    double sqreps;
    double xs;

    int nCns=problem.nPath;

    Eigen::VectorXd F1(nCns);
    Eigen::VectorXd F2(nCns);
    Eigen::VectorXd F3(nCns);

    Eigen::VectorXd dfdx_j(nCns);

    Eigen::VectorXd states(problem.nStates);
    Eigen::VectorXd controls(problem.nControls);

    Eigen::VectorXd dummy_dx(problem.nStates);

    states=x.head(problem.nStates);
    controls=x.tail(problem.nControls);

    Eigen::VectorXd xlb(problem.nStates+problem.nControls);
    Eigen::VectorXd xub(problem.nStates+problem.nControls);

    xlb<<problem.bounds.states.lower, problem.bounds.controls.lower;
    xub<<problem.bounds.states.upper, problem.bounds.controls.upper;

    double t0=0.0;
    double tf=1.0;

    F1.setZero();
    F2.setZero();
    F3.setZero();

    jacColumn.setZero(nCns);

    sqreps=sqrt(MC_EPSILON);

    bool c1,c2;

    c1=(x.array()>(xlb.array()+sqreps)).any();
    c2=(x.array()<(xub.array()-sqreps)).any();

    if(c1||c2){

        problem.dae(states,controls,t0,dummy_dx,F3,problem);
    }

    int j=jCol;

    delj=sqreps*(1+fabs(x(j)));


    xs=x(j);

    //Central difference formula
    if((xs<xub(j)-delj && xs>xlb(j)+delj) || (xub(j)==xlb(j))){

            x(j)+=delj;
            states=x.head(problem.nStates);
            controls=x.tail(problem.nControls);

            problem.dae(states,controls,t0,dummy_dx,F1,problem);

            x(j)=xs-delj;
            states=x.head(problem.nStates);
            controls=x.tail(problem.nControls);

            problem.dae(states,controls,t0,dummy_dx,F2,problem);

            dfdx_j=(F1-F2)/(2*delj);

            x(j)=xs;
    }

    //Backward difference formula (variable at upper bound)
    else if(xs>=xub(j)-delj)
    {
        x(j)=xs-delj;
        states=x.head(problem.nStates);
        controls=x.tail(problem.nControls);

        problem.dae(states,controls,t0,dummy_dx,F1,problem);

        x(j)=xs;

        dfdx_j=(F1-F3)/(-delj);
    }

    //Forward difference formula (variable at lower bound)

    else if(xs<=xlb(j)+delj){

        x(j)=xs+delj;
        states=x.head(problem.nStates);
        controls=x.tail(problem.nControls);

        problem.dae(states,controls,t0,dummy_dx,F1,problem);
        x(j)=xs;

        dfdx_j=(F1-F3)/delj;

    }

    jacColumn=dfdx_j;
}

void eventJacColumn_t0(Eigen::VectorXd &x,int jCol,Eigen::VectorXd &jacColumn,Problem &problem){


    double delj;
    double sqreps;
    double xs;

    int nCns=problem.nEvents;

    Eigen::VectorXd F1(nCns);
    Eigen::VectorXd F2(nCns);
    Eigen::VectorXd F3(nCns);

    Eigen::VectorXd dfdx_j(nCns);

    Eigen::VectorXd states(problem.nStates);

    Eigen::VectorXd dummy_states_tF(problem.nStates);
    double dummy_tf=0;

    double t0=x(0);

    states=x.tail(problem.nStates);

    Eigen::VectorXd xlb(problem.nStates+1);
    Eigen::VectorXd xub(problem.nStates+1);

    xlb<<problem.bounds.startTime.lower,problem.bounds.states.lower;
    xub<<problem.bounds.startTime.upper,problem.bounds.states.upper;

    F1.setZero();
    F2.setZero();
    F3.setZero();

    jacColumn.setZero(nCns);

    sqreps=sqrt(MC_EPSILON);

    bool c1,c2;

    c1=(x.array()>(xlb.array()+sqreps)).any();
    c2=(x.array()<(xub.array()-sqreps)).any();

    if(c1||c2){

        problem.events(states,dummy_states_tF,t0,dummy_tf,F3,problem);
    }

    int j=jCol;

    delj=sqreps*(1+fabs(x(j)));


    xs=x(j);

    //Central difference formula
    if((xs<xub(j)-delj && xs>xlb(j)+delj) || (xub(j)==xlb(j))){

            x(j)+=delj;

            t0=x(0);
            states=x.tail(problem.nStates);

            problem.events(states,dummy_states_tF,t0,dummy_tf,F1,problem);

            x(j)=xs-delj;

            t0=x(0);
            states=x.tail(problem.nStates);

            problem.events(states,dummy_states_tF,t0,dummy_tf,F2,problem);

            dfdx_j=(F1-F2)/(2*delj);

            x(j)=xs;
    }

    //Backward difference formula (variable at upper bound)
    else if(xs>=xub(j)-delj)
    {
        x(j)=xs-delj;

        t0=x(0);
        states=x.tail(problem.nStates);

        problem.events(states,dummy_states_tF,t0,dummy_tf,F1,problem);

        x(j)=xs;

        dfdx_j=(F1-F3)/(-delj);
    }

    //Forward difference formula (variable at lower bound)

    else if(xs<=xlb(j)+delj){

        x(j)=xs+delj;

        t0=x(0);
        states=x.tail(problem.nStates);

        problem.events(states,dummy_states_tF,t0,dummy_tf,F1,problem);

        x(j)=xs;

        dfdx_j=(F1-F3)/delj;

    }

    jacColumn=dfdx_j;

}

void eventJacColumn_tf(Eigen::VectorXd &x,int jCol,Eigen::VectorXd &jacColumn,Problem &problem){

    double delj;
    double sqreps;
    double xs;

    int nCns=problem.nEvents;

    Eigen::VectorXd F1(nCns);
    Eigen::VectorXd F2(nCns);
    Eigen::VectorXd F3(nCns);

    Eigen::VectorXd dfdx_j(nCns);

    Eigen::VectorXd states(problem.nStates);

    Eigen::VectorXd dummy_states_t0(problem.nStates);
    double dummy_t0=0;

    double tf=x(0);

    states=x.tail(problem.nStates);

    Eigen::VectorXd xlb(problem.nStates+1);
    Eigen::VectorXd xub(problem.nStates+1);

    xlb<<problem.bounds.finalTime.lower,problem.bounds.states.lower;
    xub<<problem.bounds.finalTime.upper,problem.bounds.states.upper;

    F1.setZero();
    F2.setZero();
    F3.setZero();

    jacColumn.setZero(nCns);

    sqreps=sqrt(MC_EPSILON);

    bool c1,c2;

    c1=(x.array()>(xlb.array()+sqreps)).any();
    c2=(x.array()<(xub.array()-sqreps)).any();

    if(c1||c2){

        problem.events(dummy_states_t0,states,dummy_t0,tf,F3,problem);
    }

    int j=jCol;

    delj=sqreps*(1+fabs(x(j)));


    xs=x(j);

    //Central difference formula
    if((xs<xub(j)-delj && xs>xlb(j)+delj) || (xub(j)==xlb(j))){

            x(j)+=delj;

            tf=x(0);
            states=x.tail(problem.nStates);

            problem.events(dummy_states_t0,states,dummy_t0,tf,F1,problem);

            x(j)=xs-delj;

            tf=x(0);
            states=x.tail(problem.nStates);

             problem.events(dummy_states_t0,states,dummy_t0,tf,F2,problem);

            dfdx_j=(F1-F2)/(2*delj);

            x(j)=xs;
    }

    //Backward difference formula (variable at upper bound)
    else if(xs>=xub(j)-delj)
    {
        x(j)=xs-delj;

        tf=x(0);
        states=x.tail(problem.nStates);

         problem.events(dummy_states_t0,states,dummy_t0,tf,F1,problem);

        x(j)=xs;

        dfdx_j=(F1-F3)/(-delj);
    }

    //Forward difference formula (variable at lower bound)

    else if(xs<=xlb(j)+delj){

        x(j)=xs+delj;

        tf=x(0);
        states=x.tail(problem.nStates);

         problem.events(dummy_states_t0,states,dummy_t0,tf,F1,problem);

        x(j)=xs;

        dfdx_j=(F1-F3)/delj;

    }

    jacColumn=dfdx_j;


}



}//END numerical namespace

}//END derivatives namespace

}//END NLP namespace

}//END OCSolver namespace
