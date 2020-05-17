/*
 *
 * Copyright (C) 2019
 * Daniel S. Cardona <ingdanielcardona@gmail.com>, Gustavo Arechavaleta <garechav@cinvestav.edu.mx>
 * CINVESTAV - Saltillo Campus
 *
 * This file is part of OpenHRC
 * OpenHRC is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * OpenHRC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.

 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301  USA
 */

/**
 *	\file src/problem.h
 *	\author Daniel S. Cardona
 *	\version 1.0
 *	\date 2019
 *
 * Problem and algorithm class implementation.
 */

#include "RobOptTraj/directTranscription.h"

namespace OCSolver
{

namespace NLP
{

namespace sparsity {

//%%%%%%%%%%%%%%%%%%%% FULL JACOBIAN FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%

void detectFullJacobianSparsity(void(*fun)(Eigen::VectorXd &x,Eigen::VectorXd &g,Problem &problem,Alg& algorithm),Eigen::VectorXd &x,Problem &problem,Alg& algorithm){


    cout<<endl;
    cout<<"-------> Detecting Full Jacobian Sparsity <------"<<endl;
    cout<<endl;

    int nDecVar=problem.nDecVar;
    int nCns=problem.nCns;

    int nzcount_G=0; //Non-zero non-constant elements (counter)

    double s=1.0e6*sqrt(MC_EPSILON);

    double tol= 1.e-16*pow(MC_EPSILON,0.8)*MAX(1.0,x.norm());

    Eigen::VectorXd jacCol1(nCns);
    Eigen::VectorXd jacCol2(nCns);
    Eigen::VectorXd jacCol3(nCns);

    Eigen::VectorXd xlb(nDecVar);
    Eigen::VectorXd xub(nDecVar);
    Eigen::VectorXd xp(nDecVar);

    Eigen::VectorXd s_v(nDecVar);

    Eigen::MatrixXd Gidx(nCns*nDecVar,2);

    Eigen::MatrixXd jacobian(nCns,nDecVar);


    xlb=problem.xlb;
    xub=problem.xub;

    s_v.setOnes();

    s_v*=s;

    //Iterate along the columns of the jacobian matrix
    for(int j=0;j<nDecVar;j++){

        xp=x;

        utils::clip_vector_given_bounds(xp,xlb,xub);

        derivatives::numerical::jacobianColumn(fun,xp,problem,algorithm,j,jacCol1);

        //Perturbate the decision vector

        xp=x.array()+(0.1*x.array().abs())+s_v.array();

        utils::clip_vector_given_bounds(xp,xlb,xub);

        derivatives::numerical::jacobianColumn(fun,xp,problem,algorithm,j,jacCol2);

        //Perturbate the decision vector(Again...)

        xp=x.array()-(0.15*x.array().abs())-(1.1*s_v.array());

        utils::clip_vector_given_bounds(xp,xlb,xub);

        derivatives::numerical::jacobianColumn(fun,xp,problem,algorithm,j,jacCol3);

        //Now iterate alog the constraints (rows)

        for(int i=0;i<nCns;i++){

            if( (fabs(jacCol1(i)) +  fabs(jacCol2(i)) + fabs(jacCol3(i)) ) >= tol){

                    Gidx(nzcount_G,0)=i;//Row
                    Gidx(nzcount_G,1)=j;//Column

                    nzcount_G++;

                }
            }

        }



    //Now save the sparsity structure

    problem.nnz=nzcount_G;

    problem.nzG.setZero(nzcount_G,2);

    problem.nzG=Gidx.topRows(nzcount_G);



    cout<<"-----Full Jacobian sparsity detected numerically-----"<<endl;
    cout<<endl;
    cout<<"+Nonzero elements             ->"<<problem.nnz<<endl;

}

void getIndexGroups(Problem &problem){

    /* This function uses the method of Curtis, Powell and Reid (1974) to find groups of variables
     * to evaluate efficiently the sparse Jacobian by perturbing simultaneously groups of variables.
     * Reference:
     * A. R. Curtis, M.J.D. Powell and J.K. Reid
     * "On the estimation of Sparse Jacobian Matrices"
     * J Inst Maths Applics (1974) 13, 117-119
     *
     */

    int nDecVar= problem.nDecVar;
    int nCns=   problem.nCns;

    Eigen::MatrixXd nnzG(problem.nnz,2);

    nnzG=problem.nzG;

    //Define a matrix with 1 in the non-zero elements

    Eigen::SparseMatrix<double, Eigen::RowMajor> J(nCns,nDecVar);
    typedef Eigen::Triplet<double,int> T;
    vector<T> entries;

      for(int k=0; k<problem.nnz;++k){
        entries.push_back( T(nnzG(k,0), nnzG(k,1), 1) );
      }

      J.setFromTriplets(entries.begin(), entries.end());

      Eigen::MatrixXd JDense(nCns,nDecVar);

      JDense=Eigen::MatrixXd(J);

      //********************************************
      std::ofstream file2("fulljacobianSparsity.txt");
        if (file2.is_open())
        {
          file2 << JDense << '\n';
        }
      //********************************************

      //cout<<JDense<<endl;

      Eigen::VectorXd C1(nCns);
      Eigen::VectorXd C2(nCns);

      Eigen::MatrixXd idxGroups(nDecVar,nDecVar);
      Eigen::VectorXd col_check(nDecVar);
      Eigen::VectorXd sizeGroups(nDecVar);


      col_check.setZero();
      idxGroups.setZero();
      sizeGroups.setZero();

        //Define the first group

      idxGroups(0,0)=0; //(Group,Index)=Column index

      col_check(0)=1; //The column 0 is already in a group...


      //Iterate along the columns of the jacobian matrix

      int gcount=1;
      int colcount=1;
      bool ok;
      double dotCols;


      for(int j=1;j<nDecVar;j++){

            ok=true;

            for(int l=0;l<gcount;l++){

                        if(j==idxGroups(0,l)){
                            ok=false;
                            break;
                        }

                         C1=JDense.col(idxGroups(0,l));
                         C2=JDense.col(j);

                         dotCols=C1.dot(C2);

                         if(dotCols>0.0){
                             ok=false;
                         }
            }

            if(ok){
                idxGroups(0,gcount)=j;
                gcount++;
                colcount++;
                col_check(j)=1;
            }
     }

      sizeGroups(0)=gcount;

      //Now lets obtain the other groups


    int group_index=0;

    int maxSize=0;

    while(colcount<nDecVar){

        group_index++;

        if(gcount>maxSize)
            maxSize=gcount;

        gcount=0;

        for(int j=1;j<nDecVar;j++){

                    ok=true;

                    if(col_check(j)==1){

                            ok=false;
                     }


                    //If the column is not in a group use it for a group

                    if(ok==true){

                          for(int l=0;l<gcount;l++){

                                    C1=JDense.col(idxGroups(group_index,l));
                                    C2=JDense.col(j);

                                    dotCols=C1.dot(C2);

                                     if(dotCols>0.0){

                                            ok=false;
                                      }

                            }
                     }

                    //If the column is not in a group AND there is
                    //no dependence with other colums use it in the gruop

                    if(ok){

                            idxGroups(group_index,gcount)=j;
                            gcount++;
                            colcount++;
                            col_check(j)=1;

                    }

        }

        //Save the size of this group

        sizeGroups(group_index)=gcount;

    }


    //Now just obtain the usable data

    Eigen::MatrixXd groups(group_index+1,maxSize);
    Eigen::VectorXd sizeG (group_index+1);

    groups=idxGroups.block(0,0,groups.rows(),groups.cols());

    sizeG=sizeGroups.head(group_index+1);


    //Return the information of the groups!

    problem.idxGroups=groups;
    problem.sizeGroups=sizeG;

    cout<<"-Number of index sets for sparse finite differences  ->"<<group_index+1<<endl;


}


//%%%%%%%%%%%%%%%%%%%% DERIVATIVE MATRIX FUNCTIONS %%%%%%%%%%%%%%%%%%%%

void detectDerivativeMatrixSparsity(void(*fun)(Eigen::VectorXd &x,Eigen::VectorXd &g,Problem &problem,Alg& algorithm),Eigen::VectorXd &x,Problem &problem,Alg& algorithm){

    int nDecVar=problem.nDecVar;
    int nCns=problem.nStates*(problem.nNodes)+problem.nEvents+(problem.nPath*problem.nNodes)+1;

    int nzcount_G=0; //Non-zero non-constant elements (counter)

    double s=1.0e6*sqrt(MC_EPSILON);


    double tol= 1.e-16*pow(MC_EPSILON,0.8)*MAX(1.0,x.norm());

    Eigen::VectorXd jacCol1(nCns);
    Eigen::VectorXd jacCol2(nCns);
    Eigen::VectorXd jacCol3(nCns);

    Eigen::VectorXd xlb(nDecVar);
    Eigen::VectorXd xub(nDecVar);
    Eigen::VectorXd xp(nDecVar);

    Eigen::VectorXd s_v(nDecVar);

    Eigen::MatrixXd Gidx(nCns*nDecVar,2);

    Eigen::MatrixXd jacobian(nCns,nDecVar);


    xlb=problem.xlb;
    xub=problem.xub;

    s_v.setOnes();

    s_v*=s;

    //Iterate along the columns of the jacobian matrix
    for(int j=0;j<nDecVar;j++){

        xp=x;

        utils::clip_vector_given_bounds(xp,xlb,xub);

        derivatives::numerical::jacobianColumn(fun,xp,problem,algorithm,j,jacCol1);

        //Perturbate the decision vector

        xp=x.array()+(0.1*x.array().abs())+s_v.array();

        utils::clip_vector_given_bounds(xp,xlb,xub);

        derivatives::numerical::jacobianColumn(fun,xp,problem,algorithm,j,jacCol2);

        //Perturbate the decision vector(Again...)

        xp=x.array()-(0.15*x.array().abs())-(1.1*s_v.array());

        utils::clip_vector_given_bounds(xp,xlb,xub);

        derivatives::numerical::jacobianColumn(fun,xp,problem,algorithm,j,jacCol3);

        //Now iterate alog the constraints (rows)

        for(int i=0;i<nCns;i++){

            if( (fabs(jacCol1(i)) +  fabs(jacCol2(i)) + fabs(jacCol3(i)) ) >= tol){

                    Gidx(nzcount_G,0)=i;//Row
                    Gidx(nzcount_G,1)=j;//Column

                    nzcount_G++;

                }
            }

        }

    //Now save the sparsity structure

    problem.nnzD=nzcount_G;

    //problem.nnzGD=nzcount_G;

    problem.nzGD.setZero(nzcount_G,2);

    problem.nzGD=Gidx.topRows(nzcount_G);


    cout<<endl;
    cout<<"---------------------------------------------------------"<<endl;
    cout<<"++++ Derivative Matrix Sparsity detected numerically ++++"<<endl;
    cout<<endl;
    cout<<"+Nonzero elements             ->"<<problem.nnzD<<endl;

}

void getIndexGroupsDerivativeMatrix(Problem &problem){

    /* This function uses the method of Curtis, Powell and Reid (1974) to find groups of variables
     * to evaluate efficiently the sparse Jacobian by perturbing simultaneously groups of variables.
     * Reference:
     * A. R. Curtis, M.J.D. Powell and J.K. Reid
     * "On the estimation of Sparse Jacobian Matrices"
     * J Inst Maths Applics (1974) 13, 117-119
     *
     */

    int nDecVar= problem.nDecVar;
    int nCns=   problem.nStates*(problem.nNodes)+problem.nEvents+(problem.nPath*problem.nNodes)+1;

    Eigen::MatrixXd nzG(problem.nnzD,2);

    nzG=problem.nzGD;

    //Define a matrix with 1 in the non-zero elements

    Eigen::SparseMatrix<double, Eigen::RowMajor> J(nCns,nDecVar);
    typedef Eigen::Triplet<double,int> T;
    vector<T> entries;

      for(int k=0; k<problem.nnzD;++k){
        entries.push_back( T(nzG(k,0), nzG(k,1), 1) );
      }

      J.setFromTriplets(entries.begin(), entries.end());

      Eigen::MatrixXd JDense(nCns,nDecVar);

      JDense=Eigen::MatrixXd(J);


      //********************************************
      std::ofstream file2("derivativeSparsity.txt");
        if (file2.is_open())
        {
          file2 << JDense << '\n';
        }
      //********************************************

      //cout<<endl;
      //cout<<"--------"<<endl;
      //cout<<JDense<<endl;

      Eigen::VectorXd C1(nCns);
      Eigen::VectorXd C2(nCns);

      Eigen::MatrixXd idxGroups(nDecVar,nDecVar);
      Eigen::VectorXd col_check(nDecVar);
      Eigen::VectorXd sizeGroups(nDecVar);


      col_check.setZero();
      idxGroups.setZero();
      sizeGroups.setZero();

        //Define the first group

      idxGroups(0,0)=0; //(Group,Index)=Column index

      col_check(0)=1; //The column 0 is already in a group...


      //Iterate along the columns of the jacobian matrix

      int gcount=1;
      int colcount=1;
      bool ok;
      double dotCols;


      for(int j=1;j<nDecVar;j++){

            ok=true;

            for(int l=0;l<gcount;l++){

                        if(j==idxGroups(0,l)){
                            ok=false;
                            break;
                        }

                         C1=JDense.col(idxGroups(0,l));
                         C2=JDense.col(j);

                         dotCols=C1.dot(C2);

                         if(dotCols>0.0){
                             ok=false;
                         }
            }

            if(ok){
                idxGroups(0,gcount)=j;
                gcount++;
                colcount++;
                col_check(j)=1;
            }
     }

      sizeGroups(0)=gcount;

      //Now lets obtain the other groups


    int group_index=0;

    int maxSize=0;

    while(colcount<nDecVar){

        group_index++;

        if(gcount>maxSize)
            maxSize=gcount;

        gcount=0;

        for(int j=1;j<nDecVar;j++){

                    ok=true;

                    if(col_check(j)==1){

                            ok=false;
                     }


                    //If the column is not in a group use it for a group

                    if(ok==true){

                          for(int l=0;l<gcount;l++){

                                    C1=JDense.col(idxGroups(group_index,l));
                                    C2=JDense.col(j);

                                    dotCols=C1.dot(C2);

                                     if(dotCols>0.0){

                                            ok=false;
                                      }

                            }
                     }

                    //If the column is not in a group AND there is
                    //no dependence with other colums use it in the gruop

                    if(ok){

                            idxGroups(group_index,gcount)=j;
                            gcount++;
                            colcount++;
                            col_check(j)=1;

                    }

        }

        //Save the size of this group

        sizeGroups(group_index)=gcount;

    }


    //Now just obtain the usable data

    Eigen::MatrixXd groups(group_index+1,maxSize);
    Eigen::VectorXd sizeG (group_index+1);

    groups=idxGroups.block(0,0,groups.rows(),groups.cols());

    sizeG=sizeGroups.head(group_index+1);


    //Return the information of the groups!

    problem.idxGroupsD=groups;
    problem.sizeGroupsD=sizeG;

    cout<<"-Number of index sets for sparse finite differences in the Derivative Matrix ->"<<group_index+1<<endl;


}

void getPerturbationMatrix(Problem &problem,Alg &algorithm){

    //This is a test function for the new version

    //Set a matrix with the perturbation vectors depending on the index groups of the jacobian

    double delj;
    int element;


    //Machine precision
    delj=sqrt(MC_EPSILON);

    //delj=1e-8;

    //Iterate along the groups

    int nGroups=problem.sizeGroupsD.size();

    Eigen::MatrixXd perturbationMatrix(problem.nDecVar,nGroups);
    Eigen::MatrixXd negativePerturbationMatrix(problem.nDecVar,nGroups);

    perturbationMatrix.setZero();
    negativePerturbationMatrix.setZero();

    for(int i=0;i<nGroups;i++){


        for(int j=0;j<problem.sizeGroupsD(i);j++){

            element=problem.idxGroupsD(i,j);

            perturbationMatrix(element,i)=delj;
            negativePerturbationMatrix(element,i)=2*delj;


        }

    }

    problem.sparsity.positivePerturbationMatrix=perturbationMatrix;
    problem.sparsity.negativePerturbationMatrix=negativePerturbationMatrix;



}

void generateSparsityTemplates(Problem &problem,Alg &algorithm){

    //This function generate the sparsity templates for the sparsity explotaiton in the jacobian of the constraints

    //By now this templates are just defined for the right hand sparsity explotaiton of the trapezoidal method (Version: 2.1)

    //c(x)=A x + B q(x)

    //d c(x)
    //----- = A + B D
    // dx

    //Variables

    int ny=problem.nStates;         //Number of states
    int nu=problem.nControls;       //Number of controls
    int M=problem.nNodes;           //Number of discrete points
    int nDecVar=problem.nDecVar;    //Number of decision variables
    int ns=M-1;                     //Number of segments
    int nPath=problem.nPath*M;      //Number of path constraints
    int nEvents=problem.nEvents;    //Number of event contraints
    int nDefects=problem.nDefectCns;//Number of defect constraints


    //-------------------------------Matrix A----------------------------

    Eigen::MatrixXd A(problem.nCns,nDecVar);
    A.setZero();

    Eigen::MatrixXd I_ns(ny,ny);
    I_ns.setIdentity();

    Eigen::MatrixXd I_u(ny,nu);
    I_u.setZero();

    Eigen::MatrixXd IInsu(ny,2*ny+2*nu);

    IInsu<<-I_ns,I_u,I_ns,I_u;

    for(int i=0;i<ns;i++){


        IInsu<<-I_ns,I_u,I_ns,I_u;
        A.block(i*ny,i*(ny+nu)+2,ny,2*ny+2*nu)=IInsu;

    }

    //-------------------------------Matrix B----------------------------

    Eigen::MatrixXd B(nDefects+(nPath)+nEvents+1,(ny*M)+(nPath)+nEvents+1);
    Eigen::MatrixXd II(ny,2*ny);
    Eigen::MatrixXd I(nEvents+nPath+1,nEvents+nPath+1);

    B.setZero();


    for(int i=0;i<ns;i++){

        II<<Eigen::MatrixXd::Identity(ny,ny)*problem.scale.defects.col(i).asDiagonal(),Eigen::MatrixXd::Identity(ny,ny)*problem.scale.defects.col(i).asDiagonal();
        B.block(i*ny,i*ny,ny,2*ny)=-0.5*II.array();

    }

    I<<Eigen::MatrixXd::Identity(nEvents+nPath+1,nEvents+nPath+1);
    B.block(nDefects,(ny*M),nEvents+nPath+1,nEvents+nPath+1)=I;



    //Return the sparsity structures as sparse matrices
    problem.sparsity.A=A.sparseView();
    problem.sparsity.B2=B.sparseView();

}

//%%%%%%%%%%%%%%%%%%% VERSION 3.0 SPARSITY FUNCTIONS %%%%%%%%%%%%%%%%%%%

void detectSparsity(Problem &problem){

    int nStates=problem.nStates;
    int nControls=problem.nControls;
    int nPath=problem.nPath;

    //Set offset of the constraints
    int eventOffset=problem.nNodes*problem.nStates;
    int pathOffset=(problem.nNodes*problem.nStates)+problem.nEvents;

    //Create the derivative matrix
    Eigen::MatrixXd derivativeMatrix(problem.derivativeMatrix.rows(),problem.derivativeMatrix.cols());
    derivativeMatrix.setZero();

    //Create the templates for the numerical derivatives of the user defined functions
    Eigen::MatrixXd dxTemplate(problem.nStates,problem.nStates+problem.nControls);
    Eigen::MatrixXd pathTemplate(problem.nPath,problem.nStates+problem.nControls);
    Eigen::MatrixXd eventTemplate_t0(problem.nEvents,problem.nStates+1);
    Eigen::MatrixXd eventTemplate_tF(problem.nEvents,problem.nStates+1);

    cout<<endl;
    cout<<"--------------> Detecting Derivative Matrix Sparsity <------------"<<endl;

    //Call the functions to obtain the sparsity template of each user defined function

    detectSparsity_dx(problem,dxTemplate);

    detectSparsity_path(problem,pathTemplate);

    detectSparsity_events(problem,eventTemplate_t0,eventTemplate_tF);

    for (int k=0; k<problem.nNodes;k++){

            derivativeMatrix.col(0).segment(k*nStates,nStates)=Eigen::VectorXd::Ones(nStates);      //Derivatives w.r.t. t0
            derivativeMatrix.col(1).segment(k*nStates,nStates)=Eigen::VectorXd::Ones(nStates);      //Derivatives w.r.t. tF

            derivativeMatrix.block(k*(nStates),2+(k*(nStates+nControls)),nStates,(nStates+nControls))=dxTemplate;

            derivativeMatrix.block(pathOffset+(k*(nPath)),2+(k*(nStates+nControls)),nPath,(nStates+nControls))=pathTemplate;

        }

    // SPARSITY OF THE EVENT CONSTRAINTS

       derivativeMatrix.block(eventOffset,0,problem.nEvents,1)=eventTemplate_t0.col(0);
       derivativeMatrix.block(eventOffset,1,problem.nEvents,1)=eventTemplate_tF.col(0);

       derivativeMatrix.block(eventOffset,2,problem.nEvents,problem.nStates)=eventTemplate_t0.rightCols(problem.nStates);
       derivativeMatrix.block(eventOffset,problem.nDecVar-(problem.nStates+problem.nControls),problem.nEvents,problem.nStates)=eventTemplate_tF.rightCols(problem.nStates);

       derivativeMatrix(derivativeMatrix.rows()-1,0)=1;
       derivativeMatrix(derivativeMatrix.rows()-1,1)=1;


   //Obtain the sparsity data of the derivative matrix

       Eigen::SparseMatrix<double> D(derivativeMatrix.rows(),derivativeMatrix.cols());
       int nnzD=0;
       int counter=0;

       D=derivativeMatrix.sparseView();
       nnzD=D.nonZeros();

       Eigen::MatrixXd nzElementsD(nnzD,2);
       nzElementsD.setZero();

       for (int k=0; k<D.outerSize(); ++k)
             for (Eigen::SparseMatrix<double>::InnerIterator it(D,k); it; ++it)
             {
               nzElementsD(counter,0)=it.row();   // row index
               nzElementsD(counter,1)=it.col();   // col index (here it is equal to k)
               counter++;

             }


       //Set the information of the sparsity

       problem.nnzD=nnzD;
       problem.nzGD.setZero(nnzD,2);
       problem.nzGD=nzElementsD;

       cout<<"+Non-zero elements in the Derivative Matrix->"<<problem.nnzD<<endl;

       cout<<endl;
       cout<<"------------------------------------------------------------------"<<endl;
       cout<<endl;
       cout<<"------------> Detecting Full Jacobian Sparsity <------------------"<<endl;
       cout<<endl;

       Eigen::SparseMatrix<double> J(problem.D.rows(),problem.D.cols());


       J=problem.sparsity.A+ (problem.sparsity.B2*D);

       int nnzJ=J.nonZeros();

       counter=0;

       Eigen::MatrixXd nzElementsJ(nnzJ,2);
       nzElementsJ.setZero();

       for (int k=0; k<J.outerSize(); ++k)
           for (Eigen::SparseMatrix<double>::InnerIterator it(J,k); it; ++it)
             {
               nzElementsJ(counter,0)=it.row();   // row index
               nzElementsJ(counter,1)=it.col();   // col index (here it is equal to k)
               counter++;
             }

       //Set the sparsity information

       problem.nnz=nnzJ;
       problem.nzG.setZero(nnzJ,2);
       problem.nzG=nzElementsJ;

       cout<<"+Nonzero elements in the Jacobian of the constrains-->"<<problem.nnz<<endl;


}

void detectSparsity_dx(Problem &problem,Eigen::MatrixXd &sparsityTemplate){


        int nDecVar=problem.nStates+problem.nControls;
        int nCns=problem.nStates;

        int nzcount_G=0; //Non-zero non-constant elements (counter)

        double s=1.0e6*sqrt(MC_EPSILON);

        Eigen::VectorXd x(nDecVar);
        Eigen::VectorXd xp(nDecVar);

        x.setRandom();

        double tol= 1.e-16*pow(MC_EPSILON,0.8)*MAX(1.0,x.norm());

        Eigen::VectorXd jacCol1(nCns);
        Eigen::VectorXd jacCol2(nCns);
        Eigen::VectorXd jacCol3(nCns);

        Eigen::VectorXd s_v(nDecVar);

        Eigen::MatrixXd Gidx(nCns*nDecVar,2);

        Eigen::MatrixXd jacobian(nCns,nDecVar);


        s_v.setOnes();

        s_v*=s;

        for(int j=0;j<nDecVar;j++){

            xp=x;

            derivatives::numerical::dxJacColumn(xp,j,jacCol1,problem);

            xp=x.array()+(0.1*x.array().abs())+s_v.array();

            derivatives::numerical::dxJacColumn(xp,j,jacCol2,problem);

            xp=x.array()-(0.15*x.array().abs())-(1.1*s_v.array());

            derivatives::numerical::dxJacColumn(xp,j,jacCol3,problem);

            for(int i=0;i<nCns;i++){

                if( (fabs(jacCol1(i)) +  fabs(jacCol2(i)) + fabs(jacCol3(i)) ) >= tol){


                        Gidx(nzcount_G,0)=i;//Row
                        Gidx(nzcount_G,1)=j;//Column

                        nzcount_G++;
                }

            }

        }


        Eigen::MatrixXd G(nzcount_G,2);

        G=Gidx.topRows(nzcount_G);

        Eigen::MatrixXd dxpattern(problem.nStates,problem.nStates+problem.nControls);
        dxpattern.setZero();

        for(int k=0; k<nzcount_G; k++){

            int row=G(k,0);
            int col=G(k,1);


            dxpattern(row,col)=1;


        }

        cout<<"+Non-zeros in the system dynamics Jacobian: "<<nzcount_G<<endl;


        //Return this!
        sparsityTemplate=dxpattern;


}

void detectSparsity_path(Problem &problem, Eigen::MatrixXd &sparsityTemplate_path){


       int nDecVar=problem.nStates+problem.nControls;
       int nCns=problem.nPath;

       int nzcount_G=0; //Non-zero non-constant elements (counter)

       double s=1.0e6*sqrt(MC_EPSILON);

       Eigen::VectorXd x(nDecVar);
       Eigen::VectorXd xp(nDecVar);

       x.setRandom();

       double tol= 1.e-16*pow(MC_EPSILON,0.8)*MAX(1.0,x.norm());

       Eigen::VectorXd jacCol1(nCns);
       Eigen::VectorXd jacCol2(nCns);
       Eigen::VectorXd jacCol3(nCns);

       Eigen::VectorXd s_v(nDecVar);

       Eigen::MatrixXd Gidx(nCns*nDecVar,2);

       Eigen::MatrixXd jacobian(nCns,nDecVar);


       s_v.setOnes();

       s_v*=s;

       for(int j=0;j<nDecVar;j++){

           xp=x;

           derivatives::numerical::pathJacColumn(xp,j,jacCol1,problem);

           xp=x.array()+(0.1*x.array().abs())+s_v.array();

           derivatives::numerical::pathJacColumn(xp,j,jacCol2,problem);

           xp=x.array()-(0.15*x.array().abs())-(1.1*s_v.array());

           derivatives::numerical::pathJacColumn(xp,j,jacCol3,problem);

           for(int i=0;i<nCns;i++){

               if( (fabs(jacCol1(i)) +  fabs(jacCol2(i)) + fabs(jacCol3(i)) ) >= tol){


                       Gidx(nzcount_G,0)=i;//Row
                       Gidx(nzcount_G,1)=j;//Column

                       nzcount_G++;
               }

           }

       }


       Eigen::MatrixXd G(nzcount_G,2);

       G=Gidx.topRows(nzcount_G);

       Eigen::MatrixXd pathPattern(problem.nPath,problem.nStates+problem.nControls);
       pathPattern.setZero();

       for(int k=0; k<nzcount_G; k++){

           int row=G(k,0);
           int col=G(k,1);


           pathPattern(row,col)=1;


       }

       cout<<"+Non-zeros in the path constraints Jacobian: "<<nzcount_G<<endl;

       //Return this!!
       sparsityTemplate_path=pathPattern;

}

void detectSparsity_events(Problem &problem,Eigen::MatrixXd &sparsityTemplate_t0,Eigen::MatrixXd &sparsityTemplate_tf){


    int nDecVar=problem.nStates+1;
    int nCns=problem.nEvents;

    int nzcount_G=0; //Non-zero non-constant elements (counter)

    double s=1.0e6*sqrt(MC_EPSILON);

    Eigen::VectorXd x(nDecVar);
    Eigen::VectorXd xp(nDecVar);

    x.setRandom();

    double tol= 1.e-16*pow(MC_EPSILON,0.8)*MAX(1.0,x.norm());

    Eigen::VectorXd jacCol1(nCns);
    Eigen::VectorXd jacCol2(nCns);
    Eigen::VectorXd jacCol3(nCns);

    Eigen::VectorXd s_v(nDecVar);

    Eigen::MatrixXd Gidx(nCns*nDecVar,2);

    Eigen::MatrixXd jacobian(nCns,nDecVar);


    s_v.setOnes();

    s_v*=s;

    for(int j=0;j<nDecVar;j++){

        xp=x;

        derivatives::numerical::eventJacColumn_t0(xp,j,jacCol1,problem);

        xp=x.array()+(0.1*x.array().abs())+s_v.array();

        derivatives::numerical::eventJacColumn_t0(xp,j,jacCol2,problem);

        xp=x.array()-(0.15*x.array().abs())-(1.1*s_v.array());

        derivatives::numerical::eventJacColumn_t0(xp,j,jacCol3,problem);

        for(int i=0;i<nCns;i++){

            if( (fabs(jacCol1(i)) +  fabs(jacCol2(i)) + fabs(jacCol3(i)) ) >= tol){

                    Gidx(nzcount_G,0)=i;//Row
                    Gidx(nzcount_G,1)=j;//Column

                    nzcount_G++;
            }

        }

    }


    Eigen::MatrixXd G(nzcount_G,2);

    G=Gidx.topRows(nzcount_G);

    Eigen::MatrixXd eventPattern(problem.nEvents,problem.nStates+1);
    eventPattern.setZero();

    for(int k=0; k<nzcount_G; k++){

        int row=G(k,0);
        int col=G(k,1);

        eventPattern(row,col)=1;

    }

    cout<<"+Non-zero elements in the initial event constraints: "<<nzcount_G<<endl;


    //%%%%%%%%%%%%%%%%%%%%%%%%%% EVENTS w.r.t. tF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nzcount_G=0;

    s_v.setOnes();

    s_v*=s;

    for(int j=0;j<nDecVar;j++){

        xp=x;

        derivatives::numerical::eventJacColumn_tf(xp,j,jacCol1,problem);

        xp=x.array()+(0.1*x.array().abs())+s_v.array();

        derivatives::numerical::eventJacColumn_tf(xp,j,jacCol2,problem);

        xp=x.array()-(0.15*x.array().abs())-(1.1*s_v.array());

        derivatives::numerical::eventJacColumn_tf(xp,j,jacCol3,problem);

        for(int i=0;i<nCns;i++){

            if( (fabs(jacCol1(i)) +  fabs(jacCol2(i)) + fabs(jacCol3(i)) ) >= tol){

                    Gidx(nzcount_G,0)=i;//Row
                    Gidx(nzcount_G,1)=j;//Column

                    nzcount_G++;
            }

        }

    }


    G=Gidx.topRows(nzcount_G);

    Eigen::MatrixXd eventPattern_tf(problem.nEvents,problem.nStates+1);
    eventPattern_tf.setZero();

    for(int k=0; k<nzcount_G; k++){

        int row=G(k,0);
        int col=G(k,1);


        eventPattern_tf(row,col)=1;


    }

    //Return this!

    cout<<"+Non-zero elements in the final event constraints: "<<nzcount_G<<endl;

    sparsityTemplate_t0=eventPattern;
    sparsityTemplate_tf=eventPattern_tf;


}

}

} //END NLP namespace

} //END OCSolver namespace


