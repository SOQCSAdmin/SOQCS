//======================================================================================================
// File qocircuit.cpp
//
// OPTICAL CIRCUIT LIBRARY.
//
// Copyright © 2023 National University of Ireland Maynooth, Maynooth University. All rights reserved.
// The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the
// root directory of this source tree. Use of the source code in this file is only permitted under the
// terms of the licence in LICENCE.TXT.
//======================================================================================================
#include "qocircuit.h"


//----------------------------------------
//
//  Create circuit
//
//----------------------------------------
qocircuit::qocircuit(int i_nch){
//  int i_nch;  // Number of channels


    create_circuit(i_nch, 1, 1, 1, -1.0 ,0, 0, false,'G');
}


//----------------------------------------
//
//  Create circuit
//
//----------------------------------------
qocircuit::qocircuit(int i_nch, int i_nm, int i_ns){
//  int i_nch;  // Number of channels
//  int i_nm;   // Number of modes
//  int i_ns;   // Number of wavepackets


    create_circuit(i_nch, i_nm, i_ns, 1, -1.0, 0 ,0, false, 'G');
}


//----------------------------------------
//
//  Create circuit
//
//----------------------------------------
qocircuit::qocircuit(int i_nch, int i_nm, int i_ns, int clock, char i_ckind){
//  int  i_nch;   // Number of channels
//  int  i_nm;    // Number of modes
//  int  i_ns;    // Number of wavepackets
//  int  clock;   // Kind of detectors
//  char i_ckind; // Kind of packets 'G'=Gaussian/'E'=Exponential


    create_circuit(i_nch, i_nm, i_ns, 1,-1.0, clock ,0, false, i_ckind);
}

//----------------------------------------
//
//  Create circuit
//
//----------------------------------------
qocircuit::qocircuit(int i_nch, int i_nm, int i_ns, int i_np, double i_dtp, int clock, int i_R, bool loss,char i_ckind){
//  int    i_nch;    // Number of channels
//  int    i_nm;     // Number of modes
//  int    i_ns;     // Number of wavepackets
//  int    i_np;     // Number of periods.
//  double dtp;      // Period length
//  int    clock;    // Kind of detectors
//  int    i_R;      // Number of iterations to calculate blinking and dark counts
//  bool   loss;     // Are losses going to be calculated explicitly? True=Yes/False=No
//  char   i_ckind;  // Kind of packets 'G'=Gaussian/'E'=Exponential


    create_circuit(i_nch, i_nm, i_ns, i_np, i_dtp, clock, i_R, loss, i_ckind);
}


//----------------------------------------
//
//  Auxiliary private method to create a circuit
//
//----------------------------------------
int qocircuit::create_circuit(int i_nch, int i_nm, int i_ns, int i_np, double i_dtp, int clock, int i_R, bool loss, char i_ckind){
//  int    i_nch;    // Number of channels
//  int    i_nm;     // Number of modes
//  int    i_ns;     // Number of wavepackets
//  int    i_np;     // Number of periods.
//  double dtp;      // Period length
//  int    clock;    // Kind of detectors
//  int    i_R;      // Number of iterations to calculate blinking and dark counts
//  bool   loss;     // Are losses going to be calculated explicitly? True=Yes/False=No
//  char   i_ckind;  // Kind of packets 'G'=Gaussian/'E'=Exponential
//  Auxiliary index
    int    i;        // Aux index
    int    j;        // Aux index
    int    k;        // Aux index
    int    l;        // Aux index


    if((i_nch<=0)||(i_nm<=0)||(i_ns<=0)||(i_np<=0)||(i_R<0)){
        cout << "Create circuit error #1: Number of photons, modes, times and periods have to be greater than zero." << endl;
        cout << "R has to be positive" << endl;
        return -1;
    }

    if((i_np>1)&&(i_dtp<=0.0)){
        cout << "Create circuit error #2: It is not possible to define more than one period with negative or zero period lengths." << endl;
        return -1;
    }

    // Initialize degrees of freedom
    nm=i_nm;
    ns=i_ns*i_np;
    np=i_np;
    nsp=i_ns;
    dtp=i_dtp;

    // If we are going to compute losses.
    // We need the double of channels.
    if(loss){
        losses=1;
        nch=2*i_nch;
    }else{
        losses=0;
        nch=i_nch;
    }

    // Initialize detector conditions
    ndetc=0;
    nignored=0;
    ncond=0;
    timed=clock;
    det_def.resize(3,i_nch);
    det_win.resize(2,i_nch);
    det_par.resize(2,i_nch);
    ch_ignored.resize(i_nch);


    // Initialize emitter/Packet definitions
    ckind=i_ckind;
    npack=0;
    pack_list.resize(5,nsp);
    emitted=new photon_mdl();
    emiss=0;
    confidence=1.0;

    // Initialize blink, dark counts and noise parameters
    R=i_R;
    dev=0.0;

    // Create dictionaries
    nlevel=nch*nm*ns;
    idx = new level[nlevel];
    i_idx= new int**[nch];
    for(i=0;i<nch;i++){
        i_idx[i]=new int*[nm];
        for(j=0;j<nm;j++){
            i_idx[i][j]=new int[ns];
        }
    }

    // Initialize dictionaries
    l=0;
    for(i=0;i<nch;i++){
        for(j=0;j<nm;j++){
            for(k=0;k<ns;k++){
                idx[l].ch=i;
                idx[l].m=j;
                idx[l].s=k;
                i_idx[i][j][k]=l;
                l++;
            }
        }
    }


    // Create and initialize circuit.
    // No elements in circuit is the identity matrix
    circmtx=matc::Identity(nlevel,nlevel);

    // Return success
    return 0;
}


//----------------------------------------
//
//  Destroy circuit
//
//----------------------------------------
qocircuit::~qocircuit(){
//  Auxiliary index
    int i; // Aux index
    int j; // Aux index


    // Destroy dictionaries/ liberate memory.
    for(i=0;i<nch;i++){
        for(j=0;j<nm;j++){
            delete[] i_idx[i][j];
        }
        delete[] i_idx[i];
    }

    // Free memory
    delete[] idx;
    delete[] i_idx;
    delete emitted;
}


//----------------------------------------------------------
//
//  Reset circuit matrix (but maintains level definitions)
//
//----------------------------------------------------------
void qocircuit::reset(){


    // Resent circuit matrix
    circmtx=matc::Identity(nlevel,nlevel);
    // Reset detectors
    ndetc=0;
    nignored=0;
    ncond=0;
    // Reset emitter/packets
    npack=0;
    emiss=0;
    confidence=1.0;
    delete emitted;
    emitted=new photon_mdl();
    // Reset noise
    dev=0.0;

}


//----------------------------------------
//
//  Copy a circuit
//
//----------------------------------------
qocircuit *qocircuit::clone(){
//  Variables
    qocircuit* newcircuit;      // New circuit recipient of the copy
//  Auxiliary index
    int i;                      // Aux index
    int j;                      // Aux index
    int k;                      // Aux index
    int l;                      // Aux index


    //Reserve memory
    newcircuit= new qocircuit(nch,nm,ns,1,dtp, timed, 0, false, ckind);
    newcircuit->np=np;
    newcircuit->nsp=nsp;
//    newcircuit->mpi=mpi;
//    newcircuit->mpf=mpf;

    newcircuit->losses=losses; // It has to be cloned that way. Otherwise the full constructor will double again the number of channels.

    // Copy dictionaries
    l=0;
    for(i=0;i<nch;i++){
        for(j=0;j<nm;j++){
            for(k=0;k<ns;k++){
                newcircuit->idx[l].ch=idx[l].ch;
                newcircuit->idx[l].m=idx[l].m;
                newcircuit->idx[l].s=idx[l].s;
                newcircuit->i_idx[i][j][k]=i_idx[i][j][k];
                l++;
            }
        }
    }

    // Copy detectors
    newcircuit->ndetc=ndetc;
    newcircuit->nignored=nignored;
    newcircuit->ncond=ncond;
    newcircuit->det_def=det_def;
    newcircuit->det_win=det_win;
    newcircuit->det_par=det_par;
    newcircuit->ch_ignored=ch_ignored;

    // Copy emitters
    newcircuit->emitted=emitted->clone();
    newcircuit->npack=npack;
    newcircuit->pack_list=pack_list;
    newcircuit->init_dmat=init_dmat;
    newcircuit->emiss=emiss;
    newcircuit->confidence=confidence;

    // Copy noise and detection effects
    newcircuit->R=R;
    newcircuit->dev=dev;

    // Copy the circuit configuration
    newcircuit->circmtx=circmtx;

    //Return a pointer to the new circuit object.
    return newcircuit;
}


//----------------------------------------
//
//  Concatenates two circuits
//
//----------------------------------------
int qocircuit::concatenate(qocircuit *qoc){
//  qocircuit *qoc;      // Circuit to be appended
//  Variables
    int nclose;          // Total number of channels


    // Check limitations
    if(nch!=qoc->nch){
        cout << "Concatenate error  #1: Circuits not compatible. Number of channels is different" << endl;
        return -1;
    }
    if(nm!=qoc->nm){
        cout << "Concatenate error  #2: Circuits not compatible. Number of polarization is different" << endl;
        return -1;
    }
    if(ns!=qoc->ns){
        cout << "Concatenate error  #3: Circuits not compatible. Number of packets is different" << endl;
        return -1;
    }
    if(qoc->npack>0){
        cout << "Concatenate error  #4: Input must be defined entirely on the first circuit! " << endl;
        return -1;
    }
    if(ncond>0){
        cout << "Concatenate error  #5: Detectors must be defined entirely on the last circuit!. " << endl;
        return -1;
    }

    if (((losses==0)&&(qoc->losses>0))||((losses>0)&&(qoc->losses==0))){
        cout << "Concatenate error  #6: Both circuits must have the same loss configuration!. " << endl;
        return -1;
    }

    // Merges circuits (in a limited way)
    circmtx=qoc->circmtx*circmtx;
    ncond=qoc->ncond;
    ndetc=qoc->ndetc;
    nignored=qoc->nignored;
    det_def=qoc->det_def;
    det_win=qoc->det_win;
    det_par=qoc->det_par;
    ch_ignored=qoc->ch_ignored;

    // Last detector operations
    // The total number of channels to put a detector
    // is different depending if we are computing losses
    if(losses==0) nclose=nch;
    else nclose=nch/2;

    // If this is the last detector compute losses if needed.
    if((losses>0)&&(ndetc==nclose)) compute_losses();

    // If this is the last detector add the emitter matrix
    // at the beginning. (Otherwise computing the losses breaks the calculation)
    if((emiss==1)&&(ndetc==nclose)) circmtx=circmtx*init_dmat;

    return 0;
}


//---------------------------------------------------------
//
//  Returns the number of levels of this circuit
//
//---------------------------------------------------------
int qocircuit::num_levels(){
    return nlevel;
}


//----------------------------------------
//
//  Adds a beamsplitter to the circuit
//
//----------------------------------------
int qocircuit:: beamsplitter(int i_ch1, int i_ch2, double d_theta, double d_phi){
//  int    i_ch1;             // Beamsplitter input channel 1
//  int    i_ch2;             // Beamsplitter input channel 2
//  double d_theta;           // Angle theta in degrees
//  double d_phi;             // Angle phi in degrees
//  Constants
    const int nbmch=2;        // U square matrix row/column dimension
//  Variables
    double theta;             // Angle theta in radians
    double phi;               // Angle phi in radians
    matc   U(nbmch,nbmch);    // Beamsplitter 2x2 matrix definitions
    mati   V(1,nbmch);        // Channels to which U refers


    // Conversion to radians
    theta=d_theta*pi/180.0;
    phi  =d_phi*pi  /180.0;

    // 2x2 Matrix initialization
    U(0,0)= cos(theta);
    U(0,1)=-exp( jm*phi)*sin(theta);
    U(1,0)= exp(-jm*phi)*sin(theta);
    U(1,1)= cos(theta);

    // Channels
    V(0,0)  = i_ch1;
    V(0,1)  = i_ch2;

    // Compute gate
    return custom_gate(V,U);
}


//----------------------------------------
//
//  Adds a dielectric to the circuit
//
//----------------------------------------
int qocircuit:: dielectric(int i_ch1, int i_ch2, cmplx t, cmplx r){
//  int   i_ch1;              // Dielectric input channel 1
//  int   i_ch2;              // Dielectric input channel 2
//  cmplx t;                  // Transmission amplitude
//  cmplx r;                  // Reflection amplitude
//  Constants
    const int nbmch=2;        // U square matrix row/column dimension
//  Variables
    cmplx rmt;                // r-t
    cmplx rpt;                // r+t
    cmplx rmt2;               // |r-t|^2
    cmplx rpt2;               // |r+t|^2
    matc  U(nbmch,nbmch);     // Dielectric 2x2 matrix definitions
    mati  V(1,nbmch);         // Channels to which U refers


    // Check physicality of the input
    rmt=r-t;
    rpt=r+t;
    rmt2=conj(rmt)*rmt;
    rpt2=conj(rpt)*rpt;
    if(real(abs(rmt2))>1.0) cout << "Dielectric warning: t-r condition broken. The result aren't physical." << endl;
    if(real(abs(rpt2))>1.0) cout << "Dielectric warning: t+r condition broken. The result aren't physical." << endl;

    // Dielectric matrix
    U(0,0)= t;
    U(0,1)= r;
    U(1,0)= r;
    U(1,1)= t;

    // Channels
    V(0,0)  = i_ch1;
    V(0,1)  = i_ch2;

    // Compute gate
    return custom_gate(V,U);

}


//-------------------------------------------------------
//
//  Adds a virtual circuit element to compute the states
//  with less photons at the input than the output.
//
//-------------------------------------------------------
void qocircuit::compute_losses(){
//  Variables
    matc M;           // Circuit non unitary matrix
    matc off;         // Off diagonal matrix Rsqrt(1-D^2)V
    matc offd;        // Off diagonal matrix sqrt(1-D^2)
    matc D;           // Diagonal matrix M=RDV
    matc R;           // Unitary matrix R
    matc V;           // Unitary matrix V
    BDCSVD<matc> svd; // Single value decomposition
//  Auxiliary index
    int  i;           // Aux index

    //Single value decomposition
    M=circmtx.block(0,0,nlevel/2,nlevel/2);
    svd.compute(M,ComputeThinU | ComputeThinV);
    D=svd.singularValues().asDiagonal();
    R=svd.matrixU();
    V=svd.matrixV().adjoint();

    // Compute off diagonal couplings
    offd.setZero(nlevel/2,nlevel/2);
    for(i=0;i<nlevel/2;i++) offd(i,i)=sqrt(abs(1.0-abs(conj(D(i,i))*D(i,i))));
    off=R*offd*V;

    // Build the extended matrix.
    // Warning! Note that is assumed that the extended levels
    // are ordered in the same way than the physical ones.
    // Therefore it is not needed to check indexes.
    circmtx.block(0,0,nlevel/2,nlevel/2)=M;
    circmtx.block(0,nlevel/2,nlevel/2,nlevel/2)=off;
    circmtx.block(nlevel/2,0,nlevel/2,nlevel/2)=off;
    circmtx.block(nlevel/2,nlevel/2,nlevel/2,nlevel/2)=-M;
}


//----------------------------------------
//
//  Creates a random circuit
//
//----------------------------------------
void qocircuit:: random_circuit(){
//  Variables
    cmplx trow;               // Row normalization constant
    HouseholderQR<matc> qr;   // Householder QR matrix;
//  Auxiliary index.
    int i;                    // Channel 1 matrix index
    int j;                    // Channel 2 matrix index


    // Create random non-unitary matrix.
    for(i=0;i<nlevel;i++){
        for(j=0;j<nlevel;j++){
            circmtx(i,j)=urand()*exp(jm*2.0*pi*urand());

        }
    }

    // We ask to the final matrix to be unitary.
    // Perform a QR decomposition and return
    // the unitary matrix Q.
    qr.compute(circmtx);
    circmtx = qr.householderQ();

}


//----------------------------------------
//
//  Adds a NSX circuit element.
//
//----------------------------------------
int qocircuit:: NSX(int i_ch1, int i_ch2, int i_ch3){
//  int   i_ch1              // NSX input channel 1
//  int   i_ch2              // NSX input channel 2
//  int   i_ch3              // NSX input channel 3
//  Constants
    const int nbmch=3;       // U square matrix row/column dimension
//  Variables
    matc  U(nbmch,nbmch);    // NSX 2x2 matrix definitions
    mati  V(1,nbmch);        // Channels to which U refers


    // 3x3 Matrix initialization
    U(0,0)= 1.0-sqrt(2.0);
    U(0,1)= 1/sqrt(sqrt(2.0));
    U(0,2)= sqrt(3.0/sqrt(2)-2.0);
    U(1,0)= 1/sqrt(sqrt(2.0));
    U(1,1)= 1.0/2.0;
    U(1,2)= 1.0/2.0-1.0/sqrt(2.0);
    U(2,0)= sqrt(3.0/sqrt(2)-2.0);
    U(2,1)= 1.0/2.0-1.0/sqrt(2.0);;
    U(2,2)= sqrt(2.0)-1.0/2.0;

    // Channels
    V(0,0)  = i_ch1;
    V(0,1)  = i_ch2;
    V(0,2)  = i_ch3;

    // Compute gate
    return custom_gate(V,U);
}


//----------------------------------------
//
//  Adds a polarization rotator to the circuit
//
//----------------------------------------
int qocircuit:: rotator(int i_ch, double d_theta, double d_phi){
//  int    i_ch               // Rotator input channel 1
//  double d_theta;           // Angle theta in degrees
//  double d_phi;             // Angle phi in degrees
//  Constants
    const  int nbmch=2;       // U square matrix row/column dimension
//  Variables
    double theta;             // Angle theta in radians
    double phi;               // Angle phi in radians
    matc   U(nbmch,nbmch);    // Rotator 2x2 matrix definitions
    mati   W(2,nbmch);        // Channels to which U refers


    // Conversion to radians
    theta=d_theta*pi/180.0;
    phi  =d_phi*pi  /180.0;

    // 2x2 Matrix initialization
    U(0,0) =  cos(theta);
    U(0,1) = -exp( jm*phi)*sin(theta);
    U(1,0) =  exp(-jm*phi)*sin(theta);
    U(1,1) =  cos(theta);

    // Channels
    W(0,0) = i_ch;
    W(0,1) = i_ch;
    W(1,0) = H;
    W(1,1) = V;

    // Compute gate
    return custom_gate(W,U);
}


//----------------------------------------
//
//  Adds a polarized beamsplitter to the circuit
//
//----------------------------------------
int qocircuit:: pol_beamsplitter(int i_ch1, int i_ch2, int pol){
//  int   i_ch1               // Polarized beamsplitter input channel 1
//  int   i_ch2               // Polarized beamsplitter input channel 2
//  int   pol;                // Polarization to be switched
//  Constants
    const int nbmch=4;        // U square matrix row/column dimension
//  Variables
    matc  U(nbmch,nbmch);     // Polarized beamsplitter 4x4 matrix definitions
    mati  W(2,nbmch);         // Channels to which U refers


    // 4x4 Matrix initialization
    U(0,0) =  0;
    U(0,1) =  0;
    U(0,2) =  1;
    U(0,3) =  0;
    U(1,0) =  0;
    U(1,1) =  1;
    U(1,2) =  0;
    U(1,3) =  0;
    U(2,0) =  1;
    U(2,1) =  0;
    U(2,2) =  0;
    U(2,3) =  0;
    U(3,0) =  0;
    U(3,1) =  0;
    U(3,2) =  0;
    U(3,3) =  1;

    // Channels
    W(0,0) =  i_ch1;
    W(0,1) =  i_ch1;
    W(0,2) =  i_ch2;
    W(0,3) =  i_ch2;
    W(1,0) =  pol;
    W(1,1) = (pol+1)%2;
    W(1,2) =  pol;
    W(1,3) = (pol+1)%2;

    // Compute gate
    return custom_gate(W,U);
}


//----------------------------------------
//
//  Adds a general waveplate to the circuit
//
//----------------------------------------
int qocircuit:: waveplate(int i_ch, double d_alpha, double d_gamma){
//  int    i_ch               // Waveplate input channel
//  double d_alpha;           // Angle of the polarization rotation
//  double d_gamma;           // Wavenumber delay of one component with respect the other
//  Constants
    const  int nbmch=2;       // U square matrix row/column dimension
//  Variables
    double alpha;             // Angle alpha in radians
    double gamma;             // Angle gamma in radians
    matc   U(nbmch,nbmch);    // Waveplate 2x2 matrix definitions
    mati   W(2,nbmch);        // Channels to which U refers


    // Conversion to radians
    alpha=d_alpha*pi/180.0;
    gamma=d_gamma*pi/180.0;

    // 2x2 Matrix initialization
    U(0,0)= jm*sin(gamma)*cos(2.0*alpha)+cos(gamma);
    U(0,1)= jm*sin(gamma)*sin(2.0*alpha);
    U(1,0)= jm*sin(gamma)*sin(2.0*alpha);
    U(1,1)=-jm*sin(gamma)*cos(2.0*alpha)+cos(gamma);

    // Channels
    W(0,0)  = i_ch;
    W(0,1)  = i_ch;
    W(1,0)  = H;
    W(1,1)  = V;

    // Compute gate
    return custom_gate(W,U);
}


//----------------------------------------
//
//  Adds a half waveplate to the circuit
//
//----------------------------------------
int qocircuit:: half(int i_ch, double alpha){
//  int    i_ch               // Waveplate input channel
//  double d_alpha;           // Angle of the polarization rotation


    return waveplate(i_ch,alpha,90.0);
}


//----------------------------------------
//
//  Adds a quarter waveplate to the circuit
//
//----------------------------------------
int qocircuit:: quarter(int i_ch, double alpha){
//  int    i_ch               // Waveplate input channel
//  double d_alpha;           // Angle of the polarization rotation

    return waveplate(i_ch,alpha,45.0);
}


//----------------------------------------
//
//  Adds a 2x2 MMI to the circuit
//
//----------------------------------------
int qocircuit:: MMI2(int i_ch1, int i_ch2){
//  int   i_ch1               // MMI2 input channel 1
//  int   i_ch2               // MMI2 input channel 2


    return dielectric(i_ch1,i_ch2,1/sqrt(2.0),jm/sqrt(2.0));
}


//------------------------------------------------
//
// Rewires / swaps two channels (of the same mode)
//
//------------------------------------------------
int qocircuit:: rewire(int i_ch1,int i_ch2){
//  int    i_ch1;             // Channel 1
//  int    i_ch2;             // Channel 2
//  Constants
    const int nbmch=2;        // U square matrix row/column dimension
//  Variables
    matc   U(nbmch,nbmch);    // Beamsplitter 2x2 matrix definitions
    mati   V(1,nbmch);        // Channels to which U refers


    // 2x2 Matrix initialization
    U(0,0)= 0.0;
    U(0,1)= 1.0;
    U(1,0)= 1.0;
    U(1,1)= 0.0;

    // Channels
    V(0,0)  = i_ch1;
    V(0,1)  = i_ch2;

    // Compute gate
    return custom_gate(V,U);
}


//------------------------------------------------
//
// Adds a gate using other circuit as the gate definition
//
//------------------------------------------------
int qocircuit:: add_gate(veci chlist, qocircuit *qoc){
//  veci chlist;     // List of channels where the gate is added
//  qocircuit *qoc;  // Quantum optical circuit
//  Variables
    int    ch;       // Channel
    int    m;        // Polarization
    int    ch1;      // Channel 1
    int    ch2;      // Channel 2
    int    m1;       // Polarization channel 1
    int    m2;       // Polarization channel 2
    int    nclose;   // number of detectors needed to close the circuit.
    matc   U;        // Gate matrix definition
    mati   W;        // Channels to which U refers
//  Auxiliary index
    int    i;        // Aux index
    int    j;        // Aux index
    int    k;        // Aux index
    int    l;        // Aux index

    // Check limitations
    if (chlist.size()!=qoc->nch){
        cout << "add_gate error (qocircuit) #1: The number of channels in the list has to be the same than in the gate circuit." << endl;
        return -1;
    }

    if (qoc->losses>0){
        cout << "add_gate error (qocircuit) #2: Losses must be configured to False in the gate circuit." << endl;
        return -1;
    }

    if ((qoc->emiss==0)&&(remdec()==qoc->ndetc)){
        cout << "add_gate error (qocircuit) #3: Photons should be emitted before adding all detectors." << endl;
        return -1;
    }

    if ((qoc->emiss==1)&&(remdec()!=qoc->ndetc)){
        cout << "add_gate error (qocircuit) #4: Photons have been already emitted." << endl;
        cout << "Note that this may happen because the definition of a delay inside the gate. Delays are forbidden in gates." << endl;
        return -1;
    }


    if (qoc->np>1){
        cout << "add_gate error (qocircuit) #5: The number of periods must be one. We can not define periods in gates." << endl;
        return -1;
    }


    // Copy gate definition
    if(qoc->ns>1){
        // From a non-ideal circuit with packets
        U.resize(qoc->nch*qoc->nm,qoc->nch*qoc->nm);
        W.resize(2,qoc->nm*qoc->nch);
        i=0;
        for(ch1=0; ch1<qoc->nch; ch1++){
        for(m1=0; m1<qoc->nm; m1++){
            j=0;
            W(0,i)=chlist(ch1);
            W(1,i)=m1;
            for(ch2=0; ch2<qoc->nch; ch2++){
            for(m2=0; m2<qoc->nm; m2++){
                k=qoc->i_idx[ch1][m1][0];
                l=qoc->i_idx[ch2][m2][0];

                U(i,j)=qoc->circmtx(k,l);
                j=j+1;
            }}
            i=i+1;
        }}
    }else{
        // From an ideal circuit
        U=qoc->circmtx;
        W.resize(2,qoc->nm*qoc->nch);
        for(i=0; i<qoc->nlevel; i++){
            ch=qoc->idx[i].ch;
            m=qoc->idx[i].m;
            W(0,i)=chlist(ch);
            W(1,i)=m;

        }
    }

    // Create gate
    custom_gate(W,U);


    // Update detectors if post-selection is defined
    ndetc=ndetc+qoc->ndetc;
    if(ndetc>nch){
        cout << "add_gate error (qocircuit) #6: More detectors than channels are being declared." << endl;
        return -1;
    }

    for(i=0;i<qoc->ncond;i++){
        det_def(0,ncond)=chlist(qoc->det_def(0,i));
        det_def(1,ncond)=qoc->det_def(1,i);
        det_def(2,ncond)=qoc->det_def(2,i);
        ncond=ncond+1;
    }

    for(i=0;i<qoc->nignored;i++){
        ch_ignored(nignored)=chlist(qoc->ch_ignored(i));
        nignored=nignored+1;
    }

    for(i=0;i<qoc->nch;i++){
        det_win(0,chlist(i))=qoc->det_win(0,i);
        det_win(1,chlist(i))=qoc->det_win(1,i);
        det_par(0,chlist(i))=qoc->det_par(0,i);
        det_par(1,chlist(i))=qoc->det_par(1,i);
    }

    // Evaluate if the circuit is closed.
    // If the circuit is closed take the same action than in detector.

    // The total number of channels to put a detector is different depending if we are computing losses
    if(losses==0) nclose=nch;
    else nclose=nch/2;

    // If this is the last detector compute losses if needed.
    if((losses==1)&&(ndetc==nclose)) compute_losses();

    // If this is the last detector add the emitter matrix
    // at the beginning. (Otherwise computing the losses breaks the calculation)
    if((emiss==1)&&(ndetc==nclose)) circmtx=circmtx*init_dmat;

    // Return success
    return 0;
}
//----------------------------------------
//
//  Adds a custom gate to the circuit
//
//----------------------------------------
int qocircuit:: custom_gate(mati iodef, matc U){
//  mati iodef;       // List of channels and polarization the define the input to the gate
//  matc U;           // Custom gate nxn matrix definitions
//  Variables
    int  nbmch;       // U square matrix row/column dimension
    int  ch1;         // Channel evaluated 1
    int  ch2;         // Channel evaluated 2
    int  km;          // Mode of the channel i_ch
    int  km1;         // Mode of the channel i_ch
    int  km2;         // Mode of the channel i_ch
    int  ks;          // Wavepacket index
    veci CH;          // Channels to which U refers
    veci P;           // Channels to which U refers
    matc oelement;    // Extended nlevelxnlevel custom gate definition
//  Auxiliary index.
    int  i;           // Aux index
    int  j;           // Aux index
    int  k;           // Aux index
    int  l;           // Aux index.


    // Definitions format adaptation
    nbmch=U.cols();
    CH.resize(nbmch);
    P.resize(nbmch);
    for(i=0;i<nbmch;i++){
        CH(i)  = iodef(0,i);
        if (CH(i)>=nch){
            cout << "Gate error: Gate declared in an undefined channel." << endl;
            return -1;
        }
        if(iodef.rows()>1){
            P(i)   = iodef(1,i);
            if (P(i)>=nm){
                cout << "Gate error: Gate declared in an undefined polarization." << endl;
                return -1;
            }
        }else{
            P(i)   = -1;
        }
    }

    //Conversion to nlevelxnlevel matrix
    oelement.resize(nlevel,nlevel);
    oelement=matc::Identity(nlevel,nlevel);
    for(k=0;k<nbmch;k++){
        ch1=CH(k);
        km1=P(k);
        for(l=0;l<nbmch;l++){
            ch2=CH(l);
            km2=P(l);
            if(iodef.rows()>1){
                for(ks=0;ks<ns;ks++){
                    i=i_idx[ch1][km1][ks];
                    j=i_idx[ch2][km2][ks];
                    oelement(i,j)=U(k,l);
                }
            }else{
                for(km=0;km<nm;km++){
                for(ks=0;ks<ns;ks++){
                    i=i_idx[ch1][km][ks];
                    j=i_idx[ch2][km][ks];
                    oelement(i,j)=U(k,l);
                }}

            }
        }
    }

    // Update circuit with new element
    circmtx=oelement*circmtx;

    // Return success
    return 0;

}


//----------------------------------------
//
//  Adds a phase to the photons in a channel
//  depending on the optical path dt and the
//  frequency of the photon packets
//
//----------------------------------------
int qocircuit:: dispersion(int ch, double dt){
//  int ch;           // Channel
//  double dt;        // Time
//  Variables
    int    km;        // Mode of the channel i_ch
    int    ks;        // Wavepacket index
    int    iw;        // Frequency index
    double w;         // Frequency
    matc oelement;    // Extended nlevelxnlevel custom gate definition
//  Auxiliary index.
    int    i;         // Aux index
    int    j;         // Aux index


    // Check
    if(ch>=nch){
        cout << "dispersion error: Dispersion declared in an undefined channel." << endl;
        return -1;
    }

    if(emiss==0){
        cout << "dispersion error: No emitter set therefore no photon packets information available to compute the phase." << endl;
        return -1;
    }

    //Conversion to nlevelxnlevel matrix
    oelement.resize(nlevel,nlevel);
    oelement=matc::Identity(nlevel,nlevel);

    for(km=0;km<nm;km++){
    for(ks=0;ks<npack;ks++){
        i=i_idx[ch][km][ks];
        j=i_idx[ch][km][ks];
        iw=emitted->pack_def(1,ks);
        w=emitted->freq(0,iw);
        oelement(i,j)=exp(jm*dt*w);
    }}

    // Update circuit with new element
    circmtx=oelement*circmtx;

    // Return success
    return 0;
}


//----------------------------------------
//
//  Adds a phase to the photons in a channel
//  depending on the optical path dt and the
//  frequency of the photon packets for the
//  photons with a selected polarization
//
//----------------------------------------
int qocircuit:: dispersion(int ch, int P, double dt){
//  int ch;           // Channel
//  in P;             // Photon polarization required to apply a phase shift
//  double dt;        // Time
//  Variables
    int    ks;        // Wavepacket index
    int    iw;        // Frequency index
    double w;         // Frequency
    matc oelement;    // Extended nlevelxnlevel custom gate definition
//  Auxiliary index.
    int    i;         // Aux index
    int    j;         // Aux index


    // Check
    if(ch>=nch){
        cout << "dispersion error (Polarized): Dispersion declared in an undefined channel." << endl;
        return -1;
    }

    if(emiss==0){
        cout << "dispersion error (Polarized): No emitter set therefore no photon packets information available to compute the phase." << endl;
        return -1;
    }

    //Conversion to nlevelxnlevel matrix
    oelement.resize(nlevel,nlevel);
    oelement=matc::Identity(nlevel,nlevel);


    for(ks=0;ks<ns;ks++){
        i=i_idx[ch][P][ks];
        j=i_idx[ch][P][ks];
        iw=emitted->pack_def(1,ks);
        w=emitted->freq(0,iw);
        oelement(i,j)=exp(jm*dt*w);
    }

    // Update circuit with new element
    circmtx=oelement*circmtx;

    // Return success
    return 0;
}


//----------------------------------------
//
//  Adds a phase shifter to the circuit
//
//----------------------------------------
int qocircuit:: phase_shifter(int i_ch, double d_phi){
//  int    i_ch               // Phase shifter input channel
//  double d_phi;             // Angle phi in degrees
//  Variables
    double phi;               // Aux index.


    // Conversion to radians
    phi=d_phi*pi/180.0;
    // Call to general phase shifter
    return phase_shifter(i_ch, exp(jm*phi));

}


//----------------------------------------
//
// General phase shifter with losses
//
//----------------------------------------
int qocircuit:: phase_shifter(int i_ch, cmplx t){
//  int   i_ch                // Phase shifter input channel
//  cmplx t                   // Phase shifter transmission amplitude of probability.
//  Constants
    const  int nbmch=1;       // U square matrix row/column dimension
//  Variables
    matc   U(nbmch,nbmch);    // Waveplate 2x2 matrix definitions
    mati   V(1,nbmch);        // Channels to which U refers

    // Matrix
    U(0,0)=t;

    // Channel
    V(0,0)=i_ch;

    // Compute gate
    return custom_gate(V,U);
}



//----------------------------------------
//
//  Adds a polarized phase shifter to the circuit
//
//----------------------------------------
int qocircuit:: pol_phase_shifter(int i_ch, int P, double d_phi){
//  int    i_ch;              // Phase shifter input channel
//  int    P;                 // Polarization to which the phase shifter is sensitive
//  double d_phi;             // Angle phi in degrees
//  Variables
    double phi;               // Aux index.


    // Conversion to radians
    phi=d_phi*pi/180.0;
    // Call to general phase shifter
    return pol_phase_shifter(i_ch, P, exp(jm*phi));

}

//----------------------------------------
//
// General polarized phase shifter with losses
//
//----------------------------------------
int qocircuit:: pol_phase_shifter(int i_ch, int P, cmplx t){
//  int   i_ch;               // Phase shifter input channel
//  int   P;                  // Polarization to which the phase shifter is sensitive
//  cmplx t;                  // Phase shifter transmission amplitude of probability.
//  Constants
    const  int nbmch=1;       // U square matrix row/column dimension
//  Variables
    matc   U(nbmch,nbmch);    // Waveplate 2x2 matrix definitions
    mati   V(2,nbmch);        // Channels to which U refers

    // Matrix
    U(0,0)=t;
    // Channel
    V(0,0)=i_ch;
    V(1,0)=P;

    // Compute gate
    return custom_gate(V,U);
}


//----------------------------------------
//
// Polarization filter
//
//----------------------------------------
int qocircuit:: pol_filter(int i_ch, int P){
//  int i_ch;     // Polarization filter input channel
//  int P;        // Polarization to be removed/filtered


    // Call to general phase shifter
    return pol_phase_shifter(i_ch, P, (cmplx)0.0);

}


//----------------------------------------
//
// Lossy medium.
// It is an alias to a general phase shifter
// but with a double as input parameter.
//
//----------------------------------------
int qocircuit:: loss(int i_ch, double l){
//  int    i_ch;          // Phase shifter input channel
//  double l;             // Loss probability


    return phase_shifter(i_ch, (cmplx)sqrt(1.0-l));
}


//--------------------------------------------------
//
// Emitter configured using a coupling model
// We provide to the routine the overlapping
// coefficients of the wave packets non orthonormal
// base.
//
//--------------------------------------------------
void qocircuit:: emitter(matc D){
//  matc D                // Overlapping coefficient matrix for the initial states.
//  Variables
    int  ch;              // Channel evaluated.
    int  m;               // Mode of the channels.
    matc DXT;             // Extended overlapping coefficient matrix to include periods.
    matc T;               // Gram-Schmidt coefficient matrix for all states (including after delays)
    matc aux;             // Auxiliary matrix.
    matc oelement;        // Extended nlevelxnlevel beamsplitter definition
//  Auxiliary index.
    int  i;               // Mode 1 matrix aux index
    int  j;               // Mode 2 matrix aux index
    int  k;               // Mode 1 "Time" aux index
    int  l;               // Mode 2 "Time" aux index


    DXT=matc::Identity(ns,ns);
    for(k=0;k<ns;k=k+nsp){
    for(i=0;i<D.rows();i++){
    for(j=0;j<D.cols();j++){
        DXT(i+k,j+k)=D(i,j);

    }}}

    // Create/adapt time base.
    // We may have more times reserved
    // than the ones we are defined in D.
    // In those cases the visibility is the
    // identity.
    aux=GSP(DXT);
    prnt_dmat=GSP(D); // Information stored just for printing.
    confidence=mat_confidence(aux);
    T.resize(ns,ns);
    T=matc::Identity(ns,ns);
    for(i=0;i<aux.rows();i++){
    for(j=0;j<aux.cols();j++){
        T(i,j)=aux(i,j);

    }}


    //Calculate circuit element
    oelement.resize(nlevel,nlevel);
    oelement=matc::Identity(nlevel,nlevel);
    for(ch=0;ch<nch;ch++){
    for(m=0;m<nm;m++){
        for(k=0;k<ns;k++){
        for(l=0;l<ns;l++){
            i=i_idx[ch][m][k];
            j=i_idx[ch][m][l];

            oelement(i,j)=T(l,k);
        }}
    }}

    // Note that in this case is a non reversible operation
    emiss=1;
    init_dmat=oelement;
}


//--------------------------------------------------
//
// Emitter configured using a packet model.
//
//--------------------------------------------------
int qocircuit:: emitter(photon_mdl *mdl){
//  photon_mdl *mdl; // Photon model/definition
//  Variables
    int    nP;       // Number of packets
    matc   c;        // Coupling coefficients between packets
//  Packet distribution variables
    double ti;       // Packet i time
    double tj;       // Packet j time
    double wi;       // Packet i frequency
    double wj;       // Packet j frequency
    double dwi;      // Packet i frequency width or characteristic time
    double dwj;      // Packet j frequency width or characteristic time
//  Model variables
    matd   P;        // Photon packet matrix (calculated from the definition)
//  Auxiliary index
    int    i;        // Aux index
    int    j;        // Aux index


    //Calculate packet matrix
    P=mdl->create_packet_mtx();
    nP=P.cols();
    if(nP>nsp){
        cout << "Emitter(photon model) error: Not enough ns degrees of freedom! Needed at least: "<< nP << endl;
        return -1;
    }


    // Calculate coupling coefficients
    c=matc::Identity(nsp,nsp);
    for(i=0;i<nP;i++){
        for(j=0;j<nP;j++){
            //<ti|tj>
            ti=  P(0,i);
            tj=  P(0,j);
            wi=  P(1,i);
            wj=  P(1,j);
            dwi= P(2,i);
            dwj= P(2,j);


            if(mdl->kind==0){
                c(i,j)=gauss_coup(ti,wi,dwi,tj,wj,dwj);
            }else{
                c(i,j)=exp_coup(ti,wi,dwi,tj,wj,dwj);
            }
        }
    }

    // Set up the emitter using the calculated Gram-Schmidt coefficients
    emitter(c);

    // Return success
    return 0;
}


//--------------------------------------------------
//
// Emitter configured using a packet description in
// in a matrix.
//
//--------------------------------------------------
veci qocircuit::emitter(int npack, matd packets){
//  int  npack;        // Number of packets defined.
//  matd packets;      // Matrix that contains the definition of the packets.
//  Variables
    veci f_conversion; // Returns the new packet numbers if they are re-arranged
    veci conversion;   // Returns the new packet numbers if they are re-arranged
    vecd longtimes;    // List with the times defined
    matd longfreq;     // List with the frequencies defined
    vecd times;        // List with the times to create the photon model. Same size as entries
    matd freq;         // List with the frequencies to create the photon model. Same number of columns as entries
    mati defs;         // List with packet definitions. (Equal to pack unless timed=1)
    mati pack;         // Packet definition matrix to create the photon model.
    int  ntimes;       // Number of times in longtimes
    int  nfreq;        // Number of frequencies in long freq
//  Aux index
    int i;             // Aux index
    int j;             // Aux index
    int k;             // Aux index
    int ip;            // Aux index


    //Check number of packets
    if(npack>nsp){
        cout << "Emitter(npack) error: Not enough ns degrees of freedom! Needed: "<< npack << endl;
        return conversion;
    }

    // Perform the format conversion
    ntimes=0;
    nfreq=0;
    longtimes.resize(npack);
    longfreq.resize(2,npack);
    defs.resize(3,npack);
    for(i=0;i<npack;i++){
        ip=packets(0,i);

        j=0;
        while((j<ntimes)&&(abs(packets(1,i)-longtimes(j))>xcut)) j++;
        if (j==ntimes){
            longtimes(ntimes)=packets(1,i);
            ntimes=ntimes+1;
        }

        k=0;
        while((k<nfreq)&&((abs(packets(2,i)-longfreq(0,k))>xcut)
                       || (abs(packets(3,i)-longfreq(1,k))>xcut))) k++;
        if (k==nfreq){
            longfreq(0,nfreq)=packets(2,i);
            longfreq(1,nfreq)=packets(3,i);
            nfreq=nfreq+1;
        }

        defs(0,ip)=(int) packets(0,i);
        defs(1,ip)=j;
        defs(2,ip)=k;
    }

    // If timed create extra packets and their translations
    conversion.resize(npack);
    if(timed==1){
        pack.resize(3,ntimes*nfreq);
        k=0;
        for(i=0;i<ntimes;i++){
            for(j=0;j<nfreq;j++){
                pack(0,k)=k;
                pack(1,k)=i;
                pack(2,k)=j;
                k=k+1;
            }
        }
        for(i=0;i<npack;i++){
                j=0;
                while((pack(1,j)!=defs(1,i))||(pack(2,j)!=defs(2,i))) j++;
                conversion(defs(0,i))=pack(0,j);
        }

    }else{
        pack=defs;
        for(i=0;i<npack;i++) conversion(defs(0,i))=defs(0,i); // It can be done easier but this way is clearer.
    }

//  Adjust sizes
    times.resize(ntimes);
    freq.resize(2,nfreq);
    for(i=0;i<ntimes;i++){
        times(i)=longtimes(i);
    }

    for(i=0;i<nfreq;i++){
        freq(0,i) = longfreq(0,i);
        freq(1,i) = longfreq(1,i);
    }

    // Create the photon model
    delete emitted;
    emitted= new photon_mdl(pack, times, freq, ckind);

    // Set up the emitter using the photon model
    emitter(emitted);

    // Extend the conversion to all periods.
    f_conversion.resize(ns);
    for(i=0;i<ns;i++) f_conversion(i)=i;
    for(i=0;i<ns;i++){
        j=i%nsp;
        k=i/nsp;
        if(j<npack) f_conversion(i)=conversion(j)+k*nsp;
    }

    // Return converted packets
    return f_conversion;
}


//--------------------------------------------------
//
// Define wave packet
//
//--------------------------------------------------
int qocircuit::def_packet(int n, double t, double f, double w){
//  int    n;    // Number of packet
//  double t;    // Emission time
//  double f;    // Emission frequency
//  double w;    // Emission width/or decay length depending on packet shape model.
//  double tc;   // Characteristic width/time of the previous photons in the cascade.
//  int cascade; // Is this packet part of a cascade (extra random emission time).
//  Variables
    int    ip;   // Period index
    double rt;   // Reduced time of the packet to its equivalent in the first period.

    if(n>=nsp){
        cout << "def_packet error!: n > nd. We need to initialize the circuit with more packets." << endl;
        return -1;
    }
    if(emiss==1){
        cout << "def_packet error!: Emitter already set!" << endl;
        return -1;
    }

    // Calculate period of the packet and the reduced
    // time in the leading period.
    if(np>1){
        ip=floor((t+0.5*dtp)/dtp);
        rt=t-ip*dtp;
    }else{
        ip=0;
        rt=t;
    }

    // Adds a new entry to the matrix of definitions of packets
    pack_list(0,npack)=(double) n;
    pack_list(1,npack)=        rt;
    pack_list(2,npack)=         f;
    pack_list(3,npack)=         w;
    npack=npack+1;

    return n+ip*nsp;
}


//--------------------------------------------------
//
// Define an emitter
//
//--------------------------------------------------
veci qocircuit:: emitter(){
//  Variables
    veci conversion;    // Returns the new packer numbers if they are re-arranged


    conversion=emitter(npack,pack_list);
    return conversion;
}


//--------------------------------------------------
//
// Adds a delay
//
//--------------------------------------------------
int qocircuit:: delay(int ch){
//  int    ch;   // Channel to be delayed
//  Variables


    // Check
    if(ch>=nch){
        cout << "delay error: Delay declared in an undefined channel." << endl;
        return -1;
    }

    // Delay one period.
    delay(ch,+1);

    // Compute the phase shift due the delay
    return dispersion(ch,dtp);
}


//--------------------------------------------------
//
// Adds a delay of n periods
//
//--------------------------------------------------
int qocircuit:: delay(int i_ch, int ip){
//  int  i_ch            // Channel where the delay is applied
//  int  ip              // Number of periods to be delayed
//  Variables
    int  m;              // Mode of the channels.
    matc invT;           // Inverse of matrix T (T^{-1})
    matc oelement;       // Extended nlevelxnlevel beamsplitter definition
//  Auxiliary index.
    int  i;              // Channel 1 matrix index
    int  j;              // Channel 2 matrix index
    int  k;              // Aux index "Times"
    int  l;              // Aux index "Times"


    // Check
    if(i_ch>=nch){
        cout << "delay error (T matrix): Delay declared in an undefined channel." << endl;
        return -1;
    }

    // Reserve memory
    oelement.resize(nlevel,nlevel);
    oelement=matc::Identity(nlevel,nlevel);

    for(m=0;m<nm;m++){
        // Initialize
        // We do not want identity for certain cases
        for(k=0;k<ns;k++){
        for(l=0;l<ns;l++){

            i=i_idx[i_ch][m][k];
            j=i_idx[i_ch][m][l];
            if(k==l+ip*nsp) oelement(i,j)=1.0;
            else oelement(i,j)=0.0;
        }}

    }

    // In this case is a non reversible operation
    circmtx=oelement*circmtx;

    // Return success
    return 0;
}


//--------------------------------------------------
//
//  Adds a virtual circuit element to define a detector
//
//--------------------------------------------------
int qocircuit:: detector(int i_ch, int cond, int pol, int mpi, int mpo, double eff, double blnk, double gamma){
//  int    i_ch;        // Chanel where detection is performed
//  int    cond;        // Post selection condition. If cond<0 there is none.
//  int    pol;         // Post selection polarization condition. If pol<0 there is none.
//  int    mpi;         // Initial detection period (if -1 takes the first one as default)
//  int    mpo;         // Last detection period (if -1 takes the last one as default)
//  double eff;         // Efficiency of the detector
//  double blnk;        // Fraction of time blinking
//  double gamma;       // Dark counts rate
//  Variables
    int nclose;         // Total number of channels


    // Adds a new detector definition entry to the corresponding matrices
    // If cond>=0 Adds a new entry to the matrix of post-selection condition.
    ndetc=ndetc+1;

    // Check configuration errors
    if(ndetc>nch){
        cout << "detector error: More detectors than channels are being declared." << endl;
        return -1;
    }
    if((i_ch>=nch)||(i_ch<0)){
        cout << "detector error: This channel does not exist." << endl;
        return -1;
    }

    // Update post-selection condition
    if(cond>=0){
        det_def(0,ncond)=i_ch;
        det_def(1,ncond)=cond;
        det_def(2,ncond)=pol;
        ncond=ncond+1;
    }
    // If cond==-1 There is no condition
    // If cond==-2 The channel is ignored in the output
    if(cond==-2){
        ch_ignored(nignored)=i_ch;
        nignored=nignored+1;
    }

    // Set up window of detection by channel
    det_win(0,i_ch)=mpi;
    det_win(1,i_ch)=mpo;

    // Adds a new entry to the matrix of detector parameters
    det_par(0,i_ch)=blnk;
    det_par(1,i_ch)=gamma;

    // Put losses in circuit to consider the efficiency
    if(eff<(1.0-xcut)){
        phase_shifter(i_ch,(cmplx)sqrt(eff));
    }

    // Last detector operations
    // The total number of channels to put a detector
    // is different depending if we are computing losses
    if(losses==0) nclose=nch;
    else nclose=nch/2;

    // If this is the last detector compute losses if needed.
    if((losses==1)&&(ndetc==nclose)) compute_losses();

    // If this is the last detector add the emitter matrix
    // at the beginning. (Otherwise computing the losses breaks the calculation)
    if((emiss==1)&&(ndetc==nclose)) circmtx=circmtx*init_dmat;

    // Return success
    return 0;
}


//--------------------------------------------------
//
//  Remaining not defined detectors.
//
//--------------------------------------------------
int qocircuit::remdec(){
//  Variables
    int nclose;
    int remaining;

    // Last detector operations
    // The total number of channels to put a detector
    // is different depending if we are computing losses
    if(losses==0) nclose=nch;
    else nclose=nch/2;

    remaining=nclose-ndetc;

    return remaining;
}


//--------------------------------------------------
//
//  Adds a virtual circuit element to define a detector
//
//--------------------------------------------------
int qocircuit:: detector(int i_ch){
//  int    i_ch;        // Chanel where detection is performed


    return detector(i_ch,-1,-1,-1,-1,1.0,0.0,0.0);
}


//--------------------------------------------------
//
//  Adds a virtual circuit element to define a detector
//
//--------------------------------------------------
int qocircuit:: detector(int i_ch, int cond){
//  int    i_ch;        // Chanel where detection is performed
//  int    cond;        // Post selection condition. If cond<0 there is none.


    return detector(i_ch,cond,-1,-1,-1,1.0,0.0,0.0);
}


//--------------------------------------------------
//
//  Adds a virtual circuit element to define a detector
//
//--------------------------------------------------
int qocircuit:: detector(int i_ch, int cond, double eff, double blnk, double gamma){
//  int    i_ch;        // Chanel where detection is performed
//  int    cond;        // Post selection condition. If cond<0 there is none.
//  double eff;         // Efficiency of the detector
//  double blnk;        // Fraction of time blinking
//  double gamma;       // Dark counts rate


    return detector(i_ch,cond,-1,-1,-1,eff,blnk,gamma);
}

//--------------------------------------------------
//
//  Adds a virtual circuit element to flag a channel to be ignored.
//
//--------------------------------------------------
int qocircuit:: ignore(int i_ch){
//  int    i_ch;        // Chanel where detection is performed


    return detector(i_ch,-2,-1,-1,-1,1.0,0.0,0.0);
}


//--------------------------------------------------
//
//  Adds a random Gaussian noise to the result.
//  This makes for different post-processing of information
//  in the p_bin class.
//
//--------------------------------------------------
void qocircuit::noise(double stdev2){
//  double stdev2;     // Dispersion of the noise


    dev=stdev2;
}


//----------------------------------------
//
//  Print circuit matrix
//
//----------------------------------------
void qocircuit:: prnt(int format){
//  int format     // Format of the output
//  Variables
    int ch1;       // Channel 1
    int ch2;       // Channel 2
    int m1;        // Mode 1
    int m2;        // Mode 2
    int s1;        // Packet 1
    int s2;        // Packet 2
    int firstline; // Is this the first line? 1=Yes/0=No
//  Auxiliary index
    int i;         // Aux index
    int j;         // Aux index


    // For all possible levels.
    for(i=0;i<nlevel;i++){
        firstline=1;
        ch1=idx[i].ch;
        m1=idx[i].m;
        s1=idx[i].s;

        //Print LHS
        cout<< i <<" :| " << ch1;
        if(format==0){
            if(nm>1) cout << ", "<< m1;
            if(ns>1) cout << ", "<< s1;
        }else{
            if(nm>1) cout << ", "<< pl[m1];
            if(ns>1) cout << ", "<< s1;
        }
        cout << " > -> ";

        // Print RHS. Local transformation rule
        for(j=0;j<nlevel;j++){
            ch2=idx[j].ch;
            m2=idx[j].m;
            s2=idx[j].s;

            if(abs(circmtx(j,i))>xcut){
                if(firstline==0){
                    cout << " + ";
                }
                firstline=0;
                cout << circmtx(j,i) << " * | "<< ch2;
                if(format==0){
                    if(nm>1) cout << ", "<< m2;
                    if(ns>1) cout << ", "<< s2;
                }else{
                    if(nm>1) cout << ", "<< pl[m2];
                    if(ns>1) cout << ", "<< s2;
                }
                cout << " >";
            }
            if(j==(nlevel-1)){
                if(firstline==1) cout << "| vac >";
                cout << endl;
            }

        }

    }
}


//----------------------------------------
//
//  Print Gram-Schmidt coefficients
//
//----------------------------------------
void qocircuit:: prntGS(){

    cout << "Gran Schmidt coefficients: " << endl;
    cout << prnt_dmat << endl;

}


//----------------------------------------
//
//  Print the visibility/overlapping probability
//  between two wavepackets of the circuit photon model.
//
//----------------------------------------
double qocircuit:: emitted_vis(int i,int j){
//  int i;      // Wave packet index 1
//  int j;      // Wave packet index 2

    return emitted->visibility(i,j,nsp);
}


//----------------------------------------
//
//  Creates an empty photon model.
//
//----------------------------------------
photon_mdl::photon_mdl(){
    // Do nothing
}


//----------------------------------------
//
//  Creates a photon model from parameters
//
//----------------------------------------
photon_mdl::photon_mdl(mati i_pack_def, vecd i_times, matd i_freq, char ckind){
//  mati pack_def       // Packet definition
//  matd times          // Times and phase times for each index/enumeration
//  matd freq           // Frequency definition for each index/enumeration
//  char ckind          // Kind of wave packets 'G': Gaussian / 'E': Exponential


    // Configure the distribution of the packets
    pack_def=create_packet_idx(i_pack_def);
    times=i_times;
    freq=i_freq;
    kind=0;
    if(ckind=='E') kind=1;
}


//----------------------------------------
//
//  Destroys photon model.
//
//----------------------------------------
photon_mdl::~photon_mdl(){
    // Do nothing.
}


//----------------------------------------
//
//  Creates a copy of a photon model.
//
//----------------------------------------
photon_mdl *photon_mdl::clone(){
//  Variables
    photon_mdl *new_mdl;    // New photon model recipient of the copy


    new_mdl=new photon_mdl();
    new_mdl->pack_def=pack_def;
    new_mdl->times=times;
    new_mdl->freq=freq;
    new_mdl->kind=kind;

    return new_mdl;
}


//---------------------------------------------------------------------------
//
//  Creates a packet matrix model from the packet definitions/specifications
//
//----------------------------------------------------------------------------
matd photon_mdl:: create_packet_mtx(){
//  Variables
    matd P;             // Packet matrix model
//  Auxiliary index
    int  i;             // Aux index


    //Reserve memory
    P.resize(3,pack_def.cols());

    // Create model
    for(i=0;i<pack_def.cols();i++){
            P(0,i)=times(pack_def(0,i));
            P(1,i)=freq(0,pack_def(1,i));
            P(2,i)=freq(1,pack_def(1,i));
    }

    //Returns the matrix model
    return P;
}


//----------------------------------------
//
//  Probability of two wave packets to overlap
//
//----------------------------------------
double photon_mdl:: visibility(int i,int j, int nsp){
//  int    i         // Packet i number
//  int    j         // Packet j number
//  int    nsp       // Number of packets by period
//  Variables
    cmplx  c;        // Coupling between two packets
//  Packet distribution variables
    int ri;          // Reduced index of packet i
    int rj;          // Rediced index of packet j
    int pi;          // Period of packet i
    int pj;          // Period of packet i
    double ti;       // Packet i time
    double tj;       // Packet j time
    double wi;       // Packet i frequency
    double wj;       // Packet j frequency
    double dwi;      // Packet i frequency width or characteristic time
    double dwj;      // Packet j frequency width or characteristic time
//  Model variables
    matd   P;        // Photon packet matrix (calculated from the definition)


    // Calculate the packet matrix from the photon model
    P=create_packet_mtx();

    // Obtain period and packet number
    pi=i/nsp;
    pj=j/nsp;
    ri=i%nsp;
    rj=j%nsp;

    // Check  validity
    if(pi!=pj) return 0.0;
    if((ri>=P.cols())||(rj>=P.cols())){
            cout << "Visibility error: One of the demanded packets does not exist!" << endl;
            return 0.0;
    }

    // Obtain the parameters from the packet matrix
    ti=  P(0,ri);
    tj=  P(0,rj);
    wi=  P(1,ri);
    wj=  P(1,rj);
    dwi= P(2,ri);
    dwj= P(2,rj);

    // Calculate the coupling
    if(kind==0){
        c=gauss_coup(ti,wi,dwi,tj,wj,dwj);
    }else{
        c=exp_coup(ti,wi,dwi,tj,wj,dwj);
    }

    // Return visibility
    return abs(conj(c)*c);
}


//----------------------------------------
//
//  Returns the packet definition in a 3-row format
//
//----------------------------------------
mati photon_mdl:: return_packet_def(){
    // Variables
    mati aux;           //Auxiliary matrix for the 3-row packet definition
    // Auxiliary index
    int  i;             // Aux index


    // Reserve memory
    aux.resize(3,pack_def.cols());

    // Rewrite packet definition
    for(i=0;i<pack_def.cols();i++){
        aux(0,i)=i;
        aux(1,i)=pack_def(0,i);
        aux(2,i)=pack_def(1,i);
    }

    // Return packet definition
    return aux;
}


//----------------------------------------
//
//  Prints the the table of times of the
//  photon model.
//
//----------------------------------------
void photon_mdl::prnt_times(){
// Auxiliary index
    int i;      // Aux index


    for(i=0;i<times.size();i++) cout << setw(2) << right << i << ": " << setw(10) << right << setprecision(5)  << times(i) << endl;
    cout << endl;
}


//----------------------------------------
//
//  Prints the the table of frequencies of the
//  photon model.
//
//----------------------------------------
void photon_mdl::prnt_freqs(){
// Auxiliary index
    int i;      // Aux index


    for(i=0;i<freq.cols();i++) cout << setw(2) << right  << i << ": " << setw(15) << right << setprecision(8) << freq(0,i) << "   "<< right  << setprecision(8) << freq(1,i) << endl;
    cout << endl;
}


//----------------------------------------
//
//  Prints the the table of packet index of the
//  photon model.
//
//----------------------------------------
void photon_mdl::prnt_packets(){
// Auxiliary index
    int i;      // Aux index


    for(i=0;i<pack_def.cols();i++) cout << setw(3) << i;
    cout << endl;
    for(i=0;i<pack_def.cols();i++) cout << setw(3) <<  pack_def(0,i);
    cout << endl;
    for(i=0;i<pack_def.cols();i++) cout << setw(3) <<  pack_def(1,i);
    cout << endl;
    cout << endl;
}


//----------------------------------------
//
//  Prints the information related with the
//  photon model.
//
//----------------------------------------
void photon_mdl::prnt(){


    cout << "----------------------------------------------------------------------------------" << endl;
    cout << "Table of times:" << endl;
    prnt_times();
    cout << "Table of frequencies:" << endl;
    prnt_freqs();
    cout << "Table of packets:" << endl;
    prnt_packets();
    cout << endl;

    if(kind==0) cout << "Kind of packet: Gaussian" << endl;
    else        cout << "Kind of packet: Exponential" << endl;
    cout << "---------------------------------------------------------------------------------" << endl;
    cout << endl;
}


//---------------------------------------------------------------------------
//
//  Auxiliary function that converts an unordered packet definition into an
//  ordered definition that consists in an indexed matrix where the column
//  index are the wave-packet numbers. It switches between two different
//  representations.
//
//----------------------------------------------------------------------------
mati create_packet_idx(mati pack_def){
//  mati pack_def    //Packet definition
//  Variables
    mati pack_idx;   //Packet index
//  Auxiliary index
    int  i;          // Aux index


    //Reserve memory
    pack_idx.resize(2,pack_def.cols());
    //Transform representation
    for(i=0;i<pack_def.cols();i++){
        pack_idx(0,pack_def(0,i))=pack_def(1,i);
        pack_idx(1,pack_def(0,i))=pack_def(2,i);
    }
    // Return matrix
    return pack_idx;
}

