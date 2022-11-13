//======================================================================================================
// File qocircuit.cpp
//
// OPTICAL CIRCUIT LIBRARY.
//
// Copyright © 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved.
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


    create_circuit(i_nch, 1, 1, 0, 0, false,'G');
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


    create_circuit(i_nch, i_nm, i_ns, 0 ,0, false, 'G');
}


//----------------------------------------
//
//  Create circuit
//
//----------------------------------------
qocircuit::qocircuit(int i_nch, int i_nm, int i_ns,int clock, char i_ckind){
//  int  i_nch;   // Number of channels
//  int  i_nm;    // Number of modes
//  int  i_ns;    // Number of wavepackets
//  int  clock;   // Kind of detectors
//  char i_ckind; // Kind of packets 'G'=Gaussian/'E'=Exponential


    create_circuit(i_nch, i_nm, i_ns, clock ,0, false, i_ckind);
}


//----------------------------------------
//
//  Create circuit
//
//----------------------------------------
qocircuit::qocircuit(int i_nch, int i_nm, int i_ns,int clock, int i_R, bool loss,char i_ckind){
//  int  i_nch;   // Number of channels
//  int  i_nm;    // Number of modes
//  int  i_ns;    // Number of wavepackets
//  int  clock;   // Kind of detectors
//  int  i_R;     // Number of iterations to calculate blinking and dark counts
//  bool loss;    // Are losses going to be calculated explicitly? True=Yes/False=No
//  char i_ckind; // Kind of packets 'G'=Gaussian/'E'=Exponential


    create_circuit(i_nch, i_nm, i_ns, clock, i_R, loss, i_ckind);
}


//----------------------------------------
//
//  Auxiliary private method to create a circuit
//
//----------------------------------------
void qocircuit::create_circuit(int i_nch, int i_nm, int i_ns,int clock, int i_R, bool loss, char i_ckind){
//  int  i_nch;   // Number of channels
//  int  i_nm;    // Number of modes
//  int  i_ns;    // Number of wavepackets
//  int  clock;   // Kind of detectors
//  int  i_R;     // Number of iterations to calculate blinking and dark counts
//  bool loss;    // Are losses going to be calculated explicitly? True=Yes/False=No
//  char i_ckind; // Kind of packets 'G'=Gaussian/'E'=Exponential
//  Auxiliary index
    int i;        // Aux index
    int j;        // Aux index
    int k;        // Aux index
    int l;        // Aux index


    // Initialize degrees of freedom
    nm=i_nm;
    ns=i_ns;
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
    det_def.resize(2,i_nch);
    det_par.resize(2,i_nch);
    ch_ignored.resize(i_nch);


    // Initialize emitter/Packet definitions
    ckind=i_ckind;
    npack=0;
    pack_list.resize(5,ns);
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
    newcircuit= new qocircuit(nch,nm,ns,timed,ckind);
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
void qocircuit:: beamsplitter(int i_ch1, int i_ch2, double d_theta, double d_phi){
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
    custom_gate(V,U);
}


//----------------------------------------
//
//  Adds a dielectric to the circuit
//
//----------------------------------------
void qocircuit:: dielectric(int i_ch1, int i_ch2, cmplx t, cmplx r){
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
    custom_gate(V,U);

}


//-------------------------------------------------------
//
//  Adds a virtual circuit element to compute the states
//  with less photons at the input than the output.
//
//-------------------------------------------------------
void qocircuit::compute_losses(){
//  Variables
    int  nchm;   // Half the number of channels
    int  ch;     // Channel
    int  m;      // Polarization
    int  s;      // Wavepacket
    matc D;      // Diagonal matrix S=U D U^*
    matc U;      // Unitary  matrix S=U D U^*
    ComplexEigenSolver<matc> es; // Eigensolver
//  Auxiliary index
    int  i;      // Aux index
    int  j;      // Aux index
    int  k;      // Aux index
    int  q;      // Aux index


    //Diagonalize
    nchm=nch/2;
    es.compute(circmtx);
    D = es.eigenvalues().asDiagonal();
    U= es.eigenvectors();

    // Compute losses eigenchannels
    for(i=0;i<nlevel;i++){
        // First we calculate the couplings
        ch=idx[i].ch;
        m=idx[i].m;
        s=idx[i].s;
        ch=ch+nchm;
        if(ch<nch){
            k=i_idx[ch][m][s];
            if((1.0-abs(conj(D(i,i))*D(i,i)))>0) D(i,k)=sqrt(1.0-abs(conj(D(i,i))*D(i,i)));
            else D(i,k)=0.0;
            D(k,i)=D(i,k);
            D(k,k)=0.0;

            // Second we write the proper unitary transformation
            for(j=0;j<nlevel;j++){
                ch=idx[j].ch;
                m=idx[j].m;
                s=idx[j].s;
                ch=ch+nchm;
                if(ch<nch){
                    q=i_idx[ch][m][s];
                    U(k,q)= U(i,j);
                }
            }
        }
    }

    // Undo diagonalization
    circmtx=U*D*U.adjoint();
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
void qocircuit:: NSX(int i_ch1, int i_ch2, int i_ch3){
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
    custom_gate(V,U);
}


//----------------------------------------
//
//  Adds a polarization rotator to the circuit
//
//----------------------------------------
void qocircuit:: rotator(int i_ch, double d_theta, double d_phi){
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
    custom_gate(W,U);
}


//----------------------------------------
//
//  Adds a polarized beamsplitter to the circuit
//
//----------------------------------------
void qocircuit:: polbeamsplitter(int i_ch1, int i_ch2, int pol){
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
    custom_gate(W,U);
}


//----------------------------------------
//
//  Adds a general waveplate to the circuit
//
//----------------------------------------
void qocircuit:: waveplate(int i_ch, double d_alpha, double d_gamma){
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
    custom_gate(W,U);
}


//----------------------------------------
//
//  Adds a half waveplate to the circuit
//
//----------------------------------------
void qocircuit:: half(int i_ch, double alpha){
//  int    i_ch               // Waveplate input channel
//  double d_alpha;           // Angle of the polarization rotation


    waveplate(i_ch,alpha,90.0);
}


//----------------------------------------
//
//  Adds a quarter waveplate to the circuit
//
//----------------------------------------
void qocircuit:: quarter(int i_ch, double alpha){
//  int    i_ch               // Waveplate input channel
//  double d_alpha;           // Angle of the polarization rotation

    waveplate(i_ch,alpha,45.0);
}


//----------------------------------------
//
//  Adds a 2x2 MMI to the circuit
//
//----------------------------------------
void qocircuit:: MMI2(int i_ch1, int i_ch2){
//  int   i_ch1               // MMI2 input channel 1
//  int   i_ch2               // MMI2 input channel 2


    dielectric(i_ch1,i_ch2,1/sqrt(2.0),jm/sqrt(2.0));
}

//------------------------------------------------
//
// Rewires / swaps two channels (of the same mode)
//
//------------------------------------------------
void qocircuit:: rewire(int i_ch1,int i_ch2){
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
    custom_gate(V,U);
}


//----------------------------------------
//
//  Adds a custom gate to the circuit
//
//----------------------------------------
void qocircuit:: custom_gate(mati iodef, matc U){
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
        if(iodef.rows()>1){
            P(i)   = iodef(1,i);
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
}


//----------------------------------------
//
//  Adds a phase to the photons in a channel
//  depending on the optical path dt and the
//  frequency of the photon packets
//
//----------------------------------------
void qocircuit:: dispersion(int ch, double dt){
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


    //Conversion to nlevelxnlevel matrix
    oelement.resize(nlevel,nlevel);
    oelement=matc::Identity(nlevel,nlevel);

    for(km=0;km<nm;km++){
    for(ks=0;ks<ns;ks++){
        i=i_idx[ch][km][ks];
        j=i_idx[ch][km][ks];
        iw=emitted->pack_def(1,ks);
        w=emitted->freq(0,iw);
        oelement(i,j)=exp(jm*dt*w);
    }}

    // Update circuit with new element
    circmtx=oelement*circmtx;
}


//----------------------------------------
//
//  Adds a phase to the photons in a channel
//  depending on the optical path dt and the
//  frequency of the photon packets for the
//  photons with a selected polarization
//
//----------------------------------------
void qocircuit:: dispersion(int ch, int P, double dt){
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
}


//----------------------------------------
//
//  Adds a phase shifter to the circuit
//
//----------------------------------------
void qocircuit:: phase_shifter(int i_ch, double d_phi){
//  int    i_ch               // Phase shifter input channel
//  double d_phi;             // Angle phi in degrees
//  Variables
    double phi;               // Aux index.


    // Conversion to radians
    phi=d_phi*pi/180.0;
    // Call to general phase shifter
    phase_shifter(i_ch, exp(jm*phi));

}


//----------------------------------------
//
// General phase shifter with losses
//
//----------------------------------------
void qocircuit:: phase_shifter(int i_ch, cmplx t){
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
    custom_gate(V,U);
}


//----------------------------------------
//
// Lossy medium.
// It is an alias to a general phase shifter
// but with a double as input parameter.
//
//----------------------------------------
void qocircuit:: loss(int i_ch, double l){
//  int    i_ch           // Phase shifter input channel
//  double l;             // Loss probability


    phase_shifter(i_ch, (cmplx)sqrt(1.0-l));
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
    int  ch;              // Channel evaluated 1
    int  m;               // Mode of the channels.
    matc T;               // Gram-Schmidt coefficient matrix for all states (including after delays)
    matc aux;             // Auxiliary matrix.
    matc oelement;        // Extended nlevelxnlevel beamsplitter definition
//  Auxiliary index.
    int  i;               // Mode 1 matrix aux index
    int  j;               // Mode 2 matrix aux index
    int  k;               // Mode 1 "Time" aux index
    int  l;               // Mode 2 "Time" aux index


    // Create/adapt time base.
    // We may have more times reserved
    // than the ones we are defined in D.
    // In those cases the visibility is the
    // identity.
    aux=GSP(D);
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
    double tri;      // Packet i phase time
    double trj;      // Packet j phase time
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
    if(nP>ns){
        cout << "Emitter(photon model) error: Not enough ns degrees of freedom! Needed: "<< nP << endl;
        return -1;
    }


    // Calculate coupling coefficients
    c=matc::Identity(ns,ns);
    for(i=0;i<nP;i++){
        for(j=0;j<nP;j++){
            //<ti|tj>
            ti=  P(0,i);
            tj=  P(0,j);
            wi=  P(1,i);
            wj=  P(1,j);
            dwi= P(2,i);
            dwj= P(2,j);
            tri=  P(3,i);
            trj=  P(3,j);

            if(mdl->kind==0){
                c(i,j)=gauss_coup(ti,wi,dwi,tj,wj,dwj,tri,trj);
            }else{
                c(i,j)=exp_coup(ti,wi,dwi,tj,wj,dwj,tri,trj);
            }
        }
    }

    // Set up the emitter using the calculated Gram-Schmidt coefficients
    emitter(c);

    return 0;
}


//--------------------------------------------------
//
// Emitter configured using a packet description in
// in a matrix.
//
//--------------------------------------------------
veci qocircuit::emitter(int npack, matd packets){
//  int  npack;       // Number of packets defined.
//  matd packets;     // Matrix that contains the definition of the packets.
//  Variables
    veci conversion;  // Returns the new packet numbers if they are re-arranged
    matd longtimes;   // List with the times defined
    matd longfreq;    // List with the frequencies defined
    matd times;       // List with the times to create the photon model. Same size as entries
    matd freq;        // List with the frequencies to create the photon model. Same number of columns as entries
    mati defs;        // List with packet definitions. (Equal to pack unless timed=1)
    mati pack;        // Packet definition matrix to create the photon model.
    int  ntimes;      // Number of times in longtimes
    int  nfreq;       // Number of frequencies in long freq
//  Aux index
    int i;            // Aux index
    int j;            // Aux index
    int k;            // Aux index
    int ip;           // Aux index


    //Check number of packets
    if(npack>ns){
        cout << "Emitter(npack) error: Not enough ns degrees of freedom! Needed: "<< npack << endl;
        return conversion;
    }

    // Perform the format conversion
    ntimes=0;
    nfreq=0;
    longtimes.resize(2,npack);
    longfreq.resize(2,npack);
    defs.resize(3,npack);
    for(i=0;i<npack;i++){
        j=0;
        ip=packets(0,i);
        while((j<ntimes)&&(abs(packets(1,ip)-longtimes(0,j))>xcut)) j++;
        if (j==ntimes){
            longtimes(0,ntimes)=packets(1,ip);
            longtimes(1,ntimes)=packets(4,ip);
            ntimes=ntimes+1;
        }

        k=0;
        while((k<nfreq)&&((abs(packets(2,ip)-longfreq(0,k))>xcut)
                       || (abs(packets(3,ip)-longfreq(1,k))>xcut))) k++;
        if (k==nfreq){
            longfreq(0,nfreq)=packets(2,ip);
            longfreq(1,nfreq)=packets(3,ip);
            nfreq=nfreq+1;
        }

        defs(0,i)=(int) packets(0,ip);
        defs(1,i)=j;
        defs(2,i)=k;
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
                conversion(i)=j;
        }
    }else{
        pack=defs;
        for(i=0;i<npack;i++) conversion(i)=i;
    }

//  Adjust sizes
    times.resize(2,ntimes);
    freq.resize(2,nfreq);
    for(i=0;i<ntimes;i++){
        times(0,i)=longtimes(0,i);
        times(1,i)=longtimes(1,i);
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

    // Return converted packets
    return conversion;
}


//--------------------------------------------------
//
// Define wave packet
//
//--------------------------------------------------
void qocircuit::def_packet(int n, double t, double f, double w){
//  int    n;    // Number of packet
//  double t;    // Emission time
//  double f;    // Emission frequency
//  double w;    // Emission width/or decay length depending on packet shape model.
//  double phi;  // Initial phase

    this->def_packet(n, t, f, w, t, 0.0, 0.0 ,0);
}


//--------------------------------------------------
//
// Define wave packet
//
//--------------------------------------------------
cfg_rnd qocircuit::def_packet(int n, double t, double f, double w, double tp, double val1, double val2 ,int rand){
//  int    n;    // Number of packet
//  double t;    // Emission time
//  double f;    // Emission frequency
//  double w;    // Emission width/or decay length depending on packet shape model.
//  doulbe val1; // Value 1: Meaning depends on rand configuration
//  doulbe val2; // Value 2: Meaning depends on rand configuration
//  int rand;    // Phase computation mode
//  Variables
    double dt1;  // "Visibility time"
    double dt2;  // "Phase time"
    cfg_rnd dt;  // Configuration of the random times. Value to be returned


    if(n>=pack_list.cols()){
        cout << "def_packet error!: n > nd. We need to initialize the circuit with more packets." << endl;
        return dt;
    }

    // Random numbers generation
    // Gaussian   : In cascade configuration val1 dwxx and val2 is T2*
    // Exponential: In cascade configuration val1  txx and val2 is T2*
    if(ckind=='G'){ // Gaussian
        dt1=erfi(2*urand()-1)/val1;
        dt2=2*w*erfi(2*urand()-1)/(val2*val1);
    }else{         // Exponential
        dt1=val1*expi(urand());
        dt2=val1*val2*expi(urand())/(2.0*w); //(val1/w)*(T*/2)*u
    }

    // Adds a new entry to the matrix of definitions of packets
    pack_list(0,npack)=(double) n;
    pack_list(2,npack)=         f;
    pack_list(3,npack)=         w;
    switch(rand) {
        case 0: // No cascade
            pack_list(1,npack)= t;
            pack_list(4,npack)= tp;
            dt.dt1=val1;
            dt.dt2=0.0;
            break;

        case 1: // Cascade *fully* automatic
            pack_list(1,npack)= t+dt1;
            pack_list(4,npack)= tp+dt1;
            dt.dt1=dt1;
            dt.dt2=dt2;
            break;
        case 2: // Cascade automatic
            pack_list(1,npack)= t+dt1;
            pack_list(4,npack)= tp+dt2;
            dt.dt1=dt1;
            dt.dt2=dt2;
            break;
        case 3: // Cascade manual
            pack_list(1,npack)= t+val1;//+val;
            pack_list(4,npack)= tp+val2;
            dt.dt1=val1;
            dt.dt2=val2;
            break;
        case 4: // Cascade semi-manual
            pack_list(1,npack)= t;
            pack_list(4,npack)= tp+dt2;
            dt.dt1=0;
            dt.dt2=dt2;
            break;
        default:
            cout << "def_packet error: rand goes form 0 to 4. Invalid configuration." << endl;
            return dt;
            break;
    }
    npack=npack+1;

    return dt;
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
int qocircuit:: delay(int ch, double dt){
//  int    ch;   // Channel to be delayed
//  double dt;   // Amount of time to be delayed
//  Variables
    photon_mdl *new_mdl;  // New photon that consider the new packet definitions created by the delay.


    if(emiss==0){
        cout << "Delay error: No emitter set therefore no photon packets available to be delayed" << endl;
        return -1;
    }

    // Apply the delay and create the new photon model
    new_mdl=delay(ch,dt, emitted);

    // Update the photon model
    delete emitted;
    emitted=new_mdl;

    // Compute the phase shift due the delay
    dispersion(ch,dt);

    return 0;
}


//--------------------------------------------------
//
// Adds a delay for a given photon model
//
//--------------------------------------------------
photon_mdl *qocircuit:: delay(int ch, double dt, photon_mdl *mdl){
//  int    ch;       // Channel where the delay is performed
//  double dt;       // Delay time
//  photon_mdl *mdl; // Photon model/definition
//  Variables
    int    nP;       // Number of packets
    matc   c;        // Coupling coefficients between packets
    matc   U;        // Gram-Schmidt coefficient matrix
//  Packet distribution variables
    double ti;       // Packet i time
    double tj;       // Packet j time
    double tri;      // Packet i phase time
    double trj;      // Packet j phase time
    double wi;       // Packet i frequency
    double wj;       // Packet j frequency
    double dwi;      // Packet i frequency width or characteristic time
    double dwj;      // Packet j frequency width or characteristic time
//  Model variables
    matd   P;        // Photon packet model (calculated from the definition)
    photon_mdl *new_mdl;  // Photon new model
//  Auxiliary index
    int    i;        // Aux index
    int    j;        // Aux index


    //Update photon model/definition
    new_mdl=mdl->clone();
    new_mdl->update(dt);

    //Calculate new packet matrix
    P=new_mdl->create_packet_mtx();
    nP=P.cols();
    if(nP>ns){
        cout << "Delay(photon model) error: Not enough ns degrees of freedom! Needed: "<< nP << endl;
        return new_mdl;
    }

    // Calculate coupling coefficients
    // It is better to recalculate everything again
    c=matc::Identity(ns,ns);
    for(i=0;i<nP;i++){
        for(j=0;j<nP;j++){
            //<ti|tj>
            ti  = P(0,i);
            tj  = P(0,j);
            wi  = P(1,i);
            wj  = P(1,j);
            dwi = P(2,i);
            dwj = P(2,j);
            tri = P(3,i);
            trj = P(3,j);

            if(new_mdl->kind==0){
                c(i,j)=gauss_coup(ti,wi,dwi,tj,wj,dwj,tri,trj);
            }else{
                c(i,j)=exp_coup(ti,wi,dwi,tj,wj,dwj,tri,trj);
            }
        }
    }

    // Set up the Gram-Schmidt coefficients for all the packets
    U=GSP(c);
    confidence=mat_confidence(U);

    // Compute the delay in the circuit matrix
    delay(ch,U,nP/2);

    // Return the updated definitions
    return new_mdl;
}


//--------------------------------------------------
//
// Adds a delay computing and applying the delay gate
// using a visibility model.
// We provide to the routine the Gram-Schmidt base
// coefficients of the new times orthonormalized with
// the previous ones and themselves.
//
//--------------------------------------------------
void qocircuit:: delay(int i_ch, matc T, int nT){
//  int  i_ch            // Channel where the delay is applied
//  matc T               // Gram-Schmidt coefficient matrix for all states (including after delays)
//  int  nT              // Number of new wave packets
//  Variables
    int  m;              // Mode of the channels.
    int  it;             // Original input time
    int  ot;             // Delayed output time
    matc invT;           // Inverse of matrix T (T^{-1})
    matc oelement;       // Extended nlevelxnlevel beamsplitter definition
//  Auxiliary index.
    int  i;              // Channel 1 matrix index
    int  j;              // Channel 2 matrix index
    int  k;              // Aux index "Times"
    int  l;              // Aux index "Times"


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
            oelement(i,j)=0.0;
        }}

        // Calculate the inverse.
        invT=T.inverse();

        // Compute the gate
        for(it=0;it<nT;it++){
            ot=nT+it;

            for(k=0;k<ns;k++){
            for(l=0;l<ns;l++){
                i=i_idx[i_ch][m][k];
                j=i_idx[i_ch][m][l];

                 oelement(i,j)=oelement(i,j)+T(ot,k)*invT(l,it);
            }}
        }
    }

    // In this case is a non reversible operation
    circmtx=oelement*circmtx;
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
//  Variables
    int nclose;         // Total number of channels


    // Adds a new detector definition entry to the corresponding matrices
    // If cond>=0 Adds a new entry to the matrix of post-selection condition.
    ndetc=ndetc+1;
    if(ndetc>nch){
        cout << "Detector error: More detectors than channels are being declared." << endl;
        return -1;
    }
    if((i_ch>=nch)||(i_ch<0)){
        cout << "Detector error: This channel does not exist." << endl;
        return -1;
    }

    if(cond>=0){
        det_def(0,ncond)=i_ch;
        det_def(1,ncond)=cond;
        ncond=ncond+1;
    }
    // If cond==-1 There is no condition
    // If cond==-2 The channel is ignored in the output
    if(cond==-2){
        ch_ignored(nignored)=i_ch;
        nignored=nignored+1;
    }

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


    return detector(i_ch,-1,1.0,0.0,0.0);
}


//--------------------------------------------------
//
//  Adds a virtual circuit element to define a detector
//
//--------------------------------------------------
int qocircuit:: detector(int i_ch, int cond){
//  int    i_ch;        // Chanel where detection is performed
//  int    cond;        // Post selection condition. If cond<0 there is none.


    return detector(i_ch,cond,1.0,0.0,0.0);
}


//--------------------------------------------------
//
//  Adds a virtual circuit element to flag a channel to be ignored.
//
//--------------------------------------------------
int qocircuit:: ignore(int i_ch){
//  int    i_ch;        // Chanel where detection is performed


    return detector(i_ch,-2,1.0,0.0,0.0);
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
//  Print the visibility/overlapping probability
//  between two wavepackets of the circuit photon model.
//
//----------------------------------------
double qocircuit:: emitted_vis(int i,int j){
//  int i;      // Wave packet index 1
//  int j;      // Wave packet index 2

    return emitted->visibility(i,j);
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
photon_mdl::photon_mdl(mati i_pack_def, matd i_times, matd i_freq, char ckind){
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
//  Updates the photon model to include definitions of packets displaced
//  a quantity dt.
//
//----------------------------------------------------------------------------
void photon_mdl::update(double dt){
//  int  dt;            // Quantity of time the new wave-packets are displaced with respect the old ones.
//  Variables
    mati new_def;       // Updated packet definition
    matd newtimes;      // Updated times
    int  nP;            // Number of packets
    int  nT;            // Number of times
//  Auxiliary index
    int  i;             // Aux index


    // Reserve memory
    nP=pack_def.cols();
    nT=times.cols();
    new_def.resize(2,2*nP);
    newtimes.resize(2,2*nT);

    // Copy elements of the old model
    for(i=0;i<nP;i++){
        new_def(0,i)=pack_def(0,i);
        new_def(1,i)=pack_def(1,i);
    }

    for(i=0;i<nT;i++){
        newtimes(0,i)=times(0,i);
        newtimes(1,i)=times(1,i);
    }

    // Add new elements to the model
    for(i=0;i<nP;i++){
        new_def(0,nP+i)=pack_def(0,i)+nT;
        new_def(1,nP+i)=pack_def(1,i);
    }

    for(i=0;i<nT;i++){
        newtimes(0,nT+i)=times(0,i)+dt;
        if(times(1,i)>xcut) newtimes(1,nT+i)=times(1,i)+dt;
        else                newtimes(1,nT+i)=times(1,i);
    }

    // Assemble the model
    pack_def=new_def;
    times=newtimes;
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
    P.resize(4,pack_def.cols());

    // Create model
    for(i=0;i<pack_def.cols();i++){
            P(0,i)=times(0,pack_def(0,i));
            P(1,i)=freq(0,pack_def(1,i));
            P(2,i)=freq(1,pack_def(1,i));
            P(3,i)=times(1,pack_def(0,i));
    }

    //Returns the matrix model
    return P;
}


//----------------------------------------
//
//  Probability of two wave packets to overlap
//
//----------------------------------------
double photon_mdl:: visibility(int i,int j){
//  int    i         // Packet i number
//  int    j         // Packet j number
//  Variables
    cmplx  c;        // Coupling between two packets
//  Packet distribution variables
    double ti;       // Packet i time
    double tj;       // Packet j time
    double tri;      // Packet i phase time
    double trj;      // Packet j phase time
    double wi;       // Packet i frequency
    double wj;       // Packet j frequency
    double dwi;      // Packet i frequency width or characteristic time
    double dwj;      // Packet j frequency width or characteristic time
//  Model variables
    matd   P;        // Photon packet matrix (calculated from the definition)


    // Calculate the packet matrix from the photon model
    P=create_packet_mtx();

    // Obtain the parameters from the packet matrix
    ti=  P(0,i);
    tj=  P(0,j);
    wi=  P(1,i);
    wj=  P(1,j);
    dwi= P(2,i);
    dwj= P(2,j);
    tri= P(3,i);
    trj= P(3,j);

    // Calculate the coupling
    if(kind==0){
        c=gauss_coup(ti,wi,dwi,tj,wj,dwj,tri,trj);
    }else{
        c=exp_coup(ti,wi,dwi,tj,wj,dwj,tri,trj);
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


    for(i=0;i<times.cols();i++) cout << setw(2) << right << i << ": " << setw(10) << right << setprecision(5)  << times(0,i) << "   " << right  << setprecision(5) << times(1,i) << endl;
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

