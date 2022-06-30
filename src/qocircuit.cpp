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


    create_circuit(i_nch, 1, 1, 0, 0, false);
}


//----------------------------------------
//
//  Create circuit
//
//----------------------------------------
qocircuit::qocircuit(int i_nch, int i_nm, int i_ns){
//  int i_nch;  // Number of channels
//  int i_nm;   // Number of modes
//  int i_ns;   // Number of spatial wavefunctions/times


    create_circuit(i_nch, i_nm, i_ns, 0 ,0, false);
}


//----------------------------------------
//
//  Create circuit
//
//----------------------------------------
qocircuit::qocircuit(int i_nch, int i_nm, int i_ns, int clock){
//  int i_nch;  // Number of channels
//  int i_nm;   // Number of modes
//  int i_ns;   // Number of spatial wavefunctions/times


    create_circuit(i_nch, i_nm, i_ns, clock ,0, false);
}


//----------------------------------------
//
//  Create circuit
//
//----------------------------------------
qocircuit::qocircuit(int i_nch, int i_nm, int i_ns, int clock, int i_R, bool loss){
//  int i_nch;  // Number of channels
//  int i_nm;   // Number of modes
//  int i_ns;   // Number of spatial wavefunctions/times
//  int i_R;    // Number of iterations to calculate blinking and dark counts
//  bool loss;  // Are losses going to be calculated explicitly? True=Yes/False=No


    create_circuit(i_nch, i_nm, i_ns, clock, i_R, loss);
}


//----------------------------------------
//
//  Auxiliary private method to create a circuit
//
//----------------------------------------
void qocircuit::create_circuit(int i_nch, int i_nm, int i_ns, int clock, int i_R, bool loss){
//  int i_nch;  // Number of channels
//  int i_nm;   // Number of modes
//  int i_ns;   // Number of spatial wavefunctions/times
//  int i_R;    // Number of iterations to calculate blinking and dark counts
//  bool loss;  //Are losses going to be calculated explicitly? True=Yes/False=No
//  Auxiliary index
    int i;      // Aux index
    int j;      // Aux index
    int k;      // Aux index
    int l;      // Aux index


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
    ncond=0;
    timed=clock;
    det_def.resize(2,i_nch);
    det_par.resize(2,i_nch);

    // Initialize emitter/Packet definitions
    npack=0;
    pack_list.resize(4,ns);
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
    newcircuit= new qocircuit(nch,nm,ns,timed);
    newcircuit->losses=losses; // It has to be cloned that way. Otherwise the full constructor will double again the number of channels.

    //Update print variable
    newcircuit->flag_prnt=flag_prnt;

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
    newcircuit->ncond=ncond;
    newcircuit->det_def=det_def;
    newcircuit->det_par=det_par;

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
    int    ch1;               // Channel evaluated 1
    int    ch2;               // Channel evaluated 2
    int    km;                // Mode of the channels
    int    ks;                // Wavefunction/"Time" of the channels
    double theta;             // Angle theta in radians
    double phi;               // Angle phi in radians
    matc   U(nbmch,nbmch);    // Beamsplitter 2x2 matrix definitions
    veci   V(nbmch);          // Channels to which U refers
    matc   oelement;          // Extended nlevelxnlevel beamsplitter definition
//  Auxiliary index.
    int    i;                 // Aux index
    int    j;                 // Aux index
    int    k;                 // Aux index
    int    l;                 // Aux index.


    // Conversion to radians
    theta=d_theta*pi/180.0;
    phi  =d_phi*pi  /180.0;

    // 2x2 Matrix initialization
    U(0,0)= cos(theta);
    U(0,1)=-exp( jm*phi)*sin(theta);
    U(1,0)= exp(-jm*phi)*sin(theta);
    U(1,1)= cos(theta);
    V(0)  = i_ch1;
    V(1)  = i_ch2;

    //Conversion to nlevelxnlevel matrix
    oelement.resize(nlevel,nlevel);
    oelement=matc::Identity(nlevel,nlevel);
    for(k=0;k<nbmch;k++){
        ch1=V(k);
        for(l=0;l<nbmch;l++){
            ch2=V(l);
            for(km=0;km<nm;km++){
            for(ks=0;ks<ns;ks++){
                if((ch1<nch)&&(ch2<nch)){
                    i=i_idx[ch1][km][ks];
                    j=i_idx[ch2][km][ks];
                    oelement(i,j)=U(k,l);
                }
            }}
        }
    }

    // Update circuit with new element
    circmtx=oelement*circmtx;
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
    int   ch1;                // Channel evaluated 1
    int   ch2;                // Channel evaluated 2
    int   km;                 // Mode of the channels
    int   ks;                 // Wavefunction/"Time" of the channels
    cmplx rmt;                // r-t
    cmplx rpt;                // r+t
    cmplx rmt2;               // |r-t|^2
    cmplx rpt2;               // |r+t|^2
    matc  U(nbmch,nbmch);     // Dielectric 2x2 matrix definitions
    veci  V(nbmch);           // Channels to which U refers
    matc  oelement;           // Extended nlevelxnlevel dielectric definition
//  Auxiliary index.
    int   i;                  // Aux index
    int   j;                  // Aux index
    int   k;                  // Aux index
    int   l;                  // Aux index.


    // Check physicality of the input
    rmt=r-t;
    rpt=r+t;
    rmt2=conj(rmt)*rmt;
    rpt2=conj(rpt)*rpt;
    if(real(abs(rmt2))>1.0) cout << "Warning! Beamsplitter: t-r condition broken. The result aren't physical." << endl;
    if(real(abs(rpt2))>1.0) cout << "Warning! Beamsplitter: t+r condition broken. The result aren't physical." << endl;

    // Dielectric matrix
    U(0,0)= t;
    U(0,1)= r;
    U(1,0)= r;
    U(1,1)= t;

    // Channels
    V(0)  = i_ch1;
    V(1)  = i_ch2;

    //Conversion to nlevelxnlevel matrix
    oelement.resize(nlevel,nlevel);
    oelement=matc::Identity(nlevel,nlevel);
    for(k=0;k<nbmch;k++){
        ch1=V(k);
        for(l=0;l<nbmch;l++){
            ch2=V(l);
            for(km=0;km<nm;km++){
            for(ks=0;ks<ns;ks++){
                if((ch1<nch)&&(ch2<nch)){
                    i=i_idx[ch1][km][ks];
                    j=i_idx[ch2][km][ks];
                    oelement(i,j)=U(k,l);
                }
            }}
        }
    }

    // Update circuit with new element
    circmtx=oelement*circmtx;
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
    int  s;      // "Time"
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
            D(i,k)=sqrt(1.0-abs(conj(D(i,i))*D(i,i)));
            D(k,i)=sqrt(1.0-abs(conj(D(i,i))*D(i,i)));
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
    int   ch1;               // Channel evaluated 1
    int   ch2;               // Channel evaluated 2
    int   km;                // Mode of the channels
    int   ks;                // Wavefunction/"Time" of the channels
    matc  U(nbmch,nbmch);    // NSX 2x2 matrix definitions
    veci  V(nbmch);          // Channels to which U refers
    matc  oelement;          // Extended nlevelxnlevel NSX definition
//  Auxiliary index.
    int   i;                 // Aux index
    int   j;                 // Aux index
    int   k;                 // Aux index
    int   l;                 // Aux index.


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
    V(0)  = i_ch1;
    V(1)  = i_ch2;
    V(2)  = i_ch3;

    //Conversion to nlevelxnlevel matrix
    oelement.resize(nlevel,nlevel);
    oelement=matc::Identity(nlevel,nlevel);
    for(k=0;k<nbmch;k++){
        ch1=V(k);
        for(l=0;l<nbmch;l++){
            ch2=V(l);
            for(km=0;km<nm;km++){
            for(ks=0;ks<ns;ks++){
                if((ch1<nch)&&(ch2<nch)){
                    i=i_idx[ch1][km][ks];
                    j=i_idx[ch2][km][ks];
                    oelement(i,j)=U(k,l);
                }
            }}
        }
    }

    // Update circuit with new element
    circmtx=oelement*circmtx;
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
    int    km1;               // Mode of the channel i_ch
    int    km2;               // Mode of the channel i_ch
    int    ks;                // Wavefunction/"Time" of the channels
    double theta;             // Angle theta in radians
    double phi;               // Angle phi in radians
    matc   U(nbmch,nbmch);    // Rotator 2x2 matrix definitions
    veci   P(nbmch);          // Channels to which U refers
    matc   oelement;          // Extended nlevelxnlevel rotator definition
//  Auxiliary index.
    int    i;                 // Aux index.
    int    j;                 // Aux index.
    int    k;                 // Aux index
    int    l;                 // Aux index.


    // Conversion to radians
    theta=d_theta*pi/180.0;
    phi  =d_phi*pi  /180.0;

    // 2x2 Matrix initialization
    U(0,0)= cos(theta);
    U(0,1)=-exp( jm*phi)*sin(theta);
    U(1,0)= exp(-jm*phi)*sin(theta);
    U(1,1)= cos(theta);
    P(0)  = H;
    P(1)  = V;

    //Conversion to nlevelxnlevel matrix
    oelement.resize(nlevel,nlevel);
    oelement=matc::Identity(nlevel,nlevel);
    for(k=0;k<nbmch;k++){
        km1=P(k);
        for(l=0;l<nbmch;l++){
            km2=P(l);
            for(ks=0;ks<ns;ks++){
                i=i_idx[i_ch][km1][ks];
                j=i_idx[i_ch][km2][ks];
                oelement(i,j)=U(k,l);

            }
        }
    }

    // Update circuit with new element
    circmtx=oelement*circmtx;
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
    int   ch1;                // Channel evaluated 1
    int   ch2;                // Channel evaluated 2
    int   km1;                // Mode of the channel i_ch1
    int   km2;                // Mode of the channel i_ch2
    int   ks;                 // Wavefunction/"Time" of the channels
    matc  U(nbmch,nbmch);     // Polarized beamsplitter 4x4 matrix definitions
    veci  CH(nbmch);          // Channels to which U refers
    veci  P(nbmch);           // Channels to which U refers
    matc  oelement;           // Extended nlevelxnlevel polarized beamsplitter definition
//  Auxiliary index.
    int   i;                  // Aux index.
    int   j;                  // Aux index.
    int   k;                  // Aux index
    int   l;                  // Aux index.


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
    CH(0)  =  i_ch1;
    CH(1)  =  i_ch1;
    CH(2)  =  i_ch2;
    CH(3)  =  i_ch2;
    P(0)   =  pol;
    P(1)   = (pol+1)%2;
    P(2)   =  pol;
    P(3)   = (pol+1)%2;

    //Conversion to nlevelxnlevel matrix
    oelement.resize(nlevel,nlevel);
    oelement=matc::Identity(nlevel,nlevel);
    for(k=0;k<nbmch;k++){
        ch1=CH(k);
        km1=P(k);
        for(l=0;l<nbmch;l++){
            ch2=CH(l);
            km2=P(l);
            for(ks=0;ks<ns;ks++){
                i=i_idx[ch1][km1][ks];
                j=i_idx[ch2][km2][ks];
                oelement(i,j)=U(k,l);

            }
        }
    }

    // Update circuit with new element
    circmtx=oelement*circmtx;
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
    int    km1;               // Mode of the channel i_ch
    int    km2;               // Mode of the channel i_ch
    int    ks;                // Wavefunction/"Time" of the channels
    double alpha;             // Angle alpha in radians
    double gamma;             // Angle gamma in radians
    matc   U(nbmch,nbmch);    // Waveplate 2x2 matrix definitions
    veci   P(nbmch);          // Channels to which U refers
    matc   oelement;          // Extended nlevelxnlevel waveplate definition
//  Auxiliary index.
    int    i;                 // Channel 1 matrix index
    int    j;                 // Channel 2 matrix index
    int    k;                 // Aux index
    int    l;                 // Aux index.


    // Conversion to radians
    alpha=d_alpha*pi/180.0;
    gamma=d_gamma*pi/180.0;

    // 2x2 Matrix initialization
    U(0,0)= jm*sin(gamma)*cos(2.0*alpha)+cos(gamma);
    U(0,1)= jm*sin(gamma)*sin(2.0*alpha);
    U(1,0)= jm*sin(gamma)*sin(2.0*alpha);
    U(1,1)=-jm*sin(gamma)*cos(2.0*alpha)+cos(gamma);
    P(0)  = H;
    P(1)  = V;

    //Conversion to nlevelxnlevel matrix
    oelement.resize(nlevel,nlevel);
    oelement=matc::Identity(nlevel,nlevel);
    for(k=0;k<nbmch;k++){
        km1=P(k);
        for(l=0;l<nbmch;l++){
            km2=P(l);
            for(ks=0;ks<ns;ks++){
                i=i_idx[i_ch][km1][ks];
                j=i_idx[i_ch][km2][ks];
                oelement(i,j)=U(k,l);

            }
        }
    }

    // Update circuit with new element
    circmtx=oelement*circmtx;
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
    int  km1;         // Mode of the channel i_ch
    int  km2;         // Mode of the channel i_ch
    int  ks;          // Wavefunction/"Time" of the channels
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
            P(i)   = 0;
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
            for(ks=0;ks<ns;ks++){
                i=i_idx[ch1][km1][ks];
                j=i_idx[ch2][km2][ks];
                oelement(i,j)=U(k,l);

            }
        }
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
//  Variables
    int   km;                 // Mode of the channels
    int   ks;                 // Wavefunction/"Time" of the channels
    matc  oelement;           // nlevelxnlevel phase shifter definition
//  Auxiliary index.
    int   i;                  // Aux index.


    // Conversion to nlevelxnlevel matrix
    oelement.resize(nlevel,nlevel);
    oelement=matc::Identity(nlevel,nlevel);
    for(km=0;km<nm;km++){
    for(ks=0;ks<ns;ks++){
        i=i_idx[i_ch][km][ks];
        oelement(i,i)=t;
    }}

    // Update circuit with new element
    circmtx=oelement*circmtx;
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
void qocircuit:: emitter(photon_mdl *mdl){
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
    if(nP>ns) cout << "Emitter(photon model): Not enough ns degrees of freedom! Needed: "<< nP << endl;


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
}


//--------------------------------------------------
//
// Emitter configured using a packet description in
// in a matrix.
//
//--------------------------------------------------
void qocircuit::emitter(int npack, matd packets, char ckind, int rand){
//  int  npack;     // Number of packets defined.
//  matd packets;   // Matrix that contains the definition of the packets.
//  matd packets;   // Kind of packets  'G': Gaussian/'E': Exponential
//  Variables
    vecd longtimes; // List with the times defined
    matd longfreq;  // List with the frequencies defined
    vecd times;     // List with the times to create the photon model. Same size as entries
    matd freq;      // List with the frequencies to create the photon model. Same number of columns as entries
    mati defs;      // List with packet definitions. (Equal to pack unless timed=1)
    mati pack;      // Packet definition matrix to create the photon model.
    int  ntimes;    // Number of times in longtimes
    int  nfreq;     // Number of frequencies in long freq
//  Aux index
    int i;          // Aux index
    int j;          // Aux index
    int k;          // Aux index


    // Perform the format conversion
    ntimes=0;
    nfreq=0;
    longtimes.resize(npack);
    longfreq.resize(2,npack);
    defs.resize(3,npack);
    for(i=0;i<npack;i++){
        j=0;
        while((j<ntimes)&&(abs(packets(1,i)-longtimes(j))>xcut)) j++;
        if (j==ntimes){
            longtimes(ntimes)=packets(1,i);
            ntimes=ntimes+1;
        }

        k=0;
        while((k<nfreq)&&((abs(packets(2,i)-longfreq(0,k))>xcut)
                       ||(abs(packets(3,i)-longfreq(1,k))>xcut))) k++;
        if (k==nfreq){
            longfreq(0,nfreq)=packets(2,i);
            longfreq(1,nfreq)=packets(3,i);
            nfreq=nfreq+1;
        }

        defs(0,i)=(int) packets(0,i);
        defs(1,i)=j;
        defs(2,i)=k;
    }

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
    }else{
        pack=defs;
    }

//  Adjust sizes
    times.resize(ntimes);
    freq.resize(2,nfreq);
    for(i=0;i<ntimes;i++) times(i)=longtimes(i);
    for(i=0;i<nfreq;i++){
        freq(0,i) = longfreq(0,i);
        freq(1,i) = longfreq(1,i);
    }

    // Create the photon model
    delete emitted;
    emitted= new photon_mdl(pack, times, freq, ckind, rand);

    // Set up the emitter using the photon model
    emitter(emitted);
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


    // Adds a new entry to the matrix of definitions of packets
    pack_list(0,npack)=(double) n;
    pack_list(1,npack)=         t;
    pack_list(2,npack)=         f;
    pack_list(3,npack)=         w;
    npack=npack+1;
}


//--------------------------------------------------
//
// Define an emitter
//
//--------------------------------------------------
void qocircuit:: emitter(char ckind, int rand){
//  char ckind;          // Kind of wave packets 'G': Gaussian / 'E': Exponential
//  int  rand;           //  How we consider the phases. 0=Not consider/1=Synchronized with the emission/2=Random (according to distribution)


    emitter(npack,pack_list,ckind,rand);
}


//--------------------------------------------------
//
// Adds a delay
//
//--------------------------------------------------
void qocircuit:: delay(int ch, double dt){
//  int    ch;   // Channel to be delayed
//  double dt;   // Amount of time to be delayed
//  Variables
    photon_mdl *new_mdl;  // New photon that consider the new packet definitions created by the delay.


    // Apply the delay and create the new photon model
    new_mdl=delay(ch,dt, emitted);

    // Update the photon model
    delete emitted;
    emitted=new_mdl;
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
    if(nP>ns) cout << "Delay(photon model): Not enough ns degrees of freedom! Needed: "<< nP << endl;

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
void qocircuit:: detector(int i_ch, int cond, double eff, double blnk, double gamma){
//  int    i_ch;        // Chanel where detection is performed
//  int    cond;        // Post selection condition. If cond<0 there is none.
//  double eff;         // Efficiency of the detector
//  double blnk;        // Fraction of time blinking
//  double ganna;       // Dark counts rate
//  Variables
    int nclose;         // Total number of channels


    // Adds a new detector definition entry to the corresponding matrices
    // If cond>=0 Adds a new entry to the matrix of post-selection condition.
    ndetc=ndetc+1;
    if(ndetc>nch)              cout << "Detector warning!: More detectors than channels are being declared." << endl;
    if((i_ch>=nch)||(i_ch<0))  cout << "Detector warning!: This channel does not exist." << endl;
    if(cond>=0){
        det_def(0,ncond)=i_ch;
        det_def(1,ncond)=cond;
        ncond=ncond+1;
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

    // If this is the last detector compute losses if needed
    if((losses==1)&&(ndetc==nclose)) compute_losses();

    // If this is the last detector add the emitter matrix
    // at the beginning. (Otherwise computing the losses breaks the calculation)
    if((emiss==1)&&(ndetc==nclose)) circmtx=circmtx*init_dmat;
}


//--------------------------------------------------
//
//  Adds a virtual circuit element to define a detector
//
//--------------------------------------------------
void qocircuit:: detector(int i_ch){
//  int    i_ch;        // Chanel where detection is performed


    detector(i_ch,-1,1.0,0.0,0.0);
}

//--------------------------------------------------
//
//  Adds a virtual circuit element to define a detector
//
//--------------------------------------------------
void qocircuit:: detector(int i_ch, int cond){
//  int    i_ch;        // Chanel where detection is performed
//  int    cond;        // Post selection condition. If cond<0 there is none.


    detector(i_ch,cond,1.0,0.0,0.0);
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
void qocircuit:: prnt(){
//  Variables
    int ch1;       // Channel 1
    int ch2;       // Channel 2
    int m1;        // Mode 1
    int m2;        // Mode 2
    int s1;        // Spatial wavefunction/"Time" 1
    int s2;        // Spatial wavefunction/"Time" 2
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
        if(flag_prnt==0){
            if(nm>1) cout << ",-log(1-urand()) "<< m1;
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
                if(flag_prnt==0){
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
//  Sets the flag that controls the print style.
//
//----------------------------------------
void qocircuit::set_prnt_flag(int flag){
 // int flag;   // Value of the flag to be adopted


    flag_prnt=flag;
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
photon_mdl::photon_mdl(mati i_pack_def, vecd i_times, matd i_freq, char ckind, int i_rand){
//  mati pack_def       // Packet definition
//  matd times          // Times definition for each index /enumeration
//  matd freq           // Frequency definition for each index/enumeration
//  char ckind          // Kind of wave packets 'G': Gaussian / 'E': Exponential
//  int  rand           //  How we consider the phases. 0=Not consider/1=Synchronized with the emission/2=Random (according to distribution)
//  Variables
    matd times_plus;    // Times plus phase times referred in the packets definitions
//  Auxiliary index
    int  i;             // Aux index


    // Extend the times definition to include phase times (in case they are needed)
    times_plus.resize(TD,i_times.size());
    for(i=0;i<i_times.size();i++){
        times_plus(0,i)=i_times(i);
        times_plus(1,i)=-log(1-urand());
    }

    // Configure the distribution of the packets
    pack_def=create_packet_idx(i_pack_def);
    times=times_plus;
    freq=i_freq;
    kind=0;
    if(ckind=='E') kind=1;
    rand=i_rand;
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
    new_mdl->rand=rand;

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
    new_def.resize(PT,2*nP);
    newtimes.resize(TD,2*nT);

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
        newtimes(1,nT+i)=times(1,i);
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
    P.resize(PD,pack_def.cols());

    // Create model
    for(i=0;i<pack_def.cols();i++){
            P(0,i)=times(0,pack_def(0,i));
            P(1,i)=freq(0,pack_def(1,i));
            P(2,i)=freq(1,pack_def(1,i));
            switch (rand){
                case 0: // Turn off
                    P(3,i)=0;
                    break;
                case 1: // Fix
                    P(3,i)=times(0,pack_def(0,i));
                    break;
                case 2: // Random
                    P(3,i)=times(0,pack_def(0,i))+times(1,pack_def(0,i));
                    break;
                default:
                    cout << "create_packet_mtx: Something went wrong." << endl;
                    break;
            }
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
    cout << "Random phase configuration: " << rand <<endl;
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
    pack_idx.resize(PT,pack_def.cols());
    //Transform representation
    for(i=0;i<pack_def.cols();i++){
        pack_idx(0,pack_def(0,i))=pack_def(1,i);
        pack_idx(1,pack_def(0,i))=pack_def(2,i);
    }
    // Return matrix
    return pack_idx;
}

