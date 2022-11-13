//======================================================================================================
// File qodev.cpp
//
// METACIRCUIT LIBRARY
//
// Copyright © 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved.
// The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the
// root directory of this source tree. Use of the source code in this file is only permitted under the
// terms of the licence in LICENCE.TXT.
//======================================================================================================
#include "qodev.h"


//-----------------------------------------------------
//
//  Creates a device
//
//-----------------------------------------------------
qodev::qodev(int i_nph, int i_nch){
//  int  i_nph;     // Maximum number of photons to be simulated
//  int  i_nch;     // Number of channels


    cfg_soqcs(i_nph);
    circ=new qocircuit(i_nch,1,1, 0, 0, false,'G');
    create_qodev(circ->num_levels(),DEFSTATEDIM,1);
}


//-----------------------------------------------------
//
//  Creates a polarized device
//
//-----------------------------------------------------
qodev::qodev(int i_nph, int i_nch, int i_nm){
//  int  i_nph;     // Maximum number of photons to be simulated
//  int  i_nch;     // Number of channels
//  int  i_nm;      // Number of modes


    cfg_soqcs(i_nph);
    circ=new qocircuit(i_nch,i_nm,1, 0, 0, false,'G');
    create_qodev(circ->num_levels(),DEFSTATEDIM,1);
}


//-----------------------------------------------------
//
//  Creates a physical device
//
//-----------------------------------------------------
qodev::qodev(int i_nph,int i_nch, int i_nm, int i_ns, int clock, char ckind){
//  int  i_nph;     // Maximum number of photons to be simulated
//  int  i_nch;     // Number of channels
//  int  i_nm;      // Number of modes
//  int  i_ns;      // Number of packets
//  int  clock;     // Kind of detectors
//  char ckind;     // Kind of packets  'G': Gaussian/'E': Exponential


    cfg_soqcs(i_nph);
    circ=new qocircuit(i_nch,i_nm,i_ns, clock, 0, false,ckind);
    create_qodev(circ->num_levels(),DEFSTATEDIM,i_ns);
}



//-----------------------------------------------------
//
//  Creates a physical device with physical detectors
//
//-----------------------------------------------------
qodev::qodev(int i_nph, int i_nch, int i_nm, int i_ns, int clock, int i_R, bool loss, char ckind, int i_maxket){
//  int  i_nph;     // Maximum number of photons to be simulated
//  int  i_nch;     // Number of channels
//  int  i_nch;     // Number of channels
//  int  i_nm;      // Number of modes
//  int  i_ns;      // Number of packets
//  int  clock;     // Kind of detectors
//  int  i_R;       // Number of iterations to calculate blinking and dark counts
//  bool loss;      // Are losses going to be calculated explicitly? True=Yes/False=No
//  char ckind;     // Kind of packets  'G': Gaussian/'E': Exponential
//  int  i_maxket;  // Maximum number of bunches in the list


    cfg_soqcs(i_nph);
    circ=new qocircuit(i_nch,i_nm,i_ns, clock, i_R, loss,ckind);
    create_qodev(circ->num_levels(),i_maxket,i_ns);
}


//----------------------------------------
//
//  Auxiliary private method to create a device
//
//----------------------------------------
void qodev::create_qodev(int i_level, int i_maxket, int i_np){
//  int  i_level;     // Number of levels to describe the state
//  int  i_maxket;    // Maximum number of bunches in the list
//  int  i_np;        // Number of packets
//  Variables
    int *occ;         // Level occupation


    npack=0;
    pack_list.resize(5,i_np);
    inpt=new state(i_level,i_maxket);
    occ=new int[i_level]();
    inpt->add_term(1.0,occ);

    delete occ;
}


//----------------------------------------
//
//  Destroys a device
//
//----------------------------------------
qodev::~qodev(){

    delete circ;
    delete inpt;

}


//----------------------------------------
//
//  Clears a device
//
//----------------------------------------
void qodev::reset(){
//  Variables
    int       *occ;       // Level occupation

    // Reset packet table
    npack=0;

    // Create a new empty initial state
    inpt->clear();
    occ=new int[inpt->nlevel]();
    inpt->add_term(1.0,occ);
    delete occ;

    // Create a new empty circuit
    circ->reset();
}


//----------------------------------------
//
//  Concatenates two devices
//
//----------------------------------------
void qodev::concatenate(qodev *dev){
//  qodev *dev;      // Meta Circuit to be appended


    if((circ->emiss==0)&&(dev->circ->remdec()==0)) send2circuit();
    circ->concatenate(dev->circ);
}


//----------------------------------------
//
//  Add photons to a device
//
//----------------------------------------
int qodev::add_photons(int N, int ch){
//  int N;      // Number of photons to be added
//  int ch;     // Chanel where the photons are added


    return add_photons(N,ch, 0, 0.0, 0.0, 0.0);
}


//----------------------------------------
//
//  Add photons to a device with physical
//  parameters about the photon shape
//
//----------------------------------------
int qodev::add_photons(int N, int ch, int P, double t, double f, double w){
//  int    N;        // Number of photons to be added
//  int   ch;        // Chanel where the photons are added
//  int    P;        // Polarization of the photons
//  double t;        // Photon emission time
//  double f;        // Photon emission frequency
//  double w;        // Photon emission  width or decay time depending on the packet shape model
//  Variables
    int    T;        // Packet number
    int    add;      // Add new packer 0=No/1=Yes. This way we check repetitions.
    int   *occ;      // Level occupation
    state *aux;      // Auxiliary state
    hterm  in_term;  // Human readable input term
    state *newinpt;  // New bunch of photons with the new photon information updated with respect the old one
// Auxiliary index
    int    i;        // Aux index
    int    j;        // Aux index


    // Check correctness
    if(circ->emiss==1){
        cout << "add_photons error #1: Photons already emitted. More photons can not be added at this stage" << endl;
        cout << "No new photons can be defined after a delay or dispersion elements " << endl;
        return -1;
    }

    // Check if the packet described for this entry already exists
    add=1;
    i=0;
    while((i<npack)&&(add==1)){
        if((pack_list(1,i)==t)&&(pack_list(2,i)==f)&&(pack_list(3,i)==w)) add=0;
        i++;
    }

    // Update the wave-packet information
    if(add==1){
        if(npack>=circ->ns) {
            cout << "add_photons error #2: Not enough ns degrees of freedom! Needed: "<< npack+1 << endl;
            return -1;
        }
        // If not exists add it
        T=npack;
        pack_list(0,npack)=(double) npack;
        pack_list(1,npack)= t;
        pack_list(2,npack)= f;
        pack_list(3,npack)= w;
        pack_list(4,npack)= t;
        npack=npack+1;
    }else{
        // If exist just take the index
        T=i-1;
    }


    // Update the input state information
    // Compute the new state
    in_term.resize(4,1);
    in_term << ch,
               P,
               T,
               N;
    aux=new state(inpt->nlevel,1);
    aux->add_term(1.0,in_term,circ);

    newinpt=new state(inpt->nlevel,inpt->maxket);
    for(j=0;j<inpt->nket;j++){
        occ=new int[inpt->nlevel]();
        for(i=0;i<inpt->nlevel;i++){
            occ[i]=inpt->ket[j][i]+aux->ket[0][i];
        }

        newinpt->add_term(1.0,occ);
        delete occ;
    }
    //Update internally the bunch
    delete inpt;
    inpt=newinpt;

    return 0;
}


//----------------------------------------
//
//  Add photons to a device as generated by a QD
//
//----------------------------------------
int qodev::add_QD(int ch1, int ch2,double i_t1, double i_f1, double i_w1, double i_t2, double i_f2, double i_w2, int PW, double S, double k, double tss, double thv){
//  int    ch1;      // Chanel where the photons are added
//  int    ch2;      // Chanel where the photons are added
//  double t1;       // Photon 1 emission time
//  double f1;       // Photon 1 emission frequency
//  double w1;       // Photon 1 emission width or decay time depending on the packet shape model
//  double t2;       // Photon 2 emission time
//  double f2;       // Photon 2 emission frequency
//  double w2;       // Photon 2 emission width or decay time depending on the packet shape model
//  int    PW;       // Do we consider frequency shift in emission because FSS?. 0='No'/1='Yes'
//  double S;        // FSS
//  double k;        // Signal to noise ration
//  double tss;      // Spin scattering characteristic time
//  double thv;      // Cross dephasing characteristic time


    return this->add_QD( ch1, ch2, i_t1, i_t1, i_f1, i_w1, i_t2, i_t2, i_f2, i_w2, PW, S, k, tss, thv, 0.0, 2);
}


//----------------------------------------
//
//  Add photons to a device as generated by a QD in a cascade
//
//----------------------------------------
int qodev::add_QD(int ch1, int ch2,double i_t1, double i_f1, double i_w1, double i_t2, double i_f2, double i_w2, int PW, double S, double k, double tss, double thv, double T2){
//  int    ch1;      // Chanel where the photons are added
//  int    ch2;      // Chanel where the photons are added
//  double t1;       // Photon 1 emission time
//  double f1;       // Photon 1 emission frequency
//  double w1;       // Photon 1 emission width or decay time depending on the packet shape model
//  double t2;       // Photon 2 emission time
//  double f2;       // Photon 2 emission frequency
//  double w2;       // Photon 2 emission width or decay time depending on the packet shape model
//  int    PW;       // Do we consider frequency shift in emission because FSS?. 0='No'/1='Yes'
//  double S;        // FSS
//  double k;        // Signal to noise ration
//  double tss;      // Spin scattering characteristic time
//  double thv;      // Cross dephasing characteristic time
//  double T2;       // Pure  dephasing characteristic time


    return this->add_QD( ch1, ch2, i_t1, i_t1, i_f1, i_w1, i_t2, i_t2, i_f2, i_w2, PW, S, k, tss, thv, T2, 1);
}


//----------------------------------------
//
//  Add photons to a device as generated by a QD
//
//----------------------------------------
int qodev::add_QD(int ch1, int ch2,double i_t1, double i_tp1, double i_f1, double i_w1, double i_t2, double i_tp2, double i_f2, double i_w2, int PW, double S, double k, double tss, double thv, double T2, int rand){
//  int    ch1;      // Chanel where the photons are added
//  int    ch2;      // Chanel where the photons are added
//  double t1;       // Photon 1 emission time
//  double tp1;      // Photon 1 emission phase time
//  double f1;       // Photon 1 emission frequency
//  double w1;       // Photon 1 emission width or decay time depending on the packet shape model
//  double t2;       // Photon 2 emission time
//  double tp2;      // Photon 2 emission phase time
//  double f2;       // Photon 2 emission frequency
//  double w2;       // Photon 2 emission width or decay time depending on the packet shape model
//  int    PW;       // Do we consider frequency shift in emission because FSS?. 0='No'/1='Yes'
//  double S;        // FSS
//  double k;        // Signal to noise ration
//  double tss;      // Spin scattering characteristic time
//  double thv;      // Cross dephasing characteristic time
//  double T2;       // Pure  dephasing characteristic time
//  int    rand      // Phase computation mode
//  Variables
    int    add;      // Add new packer 0=No/1=Yes. This way we check repetitions.
    double dt1;      // Extra emission time
    double dt2;      // Extra phase time
    veci   T;        // Packet number
    mati   ch;       // QD configuration matrix
    hterm  in_term;  // Human readable input term
    state *qdstate;  // Auxiliary state
// Auxiliary index
    int    i;        // Aux index


    // Check correctness
    if(circ->emiss==1){
        cout << "add_QD error #0: Photons already emitted. More photons can not be added at this stage" << endl;
        cout << "No new photons can be defined after a delay or dispersion elements " << endl;
        return -1;
    }

    // Reserve memory
    T.resize(4);

    // FIRST PAIR
    // FIRST PACKET
    add=1;
    i=0;
    // Check if the packet described for this entry already exists
    while((i<npack)&&(add==1)){
    if((pack_list(1,i)==i_t1)&&(pack_list(2,i)==i_f1)&&(pack_list(3,i)==i_w1)) add=0;
        i++;
    }

    // Update the wave-packet information
    if(add==1){
        if(npack>=circ->ns) {
            cout << "add_QD error #1: Not enough ns degrees of freedom! Needed: "<< npack+1 << endl;
            return -1;
        }
        T(0)=npack;
        pack_list(0,npack)=(double) npack;
        pack_list(1,npack)= i_t1;
        pack_list(2,npack)= i_f1;
        pack_list(3,npack)= i_w1;
        pack_list(4,npack)= i_tp1;
        npack=npack+1;
    }else{
        // If exist just take the index
        T(0)=i-1;
    }

    // SECOND PACKET
    if(PW==1){
        add=1;
        i=0;
        // Check if the packet described for this entry already exists
        while((i<npack)&&(add==1)){
        if((pack_list(1,i)==i_t1)&&(pack_list(2,i)==i_f1+S)&&(pack_list(3,i)==i_w1)) add=0;
            i++;
        }

        // Update the wave-packet information
        if(add==1){
            if(npack>=circ->ns) {
                cout << "add_QD error #1B: Not enough ns degrees of freedom! Needed: "<< npack+1 << endl;
                return -1;
            }
            T(1)=npack;
            pack_list(0,npack)=(double) npack;
            pack_list(1,npack)= i_t1;
            pack_list(2,npack)= i_f1+S;
            pack_list(3,npack)= i_w1;
            pack_list(4,npack)= i_tp1;
            npack=npack+1;
        }else{
            // If exist just take the index
            T(1)=i-1;
        }
    }else{
        T(1)=T(0);
    }


    // SECOND PAIR
    // FIRST PACKET
    // Check if the packet described for this entry already exists
    add=1;
    i=0;
    // Check if the packet described for this entry already exists
    while((i<npack)&&(add==1)){
    if((pack_list(1,i)==i_t2)&&(pack_list(2,i)==i_f2)&&(pack_list(3,i)==i_w2)) add=0;
        i++;
    }


    if(circ->ckind=='G'){ // Gussian
        dt1=erfi(2*urand()-1)/i_w1;
        dt2=(2*i_w2)*erfi(2*urand()-1)/(T2*i_w1);
    }else{               // Exponential
        dt1=i_w1*expi(urand());
        dt2=(T2*i_w1)*expi(urand())/(2.0*i_w2);
    }

    // Update the wave-packet information
    if(add==1){
        if(npack>=circ->ns) {
            cout << "add_QD error #2: Not enough ns degrees of freedom! Needed: "<< npack+1 << endl;
            return -1;
        }
        T(2)=npack;
        pack_list(0,npack)=(double) npack;
        pack_list(2,npack)=         i_f2;
        pack_list(3,npack)=         i_w2;
        switch(rand) {
            case 0: // Fully Automatic #1
                pack_list(1,npack)= i_t2+dt1;
                pack_list(4,npack)= i_tp2+dt1;
                break;
            case 1: // Automatic #2
                pack_list(1,npack)= i_t2+dt1;
                pack_list(4,npack)= i_tp2+dt2;
                break;
            case 2: // Semi manual #4
                pack_list(1,npack)= i_t2;
                pack_list(4,npack)= i_tp2+dt2;
                break;
            default:
                cout << "add_QD error #3: rand goes form 0 to 3. Invalid configuration." << endl;
                return -1;
                break;
        }
        npack=npack+1;

    }else{
        // If exist just take the index
        T(2)=i-1;
    }

    // SECOND PACKET
    if(PW==1){
        add=1;
        i=0;
        // Check if the packet described for this entry already exists
        while((i<npack)&&(add==1)){
        if((pack_list(1,i)==i_t2)&&(pack_list(2,i)==i_f2-S)&&(pack_list(3,i)==i_w2)) add=0;
            i++;
        }

        // Update the wave-packet information
        if(add==1){
            if(npack>=circ->ns) {
                cout << "add_QD error #2B: Not enough ns degrees of freedom! Needed: "<< npack+1 << endl;
                return -1;
            }
            T(3)=npack;
            pack_list(0,npack)=(double) npack;
            pack_list(2,npack)=         i_f2-S;
            pack_list(3,npack)=         i_w2;
            switch(rand) {
                case 0: // Fully automatic #1
                    pack_list(1,npack)= i_t2+dt1;
                    pack_list(4,npack)= i_tp2+dt1;
                    break;
                case 1: // Automatic #2
                    pack_list(1,npack)= i_t2+dt1;
                    pack_list(4,npack)= i_tp2+dt2;
                    break;
                case 2: // Semi manual #4
                    pack_list(1,npack)= i_t2;
                    pack_list(4,npack)= i_tp2+dt2;
                    break;
                default:
                    cout << "add_QD error #4: rand goes form 0 to 3. Invalid configuration." << endl;
                    return -1;
                    break;
            }
            npack=npack+1;
        }else{
            // If exist just take the index
            T(3)=i-1;
        }
    }else{
        T(3)=T(2);
    }


    // Configure table according to chosen order
    ch.resize(3,2);
    ch  << ch1,  ch2,
          T(0), T(2),
          T(1), T(3);

    // Update state
    qdstate=new state(inpt->nlevel,inpt->maxket);
    qdstate->QD(ch,k,S,i_w2,tss,thv,circ);
    inpt->dproduct(qdstate);
    delete qdstate;

    return 0;
}


//----------------------------------------
//
//  Add photons to a device as generated by
//  an unsynchronized Bell state. (Path encoding)
//
//----------------------------------------
int qodev::add_Bell(int ch1, int ch2, char kind, double phi, double t1, double f1, double w1, double t2, double f2, double w2){
//  int    ch1;      // Chanel where the photons are added
//  int    ch2;      // Chanel where the photons are added
//  char   kind;     // Kind of bell state +=Phi+/-=Phi-/p=Psi+/m=Psi-
//  double phi;      // Non-ideal emission phase
//  double t1;       // Photon 1 emission time
//  double f1;       // Photon 1 emission frequency
//  double w1;       // Photon 1 emission width or decay time depending on the packet shape model
//  double t2;       // Photon 2 emission time
//  double f2;       // Photon 2 emission frequency
//  double w2;       // Photon 2 emission width or decay time depending on the packet shape model
//  Variables
    int    add;      // Add new packer 0=No/1=Yes. This way we check repetitions.
    veci   T;        // Packet number
    mati   ch;       // QD configuration matrix
    hterm  in_term;  // Human readable input term
    state *qdstate;  // Auxiliary state
// Auxiliary index
    int    i;        // Aux index


    // Check correctness
    if(circ->emiss==1){
        cout << "add_Bell error #0: Photons already emitted. More photons can not be added at this stage" << endl;
        cout << "No new photons can be defined after a delay or dispersion elements " << endl;
        return -1;
    }

    // Reserve memory
    T.resize(4);

    // FIRST PACKET
    add=1;
    i=0;
    // Check if the packet described for this entry already exists
    while((i<npack)&&(add==1)){
    if((pack_list(1,i)==t1)&&(pack_list(2,i)==f1)&&(pack_list(3,i)==w1)) add=0;
        i++;
    }

    // Update the wave-packet information
    if(add==1){
        if(npack>=circ->ns) {
            cout << "add_Bell error #1: Not enough ns degrees of freedom! Needed: "<< npack+1 << endl;
            return -1;
        }
        T(0)=npack;
        pack_list(0,npack)=(double) npack;
        pack_list(1,npack)= t1;
        pack_list(2,npack)= f1;
        pack_list(3,npack)= w1;
        pack_list(4,npack)= t1;
        npack=npack+1;
    }else{
        // If exist just take the index
        T(0)=i-1;
    }

    // SECOND PACKET
    // Check if the packet described for this entry already exists
    add=1;
    i=0;
    // Check if the packet described for this entry already exists
    while((i<npack)&&(add==1)){
    if((pack_list(1,i)==t2)&&(pack_list(2,i)==f2)&&(pack_list(3,i)==w2)) add=0;
        i++;
    }

    // Update the wave-packet information
    if(add==1){
        if(npack>=circ->ns) {
            cout << "add_Bell error #2: Not enough ns degrees of freedom! Needed: "<< npack+1 << endl;
            return -1;
        }
        T(1)=npack;
        pack_list(0,npack)=(double) npack;
        pack_list(1,npack)= t2;
        pack_list(2,npack)= f2;
        pack_list(3,npack)= w2;
        pack_list(4,npack)= t2;
        npack=npack+1;
    }else{
        // If exist just take the index
        T(1)=i-1;
    }

    // Configure table according to chosen order
    ch.resize(2,2);
    ch  << ch1,  ch2,
          T(0), T(1);

    // Update state
    qdstate=new state(inpt->nlevel,inpt->maxket);
    qdstate->Bell(ch,kind,phi,circ);
    inpt->dproduct(qdstate);
    delete qdstate;

    return 0;
}


//----------------------------------------
//
//  Add photons to a device as generated by
//  an unsynchronized Bell state. (Polarization encoding)
//
//----------------------------------------
int qodev::add_BellP(int ch1, int ch2, char kind, double phi, double t1, double f1, double w1, double t2, double f2, double w2){
//  int    ch1;      // Chanel where the photons are added
//  int    ch2;      // Chanel where the photons are added
//  char   kind;     // Kind of bell state +=Phi+/-=Phi-/p=Psi+/m=Psi-
//  double phi;      // Non-ideal emission phase
//  double t1;       // Photon 1 emission time
//  double f1;       // Photon 1 emission frequency
//  double w1;       // Photon 1 emission width or decay time depending on the packet shape model
//  double t2;       // Photon 2 emission time
//  double f2;       // Photon 2 emission frequency
//  double w2;       // Photon 2 emission width or decay time depending on the packet shape model
//  Variables
    int    add;      // Add new packer 0=No/1=Yes. This way we check repetitions.
    veci   T;        // Packet number
    mati   ch;       // QD configuration matrix
    hterm  in_term;  // Human readable input term
    state *qdstate;  // Auxiliary state
// Auxiliary index
    int    i;        // Aux index


    // Check correctness
    if(circ->emiss==1){
        cout << "add_BellP error #0: Photons already emitted. More photons can not be added at this stage" << endl;
        cout << "No new photons can be defined after a delay or dispersion elements " << endl;
        return -1;
    }

    // Reserve memory
    T.resize(4);

    // FIRST PAIR PACKET
    add=1;
    i=0;
    // Check if the packet described for this entry already exists
    while((i<npack)&&(add==1)){
    if((pack_list(1,i)==t1)&&(pack_list(2,i)==f1)&&(pack_list(3,i)==w1)) add=0;
        i++;
    }

    // Update the wave-packet information
    if(add==1){
        if(npack>=circ->ns) {
            cout << "add_BellP error #1: Not enough ns degrees of freedom! Needed: "<< npack+1 << endl;
            return -1;
        }
        T(0)=npack;
        pack_list(0,npack)=(double) npack;
        pack_list(1,npack)= t1;
        pack_list(2,npack)= f1;
        pack_list(3,npack)= w1;
        pack_list(4,npack)= t1;
        npack=npack+1;
    }else{
        // If exist just take the index
        T(0)=i-1;
    }

    // SECOND PAIR PACKET
    // Check if the packet described for this entry already exists
    add=1;
    i=0;
    // Check if the packet described for this entry already exists
    while((i<npack)&&(add==1)){
    if((pack_list(1,i)==t2)&&(pack_list(2,i)==f2)&&(pack_list(3,i)==w2)) add=0;
        i++;
    }

    // Update the wave-packet information
    if(add==1){
        if(npack>=circ->ns) {
            cout << "add_BellP error #2: Not enough ns degrees of freedom! Needed: "<< npack+1 << endl;
            return -1;
        }
        T(1)=npack;
        pack_list(0,npack)=(double) npack;
        pack_list(1,npack)= t2;
        pack_list(2,npack)= f2;
        pack_list(3,npack)= w2;
        pack_list(4,npack)= t2;
        npack=npack+1;
    }else{
        // If exist just take the index
        T(1)=i-1;
    }

    // Configure table according to chosen order
    ch.resize(3,2);
    ch  << ch1,  ch2,
          T(0), T(1),
          T(0), T(1);

    // Update state
    qdstate=new state(inpt->nlevel,inpt->maxket);
    qdstate->BellP(ch,kind,phi,circ);
    inpt->dproduct(qdstate);
    delete qdstate;

    return 0;
}

//----------------------------------------
//
//  Add photons to a device as generated by
//  an ideal Bell state (Path encoding)
//
//----------------------------------------
int qodev::add_BellP(int ch1, int ch2, char kind){
//  int    ch1;      // Chanel where the photons are added
//  int    ch2;      // Chanel where the photons are added
//  double phi;      // Non-ideal emission phase
//  char   kind;     // Kind of bell state +=Phi+/-=Phi-/p=Psi+/m=Psi-


     return this->add_BellP(ch1, ch2, kind, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
}

//----------------------------------------
//
//  Add photons to a device as generated by
//  an ideal Bell state (Polarization encoding)
//
//----------------------------------------
int qodev::add_Bell(int ch1, int ch2, char kind){
//  int    ch1;      // Chanel where the photons are added
//  int    ch2;      // Chanel where the photons are added
//  double phi;      // Non-ideal emission phase
//  char   kind;     // Kind of bell state +=Phi+/-=Phi-/p=Psi+/m=Psi-


     return this->add_Bell(ch1, ch2, kind, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
}


//----------------------------------------
//
//  Send the photons to the circuit specifying
//  the photons shape model. Configure the emitter.
//
//----------------------------------------
void qodev::send2circuit(){
//  Variables
    veci conversion;    // Packet conversion vector
    state *newinpt;     // New converted initial state


    circ->pack_list=pack_list;
    circ->npack=npack;
    if((circ->ns>1)&&(circ->npack>0)){
        // If we only have one possible packet emitter do not perform Gram-Schmidt well
        conversion=circ->emitter();
        newinpt=inpt->convert(conversion,circ);
        delete inpt;
        inpt=newinpt;
    }
    inpt->normalize();
}


//----------------------------------------
//
//  Changes the Gram-Schmidt order before
//  sending the photons to the circuit
//
//----------------------------------------
int qodev::repack(veci ipack){
//  veci ipack;         // New packet order
//  Variables
    matd new_list;      // New updated packet list
    state* newinpt;     // New input state updated to the new packet definitions
//  Auxiliary index
    int i;              // Aux index


    if(ipack.size()<npack){
        cout << "Repack error: more indexes needed" << endl;
        return -1;
    }
    for(i=0; i<npack;i++){
        pack_list(0,i)=ipack(i);
    }

    newinpt=inpt->convert(ipack,circ);
    delete inpt;
    inpt=newinpt;

    return 0;
}


//----------------------------------------
//
//  Returns the internal photon state.
//
//----------------------------------------
state *qodev::input(){
//  Variables
    state* aux;

    aux=inpt->clone();
    return aux;
}


//----------------------------------------
//
//  Returns the internal circuit object.
//
//----------------------------------------
qocircuit *qodev::circuit(){
//  Variables
    qocircuit* aux;

    aux=circ->clone();
    return aux;
}


//----------------------------------------
//
//  Overlapping probability of two wave packets
//
//----------------------------------------
double qodev::emitted_vis(int i,int j){
//  int i;      // First packet index
//  int j;      // Second packet index


    return circ->emitted_vis(i, j);
}


//----------------------------------------
//
//  Adds a random circuit
//
//----------------------------------------
void qodev::random_circuit(){


    circ->random_circuit();
}


//----------------------------------------
//
//  Adds a built-in NSX circuit
//
//----------------------------------------
void qodev::NSX(int i_ch1, int i_ch2, int i_ch3){
//  int   i_ch1              // NSX input channel 1
//  int   i_ch2              // NSX input channel 2
//  int   i_ch3              // NSX input channel 3


    circ->NSX(i_ch1,i_ch2,i_ch3);
}


//----------------------------------------
//
//  Adds a ideal beamsplitter to the device
//
//----------------------------------------
void qodev::beamsplitter(int i_ch1, int i_ch2, double theta, double phi){
//  int    i_ch1;             // Beamsplitter input channel 1
//  int    i_ch2;             // Beamsplitter input channel 2
//  double d_theta;           // Angle theta in degrees
//  double d_phi;             // Angle phi in degrees


    circ->beamsplitter(i_ch1, i_ch2, theta, phi);
}


//----------------------------------------
//
//  Adds a physical dieletric beamsplitter
//
//----------------------------------------
void qodev::dielectric(int i_ch1, int i_ch2, cmplx t, cmplx r){
//  int   i_ch1;              // Dielectric input channel 1
//  int   i_ch2;              // Dielectric input channel 2
//  cmplx t;                  // Transmission amplitude
//  cmplx r;                  // Reflection amplitude


    circ->dielectric(i_ch1, i_ch2, t, r);
}


//----------------------------------------
//
//  Adds a 2x2 ideal MMI
//
//----------------------------------------
void qodev::MMI2(int i_ch1, int i_ch2){
//  int   i_ch1               // MMI2 input channel 1
//  int   i_ch2               // MMI2 input channel 2


    circ->MMI2(i_ch1, i_ch2);
}


//----------------------------------------
//
//  Adds a swap gate between two channels
//
//----------------------------------------
void qodev::rewire(int i_ch1,int i_ch2){
//  int    i_ch1;             // Channel 1
//  int    i_ch2;             // Channel 2


    circ->rewire( i_ch1, i_ch2);
}


//----------------------------------------
//
//  Adds a phase_shifter to the device
//
//----------------------------------------
void qodev::phase_shifter(int i_ch, double phi){
//  int    i_ch               // Phase shifter input channel
//  double d_phi;             // Angle phi in degrees


    circ->phase_shifter( i_ch, phi);
}


//----------------------------------------
//
//  Adds a lossy medium
//
//----------------------------------------
void qodev::loss(int i_ch, double l){
//  int    i_ch           // Phase shifter input channel
//  double l;             // Loss probability


    circ->loss(i_ch, l);
}


//----------------------------------------
//
//  Adds a rotator to the device.
//
//----------------------------------------
void qodev::rotator(int i_ch, double d_theta, double d_phi){
//  int    i_ch               // Rotator input channel 1
//  double d_theta;           // Angle theta in degrees
//  double d_phi;             // Angle phi in degrees


    circ->rotator(i_ch, d_theta, d_phi);
}


//----------------------------------------
//
//  Adds a polarizing beamsplitter
//
//----------------------------------------
void qodev::polbeamsplitter(int i_ch1, int i_ch2, int P){
//  int   i_ch1               // Polarized beamsplitter input channel 1
//  int   i_ch2               // Polarized beamsplitter input channel 2
//  int   pol;                // Polarization to be switched


    circ->polbeamsplitter( i_ch1, i_ch2, P);
}


//----------------------------------------
//
//  Adds a half-waveplate
//
//----------------------------------------
void qodev::half(int i_ch, double alpha){
//  int    i_ch               // Waveplate input channel
//  double d_alpha;           // Angle of the polarization rotation


    circ->half( i_ch, alpha);
}


//----------------------------------------
//
//  Adds a quarter-waveplate
//
//----------------------------------------
void qodev::quarter(int i_ch, double alpha){
//  int    i_ch               // Waveplate input channel
//  double d_alpha;           // Angle of the polarization rotation


    circ->quarter(i_ch, alpha);
}


//----------------------------------------
//
// Flags a channel to be ignored
//
//----------------------------------------
void qodev::ignore(int i_ch){
//  int    i_ch;        // Chanel where detection is performed


    detector(i_ch,-2,1.0,0.0,0.0);
}


//----------------------------------------
//
//  Adds a detector
//
//----------------------------------------
void qodev::detector(int i_ch){
//  int    i_ch;        // Chanel where detection is performed


    detector(i_ch,-1,1.0,0.0,0.0);
}


//----------------------------------------
//
//  Adds a conditional detection
//
//----------------------------------------
void qodev::detector(int i_ch, int cond){
//  int    i_ch;        // Chanel where detection is performed
//  int    cond;        // Post selection condition. If cond<0 there is none.


    detector(i_ch,cond,1.0,0.0,0.0);
}


//----------------------------------------
//
//  Adds a general physical detector
//
//----------------------------------------
void qodev::detector(int i_ch, int cond, double eff, double blnk, double gamma){
//  int    i_ch;        // Chanel where detection is performed
//  int    cond;        // Post selection condition. If cond<0 there is none.
//  double eff;         // Efficiency of the detector
//  double blnk;        // Fraction of time blinking
//  double gamma;       // Dark counts rate


    if((circ->emiss==0)&&(circ->remdec()==1)) send2circuit();  // Send to circuit the photons.
    if((npack==0)&&(circ->losses==1)) circ->losses=2;          // The circuit has losses but it does not compute the full matrix yet
    circ->detector(i_ch, cond, eff, blnk, gamma);
}


//----------------------------------------
//
//  Adds noise to the output
//
//----------------------------------------
void qodev::noise(double stdev2){
//  double stdev2;     // Dispersion of the noise


    circ->noise(stdev2);
}


//----------------------------------------
//
//  Adds a delay using the packet definition given by def_packet
//
//----------------------------------------
void qodev::delay(int ch, double dt){
//  int    ch;   // Channel to be delayed
//  double dt;   // Amount of time to be delayed

    if (circ->emiss==0) send2circuit();
    circ->delay(ch, dt);
}


//----------------------------------------
//
//  Adds a phase to the photons in a channel
//  depending on the optical path dt and the
//  frequency of the photon packets
//
//----------------------------------------
void qodev::dispersion(int ch, double dt){
//  int ch;           // Channel
//  double dt;        // Time


    if (circ->emiss==0) send2circuit();
    circ->dispersion(ch, dt);
}


//----------------------------------------
//
//  Adds a phase to the photons in a channel
//  depending on the optical path dt and the
//  frequency of the photon packets for the
//  photons with a selected polarization
//
//----------------------------------------
void qodev::dispersion(int ch, int P, double dt){
//  int ch;           // Channel
//  in P;             // Photon polarization required to apply a phase shift
//  double dt;        // Time


    if (circ->emiss==0) send2circuit();
    circ->dispersion(ch, P, dt);
}


//----------------------------------------
//
// Apply the post-selection condition defined
// by the detectors (to ideal devices nm=1 and nd=1).
//
//----------------------------------------
state* qodev:: apply_condition(state *input, bool ignore){
//  state *input;       // Input state to apply the conditions
//  bool   ignore;      // Do we also ignore the ignored channels. False=No/True=Yes. Use with caution! ket collisions may happen if True.
//  Variables
    int   newnlevels;   // New number of levels
    int  *occ;          // Occupation
    mati  cond;         // Condition (of the detectors) matrix
    state *pselected;   // State after post-selection.
    state *ignored;     // State after removing ignored channels.
    projector *prj;     // Projector
//  Auxiliary index
    int    i;           // Aux index
    int    j;           // Aux index
    int    k;           // Aux index


    // Initialize post-selection condition
    // We are limited to a single ket projector
    cond.resize(4,circ->ncond);
    for(i=0;i<circ->ncond;i++){
        cond(0,i)=circ->det_def(0,i);
        cond(1,i)=0;
        cond(2,i)=0;
        cond(3,i)=circ->det_def(1,i);;

    }

    // Perform post-selection
    prj=new projector(circ->nlevel,1);
    if(circ->ncond>0){
        prj->add_term(1.0,cond,circ);
        pselected=input->post_selection(prj);
    }else{
        pselected=input->clone();
    }


    // Remove ignored channels
    if(ignore==true){
        // Reserve memory for the new state
        newnlevels=(pselected->nlevel)-(circ->nignored);
        ignored=new state(newnlevels,pselected->maxket);

        // Update visibility
        for(i=0;i<circ->nignored;i++){
            for(j=0;j<pselected->nlevel;j++){
                if(pselected->vis[j]==circ->ch_ignored(i)) pselected->vis[j]=-2;
            }
        }

        // Store the new visibility
        k=0;
        for(j=0;j<pselected->nlevel;j++){
            if(pselected->vis[j]>=0){
                    ignored->vis[k]=pselected->vis[j];
                    k=k+1;
            }
        }

        // Copy/store the visible levels
        for(i=0; i<pselected->nket;i++){
            occ=new int[newnlevels]();
            k=0;
            for(j=0;j<pselected->nlevel;j++){
                if(pselected->vis[j]>=0){
                    occ[k]=pselected->ket[i][j];
                    k=k+1;
                }
            }
            ignored->add_term(pselected->ampl[i],occ);
            delete occ;
        }
    }else{
        ignored=pselected->clone();
    }

    // Free memory
    delete prj;
    delete pselected;

    return ignored;
}


//----------------------------------------
//
//  Prints the packet configuration of the device
//
//----------------------------------------
void qodev::prnt_packets(){


    cout << "Table of times:" << endl;
    circ->emitted->prnt_times();
    cout << "Table of frequencies:" << endl;
    circ->emitted->prnt_freqs();
    cout << "Table of packets:" << endl;
    circ->emitted->prnt_packets();
}
