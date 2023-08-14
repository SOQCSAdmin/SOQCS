//======================================================================================================
// File qodev.cpp
//
// METACIRCUIT LIBRARY
//
// Copyright © 2023 National University of Ireland Maynooth, Maynooth University. All rights reserved.
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
    circ=new qocircuit(i_nch,1,1, 1, -1.0,0,0, false,'G');
    create_qodev(i_nph,circ->num_levels(),DEFSTATEDIM,1);
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
    circ=new qocircuit(i_nch,i_nm, 1, 1, -1.0, 0, 0, false,'G');
    create_qodev(i_nph,circ->num_levels(),DEFSTATEDIM,1);
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
    circ=new qocircuit(i_nch,i_nm,i_ns, 1, -1.0, clock, 0, false, ckind);
    create_qodev(i_nph, circ->num_levels(),DEFSTATEDIM,i_ns);
}


//-----------------------------------------------------
//
//  Creates a physical device with physical detectors
//
//-----------------------------------------------------
qodev::qodev(int i_nph, int i_nch, int i_nm, int i_ns, int i_np, double i_dtp, int clock, int i_R, bool loss, char ckind){
//  int    i_nph;   // Maximum number of photons to be simulated
//  int    i_nch;   // Number of channels
//  int    i_nm;    // Number of modes
//  int    i_ns;    // Number of packets
//  int    i_np;    // Number of periods.
//  double dtp;     // Period length
//  int    clock;   // Kind of detectors
//  int    i_R;     // Number of iterations to calculate blinking and dark counts
//  bool   loss;    // Are losses going to be calculated explicitly? True=Yes/False=No
//  char   ckind;   // Kind of packets  'G': Gaussian/'E': Exponential


    cfg_soqcs(i_nph);
    circ=new qocircuit(i_nch,i_nm,i_ns, i_np, i_dtp, clock, i_R, loss,ckind);
    create_qodev(i_nph, circ->num_levels(),DEFSTATEDIM,i_ns);
}


//-----------------------------------------------------
//
//  Creates a physical device with physical detectors
//
//-----------------------------------------------------
qodev::qodev(int i_nph, int i_nch, int i_nm, int i_ns, int i_np, double i_dtp, int clock, int i_R, bool loss, char ckind, int i_maxket){
//  int    i_nph;    // Maximum number of photons to be simulated
//  int    i_nch;    // Number of channels
//  int    i_nm;     // Number of modes
//  int    i_ns;     // Number of packets
//  int    i_np;     // Number of periods.
//  double dtp;      // Period length
//  int    clock;    // Kind of detectors
//  int    i_R;      // Number of iterations to calculate blinking and dark counts
//  bool   loss;     // Are losses going to be calculated explicitly? True=Yes/False=No
//  char   ckind;    // Kind of packets  'G': Gaussian/'E': Exponential
//  int    i_maxket; // Maximum number of bunches in the list


    cfg_soqcs(i_nph);
    circ=new qocircuit(i_nch,i_nm,i_ns, i_np, i_dtp, clock, i_R, loss,ckind);
    create_qodev(i_nph, circ->num_levels(),i_maxket,i_ns);
}


//----------------------------------------
//
//  Auxiliary private method to create a device
//
//----------------------------------------
void qodev::create_qodev(int i_nph, int i_level, int i_maxket, int i_ns){
//  int  i_nph;    // Maximum number of photons to be simulated
//  int  i_level;     // Number of levels to describe the state
//  int  i_maxket;    // Maximum number of bunches in the list
//  int  i_np;        // Number of packets
//  Variables
    int *occ;         // Level occupation


    npack=0;
    pack_list.resize(4,i_ns);
    inpt=new state(i_nph, i_level,i_maxket);
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
int qodev::concatenate(qodev *dev){
//  qodev *dev;      // Meta Circuit to be appended


    if((circ->emiss==0)&&(dev->circ->remdec()==0)) send2circuit();
    return circ->concatenate(dev->circ);
}


//------------------------------------------------
//
// Adds a gate using other device as the gate definition
//
//------------------------------------------------
int qodev::add_gate(veci chlist, qodev *dev){
//  veci chlist;     // List of channels where the gate is added
//  qodev *dev;      // Meta Circuit to be appended
//  Variables
    int    ch;       // Channel
    int    P;        // Polarization
    int    S;        // Packet number
    int    error;    // Error if adding a new term
    int   *occ;      // Occupation
    double t;        // Time
    double f;        // Frequency
    double w;        // Width/Characteristic time
    hterm  in_term;  // Input term definition
    veci   T;        // List of packet numbers
    state *aux;      // Auxiliary state. ( State that contain the ancilla definitions of dev translated to this device indexes )
    state *newinpt;  // New input state. New input state of this device after adding the ancilla values of the gate.
//  Auxiliary index
    int    i;        // Aux index
    int    j;        // Aux index


    // Check limitations
    if (chlist.size()!=dev->circ->nch){
        cout << "add_gate error (qodev) #1: The number of channels in the list has to be the same than in the gate circuit." << endl;
        return -1;
    }

    // Update packets with gate packets
    T.resize(dev->npack);
    for(i=0;i<dev->npack;i++){
        t=dev->pack_list(1,i);
        f=dev->pack_list(2,i);
        w=dev->pack_list(3,i);
        T(i)=add_photons(0,0,0,t,f,w);
        if (T(i)<0){
            cout << "add_gate error (qodev) #2: The number of packets of the device has been exceeded." << endl;
            return -1;
        }
    }

    // Update the input state information with the gate ancillas.
    // First, compute the gate contribution (ancilla values).
    in_term.resize(4,dev->circ->nlevel);
    for(i=0;i<dev->circ->nlevel;i++){
        ch=dev->circ->idx[i].ch;
        P=dev->circ->idx[i].m;
        S=dev->circ->idx[i].s;

        in_term(0,i) = chlist(ch);
        in_term(1,i) = P;
        if(S<dev->npack) in_term(2,i) = T(S); // OK. Because it can not be photons in packet that do not exist.
        else in_term(2,i) = 0;
        in_term(3,i) = dev->inpt->ket[0][i];

    }

    aux=new state(inpt->nph, inpt->nlevel,1);
    error=aux->add_term(1.0,in_term,circ);
    if(error<0){
        cout << "add_gate error (qodev) #3: Photons are being created at levels not defined in the circuit." << endl;
        delete aux;
        return -1;
    }

    // Second, add the gate contribution (ancilla values) to the input state.
    newinpt=new state(inpt->nph, inpt->nlevel,inpt->maxket);
    for(j=0;j<inpt->nket;j++){
        occ=new int[inpt->nlevel]();
        for(i=0;i<inpt->nlevel;i++){
            occ[i]=inpt->ket[j][i]+aux->ket[0][i];
        }

        newinpt->add_term(1.0,occ);
        delete occ;
    }

    //Free memory.
    delete aux;
    delete inpt;

    // Update input state.
    inpt=newinpt;

    // Update circuit mtx.
    if((circ->emiss==0)&&(circ->remdec()==dev->circ->ndetc))  send2circuit();
    return circ->add_gate(chlist,dev->circ);

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
    int    ip;       // Period index
    int    T;        // Packet number
    int    add;      // Add new packer 0=No/1=Yes. This way we check repetitions.
    int    error;    // An error happened when adding the term? >=0 No/-1=Yes
    double rt;       // Reduced time of the packet to its equivalent in the first period.
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
        return -1;
    }

    // Calculate period of the packet and the reduced
    // time in the leading period.
    if(circ->np>1){
        ip=floor((t+0.5*circ->dtp)/circ->dtp);
        rt=t-ip*circ->dtp;
    }else{
        ip=0;
        rt=t;
    }

    // Check if the packet described for this entry already exists
    add=1;
    i=0;
    while((i<npack)&&(add==1)){
        if( (abs(pack_list(1,i)-rt)<xcut)&&(abs(pack_list(2,i)-f)<xcut)&&(abs(pack_list(3,i)-w)<xcut)) add=0;
        i++;
    }

    // Update the wave-packet information
    if(add==1){
        if(npack>=circ->nsp){
            cout << "add_photons error #2: Not enough ns degrees of freedom! Needed at least: "<< npack+1 << endl;
            return -1;
        }
        // If not exists add it
        T=npack;
        pack_list(0,npack)=(double) npack;
        pack_list(1,npack)= rt;
        pack_list(2,npack)= f;
        pack_list(3,npack)= w;
        npack=npack+1;
    }else{
        // If exist just take the index
        T=i-1;
    }

    T=T+ip*circ->nsp;


    // Update the input state information
    // Compute the new state
    in_term.resize(4,1);
    in_term << ch,
               P,
               T,
               N;
    aux=new state(inpt->nph, inpt->nlevel,1);
    error=aux->add_term(1.0,in_term,circ);
    if(error<0){
        cout << "add_photons error #3: Photons are being created at levels not defined in the circuit." << endl;
        return -1;
    }

    newinpt=new state(inpt->nph, inpt->nlevel,inpt->maxket);
    for(j=0;j<inpt->nket;j++){
        occ=new int[inpt->nlevel]();
        for(i=0;i<inpt->nlevel;i++){
            occ[i]=inpt->ket[j][i]+aux->ket[0][i];
        }

        newinpt->add_term(inpt->ampl[j],occ);
        delete occ;
    }

    //Update internally the bunch
    delete aux;
    delete inpt;
    inpt=newinpt;

    return T;
}


//----------------------------------------
//
//  Add photons to a device as generated by a QD
//
//----------------------------------------
int qodev::add_QD(int ch1, int ch2,double i_t1, double i_f1, double i_w1, double i_t2, double i_f2, double i_w2, double S, double k, double tss, double thv, int cascade){
//  int    ch1;      // Chanel where the photons are added
//  int    ch2;      // Chanel where the photons are added
//  double t1;       // Photon 1 emission time
//  double f1;       // Photon 1 emission frequency
//  double w1;       // Photon 1 emission width or decay time depending on the packet shape model
//  double t2;       // Photon 2 emission time
//  double f2;       // Photon 2 emission frequency
//  double w2;       // Photon 2 emission width or decay time depending on the packet shape model
//  double S;        // FSS
//  double k;        // Signal to noise ration
//  double tss;      // Spin scattering characteristic time
//  double thv;      // Cross dephasing characteristic time
//  int cascade;     // Is the second photon part of a cascade (extra random emission time).
//  Variables
    int    ip1;      // Period index packet #1
    int    ip2;      // Period index packet #2
    int    add;      // Add new packer 0=No/1=Yes. This way we check repetitions.
    int    error;    // if >=0 operation dproduct is successful otherwise there is an error.
    double rt1;      // Reduced time of the packet #1 to its equivalent in the first period.
    double rt2;      // Reduced time of the packet #2 to its equivalent in the first period.
    double dt;       // Extra emission time
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

    // Calculate periods of the packets and their reduced
    // times in the leading period.
    if(circ->np>1){
        ip1=floor((i_t1+0.5*circ->dtp)/circ->dtp);
        rt1=i_t1-ip1*circ->dtp;

        ip2=floor((i_t2+0.5*circ->dtp)/circ->dtp);
        rt2=i_t2-ip2*circ->dtp;
    }else{
        ip1=0;
        rt1=i_t1;

        ip2=0;
        rt2=i_t2;
    }

    // Calculation of extra random time for cascades.
    if(circ->ckind=='G'){ // Gussian
        dt=erfi(2*urand()-1)/i_w1;
    }else{                // Exponential
        dt=i_w1*expi(urand());
    }

    // Reserve memory
    T.resize(2);

    // FIRST PAIR
    add=1;
    i=0;
    // Check if the packet described for this entry already exists
    while((i<npack)&&(add==1)){
        if( (abs(pack_list(1,i)-rt1)<xcut)&&(abs(pack_list(2,i)-i_f1)<xcut)&&(abs(pack_list(3,i)-i_w1)<xcut)) add=0;
        i++;
    }

    // Update the wave-packet information
    if(add==1){
        if(npack>=circ->nsp){
            cout << "add_QD error #1: Not enough ns degrees of freedom! Needed at least: "<< npack+1 << endl;
            return -1;
        }
        T(0)=npack;
        pack_list(0,npack)=(double) npack;
        pack_list(1,npack)= rt1;
        pack_list(2,npack)= i_f1;
        pack_list(3,npack)= i_w1;
        npack=npack+1;
    }else{
        // If exist just take the index
        T(0)=i-1;
    }
    T(0)=T(0)+ip1*circ->nsp;

    // SECOND PAIR
    // Check if the packet described for this entry already exists
    add=1;
    i=0;
    // Check if the packet described for this entry already exists
    while((i<npack)&&(add==1)){
        if( (abs(pack_list(1,i)-rt2)<xcut)&&(abs(pack_list(2,i)-i_f2)<xcut)&&(abs(pack_list(3,i)-i_w2)<xcut)) add=0;
        i++;
    }


    // Update the wave-packet information
    if(add==1){
        if(npack>=circ->nsp){
            cout << "add_QD error #2: Not enough ns degrees of freedom! Needed at least: "<< npack+1 << endl;
            return -1;
        }
        T(1)=npack;
        pack_list(0,npack)=(double) npack;
        if(cascade==0) pack_list(1,npack)= rt2;
        else pack_list(1,npack)= rt2+dt;
        pack_list(2,npack)=         i_f2;
        pack_list(3,npack)=         i_w2;
        npack=npack+1;
    }else{
        // If exist just take the index
        T(1)=i-1;
    }
    T(1)=T(1)+ip2*circ->nsp;

    // Configure table according to chosen order
    ch.resize(3,2);
    ch  << ch1,  ch2,
          T(0), T(1),
          T(0), T(1);

    // Update state
    qdstate=new state(inpt->nph, inpt->nlevel,inpt->maxket);
    qdstate->QD(ch,k,S,i_w2,tss,thv,circ);
    error=inpt->dproduct(qdstate);
    delete qdstate;

    return error;
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
    int    ip1;      // Period index packet #1
    int    ip2;      // Period index packet #2
    int    add;      // Add new packer 0=No/1=Yes. This way we check repetitions.
    int    error;    // if >=0 operation dproduct is successful otherwise there is an error.
    double rt1;      // Reduced time of the packet #1 to its equivalent in the first period.
    double rt2;      // Reduced time of the packet #2 to its equivalent in the first period.
    veci   T;        // Packet number
    mati   ch;       // QD configuration matrix
    hterm  in_term;  // Human readable input term
    state *qdstate;  // Auxiliary state
// Auxiliary index
    int    i;        // Aux index


    // Check correctness
    if(circ->emiss==1){
        cout << "add_Bell error #0: Photons already emitted. More photons can not be added at this stage" << endl;
        return -1;
    }

    // Calculate periods of the packets and their reduced
    // times in the leading period.
    if(circ->np>1){
        ip1=floor((t1+0.5*circ->dtp)/circ->dtp);
        rt1=t1-ip1*circ->dtp;

        ip2=floor((t2+0.5*circ->dtp)/circ->dtp);
        rt2=t2-ip2*circ->dtp;
    }else{
        ip1=0;
        rt1=t1;

        ip2=0;
        rt2=t2;
    }

    // Reserve memory
    T.resize(2);

    // FIRST PACKET
    add=1;
    i=0;
    // Check if the packet described for this entry already exists
    while((i<npack)&&(add==1)){
        if( (abs(pack_list(1,i)-rt1)<xcut)&&(abs(pack_list(2,i)-f1)<xcut)&&(abs(pack_list(3,i)-w1)<xcut)) add=0;
        i++;
    }

    // Update the wave-packet information
    if(add==1){
        if(npack>=circ->nsp){
            cout << "add_Bell error #1: Not enough ns degrees of freedom! Needed at least: "<< npack+1 << endl;
            return -1;
        }
        T(0)=npack;
        pack_list(0,npack)=(double) npack;
        pack_list(1,npack)= rt1;
        pack_list(2,npack)= f1;
        pack_list(3,npack)= w1;
        npack=npack+1;
    }else{
        // If exist just take the index
        T(0)=i-1;
    }
    T(0)=T(0)+ip1*circ->nsp;

    // SECOND PACKET
    // Check if the packet described for this entry already exists
    add=1;
    i=0;
    // Check if the packet described for this entry already exists
    while((i<npack)&&(add==1)){
        if( (abs(pack_list(1,i)-rt2)<xcut)&&(abs(pack_list(2,i)-f2)<xcut)&&(abs(pack_list(3,i)-w2)<xcut)) add=0;
        i++;
    }

    // Update the wave-packet information
    if(add==1){
        if(npack>=circ->nsp){
            cout << "add_Bell error #2: Not enough ns degrees of freedom! Needed at least: "<< npack+1 << endl;
            return -1;
        }
        T(1)=npack;
        pack_list(0,npack)=(double) npack;
        pack_list(1,npack)= rt2;
        pack_list(2,npack)= f2;
        pack_list(3,npack)= w2;
        npack=npack+1;
    }else{
        // If exist just take the index
        T(1)=i-1;
    }
    T(1)=T(1)+ip2*circ->nsp;

    // Configure table according to chosen order
    ch.resize(2,2);
    ch  << ch1,  ch2,
          T(0), T(1);

    // Update state
    qdstate=new state(inpt->nph, inpt->nlevel,inpt->maxket);
    qdstate->Bell(ch,kind,phi,circ);
    error=inpt->dproduct(qdstate);
    delete qdstate;

    return error;
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
    int    ip1;      // Period index packet #1
    int    ip2;      // Period index packet #2
    int    add;      // Add new packer 0=No/1=Yes. This way we check repetitions.
    int    error;    // if >=0 operation dproduct is successful otherwise there is an error.
    double rt1;      // Reduced time of the packet #1 to its equivalent in the first period.
    double rt2;      // Reduced time of the packet #2 to its equivalent in the first period.
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

    // Calculate periods of the packets and their reduced
    // times in the leading period.
    if(circ->np>1){
        ip1=floor((t1+0.5*circ->dtp)/circ->dtp);
        rt1=t1-ip1*circ->dtp;

        ip2=floor((t2+0.5*circ->dtp)/circ->dtp);
        rt2=t2-ip2*circ->dtp;
    }else{
        ip1=0;
        rt1=t1;

        ip2=0;
        rt2=t2;
    }

    // Reserve memory
    T.resize(2);

    // FIRST PAIR PACKET
    add=1;
    i=0;
    // Check if the packet described for this entry already exists
    while((i<npack)&&(add==1)){
        if( (abs(pack_list(1,i)-rt1)<xcut)&&(abs(pack_list(2,i)-f1)<xcut)&&(abs(pack_list(3,i)-w1)<xcut)) add=0;
        i++;
    }

    // Update the wave-packet information
    if(add==1){
        if(npack>=circ->nsp){
            cout << "add_BellP error #1: Not enough ns degrees of freedom! Needed at least: "<< npack+1 << endl;
            return -1;
        }
        T(0)=npack;
        pack_list(0,npack)=(double) npack;
        pack_list(1,npack)= rt1;
        pack_list(2,npack)= f1;
        pack_list(3,npack)= w1;
        npack=npack+1;
    }else{
        // If exist just take the index
        T(0)=i-1;
    }
    T(0)=T(0)+ip1*circ->nsp;

    // SECOND PAIR PACKET
    // Check if the packet described for this entry already exists
    add=1;
    i=0;
    // Check if the packet described for this entry already exists
    while((i<npack)&&(add==1)){
        if( (abs(pack_list(1,i)-rt2)<xcut)&&(abs(pack_list(2,i)-f2)<xcut)&&(abs(pack_list(3,i)-w2)<xcut)) add=0;
        i++;
    }

    // Update the wave-packet information
    if(add==1){
        if(npack>=circ->nsp){
            cout << "add_BellP error #2: Not enough ns degrees of freedom! Needed at least: "<< npack+1 << endl;
            return -1;
        }
        T(1)=npack;
        pack_list(0,npack)=(double) npack;
        pack_list(1,npack)= rt2;
        pack_list(2,npack)= f2;
        pack_list(3,npack)= w2;
        npack=npack+1;
    }else{
        // If exist just take the index
        T(1)=i-1;
    }
    T(1)=T(1)+ip2*circ->nsp;

    // Configure table according to chosen order
    ch.resize(2,2);
    ch  << ch1,  ch2,
          T(0), T(1);

    // Update state
    qdstate=new state(inpt->nph, inpt->nlevel,inpt->maxket);
    qdstate->BellP(ch,kind,phi,circ);
//    qdstate->prnt_state(1,1,false,circ);
    error=inpt->dproduct(qdstate);
    delete qdstate;

    return error;
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


    if(ipack.size()>circ->nsp){
        cout << "Repack error: too many indexes" << endl;
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
int  qodev::NSX(int i_ch1, int i_ch2, int i_ch3){
//  int   i_ch1              // NSX input channel 1
//  int   i_ch2              // NSX input channel 2
//  int   i_ch3              // NSX input channel 3


    return circ->NSX(i_ch1,i_ch2,i_ch3);
}


//----------------------------------------
//
//  Adds a ideal beamsplitter to the device
//
//----------------------------------------
int qodev::beamsplitter(int i_ch1, int i_ch2, double theta, double phi){
//  int    i_ch1;             // Beamsplitter input channel 1
//  int    i_ch2;             // Beamsplitter input channel 2
//  double d_theta;           // Angle theta in degrees
//  double d_phi;             // Angle phi in degrees


    return circ->beamsplitter(i_ch1, i_ch2, theta, phi);
}


//----------------------------------------
//
//  Adds a physical dieletric beamsplitter
//
//----------------------------------------
int qodev::dielectric(int i_ch1, int i_ch2, cmplx t, cmplx r){
//  int   i_ch1;              // Dielectric input channel 1
//  int   i_ch2;              // Dielectric input channel 2
//  cmplx t;                  // Transmission amplitude
//  cmplx r;                  // Reflection amplitude


    return circ->dielectric(i_ch1, i_ch2, t, r);
}


//----------------------------------------
//
//  Adds a 2x2 ideal MMI
//
//----------------------------------------
int qodev::MMI2(int i_ch1, int i_ch2){
//  int   i_ch1               // MMI2 input channel 1
//  int   i_ch2               // MMI2 input channel 2


    return circ->MMI2(i_ch1, i_ch2);
}


//----------------------------------------
//
//  Adds a swap gate between two channels
//
//----------------------------------------
int qodev::rewire(int i_ch1,int i_ch2){
//  int    i_ch1;             // Channel 1
//  int    i_ch2;             // Channel 2


    return circ->rewire( i_ch1, i_ch2);
}


//----------------------------------------
//
//  Adds a phase_shifter to the device
//
//----------------------------------------
int qodev::phase_shifter(int i_ch, double phi){
//  int    i_ch               // Phase shifter input channel
//  double d_phi;             // Angle phi in degrees


    return circ->phase_shifter( i_ch, phi);
}


//----------------------------------------
//
//  Adds a lossy medium
//
//----------------------------------------
int qodev::loss(int i_ch, double l){
//  int    i_ch           // Phase shifter input channel
//  double l;             // Loss probability


    return circ->loss(i_ch, l);
}


//----------------------------------------
//
//  Adds a rotator to the device.
//
//----------------------------------------
int qodev::rotator(int i_ch, double d_theta, double d_phi){
//  int    i_ch               // Rotator input channel 1
//  double d_theta;           // Angle theta in degrees
//  double d_phi;             // Angle phi in degrees


    return circ->rotator(i_ch, d_theta, d_phi);
}


//----------------------------------------
//
//  Adds a polarizing beamsplitter
//
//----------------------------------------
int qodev::pol_beamsplitter(int i_ch1, int i_ch2, int P){
//  int   i_ch1               // Polarized beamsplitter input channel 1
//  int   i_ch2               // Polarized beamsplitter input channel 2
//  int   P;                  // Polarization to be switched


    return circ->pol_beamsplitter( i_ch1, i_ch2, P);
}


//----------------------------------------
//
//  Adds a polarizing phase shifter
//
//----------------------------------------
int qodev::pol_phase_shifter(int i_ch, int P, double phi){
//  int     i_ch             // Polarized phase shifter input channel
//  int     P;               // Polarization
//  double  phi;             // Phase


    return circ->pol_phase_shifter(i_ch, P, phi);
}


//----------------------------------------
//
//  Adds a polarization filter
//
//----------------------------------------
int qodev::pol_filter(int i_ch, int P){
//  int     i_ch             // Filter input channel
//  int     P;               // Polarization


    return circ->pol_filter(i_ch, P);
}


//----------------------------------------
//
//  Adds a half-waveplate
//
//----------------------------------------
int qodev::half(int i_ch, double alpha){
//  int    i_ch               // Waveplate input channel
//  double d_alpha;           // Angle of the polarization rotation


    return circ->half( i_ch, alpha);
}


//----------------------------------------
//
//  Adds a quarter-waveplate
//
//----------------------------------------
int qodev::quarter(int i_ch, double alpha){
//  int    i_ch               // Waveplate input channel
//  double d_alpha;           // Angle of the polarization rotation


    return circ->quarter(i_ch, alpha);
}


//----------------------------------------
//
// Flags a channel to be ignored
//
//----------------------------------------
int qodev::ignore(int i_ch){
//  int    i_ch;        // Chanel where detection is performed


    return detector(i_ch,-2,-1,-1,-1,1.0,0.0,0.0);
}


//----------------------------------------
//
//  Adds a detector
//
//----------------------------------------
int qodev::detector(int i_ch){
//  int    i_ch;        // Chanel where detection is performed


    return detector(i_ch,-1,1.0,0.0,0.0);
}


//----------------------------------------
//
//  Adds a conditional detection
//
//----------------------------------------
int qodev::detector(int i_ch, int cond){
//  int    i_ch;        // Chanel where detection is performed
//  int    cond;        // Post selection condition. If cond<0 there is none.


    return detector(i_ch,cond,-1,-1,-1,1.0,0.0,0.0);
}


//----------------------------------------
//
//  Adds a general physical detector
//
//----------------------------------------
int qodev::detector(int i_ch, int cond, double eff, double blnk, double gamma){
//  int    i_ch;        // Chanel where detection is performed
//  int    cond;        // Post selection condition. If cond<0 there is none.
//  double eff;         // Efficiency of the detector
//  double blnk;        // Fraction of time blinking
//  double gamma;       // Dark counts rate

    return detector(i_ch,cond,-1,-1,-1, eff, blnk, gamma);
}


//----------------------------------------------------------------------------------------------------------
//
//  Adds a general physical detector with a window of detection (and conditional detection by polarization)
//
//----------------------------------------------------------------------------------------------------------
int qodev::detector(int i_ch, int cond, int pol, int mpi, int mpo, double eff, double blnk, double gamma){
//  int    i_ch;        // Chanel where detection is performed
//  int    cond;        // Post selection condition. If cond<0 there is none.
//  int    pol;         // Post selection polarization condition. If pol<0 there is none.
//  int    mpi;         // Initial detection period (if -1 takes the first one as default)
//  int    mpo;         // Last detection period (if -1 takes the last one as default)
//  double eff;         // Efficiency of the detector
//  double blnk;        // Fraction of time blinking
//  double gamma;       // Dark counts rate


    if((circ->emiss==0)&&(circ->remdec()==1)) send2circuit();           // Send to circuit the photons.
    if((npack==0)&&(circ->losses==1)) circ->losses=2;                   // The circuit has losses but it does not compute the full matrix yet
    return circ->detector(i_ch, cond, pol, mpi, mpo, eff, blnk, gamma); // Adds a detector to the circuit
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
int qodev::delay(int ch){
//  int    ch;   // Channel to be delayed
//  double dt;   // Amount of time to be delayed


    if (circ->emiss==0) send2circuit();
    return circ->delay(ch);;

}


//----------------------------------------
//
//  Adds a phase to the photons in a channel
//  depending on the optical path dt and the
//  frequency of the photon packets
//
//----------------------------------------
int qodev::dispersion(int ch, double dt){
//  int ch;           // Channel
//  double dt;        // Time


    if (circ->emiss==0) send2circuit();
    return circ->dispersion(ch, dt);
}


//----------------------------------------
//
//  Adds a phase to the photons in a channel
//  depending on the optical path dt and the
//  frequency of the photon packets for the
//  photons with a selected polarization
//
//----------------------------------------
int qodev::dispersion(int ch, int P, double dt){
//  int ch;           // Channel
//  in P;             // Photon polarization required to apply a phase shift
//  double dt;        // Time


    if (circ->emiss==0) send2circuit();
    return circ->dispersion(ch, P, dt);
}


//----------------------------------------
//
// Apply the post-selection condition defined
// by the detectors (to ideal devices ns=1).
//
//----------------------------------------
state* qodev:: apply_condition(state *input, bool ignore){
//  state *input;       // Input state to apply the conditions
//  bool   ignore;      // Do we also ignore the ignored channels. False=No/True=Yes. Use with caution! ket collisions may happen if True.
//  Variables
    bool   print;       // A warning has been printed. True=Yes/False=No
    int    newnlevels;  // New number of levels
    int    nadded;      // Number of kets after removing ignored channels.
    int    idx;         // Index of the stored state.
    int   *occ;         // Occupation
    mati   cond;        // Condition (of the detectors) matrix
    state *pselected;   // State after post-selection.
    state *ignored;     // State after removing ignored channels.
    projector *prj;     // Projector
//  Auxiliary index
    int    i;           // Aux index
    int    j;           // Aux index
    int    k;           // Aux index
    int ch;


    // Initialize post-selection condition
    // We are limited to a single ket projector
    cond.resize(4,circ->nm*circ->ncond);
    k=0;
    for(i=0;i<circ->ncond;i++){
        for(j=0;j<circ->nm;j++){
            cond(0,k)=circ->det_def(0,i);
            cond(1,k)=j;
            cond(2,k)=0;
            if ((j==circ->det_def(2,i))||(circ->det_def(2,i)<0)) cond(3,k)=circ->det_def(1,i);
            else cond(3,k)=0;
            k=k+1;
        }
    }

    // Perform post-selection
    prj=new projector(input->nph,circ->nlevel,1);
    if(circ->ncond>0){
        prj->add_term(1.0,cond,circ);
        pselected=input->post_selection(prj);
    }else{
        pselected=input->clone();
    }


    // Remove ignored channels
    if(ignore==true){
        // Reserve memory for the new state
        newnlevels=(pselected->nlevel)-(circ->nm*circ->nignored);
        ignored=new state(input->nph, newnlevels,pselected->maxket);

        // Update visibility
        for(i=0;i<circ->nignored;i++){
            for(j=0;j<pselected->nlevel;j++){
                ch=circ->idx[pselected->vis[j]].ch;
                if(ch==circ->ch_ignored(i)) pselected->vis[j]=-2;
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
        nadded=0;
        print=false;
        for(i=0; i<pselected->nket;i++){
            occ=new int[newnlevels]();
            k=0;
            for(j=0;j<pselected->nlevel;j++){
                if(pselected->vis[j]>=0){
                    occ[k]=pselected->ket[i][j];
                    k=k+1;
                }
            }
            idx=ignored->add_term(pselected->ampl[i],occ);
            if((idx!=nadded)&&(print==false)){
                cout << "Apply_condition: Warning ignored lead to collision!" << endl;
                print=true;
            }
            nadded=nadded+1;
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
// Apply the post-selection condition defined
// by the detectors (to ideal devices ns=1).
//
//----------------------------------------
state* qodev:: apply_condition(state *input){
//  state *input;       // Input state to apply the conditions


    return apply_condition(input, true);
}


//----------------------------------------
//
// Initialize the device with qubit values
// ( Path encoding version)
//
//----------------------------------------
void qodev::qubits(veci qinit, mati qmap){
//  veci qinit;    // Qubit values
//  mati qmap;     // Correspondence between qubits and channels (Path encoding)
//  Variables
    int i;


    for (i=0; i<qinit.size(); i++){
        if(qinit(i)==0){
            add_photons(0, qmap(0,i));
            add_photons(1, qmap(1,i));
        }else{
            add_photons(1, qmap(0,i));
            add_photons(0, qmap(1,i));
        }
    }
}


//----------------------------------------
//
// Initialize the device with qubit values
// ( Polarization encoding version)
//
//----------------------------------------
void qodev::pol_qubits(veci qinit, veci qmap){
//  veci qinit;    // Qubit values
//  veci qmap;     // Correspondence between qubits and channels (Polarization encoding)
//  Variables
    int i;


    for (i=0; i<qinit.size(); i++){
        add_photons(1,qmap(0,i), qinit(i), 0.0, 0.0, 0.0);
    }
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
