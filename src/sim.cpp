//======================================================================================================
// File sim.cpp
//
// SIMULATOR LIBRARY
//
// Copyright © 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved.
// The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the
// root directory of this source tree. Use of the source code in this file is only permitted under the
// terms of the licence in LICENCE.TXT.
//======================================================================================================
#include "sim.h"
#include <chrono>
#include <iostream>
#include <fstream>


//----------------------------------------
//
// Create a circuit "simulator".
// The maximum quantity of memory for operation and output is set by default.
//
//----------------------------------------
simulator::simulator(){


    backend=0;
    mem=DEFSIMMEM;
}


//----------------------------------------
//
// Create a circuit "simulator".
// The cores are selected using an integer.
// Less user friendly but simpler for the Python interface.
//
//----------------------------------------
simulator::simulator(int i_backend,int i_mem){
//  int i_backend        // Method selection
//  int i_mem            // Number of memory positions reserved.


    mem=i_mem;
    backend=i_backend;
}


//----------------------------------------
//
// Create a circuit "simulator".
// The maximum quantity of memory for operation and output is set explicitly.
//
//----------------------------------------
simulator::simulator(const char *i_back,int i_mem){
//  const char *i_back   // Method selected (in compilation time. Not dynamic!)
//  int i_mem            // Number of memory positions reserved


    mem=i_mem;
    switch (str2int(i_back)){
        case str2int("Direct"):
            backend=0;
            break;
        case str2int("DirectR"):
            backend=1;
            break;
        case str2int("Glynn"):
            backend=2;
            break;
        case str2int("GlynnR"):
            backend=3;
            break;
        default:
            cout << "Simulator error: No recognized backend" << endl;
            exit(0);
            break;
    }
}


//----------------------------------------
//
// Destroy circuit "simulator"
//
//----------------------------------------
simulator::~simulator(){
}


//----------------------------------------
//
// Simulation of a device
//
//----------------------------------------
p_bin *simulator::run(qodev *circuit){
//  qodev    circuit;      // Device to be simulated.
//  Variables
    state *output;         // Output state
    p_bin *outcome;        // Device outcomes and their probabilities
    p_bin *measured;       // Measured outcomes after going through physical detectors


    // Run simulation
    output=run(circuit->inpt,circuit->circ);

    // Store the raw statistic in a probability bin
    outcome= new p_bin(output->nlevel,mem);
    outcome->add_state(output);

    // Calculate the measured outcome including possible detector errors, etc
    measured=outcome->calc_measure(circuit->circ);


    // Free memory
    delete output;
    delete outcome;

    // Return result.
    return measured;
}


//--------------------------------------------------------------
//
// Calculate output state as function of the input state
//
//---------------------------------------------------------------
state *simulator::run( state *istate, qocircuit *qoc ){
//  state     *istate;           // Input state
//  qocircuit *qoc               // Circuit to be simulated
//  Variable
    state *empty_state;          // Empty state to return in case of bad init


    switch (backend)
    {
        case 0: // DirectF
            return DirectF(istate,qoc);
            break;
        case 1: // DirectR
            return DirectR(istate,qoc);
            break;
        case 2: // GlynnF
            return GlynnF(istate,qoc);
            break;
        case 3: // GlynnR
            return GlynnR(istate,qoc);
            break;
        default:
            empty_state=new state(qoc->nlevel,mem);
            return empty_state;
            break;
    }
}


//--------------------------------------------------------------
//
// Calculate the output amplitudes for the kets in olist as a
// function of the input state
//
//---------------------------------------------------------------
state *simulator::run( state *istate, ket_list* olist, qocircuit *qoc ){
//  state     *istate;      // Input state
//  ket_list  *olist;       // Output ket list
//  qocircuit *qoc          // Circuit to be simulated
//  Variable
    state *empty_state;     // Empty state to return in case of bad init


    switch (backend)
    {
        case 0: // Direct
            return DirectS(istate,olist,qoc);
            break;
        case 2: // Glynn
            return GlynnS(istate,olist,qoc);
            break;
        default:
            empty_state=new state(qoc->nlevel,mem);
            return empty_state;
            break;
    }
}


//--------------------------------------------------------------
//
// Direct method. Full distribution
//
//---------------------------------------------------------------
state *simulator::DirectF( state *istate, qocircuit *qoc ){
//  state     *istate;      // Input state
//  qocircuit *qoc          // Circuit to be simulated
//  Variables
    int    tocc;            // Number of photons present in ket.
    int    nlevel;          // Number of levels (qoc has this information, but it is put in this variable for easy access)
    double sqfact;          // Global sqrt factor to divide to get proper normalization
    long   int ncoef;       // Number of coefficients
    cmplx  coef;            // Coefficient for the transformation of a ket.
    veci   occs;            // Occupations of the states involved in the transformation of a ket.
    veci   ilev;            // sequence input levels
    state *ostate;          // Output state
//  Index
    int    iket;            // Index of input kets elements
    int    ilin;            // Index of input levels
    int    ilout;           // Index of output levels
    int    icoef;           // Index of coefficients for a particular ket
    int    iseq;            // Index of a position in a sequence
//  Auxiliary index
    int    j;               // Aux index


    //Set up variables and reserve memory
    nlevel=qoc->nlevel;
    ostate=new state(nlevel,mem);

    // Main loop
    // For each ket of a state calculate transformation rule.
    for(iket=0;iket<istate->nket;iket++){
    if(abs(istate->ampl[iket])>xcut){
        // Calculate number of photons of the ket
        // Normalization factors
        tocc=0;
        sqfact=1.0;
        for(ilin=0;ilin<nlevel;ilin++){
            tocc=tocc+istate->ket[iket][ilin];
            sqfact=sqfact*sqrt((double)factorial(istate->ket[iket][ilin]));
        }


        //Calculate from which input constructor we obtain the output one
        ncoef=pow(nlevel,tocc);
        ilev.resize(tocc);
        iseq=0;
        for(ilin=0;ilin<nlevel;ilin++){
            for(j=0;j<istate->ket[iket][ilin];j++){
                ilev(iseq)=ilin;
                iseq++;
            }
        }

        // Translate these sequences into actual coefficients and occupations
        occs.resize(nlevel);
        //For each coefficient and its occupation that correspond with one sequence
        for(icoef=0;icoef<ncoef;icoef++){
            coef=1.0;
            occs.setZero(nlevel);
            //Transform sequences int occupations and coefficients
            iseq=0;
            while((iseq<tocc)&&(abs(coef)>xcut)){
                ilin=ilev(iseq);
                ilout=(icoef/intpow(nlevel,iseq))%nlevel;
                occs(ilout)=occs(ilout)+1;
                coef=coef*qoc->circmtx(ilout,ilin)*sqrt((double)occs(ilout));
                iseq++;
            }
            // Normalize
            coef=istate->ampl[iket]*coef/sqfact;
            if(abs(coef)>xcut) ostate->add_term(coef,(int *)(occs.data()));
        }
    }
    }

    // Return output
    return ostate;
}


//--------------------------------------------------------------
//
// Direct method. Restricted distribution.
// We consider only the output kets with occupation zero or 1.
//
//---------------------------------------------------------------

state *simulator::DirectR( state *istate, qocircuit *qoc ){
//  state     *istate;      // Input state
//  qocircuit *qoc          // Circuit to be simulated
//  Variables
    int    tocc;            // Number of photons present in ket.
    int    nlevel;          // Number of levels (qoc has this information, but it is put in this variable for easy access)
    double sqfact;          // Global sqrt factor to divide to get proper normalization
    cmplx  coef;            // Coefficient for the transformation of a ket.
    veci   occs;            // Occupations of the states involved in the transformation of a ket.
    veci   ilev;            // sequence input levels
    string bitmask;         // Bit mask
    string perm;            // Permutations
    state *ostate;          // Output state
//  Index
    int    iket;            // Index of input kets elements
    int    ilin;            // Index of input levels
    int    ilout;           // Index of output levels
    int    iseq;            // Index of a position in a sequence
//  Auxiliary index
    int    i;               // Aux index
    int    j;               // Aux index


    //Set up variables and reserve memory
    nlevel=qoc->nlevel;
    ostate=new state(nlevel,mem);

    // Main loop
    // For each ket of a state calculate transformation rule.
    for(iket=0;iket<istate->nket;iket++){
    if(abs(istate->ampl[iket])>xcut){
        // Calculate number of photons of the ket
        // Normalization factors
        tocc=0;
        sqfact=1.0;
        for(ilin=0;ilin<nlevel;ilin++){
            tocc=tocc+istate->ket[iket][ilin];
            sqfact=sqfact*sqrt((double)factorial(istate->ket[iket][ilin]));
        }


        //Calculate from which input constructor we obtain the output one
        ilev.resize(tocc);
        iseq=0;
        for(ilin=0;ilin<nlevel;ilin++){
            for(j=0;j<istate->ket[iket][ilin];j++){
                ilev(iseq)=ilin;
                iseq++;
            }
        }

        // Translate these sequences into actual coefficients and occupations
        occs.resize(nlevel);
        //Initalize the bitmask of the occupation state. Bit of level j is 1 if j is occupied.
        bitmask.resize(0,0);
        bitmask.resize(tocc,1);
        bitmask.resize(nlevel,0);
        // Initalize permutations vector
        perm.resize(tocc,0);

        // For each occupation configuration
        sort(bitmask.begin(),bitmask.end());
        do{
            //Transform into photon - level sequence
            j=0;
            for (i=0; i<nlevel; i++){
                if (bitmask[i]){
                    perm[j]=(char)i;
                    j++;
                }
            }

            // Permute all the possible photon in level configurations that give the same occupation
            do{
                coef=1.0;
                occs.setZero(nlevel);
                iseq=0;
                while((iseq<tocc)&&(abs(coef)>xcut)){
                    ilin=ilev(iseq);
                    ilout=(int)perm[iseq];
                    occs(ilout)=occs(ilout)+1;
                    coef=coef*qoc->circmtx(ilout,ilin)*sqrt((double)occs(ilout));
                    iseq++;
                }

                // Normalize
                coef=istate->ampl[iket]*coef/sqfact;

                // Store
                if(abs(coef)>xcut) ostate->add_term(coef,(int *)(occs.data()));

            }while (next_permutation(perm.begin(), perm.end()));
        }while (next_permutation(bitmask.begin(), bitmask.end()));
    }
    }

    // Return output
    return ostate;
}


//--------------------------------------------------------------
//
// Permanent calculation method (using Glynn formula). Full distribution.
//
//---------------------------------------------------------------
state *simulator::GlynnF( state *istate, qocircuit *qoc ){
//  state     *istate;           // Input state
//  qocircuit *qoc               // Circuit to be simulated
//  Variables
    int    nph;                  // Number of photons present in input ket.
    int    nlevel;               // Number of levels (qoc has this information, but it is put in this variable for easy access)
    int   *pos;                  // Level where each photon is located. "Photon position"
    int   *occ;                  // Occupation
    cmplx  coef;                 // Coefficient for the transformation of a ket.
    cmplx  s;                    // Normalization coefficient of the input ket
    cmplx  t;                    // Normalization coefficient of the output ket
    matc   Ust;                  // Matrix to calculate the permanent
    state *ostate;               // Output state
//  Index
    int    iket;                 // Index of input kets elements
    int    ilin;                 // Index of input levels
    int    ilout;                // Index of output levels
    int    irow;                 // Row index of Ust
    int    icol;                 // Col index of Ust
//  Auxiliary index
    int    i;                    // Aux index
    int    j;                    // Aux index


    //Set up variables and reserve memory
    nlevel=qoc->nlevel;
    ostate=new state(nlevel,mem);

    // Main loop
    // For each ket of a state calculate transformation rule.
    for(iket=0;iket<istate->nket;iket++){
    if(abs(istate->ampl[iket])>xcut){
        // Calculate variables from the input ket
        nph=0;
        s=1.0;
        for(i=0;i<nlevel;i++){
            nph=nph+istate->ket[iket][i];
            s=s*(cmplx)factorial(istate->ket[iket][i]);
        }


        // Check all the possible outputs
        pos=new int[nph+1]();
        occ=new int[nlevel]();
        Ust.setZero(nph,nph);

        while(pos[0] < nlevel){
             // Calculate variables from the output ket
            t=1.0;
            for(j=0;j<nlevel;j++) occ[j]=0;
            for(j=0;j<nph;j++) {
                occ[pos[j]]=occ[pos[j]]+1;
                t=t*(cmplx)occ[pos[j]]; // This is the factorial implicitly.
            }


            // If the number of photons coincide (it always should)
            if(nph>0){
                // Create Ust
                icol=0;
                for(ilin=0;ilin<nlevel;ilin++){
                for(i=0;i<istate->ket[iket][ilin];i++){
                    irow=0;
                    for(ilout=0;ilout<nlevel;ilout++){
                    for(j=0;j<occ[ilout];j++){
                        Ust(irow,icol)=qoc->circmtx(ilout,ilin);
                        irow=irow+1;
                    }}
                    icol=icol+1;
                }}

                // Calculate coefficient
                coef=istate->ampl[iket]*glynn(Ust)/(sqrt(t)*sqrt(s));
            }else{
                coef=1.0;
             }


            // Store
            if(abs(coef)>xcut) ostate->add_term(coef,occ);


            // Obtain new photon level "position"
            pos[nph-1] += 1; // xxxxN -> xxxxN+1
            for (i = nph; i > 0; i -= 1) {
                if (pos[i] > nlevel - 1) // if number spilled over: xx0(n-1)xx
                {
                    pos[i - 1] += 1; // set xx1(n-1)xx
                    for (j = i; j <= nph; j += 1)
                        pos[j] = pos[j - 1]; // set xx11..1
                }
            }
        }

        // Free memory
        delete[] pos;
        delete[] occ;

    }}
    // Return output
    return ostate;
}


//--------------------------------------------------------------
//
// Permanent calculation method (using Glynn formula). Restricted distribution.
// We consider only output states with occupation zero or one.
//
//---------------------------------------------------------------
state *simulator::GlynnR( state *istate, qocircuit *qoc ){
//  state     *istate;         // Input state
//  qocircuit *qoc             // Circuit to be simulated
//  Variables
    int    nph;                // Number of photons present in input ket.
    int    nlevel;             // Number of levels (qoc has this information, but it is put in this variable for easy access)
    int   *occ;                // Occupation
    cmplx  coef;               // Coefficient for the transformation of a ket.
    cmplx  s;                  // Normalization coefficient of the input ket
    cmplx  t;                  // Normalization coefficient of the output ket
    matc   Ust;                // Matrix to calculate the permanent
    string bitmask;            // Bit mask
    state *ostate;             // Output state
//  Index
    int    iket;               // Index of input kets elements
    int    ilin;               // Index of input levels
    int    ilout;              // Index of output levels
    int    irow;               // Row index of Ust
    int    icol;               // Col index of Ust
//  Auxiliary index
    int    i;                  // Aux idex
    int    j;                  // Aux index


    //Set up variables and reserve memory
    nlevel=qoc->nlevel;
    ostate=new state(nlevel,mem);

    // Main loop
    // For each ket of a state calculate transformation rule.
    for(iket=0;iket<istate->nket;iket++){
    if(abs(istate->ampl[iket])>xcut){
        // Calculate variables from the input ket
        nph=0;
        s=1.0;
        for(i=0;i<nlevel;i++){
            nph=nph+istate->ket[iket][i];
            s=s*(cmplx)factorial(istate->ket[iket][i]);
        }

        // Translate these sequences into actual coefficients and occupations
        occ=new int[nlevel]();
        //Initalize the bitmask of the occupation state. Bit of level j is 1 if j is occupied.
        bitmask.resize(0,0);
        bitmask.resize(nph,1);
        bitmask.resize(nlevel,0);

        // For each occupation configuration
        sort(bitmask.begin(),bitmask.end());
        // For each occupation configuration
        do{
            Ust.setZero(nph,nph);

            // Calculate variables from the output ket
            t=1.0;
            for(j=0;j<nlevel;j++) {
                t=t*(cmplx)factorial(occ[j]);
                if (bitmask[j]){
                    occ[j]=1;
                }else{
                    occ[j]=0;
                }
            }

            // If the number of photons coincide (it always should)
            if(nph>0){
                // Create Ust
                icol=0;
                for(ilin=0;ilin<nlevel;ilin++){
                for(i=0;i<istate->ket[iket][ilin];i++){
                    irow=0;
                    for(ilout=0;ilout<nlevel;ilout++){
                    for(j=0;j<occ[ilout];j++){
                        Ust(irow,icol)=qoc->circmtx(ilout,ilin);
                        irow=irow+1;
                    }}
                    icol=icol+1;
                }}

                // Calculate coefficient
                coef=istate->ampl[iket]*glynn(Ust)/(sqrt(t)*sqrt(s));
            }else{
                coef=1.0;
            }

            // Store
            if(abs(coef)>xcut) ostate->add_term(coef,occ);

        } while (std::next_permutation(bitmask.begin(), bitmask.end()));

        // Free memory
        delete occ;

    }}
    // Return output
    return ostate;
}


//--------------------------------------------------------------
//
// Direct method for a restricted list of outputs.
//
//---------------------------------------------------------------
state *simulator::DirectS( state *istate, ket_list *olist, qocircuit *qoc ){
//  state     *istate;           // Input state
//  ket_list  *olist;            // Output ket list
//  qocircuit *qoc               // Circuit to be simulated
//  Variables
    int    tocc;                 // Number of photons present in the input  ket.
    int    nph;                  // Number of photons present in the output ket.
    int    nlevel;               // Number of levels (qoc has this information, but it is put in this variable for easy access)
    double sqfacti;              // Global sqrt factor to divide to get proper normalization of the input ket
    double sqfacto;              // Global sqrt factor to divide to get proper normalization of the output ket
    cmplx  coef;                 // Coefficient for the transformation of a ket.
    veci   occs;                 // Occupations of the states involved in the transformation of a ket.
    veci   ilev;                 // sequence input levels
    string perm;                 // Permutations
    state *ostate;               // Output state
//  Index
    int    iket;                 // Index of input kets elements
    int    oket;                 // Index of output kets elements
    int    ilin;                 // Index of input levels
    int    ilout;                // Index of output levels
    int    iseq;                 // Index of a position in a sequence
//  Auxiliary index
    int    i;                    // Aux index
    int    j;                    // Aux index
    int    k;                    // Aux index


    //Set up variables and reserve memory
    nlevel=qoc->nlevel;
    ostate=new state(nlevel,olist->nket+1);

    // Main loop
    // For each ket of a state calculate transformation rule.
    for(iket=0;iket<istate->nket;iket++){
    if(abs(istate->ampl[iket])>xcut){
        // Calculate variables of the input ket
        tocc=0;
        sqfacti=1.0;
        for(ilin=0;ilin<nlevel;ilin++){
            tocc=tocc+istate->ket[iket][ilin];
            sqfacti=sqfacti*sqrt((double)factorial(istate->ket[iket][ilin]));
        }

        //Calculate from which input constructor we obtain the output one
        ilev.resize(tocc);
        iseq=0;
        for(ilin=0;ilin<nlevel;ilin++){
            for(j=0;j<istate->ket[iket][ilin];j++){
                ilev(iseq)=ilin;
                iseq++;
            }
        }

        // For each output ket of the list
        for(oket=0;oket<olist->nket;oket++){
            // Calculate variables of the output key
            nph=0;
            sqfacto=1.0;
            for(ilout=0;ilout<nlevel;ilout++){
                nph=nph+olist->ket[oket][ilout];
                sqfacto=sqfacto*sqrt((double)factorial(olist->ket[oket][ilout]));
            }

            if(tocc==nph){
                // Translate these occupation into photon level sequences
                perm.resize(tocc,0);
                k=0;
                for (i=0; i<nlevel; i++){
                    for (j=0; j<olist->ket[oket][i]; j++){
                        perm[k]=(char)i;
                        k++;
                    }
                }
                sort(perm.begin(),perm.end());

                // For all photon level sequence that gives the current occupation
                do{
                    coef=1.0;
                    occs.setZero(nlevel);
                    //Transform sequences int occupations and coefficients
                    iseq=0;
                    while((iseq<tocc)&&(abs(coef)>xcut)){
                        ilin=ilev(iseq);
                        ilout=(int)perm[iseq];
                        coef=coef*qoc->circmtx(ilout,ilin);
                        iseq++;
                    }

                // Normalize
                coef=istate->ampl[iket]*sqfacto*coef/sqfacti;

                // Store
                if(abs(coef)>xcut) ostate->add_term(coef,olist->ket[oket]);

                }while (next_permutation(perm.begin(), perm.end()));
            }
        }
    }}

    // Return output
    return ostate;
}


//--------------------------------------------------------------
//
// Permanent calculation method (using Glynn formula) for a restricted
// of outputs.
//
//---------------------------------------------------------------
state *simulator::GlynnS( state *istate, ket_list *olist, qocircuit *qoc ){
//  state     *istate;         // Input state
//  ket_list  *olist;          // Output ket list
//  qocircuit *qoc             // Circuit to be simulated
//  Variables
    int    tocc;               // Number of photons present in input ket.
    int    nph;                // Number of photons present in output ket.
    int    nlevel;             // Number of levels (qoc has this information, but it is put in this variable for easy access)
    cmplx  coef;               // Coefficient for the transformation of a ket.
    cmplx  s;                  // Normalization coefficient of the input ket
    cmplx  t;                  // Normalization coefficient of the output ket
    matc   Ust;                // Matrix to calculate the permanent
    state *ostate;             // Output state
//  Index
    int    iket;               // Index of input kets elements
    int    oket;               // Index of output kets elements
    int    ilin;               // Index of input levels
    int    ilout;              // Index of output levels
    int    irow;               // Row index of Ust
    int    icol;               // Col index of Ust
//  Auxiliary index
    int    i;                  // Aux idex
    int    j;                  // Aux index


    //Set up variables and reserve memory
    nlevel=qoc->nlevel;
    ostate=new state(nlevel,olist->nket+1);

    // Main loop
    // For each ket of a state calculate transformation rule.
    for(iket=0;iket<istate->nket;iket++){
    if(abs(istate->ampl[iket])>xcut){
        // Calculate variables of the input ket
        tocc=0;
        s=1.0;
        for(i=0;i<nlevel;i++){
            tocc=tocc+istate->ket[iket][i];
            s=s*(cmplx)factorial(istate->ket[iket][i]);
        }

        // For each output ket of the list
        for(oket=0;oket<olist->nket;oket++){
            // Calculate variables of the output key
            nph=0;
            t=1.0;
            for(i=0;i<nlevel;i++){
                nph=nph+olist->ket[oket][i];
                t=t*(cmplx)factorial(olist->ket[oket][i]);
            }

            // If the number of photons coincide
            if(nph==tocc){
                if(nph>0){
                    // Create Ust
                    Ust.setZero(nph,nph);
                    icol=0;
                    for(ilin=0;ilin<nlevel;ilin++){
                    for(i=0;i<istate->ket[iket][ilin];i++){
                        irow=0;
                        for(ilout=0;ilout<nlevel;ilout++){
                        for(j=0;j<olist->ket[oket][ilout];j++){
                            Ust(irow,icol)=qoc->circmtx(ilout,ilin);
                            irow=irow+1;
                        }}
                        icol=icol+1;
                    }}

                    // Normalize
                    coef=istate->ampl[iket]*glynn(Ust)/(sqrt(t)*sqrt(s));
                }else{
                    coef=1.0;
                }

                // Store
                if(abs(coef)>xcut) ostate->add_term(coef,olist->ket[oket]);
            }
        }
    }}

    // Return output
    return ostate;
}


//----------------------------------------
//
// Clifford A Sampling method. Sampling from a device.
//
//----------------------------------------
p_bin *simulator::sample(qodev *circuit, int N){
//  qodev    circuit;      // Device to be simulated
//  int N;                 // Number of samples
//  Variables
    p_bin *outcome;        // Device outcomes and their probabilities
    p_bin *measured;       // Measured outcomes after going through physical detectors


    // Run simulation
    outcome=sample(circuit->inpt,circuit->circ,N);

    // Calculate the measured outcome including possible detector errors, etc
    measured=outcome->calc_measure(circuit->circ);

    // Free memory
    delete outcome;

    // Return result.
    return measured;
}


//--------------------------------------------------------------
//
// Clifford A Sampling method. Sampling from a circuit.
//
//---------------------------------------------------------------
p_bin *simulator::sample( state *istate, qocircuit *qoc ,int N){
//  state     *istate;      // Input state
//  qocircuit *qoc          // Circuit to be simulated
//  int N;                  // Number of samples
//  Variables
    int     maxnph;         // Maximum number of photons in the input kets
    int     nlevel;         // Number of levels (qoc has this information, but it is put in this variable for easy access)
    int    *nph;            // Number of photons present in the input ket.
    int    *ilist;          // Sequence of input photons
    int    *occ;            // Occupation
    int    *r;              // Sequence of output photons
    double  u;              // Random number of uniform distribution [0,1)
    double *w;              // Conditioned probability distribution
    double *iw;             // Conditioned probability distribution
    cmplx   p;              // Probability auxiliary variable
    string  bitmask;        // Bit mask
    string  perm;           // Permutations
    matc    Ust;            // Matrix to calculate the permanent
    p_bin  *obin;           // Output set of bins
//  Index
    int     iket;           // Index of input kets elements
    int     isample;        // Index of sample
    int     ilin;           // Index of input levels
    int     ilout;          // Index of output levels
    int     irow;           // Row index of Ust
    int     icol;           // Col index of Ust
//  Auxiliary index
    int     i;              // Aux index
    int     j;              // Aux index
    int     k;              // Aux index
    int     m;              // Aux index


    //Set up variables and reserve memory
    nlevel=qoc->nlevel;
    obin=new p_bin(nlevel,mem);

    // Process the input ket probabilities to perform random selection
    iw=new double[istate->nket]();
    for(iket=0;iket<istate->nket;iket++) iw[iket]=pow((double)abs(istate->ampl[iket]),2);
    for(iket=1;iket<istate->nket;iket++) iw[iket]=iw[iket]+iw[iket-1];
    for(iket=0;iket<istate->nket;iket++) iw[iket]=iw[iket]/iw[istate->nket-1];


    // Calculate number of photons in input ket
    maxnph=0;
    nph=new int[istate->nket]();
    for(iket=0;iket<istate->nket;iket++){
        nph[iket]=0;
        for(i=0;i<nlevel;i++){
                nph[iket]=nph[iket]+istate->ket[iket][i];
        }
        maxnph=max(maxnph,nph[iket]);
    }

    // Reserve memory
    r=new int[maxnph]();
    w=new double[nlevel]();
    ilist=new int[maxnph]();

    // Sample N end states
    for(isample=0;isample<N;isample++){
        // Select input ket (iket)
        u=urand();
        iket=0;
        while(u>iw[iket]){
            iket=iket+1;
        }

        //Initialize sequence of photons for the selected input ket.
        k=0;
        for(i=0;i<nlevel;i++){
            for(j=0;j<istate->ket[iket][i];j++){
                ilist[k]=i;
                k=k+1;
            }
        }

        // Generate a sequence of photons
        for(k=0;k<nph[iket];k++){
            // Calculate the conditional probability of each position
            // level where the next photon can be
            for(i=0;i<nlevel;i++){
                // Initialize variables
                r[k]=i;
                w[i]=0.0;

                // Reserve memory
                Ust.resize(k+1,k+1);
                perm.resize(k+1,0);

                //All the possible subsets of columns have to be considered
                bitmask.resize(0,0);
                bitmask.resize(k+1,1);
                bitmask.resize(nph[iket],0);
                sort(bitmask.begin(),bitmask.end());
                do{
                    // Also in all possible orders
                    j=0;
                    for (m=0; m<nph[iket]; m++){
                        if(bitmask[m]==1){
                            perm[j]=(char)ilist[m];
                            j=j+1;
                        }
                    }
                    sort(perm.begin(),perm.end());
                    do{
                        // Prepare Ust
                        icol=0;
                        for(ilin=0;ilin<k+1;ilin++){
                            irow=0;
                            for(ilout=0;ilout<=k;ilout++){

                                Ust(irow,icol)=qoc->circmtx(r[ilout],(int)perm[ilin]);
                                irow=irow+1;
                            }
                            icol=icol+1;
                        }
                        p=glynn(Ust);
                        w[i]=w[i]+pow(abs(p),2);
                    }while (next_permutation(perm.begin(), perm.end()));
                }while (next_permutation(bitmask.begin(), bitmask.end()));
            }

            // Calculated accumulated probability
            for(i=1;i<nlevel;i++) w[i]=w[i]+w[i-1];
            // Normalize
            for(i=0;i<nlevel;i++) w[i]=w[i]/w[nlevel-1];

            // Obtain level i with w[i] probability
            u=urand();
            i=0;
            while(u>w[i]){
                i=i+1;
            }
            // Update sequence
            r[k]=i;
        }

        // Calculate occupation for the present sequence of photons
        occ=new int[nlevel]();
        for(i=0;i<nph[iket];i++){
            occ[r[i]]=occ[r[i]]+1;
        }
        // Store the count
        obin->add_count(occ);
        delete[] occ;
    }

    // Free memory
    delete[] r;
    delete[] w;
    delete[] iw;
    delete[] nph;
    delete[] ilist;


    // Return output
    return obin;
}
