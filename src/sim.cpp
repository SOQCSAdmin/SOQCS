//======================================================================================================
// File sim.cpp
//
// SIMULATOR LIBRARY
//
// Copyright © 2023 National University of Ireland Maynooth, Maynooth University. All rights reserved.
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


    mem=DEFSIMMEM;
}


//----------------------------------------
//
// Create a circuit "simulator".
//
//----------------------------------------
simulator::simulator(int i_mem){
//  int i_mem            // Number of memory positions reserved.


    mem=i_mem;
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
p_bin *simulator::run(qodev *circuit, int method){
//  qodev    circuit;      // Device to be simulated.
//  Variables
    state *output;         // Output state
    p_bin *outcome;        // Device outcomes and their probabilities
    p_bin *measured;       // Measured outcomes after going through physical detectors


    // Run simulation
    output=run(circuit->inpt,circuit->circ, method);

    // Store the raw statistic in a probability bin
    outcome= new p_bin(output->nph,output->nlevel,mem);
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
state *simulator::run( state *istate, qocircuit *qoc, int method ){
//  state     *istate;           // Input state
//  qocircuit *qoc               // Circuit to be simulated
//  Variable
    state *empty_state;          // Empty state to return in case of bad init


    switch (method)
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
            cout << "Run error: No recognized backend." << endl;
            empty_state=new state(istate->nph,istate->nlevel,mem);
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
state *simulator::run( state *istate, ket_list* olist, qocircuit *qoc, int method ){
//  state     *istate;      // Input state
//  ket_list  *olist;       // Output ket list
//  qocircuit *qoc          // Circuit to be simulated
//  Variable
    state *empty_state;     // Empty state to return in case of bad init


    switch (method)
    {
        case 0: // Direct
            return DirectS(istate,olist,qoc);
            break;
        case 2: // Glynn
            return GlynnS(istate,olist,qoc);
            break;
        default:
            empty_state=new state(istate->nph, istate->nlevel,mem);
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
    int    index;           // Ket list position where a new term of the output state is stored.
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
    ostate=new state(istate->nph, nlevel,mem);

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

            // Store
            if(abs(coef)>xcut){
                index=ostate->add_term(coef,(int *)(occs.data()));
                if(index<0){
                    cout << "Simulator(DirectF): Warning! Simulation canceled because the memory limit has been exceeded.  Increase *mem* for more memory." << endl;
                    return ostate;
                }
            }
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
    int    index;           // Ket list position where a new term of the output state is stored.
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
    ostate=new state(istate->nph,nlevel,mem);

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
                if(abs(coef)>xcut){
                    index=ostate->add_term(coef,(int *)(occs.data()));
                    if(index<0){
                        cout << "Simulator(DirectR): Warning! Simulation canceled because the memory limit has been exceeded.  Increase *mem* for more memory." << endl;
                        return ostate;
                    }
                }


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
    int    index;                // Ket list position where a new term of the output state is stored.
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
    ostate=new state(istate->nph,nlevel,mem);

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
            if(abs(coef)>xcut){
                index= ostate->add_term(coef,occ);
                if(index<0){
                    cout << "Simulator(GlynnF): Warning! Simulation canceled because the memory limit has been exceeded.  Increase *mem* for more memory." << endl;
                    // Free memory
                    delete[] pos;
                    delete[] occ;
                    // Return partial calculation
                    return ostate;
                }
            }


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
    int    index;              // Ket list position where a new term of the output state is stored.
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
    ostate=new state(istate->nph,nlevel,mem);

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
            if(abs(coef)>xcut){
                index= ostate->add_term(coef,occ);
                if(index<0){
                    cout << "Simulator(GlynnR): Warning! Simulation canceled because the memory limit has been exceeded.  Increase *mem* for more memory." << endl;
                    // Free memory
                    delete occ;
                    // Return partial calculation
                    return ostate;
                }
            }

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
    int    index;                // Ket list position where a new term of the output state is stored.
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
    ostate=new state(istate->nph,nlevel,olist->nket+1);

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
                if(abs(coef)>xcut){
                    index= ostate->add_term(coef,olist->ket[oket]);
                    if(index<0){
                        cout << "Simulator(DirectS): Warning! Simulation canceled because the memory limit has been exceeded.  Increase *mem* for more memory." << endl;
                        return ostate;
                    }
                }


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
    int    index;              // Ket list position where a new term of the output state is stored.
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
    ostate=new state(istate->nph,nlevel,olist->nket+1);

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
                if(abs(coef)>xcut){
                    index= ostate->add_term(coef,olist->ket[oket]);
                    if(index<0){
                        cout << "Simulator(GlynnS): Warning! Simulation canceled because the memory limit has been exceeded.  Increase *mem* for more memory." << endl;
                        return ostate;
                    }
                }

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
    int     index;          // Bin position where a new output count is stored.
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
    obin=new p_bin(maxnph, nlevel,mem);

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
        index=obin->add_count(occ);
        delete[] occ;

        if(index<0){
            cout << "Sample: Warning! Sampling canceled because the memory limit has been exceeded.  Increase *mem* for more memory." << endl;

            // Free memory
            delete[] r;
            delete[] w;
            delete[] iw;
            delete[] nph;
            delete[] ilist;

            // Return partial calculation.
            return obin;
        }
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


//--------------------------------------------------------------
//
// Metropolis sampling method. Sampling from a circuit. ( Device version )
//
//---------------------------------------------------------------
tuple<p_bin*, double> simulator::metropolis( qodev *circuit ,int method, int N, int Nburn, int Nthin){
//  qodev *circuit;     // Device to be samples.
//  int    method;      // Sampling method
                            // 0: f classical
                            // 1: g Uniform
                            // 2: g Classical
                            // 3: f classical ( Restricted )
                            // 4: g Uniform  ( Restricted )
                            // 5: g Classical  ( Restricted )
//  int N;                  // Number of samples
//  int Nburn;              // Number of initial samples before arriving to stationary distribution.
//  int Nthin;              // Number of thining samples to avoid correlation.
//  Variables
    double p;           // Success probability
    p_bin *outcome;     // Device outcomes and their probabilities
    p_bin *measured;    // Measured outcomes after going through physical detectors


    // Run simulation
    tie(outcome,p)=metropolis(circuit->inpt,circuit->circ,method,N,Nburn,Nthin);

    // Calculate the measured outcome including possible detector errors, etc
    measured=outcome->calc_measure(circuit->circ);

    // Free memory
    delete outcome;

    // Return result.
    return {measured,p};
}


//--------------------------------------------------------------
//
// Metropolis sampling method. Sampling from a circuit. ( Circuit version )
//
//---------------------------------------------------------------
tuple<p_bin*, double> simulator::metropolis( state *istate, qocircuit *qoc ,int method, int N, int Nburn, int Nthin){
//  state     *istate;      // Input state
//  qocircuit *qoc;         // Circuit to be samples.
//  int        method;      // Sampling method
                                // 0: f classical
                                // 1: g Uniform
                                // 2: g Classical
                                // 3: f classical ( Restricted )
                                // 4: g Uniform  ( Restricted )
                                // 5: g Classical  ( Restricted )
//  int N;                  // Number of samples
//  int Nburn;              // Number of initial samples before arriving to stationary distribution.
//  int Nthin;              // Number of thinning samples to avoid correlation.
//  Variables
    bool   gral;            // True='Unrestricted sampling'/'False=Restricted sampling'
    bool   uniform;         // True='Uniform sampling'/'False=Circuit classical distribution sampling'
    bool   classic;         // Classical output True='Yes'/False='No'
    int    Neff;            // Effective number of samples calculated. ( Not every sample is stored by different reasons)
    int    nph;             // Number of photons present in the input ket.
    int    nlevel;          // Number of levels (qoc has this information, but it is put in this variable for easy access)
    int    index;           // Bin position where a new output count is stored.
    int    *ilist;          // Sequence of input photons
    int    *occ;            // Occupation
    int    *r;              // Sequence of output photons
    double  T;              // Sample acceptance probability
    double  s;              // n1!n2!...nl! of the input
    double  t;              // n1!n2!...nl! of the proposed sample
    double  p=1.0;          // Probability of the sample.
    double  p_old=1.0;      // Probability of the previous sample.
    double  pc;             // Classical probability of the sample.
    double  pc_old;         // Classical probability of the previous sample.
    matc    Ust;            // Matrix to calculate the permanent
    p_bin  *obin;           // Output set of bins
//  Index
    int     isample;        // Index of accepted samples
    int     istored;        // Index of stored samples
    int     ilin;           // Index of input levels
    int     ilout;          // Index of output levels
    int     irow;           // Row index of Ust
    int     icol;           // Col index of Ust
//  Auxiliary index
    int     i;              // Aux index
    int     j;              // Aux index
    int     k;              // Aux index


    // If the input has more than one ket the metropolis method can not sample it.
    if(istate->nket>1) cout << "metropolis warning!: Multiple ket input state. All kets are ignored except the first one" << endl;

    // Set up dimensions of the probel,
    s=1.0;
    nph=0;
    nlevel=qoc->nlevel;
    for(i=0;i<nlevel;i++){
        nph=nph+istate->ket[0][i];
        s=s*(double)factorial(istate->ket[0][i]);
    }

    //Reserve memory
    ilist=new int[nph]();
    r=new int[nph]();
    Ust.resize(nph,nph);
    obin=new p_bin(nph,nlevel,mem);

    //Initialize input configuration.
    k=0;
    for(i=0;i<nlevel;i++){
        for(j=0;j<istate->ket[0][i];j++){
            ilist[k]=i;
            k=k+1;
        }
    }

    // Configure the flags of the chosen method.
    switch (method){
        case 0: // f Classical
            classic=true;
            gral=true;
            uniform=false;
            break;
        case 1: // g Uniform
            classic=false;
            gral=true;
            uniform=true;
            break;
        case 2: // g Classical
            classic=false;
            gral=true;
            uniform=false;
            break;
        case 3: // f Classical (Restricted)
            classic=true;
            gral=false;
            uniform=false;
            break;
        case 4: // g Uniform (Restricted)
            classic=false;
            gral=false;
            uniform=true;
            break;
        case 5: // g Classical  (Restricted)
            classic=false;
            gral=false;
            uniform=false;
            break;
        default:
            cout << "Metropolis error: No recognized method." << endl;
            return {obin,0.0};
            break;
    }


    // Initialize loop variables and counters
    isample=0;
    istored=0;
    Neff=0;
    p_old=1.0;
    pc_old=1.0;
    while(istored<N){
        // Generate classically distributed sample
        tie(occ,pc)=classical_sample(ilist,nph,gral,uniform,qoc);


        // Obtain acceptance probability
        if(classic==false){
            // Init variables
            for(i=0;i<nph;i++) r[i]=0;
            k=0;
            t=1.0;
            for(i=0;i<nlevel;i++){
                t=t*(double)factorial(occ[i]);
                for(j=0;j<occ[i];j++){
                    r[k]=i;
                    k=k+1;

                }
            }

            // Generate Ust
            icol=0;
            for(ilin=0;ilin<nph;ilin++){
                irow=0;
                for(ilout=0;ilout<nph;ilout++){
                    Ust(irow,icol)=qoc->circmtx(r[ilout],ilist[ilin]);
                    irow=irow+1;
                }
                icol=icol+1;
            }

            // Calculate probability
            p=pow(abs(glynn(Ust)),2)/(s*t);

            // Calculate acceptance probability
            if(uniform==false) T=min(1.0,(p*pc_old)/(pc*p_old));
            else T=min(1.0,p/p_old);
        }else{
            T=1.0;
        }


        // Accept sample with probability T
        if(urand()<T){
            p_old=p;
            pc_old=pc;
            isample=isample+1;

            // Store sample if it is not burn or thined
            if((isample>=Nburn)&&(isample%Nthin==0)){

                index=obin->add_count(occ);
                istored=istored+1;

                if(index<0){
                    cout << "Sample: Warning! Sampling canceled because the memory limit has been exceeded.  Increase *mem* for more memory." << endl;

                    // Free memory
                    delete[] r;
                    delete[] ilist;
                    delete [] occ;

                    // Return partial calculation.
                    return {obin,(double) isample/(double)Neff};
                }
            }
        }

        Neff=Neff+1;
        delete[] occ;
    }

    // Free memory
    delete[] r;
    delete[] ilist;

    // Return output
    return {obin,(double) isample/(double)Neff};
}


//-----------------------------------------------
//
//  Obtain a sample from a circuit assuming photons
//  are classical particles.
//
//-----------------------------------------------
tuple<int*, double> simulator::classical_sample(int *ilist, int nph, bool gral, bool uniform, qocircuit *qoc ){
//  int *ilist;       // List of input photons.
//  bool gral;        // Levels may have any number of photons true='Yes'/False='No'
//  bool uniform;     // Uniform or classical distribution.   True='uniform'/False='Classical'
//  qocircuit *qoc;   // Circuit being sampled.
//  Variables
    int    nlevel;    // Number of levels
    double t;         // n1!n2!..nl!
    double c;         // nph!
    double pc;        // Classical probability
    double auxpc;     // Auxiliary variable to calculate classical probability.
    int   *occ;       // Occupation
    string perm;      // Permutations of the output photon list.
//  Auxiliary index
    int    i;         // Aux index
    int    j;         // Aux index
    int    k;         // Aux index


    // Initialize variables
    pc=0.0;
    c=(double)factorial(nph);
    nlevel=qoc->nlevel;
    perm.resize(nph);


    // Generate state.
    // Accept only with probability pc
    while(urand()>=pc){
        // Generate uniform distribution
        if(gral==true) tie(occ,t)=uniform_general(nph,qoc);
        else tie(occ,t)=uniform_restricted(nph,qoc);


        // Generate classical distribution of the circuit if asked.
        if(uniform==false){
            k=0;
            for(i=0;i<nlevel;i++){
                for(j=0;j<occ[i];j++){
                    perm[k]=i;
                    k=k+1;
                }
            }

            pc=0.0;
            sort(perm.begin(),perm.end());
            do{
                auxpc=1.0;
                for(i=0;i<nph;i++){
                    auxpc=auxpc*pow(abs(qoc->circmtx(perm[i],ilist[i])),2);
                }
                pc=pc+auxpc;
            }while (next_permutation(perm.begin(), perm.end()));
            pc=pc*t/c;
        }else{
            pc=1.0;
        }
    }

    // Return sample.
    return {occ,pc};
}


//-----------------------------------------------
//
//  Generation of random nlevel states of nph photons.
//  ( Uniform distribution)
//
//-----------------------------------------------
tuple<int*, double> simulator::uniform_general(int nph, qocircuit *qoc){
//  int nph;         // Number of photons
//  int qocircuit;   // Circuit being sampled
//  Variables
    int nlevel;      // Number of levels
    double t;        // n1!n2!..nl!
    double p;        // Acceptance probability
    int *occ;        // Occupation
//  Auxiliary index
    int i;           // Aux index
    int l;           // Aux index


    // Initialize variables and reserve memory
    nlevel=qoc->nlevel;
    occ=new int[nlevel]();

    // Generate state.
    // Accept only with probability p.
    p=0.0;
    while(urand()>=p){
        // Init occupation
        for(i=0;i<nlevel;i++) occ[i]=0;
        // Generate state
        for(i=0;i<nph;i++){
            l=floor(nlevel*urand());
            occ[l]=occ[l]+1;
        }

        // Calculate acceptance probability.
        // c/t is the number of times this occupation will be generated
        // in excess over the uniform distribution. Therefore we have to
        // reject it with the inverse of this value.
        t=1.0;
        for(i=0;i<nlevel;i++) t=t*(double)factorial(occ[i]);
        p=t/(double)factorial(nph);

    }

    // Return result
    return {occ,t};
}


//-----------------------------------------------
//
//  Generation of random nlevel states of nph photons.
//  Maximum of one photon by level
//  ( Uniform distribution)
//
//-----------------------------------------------
tuple<int*, double> simulator::uniform_restricted(int nph, qocircuit *qoc){
//  int nph;         // Number of photons
//  int qocircuit;   // Circuit being sampled
//  Variables
    int iph;         // Photons allocated
    int *occ;        // Occupation
//  Auxiliary index
    int l;           // Aux index


    // Reserve memory
    occ=new int[qoc->nlevel]();

    // Generate a random state. (Max one photon by level)
    iph=0;
    while(iph<nph){
        l=floor(qoc->nlevel*urand());
        if(occ[l]==0){
            occ[l]=occ[l]+1;
            iph=iph+1;
        }
    }

    // Return result
    return {occ,(double)factorial(nph)};
}

