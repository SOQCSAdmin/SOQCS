//======================================================================================================
// File dmat.cpp
//
// DENSITY MATRIX LIBRARY
//
// Copyright © 2022 National University of Ireland Maynooth, Maynooth University. All rights reserved.
// The contents of this file are subject to the licence terms detailed in LICENCE.TXT available in the
// root directory of this source tree. Use of the source code in this file is only permitted under the
// terms of the licence in LICENCE.TXT.
//======================================================================================================
#include "dmat.h"


//----------------------------------------
//
//  Create density matrix.
//  The quantity of reserved memory for
//  the matrix is set by default.
//
//----------------------------------------
dmatrix::dmatrix(){


    create_dmtx(DEFMATDIM);
}


//----------------------------------------
//
//  Create a density matrix object.
//  The maximum dimension of the matrix
//  is set by its row number.
//
//----------------------------------------
dmatrix::dmatrix(int i_mem){
//  int i_mem;       // Number of rows of the density matrix.


    create_dmtx(i_mem);
}


//---------------------------------------------------
//
//  Auxiliary private method to create a density matrix.
//
//---------------------------------------------------
void dmatrix:: create_dmtx(int i_mem){
//  int i_mem;       // Number of rows of the density matrix.


    // Initialize configuratin
    N=0;
    mem=i_mem;

    // Reserve memory for density matrix
    dens.setZero(mem,mem);

    // Crete placeholder dictionary.
    // The good one will be created with the first element.
    // Important data extracted from the first added state.
    dicc=new ket_list(1,1);
}


//----------------------------------------
//
//  Destroy density matrix
//
//----------------------------------------
dmatrix:: ~dmatrix(){


    // Deleted dictionary
    delete dicc;
}


//----------------------------------------
//
//  Copy a density matrix
//
//----------------------------------------
dmatrix *dmatrix:: clone(){
//  Variables
    dmatrix *aux;      // New copy of this density matrix.


    // Initialize configuratin
    aux=new dmatrix(mem);
    aux->N=N;

    // Reserve memory for density matrix
    aux->dens=dens;

    // Crete placeholder dictionary.
    // The good one will be created with the first element.
    // Important data extracted from the first added state.
    aux->dicc=dicc->clone();

    return aux;
}


//---------------------------------------------------
//
//  Clear a density matrix
//
//---------------------------------------------------
void dmatrix:: clear(){


    // Initialize configuration
    N=0;

    // Reserve memory for density matrix
    dens.setZero(mem,mem);

    // Crete placeholder dictionary.
    // The good one will be created with the first element.
    // Important data extracted from the first added state.
    delete dicc;
    dicc=new ket_list(1,1);
}


//----------------------------------------
//
// Sum two compatible density matrix
//
//----------------------------------------
void dmatrix:: add(dmatrix *addm){
//  dmatrix *addm                  // Matrix to be added to the present one
//  Variables
    double Nm;                     // Normalization coefficient between matrix
//  Auxiliary index
    int   irow;                    // Index of rows
    int   icol;                    // Index of columns
    int   k;                       // Aux index
    int   l;                       // Aux index


    Nm=(double)N/(double) addm->N;
    // For each ket in the state we add it to the dictionary and
    // assign a row to it if it didn't exist before.  We obtain its
    // row otherwise.
    for(k=0;k<addm->dicc->nket;k++){
        irow=dicc->add_ket(addm->dicc->ket[k]);
        // For each ket in the state we add it to the dictionary and
        // assign a column to it if it didn't exist before.  We obtain its
        // column otherwise.
        for(l=0;l<addm->dicc->nket;l++){
            icol=dicc->add_ket(addm->dicc->ket[l]);
            dens(irow,icol)=dens(irow,icol)+addm->dens(k,l)*Nm;
        }
    }
}


//----------------------------------------
//
// Calculates trace of the density matrix
//
//----------------------------------------
double dmatrix::trace(){
//  Variables
    double tr;      // Trace
//  Auxilaty index
    int    i;       // Aux index


    tr=0.0;
    for (i=0;i<dicc->nket;i++) tr=tr+dens(i,i);
    return tr;
}


//----------------------------------------
//
// Normalize to trace=1 density matrix
//
//----------------------------------------
void dmatrix::normalize(){


    N=1;
    dens=dens/trace();
}


//----------------------------------------
//
//  Print the matrix
//  (numerically)
//
//----------------------------------------
void dmatrix::prnt_mtx(){


    aux_prnt_mtx(DEFFORMAT,xcut,nullptr);
}


//----------------------------------------
//
//  Print the matrix
//  (numerically)
//
//----------------------------------------
void dmatrix::prnt_mtx(double thresh){
//  double thresh    // Threshold value to print


    aux_prnt_mtx(DEFFORMAT,thresh,nullptr);
}


//----------------------------------------
//
//  Print the matrix
// (human readable form)
//
//----------------------------------------
void dmatrix::prnt_mtx(int format, double thresh, qocircuit *qoc){
//  int format       // Format of the index ket to be printed
//  double thresh    // Threshold value to print
//  qocircuit *qoc   // Circuit to which the density matrix is related (to translate from numeric to human form)


    aux_prnt_mtx(format, thresh, qoc);
}


//----------------------------------------
//
//  Print the matrix
// (human readable form)
//  Uses qodev instead of qocircuit
//
//----------------------------------------
void dmatrix::prnt_mtx(int format, double thresh, qodev *dev){
//  int format       // Format of the index ket to be printed
//  double thresh    // Threshold value to print
//  qodev *dev       // Device to which the density matrix is related (to translate from numeric to human form)


    if(dev > (qodev *)nullptr) aux_prnt_mtx(format, thresh, dev->circ);
    else aux_prnt_mtx(format, thresh, nullptr);
}

//----------------------------------------
//
//  Auxiliary print matrix method
//
//----------------------------------------
void dmatrix::aux_prnt_mtx(int format, double thresh, qocircuit *qoc){
//  int format       // Format of the index ket to be printed
//  double thresh    // Threshold value to print
//  qocircuit *qoc   // Circuit to which the density matrix is related (to translate from numeric to human form)
//  Variables
    int  nbase;      // nbase->dicc->nket Number of dictionary entries
    vecd rowsum;     // Greater than zero if one row element is not zero. Zero otherwise.
    vecd colsum;     // Greater than zero if one column element is not zero. Zero otherwise.
//  Auxiliary index
    int  i;          // Aux index
    int  j;          // Aux index


    nbase=dicc->nket;

    // Check which rows are non zero
    rowsum.resize(nbase);
    for(i=0;i<nbase;i++){
        rowsum(i)=0.0;
        for(j=0;j<nbase;j++) rowsum(i)=rowsum(i)+abs(dens(i,j));
    }

    // Check which columns are non zero
    colsum.resize(nbase);
    for(i=0;i<nbase;i++){
        colsum(i)=0.0;
        for(j=0;j<nbase;j++) colsum(i)=colsum(i)+abs(dens(j,i));
    }

   // Print non-zero rows and columns of the density matrix (to avoid screen cluttering)
    for(i=0;i<nbase;i++){
        if(rowsum(i)>thresh){
            dicc->prnt_ket(i,format,qoc);
            cout << " ";
            for(j=0;j<nbase;j++){
                if(colsum(j)>thresh){
                    cout  << setw(DEFWIDTH+1) << std::fixed << std::setprecision(4) << dens(i,j)/N << RESET <<" ";
                }
            }
            cout << endl;
        }
    }
    cout << endl;
}


//----------------------------------------
//
// Prints the diagonal elements of a density matrix
//
//----------------------------------------
void dmatrix::prnt_results(){

    prnt_results(DEFFORMAT,nullptr);
}


//----------------------------------------
//
// Prints the diagonal elements of a density matrix
//
//----------------------------------------
void dmatrix::prnt_results(int format,qocircuit *qoc){
//  int format       // Format of the index ket to be printed
//  qocircuit *qoc   // Circuit to which the density matrix is related (to translate from numeric to human form)
//  Auxiliary index
    int i;           // Aux index


   // Print diagonal elements of the density matrix.
    for(i=0;i<dicc->nket;i++){
        dicc->prnt_ket(i,format,qoc);
        cout << ": ";
        cout << std::fixed << std::setprecision(4) << dens(i,i)/N <<endl;

    }
}


//------------------------------------------------
//
//  Gets the probability of the event defined by def
//
//------------------------------------------------
double dmatrix::get_result(mati def,qocircuit *qoc){
//  mati       def;                 // Definition of the bra state.
//  qocircuit *qoc;                 // Circuit where the detector is placed
//  Variables
    state *bra;                     // Bra state created from the definition def
//  Dictionary variables
    int    irow;                    // Index of rows


    // Create bra state
    bra= new state(dicc->nlevel,1,dicc->vis);
    bra->add_term(1.0,def,qoc);
    irow=dicc->find_ket(bra->ket[0]);
    delete bra;

    // Return density matrix value if present.
    if(irow>=0) return dens(irow,irow);

    return 0.0;
}


//------------------------------------------------
//
//  Converts the diagonal elements of a density matrix
//  into a probability bin
//
//------------------------------------------------
p_bin *dmatrix::get_pbin(){
//  Variables
    p_bin *aux;     // Auxiliary probability bin.
//  Auxiliary index
    int    i;       // Aux index
    int    j;       // Aux index


    // Create the set of probability bin
    aux=new p_bin(dicc->nlevel,dicc->maxket,dicc->vis);

    // Store the values
    for(i=0;i<dicc->maxket;i++){
        j=aux->add_count(dicc->ket[i]);
        aux->p[j]=dens(i,i);
    }
    aux->N=N;

    // Return the bins
    return aux;
}


//----------------------------------------
//
// Calculate the Fidelity of a density matrix
// with respect to a state.
//
//----------------------------------------
double dmatrix::fidelity(state* input){
//  state *input                   // Reference state
//  Variables
    double F;                      // Fidelity
    int    irow;                   // Index of rows
    int    icol;                   // Index of columns
//  Auxiliary index
    int    k;                      // Aux index
    int    l;                      // Aux index


    F=0.0;
    // Search for the corresponding row and column of each ket in the
    // density matrix
    for(k=0;k<input->nket;k++){
        irow=dicc->find_ket(input->ket[k]);

        for(l=0;l<input->nket;l++){
            icol=dicc->find_ket(input->ket[l]);

            // IF found the we calculate the corresponding fidelity measurement.
            if((irow>=0)&&(icol>=0)){
                F=F+(real(conj(input->ampl[k])*dens(irow,icol)*input->ampl[l]));
            }
        }
    }
    return F;
}


//----------------------------------------------------
//
//  Adds conditional detection defined in the circuit.
//
//----------------------------------------------------
void dmatrix:: add_state(state* input, qocircuit *qoc){
//  state     *in_input;            // New state to be added to the density matrix
//  qocircuit *qoc;                 // Circuit where the detector is placed
//  Variables
    int  fchloss;                   // First loss channel
    int  nchloss;                   // Number of loss channels
    int  nignore;                   // Number of ignored channels
    int  nchtotal;                  // Total number of channels to be traced out
    veci chlist;                    // List of channels to trace out
//  Auxiliary index
    int  i;                         // Aux index


    // Init variables
    if(qoc->losses==1) nchloss=qoc->nch/2;
    else nchloss=0;
    nignore=qoc->nignored;
    nchtotal=nchloss+nignore;
    fchloss=qoc->nch/2;

    // Create list
    chlist.resize(nchtotal);
    for(i=0;i<nchloss;i++) chlist(i)=fchloss+i;
    for(i=0;i<nignore;i++) chlist(i+nchloss)=qoc->ch_ignored(i);

    // Obtain reduced state
    this->add_reduced_state(qoc->ncond,qoc->det_def,chlist,input,qoc);
    N++;
}


//----------------------------------------------------
//
//  Adds conditional detection defined explicitly.
//  Channels in chlist are traced out.
//
//----------------------------------------------------
void dmatrix::add_reduced_state(int ndec,mati def,veci chlist,state* input, qocircuit *qoc){
//  int        ndec;                // Number of detector <= number of channels
//  mati       def;                 // Definition of the conditions of detection
//  veci       chlist;              // List of channels to be traced out
//  state     *input;               // New state to be added to the density matrix
//  qocircuit *qoc;                 // Circuit where the detector is placed
//  Variables
    int        nph;                 // Max number of photons to be considered in the trace
    int        ntrace;              // Number of photons to be traced in the trace.
    int        nchtotal;            // Total number of channels to be traced out
    int        nsel;                // Number of post-selection conditions
    veci       chrem;               // In which channel are the traced photons
    veci       occ;                 // Occupation of the channels with losses
    mati       select;              // Post-selection condition of the losses
    mati       selectfull;          // Full post-selection condition
//  Hash table
    int       *keysel;              // Key of the projector
    long int   selvalue;            // Projector definition index
    thash      selhash;             // Definitions hash table
    thash::const_iterator vselhash; // Row hash value.
//  Auxiliary index
    int         i;                  // Aux index
    int         j;                  // Aux index
    int         k;                  // Aux index



//  Init variables
    nchtotal=chlist.size();
    if(nchtotal>0) nph=maxnph;
    else nph=0;

    // Reserve memory
    keysel=new int(nchtotal);
    // For any number of photon loss
    for(ntrace=0;ntrace<=nph;ntrace++){
        selhash.clear();
        nsel=0;
        // We initialize in which channel is each photon
        chrem.setZero(ntrace+1);

        // We post_select all the configurations of ntrace photons
        while(chrem(ntrace)==0){
            // First calculate the occupations
            occ.setZero(nchtotal);
            for(i=0;i<ntrace;i++) occ(chrem(i))=occ(chrem(i))+1;

            select.resize(2,nchtotal);
            for(k=0;k<nchtotal;k++){
                select(0,k)=chlist(k);
                select(1,k)=occ(k);
                keysel[k]=occ(k);
            }

            // Do post-selection but check if repeated
            selvalue=hashval(keysel,nchtotal,ntrace);
            vselhash=selhash.find(selvalue);
            // If the projector has not been created before
            if(vselhash==selhash.end()){
                // We updated the entry in the corresponding hash table
                selhash[selvalue]=nsel;
                nsel=nsel+1;
                selectfull.resize(2,select.cols()+ndec);
                k=0;
                for(i=0;i<ndec;i++){
                    selectfull(0,k)=def(0,i);
                    selectfull(1,k)=def(1,i);
                    k=k+1;
                }
                for(i=0;i<select.cols();i++){
                    selectfull(0,k)=select(0,i);
                    selectfull(1,k)=select(1,i);
                    k=k+1;
                }

                this->add_state_cond(selectfull.cols(),selectfull,input,qoc);
            }


            // Update in which channel is each photon to sweep them all.
            chrem(0)=chrem(0)+1;
            for(j=0;j<ntrace;j++){
                if(chrem(j)>=nchtotal){
                    chrem(j)=0;
                    chrem(j+1)=chrem(j+1)+1;
                }
            }

        }
    }

    //Free memory
    delete[] keysel;
}


//----------------------------------------
//
 // Adds a conditional detection defined explicitly
//
//----------------------------------------
void dmatrix:: add_state_cond(int ndec, mati def,state *in_state,qocircuit *qoc){
//  int    ndec;                     // Number of detector <= number of channels
//  mati   def;                      // Definition of the conditions of detection
//  state *in_state;                 // New state to be added to the density matrix
//  qocircuit *qoc;                  // Circuit where the detector is placed
//  Variables
    int    nph;                      // Number of photons
    int    maxch;                    // Larger channel number
    int    selbase;                  // Base of the projector definition hash table
    int    nempty;                   // Number of channels with zero fotons
    int    nentry;                   // Number of entries in the projector definition
    int    nprj;                     // Number of projectors
    int    iph;                      // Index of photons
    int    ich;                      // Channel index
    int    im;                       // Mode index
    int    is;                       // "Time" index
    int   *tim;                      // "Time" configuration vector
    int   *pol;                      // Polarization configuration vector
    veci   ch;                       // List of channels by photons
    hterm  select;                   // Post-selection condition
    state *newstate;                 // New post-selected state
    projector *prj;                  // Projector
//  Hash table
    int    key[3];                   // Key of the projector definition
    int   *keyprj;                   // Key of the projector
    long   int selvalue;             // Projector definition index
    long   int prjvalue;             // Projector index
    thash  selhash;                  // Definitions hash table
    thash  prjhash;                  // Projectors hash table.
    thash::const_iterator vselhash;  // Row hash value.
    thash::const_iterator vprjhash;  // Row hash value.
//  Auxiliary index
    int    i;                        // Aux index.
    int    j;                        // Aux index
    int    k;                        // Aux index
    int    l;                        // Aux index


    if(ndec>0){
    // If the number of conditional detector is larger than zero the
    // corresponding post-selections are performed
        // Calculate:
        // How many photons.
        // How many empty channels
        // Maximum number of photons in channel
        nph=0;
        nempty=0;
        for(ich=0;ich<ndec;ich++){
            nph=nph+def(1,ich);
            maxch=max(maxch,def(0,ich));
            if(def(1,ich)==0) nempty=nempty+1;
        }
        selbase=max(maxch,max(qoc->nm,qoc->ns));

        // Create a list for every photon in which channel they are.
        // Empty channels that are post-selected to zero are appended at the endl.
        ch.resize(nph+nempty);
        k=0;
        l=0;
        for(ich=0;ich<ndec;ich++){
            for(iph=0;iph<def(1,ich);iph++){
                ch[k]=def(0,ich);
                k++;
            }

            if(def(1,ich)==0){
                ch(nph+l)=def(0,ich);
                l++;
            }
        }


        // Initialize projector aux variables.
        nprj=0;
        keyprj=new int[ndec*qoc->nm*qoc->ns]();

        // From the detector definition of conditional detection create
        // all the projections necessary considering all the possible
        // polarizations and times.

        // LOOP Polarization¡
        pol=new int[nph+1]();
        while(pol[nph]==0){

            // LOOP "Time"
            tim=new int[nph+1]();
            while(tim[nph]==0){
                // Create Projector definition
                nentry=0;
                selhash.clear();
                select.setZero(4,ndec*qoc->nm*qoc->ns);
                for(iph=0;iph<(nph+nempty);iph++){
                for(im=0;im<qoc->nm;im++){
                for(is=0;is<qoc->ns;is++){
                    key[0]=ch(iph);
                    key[1]=im;
                    key[2]=is;

                    // We have to check if for those photons the entry in select
                    // already exists. If exist we update otherwise it has to be
                    // updated. Otherwise we create it.
                    selvalue=hashval(key,3,selbase);
                    vselhash=selhash.find(selvalue);
                    if(vselhash==selhash.end()){
                        k=nentry;
                        selhash[selvalue]=nentry;
                        nentry=nentry+1;
                    }else{
                        k=vselhash->second;
                    }

                    select(0,k)=ch(iph);
                    select(1,k)=im;
                    select(2,k)=is;

                    if((iph<nph)&&(im==pol[iph])&&(is==tim[iph])) select(3,k)=select(3,k)+1;
                }}}

                // Create projector.
                // Check that we have not generated that projector before.
                // Photons are indistinguishable therefore in the way this loop
                // is built repetitions may appear. If it has not been generated
                // before the projector is created
                for(i=0;i<(ndec*qoc->nm*qoc->ns);i++){
                    keyprj[i]=select(3,i);
                }
                prjvalue=hashval(keyprj,ndec*qoc->nm*qoc->ns);
                vprjhash=prjhash.find(prjvalue);
                // If the projector has not been created before
                if(vprjhash==prjhash.end()){
                    // We updated the entry in the corresponding hash table
                    prjhash[prjvalue]=nprj;
                    nprj=nprj+1;

                    // We create the projector
                    prj=new projector(qoc->num_levels(),1,in_state->vis); //<<<<Check this
                    prj->add_term(1.0,select,qoc);

                    // Apply the projector to the input state and store the result
                    newstate=in_state->post_selection(prj);
                    this->sum_state(newstate);

                    // Free memory
                    delete newstate;
                    delete prj;
                }

                // ADVANCE "Time" in LOOP
                tim[0]=tim[0]+1;
                for(j=0;j<nph;j++){
                    if(tim[j]>=qoc->ns){
                        tim[j]=0;
                        tim[j+1]=tim[j+1]+1;
                    }
                }
            }

            // ADVANCE Polarization in LOOP
            pol[0]=pol[0]+1;
            for(j=0;j<nph;j++){
                if(pol[j]>=qoc->nm){
                    pol[j]=0;
                    pol[j+1]=pol[j+1]+1;
                }
            }
            // At the end of the loop free memory
            delete[] tim;
        }
        // At the end of the loop free memory
        delete[] pol;

        // Free memory;
        delete[] keyprj;
        selhash.clear();
        prjhash.clear();

    }else{
    // If the number of conditional detector is zero
    // the the input state is just added to the matrix
        sum_state(in_state);
    }

}


//----------------------------------------
//
//  Adds a state
//
//----------------------------------------
void dmatrix:: sum_state(state *newstate){
//  state *newstate                // New state to be added to the density matrix.
//  Variables
    double P;                      // Probability
//  Auxiliary index
    int    irow;                   // Index of rows
    int    icol;                   // Index of columns
    int    k;                      // Aux index
    int    l;                      // Aux index


    // For each ket in the state we add it to the dicctionary and
    // assign a row to it if it didn't exist before.  We obtain its
    // row otherwise.
    if(dicc->nket==0){
        delete dicc;
        dicc=new ket_list(newstate->nlevel,mem,newstate->vis);
    }

    for(k=0;k<newstate->nket;k++){
        irow=dicc->add_ket(newstate->ket[k]);
        for(l=0;l<newstate->nket;l++){
            icol=dicc->add_ket(newstate->ket[l]);

            P=(real(conj(newstate->ampl[k])*newstate->ampl[l]));

            // We update the density coefficient at the proper row and column
            if(abs(P)>xcut) dens(irow,icol)=dens(irow,icol)+P;
        }
    }
}


//----------------------------------------
//
//  Calculate measurement
//
//----------------------------------------
dmatrix *dmatrix::calc_measure(qocircuit *qoc){
//  qocircuit *qoc;         // Circuit to which the measure is referred
//  Variables
    dmatrix *aux;           // Auxiliary matrix


    // Calculate measurement
    if(qoc->ns>1){
        switch(qoc->timed){
            case 0:  // Counter
                aux=get_counts(qoc->ns,qoc);
                break;
            case 1:  // Clock: Partial trace
                aux=partial_trace(qoc->emitted->pack_def,qoc);
                break;
            case 3:  // Clock: Manual mode
                aux=partial_trace(qoc->emitted->pack_def,qoc);
                break;
            case 2:  // Clock + Spectrum: Full
                aux=this->clone();
                break;
            default: // Default: Error
                cout << "Dmatrix calc_measure: Error. Timed has a non valid value" << endl;
                exit(0);
                break;
        }

    }else{
        aux=this->clone();  // If there is only one wave-packet the density matrix is just copied
    }
    return aux;
}


//-------------------------------------------------------------------------------------------------
//
// Calculates the partial trace of a matrix using a packet definition of the degrees of freedom
//
//-------------------------------------------------------------------------------------------------
dmatrix *dmatrix::partial_trace(mati pack_def,qocircuit *qoc){
//  mati pack_def       // Packet definition that relates the packet number with the numbers of the additional degrees of freedom
//  qocircuit *qoc      // Circuit to which the density matrix is referred.
//  Variables.
    mati pack_idx;      // Packet index. Degrees of freedom of each packet indexed by the packet number.
    int *rowocc;        // New dictionary ket occupation of the new partial traced density matrix (referred to rows).
    int *colocc;        // New dictionary ket occupation of the new partial traced density matrix (referred to columns).
    dmatrix *newdmat;   // New partial traced density matrix.
//  Transformation variables
    int  ch;            // Channel
    int  pol;           // Polarization
    int  ns;            // New wavepacket
    int  newk;          // Level index
    int  istore;        // Level to store index
    int  auxstore;      // Aux index with the un-traslated level to store index.
    int  irow;          // Density matrix row index
    int  icol;          // Density matrix col index
    // Note that the division between rows and columns is conceptual. There is only one hash table for both.
    // Auxiliary index
    int  i;             // Aux index
    int  j;             // Aux index
    int  k;             // Aux index


    // Reserve memory for new partial matrix
    newdmat=new dmatrix(mem);
    // Update level index to the ones of the current matrix
    delete newdmat->dicc;
    newdmat->dicc=new ket_list(dicc->nlevel,mem,dicc->vis);
    // Reorganize  the information of the pac
    if(pack_def.rows()>2){                              // Some backward compatibility
        pack_idx=create_packet_idx(pack_def);
    }else{
        pack_idx=pack_def;
    }

    // For each element of the base int rows and columns
    for(i=0;i<dicc->nket;i++){
    for(j=0;j<dicc->nket;j++){
        // The result of the traced degrees of freedom can be only 0 or 1.
        // If it is one we store the remaining degrees of freedom not traced.
        // That means to assign to the base elements a new label/number that identifies
        // uniquely those degrees of freedom.
        if(ketcompatible(i,j,pack_idx,qoc)==1){
            // A new row and column dictionary element are created
            rowocc=new int[dicc->nlevel]();
            colocc=new int[dicc->nlevel]();
            // We process each level to determine the index of the remaining degrees of freedom
            for(k=0;k<dicc->nlevel;k++){
                //Remember that the order/label of a level may have changed after post-selection.
                // We need their original ordering in the circuit index.
                newk=dicc->vis[k];

                // We recover the channel and polarization of the level.
                ch=qoc->idx[newk].ch;
                pol=qoc->idx[newk].m;
                // We change the label for a new one representing the remaining degrees of freedom.
                // Note that doing so we are re-interpreting the meaning of the levels. Therefore
                // we have new occupation vectors.

                if(qoc->idx[newk].s<pack_idx.cols()){ // This is a safety measure. Sometimes we don't use all labels.The remaining ones are zero
                    // We select the new label following the premises of our packet definitions.
                    // Each packet is labeled by a number, but is composed of different degrees of freedom also labeled by numbers.
                    // Here we change the label by the remaining one not traced out.
                    ns=pack_idx(0,qoc->idx[newk].s);
                    // Wich level correspond with that wavefunction number. Here the reinterpretation.
                    auxstore=qoc->i_idx[ch][pol][ns];
                    istore=0;
                    while(auxstore!=newdmat->dicc->vis[istore]) istore++; // Don't like it very much. But those list are not very large and this operation is not looped.
                    // We compute the new dictionary element
                    rowocc[istore]=rowocc[istore]+dicc->ket[i][k];
                    colocc[istore]=colocc[istore]+dicc->ket[j][k];
                }

            }

            // Now we store the dictionary element and the density matrix entry
            // Check rows. If present recover the index otherwise create a new one
            irow=newdmat->dicc->add_ket(rowocc);
            icol=newdmat->dicc->add_ket(colocc);
            delete[] rowocc;
            delete[] colocc;

            // Store matrix element for the new dictionary entries
            newdmat->dens(irow,icol)=newdmat->dens(irow,icol)+dens(i,j);
        }
    }}

    // Copy additional information (normalization information)
    newdmat->N=N;

    // Return new matrix
    return newdmat;
}


//-------------------------------------------------------------------------------------------------
//
// Returns a density matrix where the entries are independent of the wave packet degrees of freedom.
// This is independent of time and frequency
//
//-------------------------------------------------------------------------------------------------
dmatrix *dmatrix::get_counts(int npack,qocircuit *qoc){
//  int npack;                 // Number of packet definitions.
//  qocircuit *qoc;            // Circuit to which the density matrix is referred.
//  Variables
    dmatrix   *aux;            // Auxiliary density matrix
    ket_list  *newdicc;        // New dictionary definition
    mati       new_pack_def;   // New packet definitions
//  Auxiliary index
    int        i;              // Aux index


    // Create the packet definition
    new_pack_def.resize(3,npack);
    for(i=0;i<npack;i++){
        new_pack_def(0,i)=i;
        new_pack_def(1,i)=0;
        new_pack_def(2,i)=i;
    }

    // Perform the partial trace
    aux=partial_trace(new_pack_def,qoc);

    // Update the dictionary "labels"
    newdicc=aux->dicc->remove_time(qoc);
    delete aux->dicc;
    aux->dicc=newdicc;

    return aux;
}


//-------------------------------------------------------------------------------------------------
//
// Ket "compatibility". Do they have the same instance of the degree of freedom that is going to be traced out.
// Auxiliary private function. Not intended for external use.
//
//-------------------------------------------------------------------------------------------------
int dmatrix::ketcompatible(int A, int B,mati pack_idx,qocircuit *qoc){
//  state *A           // State A
//  state *B           // State B
//  mati pack_idx      // Packet index of the degrees of freedom
//  qocircuit *qoc     // Circuit to which the density matrix is referred.
//  Variables
    veci    flistA;    // List of the index to be traced out in ket A
    veci    flistB;    // List of the index to be traced out in ket B
    int     maxA;      // Maximum number of elements in flist A
    int     maxB;      // Maximum number of elements in flist B
    int    compatible; // Are A an B compatible 0=No/1=Yes
    int    lev;        // Level of state A
    int    s;          // Level 2 wave packet
    int    w;          // Degree of freedom to be traced of level 1. Usually frequency.
    // Auxiliary index
    int    i;          // Aux index
    int    j;          // Aux index


    // Reserve memory
    flistA.resize(maxnph);
    flistB.resize(maxnph);

    // Initialize lists an variables
    maxA=0;
    maxB=0;
    flistA.setZero(maxnph);
    flistB.setZero(maxnph);

    // Create list of indexes to compatibility comparison
    for(i=0;i<dicc->nlevel;i++){
        lev=dicc->vis[i];
        s=qoc->idx[lev].s;
        if(s<pack_idx.cols()){
            w=pack_idx(1,s);
            for(j=0;j<dicc->ket[A][i];j++){
                flistA[maxA]=w;
                maxA=maxA+1;
            }

            for(j=0;j<dicc->ket[B][i];j++){
                flistB[maxB]=w;
                maxB=maxB+1;
            }
        }
    }

    // Comparison between the lists.
    // States are compatible if the list are compatible
    i=0;
    compatible=1;
    if(maxA!=maxB) compatible=0;
    while((compatible==1)&&(i<maxA)){
        if(flistA[i]!=flistB[i]) compatible=0;
        i=i+1;
    }

    // Return compatibility.
    return compatible;
}


//----------------------------------------------------
//
//  Adds a conditional detection defined in the
//  generalized metacircuit.
//
//----------------------------------------------------
void dmatrix::add_state(state *newrun, qodev *dev){
//  state *newrun;                  // New state to be added to the density matrix
//  qodev *dev;                     // Generalized metacircuit where the detector is placed


    add_state(newrun,dev->circ);
}


//----------------------------------------
//
//  Calculate measurement
//  Uses qodev instead of qocircuit
//
//----------------------------------------
dmatrix *dmatrix::calc_measure(qodev *dev){
//  qodev *dev;                     // Generalized metacircuit to which the measurements are referred


    return this->calc_measure(dev->circ);
}
